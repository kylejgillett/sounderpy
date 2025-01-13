### IMPORT SOFTWARE ###
#########################################################################################################
import time
import sys
from urllib.request import urlopen
from urllib.error import HTTPError
import warnings
import bs4
import pandas as pd
import numpy as np
from numpy import loadtxt
import metpy.calc as mpcalc
from metpy.units import units

from .calc import *



"""
    SOUNDERPY ACARS 'GET-DATA' FUNCTIONS  

    Purpose of module: 

    House functions for loading and parsing ACARS profiles and 
    ACARS profile data. Functions here are explicitly called by
    the user.


    (C) KYLE J GILLETT, UNIVERSITY OF NORTH DAKOTA, 2024
"""






# define acars data class
class acars_data:
    """
    - NOTE: this is a Python ``Class``, not a function like the tools above.
       - This ``Class`` sets up a 'connection' to the ACARS data dataset.
       - After setting up a 'connection' to the data, you can search for available profiles using the class's function, ``.list_profiles()``
       - Then you may select one of the listed profiles and use it as an argument for the class's function, ``.get_profile()``. See below.

       :param year: observation year
       :type year: str, required
       :param month: observation month
       :type month: str, required
       :param day: observation day
       :type day: str, required
       :param hour: observation hour
       :type hour: str, required
    """

    # init
    def __init__(self, year, month, day, hour):
        self.year = year
        self.hour = hour
        self.month = month
        self.day = day
        self.hour = hour






    ### LIST PROFILES ###
    #########################################################################################################
    def list_profiles(self):
        """
        Returns
        -------

        Return a list of strings that represents ACARS profiles for a given date and hour.
        """

        st = time.time()
        print(f'> LIST ACARS PROFILES FUNCTION\n  ---------------------------------')

        # SET UP OU DIRECTORY REF
        # SEARCH FOR DIRECTORY FOR THE USER-SPECIFIED DATE
        data_dir = f'https://sharp.weather.ou.edu//soundings//acars//{self.year}//{self.month}//{self.day}//{self.hour}'

        # ACCESS THE RAW WEBSITE HTML & FIND THE ID_TIME KEYS
        print(f"> AVAILABLE ACARS PROFILES FOR {self.year}-{self.month}-{self.day} {self.hour}Z...")
        # SET UP BEAUTIFUL SOUP TO PARSE HTML
        # THIS WORKS AS A SORT OF JERRY-RIGGED WAY
        # TO REVEAL ALL AVAILABLE ACARS PROFILES
        # FOR A GIVEN DATE/TIME
        body = urlopen(data_dir).read().decode("utf-8")
        soup = bs4.BeautifulSoup(body, features="html.parser")

        # ADD PROFILES TO A LIST
        profiles_list = []
        for link in soup.select('a[href$=".txt"]'):
            profiles_list.append(link.get("href")[0:8])

        print('> COMPLETE --------')
        elapsed_time = time.time() - st
        print('> RUNTIME:', time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))

        return profiles_list


    #########################################################################################################







    ### GET PROFILE DATA ###
    #########################################################################################################
    def get_profile(self, acars_profile, hush=False, clean_it=True):

        """
        Return a ``dict`` of 'cleaned up' ACARS observation profile data. Do so by selecting one of the profile string "IDs" listed by ``list_profiles()`` and pasting it as an argument in ``get_profile()``

        :param acars_profile: profile "ID"
        :type acars_profile: str, required
        :param hush: whether to 'hush' a read-out of thermodynamic and kinematic parameters when getting a data.
        :type hush: bool, optional, default is `False`
        :param clean_it: whether to return the raw_data object or a clean_data dict.
        :type clean_it: bool, optional, default is `True`
        """

        st = time.time()
        print(f'> ACARS DATA ACCESS FUNCTION\n  ---------------------------------')

        # SET UP OU DIR REF TO THE SPECIFIC PROFILE
        profile_url = f'https://sharp.weather.ou.edu//soundings//acars//{self.year}//{self.month}//{self.day}//{self.hour}//{acars_profile}.txt'

        # SEPARATE DATA BETWEEN HEADER AND ACTUAL DATA
        try:
            raw_data = loadtxt(urlopen(profile_url).readlines()[6:-1], dtype='str', comments="%", unpack=True)
            header = loadtxt(urlopen(profile_url).readlines()[0:3], dtype='str', comments="%", unpack=True)
        except HTTPError as err:
            if err.code == 404:
                sys.exit('! ERROR ! -- Invalid profile OR profile does not exist. Try again with a valid profile (ex: BNA_2320)')
            else:
                raise

        if clean_it:
            # PARSE DATE INFO FROM OU FILE
            year = f'20{header[1][0:2]}'
            month = header[1][2:4]
            day = header[1][4:6]
            hour = f'{header[1][7:9]}:{header[1][9:11]}'

            # PARSE PROFILE DATA FROM OU FILE IN DICT
            new_keys = ['p', 'z', 'T', 'Td', 'u', 'v']
            units_list = ['hPa', 'meter', 'degC', 'degC']
            clean_data = {}
            for new_key, idx, unit in zip(new_keys, range(0, 4), units_list):
                clean_data[new_key] = np.array([float(ele) for ele in [ele[0:-1] for ele in raw_data[idx]]]) * units(unit)

            clean_data['u'], clean_data['v'] = mpcalc.wind_components(
                [float(ele) for ele in [ele[0:-1] for ele in raw_data[5]]] * units.kts,
                [float(ele) for ele in [ele[0:-1] for ele in raw_data[4]]] * units.deg)

            # GET AIRPORT INFO FROM GITHUB AIRPORTS.CSV
            airports_csv = pd.read_csv(f'https://raw.githubusercontent.com/kylejgillett/sounderpy/main/src/AIRPORTS.csv',
                                       skiprows=7, skipinitialspace=True)
            where = [np.where(airports_csv['IATA'].str.contains(header[0], na=False, case=True))[0]][0][0]

            # ADD AIRPORT DATA INTO DICT
            keys = ['Name', 'City', 'Country', 'Latitude', 'Longitude', 'Altitude']
            airport_info = []
            for key in keys:
                airport_info.append(airports_csv[key][where])

            # add clean_data site_info dict
            clean_data['site_info'] = {
                'site-id': header[0],
                'site-name': airport_info[0],
                'site-lctn': airport_info[2],
                'site-latlon': [np.round(airport_info[3], 2), np.round(airport_info[4], 2)],
                'site-elv': str(int(airport_info[5])),
                'source': 'ACARS OBSERVED AIRCRAFT PROFILE',
                'model': 'no-model',
                'fcst-hour': 'no-fcst-hour',
                'run-time': ['none', 'none', 'none', 'none'],
                'valid-time': [f'20{header[1][0:2]}', header[1][2:4], header[1][4:6], f'{header[1][7:9]}:{header[1][9:11]}']
            }

            # add clean_data title dict
            clean_data['titles'] = {

                'top_title': f"ACARS AIRCRAFT OBSERVATION VERTICAL PROFILE",
                'left_title': f"VALID: {clean_data['site_info']['valid-time'][1]}-{clean_data['site_info']['valid-time'][2]}-{clean_data['site_info']['valid-time'][0]} {clean_data['site_info']['valid-time'][3]}Z",
                'right_title': f"{clean_data['site_info']['site-id']} - {clean_data['site_info']['site-name']}, {clean_data['site_info']['site-lctn']} | {clean_data['site_info']['site-latlon'][0]}, {clean_data['site_info']['site-latlon'][1]}    "
            }


            print('    > COMPLETE --------')
            elapsed_time = time.time() - st
            print('    > RUNTIME:', time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))

            if not hush:
                print(
                    f"    > SUMMARY: {clean_data['site_info']['valid-time'][3]}Z Flight from {clean_data['site_info']['site-id']}, {clean_data['site_info']['site-name']} at {clean_data['site_info']['valid-time'][1]}-{clean_data['site_info']['valid-time'][2]}-{clean_data['site_info']['valid-time'][0]}-{clean_data['site_info']['valid-time'][3]}Z")
                warnings.filterwarnings("ignore")
                sounding_params(clean_data).print_vals()

            return clean_data

        else:
            return raw_data

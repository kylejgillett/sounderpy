import time
import warnings
from datetime import datetime

import pandas as pd
import numpy as np

from metpy.units import units

from siphon.simplewebservice.wyoming import WyomingUpperAir
from siphon.simplewebservice.igra2 import IGRAUpperAir

from .utils import get_latlon
from .calc import sounding_params







"""
    SOUNDERPY RAOB OBSERVATIONS 'GET-DATA' FUNCTIONS  

    Purpose of module: 

    House function for loading and parsing RAOB & IGRAv2 observation
    data. Functions here are referenced by sounderpy.py


    (C) KYLE J GILLETT, UNIVERSITY OF NORTH DAKOTA, 2024
"""






#######################
# RAOB AND IRGAv2 OBSERVATIONS DATA
#########################################################################

def fetch_obs(station, year, month, day, hour, hush, clean_it):

    # record process time
    st = time.time()

    r"""
       Return a ``dict`` of 'cleaned up' observed profile data
    
       :param station: a three digit RAOB identifier (such as: 'DTX') or 11 digit IGRAv2 identifier (such as: 'GMM00010393')
       :type station: str, required
       :param year: launch year
       :type year: str, required
       :param month: launch month
       :type month: str, required
       :param day: launch day
       :type day: str, required
       :param hour: launch hour
       :type hour: str, required
       :param clean_it: whether to return the raw_data object or a clean_data dict.
       :type clean_it: bool, optional, default is `True`
       :return: clean_data, a dict of ready-to-use vertical profile data including pressure, height, temperature, dewpoint, u-wind, v-wind, & model information
       :rtype: dict
    """

    print(f'> OBSERVED DATA ACCESS FUNCTION\n  -----------------------------------')

    station = str.upper(station)
    # get station lists from SounderPy GitHub Repo
    RAOB_STATIONS = pd.read_csv(f'https://raw.githubusercontent.com/kylejgillett/sounderpy/main/src/RAOB-STATIONS.txt',
                                skiprows=7, skipinitialspace=True)
    IGRA_STATIONS = pd.read_csv(f'https://raw.githubusercontent.com/kylejgillett/sounderpy/main/src/IGRA-STATIONS.txt',
                                skiprows=7, skipinitialspace=True)

    got_data = False

    # set up siphon API call for raob data -- if station ID is found in RAOB_STATIONS, it is
    # a RAOB ID and siphon UW or ISU must be used to get data
    if len(station) == 11:
        search_for = 'igra'
    else:
        search_for = 'raob'






    ### RAOB OBSERVATIONS ###
    #########################################################################################################
    if search_for == 'raob':
        # try this process 10 times, sometimes requests fail due to temporary 404 errors
        for i in range(1, 11):
            try:
                # try UW data request
                df = WyomingUpperAir.request_data(datetime(int(year), int(month), int(day), int(hour)), station)
                got_data = True
                if got_data:
                    print(f'    > PROFILE FOUND: {station} on {month}/{day}/{year} at {hour}z | From UW')
                    break
            except:
                got_data = False
                pass

        # search through RAOB sites list with provided RAOB ID, first try ICAO ID, then WMO ID
        if got_data:
            if clean_it:
                try:
                    station = \
                    RAOB_STATIONS['ICAO'][np.where(RAOB_STATIONS['ICAO'].str.contains(station, na=False, case=True))[0]].values[
                        0].strip()
                    name_idx = 'ICAO'
                except:
                    try:
                        station = RAOB_STATIONS['WMO'][np.where(RAOB_STATIONS['WMO'] == int(station))[0]].values[0]
                        name_idx = 'WMO'
                    except:
                        raise ValueError(
                            f'ICAO or WMO identifier not found, please make sure you provided the correct RAOB ID. If you think this is an error' +
                            'contact the author: https://kylejgillett.github.io/sounderpy/about.html#about-the-author')
                        pass
                    pass

                # begin loading data
                # create dict of data
                new_keys = ['p', 'z', 'T', 'Td', 'u', 'v']
                old_keys = ['pressure', 'height', 'temperature', 'dewpoint', 'u_wind', 'v_wind']  # 'latitude', 'longitude']
                units_list = ['hPa', 'meter', 'degC', 'degC', 'kt', 'kt']
                clean_data = {}
                non_dups = np.concatenate(([True], np.diff(df.to_dict('list')['pressure']) != 0))
                for old_key, new_key, unit in zip(old_keys, new_keys, units_list):
                    clean_data[new_key] = np.array(df.to_dict('list')[old_key])[non_dups] * units(unit)
                clean_data['site_info'] = {
                    'site-id': RAOB_STATIONS[RAOB_STATIONS[name_idx] == station][name_idx].values[0],
                    'site-name': RAOB_STATIONS[RAOB_STATIONS[name_idx] == station]['NAME'].values[0],
                    'site-lctn': RAOB_STATIONS[RAOB_STATIONS[name_idx] == station]['LOC'].values[0],
                    'site-latlon': get_latlon('raob', str(station)),
                    'site-elv': RAOB_STATIONS[RAOB_STATIONS[name_idx] == station]['EL(m)'].values[0],
                    'source': 'RAOB OBSERVED PROFILE',
                    'model': 'no-model',
                    'fcst-hour': 'no-fcst-hour',
                    'run-time': ['none', 'none', 'none', 'none'],
                    'valid-time': [year, month, day, hour]}

                try:
                    # trim data to 98hPa and below for less process time
                    slc = (len(clean_data['p']) - np.where(clean_data['p'] <= 98. * units('hPa'))[0][0])
                    for key in new_keys:
                        clean_data[key] = clean_data[key][:-slc]
                except:
                    pass
        else:
            raise ValueError(
                f'Wyoming Upper Air Archive connection failed -- ensure you have the correct dates and corresponding station identifier\n' +
                f'There is likely no available data for station {station} on {month}/{day}/{year} at {hour}z')
    #########################################################################################################







    ### IGRAv2 OBSERVATIONS ###
    #########################################################################################################
    elif search_for == 'igra':
        for i in range(1, 3):
            try:
                # try siphon IGRA request
                df = IGRAUpperAir.request_data(datetime(int(year), int(month), int(day), int(hour)), station)
                got_data = True
                if got_data:
                    print(f'    > PROFILE FOUND: {station} on {month}/{day}/{year} at {hour}z | From IGRAv2')
                    break
            except:
                got_data = False
                pass

        # if data is found, parse data and create a dict of clean data
        if got_data:
            if clean_it:
                station = \
                IGRA_STATIONS['ID'][np.where(IGRA_STATIONS['ID'].str.contains(station, na=False, case=True))[0]].values[
                    0].strip()

                # create dict of data
                head = df[1]
                df = df[0]
                new_keys = ['p', 'z', 'T', 'Td', 'u', 'v']
                old_keys = ['pressure', 'height', 'temperature', 'dewpoint', 'u_wind', 'v_wind']  # 'latitude', 'longitude']
                units_list = ['hPa', 'meter', 'degC', 'degC', 'kt', 'kt']
                clean_data = {}
                zflag = np.array(df['zflag'])
                pflag = np.array(df['pflag'])
                tflag = np.array(df['tflag'])
                for old_key, new_key, unit in zip(old_keys, new_keys, units_list):
                    clean_data[new_key] = np.array(df.to_dict('list')[old_key])[zflag + pflag + tflag >= 4] * units(unit)
                clean_data['site_info'] = {
                    'site-id': IGRA_STATIONS[IGRA_STATIONS['ID'] == station]['ID'].str.strip().values[0],
                    'site-name': IGRA_STATIONS[IGRA_STATIONS['ID'] == station]['NAME'].str.strip().values[0],
                    'site-lctn': '',
                    'site-latlon': get_latlon('igra', station),
                    'site-elv': IGRA_STATIONS[IGRA_STATIONS['ID'] == station]['EL(m)'].values[0],
                    'source': 'RAOB OBSERVED PROFILE',
                    'model': 'no-model',
                    'fcst-hour': 'no-fcst-hour',
                    'run-time': ['none', 'none', 'none', 'none'],
                    'valid-time': [year, month, day, hour]}
                # correct u & v units
                clean_data['u'] = clean_data['u'] * 1.94384
                clean_data['v'] = clean_data['v'] * 1.94384

        else:
            raise ValueError(
                f'IGRAv2 Dataset connection failed -- ensure you have the correct dates and corresponding station identifier\n' +
                f'There is likely no available data for station {station} on {month}/{day}/{year} at {hour}z')
    #########################################################################################################




    if clean_it:

        ### CHECK AND FIX BAD HEIGHT DATA ###
        ######################################
        def are_heights_bad(arr):
            # check if heights are bad

            # inc flag
            increasing = False

            for i in range(1, len(arr)):
                if arr[i] > arr[i - 1]:
                    increasing = True
                elif increasing:  # are_heights_bad == True
                    return True

            return False


        def fix_bad_heights(heights):
            # fix bad heights

            increasing = False

            last_increasing_index = 0
            # search through height array to find the last
            # value in the sequence to 'increase' correctly
            for i in range(1, len(heights)):
                if heights[i] > heights[i - 1]:  # do so by checking if i is > than i-1
                    # record its index
                    last_increasing_index = i
                else:
                    # the last increasing value is found
                    # return its index
                    break
            last_increasing_value = heights[last_increasing_index]

            # for every 'bad' value afterward, add the last_increasing_index value
            # to it, to correct the height data
            for i in range(last_increasing_index + 1, len(heights)):
                heights[i] += last_increasing_value

                fixed_heights = heights

            return fixed_heights


        bad_heights = are_heights_bad(clean_data['z'])

        if bad_heights:
            clean_data['z'] = fix_bad_heights(clean_data['z'])



        # add plot title to clean_data dict
        clean_data['titles'] = {
            # title specific to RAOB Obs
            'top_title': "RAOB OBSERVED VERTICAL PROFILE",
            'left_title': f"VALID: {clean_data['site_info']['valid-time'][1]}-{clean_data['site_info']['valid-time'][2]}-{clean_data['site_info']['valid-time'][0]} {clean_data['site_info']['valid-time'][3]}Z",
            'right_title': f"{clean_data['site_info']['site-id']} - {clean_data['site_info']['site-name']}, {clean_data['site_info']['site-lctn']} | {clean_data['site_info']['site-latlon'][0]}, {clean_data['site_info']['site-latlon'][1]}    "
        }




        print('    > COMPLETE --------')
        elapsed_time = time.time() - st
        print('    > RUNTIME:', time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))

        if not hush:
            print(
                f"    > SUMMARY: {clean_data['site_info']['valid-time'][3]}Z Launch for {clean_data['site_info']['site-id']}, {clean_data['site_info']['site-name']} at {clean_data['site_info']['valid-time'][1]}"
                f"-{clean_data['site_info']['valid-time'][2]}-{clean_data['site_info']['valid-time'][0]}-{clean_data['site_info']['valid-time'][3]}Z")

            warnings.filterwarnings("ignore")

            sounding_params(clean_data).print_vals()

        return clean_data


    else:
        # if user sets `clean_it=False` return raw data ('dirty') format
        return df

    #########################################################################################################
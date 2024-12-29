from datetime import datetime, timedelta
from urllib.request import urlopen
import warnings
import time
import pandas as pd
import numpy as np
import metpy.calc as mpcalc
from metpy.units import units

from .calc import *




"""
    SOUNDERPY BUFKIT 'GET-DATA' FUNCTIONS  

    Purpose of module: 

    House functions for loading and parsing BUFKIT model forecast data.
    Functions here are referenced by sounderpy.py


    (C) KYLE J GILLETT, UNIVERSITY OF NORTH DAKOTA, 2024
"""



#######################
# BUFKIT MODEL FORECAST DATA
#########################################################################

def fetch_bufkit(model, station, fcst_hour, run_year, run_month, run_day, run_hour,
                   hush, clean_it):

    # record process time
    st = time.time()

    r"""Get BUFKIT forecast model vertical profile data
       Return a ``dict`` of 'cleaned up' model forecast data from a given model, for a given BUFKIT site identifier, forecast hour, & model-run-date

       :param model: the requested model to use (such as hrrr, nam, gfs, etc)
       :type model: str, required
       :param station: a 3-4 digit BUFKIT site identifier
       :type station: str, required
       :param fcst_hour: valid forecast hour
       :type fcst_hour: int, required
       :param year: valid year
       :type year: str, required
       :param month: valid month
       :type month: str, required
       :param day: valid day
       :type day: str, required
       :param hour: valid hour
       :type hour: str, required
       :param hush: whether to 'hush' a read-out of thermodynamic and kinematic parameters when getting a data.
       :type hush: bool, optional, default is `False`
       :return: clean_data, a dict of ready-to-use vertical profile data including pressure, height, temperature, dewpoint, u-wind, v-wind, & model information
       :rtype: dict
    """

    print(f'> BUFKIT DATA ACCESS FUNCTION\n   ---------------------------------')

    # make sure variables are in the correct case
    model = str.lower(model)
    station = str.upper(station)

    # remove '#' for URL
    if '#' in station:
        url_station = station.replace('#', '%23')
    else:
        url_station = station



    # GET MOST-RECENT RUNS FROM PSU SERVERS
    # if date variables (year, month, day) are not given, the user has 'selected' a most
    # recent forecast run, get that from PSU
    if run_year == None:
        if model not in ['gfs', 'nam', 'namnest', 'rap', 'hrrr', 'sref', 'hiresw']:
            raise ValueError(
                f"{model} is not a valid model option. Valid models include ['GFS', 'NAM', 'NAMNEST', 'RAP', 'HRRR', 'SREF', 'HIRESW']")
        if model == 'gfs':
            model3 = 'gfs3'
        else:
            model3 = model
        data_conn = f'https://www.meteo.psu.edu/bufkit/data/{model.upper()}/{model3}_{url_station.lower()}.buf'




    # GET ARCHIVE DATA FROM THE IEM SERVERS. CORRECT GFS & NAM MODEL NAMES
    # if date variables (year, month, day) are given, the user has 'selected' a
    # archived forecast for the given date
    else:
        if model not in ['gfs', 'nam', 'namnest', 'rap', 'hrrr']:
            raise ValueError(
                f"{model} is not a valid model option. Valid models include ['GFS', 'NAM', 'NAMNEST', 'RAP', 'HRRR']")
        if model == 'namnest':
            model = 'nam4km'
        if model == 'gfs':
            model3 = 'gfs3'
        else:
            model3 = model
        data_conn = f'https://mtarchive.geol.iastate.edu/{run_year}/{run_month}/{run_day}/bufkit/{run_hour}/{model}/{model3}_{url_station.lower()}.buf'



    # Check to make sure the user-defined site ID is a valid bufkit site before continuing
    try:
        # GET BUFKIT STATIONS LISTING FROM SOUNDERPY GITHUB REPO
        BUFKIT_STATIONS = pd.read_csv(
            f'https://raw.githubusercontent.com/kylejgillett/sounderpy/main/src/BUFKIT-STATIONS-MASTER.txt',
            skiprows=7, skipinitialspace=True)
        # ATTEMPT TO FIND THE STATION
        station = \
        BUFKIT_STATIONS['ID'][np.where(BUFKIT_STATIONS['ID'].str.contains(station, na=False, case=True))[0]].values[0]
        worked = True
    except:
        worked = False
        pass

    if not worked:
        raise ValueError(
            f"{station} does not appear to be a valid BUFKIT site identifier. A map of valid BUFKIT stations can be found here from Penn State: http://www.meteo.psu.edu/bufkit/CONUS_RAP_00.html")




    # GET BUFKIT FILE
    # CONVERT LINES OF BYTES TO STRINGS
    buf_file = urlopen(data_conn)
    buf_file = [str(line).replace("b'", "").replace("\\r\\n'", "") for line in buf_file]


    if clean_it:
        # CREATE TEMP DATA
        tmp_data, sounding_headers, derived_headers = [], '', ''
        recordSounding = False



        # SET UP DATE / TIME OBJECTS FROM THE BUFKIT FILE
        run_time = buf_file[4][buf_file[4].index('TIME') + 7:(buf_file[4].index('TIME') + 9) + 9]
        run_dt = datetime(int(f'20{run_time[0:2]}'), int(run_time[2:4]), int(run_time[4:6]), int(run_time[7:9]))
        fct_dt = run_dt + timedelta(hours=fcst_hour)
        hr_deltas = {
            'gfs': [1, 180], 'hrrr': [1, 48],
            'rap': [1, 51], 'nam': [1, 48],
            'namnest': [1, 60], 'nam4km': [1, 60],
            'sref': [1, 84], 'hiresw': [1, 48]}
        stp_dt = fct_dt + timedelta(hours=hr_deltas[model][0])



        if hr_deltas[model][1] < fcst_hour:
            raise ValueError(
                f'Invalid forecast hour -- BUFKIT only stores up to F0{hr_deltas[model][1]} for the {str.upper(model)}')



        # Loop over each line in data file
        for line in buf_file:
            # Find start of sounding data
            if f'TIME = {fct_dt.strftime("%Y")[2:4]}{fct_dt.strftime("%m")}{fct_dt.strftime("%d")}/{fct_dt.strftime("%H")}00' in line:
                recordSounding = True
            if 'SNPARM' in line:
                sounding_headers = line[line.index('=') + 2:].replace(' ', '').split(';')
            if 'STNPRM' in line:
                derived_headers = line[line.index('=') + 2:].replace(' ', '').split(';')
            # Append data line to temp data list
            if recordSounding:
                tmp_data.append(line)
            # Break out of loop when end key reached
            if f'TIME = {stp_dt.strftime("%Y")[2:4]}{stp_dt.strftime("%m")}{stp_dt.strftime("%d")}/{stp_dt.strftime("%H")}00' in line:
                tmp_data.pop(-1)
                break
            elif 'YYMMDD/HHMM' in line:
                tmp_data.pop(-1)
                break



        # SET UP UTILS
        station_headers = ['STID', 'STNM', 'TIME', 'SLAT', 'SLON', 'SELV', 'STIM']
        tmp_str = ''
        recordStationInfo, recordDerivedQty, recordSoundingQty = False, False, True
        station_metadata, derived_data, sounding_data = [], [], []



        # Check if the last line is blank
        if tmp_data[-1].strip() != '':
            # Add a blank line at the end
            tmp_data.append('')



        # PARSE THROUGH FILE, SPLIT LINES AND RECORD DATA WE WANT TO KEEP
        for line in tmp_data:
            # Check for station information
            if recordStationInfo and line == '':
                # Break values up to only be seperated by one whitespace
                station_info = (tmp_str.replace(' = ', ' '))
                # Split values into list
                station_info = station_info.split(' ')
                # Remove label values
                station_info = [x for x in station_info if x not in station_headers]
                while '' in station_info:
                    station_info.remove('')
                # Add to main list
                station_metadata.append(station_info)
                # Reset temp vars
                tmp_str = ''
                recordStationInfo = False
            if any(var in line for var in station_headers):
                recordStationInfo = True
                tmp_str += (' ' + line)
            # Check for derived sounding quantities
            if recordDerivedQty == True and line == '':
                # Break values up to only be seperated by one whitespace
                derived_qty = (tmp_str.replace(' = ', ' '))
                # Split values into list
                derived_qty = derived_qty.split(' ')
                # Remove non-numeric values
                derived_qty = [x for x in derived_qty if x not in derived_headers]
                while '' in derived_qty:
                    derived_qty.remove('')
                # Add to main list
                derived_data.append(derived_qty)
                # Reset temp vars
                tmp_str = ''
                recordDerivedQty = False
            if any(var in line for var in derived_headers):
                recordDerivedQty = True
                tmp_str += (' ' + line)
            # Check for sounding quantities
            if any(var in line for var in sounding_headers):
                recordSoundingQty = True
            if recordSoundingQty and line == '':
                level_list = []
                # Split data string into values
                data_list = tmp_str.split(' ')
                # Remove empty indices
                while '' in data_list:
                    data_list.remove('')
                # Break data up into pressure levels
                for i in range(0, len(data_list), len(sounding_headers)):
                    level_list.append(data_list[i:len(sounding_headers) + i])
            elif recordSoundingQty:
                if any(var in line for var in sounding_headers) == False:
                    tmp_str += (' ' + line)



        if 'level_list' not in locals():
            raise ValueError(
                f"The data for the model and forecast hour you requested, from the model-run date you requested, at the BUFKIT site you requested, does not appear to exist\n" +
                f"Please try a different model, forecast hour, run date or BUFKIT site")


        # CREATE BLANK LISTS
        p = []
        z = []
        T = []
        Td = []
        ws = []
        wd = []
        omega = []

        # APPEND LISTS WITH DATA FROM BUFKIT FILES
        if model in ['gfs']:
            for i in range(0, len(level_list)):
                p.append(float(level_list[i][0]))
                z.append(float(level_list[i][8]))
                T.append(float(level_list[i][1]))
                Td.append(float(level_list[i][3]))
                ws.append(float(level_list[i][6]))
                wd.append(float(level_list[i][5]))
                omega.append(float(level_list[i][7]))
        else:
            for i in range(0, len(level_list)):
                p.append(float(level_list[i][0]))
                z.append(float(level_list[i][9]))
                T.append(float(level_list[i][1]))
                Td.append(float(level_list[i][3]))
                ws.append(float(level_list[i][6]))
                wd.append(float(level_list[i][5]))
                omega.append(float(level_list[i][7]))


        # CALCULATE U AND V COMPONENTS
        u = list(mpcalc.wind_components(ws * units.kts, wd * units.degrees)[0].m)
        v = list(mpcalc.wind_components(ws * units.kts, wd * units.degrees)[1].m)



        # DEFINE find_nearest() FUNCTION
        def find_nearest(array, value):
            array = np.asarray(array)
            nearest_idx = (np.abs(array - value)).argmin()
            return nearest_idx



        # FIND P LEVEL AT 50hPa
        hPa50 = find_nearest(p, 50)



        # ARRANGE DATA IN CLEAN_DATA DICT
        clean_data = {}
        lists = [p[0:hPa50], z[0:hPa50], T[0:hPa50], Td[0:hPa50], u[0:hPa50], v[0:hPa50], omega[0:hPa50]]
        keys = ['p', 'z', 'T', 'Td', 'u', 'v', 'omega']
        units_list = ['hPa', 'meter', 'degC', 'degC', 'kt', 'kt', 'Pa/sec']
        for key, lst, unit in zip(keys, lists, units_list):
            clean_data[key] = lst * units(unit)



        # create dict of data
        clean_data['site_info'] = {
            'site-id': BUFKIT_STATIONS[BUFKIT_STATIONS['ID'] == station]['ID'].str.strip().values[0],
            'site-name': BUFKIT_STATIONS[BUFKIT_STATIONS['ID'] == station]['NAME'].str.strip().values[0],
            'site-lctn': BUFKIT_STATIONS[BUFKIT_STATIONS['ID'] == station]['LOC'].str.strip().values[0],
            'site-latlon': [BUFKIT_STATIONS[BUFKIT_STATIONS['ID'] == station]['LAT'].values[0],
                            BUFKIT_STATIONS[BUFKIT_STATIONS['ID'] == station]['LON'].values[0]],
            'site-elv': BUFKIT_STATIONS[BUFKIT_STATIONS['ID'] == station]['EL(m)'].values[0],
            'source': 'BUFKIT FORECAST PROFILE',
            'model': str.upper(model),
            'fcst-hour': f'F0{fcst_hour}',
            'run-time': [run_dt.strftime("%Y"), run_dt.strftime("%m"), run_dt.strftime("%d"), run_dt.strftime("%H")],
            'valid-time': [fct_dt.strftime("%Y"), fct_dt.strftime("%m"), fct_dt.strftime("%d"), fct_dt.strftime("%H")]
        }

        clean_data['titles'] = {
            'top_title': f"BUFKIT MODEL FORECAST PROFILE | {clean_data['site_info']['run-time'][3]}Z {clean_data['site_info']['model']} {clean_data['site_info']['fcst-hour']}",
            'left_title': f" RUN: {clean_data['site_info']['run-time'][1]}/{clean_data['site_info']['run-time'][2]}/{clean_data['site_info']['run-time'][0]} {clean_data['site_info']['run-time'][3]}Z  |  VALID: {clean_data['site_info']['valid-time'][1]}/{clean_data['site_info']['valid-time'][2]}/{clean_data['site_info']['valid-time'][0]} {clean_data['site_info']['valid-time'][3]}Z",
            'right_title': f"{clean_data['site_info']['site-id']} - {clean_data['site_info']['site-name']}, {clean_data['site_info']['site-lctn']} | {clean_data['site_info']['site-latlon'][0]}, {clean_data['site_info']['site-latlon'][1]}    "
        }


        def clean_dewpoints(td):
            for i in range(len(td)):
                if td[i] < -130 * td[i].units: td[i] = -130 * td[i].units

        clean_dewpoints(clean_data['Td'])



        print('    > COMPLETE --------')
        elapsed_time = time.time() - st
        print('    > RUNTIME:', time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))


        if not hush:
            print(
                f"    > SUMMARY: SUMMARY: {clean_data['site_info']['run-time'][3]}Z {clean_data['site_info']['model']} {clean_data['site_info']['fcst-hour']} for {clean_data['site_info']['site-id']},"
                f"{clean_data['site_info']['site-name']} at {clean_data['site_info']['valid-time'][1]}-{clean_data['site_info']['valid-time'][2]}-{clean_data['site_info']['valid-time'][0]}-{clean_data['site_info']['valid-time'][3]}Z")

            warnings.filterwarnings("ignore")

            sounding_params(clean_data).print_vals()


        return clean_data

    else:
        return buf_file

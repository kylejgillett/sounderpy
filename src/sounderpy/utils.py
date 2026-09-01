from datetime import datetime
import csv
import requests
import pandas as pd
import numpy as np
import numpy.ma as ma
import metpy.calc as mpcalc
from metpy.units import units
import copy
import math






"""
    SOUNDERPY UTILITY FUNCTIONS 

    Purpose of module: 

     A collection of helper functions that users may find useful for processing 
     vertical profile data but are not necessary to use the basic functions of 
     SounderPy. May be called by other modules or explicitly called by the user.


    (C) KYLE J GILLETT, UNIVERSITY OF NORTH DAKOTA, 2024
"""





#########################
# MODIFY SURFACE VALUES OF SOUNDING ARRAYS
#########################################################################
def modify_surface(clean_data, modify_sfc):
    """
    where,
    modify_sfc = {'T': 20, 'Td': 12, 'ws': 10, 'wd': 270}
    and units = degC, degC, kts, degrees
    """

    sounding_data = copy.deepcopy(clean_data)

    # Modify arrays based on the provided modifications
    if 'T' in modify_sfc:
        sounding_data['T'][0] = modify_sfc['T']*units.degC

    if 'Td' in modify_sfc:
        sounding_data['Td'][0] = modify_sfc['Td']*units.degC

    if 'ws' in modify_sfc:
        if 'wd' in modify_sfc:
            new_u, new_v = mpcalc.wind_components(modify_sfc['ws'] * units.kts, modify_sfc['wd'] * units.deg)

            sounding_data['u'][0] = new_u
            sounding_data['v'][0] = new_v

    return sounding_data





#########################
# BARNES INTERP VERTICAL ARRAYS FOR SFC MODIFICATION
#########################################################################
def barnes_interp(new_val, xarr, yarr, kappa=0.5, num_points=5):
    """
    conduct barnes interp on a vertical array of T,Td,u, or v arrays

    - new_val: The new value for the surface point (index 0).
    - xarr: The original array of data (temperature, dewpoint, u, or v).
    - yarr: The pressure levels corresponding to the data.
    - kappa: Smoothing parameter for Barnes interpolation.
    - num_points: The number of lowest points to apply the interpolation to (default is 5).

    Currently not being used, considered a beta feature.
    """

    xarr = np.array(xarr)
    pressure_levels = np.array(yarr)

    # number of points to interp through (lowest 5)
    num_points = min(num_points, len(xarr))

    # replace first val in xarr with new sfc val
    xarr[0] = new_val

    # interpolate the xarr with the new val against the yarr
    interpolated_data = np.copy(xarr)
    for i in range(1, num_points):
        influence = np.exp(-((yarr[i] - yarr[0]) ** 2) / (kappa * (yarr[num_points - 1] - yarr[0]) ** 2))
        interpolated_data[i] = influence * new_val + (1 - influence) * xarr[i]

    # return interp'd xarr
    return interpolated_data





#########################
# INTERPOLATE DATA
#########################################################################
def interp_data(variable, heights, step=100):
    '''
    Interpolate a 1D array of data (such as a temperature profile) over a given interval (step) based on a corresponding array of height values.

    :param variable: an array of data to be interpolated. Must be same length as height array.
    :type variable: arr, required
    :param heights: heights corresponding to the vertical profile used to interpolate. Must be same length as variable array.
    :type heights: arr, required
    :param step: the resolution of interpolation. Default is 100 (recommended value is 100)
    :type step: int, optional
    :return: interp_var, an array of interpolated data.
    :rtype: arr
    '''

    try:
        variable.units
        variable = variable.m
    except:
        variable = variable
    try:
        heights.units
        heights = heights.m
    except:
        heights = heights

    levels = np.arange(0, np.max(heights), step)
    varinterp = np.zeros(len(levels))
    for i in range(0, len(levels)):
        lower = np.where(heights - levels[i] <= 0, heights - levels[i], -np.inf).argmax()
        varinterp[i] = (((variable[lower + 1] - variable[lower]) / (heights[lower + 1] - heights[lower])) * (
                    levels[i] - heights[lower]) + variable[lower])
    return varinterp






#########################
# FIND NEAREST
#########################################################################
def find_nearest(array, value):
    """
    search through an array to find the index of the value nearest to a given value
    """
    array = np.asarray(array)
    nearest_idx = (np.abs(array - value)).argmin()
    return nearest_idx








#########################
# TXT MANIPULATION | mag and mag_round functions
#########################################################################
# make text prettier
# magnitude number text
def mag(param):
    if ma.is_masked(param):
        fixed = '---'
    else:
        try:
            fixed = int(param.m)
        except:
            try:
                fixed = int(param)
            except:
                fixed = param
    return fixed


# make text prettier
# magnitude and round number text
# kawrgs - 'param': the value to round; 'dec': the decimal places to round to;
# 'mag': bool, remove units (magintude of the value)
def mag_round(param, dec, mag=False):
    if ma.is_masked(param):
        fixed = '---'
    elif mag == True:
        fixed = np.round(param.m, dec)
    else:
        fixed = np.round(param, dec)
    return fixed




#########################
# GET SFC INDEX
#########################################################################
def get_sfc_index(height_arr):
    """
    Return a value of an index of an array who's value is closest to a define value.

    :param array: an array of data to be searched through
    :type array: arr, required
    :param heights: the value used to compare against the array of data
    :type heights: int or float, required
    :return: nearest_idx, index of the data array that corresponds with the nearest value to the given value
    :rtype: int
    """

    i = 0
    # Search the array for the sfc index
    while i < len(height_arr):
        if height_arr[i] >= 0:
            return i
        else:
            i += 1
    # Did not find a positive index
    return -1







####################
# MAKE SURFACE BASED
#########################################################################
def make_sfc_based(arr, sfc_val, sfc_index):
    """
    takes an array and a valid index in that array, then returns a copy of the
    array beginning at the provided index i.e., chops off below-ground values
    """

    # Initialize an empty numpy array as the modified array
    mod_arr = np.empty(len(arr) - sfc_index)
    # Insert the sfc and higher values into the modified array
    i = 0
    while i < len(mod_arr):
        mod_arr[i] = arr[i + sfc_index]
        i = i + 1
    # Inserts surface values
    mod_arr = np.insert(mod_arr, 0, sfc_val)
    return mod_arr








#########################
# MAKE SURFACE-BASED 3D
#########################################################################
def make_sfc_based_3D(arr, sfc_arr):
    '''
    takes a 3D array of mandatory level data and a 2D array of surface data,
    appends the surface data onto the mandatory level array, and returns a single array
    of both surface and mandatory level data
    '''

    mod_arr = np.zeros((np.shape(arr)[0] + 1, np.shape(arr)[1], np.shape(arr)[2]))
    for j in range(np.shape(arr)[1]):
        for k in range(np.shape(arr)[2]):
            mod_arr[0, j, k] = sfc_arr[j, k]
            mod_arr[1:, j, k] = arr[:, j, k]
    return mod_arr




#########################
# LAT-LON FINDER FUNCTION
#########################################################################

def get_latlon(station_type, station_id):
    '''
    Return a latitude-longitude float pair in a ``list``

    :param station_type: the station 'type' that corresponds with the given station ID
    :type station_type: str, required
    :param station_id: the station ID for the given station type
    :type station_id: str, required
    :return: lat/lon float pair
    :rtype: list
    '''

    station_id = str.upper(station_id)

    # DMS to decimal degrees
    def dms2dd_min(degrees, minutes, direction):
        dd = float(degrees) + float(minutes) / 60
        if direction == "S" or direction == "W":
            dd *= -1
        return dd

    def dms2dd(degrees, direction):
        dd = float(degrees)
        if direction == "S" or direction == "W":
            dd *= -1
        return dd



    ############################################### METAR #####################################################
    if station_type.casefold() == 'metar':
        '''
        takes a METAR site id such as 'KMBS', searches over 9000 station IDs  and returns
        a list including the lat/lon for the METAR site
        '''
        # get METAR stations list
        request = requests.get("https://raw.githubusercontent.com/kylejgillett/sounderpy/main/src/METAR-STATIONS.txt",
                               stream=True)
        stations = {}
        for line in request.iter_lines():
            data = line.decode("ascii")
            if data:
                if data[0] == "!" or len(data) != 83:
                    continue
                province = data[0:2]
                station = data[3:19].strip()
                icao = data[20:24].strip()
                lat = dms2dd_min(data[39:41], data[42:44], data[44:45])
                lon = dms2dd_min(data[47:50], data[51:53], data[53:54])
                altitude = int(data[55:59])
                country = data[81:83]

                if icao:
                    stations[icao] = {"name": station, "lat": lat, "lon": lon, "altitude": altitude, "country": country}
        try:
            latlon = [np.round(stations.get(station_id)['lat'], 2), np.round(stations.get(station_id)['lon'], 2)]
            return latlon
        except:
            raise ValueError(f"The station you requested ({station_id}) doesn't seem to exist\n" +
                             "TIP: most METAR IDs include a 'K' in front, such as 'KMOP'")



    ############################################### BUFKIT #####################################################
    elif station_type.casefold() == 'bufkit':
        '''
        takes a BUFKIT site id such as 'KMOP', searches over 1200 station IDs and returns
        a list including the lat/lon for the BUFKIT site
        '''
        # get BUFKIT stations list
        BUFKIT_STATIONS = pd.read_csv(
            f'https://raw.githubusercontent.com/kylejgillett/sounderpy/main/src/BUFKIT-STATIONS-MASTER.txt',
            skiprows=7, skipinitialspace=True)

        try:
            station = BUFKIT_STATIONS['ID'][
                np.where(BUFKIT_STATIONS['ID'].str.contains(station_id, na=False, case=True))[0]].values[0]
            lat = (BUFKIT_STATIONS[BUFKIT_STATIONS['ID'] == station]['LAT'].values[0])
            lon = (BUFKIT_STATIONS[BUFKIT_STATIONS['ID'] == station]['LON'].values[0])
            return [lat, lon]
        except:
            raise ValueError(f"The station you requested ({station_id}) doesn't seem to exist\n" +
                             "TIP: some IDs include a 'K' in front, such as 'KMOP', others are 3 digits, such a 'DTX'")



    ############################################### RAOB #####################################################
    elif station_type.casefold() == 'raob':

        '''
        takes a RAOB site id such as 'DTX', searches over 9000 station IDs  and returns
        a list including the lat/lon for the RAOB site
        '''
        # get RAOB stations list
        RAOB_STATIONS = pd.read_csv(
            f'https://raw.githubusercontent.com/kylejgillett/sounderpy/main/src/RAOB-STATIONS.txt',
            skiprows=7, skipinitialspace=True)

        # find lat-lon from stations list if it exists
        try:
            station = RAOB_STATIONS['ICAO'][
                np.where(RAOB_STATIONS['ICAO'].str.contains(station_id, na=False, case=True))[0]].values[0]
            lat = dms2dd(RAOB_STATIONS[RAOB_STATIONS['ICAO'] == station]['LAT'].values[0],
                         RAOB_STATIONS[RAOB_STATIONS['ICAO'] == station]['A'].values[0])
            lon = dms2dd(RAOB_STATIONS[RAOB_STATIONS['ICAO'] == station]['LON'].values[0],
                         RAOB_STATIONS[RAOB_STATIONS['ICAO'] == station]['B'].values[0])
            return [lat, lon]
        except:
            try:
                station = RAOB_STATIONS['WMO'][RAOB_STATIONS[RAOB_STATIONS['WMO'] == int(station_id)].index[0]]
                lat = dms2dd(RAOB_STATIONS[RAOB_STATIONS['WMO'] == station]['LAT'].values[0],
                             RAOB_STATIONS[RAOB_STATIONS['WMO'] == station]['A'].values[0])
                lon = dms2dd(RAOB_STATIONS[RAOB_STATIONS['WMO'] == station]['LON'].values[0],
                             RAOB_STATIONS[RAOB_STATIONS['WMO'] == station]['B'].values[0])
                return [lat, lon]
            except:
                raise ValueError(f"The station you requested ({station_id}) doesn't seem to exist")



    ############################################### IGRA #####################################################
    elif station_type.casefold() == 'igra':
        '''
        takes a IGRA2 site id such as 'GMM00010393', searches nearly 3000 station IDs  and returns
        a list including the lat/lon for the IGRA2 site
        '''

        # get IGRA stations list
        IGRA_STATIONS = pd.read_csv(
            f'https://raw.githubusercontent.com/kylejgillett/sounderpy/main/src/IGRA-STATIONS.txt',
            skiprows=7, skipinitialspace=True)

        # find lat-lon from stations list if it exists
        try:
            station = \
            IGRA_STATIONS['ID'][np.where(IGRA_STATIONS['ID'].str.contains(station_id, na=False, case=True))[0]].values[
                0]

            lat = np.round(IGRA_STATIONS[IGRA_STATIONS['ID'] == station]['LAT'].values[0], 2)
            lon = np.round(IGRA_STATIONS[IGRA_STATIONS['ID'] == station]['LON'].values[0], 2)
            return [lat, lon]
        except:
            raise ValueError(f"The station you requested ({station_id}) doesn't seem to exist")



    ############################################### BUOY #####################################################
    elif station_type.casefold() == 'buoy':
        '''
        takes a BUOY/CMAN site id such as '41001', searches through a number of station IDs and returns
        a list including the lat/lon for the BUOY/CMAN  site
        '''

        # get buoy stations list
        BUOY_STATIONS = pd.read_csv(
            'https://raw.githubusercontent.com/kylejgillett/sounderpy/main/src/BUOY-STATIONS.txt',
            skiprows=7, skipinitialspace=True)

        # find lat-lon from stations list if it exists
        try:
            station = \
            BUOY_STATIONS['ID'][np.where(BUOY_STATIONS['ID'].str.contains(station_id, na=False, case=True))[0]].values[
                0]

            lat = dms2dd(BUOY_STATIONS[BUOY_STATIONS['ID'] == station]['LAT'].values[0],
                         BUOY_STATIONS[BUOY_STATIONS['ID'] == station]['A'].values[0])
            lon = dms2dd(BUOY_STATIONS[BUOY_STATIONS['ID'] == station]['LON'].values[0],
                         BUOY_STATIONS[BUOY_STATIONS['ID'] == station]['B'].values[0])
            return [lat, lon]
        except:
            raise ValueError(f"The station you requested ({station_id}) doesn't seem to exist\n" +
                             "TIP: buoy IDs typically look like this: '41001'")
    else:
        raise ValueError(
            f"Incorrect station_type argument. Valid station_type-s are 'metar', 'raob', 'igra', 'bufkit', 'buoy'")



#########################
# U/V WIND TO DIR STRING
#########################################################################
def wind_to_dir(u, v):

    wind_dir_deg = (math.degrees(math.atan2(u, v)) + 360) % 360

    dirs = ["N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE",
            "S", "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW"]

    # Each sector covers 360/16 = 22.5 degrees
    idx = int((wind_dir_deg + 11.25) // 22.5) % 16

    return dirs[idx]
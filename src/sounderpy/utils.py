from datetime import datetime
import csv
import requests
import pandas as pd
import numpy as np
import numpy.ma as ma
import metpy.calc as mpcalc
from metpy.units import units
import copy
import difflib
import functools
import importlib.resources
import math






"""
    SOUNDERPY UTILITY FUNCTIONS 

    Purpose of module: 

     A collection of helper functions that users may find useful for processing 
     vertical profile data but are not necessary to use the basic functions of 
     SounderPy. May be called by other modules or explicitly called by the user.


    (C) KYLE J GILLETT, UNIVERSITY OF NORTH DAKOTA, 2024
"""


class SounderPyWarning(UserWarning):
    """A warning about the user's figure or data, such as missing radar; the CLI shows these even without --verbose."""


class SounderPyNotice(UserWarning):
    """Information rather than a problem, such as the site meta-data a file reader couldn't fill in.

    Filter these on their own (``warnings.filterwarnings("ignore", category=SounderPyNotice)``) to keep a
    log quiet while still seeing ``SounderPyWarning``s about the data. The CLI hides them.
    """





#########################
# MODIFY SURFACE VALUES OF SOUNDING ARRAYS
#########################################################################
def modify_surface(clean_data, modify_sfc):
    """
    where,
    modify_sfc = {'T': 20, 'Td': 12, 'ws': 10, 'wd': 270, 'method': 'mix'}
    and units = degC, degC, kts, degrees

    method 'mix' (the default):
    A warmer surface temperature is mixed upward along its dry adiabat until it meets the
    original profile (Holzworth 1964), the depth a heated surface mixes to. Through that
    layer, a new surface dewpoint is held at a constant mixing ratio, as in a well-mixed
    boundary layer (Betts 1973; Carson 1973; Tennekes 1973). A cooler surface temperature,
    a dewpoint-only change, and the wind change only the surface level, since a stable
    surface layer and surface-layer wind shear are realistic. Dewpoint is capped at the
    temperature.

    method 'replace':
    Only the surface level changes, whatever the new values, so a superadiabatic surface
    layer is kept (useful when the lowest few hundred meters matter, e.g. blowing snow).

    method 'blend' (with 'depth', m AGL, default 500):
    The surface temperature and dewpoint changes taper linearly to zero at 'depth', so a
    warmer surface gives a superadiabatic layer of that depth and the profile above it is
    unchanged. Dewpoint is capped at the temperature.

    Holzworth, G. C., 1964: Estimates of mean maximum mixing depths in the contiguous
        United States. Mon. Wea. Rev., 92, 235-242.
    Betts, A. K., 1973: Non-precipitating cumulus convection and its parameterization.
        Quart. J. Roy. Meteor. Soc., 99, 178-196.
    Carson, D. J., 1973: The development of a dry inversion-capped convectively unstable
        boundary layer. Quart. J. Roy. Meteor. Soc., 99, 450-467.
    Tennekes, H., 1973: A model for the dynamics of the inversion above a convective
        boundary layer. J. Atmos. Sci., 30, 558-567.
    """

    method = modify_sfc.get('method', 'mix')
    if method not in ('mix', 'replace', 'blend'):
        raise ValueError(f"modify_sfc 'method' must be 'mix', 'replace' or 'blend'. Got: {method!r}")
    depth = modify_sfc.get('depth', 500)
    if method == 'blend' and not depth > 0:
        raise ValueError(f"modify_sfc 'depth' must be a height in m above 0. Got: {depth!r}")

    sounding_data = copy.deepcopy(clean_data)
    p = sounding_data['p'].to(units.hPa)
    T = sounding_data['T'].to(units.degC)
    Td = sounding_data['Td'].to(units.degC)

    T_sfc = modify_sfc.get('T', T[0].m) * units.degC
    Td_sfc = min(modify_sfc.get('Td', Td[0].m), T_sfc.m) * units.degC

    # levels the new surface temperature's dry adiabat is warmer than (the mixed layer)
    top = 0
    if 'T' in modify_sfc and method == 'mix':
        theta = mpcalc.potential_temperature(p, T)
        theta_sfc = mpcalc.potential_temperature(p[0], T_sfc)
        while top + 1 < len(p) and theta[top + 1] < theta_sfc:
            top += 1

    # blend: the surface changes taper linearly to zero at `depth` m AGL
    if method == 'blend':
        z = sounding_data['z'].to(units.m).m
        weight = np.clip(1 - (z - z[0]) / depth, 0, None)
        T = T + (T_sfc - T[0]) * weight
        Td = np.minimum((Td + (Td_sfc - Td[0]) * weight).m, T.m) * units.degC

    T[0] = T_sfc
    Td[0] = Td_sfc
    if top > 0:
        mixed = slice(0, top + 1)
        T[mixed] = mpcalc.temperature_from_potential_temperature(p[mixed], theta_sfc * np.ones(top + 1)).to(units.degC)
        if 'Td' in modify_sfc:
            w_sfc = mpcalc.saturation_mixing_ratio(p[0], Td_sfc)
            Td[mixed] = mpcalc.dewpoint(mpcalc.vapor_pressure(p[mixed], w_sfc * np.ones(top + 1)))
        Td[mixed] = np.minimum(Td[mixed].m, T[mixed].m) * units.degC

    sounding_data['T'] = T.to(clean_data['T'].units)
    sounding_data['Td'] = Td.to(clean_data['Td'].units)

    if 'ws' in modify_sfc:
        if 'wd' in modify_sfc:
            new_u, new_v = mpcalc.wind_components(modify_sfc['ws'] * units.kts, modify_sfc['wd'] * units.deg)

            sounding_data['u'][0] = new_u
            sounding_data['v'][0] = new_v

    return sounding_data





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
# STATION ID LOOKUP
#########################################################################
def find_station(ids, station_id, k_prefix=False, required=True, hint=''):
    '''
    Return the row label of a station ID in a column of station IDs, matching the whole ID
    (ignoring case & padding) instead of any ID that contains it

    :param ids: the column of station IDs to search, such as ``BUFKIT_STATIONS['ID']``
    :type ids: pandas.Series, required
    :param station_id: the station ID to find
    :type station_id: str, required
    :param k_prefix: if a 3-letter ID isn't found, also try it with a 'K' in front ('OAX' -> 'KOAX')
    :type k_prefix: bool, optional, default is `False`
    :param required: raise an error if the station isn't found; otherwise return `None`
    :type required: bool, optional, default is `True`
    :param hint: extra text added to the error message
    :type hint: str, optional
    :return: the row label of the first matching station
    :rtype: int or None
    '''
    # the buoy list marks most IDs with a trailing '*', which isn't part of the ID
    clean_ids = ids.astype(str).str.strip().str.rstrip('*').str.upper()
    target = str(station_id).strip().upper()

    candidates = [target]
    if k_prefix and len(target) == 3 and target.isalpha():
        candidates.append('K' + target)

    for candidate in candidates:
        matches = clean_ids.index[clean_ids == candidate]
        if len(matches):
            return matches[0]

    if not required:
        return None

    close = difflib.get_close_matches(target, clean_ids[clean_ids != 'NAN'].unique(), n=3, cutoff=0.6)
    suggestion = f" Did you mean {', '.join(repr(c) for c in close)}?" if close else ''
    raise ValueError(f"The station you requested ({station_id}) doesn't seem to exist.{suggestion}" + (f"\n{hint}" if hint else ''))





#########################
# RAOB / IGRA STATION LIST
#########################################################################
@functools.lru_cache(maxsize=1)
def load_raob_stations():
    '''
    Return SounderPy's radiosonde station list (WMO, ICAO & IGRAv2 IDs) as a ``pandas.DataFrame``

    The list ships with SounderPy (``sounderpy/RAOB-STATIONS-MASTER.csv``) and is read once per session.
    It is built by ``tools/build_station_list.py`` (which keeps its own copy in ``tools/station_list/``)
    and copied into the package with ``build ... --install``. Every column is text, so WMO
    numbers keep their leading zeros (``'01001'``). The file's columns are lined up with spaces after
    each comma, which are skipped when reading.

    :return: one row per station
    :rtype: pandas.DataFrame
    '''
    path = importlib.resources.files('sounderpy').joinpath('RAOB-STATIONS-MASTER.csv')
    with path.open('r', encoding='utf-8') as f:
        df = pd.read_csv(f, dtype=str, keep_default_na=False, skipinitialspace=True)
    df.columns = df.columns.str.strip()
    return df


def find_raob_station(station_id):
    '''
    Return the station list row for a radiosonde station ID, as a ``dict``

    Accepts an ICAO code (``'KOUN'``, or ``'OUN'`` for U.S. sites), a WMO number (``'72357'``,
    ``'1001'`` for ``'01001'``), an IGRAv2 ID (``'USM00072357'``) or an alternate ICAO code (``'OTX'``).

    :param station_id: the station ID
    :type station_id: str, required
    :return: the station's row (all values are text; blank when unknown)
    :rtype: dict
    '''
    df = load_raob_stations()
    target = str(station_id).strip().upper()

    def first(mask):
        hits = df.index[mask]
        if not len(hits):
            return None
        active = hits[df.loc[hits, 'ACTIVE'] == 'Y']
        return df.loc[(active if len(active) else hits)[0]].to_dict()

    k_target = 'K' + target if len(target) == 3 and target.isalpha() else None
    checks = [df['ICAO'] == target]
    if k_target:
        checks.append(df['ICAO'] == k_target)
    if target.isdigit() and len(target) <= 5:
        checks.append(df['WMO'] == target.zfill(5))
    checks.append(df['IGRA_ID'] == target)
    alt = df['ICAO_ALT'].str.split()
    checks.append(alt.apply(lambda codes: target in codes))
    if k_target:
        checks.append(alt.apply(lambda codes: k_target in codes))

    # checks run in order of preference, but an active station wins over an inactive one: a code that is
    # the main ICAO of a closed station and an alternate of a working one finds the working one
    hits = [row for row in map(first, checks) if row is not None]
    if hits:
        return next((row for row in hits if row['ACTIVE'] == 'Y'), hits[0])

    ids = pd.concat([df['ICAO'], df['WMO'], df['IGRA_ID']])
    close = difflib.get_close_matches(target, ids[ids != ''].unique(), n=3, cutoff=0.6)
    suggestion = f" Did you mean {', '.join(repr(c) for c in close)}?" if close else ''
    raise ValueError(f"The station you requested ({station_id}) isn't in SounderPy's RAOB/IGRA station list.{suggestion}\n"
                     "TIP: use an ICAO code ('KOUN' or 'OUN'), a WMO number ('72357') or an IGRAv2 ID ('USM00072357')")


def raob_station_has_no_known_data(row):
    '''
    ``True`` when a station list row (from ``find_raob_station``) has no IGRAv2 record, wasn't seen in the
    recent Wyoming files and isn't marked active: mostly closed military or old city sites from the legacy
    RAOB list. Wyoming's archive may still hold old soundings for them, so they stay in the list.
    '''
    return not (row.get('IGRA_ID') or row.get('WYOMING_SEEN') or row.get('ACTIVE') == 'Y')


def split_display_name(display_name):
    '''
    Split a station display name into (name, region): ``'Norman, OK'`` -> ``('Norman', 'OK')``,
    ``'Guam'`` -> ``('Guam', '')``
    '''
    name, sep, region = str(display_name).rpartition(', ')
    return (name, region) if sep else (region, '')





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

        idx = find_station(BUFKIT_STATIONS['ID'], station_id, k_prefix=True,
                           hint="TIP: some IDs include a 'K' in front, such as 'KMOP', others are 3 digits, such a 'DTX'")
        return [BUFKIT_STATIONS['LAT'][idx], BUFKIT_STATIONS['LON'][idx]]



    ############################################### RAOB / IGRA #####################################################
    elif station_type.casefold() in ('raob', 'igra'):
        '''
        takes a RAOB or IGRA2 site id such as 'DTX', '72632' or 'USM00072632' and returns
        a list including the lat/lon for the site, from SounderPy's packaged station list
        '''
        row = find_raob_station(station_id)
        if not row['LAT'] or not row['LON']:
            raise ValueError(f"The station you requested ({station_id}) has no known location")
        lat, lon = float(row['LAT']), float(row['LON'])
        if station_type.casefold() == 'igra':
            return [np.round(lat, 2), np.round(lon, 2)]
        return [lat, lon]



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
        idx = find_station(BUOY_STATIONS['ID'], station_id, hint="TIP: buoy IDs typically look like this: '41001'")
        return [dms2dd(BUOY_STATIONS['LAT'][idx], BUOY_STATIONS['A'][idx]),
                dms2dd(BUOY_STATIONS['LON'][idx], BUOY_STATIONS['B'][idx])]
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
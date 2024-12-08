import numpy as np
from datetime import datetime
from scipy import interpolate
from metpy.units import units

"""
    SOUNDERPY WRF UTILITIES

    Purpose of module:

        Parse WRF-ARW simulation output and build a SounderPy
        `clean_data` dictionary for a given latitude and longitude
        point to plot a sounding.

    Input type: netcdf4 dataset obj of a wrfout___.nc file

    (C) KYLE J GILLETT, UNIVERSITY OF NORTH DAKOTA, 2024
"""



#############################################
# FIND NEAREST VALUE IN ARRAY
#############################################

def find_nearest(array, value):
    """
    search through an array to find the index of the value nearest to a given value
    """
    array = np.asarray(array)
    nearest_idx = (np.abs(array - value)).argmin()
    return nearest_idx




#############################################
# PARSE DATASET TIME OBJECT
#############################################
def parse_time(ds):
    # convert byte data to string
    byte_data = b''.join(ds['Times'][:].data[0])
    time_str = byte_data.decode('utf-8')

    # replace underscores with spaces to match datetime format
    time_str = time_str.replace('_', ' ')

    # convert to datetime object
    time_obj = datetime.strptime(time_str, '%Y-%m-%d %H:%M:%S')

    return [str(time_obj)[0:4], str(time_obj)[5:7], str(time_obj)[8:10], str(time_obj)[11:16]]




#############################################
# COMPUTE TEMPERATURE FROM WRF OUTPUT
#############################################
def calc_temp(theta_pert, base_pres):
    '''
        Using theta perturbation from reference theta and pressure output from WRF,
        use a rearanged possion's theta eqn to find isobaric temperature
        where,
         > theta = T(P0/P)*Rd/Cp

        solve for T,
        > T = theta/(P0/P)*Rd/cp
    '''

    # CONSTANTS -------------------
    # specific heat of dry air
    cp = 1006  # J/kg/K

    # gas constant for dry air
    Rd = 287  # J/kg/K

    # reference values
    ref_theta = 300  # K
    ref_press = 100000  # Pa
    # ------------------------------

    # BUILD EQUATION --------------
    # find theta values using theta
    # perturbation from reference theta (300K)
    numerator = (theta_pert + ref_theta)  # K

    denominator = (ref_press / (base_pres)) ** (Rd / cp)  # unitless
    # -------------------------------

    return (numerator / denominator)






#############################################
# COMPUTE DEWPOINT FROM WRF OUTPUT
#############################################
def calc_dwpt(q, p):
    '''
        Using specific humidity and pressure output from WRF,

        1. Compute mixing ratio q/(1-q)
        2. Compute vapor pressure (w/(epsilon + w))*p
        3. Compute dewpoint using Td form of Clausius-Clapeyron

        Td form of Clausius-Clapetron derived from...
        Bolton, D., 1980: The computation of equivalent potential temperature.
            Mon. Wea. Rev., 108, 1046-1053, doi:10.1175/1520-0493(1980)108%3C1046:TCOEPT%3E2.0.CO;2.

    '''

    # CONSTANTS -------------------
    epsilon = 0.622
    # ------------------------------

    # compute mixing ratio
    w = (q) / (1 - q)

    # compute vapor pressure
    e = (w / (epsilon + w)) * p / 100

    # BUILD CLAUSIUS-CLAPEYRON EQUATION --------------
    numerator = (243.5 * np.log(e / 6.112))

    denominator = (17.67 - np.log(e / 6.112))
    # -------------------------------------------------

    return (numerator / denominator)






#############################################
# CREATE SMALL "AVG" BOX TO CREATE PROFILE
#############################################
def profile_location(lat_pt, lon_pt, lat_array, lon_array):
    # find latitude max and min from requested center pt
    max_lat_idx = find_nearest(lat_array, lat_pt + 0.1)
    min_lat_idx = find_nearest(lat_array, lat_pt - 0.1)

    # find longitude mx and min from requested center pt
    min_lon_idx = find_nearest(lon_array, lon_pt - 0.1)
    max_lon_idx = find_nearest(lon_array, lon_pt + 0.1)

    # return a list of min, max lats & min, max lons
    return [min_lat_idx, max_lat_idx, min_lon_idx, max_lon_idx]







#############################################
# BUILD AND RETURN 'clean_data' DICTIONARY
#############################################
def make_wrf_profile(ds, latlon, model_name='WRF-ARW', run_name='CONTROL RUN'):
    '''
        params
            > ds: netcdf file of WRF output
            > latlon: list of latitude and longitude `[45.56, -100.89]`
            > run_name: str, name of the run for title plotting 'CONTROL'
    '''

    # define latitude and longitude max/min index values for "bounding box"
    latmn, latmx, lonmn, lonmx = profile_location(latlon[0], latlon[1], ds['XLAT'][0, :, 0], ds['XLONG'][0, 0, :])

    # build vertical data dictionary
    vert_data = {
        'vert_T': np.mean(calc_temp(ds['T'][0, :, latmn:latmx, lonmn:lonmx], ds['PB'][0, :, latmn:latmx, lonmn:lonmx]),
                          axis=(1, 2)) - 273.15,
        'vert_p': np.mean((ds['PB'][:] + ds['P'][:])[0, :, latmn:latmx, lonmn:lonmx], axis=(1, 2)) / 100,
        'vert_z': np.mean((ds['PHB'][:] + ds['PH'][:])[0, 1:, latmn:latmx, lonmn:lonmx], axis=(1, 2)) / 9.81,
        'vert_u': np.mean(ds['U'][0, :, latmn:latmx, lonmn:lonmx], axis=(1, 2)) * 1.94384,
        'vert_v': np.mean(ds['V'][0, :, latmn:latmx, lonmn:lonmx], axis=(1, 2)) * 1.94384,
        'vert_Td': np.mean(calc_dwpt(ds['QVAPOR'][0, :, latmn:latmx, lonmn:lonmx],
                                     (ds['PB'][:] + ds['P'][:])[0, :, latmn:latmx, lonmn:lonmx]), axis=(1, 2))
    }

    # build sirface data dictionary
    sfc_data = {
        'sfc_T': np.mean(ds['T2'][0, latmn:latmx, lonmn:lonmx]) - 273.15,
        'sfc_p': np.mean(ds['PSFC'][0, latmn:latmx, lonmn:lonmx]) / 100,
        'sfc_z': np.mean(ds['HGT'][0, latmn:latmx, lonmn:lonmx]),
        'sfc_u': np.mean(ds['U10'][0, latmn:latmx, lonmn:lonmx]) * 1.94384,
        'sfc_v': np.mean(ds['V10'][0, latmn:latmx, lonmn:lonmx]) * 1.94384,
        'sfc_Td': np.mean(calc_dwpt(ds['Q2'][0, latmn:latmx, lonmn:lonmx], ds['PSFC'][0, latmn:latmx, lonmn:lonmx]))

    }

    # build latlon metadata dictionary
    latlon_data = {
        'data_lat': (ds['XLAT'][0, latmn:latmx, 0]),
        'data_lon': (ds['XLONG'][0, 0, lonmn:lonmx]),
        'data_latnum': (ds['XLAT'][0, latmn:latmx, 0]).shape[0],
        'data_lonnum': (ds['XLONG'][0, 0, lonmn:lonmx]).shape[0],
        'data_time': parse_time(ds)
    }

    # perform surface base-ing of isobaric values with sfc values
    sb_dict = {}
    new_keys = ['T', 'Td', 'u', 'v', 'z', 'p', ]
    sfc_keys = ['sfc_T', 'sfc_Td', 'sfc_u', 'sfc_v', 'sfc_z', 'sfc_p']
    vert_keys = ['vert_T', 'vert_Td', 'vert_u', 'vert_v', 'vert_z', 'vert_p']

    # create a dict of surface-based data
    for vert_key, sfc_key, new_key in zip(vert_keys, sfc_keys, new_keys):
        sb_dict[new_key] = np.insert(vert_data[vert_key][vert_data['vert_z'] >= sfc_data['sfc_z']], 0,
                                     sfc_data[sfc_key])
    sb_dict['z'] = sb_dict['z'] - sb_dict['z'][0]

    # interpolates data
    dz = 250
    soundingtop_hght = sb_dict['z'][-1]
    toplvl = int(soundingtop_hght / dz) * dz
    numlvls = int(toplvl / dz)
    interp_lvls = np.linspace(0, toplvl, numlvls + 1)

    # prepare new dicts
    keys = ['T', 'Td', 'u', 'v', 'z', 'p', ]
    units_list = ['degC', 'degC', 'kt', 'kt', 'm', 'hPa']
    interp_dict = {}
    zeros_dict = {}
    clean_data = {}

    surface_height = sfc_data['sfc_z']

    # create dict of clean data
    for key in keys:
        interp_dict[key] = (interpolate.interp1d(sb_dict['z'], sb_dict[key]))
        zeros_dict[key] = np.zeros((len(interp_lvls)))
        for zeros_arr in zeros_dict.values():
            zeros_dict[key][0] = sb_dict[key][0]
        for i in range(1, len(zeros_dict[key]), 1):
            zeros_dict[key][i] = interp_dict[key](dz * i)

    for i, unit, key in zip(range(0, len(units_list)), units_list, keys):
        clean_data[key] = zeros_dict[key] * units(unit)

    # add site & profile info to clean_data dict
    clean_data['site_info'] = {
        'site-id': 'no-site-id',
        'site-name': 'no-site-name',
        'site-lctn': 'no-site-location',
        'site-latlon': [latlon[0], latlon[1]],
        'site-elv': surface_height,
        'source': f'MODEL REANALYSIS',
        'model': f'{str.upper(model_name)}',
        'fcst-hour': f'{str.upper(run_name)}',
        'run-time': latlon_data['data_time'],
        'valid-time': latlon_data['data_time'],
        'box_area': f'0.1° BOX AVG'}

    # add plot title to clean_data dict
    clean_data['titles'] = {
        'top_title': f"{str.upper(model_name)} MODEL SIMULATION VERTICAL PROFILE | {str.upper(run_name)}",
        'left_title': f"VALID: {clean_data['site_info']['valid-time'][1]}/{clean_data['site_info']['valid-time'][2]}/{clean_data['site_info']['valid-time'][0]} {clean_data['site_info']['valid-time'][3]}Z",
        'right_title': f"{clean_data['site_info']['site-latlon'][0]}, {clean_data['site_info']['site-latlon'][1]} | {clean_data['site_info']['box_area']}    "
    }

    return clean_data
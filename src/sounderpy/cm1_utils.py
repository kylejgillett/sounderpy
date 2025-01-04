import numpy as np
import metpy.calc as mpcalc
from metpy.units import units

"""
    SOUNDERPY CM1 UTILITIES  

    Purpose of module: 

     A collection of functions that parse and process CM1 input files
     into SounderPy `clean_data` dictionaries for easy analysis of 
     input_sounding files.


    (C) KYLE J GILLETT, UNIVERSITY OF NORTH DAKOTA, 2024, 2025
"""



#########################
# EXCEPTION CHECKS
#########################################################################
def check_sfc_hgt(full_z):
    if full_z[1] <= full_z[0]:
        raise ValueError(
            f"Invalid Surface Height Value | Check or add an elevation value in `meta_data_dict`. Value should be in meters (m).\nYour current elevation is {full_z[0]}m which is larger than the next height value of {full_z[1]}m")


def check_latlon(meta_data_dict):
    if len(meta_data_dict['latlon']) < 2:
        raise ValueError(
            f"Invalid latitude-longitude pair | Check or add a lat-lon pair to `meta_data_dict`. Values should be in a list, ex: [48.57, -100.98]")




#########################
# MAKE CM1 PROFILE
#########################################################################
def make_cm1_profile(filename, meta_data_dict):

    '''
        process a cm1 input_sounding file for sounderpy
    '''


    #############################################
    # CONSTANTS
    #############################################
    Rd = 287.05  # specific gas constant, dry air (J/kg/K)
    cp = 1005.7  # specific heat at constant pressure (J/kg/K)
    g = 9.80665  # gravitational accel. (m/s^2)
    P0 = 100000  # ref pres (hPa)

    #############################################
    # READ CM1 INPUT_SOUNDING VALUES
    #############################################
    with open(filename, 'r') as f:
        # extract surface values from the first line
        sfc_p, sfc_th, sfc_qv = map(float, f.readline().split())

        # extract vertical arrays from the remaining lines
        vert_z, vert_th, vert_qv, vert_u, vert_v = np.loadtxt(f, unpack=True)

    #############################################
    # SORT AND PREPARE DATA
    #############################################
    # determine length of vertical data arrays
    num_vals = len(vert_z) + 1

    # initialize empty arrays of len = num_vals +1 (for sfc val)
    full_z = np.zeros(num_vals, np.float32)
    full_th = np.zeros(num_vals, np.float32)
    full_qv = np.zeros(num_vals, np.float32)
    full_u = np.zeros(num_vals, np.float32)
    full_v = np.zeros(num_vals, np.float32)
    full_T = np.zeros(num_vals, np.float32)
    full_vp = np.zeros(num_vals, np.float32)
    full_Td = np.zeros(num_vals, np.float32)
    full_p = np.zeros(num_vals, np.float32)

    # add vertical array values at and beyond idx 1 in full arrays
    full_z[1:num_vals] = vert_z
    full_th[1:num_vals] = vert_th
    full_qv[1:num_vals] = vert_qv / 1000.
    full_u[1:num_vals] = vert_u
    full_v[1:num_vals] = vert_v

    # add surface values to index 0 in full arrays
    # if elev is not provided by the user, assume hydrostatic balance and
    # compute surface elevation using the hydrostatic equation
    if meta_data_dict['elev'] is None:
        elev = 1.70 * full_z[1] - full_z[2] + 0.25 * full_z[
            3]  # cm1 input has no surface hgt, use quadratic extrapolation to estimate sfc z
    else:
        elev = meta_data_dict['elev']

    full_z[0] = elev

    full_th[0] = float(sfc_th)
    full_qv[0] = float(sfc_qv) / 1000.  # convert to g/kg
    full_u[0] = 1.75 * full_u[1] - full_u[2] + 0.25 * full_u[
        3]  # cm1 input has no surface wind, use quadratic extrapolation to estimate sfc u & v
    full_v[0] = 1.75 * full_v[1] - full_v[2] + 0.25 * full_v[
        3]  # cm1 input has no surface wind, use quadratic extrapolation to estimate sfc u & v
    full_p[0] = float(sfc_p) * 100.  # convert to Pa for later calcs

    #############################################
    # COMPUTE NECESSARY DATA
    #############################################

    # convert hgt values to include elevation
    full_z[1:] = full_z[1:] + elev
    
    # compute remaining pressure values above the surface assuming hydrostatic balance
    # where dp = -rho*g*dz
    for k in np.arange(1, num_vals):
        full_p[k] = full_p[k - 1] * np.exp((g * (full_z[k - 1] - full_z[k])) / (
                    Rd * (((full_th[k] + full_th[k - 1]) / 2.) * (P0 / full_p[k - 1]) ** -(Rd / cp))))

    for k in np.arange(num_vals):
        # compute temperature from theta
        full_T[k] = mpcalc.temperature_from_potential_temperature(full_p[k] * units.Pa, full_th[k] * units.K).m
        # compute vapor pressure using the mixing ratio and pressure
        full_vp[k] = mpcalc.vapor_pressure(full_p[k] * units.Pa, full_qv[k]).m
        # compute dewpoint from vapor pressure
        full_Td[k] = mpcalc.dewpoint(full_vp[k] * units.Pa).m

    #check_sfc_hgt(full_z)
    check_latlon(meta_data_dict)

    #############################################
    # CONSTRUCT CLEAN_DATA DICTIONARY
    #############################################
    clean_data = {
        'T': np.array(full_T[0:-2] - 273.15) * units.degC,
        'Td': np.array(np.nan_to_num(full_Td, nan=-80)[0:-2]) * units.degC,
        'z': np.array(full_z[0:-2]) * units.m,
        'p': np.array(full_p[0:-2] / 100) * units.hPa,
        'u': np.array(full_u[0:-2] * 1.94384) * units.kts,
        'v': np.array(full_v[0:-2] * 1.94384) * units.kts}

    # add site & profile info to clean_data dict
    clean_data['site_info'] = {
        'site-id': 'no-site-id',
        'site-name': 'no-site-name',
        'site-lctn': 'no-site-location',
        'site-latlon': meta_data_dict['latlon'],
        'site-elv': meta_data_dict['elev'],
        'source': 'none',
        'model': 'custom',
        'fcst-hour': 'none',
        'run-time': 'none',
        'valid-time': 'none',
        'box_area': f'none'}

    # add plot title to clean_data dict
    clean_data['titles'] = {
        'top_title': meta_data_dict['top_title'],
        'left_title': meta_data_dict['left_title'],
        'right_title': f"{meta_data_dict['right_title']}    "
    }

    return clean_data
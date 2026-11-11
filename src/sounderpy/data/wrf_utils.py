import warnings
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
def parse_time(ds, time_idx=0):
    # convert byte data to string
    byte_data = b''.join(ds['Times'][:].data[time_idx])
    time_str = byte_data.decode('utf-8')

    # replace underscores with spaces to match datetime format
    time_str = time_str.replace('_', ' ')

    # convert to datetime object
    time_obj = datetime.strptime(time_str, '%Y-%m-%d %H:%M:%S')

    return [str(time_obj)[0:4], str(time_obj)[5:7], str(time_obj)[8:10], str(time_obj)[11:16]]




#############################################
# COMPUTE TEMPERATURE FROM WRF OUTPUT
#############################################
def calc_temp(theta_pert, pres):
    '''
        Using theta perturbation from reference theta and full pressure (PB + P) output from WRF,
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

    denominator = (ref_press / (pres)) ** (Rd / cp)  # unitless
    # -------------------------------

    return (numerator / denominator)






#############################################
# COMPUTE DEWPOINT FROM WRF OUTPUT
#############################################
def calc_dwpt(w, p):
    '''
        Using water vapor mixing ratio (QVAPOR, Q2) and pressure output from WRF,

        1. Compute vapor pressure (w/(epsilon + w))*p
        2. Compute dewpoint using Td form of Clausius-Clapeyron

        Td form of Clausius-Clapetron derived from...
        Bolton, D., 1980: The computation of equivalent potential temperature.
            Mon. Wea. Rev., 108, 1046-1053, doi:10.1175/1520-0493(1980)108%3C1046:TCOEPT%3E2.0.CO;2.

    '''

    # CONSTANTS -------------------
    epsilon = 0.622
    # ------------------------------

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
def profile_location(lat_pt, lon_pt, lat_array, lon_array, box_avg_size=0.10):
    # use the full 2-D XLAT/XLONG arrays, since lines of latitude and
    # longitude curve across most WRF map projections
    lat_array = np.asarray(lat_array)
    lon_array = np.asarray(lon_array)

    # nearest grid point to the requested center pt. a nearest pt on the
    # domain edge means the requested pt is outside the domain
    dist = (lat_array - lat_pt) ** 2 + ((lon_array - lon_pt) * np.cos(np.radians(lat_pt))) ** 2
    near_lat_idx, near_lon_idx = np.unravel_index(np.argmin(dist), dist.shape)
    if near_lat_idx in (0, dist.shape[0] - 1) or near_lon_idx in (0, dist.shape[1] - 1):
        raise ValueError(f"{lat_pt}, {lon_pt} is outside of the WRF domain")

    # grid pts within the box, or the nearest pt if the box is smaller than one grid cell
    half = box_avg_size / 2
    mask = (np.abs(lat_array - lat_pt) <= half) & (np.abs(lon_array - lon_pt) <= half)
    if not mask.any():
        mask[near_lat_idx, near_lon_idx] = True
    lat_idx, lon_idx = np.where(mask)

    # return a list of min, max lats & min, max lons, with max as the slice end
    return [lat_idx.min(), lat_idx.max() + 1, lon_idx.min(), lon_idx.max() + 1]







#############################################
# BUILD AND RETURN 'clean_data' DICTIONARY
#############################################
def from_wrf(ds, latlon, model_name='WRF-ARW', run_name='CONTROL RUN', box_avg_size=0.10, time_idx=0):
    '''
        make a SounderPy clean_data profile from WRF output at a location (box-averaged)

        params
            > ds: an opened wrfout file, e.g. netCDF4.Dataset('wrfout_d01_...'); a dataset, not a file path
            > latlon: list of latitude and longitude `[45.56, -100.89]`
            > run_name: str, name of the run for title plotting 'CONTROL'
            > box_avg_size: float, width of the box-average domain in degrees
            > time_idx: int, index of the output time to use, for files with several times
    '''
    t = time_idx

    # define latitude and longitude max/min index values for "bounding box"
    latmn, latmx, lonmn, lonmx = profile_location(latlon[0], latlon[1], ds['XLAT'][t], ds['XLONG'][t], box_avg_size)

    # full pressure is the base state plus the perturbation
    pres = (ds['PB'][t, :, latmn:latmx, lonmn:lonmx] + ds['P'][t, :, latmn:latmx, lonmn:lonmx])

    # geopotential sits on the layer edges; average neighboring edges to the
    # layer middles, where T, moisture, and winds sit
    z_stag = (ds['PHB'][t, :, latmn:latmx, lonmn:lonmx] + ds['PH'][t, :, latmn:latmx, lonmn:lonmx]) / 9.81
    z_mass = 0.5 * (z_stag[:-1] + z_stag[1:])

    # U and V sit on the cell edges; average neighboring edges to the cell centers
    u = 0.5 * (ds['U'][t, :, latmn:latmx, lonmn:lonmx] + ds['U'][t, :, latmn:latmx, lonmn + 1:lonmx + 1])
    v = 0.5 * (ds['V'][t, :, latmn:latmx, lonmn:lonmx] + ds['V'][t, :, latmn + 1:latmx + 1, lonmn:lonmx])
    u10 = ds['U10'][t, latmn:latmx, lonmn:lonmx]
    v10 = ds['V10'][t, latmn:latmx, lonmn:lonmx]

    # WRF winds are grid relative; rotate them to earth relative
    if 'COSALPHA' in ds.variables and 'SINALPHA' in ds.variables:
        cosa = ds['COSALPHA'][t, latmn:latmx, lonmn:lonmx]
        sina = ds['SINALPHA'][t, latmn:latmx, lonmn:lonmx]
        u, v = u * cosa - v * sina, v * cosa + u * sina
        u10, v10 = u10 * cosa - v10 * sina, v10 * cosa + u10 * sina
    else:
        warnings.warn("COSALPHA/SINALPHA not found in this WRF file; winds are left grid relative", Warning)

    # build vertical data dictionary
    vert_data = {
        'vert_T': np.mean(calc_temp(ds['T'][t, :, latmn:latmx, lonmn:lonmx], pres), axis=(1, 2)) - 273.15,
        'vert_p': np.mean(pres, axis=(1, 2)) / 100,
        'vert_z': np.mean(z_mass, axis=(1, 2)),
        'vert_u': np.mean(u, axis=(1, 2)) * 1.94384,
        'vert_v': np.mean(v, axis=(1, 2)) * 1.94384,
        'vert_Td': np.mean(calc_dwpt(ds['QVAPOR'][t, :, latmn:latmx, lonmn:lonmx], pres), axis=(1, 2))
    }

    # build sirface data dictionary
    sfc_data = {
        'sfc_T': np.mean(ds['T2'][t, latmn:latmx, lonmn:lonmx]) - 273.15,
        'sfc_p': np.mean(ds['PSFC'][t, latmn:latmx, lonmn:lonmx]) / 100,
        'sfc_z': np.mean(ds['HGT'][t, latmn:latmx, lonmn:lonmx]),
        'sfc_u': np.mean(u10) * 1.94384,
        'sfc_v': np.mean(v10) * 1.94384,
        'sfc_Td': np.mean(calc_dwpt(ds['Q2'][t, latmn:latmx, lonmn:lonmx], ds['PSFC'][t, latmn:latmx, lonmn:lonmx]))

    }

    # build latlon metadata dictionary
    latlon_data = {
        'data_lat': (ds['XLAT'][t, latmn:latmx, lonmn:lonmx]),
        'data_lon': (ds['XLONG'][t, latmn:latmx, lonmn:lonmx]),
        'data_latnum': latmx - latmn,
        'data_lonnum': lonmx - lonmn,
        'data_time': parse_time(ds, t)
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
        'source': f'MODEL SIMULATION',
        'model': f'{str.upper(model_name)}',
        'fcst-hour': f'{str.upper(run_name)}',
        'run-time': latlon_data['data_time'],
        'valid-time': latlon_data['data_time'],
        'box_area': f'{box_avg_size}° BOX AVG'}

    # add plot title to clean_data dict
    clean_data['titles'] = {
        'top_title': f"{str.upper(model_name)} MODEL SIMULATION VERTICAL PROFILE | {str.upper(run_name)}",
        'left_title': f"VALID: {clean_data['site_info']['valid-time'][1]}/{clean_data['site_info']['valid-time'][2]}/{clean_data['site_info']['valid-time'][0]} {clean_data['site_info']['valid-time'][3]}Z",
        'right_title': f"{clean_data['site_info']['site-latlon'][0]}, {clean_data['site_info']['site-latlon'][1]} | {clean_data['site_info']['box_area']}    "
    }

    return clean_data




# OLD NAME
# make_wrf_profile() was renamed from_wrf() to match SounderPy's other readers (from_sharppy, from_graw, ...). The old
# name still works, with a FutureWarning; remove it in a later release.
def make_wrf_profile(*args, **kwargs):
    """Deprecated: use from_wrf(), which takes the same arguments (ds, latlon, model_name='WRF-ARW', run_name='CONTROL RUN', box_avg_size=0.10, time_idx=0)."""
    warnings.warn("make_wrf_profile() is now from_wrf(); the old name will be removed in a later SounderPy release",
                  FutureWarning, stacklevel=2)
    return from_wrf(*args, **kwargs)

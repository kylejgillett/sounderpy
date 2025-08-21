# IMPORTS SOFTWARE
import numpy as np
from scipy.interpolate import interp1d
from metpy.units import units
from scipy.ndimage import gaussian_filter
from metpy.units import units as mpunits


"""
    SOUNDERPY MERGE & SMOOTH PROFILE FUNCTIONS

    Purpose of module: 

     A function to merge two profiles together using weighted averaging of values along a 
     common height axis to return a single profile.


    (C) KYLE J GILLETT, UNIVERSITY OF NORTH DAKOTA, 2024, 2025
"""


# PROFILE SMOOTHING SCHEME
def smooth_profile(clean_data, sigma=0.5):
    for key in clean_data.keys():
        try:
            unit = clean_data[key].units
            clean_data[key] = gaussian_filter(clean_data[key], sigma)*mpunits(str(unit))
        except:
            pass
    return clean_data



# INTERPOLATION AND AVERAGING SCHEME
def weighted_average_interp(z_a, z_b, var_a, var_b, weight_a=0.5):
    """
    Interpolates two profiles onto a common height grid and computes a weighted average profile.

    Parameters:
    - z_a: Height array for profile A.
    - z_b: Height array for profile B.
    - var_a: Variable array (e.g., temperature, pressure) for profile A.
    - var_b: Variable array for profile B.
    - weight_a: Weight for profile A in the averaging.

    Returns:
    - avg: A weighted average array of the interpolated variables.
    """

    # build a common height grid between both profiles
    z_common = np.linspace(
        max(np.min(z_a), np.min(z_b)),
        min(np.max(z_a), np.max(z_b)),
        num=100
    )

    # interpolate height vars
    interp_func_a = interp1d(z_a, var_a, bounds_error=False, fill_value='extrapolate')
    interp_func_b = interp1d(z_b, var_b, bounds_error=False, fill_value='extrapolate')

    interp_a = interp_func_a(z_common)
    interp_b = interp_func_b(z_common)

    # compute the weighted avgs
    avg = (weight_a * interp_a) + ((1 - weight_a) * interp_b)

    return avg, z_common



def merge_profiles(profile_a, profile_b, weight_a=0.5):
    """
    Averages two profiles using height coordinates.

    Parameters:
    - profile_a: Dictionary containing height and variable arrays for profile A.
    - profile_b: Dictionary containing height and variable arrays for profile B.
    - weight_a: Weight for profile A in the averaging.

    Returns:
    - averaged_profile: Dictionary with weighted-averaged values for z, T, Td, u, v, and p.
    """

    averaged_profile = {}

    unit_list = ['m', 'degC', 'degC', 'kts', 'kts']
    key_list  = ['z', 'T', 'Td', 'u', 'v']

    for key, unit in zip(key_list, unit_list):
        avg, z_common = weighted_average_interp(
            profile_a['z'], profile_b['z'],
            profile_a[key], profile_b[key],
            weight_a
        )
        averaged_profile[key] = avg * units(unit)

    interp_p_a = interp1d(profile_a['z'], profile_a['p'], bounds_error=False, fill_value='extrapolate')
    interp_p_b = interp1d(profile_b['z'], profile_b['p'], bounds_error=False, fill_value='extrapolate')
    interp_p = (weight_a * interp_p_a(z_common)) + ((1 - weight_a) * interp_p_b(z_common))
    averaged_profile['p'] = interp_p * units.hPa

    return averaged_profile

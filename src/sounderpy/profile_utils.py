"""
    SOUNDERPY PROFILE UTILITIES

    Purpose of module:

    Merge and smooth profiles: `merge_profiles()` blends two `clean_data` profiles by weighted
    averaging on a common height axis, and `smooth_profile()` smooths a profile with a Savitzky-Golay
    or Gaussian filter.

    Docs:
        merging profiles:     https://sounderpy.app/docs/reference/helpertools.html#merging-two-profiles-via-weighted-averaging
        smoothing a profile:  https://sounderpy.app/docs/reference/helpertools.html#smooth-a-profile

    References:
        Savitzky, A., and M. J. E. Golay, 1964: Smoothing and differentiation of data by simplified least
            squares procedures. Anal. Chem., 36, 1627-1639, https://doi.org/10.1021/ac60214a047
            (Savitzky-Golay smoothing)
        Virtanen, P., R. Gommers, T. E. Oliphant, and Coauthors, 2020: SciPy 1.0: Fundamental algorithms for
            scientific computing in Python. Nature Methods, 17 (3), 261-272 (Gaussian and Savitzky-Golay
            filters)

    (C) KYLE J GILLETT, UNIVERSITY OF NORTH DAKOTA, 2024, 2025, 2026
"""

# IMPORTS SOFTWARE
import numpy as np
from scipy.interpolate import interp1d
from metpy.units import units
from scipy.ndimage import gaussian_filter
from scipy.signal import savgol_filter
from metpy.units import units as mpunits


# PROFILE SMOOTHING SCHEME
# kwargs - 'method': 'savgol' (Savitzky-Golay, default) or 'gaussian';
# 'sigma': the gaussian kernel's standard deviation (passing it alone picks 'gaussian', as before 3.3.0);
# 'window_length' & 'polyorder': the Savitzky-Golay window (levels, odd) and polynomial order
def smooth_profile(clean_data, sigma=None, method=None, window_length=5, polyorder=2):
    r"""Smooth a profile's data arrays

       Every array with units in ``clean_data`` (p, z, T, Td, u, v, and omega when present) is smoothed; metadata
       is left as is. The smoothed arrays replace the ones in ``clean_data``, which is also returned.

       :param clean_data: a SounderPy ``clean_data`` dictionary
       :type clean_data: dict, required
       :param sigma: standard deviation of the gaussian kernel; larger values smooth more. Passing ``sigma`` without a ``method`` uses the gaussian smoother, as before SounderPy 3.3.0. Default is ``0.5``
       :type sigma: float, optional
       :param method: ``'savgol'`` (or ``'sg'``, a Savitzky-Golay filter, which keeps peaks and inversions sharper) or ``'gaussian'``. Default is ``'savgol'`` (``'gaussian'`` when only ``sigma`` is given)
       :type method: str, optional
       :param window_length: the number of levels in the Savitzky-Golay window (odd); larger values smooth more. Default is ``5``
       :type window_length: int, optional
       :param polyorder: the order of the polynomial fit in each Savitzky-Golay window; lower values smooth more. Default is ``2``
       :type polyorder: int, optional
       :return: the smoothed ``clean_data`` dictionary
       :rtype: dict
    """
    if method is None:
        method = 'gaussian' if sigma is not None else 'savgol'
    if method not in ('savgol', 'sg', 'gaussian'):
        raise ValueError(f"method must be 'savgol' ('sg') or 'gaussian', not '{method}'")
    if sigma is None:
        sigma = 0.5

    for key in clean_data.keys():
        try:
            unit = clean_data[key].units
            values = np.asarray(clean_data[key].m)
            if method in ['savgol', 'sg']:
                # the window can't be longer than the profile, and has to be odd
                window = min(window_length, values.size if values.size % 2 else values.size - 1)
                values = savgol_filter(values, window, min(polyorder, window - 1))
            else:
                values = gaussian_filter(values, sigma)
            clean_data[key] = values*mpunits(str(unit))
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

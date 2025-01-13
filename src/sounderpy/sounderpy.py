### IMPORT SOFTWARE ###
#########################################################################################################
# SOUNDERPY
from typing import List, Tuple, Optional

from .plot import __full_sounding, __full_hodograph, __composite_sounding, __vad_hodograph
from .calc import sounding_params
from .model_data import fetch_model
from .bufkit_data import fetch_bufkit
from .obs_data import fetch_obs
from .acars_data import *
from .wrf_utils import make_wrf_profile
from .cm1_utils import make_cm1_profile
#########################################################################################################





'''
    SOUNDERPY | Vertical Profile Data Retrieval and Analysis Tool For Python
    -------------------------------------------------------------------------
    SounderPy is an open-source atmospheric science Python package for vertical profile analysis. 
    This tool is designed to get data, ‘clean it up’ for simple use, and plot the data on advanced-sounding 
    plots. SounderPy was developed with the goal in mind to keep the code simple and efficient for users of 
    all experience levels and for reliability in all use cases. 


    THIS RELEASE
    -------
    Version: 3.0.8 | January 2025


    DOCUMENTATION
    -------
    Docs: https://kylejgillett.github.io/sounderpy/
    Code: https://github.com/kylejgillett/sounderpy
    PyPi: https://pypi.org/project/sounderpy/
    Operational Site: https://sounderpysoundings.anvil.app/

    COPYRIGHT
    ---------
    Created & maintained by Kyle J Gillett (@wxkylegillett) 2023, 2024, 2025
    
'''




citation_text = f"""
## ---------------------------------- SOUNDERPY ----------------------------------- ##
##          Vertical Profile Data Retrieval and Analysis Tool For Python            ##
##                      v3.0.8 | Jan 2025 | (C) Kyle J Gillett                      ##
##                 Docs: https://kylejgillett.github.io/sounderpy/                  ##
## --------------------- THANK YOU FOR USING THIS PACKAGE! ------------------------ ##
"""
print(citation_text)

#########################################################################################################




#########################################################################
##### GET DATA FUNCTIONS  #####
#########################################################################


###############
# MODEL REANALYSIS DATA
#########################################################################
def get_model_data(model: str, latlon: list, year: str, month: str, day: str, hour: str,
                   dataset: Optional[str] = None, box_avg_size: float = 0.10, hush: bool = False, clean_it: bool = True):

    r"""Get model reanalysis vertical profile data

       Return a ``dict`` of 'cleaned up' model reanalysis data from a given model, for a given location, date, and time

       :param model: the requested model to use (rap-ruc, era5, ncep)
       :type model: str, required
       :param latlon: the latitude & longitude pair for sounding ([44.92, -84.72])
       :type latlon: list, required
       :param year: valid year
       :type year: str, required
       :param month: valid month
       :type month: str, required
       :param day: valid day
       :type day: str, required
       :param hour: required, valid hour
       :type hour: str, required
       :param dataset: optional, target a specific dataset instead of searching for the first one with data.
       :type dataset: str, optional
       :param box_avg_size: optional, determine an area-averaged box size in degrees, default is 0.10 degrees.
       :type box_avg_size: int, optional
       :param hush: whether to 'hush' a read-out of thermodynamic and kinematic parameters when getting a data.
       :type hush: bool, optional, default is `False`
       :param clean_it: whether to return the raw_data object or a clean_data dict.
       :type clean_it: bool, optional, default is `True`


       :return: clean_data, a dict of ready-to-use vertical profile data including pressure, height, temperature, dewpoint, u-wind, v-wind, & model information
       :rtype: dict
    """

    data = fetch_model(model, latlon, year, month, day, hour, dataset, box_avg_size, hush, clean_it)

    return data



###############
# OBSERVED DATA
#########################################################################
def get_obs_data(station, year, month, day, hour, hush=False, clean_it=True):

    r"""Get observed vertical profile data
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

    data = fetch_obs(station, year, month, day, hour, hush, clean_it)

    return data



###############
# BUFKIT DATA
#########################################################################
def get_bufkit_data(model, station, fcst_hour, run_year=None, run_month=None, run_day=None, run_hour=None,
                   hush=False, clean_it=True):

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
       :param clean_it: whether to return the raw_data object or a clean_data dict.
       :type clean_it: bool, optional, default is `True`


       :return: clean_data, a dict of ready-to-use vertical profile data including pressure, height, temperature, dewpoint, u-wind, v-wind, & model information
       :rtype: dict
    """

    data = fetch_bufkit(model, station, fcst_hour, run_year, run_month, run_day, run_hour,
                   hush, clean_it)

    return data



    

############
# ACARS DATA 
#########################################################################













#########################################################################
##### PLOTTING FUNCTIONS  #####
#########################################################################


############
# FULL SOUNDING
#########################################################################
def build_sounding(clean_data, style='full', color_blind=False, dark_mode=False, storm_motion='right_moving',
                   special_parcels=None, show_radar=True, radar_time='sounding', map_zoom=2, modify_sfc=None,
                   show_theta=False, save=False, filename='sounderpy_sounding'):
    
    '''
    Return a full sounding plot of SounderPy data, ``plt``

    :param clean_data: the dictionary of data to be plotted (see :doc:`gettingdata`)
    :type clean_data: dict, required
    :param style: may be `simple` or `full`. Default is `full`.
    :type style: str, optional
    :param color_blind: whether to change the dewpoint trace line from green to blue for improved readability for color deficient users/readers. Default is ``False``
    :type color_blind: bool, optional
    :param dark_mode: ``True`` will invert the color scheme for a 'dark-mode' sounding. Default is ``False``.
    :type dark_mode: bool, optional
    :param storm_motion: the storm motion used for plotting and calculations. Default is 'right_moving'. Custom storm motions are accepted as a `list` of `floats` representing direction and speed. Ex: ``[270.0, 25.0]`` where '270.0' is the *direction in degrees* and '25.0' is the *speed in kts*. See the :ref:`storm_motions` section for more details.
    :type storm_motion: str or list of floats, optional
    :param special_parcels: a nested list of special parcels from the ``ecape_parcels`` library. The nested list should be a list of two lists (`[[a, b], [c, d]]`) where the first list should include 'highlight parcels' and second list should include 'background parcels'. For more details, see the :ref:`parcels_logic` section.
    :type special_parcels: nested `list` of two `lists`, optional
    :param show_theta: bool, optional
    :param save: whether to show the plot inline or save to a file. Default is ``False`` which displays the file inline.
    :type save: bool, optional
    :param filename: the filename by which a file should be saved to if ``save = True``. Default is `sounderpy_sounding`.
    :type filename: str, optional
    :return: plt, a SounderPy sounding built with Matplotlib, MetPy, SharpPy, & SounderPy.
    :rtype: plt
    '''
    
    print(f'> SOUNDING PLOTTER FUNCTION\n  ---------------------------------')

    plt = __full_sounding(clean_data, color_blind, dark_mode, storm_motion, special_parcels, show_radar, radar_time, map_zoom, modify_sfc, show_theta)
    if save:
        plt.savefig(filename, bbox_inches='tight')
    else:
        plt.show()






############
# FULL HODOGRAPH
#########################################################################
def build_hodograph(clean_data, save=False, dark_mode=False, storm_motion='right_moving', sr_hodo=False, modify_sfc=None, filename='sounderpy_hodograph'):
    
    '''
       Return a full sounding plot of SounderPy data, ``plt`` 

       :param clean_data: the dictionary of data to be plotted (see :doc:`gettingdata`)
       :type clean_data: dict, required
       :param save: whether to show the plot inline or save to a file. Default is ``False`` which displays the file inline.
       :type save: bool, optional
       :param filename: the filename by which a file should be saved to if ``save = True``. Default is `sounderpy_sounding`.
       :type filename: str, optional
       :param dark_mode: ``True`` will invert the color scheme for a 'dark-mode' sounding. Default is ``False``.
       :type dark_mode: bool, optional
       :param storm_motion: the storm motion used for plotting and calculations. Default is 'right_moving'. Custom storm motions are accepted as a `list` of `floats` representing direction and speed. Ex: ``[270.0, 25.0]`` where '270.0' is the *direction in degrees* and '25.0' is the *speed in kts*. See the :ref:`storm_motions` section for more details.
       :type storm_motion: str or list of floats, optional
       :param sr_hodo: transform the hodograph from ground relative to storm relative 
       :type sr_hodo: bool, optional, default is ``False``
       :return: plt, a SounderPy sounding built with Matplotlib, MetPy, SharpPy, & SounderPy.
       :rtype: plt
    '''
    
    print(f'> HODOGRAPH PLOTTER FUNCTION --\n-------------------------------')
    
    if save:
        __full_hodograph(clean_data, dark_mode, storm_motion, sr_hodo, modify_sfc).savefig(filename, bbox_inches='tight')
    else:
        __full_hodograph(clean_data, dark_mode, storm_motion, sr_hodo, modify_sfc).show()




############
# COMPOSITE SOUNDING
#########################################################################
def build_composite(data_list, shade_between=True, cmap='viridis', colors_to_use='none', ls_to_use='none', alphas_to_use='none',
                    lw_to_use='none', dark_mode=False, save=False, filename='sounderpy_composite'):

    '''
    Return a composite sounding plot of multiple profiles, ``plt``

    :param data_list: a list of data dictionaries for each profile to be plotted
    :type data_list: list of dicts, required
    :param shade_between: Lightly shade between the dewpoint & temperature trace. In many cases, this improves readability. Default is ``True``.
    :type shade_between: bool, optional
    :param cmap: a linear colormap, may be any custom or matplotlib cmap. Default is 'viridis'. If `colors_to_use` kwarg is provided, `colors_to_use` will be used instead.
    :type cmap: `matplotlib.colors.LinearSegmentedColormap` or `str` representing the name of a matplotlib cmap, optional
    :param colors_to_use: A list of custom matplotlib color name stings. List length must match the number of profiles listed in ``data_list``. Default is 'none'.
    :type colors_to_use: list of strings, optional
    :param alphas_to_use: A list of custom alphas (0.0-1.0). List length must match the number of profiles listed in ``data_list``. Default is 'none'. Default alpha is 1.
    :type alphas_to_use: list of floats, optional
    :param ls_to_use: A list of custom matplotlib line styles. List length must match the number of profiles listed in ``data_list``. Default is 'none'. Default line style is '-'.
    :type ls_to_use: list of stings, optional
    :param lw_to_use: A list of custom line widths. List length must match the number of profiles listed in ``data_list``. Default is 'none'. Default line width is 3.
    :type lw_to_use: list of floats, optional
    :param dark_mode: ``True`` will invert the color scheme for a 'dark-mode' sounding. Default is ``False``.
    :type dark_mode: bool, optional
    :param save: whether to show the plot inline or save to a file. Default is ``False`` which displays the file inline.
    :type save: bool, optional
    :param filename: the filename by which a file should be saved to if ``save = True``. Default is `sounderpy_sounding`.
    :type filename: str, optional
    :return: plt, a SounderPy composite sounding built with Matplotlib, MetPy, SharpPy, & SounderPy.
    :rtype: plt
    '''
    
    print(f'> COMPOSITE SOUNDING FUNCTION\n  -------------------------------')
    
    if save:
        __composite_sounding(data_list, shade_between, cmap, colors_to_use,
            ls_to_use, alphas_to_use, lw_to_use, dark_mode).savefig(filename, bbox_inches='tight')
    else:
        __composite_sounding(data_list, shade_between, cmap, colors_to_use, 
            ls_to_use, alphas_to_use, lw_to_use, dark_mode).show()












############
# VAD HODOGRAPH
#########################################################################
def build_vad_hodograph(vad_data, save=False, dark_mode=False, storm_motion='right_moving', sr_hodo=False,
                        filename='vad_hodograph'):
    '''
       Return a VAD hodograph plot of SounderPy VAD data, ``plt``

       :param vad_data: the dictionary of VAD data to be plotted
       :type vad_data: dict, required
       :param save: whether to show the plot inline or save to a file. Default is ``False`` which displays the file inline.
       :type save: bool, optional
       :param filename: the filename by which a file should be saved to if ``save = True``. Default is `sounderpy_sounding`.
       :type filename: str, optional
       :param dark_mode: ``True`` will invert the color scheme for a 'dark-mode' sounding. Default is ``False``.
       :type dark_mode: bool, optional
       :param storm_motion: the storm motion used for plotting and calculations. Default is 'right_moving'. Custom storm motions are accepted as a `list` of `floats` representing direction and speed. Ex: ``[270.0, 25.0]`` where '270.0' is the *direction in degrees* and '25.0' is the *speed in kts*. See the :ref:`storm_motions` section for more details.
       :type storm_motion: str or list of floats, optional
       :param sr_hodo: transform the hodograph from ground relative to storm relative
       :type sr_hodo: bool, optional, default is ``False``
       :return: plt, a SounderPy sounding built with Matplotlib, MetPy, SharpPy, & SounderPy.
       :rtype: plt
    '''

    print(f'> VAD HODOGRAPH PLOTTER FUNCTION --\n------------------------------------')

    if save:
        __vad_hodograph(vad_data, dark_mode, storm_motion, sr_hodo).savefig(filename, bbox_inches='tight')
    else:
        __vad_hodograph(vad_data, dark_mode, storm_motion, sr_hodo).show()












############
# PRINT DATA TO CONSOLE
#########################################################################
def print_variables(clean_data):
    sounding_params(clean_data, storm_motion='right_moving', modify_sfc=None).print_vals()




    
    

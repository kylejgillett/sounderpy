### IMPORT SOFTWARE ###
#########################################################################################################
# SOUNDERPY
import warnings
from typing import List, Tuple, Optional

from .plot import __full_sounding, __winter_sounding, _resolve_sounding_mode, __full_hodograph, __composite_sounding, __vad_hodograph
from .calc import sounding_params
from .data.model_data import fetch_model
from .data.bufkit_data import fetch_bufkit
from .data.obs_data import fetch_obs
from .data.acars_data import *
from .data.wrf_utils import make_wrf_profile
from .data.cm1_utils import make_cm1_profile
from .data.rapruc_grib_retriever import get_rapruc_data
#########################################################################################################





'''
    SOUNDERPY | Vertical Profile Data Retrieval and Analysis Tool For Python
    -------------------------------------------------------------------------
    SounderPy is an open-source atmospheric science Python package for vertical profile analysis. 
    This tool is designed to get data, 'clean it up' for simple use, and plot the data on advanced-sounding 
    plots. SounderPy was developed with the goal in mind to keep the code simple and efficient for users of 
    all experience levels and for reliability in all use cases. 


    THIS RELEASE
    -------
    Version: 3.3.0 | November 2026


    DOCUMENTATION
    -------
    Docs: https://kylejgillett.github.io/sounderpy/
    Code: https://github.com/kylejgillett/sounderpy
    PyPi: https://pypi.org/project/sounderpy/
    Operational Site: https://sounderpy.app/

    COPYRIGHT
    ---------
    Created & maintained by Kyle J Gillett (@wxkylegillett) 2023, 2024, 2025, 2026
    
'''




citation_text = f"""
## ---------------------------------- SOUNDERPY ----------------------------------- ##
##          Vertical Profile Data Retrieval and Analysis Tool For Python            ##
##                   v3.3.0 | November 2026 | (C) Kyle J Gillett                    ##
##                 Docs: https://kylejgillett.github.io/sounderpy/                  ##
## --------------------- THANK YOU FOR USING THIS PACKAGE! ------------------------ ##
"""
#print(citation_text)

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
       :param clean_it: whether to return the raw_data object or a clean_data dict. With ``False``, RAP/RUC returns a dict of xarray datasets ('pressure', '2m', '10m', 'surface') cut to the box-average domain, ERA5 returns an xarray dataset, and NCEP-FNL and RAP-NOW return a netCDF4 dataset.
       :type clean_it: bool, optional, default is `True`


       :return: clean_data, a dict of ready-to-use vertical profile data including pressure, height, temperature, dewpoint, u-wind, v-wind, & model information
       :rtype: dict
    """

    if model.casefold() in ['rap', 'ruc', 'rap-ruc']:
        data = get_rapruc_data(latlon, year, month, day, hour, dataset, box_avg_size, hush, clean_it)
    else:
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
                   hush=False, clean_it=True, all_fcst_hours=False):

    r"""Get BUFKIT forecast model vertical profile data
       Return a ``dict`` of 'cleaned up' model forecast data from a given model, for a given BUFKIT site identifier, forecast hour, & model-run-date

       :param model: the requested model to use (such as hrrr, rap, gfs, etc)
       :type model: str, required
       :param station: a 3-4 digit BUFKIT site identifier
       :type station: str, required
       :param fcst_hour: the forecast hour; the profile is valid at the model run time plus this many hours
       :type fcst_hour: int, required
       :param run_year: model run year. Leave the run date as ``None`` for the latest run.
       :type run_year: str, optional, default is `None`
       :param run_month: model run month
       :type run_month: str, optional, default is `None`
       :param run_day: model run day
       :type run_day: str, optional, default is `None`
       :param run_hour: model run hour (the initialization time, not the valid time)
       :type run_hour: str, optional, default is `None`
       :param hush: whether to 'hush' a read-out of thermodynamic and kinematic parameters when getting a data.
       :type hush: bool, optional, default is `False`
       :param clean_it: whether to return the raw_data object or a clean_data dict.
       :type clean_it: bool, optional, default is `True`
       :param all_fcst_hours: return a list with every forecast hour in the BUFKIT file, from one download. ``fcst_hour`` is then ignored. Not available for ``sref`` and ``hiresw``, whose files hold several ensemble members.
       :type all_fcst_hours: bool, optional, default is `False`


       :return: clean_data, a dict of ready-to-use vertical profile data including pressure, height, temperature, dewpoint, u-wind, v-wind, & model information, or a list of them, one per forecast hour, with ``all_fcst_hours``
       :rtype: dict or list of dicts
    """

    data = fetch_bufkit(model, station, fcst_hour, run_year, run_month, run_day, run_hour,
                   hush, clean_it, all_fcst_hours)

    return data




#########################################################################
##### PLOTTING FUNCTIONS  #####
#########################################################################


############
# FULL SOUNDING
#########################################################################
def build_sounding(clean_data, color_blind=False, dark_mode=False, storm_motion='right_moving',
                   special_parcels=None, radar='mosaic', radar_time='sounding', map_zoom=2, modify_sfc=None,
                   show_theta=False, hodo_boundary=None, dpi=100, save=False, filename='sounderpy_sounding', max_levels=250, ylim=None, theme=None,
                   sr_hodo=False, mode='severe', show_obs=False):
    
    '''
    Return a full sounding plot of SounderPy data, ``plt``

    :param clean_data: the dictionary of data to be plotted (see :doc:`gettingdata`)
    :type clean_data: dict, required
    :param color_blind: whether to change the dewpoint trace line from green to blue for improved readability for color deficient users/readers. Default is ``False``
    :type color_blind: bool, optional
    :param dark_mode: ``True`` will invert the color scheme for a 'dark-mode' sounding. Default is ``False``.
    :type dark_mode: bool, optional
    :param storm_motion: the storm motion used for plotting and calculations. Default is 'right_moving'. Custom storm motions are accepted as a `list` of `floats` representing direction and speed. Ex: ``[270.0, 25.0]`` where '270.0' is the *direction in degrees* and '25.0' is the *speed in kts*. See the :ref:`storm_motions` section for more details.
    :type storm_motion: str or list of floats, optional
    :param special_parcels: a nested list of special parcels from the ``ecape_parcels`` library. The nested list should be a list of two lists (`[[a, b], [c, d]]`) where the first list should include 'highlight parcels' and second list should include 'background parcels'. For more details, see the :ref:`parcels_logic` section.
    :type special_parcels: nested `list` of two `lists`, optional
    :param radar: type of radar data to retrieve for radar inset, ``mosaic`` will plot a CONUS reflectivity mosaic (data only available within one month of current time), ``single`` will plot reflectivity from the nearest NEXRAD WRS-88D radar site (data availability depends on site). Default is ``mosaic``.
    :type radar: str, optional
    :param radar_time: radar mosaic data valid time. May be ``sounding`` (uses the valid time of the sounding data), or ``now`` (current time/date). **Note: radar mosaic data only goes back 1 month from current date**
    :type radar_time: str, optional, Default is ``sounding``.
    :param map_zoom: a 'zoom' level for the map inset as an `int`. **Note: Setting ``map_zoom=0`` will hide the map**
    :type map_zoom: int, optional, Default is ``2``.
    :param modify_sfc: a `dict` in the format ``{'T': 25, 'Td': 21, 'ws': 20, 'wd': 270}`` to modify the surface values of the ``clean_data`` dict. See the :ref:`modify_sfc` section for more details.
    :type modify_sfc: None or dict, optional, default is None
    :param show_theta: whether to replace the piecewise CAPE plot with theta & theta-e. Default is ``False``
    :type show_theta: bool, optional
    :param hodo_boundary: plot a "boundary" (straight line) on the hodograph axis. Do so by providing an angle and a color. The angle represents the angle between zero-degrees (north) and the upper-half of the boundary line, to the right. Multiple boundaries can be plotted. Default is ``None``
    :type hodo_boundary: dict of lists -- structured as ``{'angle':[], 'color':[]}``, optional
    :param dpi: Dots per inch resolution of the figure, default is 100.
    :type dpi: int, optional
    :param save: whether to show the plot inline or save to a file. Default is ``False`` which displays the file inline.
    :type save: bool, optional
    :param filename: the filename by which a file should be saved to if ``save = True``. Default is `sounderpy_sounding`.
    :type filename: str, optional
    :param max_levels: the most vertical levels used to calculate and plot the sounding. Profiles with more levels are thinned to this many evenly spaced levels. Set to ``None`` to use every level. Default is ``250``.
    :type max_levels: int or None, optional
    :param ylim: the skew-t's pressure range as (bottom, top) in hPa, such as ``(1000, 700)``, to look at one part of the atmosphere. The temperature axis and pressure ticks fit to that range, while the rest of the figure still uses the full profile. Default is ``None``, which shows the full profile.
    :type ylim: tuple of two numbers, optional
    :param theme: changes to the figure's colors, such as ``{'border': 'k', 'temperature': 'navy', 'parcel': {'mu': 'black'}}``. An invalid name lists the valid ones. Default is ``None``.
    :type theme: dict, optional
    :param sr_hodo: draw the hodograph storm relative, shifting the winds by the storm motion from ``storm_motion``. If the profile has no storm motion, the hodograph stays ground relative. Default is ``False``.
    :type sr_hodo: bool, optional
    :param mode: the kind of sounding to draw: ``'severe'``, ``'winter'`` (precipitation type, snow growth, and winter hazards), or ``'auto'``, which draws the winter sounding when the surface wet-bulb temperature is below 35 degrees F and the severe sounding otherwise. The winter sounding ignores ``special_parcels``, ``show_theta``, ``sr_hodo`` and ``hodo_boundary``, with a warning. Default is ``'severe'``.
    :type mode: str, optional
    :param show_obs: winter sounding only: draw the surface observation (METAR) nearest the sounding in space and time as a station model in the hodograph's corner, from the Iowa Environmental Mesonet archive. Needs an internet connection; nothing is drawn when no station within 75 km has a report within an hour of the valid time. Default is ``False``.
    :type show_obs: bool, optional
    :return: plt, a SounderPy sounding built with Matplotlib, MetPy, SharpPy, & SounderPy.
    :rtype: plt
    '''
    
    print(f'> SOUNDING PLOTTER FUNCTION\n  ---------------------------------')

    mode = _resolve_sounding_mode(mode, clean_data, modify_sfc, max_levels)
    if mode == 'winter':
        # severe-only settings changed from their defaults are ignored, with a warning
        ignored = [name for name, value, default in (('special_parcels', special_parcels, None), ('show_theta', show_theta, False),
                                                     ('sr_hodo', sr_hodo, False), ('hodo_boundary', hodo_boundary, None))
                   if value != default]
        if ignored:
            warnings.warn(f"The winter sounding doesn't use {', '.join(ignored)}; ignoring {'it' if len(ignored) == 1 else 'them'}.", Warning)
        plt = __winter_sounding(clean_data, color_blind, dark_mode, storm_motion, radar, radar_time, map_zoom, modify_sfc, dpi, max_levels, ylim, theme, show_obs)
    else:
        if show_obs:
            warnings.warn("show_obs is only used by the winter sounding (mode='winter'); ignoring it.", Warning)
        plt = __full_sounding(clean_data, color_blind, dark_mode, storm_motion, special_parcels, radar, radar_time, map_zoom, modify_sfc, show_theta, hodo_boundary, dpi, max_levels, ylim, theme, sr_hodo)
    if save:
        plt.savefig(filename, bbox_inches='tight')
    else:
        plt.show()






############
# FULL HODOGRAPH
#########################################################################
def build_hodograph(clean_data, dark_mode=False, storm_motion='right_moving', sr_hodo=False, modify_sfc=None,
                    radar='mosaic', radar_time='sounding', map_zoom=2, hodo_boundary=None, dpi=100, save=False, filename='sounderpy_hodograph', max_levels=250, theme=None):

    '''
       Return a hodograph plot of SounderPy data, ``plt`` 

       :param clean_data: the dictionary of data to be plotted (see :doc:`gettingdata`)
       :type clean_data: dict, required
       :param dark_mode: ``True`` will invert the color scheme for a 'dark-mode' sounding. Default is ``False``.
       :type dark_mode: bool, optional
       :param storm_motion: the storm motion used for plotting and calculations. Default is 'right_moving'. Custom storm motions are accepted as a `list` of `floats` representing direction and speed. Ex: ``[270.0, 25.0]`` where '270.0' is the *direction in degrees* and '25.0' is the *speed in kts*. See the :ref:`storm_motions` section for more details.
       :type storm_motion: str or list of floats, optional
       :param sr_hodo: transform the hodograph from ground relative to storm relative 
       :type sr_hodo: bool, optional, default is ``False``
       :param radar: type of radar data to retrieve for radar inset, ``mosaic`` will plot a CONUS reflectivity mosaic (data only available within one month of current time), ``single`` will plot reflectivity from the nearest NEXRAD WRS-88D radar site (data availability depends on site). Default is ``mosaic``.
       :type radar: str, optional
       :param radar_time: radar mosaic data valid time. May be ``sounding`` (uses the valid time of the sounding data), or ``now`` (current time/date). **Note: radar mosaic data only goes back 1 month from current date**
       :type radar_time: str, optional, Default is ``sounding``.
       :param map_zoom: a 'zoom' level for the map inset as an `int`. **Note: Setting ``map_zoom=0`` will hide the map**
       :type map_zoom: int, optional, Default is ``2``.
       :param modify_sfc: a `dict` in the format ``{'T': 25, 'Td': 21, 'ws': 20, 'wd': 270}`` to modify the surface values of the ``clean_data`` dict. See the :ref:`modify_sfc` section for more details.
       :type modify_sfc: None or dict, optional, default is None
       :param hodo_boundary: plot a "boundary" (straight line) on the hodograph axis. Do so by providing an angle and a color. The angle represents the angle between zero-degrees (north) and the upper-half of the boundary line, to the right. Multiple boundaries can be plotted. Default is ``None``
       :type hodo_boundary: dict of lists -- structured as ``{'angle':[], 'color':[]}``, optional
       :param dpi: Dots per inch resolution of the figure, default is 100.
       :type dpi: int, optional
       :param save: whether to show the plot inline or save to a file. Default is ``False`` which displays the file inline.
       :type save: bool, optional
       :param filename: the filename by which a file should be saved to if ``save = True``. Default is `sounderpy_hodograph`.
       :type filename: str, optional
       :param max_levels: the most vertical levels used to calculate and plot the hodograph. Profiles with more levels are thinned to this many evenly spaced levels. Set to ``None`` to use every level. Default is ``250``.
       :type max_levels: int or None, optional
       :param theme: changes to the figure's colors, such as ``{'border': 'k', 'parcel': {'ml': 'black'}}``. An invalid name lists the valid ones. Default is ``None``.
       :type theme: dict, optional
       :return: plt, a SounderPy hodograph built with Matplotlib, MetPy, SharpPy, & SounderPy.
       :rtype: plt
    '''


    print(f'> HODOGRAPH PLOTTER FUNCTION --\n-------------------------------')

    plt = __full_hodograph(clean_data, dark_mode, storm_motion, sr_hodo, modify_sfc, radar, radar_time, map_zoom, hodo_boundary, dpi, max_levels, theme)
    if save:
        plt.savefig(filename, bbox_inches='tight')
    else:
        plt.show()




############
# COMPOSITE SOUNDING
#########################################################################
def build_composite(data_list, shade_between=True, cmap='viridis', colors=None, linestyles=None, alphas=None,
                    linewidths=None, dark_mode=False, save=False, filename='sounderpy_composite',
                    show_params=True, storm_motion='right_moving', max_levels=250, dpi=100, ylim=None,
                    colors_to_use=None, ls_to_use=None, alphas_to_use=None, lw_to_use=None, theme=None, style='full', sr_hodo=False):

    '''
    Return a composite sounding plot of multiple profiles, ``plt``

    :param data_list: a list of data dictionaries for each profile to be plotted. Profiles may come from any mix of sources, such as observed, model, BUFKIT, ACARS, and custom data.
    :type data_list: list of dicts, required
    :param shade_between: Lightly shade between the dewpoint & temperature trace. In many cases, this improves readability. Default is ``True``.
    :type shade_between: bool, optional
    :param cmap: a linear colormap, may be any custom or matplotlib cmap. Default is 'viridis', using only the part of it that stands out from the figure's background. If `colors` kwarg is provided, `colors` will be used instead.
    :type cmap: `matplotlib.colors.LinearSegmentedColormap` or `str` representing the name of a matplotlib cmap, optional
    :param colors: A list of custom matplotlib color name stings. List length must match the number of profiles listed in ``data_list``. Default is 'none'.
    :type colors: list of strings, optional
    :param alphas: A list of custom alphas (0.0-1.0). List length must match the number of profiles listed in ``data_list``. Default is 'none'. Default alpha is 1.
    :type alphas: list of floats, optional
    :param linestyles: A list of custom matplotlib line styles. List length must match the number of profiles listed in ``data_list``. Default is 'none'. Default line style is '-'.
    :type linestyles: list of stings, optional
    :param linewidths: A list of custom line widths. List length must match the number of profiles listed in ``data_list``. Default is 'none'. Default line width is 3.
    :type linewidths: list of floats, optional
    :param dark_mode: ``True`` will invert the color scheme for a 'dark-mode' sounding. Default is ``False``.
    :type dark_mode: bool, optional
    :param save: whether to show the plot inline or save to a file. Default is ``False`` which displays the file inline.
    :type save: bool, optional
    :param filename: the filename by which a file should be saved to if ``save = True``. Default is `sounderpy_sounding`.
    :type filename: str, optional
    :param show_params: whether to compute each profile's parameters and show them in the table. Setting this to ``False`` is faster. Default is ``True``.
    :type show_params: bool, optional
    :param storm_motion: the storm motion used to compute every profile's parameters. Default is ``right_moving``.
    :type storm_motion: str or list of floats, optional
    :param max_levels: the most vertical levels used to calculate and plot each profile. Profiles with more levels are thinned to this many evenly spaced levels. Set to ``None`` to use every level. Default is ``250``.
    :type max_levels: int or None, optional
    :param dpi: dots per inch resolution of the figure. Default is ``100``.
    :type dpi: int, optional
    :param ylim: the skew-t's pressure range as (bottom, top) in hPa, such as ``(1000, 700)``, to look at one part of the atmosphere. The temperature axis and pressure ticks fit to that range, while the rest of the figure still uses the full profile. Default is ``None``, which shows the full profile.
    :type ylim: tuple of two numbers, optional
    :param colors_to_use: deprecated, use ``colors``
    :param ls_to_use: deprecated, use ``linestyles``
    :param alphas_to_use: deprecated, use ``alphas``
    :param lw_to_use: deprecated, use ``linewidths``
    :param theme: changes to the figure's colors, such as ``{'border': 'k', 'temperature': 'navy', 'parcel': {'mu': 'black'}}``. An invalid name lists the valid ones. Default is ``None``.
    :type theme: dict, optional
    :param style: the panels to show. ``'full'`` shows the skew-t and hodograph, ``'sounding'`` shows only the skew-t, and ``'hodograph'`` shows only the hodograph. Each style's table lists the parameters that go with its panels. Default is ``'full'``.
    :type style: str, optional
    :param sr_hodo: draw the hodograph storm relative, shifting each profile's winds by its own storm motion from ``storm_motion``. Not used with ``style='sounding'``. If any profile has no storm motion, the hodograph stays ground relative. Default is ``False``.
    :type sr_hodo: bool, optional
    :return: plt, a SounderPy composite sounding built with Matplotlib, MetPy, SharpPy, & SounderPy.
    :rtype: plt
    '''
    
    print(f'> COMPOSITE SOUNDING FUNCTION\n  -------------------------------')

    # the old argument names still work for now
    def _use_new_name(old_name, old_value, new_name, new_value):
        if old_value is None:
            return new_value
        warnings.warn(f"'{old_name}' is deprecated and will be removed in a future release. Use '{new_name}' instead.", FutureWarning, stacklevel=3)
        return old_value if new_value is None else new_value

    colors = _use_new_name('colors_to_use', colors_to_use, 'colors', colors)
    linestyles = _use_new_name('ls_to_use', ls_to_use, 'linestyles', linestyles)
    alphas = _use_new_name('alphas_to_use', alphas_to_use, 'alphas', alphas)
    linewidths = _use_new_name('lw_to_use', lw_to_use, 'linewidths', linewidths)

    plt =  __composite_sounding(data_list, shade_between, cmap, colors,
            linestyles, alphas, linewidths, dark_mode, show_params, storm_motion, max_levels, dpi, ylim, theme, style, sr_hodo)
    if save:
        plt.savefig(filename, bbox_inches='tight')
    else:
        plt.show()












############
# VAD HODOGRAPH
#########################################################################
def build_vad_hodograph(vad_data, save=False, dark_mode=False, storm_motion='right_moving', sr_hodo=False,
                        filename='vad_hodograph', radar='mosaic', radar_time='sounding', map_zoom=2, hodo_boundary=None, dpi=100, theme=None):
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
       :param radar: the type of radar on the map: ``'mosaic'``, ``'single-site'``, or ``None``. Default is ``'mosaic'``.
       :type radar: str or None, optional
       :param radar_time: the time of radar to show: ``'sounding'`` for the data's valid time, or ``'now'``. Default is ``'sounding'``.
       :type radar_time: str, optional
       :param map_zoom: a zoom level for the map inset. Setting ``map_zoom=0`` hides the map. Default is ``2``.
       :type map_zoom: int, optional
       :param hodo_boundary: plot a boundary on the hodograph, structured as ``{'angle':[], 'color':[]}``. Default is ``None``.
       :type hodo_boundary: dict of lists, optional
       :param dpi: dots per inch resolution of the figure. Default is ``100``.
       :type dpi: int, optional
       :param theme: changes to the figure's colors, such as ``{'border': 'k', 'temperature': 'navy', 'parcel': {'mu': 'black'}}``. An invalid name lists the valid ones. Default is ``None``.
       :type theme: dict, optional
       :return: plt, a SounderPy sounding built with Matplotlib, MetPy, SharpPy, & SounderPy.
       :rtype: plt
    '''

    print(f'> VAD HODOGRAPH PLOTTER FUNCTION --\n------------------------------------')

    plt =  __vad_hodograph(vad_data, dark_mode, storm_motion, sr_hodo, radar, radar_time, map_zoom, hodo_boundary, dpi, theme)
    if save:
        plt.savefig(filename, bbox_inches='tight')
    else:
        plt.show()












############
# PRINT DATA TO CONSOLE
#########################################################################
def print_variables(clean_data):
    sounding_params(clean_data, storm_motion='right_moving', modify_sfc=None).print_vals()




    
    

"""
    SOUNDERPY OPERATIONAL FIGURES

    Purpose of module:

    House tools for plotting sounding, hodograph, and composite sounding
    figures from data provided by the user. Functions here are referenced
    by sounderpy.py.

    (C) KYLE J GILLETT, UNIVERSITY OF NORTH DAKOTA, 2024, 2026
"""

# METPY OPERATIONS

# MATPLOTLIB OPERATIONS
import matplotlib.pyplot as plt

# NUMPY OPERATIONS
import numpy as np
import numpy.ma as ma

# STANDARD PYTHON
import warnings
import time

# FOR MAPPING
import cartopy.crs as ccrs

# SOUNDERPY MODULES
from ..calc import sounding_params, vad_params
from ..radar_utils import *

from .style import (FONT_SIZES,
                    COMPOSITE_LAYOUTS,
                    FULL_HODOGRAPH_LAYOUT,
                    FULL_SOUNDING_LAYOUT,
                    _get_plot_theme)
from .components import (_clip_skewt_text,
                         _get_height_markers,
                         _setup_zoomed_skewt,
                         _validate_pressure_limits,
                         STABILITY_PANEL_DEPTH,
                         _add_figure_margin,
                         _calculate_special_parcel,
                         _calculate_special_parcels,
                         _get_composite_titles,
                         _get_hodograph_bounds,
                         _get_missing_storm_motions,
                         _get_profile_colors,
                         _get_profile_letter,
                         _get_profile_styles,
                         _get_profile_winds,
                         _get_skewt_temperature_bounds,
                         _make_storm_relative,
                         _make_winds_storm_relative,
                         _plot_basic_parcels,
                         _plot_cape_cin_panel,
                         _plot_composite_hodograph_key,
                         _plot_composite_hodograph_traces,
                         _plot_composite_skewt_profiles,
                         _plot_composite_table,
                         _plot_credits,
                         _plot_downdraft_parcel,
                         _plot_figure_titles,
                         _plot_hodograph_boundary,
                         _plot_hodograph_height_markers,
                         _plot_hodograph_motion_labels,
                         _plot_hodograph_motion_readout,
                         _plot_hodograph_parameter_box,
                         _plot_hodograph_sm_vector,
                         _plot_hodograph_srh,
                         _plot_hodograph_trace,
                         _plot_hodograph_velocity_labels,
                         _plot_indices_box,
                         _plot_kinematics_table,
                         _plot_logo,
                         _plot_map_mosaic_radar,
                         _plot_map_profile_location,
                         _plot_map_single_site_radar,
                         _plot_moisture_table,
                         _plot_omega,
                         _plot_skewt_accessories,
                         _plot_skewt_cape_annotations,
                         _plot_skewt_freezing_annotations,
                         _plot_skewt_height_annotations,
                         _plot_skewt_layer_annotations,
                         _plot_skewt_legend,
                         _plot_skewt_parcel_level_annotations,
                         _plot_skewt_profiles,
                         _plot_skewt_rh_annotations,
                         _plot_skewt_surface_annotations,
                         _plot_skewt_wind_barbs,
                         _plot_special_parcels,
                         _plot_sr_wind_panel,
                         _plot_streamwiseness_panel,
                         _plot_temperature_advection,
                         _plot_thermo_table,
                         _plot_theta_panel,
                         _plot_vorticity_panel,
                         _plot_wind_component_panel,
                         _plot_wind_parameter_box,
                         _prepare_hodograph_data,
                         _prepare_sounding_data,
                         _setup_hodograph,
                         _setup_map_inset,
                         _setup_skewt)


#########################################################################
############################## FIGURES ##################################
#########################################################################




#########################################################################
########################## FULL SOUNDING ################################

def __full_sounding(clean_data, color_blind, dark_mode, storm_motion, special_parcels=None, 
                    radar='mosaic', radar_time='sounding', map_zoom=2, modify_sfc=None,
                    show_theta=False, hodo_boundary=None, dpi=100, max_levels=250, ylim=None, theme_overrides=None):
    
    # record process time 
    st = time.time()

    #################################################################
    ### SET UP PLOT THEME ###
    #################################################################
    theme = _get_plot_theme(dark_mode=dark_mode, colorblind=color_blind, overrides=theme_overrides)
    layout = FULL_SOUNDING_LAYOUT
    

    #################################################################
    ### SET UP THE DATA ###
    #################################################################
    print(f'    > PREPARING DATA')
    sounding_data = _prepare_sounding_data(clean_data, modify_sfc, max_levels)
    surface_modified = modify_sfc is not None
    
    # calculate other sounding parameters using SounderPy Calc
    snd_params = sounding_params(sounding_data, storm_motion, include_all_parcels=not show_theta, parcel_depth=STABILITY_PANEL_DEPTH).calc()
    # snd_params = {"general": general, "thermo": thermo, "kinem": kinem,"intrp": intrp}
    #################################################################
    
    

    #########################################################################
    ################################ SKEW-T ################################# 
    #########################################################################
    
    
    #################################################################
    ### BUILD FIGURE AND SKEW-T ###
    #################################################################
    print(f'    > BUILDING SKEW-T LOG-P DIAGRAM')
    fig = plt.figure(figsize=(22, 13), linewidth=10, edgecolor=theme["border"], facecolor=theme["background"], dpi=dpi)        
    # create skew-t 
    ylim = _validate_pressure_limits(ylim)
    if ylim is None:
        skew = _setup_skewt(fig, sounding_data["T"], theme, rect=layout["skewt"])
    else:
        # leave room on the right for the level labels and wind barbs
        # keep the parcel traces in view too
        thermo = snd_params["thermo"]
        parcel_traces = [(thermo[f"{key}P_trace"], thermo[f"{key}T_trace"]) for key in ("sb", "ml", "mu") if not ma.is_masked(thermo[f"{key}cape"]) and thermo[f"{key}cape"] > 0]
        skew = _setup_zoomed_skewt(fig, [sounding_data], theme, ylim, layout["skewt"], right_margin=0.125, extra_traces=parcel_traces)
    # add skew-t extras (adiabats, mixing ratio, highlighted isotherms, shading)
    _plot_skewt_accessories(skew, theme)


    #################################################################
    ### PLOT SKEW-T DATA ###
    #################################################################
    # plot skew-t profiles (T, Td, Tv, Tw)
    _plot_skewt_profiles(skew, sounding_data, snd_params, theme)


    #################################################################
    ### PARCELS LOGIC ###
    #################################################################
    '''
    SounderPy 'Special Parcels Logic' using Amelia Urquhart's ecape_parcels library: https://github.com/a-urq/ecape-parcel-py
    
    Docs: https://kylejgillett.github.io/sounderpy/reference/plottingdata.html#parcel-logic
    '''

    # user-requested special parcels
    if isinstance(special_parcels, list):

        print("    > COMPUTING & PLOTTING SPECIAL PARCELS -- this may take a moment")

        parcel_data = _calculate_special_parcels(special_parcels, sounding_data, snd_params)
        _plot_special_parcels(skew, special_parcels, parcel_data, theme)

    # default/simple parcels
    elif special_parcels in (None, "simple"):
        muecape_trace = None
        if special_parcels is None and snd_params["thermo"]["mu_ecape"] > 0:
            muecape_trace = _calculate_special_parcel("mu_ps_ecape", sounding_data, snd_params)
        _plot_basic_parcels(skew, snd_params, theme, muecape_trace=muecape_trace)

    # downdraft parcel is always plotted
    _plot_downdraft_parcel(skew, snd_params, theme)
    
    
    #################################################################
    ### PLOT SKEW-T ANNOTATIONS ###
    #################################################################
    # height labels and surface elevation
    if ylim is None:
        _plot_skewt_height_annotations(skew, sounding_data, snd_params, theme)
    else:
        _plot_skewt_height_annotations(skew, sounding_data, snd_params, theme, heights=_get_height_markers(snd_params, ylim))

    # surface T/Td and modified-surface flag
    _plot_skewt_surface_annotations(skew, sounding_data, theme, surface_modified=surface_modified,)

    # LCL, LFC, EL, MPL
    _plot_skewt_parcel_level_annotations(skew, sounding_data, snd_params, theme, parcel="ml", x=0.83)
    _plot_skewt_parcel_level_annotations(skew, sounding_data, snd_params, theme, parcel="mu", levels=("lcl", "lfc", "el", "mpl"), x=0.90)

    # freezing and wet-bulb-zero levels
    _plot_skewt_freezing_annotations(skew, snd_params, theme)

    # CAPE labels
    _plot_skewt_cape_annotations(skew, snd_params, theme)

    # PBL, lapse rate, EIL, DGZ/HGZ
    _plot_skewt_layer_annotations(skew, sounding_data, snd_params, theme)

    # RH labels
    _plot_skewt_rh_annotations(skew, snd_params,theme)


    #################################################################
    ### PLOT SKEW-T WIND BARBS ###
    #################################################################
    _plot_skewt_wind_barbs(skew,sounding_data, snd_params, theme, add_sr_barbs=True, ylim=ylim)
    if ylim is not None:
        _clip_skewt_text(skew)

    


    #################################################################
    ### ADD HODOGRAPH ###
    #################################################################
    print("    > BUILDING HODOGRAPH")
    
    hodo_data = _prepare_hodograph_data(snd_params, max_height=12000,)
    hodo_xlim, hodo_ylim = _get_hodograph_bounds(hodo_data)
    hodo = _setup_hodograph(fig.add_axes(layout["hodograph"]), theme, hodo_xlim, hodo_ylim)

    _plot_hodograph_velocity_labels(hodo, theme)
    _plot_hodograph_height_markers(hodo, snd_params, theme)
    _plot_hodograph_trace(hodo, snd_params, theme)

    _plot_hodograph_motion_labels(hodo, snd_params, storm_motion, theme)
    _plot_hodograph_sm_vector(hodo, snd_params, theme)
    _plot_hodograph_motion_readout(hodo, snd_params, storm_motion, theme)
    _plot_hodograph_srh(hodo, snd_params, theme, layer="auto")
    
    if hodo_boundary is not None:
        _plot_hodograph_boundary(hodo, hodo_boundary,snd_params, theme, show_mw_vector=True)
        
    
    
    #########################################################################
    ############################ OTHER AXES ################################# 
    #########################################################################

    

    #################################################################
    ### SOUNDING MAP INSET ###
    #################################################################

    if map_zoom > 0:
        print("    > BUILDING MAP INSET")
        latlon = clean_data["site_info"]["site-latlon"]
        map_ax = fig.add_axes(layout["map"], projection=ccrs.PlateCarree(), zorder=12)
        _setup_map_inset(map_ax, latlon, theme, map_zoom=map_zoom)
        
        if radar == "mosaic":
            _plot_map_mosaic_radar(map_ax, sounding_data, map_zoom, radar_time, theme)
        
        elif radar in ("single-site", "single"):
            _plot_map_single_site_radar(map_ax, clean_data, radar_time, theme)

        _plot_map_profile_location(map_ax, latlon, theme)
        


    #################################################################
    ### TEMPERATURE ADVECTION ###
    #################################################################
    print(f'    > BUILDING ACCESSORY PLOTS')
    _plot_temperature_advection(fig.add_axes(layout["temp_adv"]), snd_params, theme, pressure_limits=ylim if ylim is not None else (1050, 100))


    #################################################################
    ### OMEGA ###
    #################################################################
    if 'omega' in sounding_data:
        _plot_omega(fig.add_axes(layout["omega"], zorder=10), sounding_data, theme, pressure_limits=ylim if ylim is not None else (1050, 100))


    #################################################################
    ### LOW-LEVEL PANELS ###
    #################################################################
    _plot_streamwiseness_panel(fig.add_axes(layout["streamwiseness"]), snd_params, theme)
    _plot_vorticity_panel(fig.add_axes(layout["vorticity"]), snd_params, theme)
    _plot_sr_wind_panel(fig.add_axes(layout["sr_wind"]), snd_params, theme)

    stability_ax = fig.add_axes(layout["stability"])
    if show_theta:
        _plot_theta_panel(stability_ax, snd_params, theme)
    else:
        _plot_cape_cin_panel(stability_ax, snd_params, theme)



    #########################################################################
    ############################## TEXT PLOTS ###############################
    #########################################################################
    print(f'    > PLOTTING SOUNDING PARAMETERS')
    _plot_thermo_table(fig.add_axes(layout["thermo_table"], zorder=3), snd_params, theme)
    _plot_moisture_table(fig.add_axes(layout["moisture_table"], zorder=3), snd_params, theme)
    _plot_kinematics_table(fig.add_axes(layout["kinematics_table"], zorder=3), snd_params, theme)
    _plot_indices_box(fig.add_axes(layout["indices_box"], zorder=3), snd_params, theme)




    #########################################################################
    ############################# PLOT EXTRAS ############################### 
    ######################################################################### 
    print(f'    > FIGURE LOADING...')
    _plot_figure_titles(fig, clean_data["titles"], theme)
    _plot_skewt_legend(skew, theme)
    _plot_credits(fig, theme)
    _plot_logo(fig.add_axes(layout["logo"], anchor="SE", zorder=3))
    _add_figure_margin(fig, layout["margin"])

    
    print('    > COMPLETE --------')
    elapsed_time = time.time() - st
    print('    > RUNTIME:', time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
    
    return plt
#########################################################################




















#########################################################################
############################ HODOGRAPH ##################################

def _build_hodograph_figure(snd_params, site_data, titles, dark_mode, storm_motion, sr_hodo,
                            radar, radar_time, map_zoom, hodo_boundary, dpi, theme_overrides=None):
    """Build a hodograph figure from calculated parameters, shared by sounding and VAD data.

    Data with thermodynamics gets the stepwise CAPE/CIN panel and parcel values. Wind-only data gets
    the wind component panel and kinematic values.
    """
    # record process time
    st = time.time()

    #################################################################
    ### SET UP PLOT THEME ###
    #################################################################
    theme = _get_plot_theme(dark_mode=dark_mode, overrides=theme_overrides)
    layout = FULL_HODOGRAPH_LAYOUT
    has_thermo = "thermo" in snd_params
    has_map = map_zoom > 0


    #################################################################
    ### STORM-RELATIVE HODO LOGIC ###
    #################################################################
    # readouts always list ground-relative motions, while the hodograph can be drawn storm-relative
    hodo_params = snd_params
    hodo_title = 'HODOGRAPH'
    if sr_hodo:
        if ma.is_masked(snd_params["kinem"]["sm_u"]) or np.isnan(snd_params["kinem"]["sm_u"]):
            sr_hodo = False
            warnings.warn("This profile can not be plotted as storm relative because a storm-motion does not exist for this data.\n"
                          "This plot will feature a ground relative instead.", Warning)
        else:
            hodo_params = _make_storm_relative(snd_params)
            hodo_title = 'STORM RELATIVE HODOGRAPH'


    #################################################################
    ### BUILD FIGURE AND HODOGRAPH ###
    #################################################################
    print("    > BUILDING HODOGRAPH")
    fig = plt.figure(figsize=(16, 12), linewidth=1.5, edgecolor=theme["border"], facecolor=theme["background"], dpi=dpi)

    hodo = _setup_hodograph(fig.add_axes(layout["hodograph"]), theme, (-90, 90), (-90, 90), component_range=180)
    _plot_hodograph_velocity_labels(hodo, theme, fontsize=FONT_SIZES["small"])
    _plot_hodograph_height_markers(hodo, hodo_params, theme, marker_size=20, fontsize=FONT_SIZES["tiny"])
    _plot_hodograph_trace(hodo, hodo_params, theme, linewidth=5)

    _plot_hodograph_motion_labels(hodo, hodo_params, storm_motion, theme)
    if not sr_hodo:
        _plot_hodograph_sm_vector(hodo, hodo_params, theme)
    _plot_hodograph_motion_readout(hodo, snd_params, storm_motion, theme, x=0.015, y=0.985, fontsize=FONT_SIZES["note"])
    _plot_hodograph_srh(hodo, hodo_params, theme, layer="auto")

    if hodo_boundary is not None:
        _plot_hodograph_boundary(hodo, hodo_boundary, hodo_params, theme, show_mw_vector=True, show_lm_vector=True)


    #################################################################
    ### MAP INSET ###
    #################################################################
    if has_map:
        print("    > BUILDING MAP INSET")
        latlon = site_data["site_info"]["site-latlon"]
        map_ax = fig.add_axes(layout["map"], projection=ccrs.PlateCarree(), zorder=12)
        _setup_map_inset(map_ax, latlon, theme, map_zoom=map_zoom)

        if radar == "mosaic":
            _plot_map_mosaic_radar(map_ax, site_data, map_zoom, radar_time, theme)

        elif radar in ("single-site", "single"):
            _plot_map_single_site_radar(map_ax, site_data, radar_time, theme)

        _plot_map_profile_location(map_ax, latlon, theme)


    #################################################################
    ### LOW-LEVEL PANELS ###
    #################################################################
    print(f'    > BUILDING ACCESSORY PLOTS')
    _plot_streamwiseness_panel(fig.add_axes(layout["streamwiseness"]), snd_params, theme)
    _plot_vorticity_panel(fig.add_axes(layout["vorticity"]), snd_params, theme)
    _plot_sr_wind_panel(fig.add_axes(layout["sr_wind"]), snd_params, theme)
    plot_last_panel = _plot_cape_cin_panel if has_thermo else _plot_wind_component_panel
    plot_last_panel(fig.add_axes(layout["stability"]), snd_params, theme)


    #################################################################
    ### TEXT PLOTS ###
    #################################################################
    print(f'    > PLOTTING PARAMETERS')
    plot_parameters = _plot_hodograph_parameter_box if has_thermo else _plot_wind_parameter_box
    plot_parameters(fig.add_axes(layout["parameter_box"], zorder=3), snd_params, theme)


    #################################################################
    ### PLOT EXTRAS ###
    #################################################################
    print(f'    > FIGURE LOADING...')
    titles = dict(titles, top_title=f"{titles['top_title']} {hodo_title}")
    _plot_figure_titles(fig, titles, theme, left_x=0.2262, right_x=1.1815, top_y=0.92, sub_y=0.89, top_fontsize=FONT_SIZES["title_small"], sub_fontsize=FONT_SIZES["table_large"])
    # the watermark sits beside the map, or at the edge of the axis without one
    _plot_credits(hodo.ax, theme, x=0.298 if has_map else 0.0873, title_y=0.045, subtitle_y=0.017, style="compact", title_fontsize=FONT_SIZES["table_large"], subtitle_fontsize=FONT_SIZES["label"], max_x=0.99)
    _plot_logo(fig.add_axes(layout["logo" if has_map else "logo_no_map"], anchor="SW"))
    _add_figure_margin(fig, layout["margin"])

    print('    > COMPLETE --------')
    elapsed_time = time.time() - st
    print('    > RUNTIME:', time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))

    return plt



def __full_hodograph(clean_data, dark_mode, storm_motion, sr_hodo, modify_sfc,
                     radar, radar_time, map_zoom, hodo_boundary, dpi, max_levels=250, theme_overrides=None):

    print(f'    > PREPARING DATA')
    sounding_data = _prepare_sounding_data(clean_data, modify_sfc, max_levels)

    # calculate other sounding parameters using SounderPy Calc
    snd_params = sounding_params(sounding_data, storm_motion, include_all_parcels=True, parcel_depth=STABILITY_PANEL_DEPTH).calc()

    return _build_hodograph_figure(snd_params, sounding_data, clean_data["titles"], dark_mode, storm_motion, sr_hodo,
                                   radar, radar_time, map_zoom, hodo_boundary, dpi, theme_overrides)




















#########################################################################
######################## COMPOSITE SOUNDING #############################

def __composite_sounding(data_list, shade_between, cmap, colors, linestyles, alphas, linewidths, dark_mode,
                         show_params=True, storm_motion='right_moving', max_levels=250, dpi=100, ylim=None, theme_overrides=None,
                         style="full", sr_hodo=False):

    # record process time
    st = time.time()

    if style not in COMPOSITE_LAYOUTS:
        raise ValueError(f"'style' must be one of {', '.join(repr(name) for name in COMPOSITE_LAYOUTS)}. Got: {style!r}")
    if style == "hodograph" and ylim is not None:
        warnings.warn("'ylim' only applies to the skew-t, so it is not used with style='hodograph'.", Warning)
    show_skewt = style in ("full", "sounding")
    show_hodograph = style in ("full", "hodograph")
    if sr_hodo and not show_hodograph:
        warnings.warn("'sr_hodo' only applies to the hodograph, so it is not used with style='sounding'.", Warning)
        sr_hodo = False

    #################################################################
    ### SET UP PLOT THEME AND STYLES ###
    #################################################################
    theme = _get_plot_theme(dark_mode=dark_mode, overrides=theme_overrides)
    layout = COMPOSITE_LAYOUTS[style]
    n_profiles = len(data_list)
    colors = _get_profile_colors(colors, cmap, n_profiles, dark_mode)
    styles = _get_profile_styles(linestyles, "-", n_profiles)
    widths = _get_profile_styles(linewidths, 3, n_profiles)
    alphas = _get_profile_styles(alphas, 1, n_profiles)


    #################################################################
    ### SET UP THE DATA ###
    #################################################################
    print(f'    > PREPARING DATA')
    profiles = [_prepare_sounding_data(profile, None, max_levels) for profile in data_list]
    winds = [_get_profile_winds(profile) for profile in profiles]

    # a storm-relative hodograph needs each profile's storm motion, even without the table's parameters
    snd_params = None
    if show_params or sr_hodo:
        print(f'    > COMPUTING PARAMETERS')
        snd_params = []
        for i, profile in enumerate(profiles):
            try:
                snd_params.append(sounding_params(profile, storm_motion).calc())
            except Exception as exc:
                warnings.warn(f"Parameters could not be computed for profile {_get_profile_letter(i)}: {exc}", Warning)
                snd_params.append(None)

    if sr_hodo:
        missing = _get_missing_storm_motions(snd_params)
        if missing:
            sr_hodo = False
            warnings.warn(f"Profile{'s' if len(missing) > 1 else ''} {', '.join(missing)} {'have' if len(missing) > 1 else 'has'} no storm motion, "
                          "so the hodograph is ground relative instead.", Warning)
        else:
            winds = _make_winds_storm_relative(winds, snd_params)
    table_params = snd_params if show_params else None


    #################################################################
    ### BUILD FIGURE AND SKEW-T ###
    #################################################################
    fig = plt.figure(figsize=layout["figsize"], linewidth=10, edgecolor=theme["border"], facecolor=theme["background"], dpi=dpi)

    skew = None
    if show_skewt:
        print(f'    > BUILDING SKEW-T LOG-P DIAGRAM')
        # the temperature axis fits both the coldest and warmest surfaces
        surface_temperatures = [profile["T"][0] for profile in profiles]
        temperature_limits = (_get_skewt_temperature_bounds(min(surface_temperatures))[0], _get_skewt_temperature_bounds(max(surface_temperatures))[1])
        ylim = _validate_pressure_limits(ylim)
        if ylim is None:
            skew = _setup_skewt(fig, profiles[0]["T"], theme, rect=layout["skewt"], temperature_limits=temperature_limits)
        else:
            # a zoomed view fits its temperature range and pressure ticks to that part of the atmosphere
            skew = _setup_zoomed_skewt(fig, profiles, theme, ylim, layout["skewt"])
        _plot_skewt_accessories(skew, theme)
        _plot_composite_skewt_profiles(skew, profiles, colors, styles, widths, alphas, shade_between)


    #################################################################
    ### ADD HODOGRAPH ###
    #################################################################
    if show_hodograph:
        print("    > BUILDING HODOGRAPH")
        all_winds = {"u": np.concatenate([w["u"] for w in winds]), "v": np.concatenate([w["v"] for w in winds])}
        hodo_xlim, hodo_ylim = _get_hodograph_bounds(all_winds)
        hodo = _setup_hodograph(fig.add_axes(layout["hodograph"]), theme, hodo_xlim, hodo_ylim)
        _plot_hodograph_velocity_labels(hodo, theme)
        _plot_composite_hodograph_traces(hodo, winds, colors)
        _plot_composite_hodograph_key(hodo, theme)
        # the hodograph-only figure says this in its title
        if sr_hodo and style == "full":
            hodo.ax.text(0.02, 0.98, "STORM RELATIVE", transform=hodo.ax.transAxes, ha="left", va="top",
                         color=theme["text"], weight="bold", fontsize=FONT_SIZES["table"], zorder=8)


    #################################################################
    ### TEXT PLOTS ###
    #################################################################
    print(f'    > PLOTTING PROFILE TABLE')
    # the table and everything below it grow with the number of profiles
    table_height = layout["table_row_height"] * (n_profiles + 1.2)
    table_bottom = layout["table_top"] - table_height
    _plot_composite_table(fig.add_axes((layout["table_left"], table_bottom, layout["table_width"], table_height), zorder=3),
                          profiles, table_params, colors, theme, style)


    #################################################################
    ### PLOT EXTRAS ###
    #################################################################
    print(f'    > FIGURE LOADING...')
    left_x, right_x = layout["title_x"]
    _plot_figure_titles(fig, _get_composite_titles(profiles, style, sr_hodo), theme, left_x=left_x, right_x=right_x, sub_fontsize=layout["subtitle_fontsize"])
    if skew is not None:
        _plot_skewt_legend(skew, theme)
    credits_x, credits_max_x = layout["credits_x"]
    _plot_credits(fig, theme, x=credits_x, title_y=table_bottom - 0.025, subtitle_y=table_bottom - 0.045, max_x=credits_max_x)
    logo_x, logo_width = layout["logo"]
    _plot_logo(fig.add_axes((logo_x, table_bottom - 0.055, logo_width, 0.05), anchor="SE", zorder=3))
    margin_left, margin_right = layout["margin_x"]
    margin_bottom = table_bottom - 0.065
    _add_figure_margin(fig, (margin_left, margin_bottom, margin_right - margin_left, 1.0185 - margin_bottom))

    print('    > COMPLETE --------')
    elapsed_time = time.time() - st
    print('    > RUNTIME:', time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))

    return plt





#########################################################################
########################### VAD HODOGRAPH ###############################

def __vad_hodograph(vad_data, dark_mode, storm_motion, sr_hodo, radar='mosaic', radar_time='sounding', map_zoom=2, hodo_boundary=None, dpi=100, theme_overrides=None):

    print(f'    > PREPARING DATA')
    kinem, intrp = vad_params(vad_data, storm_motion).calc()

    site = vad_data["site_info"]
    valid = site["valid-time"]
    titles = {"top_title": site["source"],
              "left_title": f"VALID: {valid[1]}-{valid[2]}-{valid[0]} {valid[3]}Z",
              "right_title": f"{site['site-id']} - {site['site-name']} | {site['site-latlon'][0]}, {site['site-latlon'][1]}"}

    return _build_hodograph_figure({"kinem": kinem, "intrp": intrp}, vad_data, titles, dark_mode, storm_motion, sr_hodo,
                                   radar, radar_time, map_zoom, hodo_boundary, dpi, theme_overrides)

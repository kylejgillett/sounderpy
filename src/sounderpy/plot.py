# METPY OPERATIONS 
import metpy.calc as mpcalc
from metpy.units import units
from metpy.plots import SkewT, Hodograph, USCOUNTIES

#MATPLOTLIB OPERATIONS
import matplotlib.lines    as mlines
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import matplotlib.colors as mcolors

# NUMPY OPERATIONS
import numpy as np
import numpy.ma as ma
import copy

# STANDARD PYTHON
import warnings
import time
from datetime import datetime

# ECAPE
from ecape_parcel.calc import calc_ecape_parcel, density_temperature

# FOR MAPPING
from urllib.request import urlopen
from PIL import Image
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt

# DATA RETRIEVAL
from xarray.backends import NetCDF4DataStore
from xarray import open_dataset
from siphon.catalog import TDSCatalog
from netCDF4 import Dataset
from matplotlib.colors import LinearSegmentedColormap

# SOUNDERPY MODULES
from .calc import sounding_params, vad_params
from .utils import modify_surface, find_nearest, mag, mag_round




"""
    SOUNDERPY SPYPLOT FUNCTIONS  

    Purpose of module: 

    House tools for plotting sounding, hodograph, and composite sounding
    figures from data provided by the user. Functions here are referenced
    by sounderpy.py.     

    (C) KYLE J GILLETT, UNIVERSITY OF NORTH DAKOTA, 2024
"""



#########################################################################
########################## FULL SOUNDING ################################

def __full_sounding(clean_data, color_blind, dark_mode, storm_motion, special_parcels=None, 
                    show_radar=True, radar_time='sounding', map_zoom=2, modify_sfc=None,
                    show_theta=False):
    
    
    # record process time 
    st = time.time()  


    
    
    # get radar mosaic data for map inset
    def get_radar_mosaic(clean_data, map_zoom, radar_time):
        '''
        GET RADAR MOSAIC DATA FOR SOUNDING PLOT
        
        - gets data from NCEI for 'sounding time' or now
        - data only goes back 1 month from current time
        '''
    
        # search the NCEI database for the file closest to 
        # the requested radar date/time
        def find_file(target_time, file_names):
            closest_index = None
            closest_difference = float('inf')

            for idx, file_name in enumerate(file_names):
                # Extract the final component (ex: "1220") from the file name
                components = file_name.split('_')
                if len(components) > 1:
                    # Remove file extension if present
                    last_component = components[-1].split('.')[0]
                    try:
                        file_time = int(last_component)
                        # Convert the target_time to an integer for comparison
                        target_time_int = int(target_time)
                        difference = abs(file_time - target_time_int)
                        if difference < closest_difference:
                            closest_difference = difference
                            closest_index = idx
                    except ValueError:
                        pass  # Skip if the last component is not a valid integer

            return closest_index


        # PREFORM DATETIME OPERATIONS------------------------------------------------------------------------
        # get sounding data from from clean_data
        data_dt = datetime(int(clean_data['site_info']['valid-time'][0]), int(clean_data['site_info']['valid-time'][1]), 
             int(clean_data['site_info']['valid-time'][2]), int(clean_data['site_info']['valid-time'][3][0:2]))

        # get current time
        curr_dt = datetime.utcnow()

        # if data time is after current time, default to now
        if data_dt > curr_dt:
            time = curr_dt.strftime('%H%M')
            datestr = curr_dt.strftime('%Y%m%d')
        
        # if radar time is given as now, time is now
        elif radar_time == 'now':
            time = curr_dt.strftime('%H%M')
            datestr = curr_dt.strftime('%Y%m%d')

        # if radar time is given as sounding, use the valid-time info from clean_data
        elif radar_time == 'sounding':
            if len(clean_data['site_info']['valid-time'][3]) == 2:
                time = f"{clean_data['site_info']['valid-time'][3]}00"
            else:
                time = clean_data['site_info']['valid-time'][3]
            datestr = data_dt.strftime('%Y%m%d')


        # PULL SOUNDING DATA LAT/LON ----------------------------------------------------------------------
        data_lat = clean_data['site_info']['site-latlon'][0]
        data_lon = clean_data['site_info']['site-latlon'][1]


        # GET RADAR DATA ----------------------------------------------------------------------------------
        composite_url = 'https://thredds.ucar.edu/thredds/catalog/nexrad/composite/gini/dhr/1km/'+datestr+'/catalog.xml'
        try:
            # try accessing the NCEI TDS database
            composite_catalog = TDSCatalog(composite_url).datasets
            data_found = True
        except:
            # no data is available, return no data
            print('- no radar data available -')
            data_found = False
            pass
        
        if data_found == True:
            # if data is foudn, parse it for plotting.
            filename = composite_catalog[find_file(time, composite_catalog)].subset()
            composite_query = filename.query()
            composite_query.lonlat_box(north = data_lat + (map_zoom + 1),
                             south = data_lat - (map_zoom + 1),
                             east  = data_lon + (map_zoom + 1),
                             west  = data_lon - (map_zoom + 1))
            composite_query.add_lonlat(value=True)
            composite_query.accept('netcdf4')
            composite_query.variables('Reflectivity')
            radar_data = filename.get_data(composite_query)
            radar_data = open_dataset(NetCDF4DataStore(radar_data))

            return radar_data
        
    # define display mode colors and alphas    
    if dark_mode == True:
        gen_txt_clr = 'white'
        bckgrnd_clr = 'black'
        brdr_clr    = 'white'
        barb_clr    = 'white'
        shade_alpha = 0.06
        skw_ln_clr = 'white'
        marker_clr = 'white'
    else: 
        gen_txt_clr = 'black'
        bckgrnd_clr = 'white'
        brdr_clr    = 'black'
        barb_clr    = 'black'
        shade_alpha = 0.02
        skw_ln_clr = 'black'
        marker_clr = 'black'


    # define colorblindness colors
    if color_blind == True:
        td_color = 'cornflowerblue'
    else:
        td_color = 'green'
        
    

    #################################################################
    ### SET UP THE DATA ###
    #################################################################
    # SFC CORRECTION
    if str(type(modify_sfc)) == "<class 'dict'>":
        sounding_data = modify_surface(clean_data, modify_sfc)
        modify_sfc_txt = True
    else:
        sounding_data = clean_data
        modify_sfc_txt = False
    
    # declare easy variable names for reuse from `clean_data` 
    T  = sounding_data['T']
    Td = sounding_data['Td']
    p  = sounding_data['p']
    z  = sounding_data['z']
    u  = sounding_data['u']
    v  = sounding_data['v']
    wd = mpcalc.wind_direction(u, v)
    ws = mpcalc.wind_speed(u, v) 
    
    # calculate other sounding parameters using SounderPy Calc
    general, thermo, kinem, intrp = sounding_params(sounding_data, storm_motion, include_all_parcels=not show_theta).calc()
    #################################################################
    
    
    
    #################################################################
    ### DECLARE PLOT TITLES FROM CLEAN_DATA ###
    #################################################################
    
    top_title   = clean_data['titles']['top_title']
    left_title  = clean_data['titles']['left_title']
    right_title = clean_data['titles']['right_title']
    ################################################################    
        
    
    #########################################################################
    ################################ SKEW-T ################################# 
    #########################################################################
    
    
    
    #################################################################
    ### CREATE FIGURE ###
    #################################################################
    #################################################################
    # create master figure
    fig = plt.figure(figsize=(22,13), linewidth=10, edgecolor=brdr_clr)         
    # create skew-t obj and axes params
    skew = SkewT(fig, rotation=47, rect=(0.1124, 0.1005, 0.60, 0.85))  
    skew.ax.set_box_aspect(0.87)
    skew.ax.zorder = 5
    # Define axis bounds, temp axis based on data
    skew.ax.set_adjustable('datalim')
    skew.ax.set_ylim(1050, 100)    
    if T[0].m <= -10:
        bounds = [-42, 42]
        if T[0].m <= -20:
            bounds = [-52, 32]
    elif T[0].m >= 20:
        bounds = [-22, 62]
    else: 
        bounds = [-32, 52]
    skew.ax.set_xlim(bounds[0], bounds[1])   
    # Define and customize axis labels 
    plt.xlabel("  ", fontsize=12)
    plt.ylabel("  ", fontsize=12) 
    plt.xticks(fontsize=13, weight='bold')  
    plt.yticks(fontsize=13, ha='left', weight='bold')
    plt.tick_params(axis="x", direction="in", pad=-12, colors=gen_txt_clr)
    plt.tick_params(axis="y", direction="in", pad=-7, colors=gen_txt_clr)
    skew.ax.set_yticks([1000, 900, 800, 700, 600, 500, 400, 300, 200])
    skew.ax.set_yticklabels([1000, 900, 800, 700, 600, 500, 400, 300, 200], color=gen_txt_clr)
    # add skew-t axis spines
    skew.ax.spines["top"].set_color(brdr_clr)
    skew.ax.spines["left"].set_color(brdr_clr)
    skew.ax.spines["right"].set_color(brdr_clr)
    skew.ax.spines["bottom"].set_color(brdr_clr)
    skew.ax.spines["bottom"].set_color(brdr_clr)
    # Define background colors
    fig.set_facecolor(bckgrnd_clr)         
    skew.ax.set_facecolor(bckgrnd_clr)    
    # Add shaded isotherms
    x1 = np.linspace(-100, 40, 8)                                                          
    x2 = np.linspace(-90, 50, 8)                                                         
    y = [1200, 50]                                                                      
    for i in range(0, 8):              
        skew.shade_area(y=y, x1=x1[i], x2=x2[i], color='gray', alpha=shade_alpha, zorder=1)   
    #################################################################
    
    
    
    #################################################################
    ### PLOT SKEW T LINES ###
    #################################################################
    # Plot relevent Skew-T lines
    # 0C and -20C highlights
    skew.ax.axvline(0 * units.degC, linestyle='--', color='blue', alpha=0.3)
    skew.ax.axvline(-20 * units.degC, linestyle='--', color='blue', alpha=0.3)
    # dry adiabats
    skew.plot_dry_adiabats(color='cornflowerblue', linewidth=0.2, alpha=0.7) 
    # moist adiabats
    skew.plot_moist_adiabats(color='cornflowerblue', linewidth=0.2, alpha=0.7)
    # mixing ratio lines
    skew.plot_mixing_lines(color=skw_ln_clr, linewidth=0.2, alpha=0.7)
    # wet bulb temp (Tw)
    twline = skew.plot(p, general['wet_bulb'], '#3d8aff', linewidth=2, label='WETBULB TEMP', alpha=0.6)
    # virtual temp (Tv)
    tvline = skew.plot(p, general['virt_temp'], 'darkred', linestyle=(0, (1, 1)), linewidth=4, label='VIRTUAL TEMP', alpha=0.7)  
    # dewpoint temp (Td)
    tdline = skew.plot(p, Td, td_color, linewidth=5, label='DEWPOINT')
    # temp (T)
    tline1 = skew.plot(p, T, 'red', linewidth=5, label='TEMPERATURE')  

    
    
    #################################################################
    ### PARCELS LOGIC ###
    #################################################################
    '''
    SounderPy 'Special Parcels Logic' using Amelia Urquhart’s ecape_parcels library: https://github.com/a-urq/ecape-parcel-py
    
    Docs: https://kylejgillett.github.io/sounderpy/plottingdata.html#parcel-logic
    '''
    parcel_details = {

            'term1': {'sb': 'surface_based', 
                      'ml': 'mixed_layer', 
                      'mu': 'most_unstable',
                     },
            'term2': {'ps': [True,  'pseudoadiabatic'],
                      'ia': [False, 'irrev-adiabatic']
                     },
            'term3': {'ecape': True,
                       'cape': False
                     },
            'sm' :   {'rm': 'right_moving',
                      'lm': 'left_moving',
                      'mw': 'mean_wind',
                     }
    }
    
    valid_parcel_codes = ['mu_ps_cape', 'mu_ia_cape', 'mu_ps_ecape', 'mu_ia_ecape',
                          'ml_ps_cape', 'ml_ia_cape', 'ml_ps_ecape', 'ml_ia_ecape',
                          'sb_ps_cape', 'sb_ia_cape', 'sb_ps_ecape', 'sb_ia_ecape']
    
    # CREATE AND ADD SPECIAL PARCELS IF THE SPECIAL_PARCELS
    # VARIABLE CONTAINS PARCEL NAMES
    if str(type(special_parcels)) == "<class 'list'>":

        print('COMPUTING & PLOTTING SPECIAL PARCELS -- this may take a moment')
            
        # combine the two parcel lists to a single list
        all_parcels = sum(special_parcels, [])

        # create empty dict of parcel data
        parcel_dict = {}

        for parcel in all_parcels:
            
            if parcel not in valid_parcel_codes:
                    sys.exit(f"Invalid parcel code provided in 'special_parcels' kwarg. Please make sure parcel codes are one of these:\n{valid_parcel_codes}\n" +
                            "Valid parcel codes include the parcel type ('mu', 'sb', 'ml'), adiabatic scheme ('ps', 'ia'), and cape type ('cape', 'ecape').\n" +
                            "See 'https://kylejgillett.github.io/sounderpy/plottingdata.html#parcel-logic' for more info.")

            # set vars to None 
            trace = None 
            trace_Trho = [None]

            # build parcel trace obj
            trace = calc_ecape_parcel(p, z, T, Td, u, v, 
                              False, 
                              cape_type=parcel_details['term1'][parcel.split('_')[0]],
                              pseudoadiabatic_switch=parcel_details['term2'][parcel.split('_')[1]][0], 
                              entrainment_switch=parcel_details['term3'][parcel.split('_')[2]], 
                              storm_motion_type='user_defined', 
                              storm_motion_u=kinem['sm_u']*units.kts, storm_motion_v=kinem['sm_v']*units.kts)

            # calc parcel density temperature
            if (trace[0][0] != None):
                 trace_Trho = density_temperature(trace[2], trace[3], trace[4])
                 # add Trho to the end of the trace tuple
                 trace = trace + (trace_Trho, ) 


            # add parcel data to the parcel_dict
            parcel_dict[parcel] = trace


        # loop through the parcels again to plot them  
        for parcel in parcel_dict.keys():

            if len(parcel_dict[parcel]) == 6:

                parcel_name = str.upper(f"{parcel.split('_')[0]} {parcel_details['term2'][parcel.split('_')[1]][1]} {parcel.split('_')[2]}")
                # if parcel exists, plot it
                if parcel in special_parcels[0]:
                            skew.plot(parcel_dict[parcel][0], parcel_dict[parcel][-1], color='red', linestyle='--',  
                                             linewidth=2, alpha=1, label=parcel_name)
                # if parcel exists, plot it
                elif parcel in special_parcels[1]:
                            skew.plot(parcel_dict[parcel][0], parcel_dict[parcel][-1], color='#808080', linestyle='--',  
                                             linewidth=2, alpha=0.5, label=parcel_name)
    
    
    
    
    # ADD BASIC PARCEL TRACES IF NO SPECIAL PARCFELS ARE 
    # PROVIDED BY THE USER.
    elif special_parcels in [None, 'simple']:
        if (special_parcels != 'simple'):
            if thermo['mu_ecape'] > 0:
                trace = calc_ecape_parcel(p, z, T, Td, u, v, 
                                  False, 
                                  cape_type='most_unstable',
                                  pseudoadiabatic_switch='true', 
                                  entrainment_switch=True, 
                                  storm_motion_type='user_defined', 
                                  storm_motion_u=kinem['sm_u']*units.kts, storm_motion_v=kinem['sm_v']*units.kts)

                trace_Trho = density_temperature(trace[2], trace[3], trace[4])

                muecapeline = skew.plot(trace[0], trace_Trho, 
                                        linestyle='--', linewidth=3, alpha=1, color='red',
                                        label='MUECAPE PARCEL')
            
#                 print(trace_Trho[0:len(T)].m-273.15)
#                 return trace[0]
                
#                 skew.ax.fill_betweenx(p, T, trace_Trho[0:len(T)].m-273.15,
#                                       color='red', alpha=0.1)
                
                #skew.shade_cape(p, T, np.linspace(trace_Trho[0], trace_Trho[-1], len(T)), color='red', alpha=0.4) 
                #skew.shade_cape(p, T, mu_parcel_path, color='orange') 
        
        # if CAPE for SB, MU, ML parcels is >0, plotm them
        if thermo['sbcape'] > 0:
            sbparcelline = skew.plot(thermo['sbP_trace'], thermo['sbT_trace'], 
                                     linestyle='--', linewidth=2, alpha=0.5, color='#808080', 
                                     label='SBCAPE PARCEL')    
        if thermo['mlcape'] > 0:
            mlparcelline = skew.plot(thermo['mlP_trace'], thermo['mlT_trace'],
                                     linestyle='--', linewidth=2, alpha=0.5, color='#808080',
                                     label='MLCAPE PARCEL') 
        if thermo['mucape'] > 0:
            muparcelline = skew.plot(thermo['muP_trace'], thermo['muT_trace'],
                                     linestyle='--', linewidth=2, alpha=0.5, color='#808080',
                                     label='MUCAPE PARCEL')

            
    # ADD DOWNDRAFT PARCEL TRACE
    skew.plot(thermo['dparcel_p'], thermo['dparcel_T'], linestyle='--',linewidth=0.7, color='purple', 
              alpha=0.8, label='DWNDRFT PARCEL')
    
    
    #################################################################
    ### PLOT SKEW T ANNOTATIONS ###
    #################################################################  
    #plot left side hgt level markers
    hgt_lvls =[]
    for key in intrp['hgt_lvls'].keys():
        hgt_lvls.append(intrp['hgt_lvls'][key])
    hgt_lvls.pop(0) 

    for key in hgt_lvls[1::4]:
        # for a hgt_lvl key, find the pressure at that level's index, plot the hgt at the pressure
        trans, _, _ = skew.ax.get_yaxis_text1_transform(0)
        skew.ax.text(0.048, intrp['pINTRP'][key], f"{int(intrp['zINTRP'][key]/1000)}km", 
                     fontsize=11, transform=trans, alpha=0.6, weight='bold', color=gen_txt_clr) 
        
    # add a '-sfc-' marker with an elevation of the sfc in meters
    trans, _, _ = skew.ax.get_yaxis_text1_transform(0)
    sfc = mpcalc.height_to_pressure_std(general['elevation']*units.m)
    skew.ax.text(0.048, p[0], f"-SFC ({mag(general['elevation'])}m) -", 
                 fontsize=11, transform=trans, alpha=0.6, weight='bold', color=gen_txt_clr)
    
    
    # SFC TEMPERATURE AND DEWPOINT ANNOTATIONS---------------------------------------------
    # plot sfc T value in degF
    T_degF = np.round(T.to(units.degF), 1)
    T_degF_label = '{}°F'.format(int(T_degF[0].magnitude))                             
    plt.annotate(T_degF_label, (T[0], p[0]), textcoords="offset points", xytext=(16,-15),
                     fontsize=12, color='red', weight='bold', alpha=0.7, ha='center') 
    
    # plot sfc Td value in DegF
    Td_degF = np.round(Td.to(units.degF), 1) 
    Td_degF_label = '{}°F'.format(int(Td_degF[0].magnitude))                             
    plt.annotate(Td_degF_label,(Td[0], p[0]),textcoords="offset points",xytext=(-16,-15), 
                     fontsize=12, color=td_color, weight='bold', alpha=0.7, ha='center') 
    
    # if the sfc was modified by the user, plot a modified sfc flag
    if modify_sfc_txt == True: 
        plt.annotate("USER MODIFIED SFC →", (Td[0], p[0]),textcoords="offset points",xytext=(-50,0), 
                     fontsize=10, color=gen_txt_clr, weight='bold', alpha=0.5, ha='right')
        
        
    # PARCEL HEIGHT ANNOTATIONS------------------------------------------------------------- 
    # plot the SBLCL and draw a line between the sfc and the lcl
    plt.text((0.82), (thermo['sb_lcl_p']), "←SBLCL", weight='bold',color='gray',
             alpha=0.9, fontsize=15, transform=skew.ax.get_yaxis_transform(), clip_on=True)
    lcl_line = plt.Line2D([0.86, 0.86], (p[0], thermo['sb_lcl_p']), 
                          color='gray', linestyle=':', alpha=1, transform=skew.ax.get_yaxis_transform())
    skew.ax.add_artist(lcl_line)
    
    # if a mulfc doesn't exist,
    # plot the sblfc, and sbel and draw a line between the lcl and lfc and the lcl to the el
    if ma.is_masked(thermo['mu_lfc_z']) == True:
        plt.text((0.82), (thermo['sb_lfc_p']), "←SBLFC", weight='bold',color='gray',
                 alpha=0.9, fontsize=15, transform=skew.ax.get_yaxis_transform(), clip_on=True)
        plt.text((0.82), (thermo['sb_el_p']), "←SBEL", weight='bold',color='gray', 
                 alpha=0.9, fontsize=15, transform=skew.ax.get_yaxis_transform(), clip_on=True)
        el_line = plt.Line2D([0.86, 0.86], (p[0], thermo['sb_el_p']), 
                             color='gray', linestyle=':', alpha=0.3, transform=skew.ax.get_yaxis_transform())
        skew.ax.add_artist(el_line)
        lfc_line = plt.Line2D([0.86, 0.86], (p[0], thermo['sb_lfc_p']), 
                              color='gray', linestyle=':', alpha=0.7, transform=skew.ax.get_yaxis_transform())
        skew.ax.add_artist(lfc_line)
    
    else: 
        # if not, plot mulfc and muel and draw a line between the lcl and lfc and the lcl to the el
        plt.text((0.82), (thermo['mu_lfc_p']), "←MULFC", weight='bold',color='gray',
                 alpha=0.9, fontsize=15, transform=skew.ax.get_yaxis_transform(), clip_on=True)
        plt.text((0.82), (thermo['mu_el_p']), "←MUEL", weight='bold',color='gray',
                 alpha=0.9, fontsize=15, transform=skew.ax.get_yaxis_transform(), clip_on=True)
        el_line = plt.Line2D([0.86, 0.86], (p[0], thermo['mu_el_p']), 
                             color='gray', linestyle=':', alpha=0.3, transform=skew.ax.get_yaxis_transform())
        skew.ax.add_artist(el_line)
        lfc_line = plt.Line2D([0.86, 0.86], (p[0], thermo['mu_lfc_p']), 
                              color='gray', linestyle=':', alpha=0.7, transform=skew.ax.get_yaxis_transform())
        skew.ax.add_artist(lfc_line)
        
        # plot mumpl, but if its above 110 hPa, just plot whatever the value is at 110 hPa
        if thermo['mu_mpl_p'] < 110:
            plt.text((0.82), (110), f"↑MUMPL: {mag(thermo['mu_mpl_p'])}hPa", weight='bold',color='gray',
                     alpha=0.9, fontsize=15, transform=skew.ax.get_yaxis_transform(), clip_on=True)
        else: 
            plt.text((0.82), (thermo['mu_mpl_p']), "←MUMPL", weight='bold',color='gray',
                     alpha=0.9, fontsize=15, transform=skew.ax.get_yaxis_transform(), clip_on=True)
        
        
        
        
    # T & WB FRREZING POINT ANNOTATION--------------------------------------------------------------
    if ma.is_masked(general['frz_pt_z']) == False:
        if general['frz_pt_z'] >= 50*units.m:
            plt.text((0.765), (general['frz_pt_p']), "←FRZ", weight='bold',color='cornflowerblue',  
                     alpha=0.6, fontsize=12, transform=skew.ax.get_yaxis_transform(), clip_on=True)
    
    if ma.is_masked(general['wb_frz_pt_z']) == False:
        if general['wb_frz_pt_z'] >= 50*units.m:
            plt.text((0.765), (general['wb_frz_pt_p']), "←WB0", weight='bold',color='cornflowerblue',  
                     alpha=0.6, fontsize=12, transform=skew.ax.get_yaxis_transform(), clip_on=True)

            
    # PBL TOP POINT ANNOTATION---------------------------------------------------------------                   
    plt.text((0.78), (thermo['pbl_top']), "←PBL", weight='bold',color='gray', 
             alpha=0.9, fontsize=10, transform=skew.ax.get_yaxis_transform(), clip_on=True)
    pbl_line = plt.Line2D([0.80, 0.80], (p[0], thermo['pbl_top']), color='gray', 
                          linestyle=':', alpha=0.7, transform=skew.ax.get_yaxis_transform())
    skew.ax.add_artist(pbl_line)

    
    # 0-3km & 0-6km CAPE ANNOTATIONS--------------------------------------------------------
    if thermo['mu3cape'] > 10:
        idx = find_nearest(thermo['muZ_trace'], 3000)
        cape03_label = " ←{}J/kg".format(mag(thermo['mu3cape']))     
        plt.annotate(cape03_label,  ((thermo['muT_trace'][idx]+7), intrp['pINTRP'][intrp['hgt_lvls']['h3']]),  textcoords="offset points",  xytext=(20, 0), 
                     color='red', alpha=0.6, fontsize=13.5, ha='right')  
            
            
    if thermo['mu6cape'] > thermo['mu3cape']:
        idx = find_nearest(thermo['muZ_trace'], 6000)
        cape06_label = " ←{}J/kg".format(mag(thermo['mu6cape']))                                            
        plt.annotate(cape06_label, ((thermo['muT_trace'][idx]+7), intrp['pINTRP'][intrp['hgt_lvls']['h6']]), textcoords="offset points",  xytext=(10, 0), 
                         color='red', alpha=0.6, fontsize=13.5, ha='right') 
            
            
    if thermo['mucape'] > thermo['mu6cape']:
            cape_label = "←{}J/kg".format(mag(thermo['mucape']))                                            
            plt.annotate(cape_label,((thermo['mu_el_T']), thermo['mu_el_p']), textcoords="offset points",  xytext=(5, 0), 
                         color='red', alpha=0.6, fontsize=13.5, ha='left') 
            
            
            
    # MAX LAPSE RATE ANNOTATION---------------------------------
    x_start, x_end = 0.15, 0.17
    x_mid = (x_start + x_end)/2
    Lapse_line = plt.Line2D([x_mid, x_mid], (thermo['lr_max'][1], thermo['lr_max'][2]), color='firebrick', linewidth=3, alpha=0.3, transform=skew.ax.get_yaxis_transform())
    plt.text((x_start-0.005), (thermo['lr_max'][2]-8), f"{np.round(thermo['lr_max'][0],1)}", 
             color='firebrick', weight='bold', fontsize=13, alpha=0.3, transform=skew.ax.get_yaxis_transform())
    skew.ax.add_artist(Lapse_line)

    # EFFECTIVE INFLOW LAYER ANNOTATION--------------------------
    x_start, x_end = 0.18, 0.20
    x_mid = (x_start + x_end)/2
    plt.text((x_start+0.01), (kinem['eil'][1]-8), "EIL", weight='bold',color='lightblue', alpha=0.4, ha='center', fontsize=13, transform=skew.ax.get_yaxis_transform())
    EIL_line = plt.Line2D([x_mid, x_mid], (kinem['eil'][0], kinem['eil'][1]), color='lightblue', linewidth=3, alpha=0.4, transform=skew.ax.get_yaxis_transform())
    skew.ax.add_artist(EIL_line)

    # DGZ ANNOTATION-------------------------------
    if T[0].m < 5:
        if thermo['dgz'][1] < p[0].m:
            x_start, x_end = 0.12, 0.14
            x_mid = (x_start + x_end)/2
            plt.text((x_start+0.01), (thermo['dgz'][1]-8), "DGZ", weight='bold',color='blue', alpha=0.4, ha='center', fontsize=13, transform=skew.ax.get_yaxis_transform())
            dgz_line = plt.Line2D([x_mid, x_mid], (thermo['dgz'][0], thermo['dgz'][1]), color='blue', linewidth=3, alpha=0.4, transform=skew.ax.get_yaxis_transform())
            skew.ax.add_artist(dgz_line)
      
    # HGZ ANNOTATION-------------------------------
    else:
        if thermo['hgz'][1] < p[0].m:
            x_start, x_end = 0.12, 0.14
            x_mid = (x_start + x_end)/2
            plt.text((x_start+0.01), (thermo['hgz'][1]-8), "HGZ", weight='bold',color='green', alpha=0.4, ha='center', fontsize=13, transform=skew.ax.get_yaxis_transform())
            hgz_line = plt.Line2D([x_mid, x_mid], (thermo['hgz'][0], thermo['hgz'][1]), color='green', alpha=0.4, linewidth=3, transform=skew.ax.get_yaxis_transform())
            skew.ax.add_artist(hgz_line)
            
            
    # RELATIVE HUMIDITY ANNOTATIONS --------------------------      
    hgts = [10, 30, 60, 90, 120]

    try:
        for hgt_idx in hgts:
            label = "{}% →".format(int(intrp['rhINTRP'][hgt_idx]))
            plt.annotate(label, (mpcalc.dewpoint_from_relative_humidity(intrp['tINTRP'][hgt_idx]*units.degC, 
                                                                    intrp['rhINTRP'][hgt_idx]/100),
                         intrp['pINTRP'][hgt_idx]*units.hPa), textcoords="offset points", xytext=(-40, 0), 
                         color='green', alpha=0.6, fontsize=13.5, ha='center', clip_on=True) 
    except:
        pass
    #################################################################
    
    
    
    
    #################################################################
    ### PLOT SKEW T WIND BARBS ###
    #################################################################
    # Arrange wind barbs for best fit & resample wind barbs for best fit
    interval = np.logspace(2.113, 3, 30) *units.hPa
    idx = mpcalc.resample_nn_1d(p, interval) 

    # create blank barbs for small dot at the start of each actual barb
    blank_len = len(u[idx])     
    blank = np.zeros(blank_len)  
    skew.plot_barbs(pressure=p[idx], u=blank, v=blank, xloc=0.955, fill_empty=True, color=barb_clr,
                    sizes=dict(emptybarb=0.075, width=0.18, height=0.4))

    # plot actual wind barbs
    skew.plot_barbs(pressure=p[idx], u=u[idx], v=v[idx], xloc=0.955, fill_empty=True, color=barb_clr,
                    sizes=dict(emptybarb=0.075, width=0.18, height=0.4), length=8)

    # Draw line underneath wind barbs
    line = mlines.Line2D([0.955, 0.955], [0.01,0.95],color=barb_clr,linewidth=0.5,
                         transform=skew.ax.transAxes,clip_on=False,zorder=1)
    skew.ax.add_line(line)  
    ################################################################
    
    
    #########################################################################
    ############################## HODOGRAPH ################################
    #########################################################################
    
    
    #################################################################
    ### DEFINE HODOGRAPH BOUNDS ###
    #################################################################
    # restructure u and v, p, ws and z data arrays based on corrected u and v arrays and hodo_layer depth
    if (z.max().m - z[0].m) > 9001: 
        hodo_hgt = 9000*units.m
        p_hodo, u_hodo, v_hodo, z_hodo = mpcalc.get_layer(p, u, v, z, depth=hodo_hgt)
    else:
        p_hodo = p
        u_hodo = u
        v_hodo = v
        z_hodo = z
    # determine max height of wind data to plot on hodograph in km (if hodo_layer = 9, 0-9km u and v are plotted)
    # remove nan values from base wind u and v component arrays to find min & max values.
    u_clean = u.magnitude[np.logical_not(np.isnan(u.magnitude))]
    v_clean = v.magnitude[np.logical_not(np.isnan(v.magnitude))]
    # restructure u and v, p, ws and z data arrays based on corrected u and v arrays and hodo_layer depth
    # define x and y min/max values from 'cleaned' and restructured u and v arrays
    x_min = u_hodo.min().m
    y_min = v_hodo.min().m
    x_max = u_hodo.max().m
    y_max = v_hodo.max().m
    
    # if statements to determine approprate x axis and y axis limits (to change dynamically with the data)
    if y_max >= 0:
        y_Maxlimit = (y_max + 30)
    if y_max < 0:
        y_Maxlimit = (y_max + 30)

    if x_max >= 0:
        x_Maxlimit = (x_max + 30)
    if x_max < 0:
        x_Maxlimit = (x_max + 30)

    if y_min >= 0:
        y_Minlimit = (y_min - 45)
    if y_min < 0:
        y_Minlimit = (y_min - 45)

    if x_min >= 0:
        x_Minlimit = (x_min - 45)
    if x_min < 0:
        x_Minlimit = (x_min - 45)
    #################################################################
        
        
        
    #################################################################
    ### CREATE HODOGRAPH OBJECT ###
    #################################################################
    # create hodograph object
    hod_ax = plt.axes((0.53, 0.1003, 0.85, 0.85))
    h = Hodograph(hod_ax, component_range=160.)
    # dynamically set hodograph axis bounds
    try:
        h.ax.set_xlim(x_Minlimit, x_Maxlimit)                                  
        h.ax.set_ylim(y_Minlimit, y_Maxlimit)                             
    except:
        h.ax.set_xlim(-65,65)
        h.ax.set_ylim(-65,65)
        pass
    # add hodograph grid lines
    h.add_grid(increment=20, color=gen_txt_clr, linestyle='-', linewidth=1.5, alpha=0.2) 
    h.add_grid(increment=10, color=gen_txt_clr, linewidth=1, linestyle='--', alpha=0.2) 
    # add background color
    h.ax.set_facecolor(bckgrnd_clr)
    # add spines 
    h.ax.spines["top"].set_color(brdr_clr)
    h.ax.spines["left"].set_color(brdr_clr)
    h.ax.spines["right"].set_color(brdr_clr)
    h.ax.spines["bottom"].set_color(brdr_clr)
    h.ax.spines["bottom"].set_color(brdr_clr)
    # define box aspect
    h.ax.set_box_aspect(1)
    # define tick label params (remove them)
    h.ax.set_yticklabels([])
    h.ax.set_xticklabels([])
    h.ax.set_xticks([])
    h.ax.set_yticks([])
    h.ax.set_xlabel(' ')
    h.ax.set_ylabel(' ')
    #################################################################
    
    
    
    #################################################################
    ### PLOT VELOCITY AND HEIGHT MARKERS ###
    #################################################################
    # add small velocity markers on the rings 
    plt.xticks(np.arange(0,0,1))
    plt.yticks(np.arange(0,0,1))
    for i in range(10,130,20):
        h.ax.annotate(str(i),(i,0),xytext=(0,2),textcoords='offset pixels',
                      clip_on=True,fontsize=12,weight='bold',alpha=0.2,zorder=0, color=gen_txt_clr)
    for i in range(10,130,20):
        h.ax.annotate(str(i),(0,i),xytext=(0,2),textcoords='offset pixels',
                      clip_on=True,fontsize=12,weight='bold',alpha=0.2,zorder=0, color=gen_txt_clr)
    for i in range(10,130,20):
        h.ax.annotate(str(i),(-i,0),xytext=(0,2),textcoords='offset pixels',
                      clip_on=True,fontsize=12,weight='bold',alpha=0.2,zorder=0, color=gen_txt_clr)
    for i in range(10,130,20):
        h.ax.annotate(str(i),(0,-i),xytext=(0,2),textcoords='offset pixels',
                      clip_on=True,fontsize=12,weight='bold',alpha=0.2,zorder=0, color=gen_txt_clr)

    
    # add 0.5km marker
    h.plot(intrp['uINTRP'][intrp['hgt_lvls']['h05']],intrp['vINTRP'][intrp['hgt_lvls']['h05']],marker='.', markeredgecolor='black',
           color='white', alpha=1, markersize=30, clip_on=True, zorder=5)
    h.ax.annotate(str('.5'),(intrp['uINTRP'][intrp['hgt_lvls']['h05']],intrp['vINTRP'][intrp['hgt_lvls']['h05']]),
                  weight='bold', fontsize=11, color='black',xytext=(0.02,-5),
                  textcoords='offset pixels',horizontalalignment='center',clip_on=True, zorder=6) 

    # add hgt markers 
    hgt_lvls = [] 
    for key in intrp['hgt_lvls'].keys():
        hgt_lvls.append(intrp['hgt_lvls'][key])
    hgt_lvls.pop(0) 

    for lvl in hgt_lvls[1::2]:
        if lvl < 130:
            h.plot(intrp['uINTRP'][lvl],intrp['vINTRP'][lvl], marker='.', color='white', 
                   markeredgecolor='black', alpha=1, markersize=30, zorder=5)
            h.ax.annotate(str(int(round(intrp['zINTRP'][lvl]/1000,0))),(intrp['uINTRP'][lvl],intrp['vINTRP'][lvl]), 
                          weight='bold', fontsize=11, color='black',xytext=(0.02,-5),
                          textcoords='offset pixels',horizontalalignment='center',clip_on=True, zorder=5.1) 
    #################################################################
    
    
    
    
    
    
    #################################################################
    ### PLOT HODOGRAPH LINE ###
    #################################################################
    hodo_color = ['purple','red','darkorange','gold','#fff09f']

    h.ax.plot(intrp['uINTRP'][0:10+1],   intrp['vINTRP'][0:10+1],   color=hodo_color[0], linewidth=7, zorder=4, clip_on=True)
    h.ax.plot(intrp['uINTRP'][10:30+1],  intrp['vINTRP'][10:30+1],  color=hodo_color[1], linewidth=7, zorder=4, clip_on=True)
    h.ax.plot(intrp['uINTRP'][30:60+1],  intrp['vINTRP'][30:60+1],  color=hodo_color[2], linewidth=7, zorder=4, clip_on=True)
    h.ax.plot(intrp['uINTRP'][60:90+1],  intrp['vINTRP'][60:90+1],  color=hodo_color[3], linewidth=7, zorder=4, clip_on=True)
    h.ax.plot(intrp['uINTRP'][90:120+1], intrp['vINTRP'][90:120+1], color=hodo_color[4], linewidth=7, zorder=4, clip_on=True) 

    
    
    
    #################################################################
    ### ADD HODOGRAPH ANNOTATIONS ###
    #################################################################
    # BUNKERS STORM MOTION
    if ma.is_masked(kinem['sm_rm']) == False:
        h.ax.text((kinem['sm_rm'][0]+0.5), (kinem['sm_rm'][1]-0.5), 'RM', 
                  weight='bold', ha='left', fontsize=14, zorder=7, alpha=0.9, color=gen_txt_clr)
        h.ax.text((kinem['sm_lm'][0]+0.5), (kinem['sm_lm'][1]-0.5), 'LM', 
                  weight='bold', ha='left', fontsize=14, zorder=7, alpha=0.9, color=gen_txt_clr)
        h.ax.text((kinem['sm_mw'][0]+0.5), (kinem['sm_mw'][1]-0.5), 'MW', 
                  weight='bold', ha='left', fontsize=14, zorder=7, alpha=0.9, color=gen_txt_clr)
    
    # DEVIANT TORNADO MOTION
    if ma.is_masked(kinem['sm_u']) == False:
        h.ax.text(kinem['dtm'][0], (kinem['dtm'][1] + 2), 'DTM', 
                  weight='bold', fontsize=10, zorder=7, color='brown', ha='center')
        h.plot(kinem['dtm'][0], kinem['dtm'][1], marker='v', color='brown', 
               markersize=8, zorder=7, alpha=0.8, ls='', label='DEVIANT TORNADO MOTION')
    
    # ADD SM POINT IF ITS A CUSTOM STORM MOTION
        if str(type(storm_motion)) == "<class 'list'>":
            h.ax.text((kinem['sm_u']+0.5), (kinem['sm_v']-0.5), 'SM', 
                      weight='bold', ha='left', fontsize=14, alpha=0.9, zorder=7, color=gen_txt_clr)
        
        h.ax.arrow(0,0,kinem['sm_u']-0.3, kinem['sm_v']-0.3, linewidth=3, color=gen_txt_clr, alpha=0.2, 
                label='SM Vector', length_includes_head=True, head_width=0.6)
    
    
    # EFFECTIVE INFLOW LAYER SRH FILL
    if ma.is_masked(kinem['eil_z'][0]) == False:
        
        ebot = h.ax.plot((kinem['sm_u'], intrp['uINTRP'][find_nearest(intrp['zINTRP'],kinem['eil_z'][0])]), 
                         (kinem['sm_v'], intrp['vINTRP'][find_nearest(intrp['zINTRP'],kinem['eil_z'][0])]),  
                 linestyle='-', linewidth=2.3, alpha=0.5, zorder=3, color='lightblue', label='Effective Inflow Layer')
        etop = h.ax.plot((kinem['sm_u'], intrp['uINTRP'][find_nearest(intrp['zINTRP'],kinem['eil_z'][1])]), 
                         (kinem['sm_v'], intrp['vINTRP'][find_nearest(intrp['zINTRP'],kinem['eil_z'][1])]),  
                 linestyle='-', linewidth=2.3, alpha=0.5, zorder=3, color='lightblue')
        fill_srh = h.ax.fill(np.append(intrp['uINTRP'][find_nearest(intrp['zINTRP'], kinem['eil_z'][0]):find_nearest(intrp['zINTRP'], kinem['eil_z'][1])+1], kinem['sm_u']), np.append(intrp['vINTRP'][find_nearest(intrp['zINTRP'], kinem['eil_z'][0]):find_nearest(intrp['zINTRP'], kinem['eil_z'][1])+1], kinem['sm_v']),'lightblue',alpha=0.3, zorder=2, label='EIL SRH')
    else:
        fill_srh = h.ax.fill(np.append(intrp['uINTRP'][find_nearest(intrp['zINTRP'], 0):find_nearest(intrp['zINTRP'], 3000)], kinem['sm_u']), 
                     np.append(intrp['vINTRP'][find_nearest(intrp['zINTRP'], 0):find_nearest(intrp['zINTRP'], 3000)], kinem['sm_v']),
                     'lightblue',alpha=0.1, zorder=2, label='0-3 SRH')
    
    
    # MCS MOTION
    h.ax.text(kinem['mcs'][0], kinem['mcs'][1], 'UP', 
              weight='bold', fontsize=12, color='orange', zorder=7, ha='center', alpha=0.5, clip_on=True)
    h.ax.text(kinem['mcs'][2], kinem['mcs'][3], 'DN', 
              weight='bold', fontsize=12, color='orange', zorder=7, ha='center', alpha=0.5, clip_on=True)
    
    
    # STORM MOTION PRINTOUT
    if ma.is_masked(kinem['sm_u']) == False:
        try:
            keys = ['sm_rm', 'sm_lm', 'sm_mw', 'dtm'] 

            speeds = []
            directions = []
            
            speeds.append(mpcalc.wind_speed(kinem['sm_u']*units.kts,kinem['sm_v']*units.kts).m) 
            directions.append(mpcalc.angle_to_direction((np.degrees(np.arctan2(kinem['sm_u'],kinem['sm_v']))), 
                                                        full=False, level=3)) 

            for key in keys:
                speeds.append(mpcalc.wind_speed(kinem[key][0]*units.kts,kinem[key][1]*units.kts).m) 
                directions.append(mpcalc.angle_to_direction((np.degrees(np.arctan2(kinem[key][0],kinem[key][1]))), 
                                                            full=False, level=3))

            speeds.append(mpcalc.wind_speed(kinem['mcs'][0]*units.kts,kinem['mcs'][1]*units.kts).m) 
            directions.append(mpcalc.angle_to_direction((np.degrees(np.arctan2(kinem['mcs'][0],kinem['mcs'][1]))), 
                                                        full=False, level=3))   

            speeds.append(mpcalc.wind_speed(kinem['mcs'][2]*units.kts,kinem['mcs'][3]*units.kts).m) 
            directions.append(mpcalc.angle_to_direction((np.degrees(np.arctan2(kinem['mcs'][2],kinem['mcs'][3]))), 
                                                        full=False, level=3)) 
                    
                    
            #plot Bunkers Storm Motion & DTM Data in box on Hodograph 
            plt.figtext(0.748, 0.919, 
                        f' RM: {directions[1]} @ {mag(speeds[1])} kts\n LM: {directions[2]} @ {mag(speeds[2])} kts\n' + 
                        f' MW: {directions[3]} @ {mag(speeds[3])} kts\nDTM: {directions[4]} @ {mag(speeds[4])} kts\n US: {directions[5]} @ {mag(speeds[5])} kts\n' + 
                        f' DS: {directions[6]} @ {mag(speeds[6])} kts\n',
                        color=gen_txt_clr, fontsize=15, verticalalignment='top', linespacing=2.2, alpha=0.6)  
        except IndexError:
            pass
        
        def sm_str(storm_motion, speeds, directions):
            
            if storm_motion in ['right_moving',
                                'left_moving',
                                'mean_wind']:
                return f"{str.upper(storm_motion.replace('_', ' '))} | {directions[0]} @ {mag(speeds[0])}"

            else:
                return f"USER DEFINED | {directions[0]} @ {mag(speeds[0])}"
            
        plt.figtext(0.748, 0.942, 
                f' SM: {sm_str(storm_motion, speeds, directions)} kts',
                color=gen_txt_clr, weight='bold', fontsize=15, verticalalignment='top', linespacing=2.2, alpha=0.6) 
    ################################################################
    
    
    #########################################################################
    ############################ OTHER AXES ################################# 
    #########################################################################

    
    
    #################################################################
    ### TEMPERATURE ADVECTION ###
    #################################################################
    # GET TOP AND BOTTOM BOUND FOR EACH 'BIN'
    bot_arr, top_arr = np.hsplit(thermo['temp_adv'][1],2)
    bot_arr, top_arr = bot_arr.flatten(), top_arr.flatten()
    #CREATE FIGURE
    temp_adv_ax = plt.axes((0.701,0.1003,0.04,0.85))
    temp_adv_ax.spines["top"].set_color(brdr_clr)
    temp_adv_ax.spines["left"].set_color(brdr_clr)
    temp_adv_ax.spines["right"].set_color(brdr_clr)
    temp_adv_ax.spines["bottom"].set_color(brdr_clr)
    temp_adv_ax.spines["bottom"].set_color(brdr_clr)   
    temp_adv_ax.set_facecolor(bckgrnd_clr)    
    plt.yscale('log')
    temp_adv_ax.set_ylim(1050, 100)
    temp_adv_ax.set_xlim(np.nanmin(thermo['temp_adv'][0])-4, np.nanmax(thermo['temp_adv'][0])+4)
    plt.ylabel(' '), plt.xlabel(' ')
    temp_adv_ax.set_yticklabels([]), temp_adv_ax.set_xticklabels([])
    temp_adv_ax.tick_params(axis='y', length = 0), temp_adv_ax.tick_params(axis='x', length = 0)
    # ADD LINES
    lvls = [1000, 900, 800, 700, 600, 500, 400, 300, 200]
    for lvl in lvls:
        plt.plot((-20,20), (lvl,lvl), color='gray', alpha=0.8, linewidth=1, linestyle='-', clip_on=True)
    # PLOT TEMP ADV BINS 
    for i in range(len(thermo['temp_adv'][0])):
                        if thermo['temp_adv'][0][i] <= 0:
                            temp_adv_bxclr = 'cornflowerblue'
                        elif thermo['temp_adv'][0][i]  > 0:
                            temp_adv_bxclr = 'red'
                        temp_adv_ax.barh(top_arr[i], thermo['temp_adv'][0][i], align='center', 
                                         height=bot_arr[i]-top_arr[i], edgecolor='black', alpha=0.3, color=temp_adv_bxclr)
                        if thermo['temp_adv'][0][i] > 0:
                            temp_adv_ax.annotate((np.round(thermo['temp_adv'][0][i],1)), 
                                                 xy=(0.3, top_arr[i]+10), color=gen_txt_clr, textcoords='data', 
                                                 ha='left', weight='bold')
                        if thermo['temp_adv'][0][i] < 0:
                            temp_adv_ax.annotate((np.round(thermo['temp_adv'][0][i],1)), 
                                                 xy=(-0.3, top_arr[i]+10), color=gen_txt_clr, textcoords='data', 
                                                 ha='right', weight='bold')
    temp_adv_ax.axvline(x=0, color=gen_txt_clr, linewidth=1, linestyle='--', clip_on=True)
    #################################################################





    #################################################################
    ### OMEGA PLOT ###
    #################################################################
    if 'omega' in sounding_data.keys():

        # convert omega (Pa/sec) to Pa/sec raised by an order of magnitude
        # for plotting simplicity
        omega = sounding_data['omega'] * 10

        # create omega axis
        omega_ax = plt.axes((0.127,0.1003,0.3,0.85))
        omega_ax.set_zorder(10)

        # set axis ticks
        plt.yscale('log')
        omega_ax.set_ylim(1050, 100)
        omega_ax.set_xlim(40, -300)

        omega_pos = np.where(omega > 0, omega, np.nan)
        omega_neg = np.where(omega < 0, omega, np.nan)

        # Plot the line for all omega values
        omega_ax.plot(omega, p, linewidth=2, color='black', alpha=0.1)

        # Fill the positive (red) and negative (blue) regions
        omega_ax.fill_betweenx(p, 0, omega_neg, where=omega_neg < 0, color='firebrick', alpha=0.08, interpolate=True)
        omega_ax.fill_betweenx(p, 0, omega_pos, where=omega_pos > 0, color='cornflowerblue', alpha=0.08, interpolate=True)

        omega_ax.axvline(x=0, color='b', alpha=0.2)
        omega_ax.axvline(x=10, color='purple', linestyle='-.', alpha=0.2)
        omega_ax.axvline(x=-10, color='purple', linestyle='-.', alpha=0.2)
        omega_ax.text(10, 150, "+10", weight="bold", alpha=0.4, color='k', fontsize=12, ha='center')
        omega_ax.text(-10, 150, "-10", weight="bold", alpha=0.4, color='k', fontsize=12, ha='center')
        omega_ax.set_axis_off()
        omega_ax.patch.set_alpha(0)
    #################################################################




    
    
    #################################################################
    ### SOUNDING MAP INSET -- TEST ###
    #################################################################
    # BUILD SIMPLE MAP -----------------------------------------------------------------------------------------------------
    if map_zoom > 0:
        proj = ccrs.PlateCarree()
        map_ax = plt.axes((0.7005, 0.1, 0.20, 0.20), projection=proj, zorder=12) # 1, 0.1 | 0.542, 0.75, 0.20, 0.20
        map_ax.spines['geo'].set_color(gen_txt_clr)

        zoom_factor = map_zoom
        map_ax.set_extent([clean_data['site_info']['site-latlon'][1] - zoom_factor,
                           clean_data['site_info']['site-latlon'][1] + zoom_factor,
                           clean_data['site_info']['site-latlon'][0] - zoom_factor,
                           clean_data['site_info']['site-latlon'][0] + zoom_factor*1.6])
        map_ax.set_box_aspect(1) 

        map_ax.add_feature(cfeature.STATES, color=bckgrnd_clr, zorder=1)
        map_ax.add_feature(cfeature.STATES, edgecolor=gen_txt_clr, alpha=0.7, 
                       linestyle='-', linewidth=2, zorder=6)
        map_ax.add_feature(USCOUNTIES, alpha=0.3, edgecolor=gen_txt_clr, 
                           linestyle='-', linewidth=0.2, zorder=5)
        map_ax.add_feature(cfeature.LAKES, color=bckgrnd_clr)
        map_ax.add_feature(cfeature.OCEAN, color=bckgrnd_clr)

        if show_radar == True:
            # BUILD COLOR MAP -----------------------------------------------------------------------------------------------------
            radar_cmap = LinearSegmentedColormap.from_list('custom_cmap',
                        ['#B0C4DE', '#4682B4', '#90EE90', '#228B22', '#FFFF4D', '#E6E600',
                        '#FFC34D', '#E69900', '#FF4D4D', '#E60000', '#FFCCEE', '#FF198C',
                        '#D400FF', '#550080'], N=256)


            # PLOT RADAR DATA -----------------------------------------------------------------------------------------------------
            try:
                radar_data = get_radar_mosaic(sounding_data, map_zoom, radar_time)

                reflectivity = np.ma.masked_array(np.array(radar_data['Reflectivity'])[0,:,:],
                                                  np.array(radar_data['Reflectivity'])[0,:,:]<10)

                map_ax.pcolormesh(radar_data['lon']+0.05, radar_data['lat']+0.05, reflectivity, 
                                  vmin=10, vmax=80, cmap=radar_cmap, alpha=1, zorder=4)

                radar_timestamp = f"{str(radar_data['time'].values[0])[5:10]} | {str(radar_data['time'].values[0])[11:16]}"


                plt.figtext(0.80, 0.11, f'Valid: {radar_timestamp}', 
                        ha='center', alpha=0.9, weight='bold', fontsize=13, zorder=13, color=gen_txt_clr,
                        bbox=dict(facecolor=bckgrnd_clr, alpha=0.8, edgecolor=gen_txt_clr, pad=3))
            except: 
                pass

        # ADD PROFILE LOCATION ------------------------------------------------------------------------------------------------
        map_ax.plot(clean_data['site_info']['site-latlon'][1],clean_data['site_info']['site-latlon'][0], 
                    marker='o', mfc='none', markeredgewidth=1.5, 
                    color=marker_clr, markersize='11', transform=ccrs.PlateCarree(), zorder=10, clip_on=True)
        map_ax.plot(clean_data['site_info']['site-latlon'][1],clean_data['site_info']['site-latlon'][0], 
                    marker='+', color=marker_clr, markeredgewidth=1.5, 
                    markersize='18', transform=ccrs.PlateCarree(), zorder=10, clip_on=True)

    #################################################################
    
    
    #################################################################
    ### STREAMWISENESS AND RH W/HGT ###
    #################################################################
    # PLOT AXIS/LOC
    strmws_ax = plt.axes((0.945, -0.13, 0.065, 0.23))
    plt.figtext(0.978, 0.07, f'Streamwiseness\n of ζ (%)', color=gen_txt_clr, weight='bold', fontsize=10, ha='center', alpha=0.9)
    strmws_ax.spines["top"].set_color(brdr_clr)
    strmws_ax.spines["left"].set_color(brdr_clr)
    strmws_ax.spines["right"].set_color(brdr_clr)
    strmws_ax.spines["bottom"].set_color(brdr_clr)
    strmws_ax.spines["bottom"].set_color(brdr_clr)   
    strmws_ax.set_facecolor(bckgrnd_clr) 

    #YTICKS
    strmws_ax.set_ylim(0, 3000)
    strmws_ax.set_yticklabels([])
    strmws_ax.set_ylabel(' ')
    strmws_ax.tick_params(axis='y', length = 0)
    strmws_ax.grid(True, axis='y')
    strmws_ax.tick_params(axis="x",direction="in", pad=-12)

    #XTICKS
    strmws_ax.set_xlim(40, 102)
    strmws_ax.set_xticks([50, 70, 90])
    strmws_ax.set_xticklabels([50, 70, 90], weight='bold', alpha=0.5, fontstyle='italic', color=gen_txt_clr)
    strmws_ax.set_xlabel(' ')

    #HGT LABLES 
    strmws_ax.text(47, 502 , '.5 km', fontsize=10, weight='bold', alpha=1, color=gen_txt_clr)
    strmws_ax.text(47, 1002, '1 km', fontsize=10, weight='bold', alpha=1, color=gen_txt_clr)
    strmws_ax.text(47, 1502, '1.5 km', fontsize=10, weight='bold', alpha=1, color=gen_txt_clr)
    strmws_ax.text(47, 2002, '2 km', fontsize=10, weight='bold', alpha=1, color=gen_txt_clr)
    #strmws_ax.text(47, 2502, '2.5 km', fontsize=9, weight='bold', alpha=0.9, color=gen_txt_clr)

    if ma.is_masked(kinem['sm_u']) == False:
        strmws_ax.plot(kinem['swv_perc'][0:11],  intrp['zINTRP'][0:11],  color=hodo_color[0], lw=3, clip_on=True)
        strmws_ax.plot(kinem['swv_perc'][10:31], intrp['zINTRP'][10:31], color=hodo_color[1], lw=3, clip_on=True)
        strmws_ax.fill_betweenx(intrp['zINTRP'][0:31], kinem['swv_perc'][0:31],
                                color='cornflowerblue', linewidth=0, alpha=0.2, clip_on=True)

    
        if ma.is_masked(kinem['eil_z'][0]) == False:
            strmws_ax.fill_between(x=(40,102), y1=kinem['eil_z'][0], y2=kinem['eil_z'][1], color='lightblue', alpha=0.2)
    else:
        warnings.warn("Streamwiseness could not be plotted (no valid storm motion/not enough data)", Warning)
    #################################################################

    
    
    
    #################################################################
    ### VORTICITY W/HGT ###
    #################################################################
    vort_ax = plt.axes((1.01, -0.13, 0.066, 0.23))
    plt.figtext(1.043, 0.05, f'Total ζ &\n Streamwise ζ \n(/sec)', color=gen_txt_clr, weight='bold',
                fontsize=12, ha='center', alpha=0.9)
    vort_ax.spines["top"].set_color(brdr_clr)
    vort_ax.spines["left"].set_color(brdr_clr)
    vort_ax.spines["right"].set_color(brdr_clr)
    vort_ax.spines["bottom"].set_color(brdr_clr)
    vort_ax.spines["bottom"].set_color(brdr_clr)   
    vort_ax.set_facecolor(bckgrnd_clr) 

    #YTICKS
    vort_ax.tick_params(axis='y', length = 0)
    vort_ax.grid(True, axis='y')
    vort_ax.set_ylim(0, 3000)
    vort_ax.set_yticklabels([])
    vort_ax.set_ylabel(' ')
    vort_ax.tick_params(axis="x",direction="in", pad=-12)
        
    if ma.is_masked(kinem['sm_u']) == False: 
        #XTICKS
        vort_ax.set_xlabel(' ')
    
        vort_ax.set_xlim(0, 0.06)
        vort_ax.set_xticks([.01, .03, .05])
        vort_ax.set_xticklabels([".01", ".03", ".05"],
                                weight='bold', alpha=0.5, fontstyle='italic', color=gen_txt_clr)
 
        vort_ax.plot(kinem['swv'][0:31],  intrp['zINTRP'][0:31], color='orange', linewidth=3, alpha=0.8, label='SW ζ')
        vort_ax.plot(kinem['vort'][0:31], intrp['zINTRP'][0:31], color=gen_txt_clr,  linewidth=4, alpha=0.4, label='Total ζ')
        vort_ax.fill_betweenx(intrp['zINTRP'][0:31], kinem['vort'][0:31], where=(kinem['vort'][0:31] > kinem['swv'][0:31]),
                                color=gen_txt_clr, linewidth=0, alpha=0.1, clip_on=True)
        vort_ax.fill_betweenx(intrp['zINTRP'][0:31], kinem['swv'][0:31],
                                color='orange', linewidth=0, alpha=0.2, clip_on=True)
        
        if ma.is_masked(kinem['eil_z'][0]) == False:
            vort_ax.fill_between(x=(0, 0.065), y1=kinem['eil_z'][0],
                                 y2=kinem['eil_z'][1], color='lightblue', alpha=0.2)
    else:
        warnings.warn("Total Vorticity could not be plotted (no valid storm motion/not enough data)", Warning)
    #################################################################
    
    
    
    
    
    #################################################################
    ### SRW W/HGT ###
    #################################################################
    wind_ax = plt.axes((1.0755, -0.13, 0.066, 0.23))
    plt.figtext(1.108, 0.06, f'Storm Relative\nWind (kts)', weight='bold', color=gen_txt_clr, fontsize=12, ha='center', alpha=0.9)
    wind_ax.spines["top"].set_color(brdr_clr)
    wind_ax.spines["left"].set_color(brdr_clr)
    wind_ax.spines["right"].set_color(brdr_clr)
    wind_ax.spines["bottom"].set_color(brdr_clr)
    wind_ax.spines["bottom"].set_color(brdr_clr)   
    wind_ax.set_facecolor(bckgrnd_clr) 
    plt.ylabel(' ')
    plt.xlabel(' ')

    #YTICKS
    wind_ax.set_ylim(0, 3000)
    wind_ax.grid(True, axis='y')
    wind_ax.set_yticklabels([])
    wind_ax.tick_params(axis='y', length = 0)
    wind_ax.tick_params(axis="x",direction="in", pad=-12)

    wind_ax.text(40, thermo['sb_lcl_z'], '-LCL-', fontsize=10, weight='bold',
                 alpha=1, color=gen_txt_clr, clip_on=True)
    if thermo['mu_lfc_z'] < 2500:
        wind_ax.text(40, thermo['mu_lfc_z'], '-LFC-', fontsize=10, weight='bold',
                     alpha=1, color=gen_txt_clr, clip_on=True)


    if ma.is_masked(kinem['sm_u']) == False:
        #XTICKS
        wind_ax.set_xlim(10, 50)
        wind_ax.set_xticks([20, 30, 40])
        wind_ax.set_xticklabels([20, 30, 40], weight='bold',
                                alpha=0.5, fontstyle='italic', color=gen_txt_clr)

        #PLOT SR WIND  
        wind_ax.plot(kinem['srw'][0:11],  intrp['zINTRP'][0:11],  color=hodo_color[0], clip_on=True, 
                     linewidth=3, alpha=0.8, label='0-1 SR Wind')
        wind_ax.plot(kinem['srw'][10:31], intrp['zINTRP'][10:31], color=hodo_color[1], clip_on=True,
                     linewidth=3, alpha=0.8, label='1-3 SR Wind')
        wind_ax.fill_betweenx(intrp['zINTRP'][0:31], kinem['srw'][0:31],
                                color='cornflowerblue', linewidth=0, alpha=0.2, clip_on=True)
        
        if ma.is_masked(kinem['eil_z'][0]) == False:
            wind_ax.fill_between(x=(5, 55), y1=kinem['eil_z'][0], y2=kinem['eil_z'][1],
                                 color='lightblue', alpha=0.2)
    else:
        warnings.warn("Storm Relative Wind could not be plotted (no valid storm motion/not enough data)", Warning)
        
    #################################################################
    
    
    
    
    
    if show_theta:
        #############################################################
        ### THETA & THETA E W/HGT ###
        #############################################################
        #PLOT AXES/LOC
        theta_ax = plt.axes((1.141, -0.13, 0.0653, 0.23))
        plt.figtext(1.175, 0.06, f'Theta-e &\nTheta (K)', weight='bold', color=gen_txt_clr, fontsize=12, ha='center', alpha=0.9)
        theta_ax.spines["top"].set_color(brdr_clr)
        theta_ax.spines["left"].set_color(brdr_clr)
        theta_ax.spines["right"].set_color(brdr_clr)
        theta_ax.spines["bottom"].set_color(brdr_clr)
        theta_ax.spines["bottom"].set_color(brdr_clr)   
        theta_ax.set_facecolor(bckgrnd_clr) 

        maxtheta = intrp['thetaeINTRP'][0:31].max()
        mintheta = intrp['thetaINTRP'][0:31].min()

        #YTICKS
        theta_ax.set_ylim(intrp['zINTRP'][0], 3000)
        theta_ax.set_yticklabels([])
        plt.ylabel(' ')
        theta_ax.tick_params(axis='y', length = 0)
        theta_ax.grid(True, axis='y')

        #XTICKS
        theta_ax.set_xlim(mintheta - 5, maxtheta + 5)
        theta_ax.set_xticks([(mintheta), (maxtheta)])
        theta_ax.set_xticklabels([int(mintheta), int(maxtheta)], weight='bold', alpha=0.5, fontstyle='italic', color=gen_txt_clr)

        theta_ax.tick_params(axis="x", direction="in", pad=-12)
        theta_ax.set_xlabel(' ')

        #PLOT THETA VS HGT
        plt.plot(intrp['thetaINTRP'], intrp['zINTRP'], color='purple', linewidth=3.5, alpha=0.5, clip_on=True)
        plt.plot(intrp['thetaeINTRP'], intrp['zINTRP'], color='purple', linewidth=3.5, alpha=0.8, clip_on=True)

        if ma.is_masked(kinem['eil_z'][0]) == False:
            theta_ax.fill_between(x=(mintheta - 5, maxtheta + 5), y1=kinem['eil_z'][0], 
                                y2=kinem['eil_z'][1], color='lightblue', alpha=0.2)
    else:
        #############################################################
        ### CIN & 3CAPE W/HGT ###
        #############################################################
        cin_color, cape_color = 'teal', 'tab:orange'

        mincin = intrp['cin_profileINTRP'][0:31].min()
        max3cape = intrp['3cape_profileINTRP'][0:31].max()

        # PLOT AXES/LOC
        inflow_ax = plt.axes((1.141, -0.13, 0.0653, 0.23))
        plt.figtext(1.173, 0.05, f'Stepwise\n CIN & CAPE\n (J/kg)', weight='bold',
                    color=gen_txt_clr, fontsize=12, ha='center', alpha=0.9)
        inflow_ax.spines["top"].set_color(brdr_clr)
        inflow_ax.spines["left"].set_color(brdr_clr)
        inflow_ax.spines["right"].set_color(brdr_clr)
        inflow_ax.spines["bottom"].set_color(brdr_clr)
        inflow_ax.set_facecolor(bckgrnd_clr)

        ### AXIS 1 -- CIN ###

        # YAXIS PARAMS AX1 & AX2
        inflow_ax.set_ylim(intrp['zINTRP'][0], 3000)
        inflow_ax.set_yticklabels([])
        inflow_ax.set_ylabel(' ')
        inflow_ax.tick_params(axis='y', length=0)
        inflow_ax.grid(True, axis='y')
        inflow_ax.axvline(0, linestyle=':', alpha=.8)

        # XAXIS PARAMS AX1
        inflow_ax.set_xlim(-300, 300)
        inflow_ax.set_xticks([-200, -100])
        inflow_ax.set_xticklabels([-200, -100], weight='bold', alpha=0.5, rotation=60,
                                  fontstyle='italic', color=gen_txt_clr, zorder=10)
        inflow_ax.tick_params(axis="x", direction="in", pad=-25)
        inflow_ax.set_xlabel(' ')
        plt.xticks(rotation=45)

        # PLOT CIN VS HGT
        inflow_ax.fill_betweenx(intrp['zINTRP'], intrp['cin_profileINTRP'],
                                color=cin_color, linewidth=0, alpha=0.5, clip_on=True)

        if ma.is_masked(kinem['eil_z'][0]) == False:
            inflow_ax.fill_between(x=(-305, 305), y1=kinem['eil_z'][0],
                                   y2=kinem['eil_z'][1], color='lightblue', alpha=0.2)

        ### AXIS 2 -- CAPE ###
        inflow_ax2 = inflow_ax.twiny()
        inflow_ax2.spines["top"].set_color(brdr_clr)
        inflow_ax2.spines["left"].set_color(brdr_clr)
        inflow_ax2.spines["right"].set_color(brdr_clr)
        inflow_ax2.spines["bottom"].set_color(brdr_clr)
        inflow_ax2.set_facecolor(bckgrnd_clr)

        # XAXIS PARAMS AX2
        inflow_ax2.set_xlim(-3000, 3000)
        inflow_ax2.set_xticks([0, 1000, 2000])
        inflow_ax2.set_xticklabels(['0 ', '1k', '2k'], weight='bold', alpha=0.5,
                                  fontstyle='italic', color=gen_txt_clr, zorder=10)
        inflow_ax2.xaxis.set_ticks_position('bottom')
        inflow_ax2.tick_params(axis="x", direction="in", pad=-20)
        inflow_ax2.set_xlabel(' ')
        plt.xticks(rotation=45)

        # YAXIS PARAMS AX2
        inflow_ax2.set_ylim(intrp['zINTRP'][0], 3000)
        inflow_ax2.set_yticklabels([])
        inflow_ax2.set_ylabel(' ')
        inflow_ax2.tick_params(axis="y", length=0)
        inflow_ax2.grid(True, axis='y')

        # PLOT 3CAPE VS HGT
        inflow_ax2.fill_betweenx(intrp['zINTRP'], intrp['cape_profileINTRP'],
                                 color=cape_color, linewidth=0, alpha=0.8, clip_on=True)
        #############################################################

    #########################################################################
    ############################## TEXT PLOTS ###############################
    #########################################################################

    
    
    #################################################################
    ### BOXES FOR TEXT ###
    #################################################################
    #THERMO BOX
    fig.patches.extend([plt.Rectangle((0.124, -0.08),0.556,0.18,
                                      edgecolor=brdr_clr, facecolor=bckgrnd_clr, linewidth=1, alpha=1,
                                      transform=fig.transFigure, figure=fig)])
    # MOISTURE BOX
    fig.patches.extend([plt.Rectangle((0.53, -0.13),0.15,0.23,
                                      edgecolor=brdr_clr, facecolor=bckgrnd_clr, linewidth=1, alpha=1,
                                      transform=fig.transFigure, figure=fig)])
    # KINEMATICS BOX
    fig.patches.extend([plt.Rectangle((0.68, -0.13),0.265,0.23,
                                      edgecolor=brdr_clr, facecolor=bckgrnd_clr, linewidth=1, alpha=1,
                                      transform=fig.transFigure, figure=fig)])

    # OTHER KINEMATICS BOX
    fig.patches.extend([plt.Rectangle((1.106, 0.1),0.10,0.14,
                                      edgecolor=brdr_clr, facecolor=bckgrnd_clr, linewidth=1, alpha=1,
                                      transform=fig.transFigure, figure=fig)])
    #################################################################
    
    
    
    #################################################################
    ### THERMODYNAMICS ###
    #################################################################
    plt.figtext( 0.17, 0.07, 'SR-ECAPE       CAPE         6CAPE         3CAPE          CIN            LCL', 
                color=gen_txt_clr, weight='bold', fontsize=15)
    #SBCAPE
    plt.figtext( 0.13,  0.04, f"SB:", weight='bold',   fontsize=15, color=gen_txt_clr)
    plt.figtext( 0.17,  0.04, f"{mag(thermo['sb_ecape'])} J/kg",  fontsize=15, color='firebrick', weight='bold')
    plt.figtext( 0.235,  0.04, f"{mag(thermo['sbcape'])} J/kg",  fontsize=15, color='orangered', weight='bold')
    plt.figtext( 0.30,  0.04, f"{mag(thermo['sb6cape'])} J/kg", fontsize=15, color='orangered', weight='bold')
    plt.figtext( 0.364, 0.04, f"{mag(thermo['sb3cape'])} J/kg", fontsize=15, color='orangered', weight='bold')
    plt.figtext( 0.425, 0.04, f"{mag(thermo['sbcin'])} J/kg",   fontsize=15, color='orangered', weight='bold')
    plt.figtext( 0.485, 0.04, f"{mag(thermo['sb_lcl_z'])} m",   fontsize=15, color='orangered', weight='bold')
    #MUCAPE
    plt.figtext( 0.13,  0.01, f"MU:", weight='bold',   fontsize=15, color=gen_txt_clr)
    plt.figtext( 0.17,  0.01, f"{mag(thermo['mu_ecape'])} J/kg",  fontsize=15, color='firebrick', weight='bold')
    plt.figtext( 0.235,  0.01, f"{mag(thermo['mucape'])} J/kg",  fontsize=15, color='orangered', weight='bold')
    plt.figtext( 0.30,  0.01, f"{mag(thermo['mu6cape'])} J/kg", fontsize=15, color='orangered', weight='bold')
    plt.figtext( 0.364, 0.01, f"{mag(thermo['mu3cape'])} J/kg", fontsize=15, color='orangered', weight='bold')
    plt.figtext( 0.425, 0.01, f"{mag(thermo['mucin'])} J/kg",   fontsize=15, color='orangered', weight='bold')
    plt.figtext( 0.485, 0.01, f"{mag(thermo['mu_lcl_z'])} m",   fontsize=15, color='orangered', weight='bold')
    #MLCAPE
    plt.figtext( 0.13,  -0.02, f"ML:", weight='bold',   fontsize=15, color=gen_txt_clr)
    plt.figtext( 0.17,  -0.02, f"{mag(thermo['ml_ecape'])} J/kg",  fontsize=15, color='firebrick', weight='bold')
    plt.figtext( 0.235,  -0.02, f"{mag(thermo['mlcape'])} J/kg",  fontsize=15, color='orangered', weight='bold')
    plt.figtext( 0.30,  -0.02, f"{mag(thermo['ml6cape'])} J/kg", fontsize=15, color='orangered', weight='bold')
    plt.figtext( 0.364, -0.02, f"{mag(thermo['ml3cape'])} J/kg", fontsize=15, color='orangered', weight='bold')
    plt.figtext( 0.425, -0.02, f"{mag(thermo['mlcin'])} J/kg",   fontsize=15, color='orangered', weight='bold')
    plt.figtext( 0.485, -0.02, f"{mag(thermo['ml_lcl_z'])} m",   fontsize=15, color='orangered', weight='bold')
    #DCAPE
    plt.figtext( 0.13,  -0.046, f"DCAPE:", weight='bold', fontsize=15, color=gen_txt_clr)
    plt.figtext( 0.175, -0.046, f"{mag(thermo['dcape'])} J/kg", fontsize=15, color='cornflowerblue', alpha=0.7, weight='bold')
    plt.figtext( 0.13,  -0.066, f"DCIN:", weight='bold', fontsize=15, color=gen_txt_clr)
    plt.figtext( 0.175, -0.066, f"{mag(thermo['dcin'])} J/kg", fontsize=15, color='cornflowerblue', alpha=0.7, weight='bold')
    # MUNCAPE
    plt.figtext( 0.24,  -0.061, f"MUNCAPE:", weight='bold', fontsize=15, color=gen_txt_clr)
    plt.figtext( 0.295, -0.061, f" {(thermo['mu_ncape'])}", fontsize=15, color='brown', alpha=0.7, weight='bold')
    #0-3KM LR
    plt.figtext( 0.335, -0.061, f"Γ₀₋₃:", weight='bold', fontsize=16, color=gen_txt_clr)
    plt.figtext( 0.365, -0.061, f"{mag(thermo['lr_03km'])} Δ°C/km", fontsize=15, color='saddlebrown', weight='bold')
    #3-6km LR
    plt.figtext( 0.425, -0.061, f"Γ₃₋₆:", weight='bold', fontsize=16, color=gen_txt_clr)
    plt.figtext( 0.458, -0.061, f"{mag(thermo['lr_36km'])} Δ°C/km", fontsize=15, color='saddlebrown', weight='bold')
    #################################################################
    

    
    #################################################################
    ### MOSITURE ###
    ################################################################# 
    plt.figtext( 0.55,  0.07, '          RH(%)       ω', color=gen_txt_clr, weight='bold', fontsize=15)
    plt.figtext( 0.54,  0.04, f"0-.5ₖₘ:", fontsize=15, color=gen_txt_clr, weight='bold')
    plt.figtext( 0.541, 0.06, f"___", fontsize=15, color='black' , alpha=1, weight='bold')
    plt.figtext( 0.585, 0.04, f"{mag(general['rh_0_500'])} %", fontsize=15, color='forestgreen', alpha=1, weight='bold')
    plt.figtext( 0.62,  0.04, f"{mag_round(general['w_0_500'],1)} g/kg", fontsize=15, color='forestgreen', alpha=1, weight='bold')

    plt.figtext( 0.54,  0.01, f"0-1ₖₘ:", fontsize=15, color=gen_txt_clr, weight='bold',)
    plt.figtext( 0.541, 0.03, f"___", fontsize=15, color='black' , alpha=0.9, weight='bold')
    plt.figtext( 0.585, 0.01, f"{mag(general['rh_0_1000'])} %", fontsize=15, color='forestgreen', alpha=0.9, weight='bold')
    plt.figtext( 0.62,  0.01, f"{mag_round(general['w_0_1000'],1)} g/kg", fontsize=15, color='forestgreen', alpha=0.9, weight='bold')

    plt.figtext( 0.54,  -0.02, f"1-3ₖₘ:", fontsize=15, color=gen_txt_clr, weight='bold',)
    plt.figtext( 0.541, -0.00, f"___", fontsize=15, color='black' , alpha=0.8, weight='bold')
    plt.figtext( 0.585, -0.02, f"{mag(general['rh_1_3000'])} %", fontsize=15, color='forestgreen', alpha=0.8, weight='bold')
    plt.figtext( 0.62,  -0.02, f"{mag_round(general['w_1_3000'],1)} g/kg", fontsize=15, color='forestgreen', alpha=0.8, weight='bold')

    plt.figtext( 0.54,  -0.05, f"3-6ₖₘ:", fontsize=15, color=gen_txt_clr, weight='bold',)
    plt.figtext( 0.541, -0.03, f"___", fontsize=15, color='black' , alpha=0.7, weight='bold')
    plt.figtext( 0.585, -0.05, f"{mag(general['rh_3_6000'])} %", fontsize=15, color='forestgreen', alpha=0.7, weight='bold')
    plt.figtext( 0.62,  -0.05, f"{mag_round(general['w_3_6000'],1)} g/kg", fontsize=15, color='forestgreen', alpha=0.7, weight='bold')

    plt.figtext( 0.54,  -0.08, f"6-9ₖₘ:", fontsize=15, color=gen_txt_clr, weight='bold',)
    plt.figtext( 0.541, -0.06, f"___", fontsize=15, color='black', alpha=0.64, weight='bold')
    plt.figtext( 0.585, -0.08, f"{mag(general['rh_6_9000'])} %", fontsize=15, color='forestgreen', alpha=0.64, weight='bold')
    plt.figtext( 0.62,  -0.08, f"{mag_round(general['w_6_9000'],1)} g/kg", fontsize=15, color='forestgreen',alpha=0.64, weight='bold')

    plt.figtext( 0.54,  -0.102, f"PWAT:", fontsize=12, weight='bold', color=gen_txt_clr)
    plt.figtext( 0.59, -0.102, f"{str(general['pwat'])} in", fontsize=12, color='darkgreen', ha='center', weight='bold')

    plt.figtext( 0.625, -0.102, f"WB:", fontsize=12, weight='bold', color=gen_txt_clr, ha='center')
    plt.figtext( 0.655,  -0.102, f"{mag(general['wet_bulb'][0])} °C", fontsize=12, color='darkgreen', ha='center', weight='bold')
    
    plt.figtext( 0.54,  -0.122, f"FRZ:", fontsize=12, weight='bold', color=gen_txt_clr)
    plt.figtext( 0.59, -0.122, f"{mag(general['frz_pt_z'])}m", fontsize=12, color='cornflowerblue', ha='center', weight='bold')
    
    plt.figtext( 0.625, -0.122, f"WB0:", fontsize=12, weight='bold', color=gen_txt_clr, ha='center')
    plt.figtext( 0.656,  -0.122, f"{mag(general['wb_frz_pt_z'])}m", fontsize=12, color='cornflowerblue', ha='center', weight='bold')
    #################################################################
    
    
    
    
    #################################################################
    ### KINEMATICS ###
    #################################################################   
    met_per_sec = (units.m*units.m)/(units.sec*units.sec)

    plt.figtext( 0.735, 0.07, 'BWD       SRH       SRW   SWζ%    SWζ', color=gen_txt_clr, weight='bold', fontsize=15, alpha=0.8)

    plt.figtext( 0.689, 0.04, f"0-.5ₖₘ:", weight='bold', fontsize=15, color=gen_txt_clr, alpha=0.9)
    plt.figtext( 0.732, 0.04, f"{mag(kinem['shear_0_to_500'])} kt", fontsize=15, color='#6495ED', weight='bold')
    plt.figtext( 0.769, 0.04, f"{mag(kinem['srh_0_to_500'])* met_per_sec:~P}", fontsize=15, color='#6495ED', weight='bold')
    plt.figtext( 0.828, 0.04, f"{mag(kinem['srw_0_to_500'])} kt", fontsize=15, color='#6495ED', weight='bold')
    plt.figtext( 0.870, 0.04, f"{mag(kinem['swv_perc_0_to_500'])}", fontsize=15, color='#6495ED', weight='bold')
    plt.figtext( 0.905, 0.04, f"{mag_round(kinem['swv_0_to_500'], 3)}", fontsize=15, color='#6495ED', weight='bold')

    plt.figtext( 0.689, 0.01, f"0-1ₖₘ:", weight='bold', fontsize=15, color=gen_txt_clr, alpha=0.9)
    plt.figtext( 0.732, 0.01, f"{mag(kinem['shear_0_to_1000'])} kt", fontsize=15, color='#7B68EE', weight='bold')
    plt.figtext( 0.769, 0.01, f"{mag(kinem['srh_0_to_1000'])* met_per_sec:~P}", fontsize=15, color='#7B68EE', weight='bold')
    plt.figtext( 0.828, 0.01, f"{mag(kinem['srw_0_to_1000'])} kt", fontsize=15, color='#7B68EE', weight='bold')
    plt.figtext( 0.870, 0.01, f"{mag(kinem['swv_perc_0_to_1000'])}", fontsize=15, color='#7B68EE', weight='bold')
    plt.figtext( 0.905, 0.01, f"{mag_round(kinem['swv_0_to_1000'], 3)}", fontsize=15, color='#7B68EE', weight='bold')

    plt.figtext( 0.689, -0.02, f"1-3ₖₘ:", weight='bold', fontsize=15, color=gen_txt_clr, alpha=0.9)
    plt.figtext( 0.732, -0.02, f"{mag(kinem['shear_1_to_3000'])} kt", fontsize=15, color='#4169E1', weight='bold')
    plt.figtext( 0.769, -0.02, f"{mag(kinem['srh_1_to_3000'])* met_per_sec:~P}", fontsize=15, color='#4169E1', weight='bold')
    plt.figtext( 0.828, -0.02, f"{mag(kinem['srw_1_to_3000'])} kt", fontsize=15, color='#4169E1', weight='bold')
    plt.figtext( 0.870, -0.02, f"{mag(kinem['swv_perc_1_to_3000'])}", fontsize=15, color='#4169E1', weight='bold')
    plt.figtext( 0.905, -0.02, f"{mag_round(kinem['swv_1_to_3000'], 3)}", fontsize=15, color='#4169E1', weight='bold')

    plt.figtext( 0.689, -0.05, f"3-6ₖₘ:", weight='bold', fontsize=15, color=gen_txt_clr, alpha=0.9)
    plt.figtext( 0.732, -0.05, f"{mag(kinem['shear_3_to_6000'])} kt", fontsize=15, color='#0000CD', weight='bold')
    plt.figtext( 0.769, -0.05, f"{mag(kinem['srh_3_to_6000'])* met_per_sec:~P}", fontsize=15, color='#0000CD', weight='bold')
    plt.figtext( 0.828, -0.05, f"{mag(kinem['srw_3_to_6000'])} kt", fontsize=15, color='#0000CD', weight='bold')
    plt.figtext( 0.870, -0.05, f"{mag(kinem['swv_perc_3_to_6000'])}", fontsize=15, color='#0000CD', weight='bold')
    plt.figtext( 0.905, -0.05, f"{mag_round(kinem['swv_3_to_6000'], 3)}", fontsize=15, color='#0000CD', weight='bold')
    
    plt.figtext( 0.689, -0.08, f"6-9ₖₘ:", weight='bold', fontsize=15, color=gen_txt_clr, alpha=0.9)
    plt.figtext( 0.732, -0.08, f"{mag(kinem['shear_6_to_9000'])} kt", fontsize=15, color='#00008B', weight='bold')
    plt.figtext( 0.769, -0.08, f"{mag(kinem['srh_6_to_9000'])* met_per_sec:~P}", fontsize=15, color='#00008B', weight='bold')
    plt.figtext( 0.828, -0.08, f"{mag(kinem['srw_6_to_9000'])} kt", fontsize=15, color='#00008B', weight='bold')
    plt.figtext( 0.870, -0.08, f"{mag(kinem['swv_perc_6_to_9000'])}", fontsize=15, color='#00008B', weight='bold')
    plt.figtext( 0.905, -0.08, f"{mag_round(kinem['swv_6_to_9000'], 3)}", fontsize=15, color='#00008B', weight='bold')
    
    plt.figtext( 0.689, -0.11, f"EIL:", weight='bold', fontsize=15, color=gen_txt_clr, alpha=0.9)
    plt.figtext( 0.732, -0.11, f"{mag(kinem['shear_eil'])} kt", fontsize=15, color='#000080', weight='bold')
    plt.figtext( 0.769, -0.11, f"{mag(kinem['srh_eil'])* met_per_sec:~P}", fontsize=15, color='#000080', weight='bold')
    plt.figtext( 0.828, -0.11, f"{mag(kinem['srw_eil'])} kt", fontsize=15, color='#000080', weight='bold')
    plt.figtext( 0.870, -0.11, f"{mag(kinem['swv_perc_eil'])}", fontsize=15, color='#000080', weight='bold')
    plt.figtext( 0.905, -0.11, f"{mag_round(kinem['swv_eil'], 3)}", fontsize=15, color='#000080', weight='bold')
    #################################################################




    #################################################################
    ### OTHER KINEMATICS ###
    #################################################################
    plt.figtext(1.136, 0.219, f"  0-3ₖₘ SRH,", weight='bold', fontsize=13, color=gen_txt_clr, alpha=0.9, ha='center')
    plt.figtext(1.136, 0.196, f"{mag(kinem['srh_0_to_3000']) * met_per_sec:~P}", weight='bold', fontsize=13,
                color='#4169E1', alpha=0.9, ha='center')

    plt.figtext(1.183, 0.219, f"BWD", weight='bold', fontsize=13, color=gen_txt_clr, alpha=0.9, ha='center')
    plt.figtext(1.183, 0.196, f"{mag(kinem['shear_0_to_3000'])} kt", weight='bold', fontsize=13, color='#4169E1',
                alpha=0.9, ha='center')

    plt.figtext(1.136, 0.173, f"  0-6ₖₘ SRH,", weight='bold', fontsize=13, color=gen_txt_clr, alpha=0.9, ha='center')
    plt.figtext(1.136, 0.152, f"{mag(kinem['srh_0_to_6000']) * met_per_sec:~P}", weight='bold', fontsize=13,
                color='#4169E1', alpha=0.9, ha='center')

    plt.figtext(1.183, 0.173, f"BWD", weight='bold', fontsize=13, color=gen_txt_clr, alpha=0.9, ha='center')
    plt.figtext(1.183, 0.152, f"{mag(kinem['shear_0_to_6000'])} kt", weight='bold', fontsize=13, color='#4169E1',
                alpha=0.9, ha='center')

    plt.figtext(1.12, 0.128, f"SCP", weight='bold', fontsize=13, color=gen_txt_clr, alpha=0.9, ha='center')
    plt.figtext(1.12, 0.1085, f"{mag(kinem['eil_scp'])}", weight='bold', fontsize=13, color='firebrick', alpha=0.9,
                ha='center')

    plt.figtext(1.145, 0.128, f"STP", weight='bold', fontsize=13, color=gen_txt_clr, alpha=0.9, ha='center')
    plt.figtext(1.145, 0.1085, f"{mag(kinem['eil_stp'])}", weight='bold', fontsize=13, color='firebrick', alpha=0.9,
                ha='center')

    plt.figtext(1.17, 0.129, f"EHI", weight='bold', fontsize=13, color=gen_txt_clr, alpha=0.9, ha='center')
    plt.figtext(1.17, 0.121, f"₀₋₁ₖₘ", weight='bold', fontsize=12, color=gen_txt_clr, alpha=0.9, ha='center')
    plt.figtext(1.17, 0.108, f"{mag(kinem['ehi_0_to_1000'])}", weight='bold', fontsize=13, color='orangered', alpha=0.9,
                ha='center')

    plt.figtext(1.195, 0.129, f"EHI", weight='bold', fontsize=13, color=gen_txt_clr, alpha=0.9, ha='center')
    plt.figtext(1.195, 0.121, f"₀₋₃ₖₘ", weight='bold', fontsize=12, color=gen_txt_clr, alpha=0.9, ha='center')
    plt.figtext(1.195, 0.108, f"{mag(kinem['ehi_0_to_3000'])}", weight='bold', fontsize=13, color='orangered',
                alpha=0.9, ha='center')
    #################################################################




    #########################################################################
    ############################# PLOT EXTRAS ############################### 
    #########################################################################   
    plt.figtext( 0.125, 0.985, top_title, weight='bold', ha='left', fontsize=30, color=gen_txt_clr)
    plt.figtext( 0.125, 0.959, left_title, ha='left', fontsize=23, color=gen_txt_clr)
    plt.figtext( 1.22, 0.959, right_title, ha='right', fontsize=23, color=gen_txt_clr)
    skewleg1 = skew.ax.legend(loc='upper left', framealpha=0.5, labelcolor=gen_txt_clr, facecolor=bckgrnd_clr)
    
    # plot author credit information 
    plt.figtext( 0.34, -0.105, 'SOUNDERPY VERTICAL PROFILE ANALYSIS TOOL', 
                fontsize=20, ha='center', color='cornflowerblue', weight='bold', alpha=0.6)
    plt.figtext( 0.34, -0.125, '(C) KYLE J GILLETT | sounderpysoundings.anvil.app', 
                fontsize=15, ha='center', color='cornflowerblue', weight='bold', alpha=0.6)

    # sounderpy logo
    img = Image.open(urlopen('https://user-images.githubusercontent.com/100786530/251580013-2e9477c9-e36a-4163-accb-fe46780058dd.png'))
    #                  side-side  up-down  size   size
    imgax = fig.add_axes([0.115, -0.135, 0.05, 0.05], anchor='SE', zorder=3)
    imgax.imshow(img)
    imgax.axis('off')
    
    plt.tight_layout()
    
    print('    > COMPLETE --------')
    elapsed_time = time.time() - st
    print('    > RUNTIME:', time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
    
    return plt
#########################################################################




















#########################################################################
############################ HODOGRAPH ##################################

def __full_hodograph(clean_data, dark_mode, storm_motion, sr_hodo, modify_sfc):


    if dark_mode == True:
        gen_txt_clr = 'white'
        bckgrnd_clr = 'black'
        brdr_clr    = 'white'
        barb_clr    = 'white'
        shade_alpha = 0.06
    else: 
        gen_txt_clr = 'black'
        bckgrnd_clr = 'white'
        brdr_clr    = 'black'
        barb_clr    = 'black'
        shade_alpha = 0.02
    
    
    # record process time 
    st = time.time()  


        
    hodo_color = ['purple','red','darkorange','gold','#fff09f']

    #################################################################
    ### SET UP THE DATA ###
    #################################################################
    # SFC CORRECTION
    if str(type(modify_sfc)) == "<class 'dict'>":
        sounding_data = modify_surface(clean_data, modify_sfc)
    else:
        sounding_data = clean_data

    # declare easy variable names for reuse from `clean_data`
    T = sounding_data['T']
    Td = sounding_data['Td']
    p = sounding_data['p']
    z = sounding_data['z']
    u = sounding_data['u']
    v = sounding_data['v']
    wd = mpcalc.wind_direction(u, v)
    ws = mpcalc.wind_speed(u, v)

    # calculate other sounding parameters using SounderPy Calc
    general, thermo, kinem, intrp = sounding_params(sounding_data, storm_motion).calc()
    #################################################################




    #################################################################
    ### STORM-RELATIVE HODO LOGIC ###
    #################################################################
    hodo_title = 'HODOGRAPH'
    
    if sr_hodo == True:
        if ma.is_masked(kinem['sm_u']) == False:
            if np.isnan(kinem['sm_u']) == False:
                u = u - (kinem['sm_u']*units.kts)
                v = v - (kinem['sm_v']*units.kts)
                intrp['uINTRP'] = intrp['uINTRP'] - kinem['sm_u']
                intrp['vINTRP'] = intrp['vINTRP'] - kinem['sm_v']
                hodo_title = 'STORM RELATIVE HODOGRAPH'
            else:
                sr_hodo = False
                warnings.warn("This profile can not be plotted as storm relative because a storm-motion does not exist for this data"+
                              "This plot will feature a ground relative instead.", Warning)
        else:
            sr_hodo = False
            warnings.warn("This profile can not be plotted as storm relative because a storm-motion does not exist for this data"+
                  "This plot will feature a ground relative instead.", Warning)

    #################################################################
    ### DECLARE PLOT TITLES FROM CLEAN_DATA ###
    #################################################################

    top_title = f"{clean_data['titles']['top_title']} {hodo_title}"
    left_title = clean_data['titles']['left_title']
    right_title = clean_data['titles']['right_title']
    ################################################################
    
    
    
    #################################################################
    ### DEFINE HODOGRAPH BOUNDS ###
    #################################################################
    if (z.max().m - z[0].m) > 9001: 
        hodo_hgt = 9000*units.m
        # restructure u and v, p, ws and z data arrays based on corrected u and v arrays and hodo_layer depth
        p_hodo, u_hodo, v_hodo, z_hodo = mpcalc.get_layer(p, u, v, z, depth=hodo_hgt)
    else:
        p_hodo = p
        u_hodo = u
        v_hodo = v
        z_hodo = z
    # determine max height of wind data to plot on hodograph in km (if hodo_layer = 9, 0-9km u and v are plotted)
    # remove nan values from base wind u and v component arrays to find min & max values.
    u_clean = u.magnitude[np.logical_not(np.isnan(u.magnitude))]
    v_clean = v.magnitude[np.logical_not(np.isnan(v.magnitude))]
    # restructure u and v, p, ws and z data arrays based on corrected u and v arrays and hodo_layer depth
    # define x and y min/max values from 'cleaned' and restructured u and v arrays
    x_min = u_hodo.min().m
    y_min = v_hodo.min().m
    x_max = u_hodo.max().m
    y_max = v_hodo.max().m
    # if statements to determine approprate x axis and y axis limits (to change dynamically with the data)
    if y_max >= 0:
        y_Maxlimit = (y_max + 15)
    if y_max < 0:
        y_Maxlimit = (y_max + 15)

    if x_max >= 0:
        x_Maxlimit = (x_max + 15)
    if x_max < 0:
        x_Maxlimit = (x_max + 15)

    if y_min >= 0:
        y_Minlimit = (y_min - 40)
    if y_min < 0:
        y_Minlimit = (y_min - 40)

    if x_min >= 0:
        x_Minlimit = (x_min - 40)
    if x_min < 0:
        x_Minlimit = (x_min - 40)
    #################################################################
        
        
    #################################################################
    ### CREATE HODOGRAPH OBJECT ###
    #################################################################
    fig = plt.figure(figsize=(16, 12), linewidth=1.5, edgecolor=brdr_clr)
    fig.set_facecolor(bckgrnd_clr)  
    hod_ax = plt.axes((0.13, 0.11, 0.77, 0.77))
    h = Hodograph(hod_ax, component_range=150.)
    try:
        h.ax.set_xlim(x_Minlimit, x_Maxlimit)                                  
        h.ax.set_ylim(y_Minlimit, y_Maxlimit)                             
    except:
        h.ax.set_xlim(-65,65)
        h.ax.set_ylim(-65,65)
        pass                                                                         
    h.add_grid(increment=20, color=gen_txt_clr, linestyle='-', linewidth=1.5, alpha=0.4) 
    h.add_grid(increment=10, color=gen_txt_clr, linewidth=1, linestyle='--', alpha=0.4) 
    h.ax.set_facecolor(bckgrnd_clr)
    h.ax.spines["top"].set_color(brdr_clr)
    h.ax.spines["left"].set_color(brdr_clr)
    h.ax.spines["right"].set_color(brdr_clr)
    h.ax.spines["bottom"].set_color(brdr_clr)
    h.ax.spines["bottom"].set_color(brdr_clr)
    h.ax.set_box_aspect(1) 
    h.ax.set_yticklabels([])
    h.ax.set_xticklabels([])
    h.ax.set_xticks([])
    h.ax.set_yticks([])
    h.ax.set_xlabel(' ')
    h.ax.set_ylabel(' ')
    #################################################################
    
    
    
    #################################################################
    ### PLOT HEIGHT MARKERS ###
    #################################################################
    plt.xticks(np.arange(0,0,1))
    plt.yticks(np.arange(0,0,1))
    for i in range(10,130,20):
        h.ax.annotate(str(i),(i,0),xytext=(0,2),textcoords='offset pixels',clip_on=True,fontsize=12,weight='bold',alpha=0.2,zorder=0, color=gen_txt_clr)
    for i in range(10,130,20):
        h.ax.annotate(str(i),(0,i),xytext=(0,2),textcoords='offset pixels',clip_on=True,fontsize=12,weight='bold',alpha=0.2,zorder=0, color=gen_txt_clr)
    for i in range(10,130,20):
        h.ax.annotate(str(i),(-i,0),xytext=(0,2),textcoords='offset pixels',clip_on=True,fontsize=12,weight='bold',alpha=0.2,zorder=0, color=gen_txt_clr)
    for i in range(10,130,20):
        h.ax.annotate(str(i),(0,-i),xytext=(0,2),textcoords='offset pixels',clip_on=True,fontsize=12,weight='bold',alpha=0.2,zorder=0, color=gen_txt_clr)

        
        
    h.plot(intrp['uINTRP'][intrp['hgt_lvls']['h05']],intrp['vINTRP'][intrp['hgt_lvls']['h05']],marker='.', markeredgecolor='black',
           color='white', alpha=1, markersize=30, clip_on=True, zorder=5)
    h.ax.annotate(str('.5'),(intrp['uINTRP'][intrp['hgt_lvls']['h05']],intrp['vINTRP'][intrp['hgt_lvls']['h05']]),
                  weight='bold', fontsize=11, color='black',xytext=(0.02,-5),textcoords='offset pixels',horizontalalignment='center',clip_on=True, zorder=6) 

    hgt_lvls = [] 
    for key in intrp['hgt_lvls'].keys():
        hgt_lvls.append(intrp['hgt_lvls'][key])
    hgt_lvls.pop(0) 

    for lvl in hgt_lvls[1::2]:
        if lvl < 130:
            h.plot(intrp['uINTRP'][lvl],intrp['vINTRP'][lvl], marker='.', color='white', markeredgecolor='black', alpha=1, markersize=30, zorder=5)
            h.ax.annotate(str(int(round(intrp['zINTRP'][lvl]/1000,0))),(intrp['uINTRP'][lvl],intrp['vINTRP'][lvl]), 
                          weight='bold', fontsize=11, color='black',xytext=(0.02,-5),textcoords='offset pixels',horizontalalignment='center',clip_on=True, zorder=5.1) 
    #################################################################
    
    
    
    
    
    #################################################################
    ### PLOT HODOGRAPH LINE ###
    #################################################################
    hodo_color = ['purple','red','darkorange','gold','#fff09f']

    h.ax.plot(intrp['uINTRP'][0:10+1],   intrp['vINTRP'][0:10+1],   color=hodo_color[0], linewidth=7, clip_on=True)
    h.ax.plot(intrp['uINTRP'][10:30+1],  intrp['vINTRP'][10:30+1],  color=hodo_color[1], linewidth=7, clip_on=True)
    h.ax.plot(intrp['uINTRP'][30:60+1],  intrp['vINTRP'][30:60+1],  color=hodo_color[2], linewidth=7, clip_on=True)
    h.ax.plot(intrp['uINTRP'][60:90+1],  intrp['vINTRP'][60:90+1],  color=hodo_color[3], linewidth=7, clip_on=True)
    h.ax.plot(intrp['uINTRP'][90:120+1], intrp['vINTRP'][90:120+1], color=hodo_color[4], linewidth=7, clip_on=True) 

    
    
    #################################################################
    ### ADD HODOGRAPH ANNOTATION ###
    #################################################################
    if ma.is_masked(kinem['sm_rm']) == False:
        # BUNKERS STORM MOTION
        if sr_hodo == False:
            h.ax.text((kinem['sm_rm'][0]+0.5), (kinem['sm_rm'][1]-0.5), 'RM', weight='bold', ha='left', fontsize=14, alpha=0.9, color=gen_txt_clr)
            h.ax.text((kinem['sm_lm'][0]+0.5), (kinem['sm_lm'][1]-0.5), 'LM', weight='bold', ha='left', fontsize=14, alpha=0.9, color=gen_txt_clr)
            h.ax.text((kinem['sm_mw'][0]+0.5), (kinem['sm_mw'][1]-0.5), 'MW', weight='bold', ha='left', fontsize=14, alpha=0.9, color=gen_txt_clr)
        elif sr_hodo == True:
            h.ax.text((kinem['sm_lm'][0] - kinem['sm_u'] +0.5), (kinem['sm_lm'][1] - kinem['sm_v'] -0.5), 'LM', weight='bold', ha='left', fontsize=14, alpha=0.9, color=gen_txt_clr)
            h.ax.text((kinem['sm_mw'][0] - kinem['sm_u'] +0.5), (kinem['sm_mw'][1] - kinem['sm_v'] -0.5), 'MW', weight='bold', ha='left', fontsize=14, alpha=0.9, color=gen_txt_clr)
    

    if ma.is_masked(kinem['sm_u']) == False:    
        # ADD SM POINT IF ITS A CUSTOM STORM MOTION
        if str(type(storm_motion)) == "<class 'list'>":
            h.ax.text((kinem['sm_u']+0.5), (kinem['sm_v']-0.5), 'SM', weight='bold', ha='left', fontsize=14, alpha=0.9, color=gen_txt_clr)
    
        h.ax.arrow(0,0,kinem['sm_u']-0.3, kinem['sm_v']-0.3, linewidth=3, color=gen_txt_clr, alpha=0.2, 
                label='SM Vector', length_includes_head=True, head_width=0.5)
        
        if sr_hodo == False:
            # DEVIANT TORNADO MOTION
            h.ax.text(kinem['dtm'][0], (kinem['dtm'][1] + 2), 'DTM', weight='bold', fontsize=10, color='brown', ha='center')
            h.plot(kinem['dtm'][0], kinem['dtm'][1], marker='v', color='brown', markersize=8, alpha=0.8, ls='', label='DEVIANT TORNADO MOTION')
            
        elif sr_hodo == True:
            # DEVIANT TORNADO MOTION
            h.ax.text(kinem['dtm'][0] - kinem['sm_u'], (kinem['dtm'][1] - kinem['sm_v'] + 2), 'DTM', weight='bold', fontsize=10, color='brown', ha='center')
            h.plot(kinem['dtm'][0] - kinem['sm_u'], kinem['dtm'][1] - kinem['sm_v'], marker='v', color='brown', markersize=8, alpha=0.8, ls='', label='DEVIANT TORNADO MOTION')
          
        
        
    # EFFECTIVE INFLOW LAYER SRH FILL
    if sr_hodo == False:
        if ma.is_masked(kinem['eil_z'][0]) == False:

            ebot = h.ax.plot((kinem['sm_u'], intrp['uINTRP'][find_nearest(intrp['zINTRP'],kinem['eil_z'][0])]), (kinem['sm_v'], intrp['vINTRP'][find_nearest(intrp['zINTRP'],kinem['eil_z'][0])]),  
                     linestyle='--', linewidth=2.3, alpha=0.5, color='lightblue', label='Effective Inflow Layer')
            etop = h.ax.plot((kinem['sm_u'], intrp['uINTRP'][find_nearest(intrp['zINTRP'],kinem['eil_z'][1])]), (kinem['sm_v'], intrp['vINTRP'][find_nearest(intrp['zINTRP'],kinem['eil_z'][1])]),  
                     linestyle='--', linewidth=2.3, alpha=0.5, color='lightblue')
            fill_srh = h.ax.fill(np.append(intrp['uINTRP'][find_nearest(intrp['zINTRP'], kinem['eil_z'][0]):find_nearest(intrp['zINTRP'], kinem['eil_z'][1])+1], kinem['sm_u']), 
                                 np.append(intrp['vINTRP'][find_nearest(intrp['zINTRP'], kinem['eil_z'][0]):find_nearest(intrp['zINTRP'], kinem['eil_z'][1])+1], kinem['sm_v']),
                                 'lightblue',alpha=0.1, label='EIL SRH')
        else:
            fill_srh = h.ax.fill(np.append(intrp['uINTRP'][find_nearest(intrp['zINTRP'], 0):find_nearest(intrp['zINTRP'], 3000)], kinem['sm_u']), 
                         np.append(intrp['vINTRP'][find_nearest(intrp['zINTRP'], 0):find_nearest(intrp['zINTRP'], 3000)], kinem['sm_v']),
                         'lightblue',alpha=0.1, label='0-3 SRH')
        
    elif sr_hodo == True: 
        if ma.is_masked(kinem['eil_z'][0]) == False:

            ebot = h.ax.plot((0, intrp['uINTRP'][find_nearest(intrp['zINTRP'],kinem['eil_z'][0])]), (0, intrp['vINTRP'][find_nearest(intrp['zINTRP'],kinem['eil_z'][0])]),  
                     linestyle='--', linewidth=2.3, alpha=0.5, color='lightblue', label='Effective Inflow Layer')
            etop = h.ax.plot((0, intrp['uINTRP'][find_nearest(intrp['zINTRP'],kinem['eil_z'][1])]), (0, intrp['vINTRP'][find_nearest(intrp['zINTRP'],kinem['eil_z'][1])]),  
                     linestyle='--', linewidth=2.3, alpha=0.5, color='lightblue')
            fill_srh = h.ax.fill(np.append(intrp['uINTRP'][find_nearest(intrp['zINTRP'], kinem['eil_z'][0]):find_nearest(intrp['zINTRP'], kinem['eil_z'][1])+1], 0), 
                                 np.append(intrp['vINTRP'][find_nearest(intrp['zINTRP'], kinem['eil_z'][0]):find_nearest(intrp['zINTRP'], kinem['eil_z'][1])+1], 0),
                                 'lightblue',alpha=0.1, label='EIL SRH')
        else:
            fill_srh = h.ax.fill(np.append(intrp['uINTRP'][find_nearest(intrp['zINTRP'], 0):find_nearest(intrp['zINTRP'], 3000)], 0), 
                         np.append(intrp['vINTRP'][find_nearest(intrp['zINTRP'], 0):find_nearest(intrp['zINTRP'], 3000)], 0),
                         'lightblue',alpha=0.1, label='0-3 SRH')

    if sr_hodo == False:
        h.ax.text(kinem['mcs'][0], kinem['mcs'][1], 'UP', weight='bold', fontsize=12, color='orange', ha='center', alpha=0.5, clip_on=True)
        h.ax.text(kinem['mcs'][2], kinem['mcs'][3], 'DN', weight='bold', fontsize=12, color='orange', ha='center', alpha=0.5, clip_on=True)
    
    elif sr_hodo == True:
        h.ax.text(kinem['mcs'][0]- kinem['sm_u'], kinem['mcs'][1]- kinem['sm_v'], 'UP', weight='bold', fontsize=12, color='orange', ha='center', alpha=0.5, clip_on=True)
        h.ax.text(kinem['mcs'][2]- kinem['sm_u'], kinem['mcs'][3]- kinem['sm_v'], 'DN', weight='bold', fontsize=12, color='orange', ha='center', alpha=0.5, clip_on=True)
    
    
    
    # STORM MOTION PRINTOUT
    if ma.is_masked(kinem['sm_u']) == False:
        try:
            keys = ['sm_rm', 'sm_lm', 'sm_mw', 'dtm'] 

            speeds = []
            directions = []
            
            speeds.append(mpcalc.wind_speed(kinem['sm_u']*units.kts,kinem['sm_v']*units.kts).m) 
            directions.append(mpcalc.angle_to_direction((np.degrees(np.arctan2(kinem['sm_u'],kinem['sm_v']))), full=False, level=3)) 

            for key in keys:
                speeds.append(mpcalc.wind_speed(kinem[key][0]*units.kts,kinem[key][1]*units.kts).m) 
                directions.append(mpcalc.angle_to_direction((np.degrees(np.arctan2(kinem[key][0],kinem[key][1]))), full=False, level=3))

            speeds.append(mpcalc.wind_speed(kinem['mcs'][0]*units.kts,kinem['mcs'][1]*units.kts).m) 
            directions.append(mpcalc.angle_to_direction((np.degrees(np.arctan2(kinem['mcs'][0],kinem['mcs'][1]))), full=False, level=3))   

            speeds.append(mpcalc.wind_speed(kinem['mcs'][2]*units.kts,kinem['mcs'][3]*units.kts).m) 
            directions.append(mpcalc.angle_to_direction((np.degrees(np.arctan2(kinem['mcs'][2],kinem['mcs'][3]))), full=False, level=3)) 
                    
                    
            #plot Bunkers Storm Motion & DTM Data in box on Hodograph 
            plt.figtext(0.228, 0.845, 
                        f' RM: {directions[1]} @ {mag(speeds[1])} kts\n LM: {directions[2]} @ {mag(speeds[2])} kts\n' + 
                        f' MW: {directions[3]} @ {mag(speeds[3])} kts\nDTM: {directions[4]} @ {mag(speeds[4])} kts\n US: {directions[5]} @ {mag(speeds[5])} kts\n' + 
                        f' DS: {directions[6]} @ {mag(speeds[6])} kts\n',
                        color=gen_txt_clr, fontsize=10, verticalalignment='top', linespacing=2.2, alpha=0.6)  
        except IndexError:
            pass
        
        def sm_str(storm_motion, speeds, directions):
            
            if storm_motion in ['right_moving',
                                'left_moving',
                                'mean_wind']:
                return f"{str.upper(storm_motion.replace('_', ' '))} | {directions[0]} @ {mag(speeds[0])}"

            else:
                return f"USER DEFINED | {directions[0]} @ {mag(speeds[0])}"
            
        plt.figtext(0.228, 0.87, 
                f' SM: {sm_str(storm_motion, speeds, directions)} kts',
                color=gen_txt_clr, weight='bold', fontsize=10, verticalalignment='top', linespacing=2.2, alpha=0.6) 
    ################################################################
    
    
    
    
    #########################################################################
    ############################ OTHER AXES ################################# 
    #########################################################################
    
    
    
    #################################################################
    ### STREAMWISENESS AND RH W/HGT ###
    #################################################################
    # PLOT AXIS/LOC
    strmws_ax = plt.axes((0.8015, 0.11, 0.095, 0.32))
    plt.figtext(0.85, 0.40, f'SW ζ (%)', color=gen_txt_clr, weight='bold', fontsize=12, ha='center', alpha=0.7)
    strmws_ax.spines["top"].set_color(brdr_clr)
    strmws_ax.spines["left"].set_color(brdr_clr)
    strmws_ax.spines["right"].set_color(brdr_clr)
    strmws_ax.spines["bottom"].set_color(brdr_clr)
    strmws_ax.spines["bottom"].set_color(brdr_clr)   
    strmws_ax.set_facecolor(bckgrnd_clr) 

    #YTICKS
    strmws_ax.set_ylim(0, 3000)
    strmws_ax.set_yticklabels([])
    strmws_ax.set_ylabel(' ')
    strmws_ax.tick_params(axis='y', length = 0)
    strmws_ax.grid(True, axis='y')
    strmws_ax.tick_params(axis="x",direction="in", pad=-12)

    #XTICKS
    strmws_ax.set_xlim(40, 102)
    strmws_ax.set_xticks([50, 90])
    strmws_ax.set_xticklabels([50, 90], weight='bold', alpha=0.5, fontstyle='italic', color=gen_txt_clr)
    strmws_ax.set_xlabel(' ')

    #HGT LABLES 
    strmws_ax.text(47, 502 , '0.5km', fontsize=8, alpha=0.6, color=gen_txt_clr)
    strmws_ax.text(47, 1002, '1.0km', fontsize=8, alpha=0.6, color=gen_txt_clr)
    strmws_ax.text(47, 1502, '1.5km', fontsize=8, alpha=0.6, color=gen_txt_clr)
    strmws_ax.text(47, 2002, '2.0km', fontsize=8, alpha=0.6, color=gen_txt_clr)
    strmws_ax.text(47, 2502, '2.5km', fontsize=8, alpha=0.6, color=gen_txt_clr)

    if ma.is_masked(kinem['sm_u']) == False:
        plt.plot(kinem['swv_perc'][0:11],  intrp['zINTRP'][0:11],  color=hodo_color[0], lw=3, clip_on=True)
        plt.plot(kinem['swv_perc'][10:30], intrp['zINTRP'][10:30], color=hodo_color[1], lw=3, clip_on=True)
    
        if ma.is_masked(kinem['eil_z'][0]) == False:
            strmws_ax.fill_between(x=(40,102), y1=kinem['eil_z'][0], y2=kinem['eil_z'][1], color='lightblue', alpha=0.2)
    else:
        warnings.warn("Streamwiseness could not be plotted (no valid storm motion/not enough data)", Warning)
    #################################################################

    
    
    
    #################################################################
    ### VORTICITY W/HGT ###
    #################################################################
    vort_ax = plt.axes((0.8965, 0.11, 0.095, 0.32))
    plt.figtext(0.945, 0.39, f'ζₜₒₜ & ζSW\n(/sec)', color=gen_txt_clr, weight='bold', fontsize=12, ha='center', alpha=0.7)
    vort_ax.spines["top"].set_color(brdr_clr)
    vort_ax.spines["left"].set_color(brdr_clr)
    vort_ax.spines["right"].set_color(brdr_clr)
    vort_ax.spines["bottom"].set_color(brdr_clr)
    vort_ax.spines["bottom"].set_color(brdr_clr)   
    vort_ax.set_facecolor(bckgrnd_clr) 

    #YTICKS
    vort_ax.tick_params(axis='y', length = 0)
    vort_ax.grid(True, axis='y')
    vort_ax.set_ylim(0, 3000)
    vort_ax.set_yticklabels([])
    vort_ax.set_ylabel(' ')
    vort_ax.tick_params(axis="x",direction="in", pad=-12)

        
    if ma.is_masked(kinem['sm_u']) == False: 
        #XTICKS
        vort_ax.set_xlabel(' ')
        vort_max = kinem['vort'][0:30].max()+0.005
        vort_min = kinem['vort'][0:30].min()-0.005
    
        vort_ax.set_xlim(vort_min-0.002, vort_max+0.002)
        vort_ax.set_xticks([(vort_min+0.005),(vort_max-0.005)])
        vort_ax.set_xticklabels([(np.round(vort_min+0.002,2)),(np.round(vort_max-0.002,2))], weight='bold', alpha=0.5, fontstyle='italic', color=gen_txt_clr)
 
        vort_ax.plot(kinem['swv'][0:30],  intrp['zINTRP'][0:30], color='orange', linewidth=3, alpha=0.8, label='SW ζ')
        vort_ax.plot(kinem['vort'][0:30], intrp['zINTRP'][0:30], color=gen_txt_clr,  linewidth=4, alpha=0.4, label='Total ζ')
        
        if ma.is_masked(kinem['eil_z'][0]) == False:
            vort_ax.fill_between(x=(vort_min-0.002, vort_max+0.002), y1=kinem['eil_z'][0], y2=kinem['eil_z'][1], color='lightblue', alpha=0.2)
    else:
        warnings.warn("Total Vorticity could not be plotted (no valid storm motion/not enough data)", Warning)
    #################################################################
    
    
    
    
    
    #################################################################
    ### SRW W/HGT ###
    #################################################################
    wind_ax = plt.axes((0.9915, 0.11, 0.095, 0.32))
    plt.figtext(1.037, 0.39, f'SR Wind\n(kts)', weight='bold', color=gen_txt_clr, fontsize=12, ha='center', alpha=0.7)
    wind_ax.spines["top"].set_color(brdr_clr)
    wind_ax.spines["left"].set_color(brdr_clr)
    wind_ax.spines["right"].set_color(brdr_clr)
    wind_ax.spines["bottom"].set_color(brdr_clr)
    wind_ax.spines["bottom"].set_color(brdr_clr)   
    wind_ax.set_facecolor(bckgrnd_clr) 
    plt.ylabel(' ')
    plt.xlabel(' ')

    #YTICKS
    wind_ax.set_ylim(0, 3000)
    wind_ax.grid(True, axis='y')
    wind_ax.set_yticklabels([])
    wind_ax.tick_params(axis='y', length = 0)
    wind_ax.tick_params(axis="x",direction="in", pad=-12)

    if ma.is_masked(kinem['sm_u']) == False:
        #XTICKS
        wind_max = kinem['srw'][0:30].max()+1
        wind_min = kinem['srw'][0:30].min()-1
        wind_ax.set_xlim(wind_min-5, wind_max+5)
        wind_ax.set_xticks([(wind_min)+2, (wind_max)-2])
        wind_ax.set_xticklabels([(int(wind_min)+2), (int(wind_max)-2)], weight='bold', alpha=0.5, fontstyle='italic', color=gen_txt_clr)

        #PLOT SR WIND  
        wind_ax.plot(kinem['srw'][0:11],  intrp['zINTRP'][0:11],  color=hodo_color[0], clip_on=True, linewidth=3, alpha=0.8, label='0-1 SR Wind')
        wind_ax.plot(kinem['srw'][10:30], intrp['zINTRP'][10:30], color=hodo_color[1], clip_on=True, linewidth=3, alpha=0.8, label='1-3 SR Wind')
        
        if ma.is_masked(kinem['eil_z'][0]) == False:
            wind_ax.fill_between(x=(wind_min-5, wind_max+5), y1=kinem['eil_z'][0], y2=kinem['eil_z'][1], color='lightblue', alpha=0.2)
    else:
        warnings.warn("Storm Relative Wind could not be plotted (no valid storm motion/not enough data)", Warning)
        
    #################################################################
    
    
    
    
    
    #################################################################
    ### THETA & THETA E W/HGT ###
    #################################################################
    #PLOT AXES/LOC
    theta_ax = plt.axes((1.0865, 0.11, 0.095, 0.32))
    plt.figtext(1.135, 0.39, f'Theta-e &\nTheta (K)', weight='bold', color=gen_txt_clr, fontsize=12, ha='center', alpha=0.7)
    theta_ax.spines["top"].set_color(brdr_clr)
    theta_ax.spines["left"].set_color(brdr_clr)
    theta_ax.spines["right"].set_color(brdr_clr)
    theta_ax.spines["bottom"].set_color(brdr_clr)
    theta_ax.spines["bottom"].set_color(brdr_clr)   
    theta_ax.set_facecolor(bckgrnd_clr) 

    maxtheta = intrp['thetaeINTRP'][0:30].max()
    mintheta = intrp['thetaINTRP'][0:30].min()

    #YTICKS
    theta_ax.set_ylim(intrp['zINTRP'][0], 3000)
    theta_ax.set_yticklabels([])
    plt.ylabel(' ')
    theta_ax.tick_params(axis='y', length = 0)
    theta_ax.grid(True, axis='y')

    #XTICKS
    theta_ax.set_xlim(mintheta - 5, maxtheta + 5)
    theta_ax.set_xticks([(mintheta), (maxtheta)])
    theta_ax.set_xticklabels([int(mintheta), int(maxtheta)], weight='bold', alpha=0.5, fontstyle='italic', color=gen_txt_clr)

    theta_ax.tick_params(axis="x", direction="in", pad=-12)
    theta_ax.set_xlabel(' ')

    #PLOT THETA VS HGT
    plt.plot(intrp['thetaINTRP'], intrp['zINTRP'], color='purple', linewidth=3.5, alpha=0.5, clip_on=True)
    plt.plot(intrp['thetaeINTRP'], intrp['zINTRP'], color='purple', linewidth=3.5, alpha=0.8, clip_on=True)

    if ma.is_masked(kinem['eil_z'][0]) == False:
        theta_ax.fill_between(x=(mintheta - 5, maxtheta + 5), y1=kinem['eil_z'][0], y2=kinem['eil_z'][1], color='lightblue', alpha=0.2)
    #################################################################

    
    
    
    
    #########################################################################
    ############################## TEXT PLOTS ###############################
    #########################################################################

    
    
    #################################################################
    ### BOXES FOR TEXT ###
    #################################################################
    #                                 xloc   yloc   xsize  ysize
    fig.patches.extend([plt.Rectangle((0.8015, 0.43), 0.38, 0.45,
                                      edgecolor=brdr_clr, facecolor=bckgrnd_clr, linewidth=1, alpha=1,
                                      transform=fig.transFigure, figure=fig)])
    
    #################################################################
    ### KINEMATICS ###
    #################################################################
    # now some kinematic parameters
    met_per_sec = (units.m*units.m)/(units.sec*units.sec)
    plt.figtext( 0.88, 0.84, 'BS         SRH        SRW     SWζ%     SWζ', color=gen_txt_clr, weight='bold', fontsize=15, alpha=0.8)
    
    plt.figtext( 0.815, 0.81, f"0-.5ₖₘ:", weight='bold', fontsize=15, color=gen_txt_clr, alpha=0.9)
    plt.figtext( 0.875, 0.81, f"{mag(kinem['shear_0_to_500'])} kt", fontsize=15, color='deepskyblue', weight='bold')
    plt.figtext( 0.929, 0.81,  f"{mag(kinem['srh_0_to_500'])* met_per_sec:~P}", fontsize=15, color='deepskyblue', weight='bold')
    plt.figtext( 1.006,   0.81, f"{mag(kinem['srw_0_to_500'])} kt", fontsize=15, color='deepskyblue', weight='bold')
    plt.figtext( 1.070, 0.81, f"{mag(kinem['swv_perc_0_to_500'])}", fontsize=15, color='deepskyblue', weight='bold')
    plt.figtext( 1.12, 0.81,  f"{mag_round(kinem['swv_0_to_500'], 3)}", fontsize=15, color='deepskyblue', weight='bold')

    plt.figtext( 0.815, 0.78, f"0-1ₖₘ:", weight='bold', fontsize=15, color=gen_txt_clr, alpha=0.9)
    plt.figtext( 0.875, 0.78, f"{mag(kinem['shear_0_to_1000'])} kt", fontsize=15, color='mediumslateblue', weight='bold')
    plt.figtext( 0.929, 0.78,  f"{mag(kinem['srh_0_to_1000'])* met_per_sec:~P}", fontsize=15, color='mediumslateblue', weight='bold')
    plt.figtext( 1.006,   0.78, f"{mag(kinem['srw_0_to_1000'])} kt", fontsize=15, color='mediumslateblue', weight='bold')
    plt.figtext( 1.070, 0.78, f"{mag(kinem['swv_perc_0_to_1000'])}", fontsize=15, color='mediumslateblue', weight='bold')
    plt.figtext( 1.12, 0.78,  f"{mag_round(kinem['swv_0_to_1000'], 3)}", fontsize=15, color='mediumslateblue', weight='bold')

    plt.figtext( 0.815, 0.75, f"1-3ₖₘ:", weight='bold', fontsize=15, color=gen_txt_clr, alpha=0.9)
    plt.figtext( 0.875, 0.75, f"{mag(kinem['shear_1_to_3000'])} kt", fontsize=15, color='slateblue', weight='bold')
    plt.figtext( 0.929, 0.75,  f"{mag(kinem['srh_1_to_3000'])* met_per_sec:~P}", fontsize=15, color='slateblue', weight='bold')
    plt.figtext( 1.006,   0.75, f"{mag(kinem['srw_1_to_3000'])} kt", fontsize=15, color='slateblue', weight='bold')
    plt.figtext( 1.070, 0.75, f"{mag(kinem['swv_perc_1_to_3000'])}", fontsize=15, color='slateblue', weight='bold')
    plt.figtext( 1.12, 0.75,  f"{mag_round(kinem['swv_1_to_3000'], 3)}", fontsize=15, color='slateblue', weight='bold')

    plt.figtext( 0.815, 0.72, f"3-6ₖₘ:", weight='bold', fontsize=15, color=gen_txt_clr, alpha=0.9)
    plt.figtext( 0.875, 0.72, f"{mag(kinem['shear_3_to_6000'])} kt", fontsize=15, color='darkslateblue', weight='bold')
    plt.figtext( 0.929, 0.72,  f"{mag(kinem['srh_3_to_6000'])* met_per_sec:~P}", fontsize=15, color='darkslateblue', weight='bold')
    plt.figtext( 1.006, 0.72, f"{mag(kinem['srw_3_to_6000'])} kt", fontsize=15, color='darkslateblue', weight='bold')
    plt.figtext( 1.070, 0.72, f"{mag(kinem['swv_perc_3_to_6000'])}", fontsize=15, color='darkslateblue', weight='bold')
    plt.figtext( 1.12,  0.72,  f"{mag_round(kinem['swv_3_to_6000'], 3)}", fontsize=15, color='darkslateblue', weight='bold')

    plt.figtext( 0.815, 0.69,  f"6-9ₖₘ:", weight='bold', fontsize=15, color=gen_txt_clr, alpha=0.9)
    plt.figtext( 0.875, 0.69,  f"{mag(kinem['shear_6_to_9000'])} kt", fontsize=15, color='darkslateblue', weight='bold')
    plt.figtext( 0.929, 0.69,  f"{mag(kinem['srh_6_to_9000'])* met_per_sec:~P}", fontsize=15, color='darkslateblue', weight='bold')
    plt.figtext( 1.006, 0.69,  f"{mag(kinem['srw_6_to_9000'])} kt", fontsize=15, color='darkslateblue', weight='bold')
    plt.figtext( 1.070, 0.69,  f"{mag(kinem['swv_perc_6_to_9000'])}", fontsize=15, color='darkslateblue', weight='bold')
    plt.figtext( 1.12,  0.69,  f"{mag_round(kinem['swv_6_to_9000'], 3)}", fontsize=15, color='darkslateblue', weight='bold')
    #################################################################
    
    
    #################################################################
    ### THERMODYNAMICS ###
    #################################################################
    plt.figtext( 0.86, 0.65, 'SR-ECAPE       CAPE         6CAPE       3CAPE', color=gen_txt_clr, weight='bold', fontsize=15)
    #SBCAPE
    plt.figtext( 0.815, 0.62,  f"SB:", weight='bold',   fontsize=15, color=gen_txt_clr)
    plt.figtext( 0.86, 0.62,   f"{mag(thermo['sb_ecape'])} J/kg",  fontsize=15, color='orangered', weight='bold')
    plt.figtext( 0.95, 0.62,   f"{mag(thermo['sbcape'])} J/kg",   fontsize=15, color='orangered', weight='bold')
    plt.figtext( 1.035, 0.62,  f"{mag(thermo['sb6cape'])} J/kg",  fontsize=15, color='orangered', weight='bold')
    plt.figtext( 1.116, 0.62,  f"{mag(thermo['sb3cape'])} J/kg",   fontsize=15, color='orangered', weight='bold')
    #MUCAPE
    plt.figtext( 0.815, 0.59,  f"MU:", weight='bold',   fontsize=15, color=gen_txt_clr)
    plt.figtext( 0.86, 0.59,   f"{mag(thermo['mu_ecape'])} J/kg", fontsize=15, color='red', weight='bold')
    plt.figtext( 0.95, 0.59,    f"{mag(thermo['mucape'])} J/kg",  fontsize=15, color='red', weight='bold')
    plt.figtext( 1.035, 0.59,   f"{mag(thermo['mu6cape'])} J/kg", fontsize=15, color='red', weight='bold')
    plt.figtext( 1.116, 0.59,   f"{mag(thermo['mu3cape'])} J/kg", fontsize=15, color='red', weight='bold')
    #MLCAPE
    plt.figtext( 0.815, 0.56,  f'ML:', weight='bold', fontsize=15, color=gen_txt_clr)
    plt.figtext( 0.86, 0.56,   f"{mag(thermo['ml_ecape'])} J/kg", fontsize=15, color='darkred', weight='bold')
    plt.figtext( 0.95, 0.56,    f"{mag(thermo['mlcape'])} J/kg",  fontsize=15, color='darkred', weight='bold')
    plt.figtext( 1.035, 0.56,   f"{mag(thermo['ml6cape'])} J/kg", fontsize=15, color='darkred', weight='bold')
    plt.figtext( 1.116, 0.56,   f"{mag(thermo['ml3cape'])} J/kg", fontsize=15, color='darkred', weight='bold')
    # NCAPE
    plt.figtext( 0.86, 0.49, f"MUNCAPE:", weight='bold', fontsize=15, color=gen_txt_clr)
    plt.figtext( 0.943, 0.49, f" {mag(thermo['mu_ncape'])} J/kg", fontsize=15, color='brown', alpha=0.7, weight='bold')
    #DCAPE
    plt.figtext( 1.01, 0.49, f"DCAPE:", weight='bold', fontsize=15, color=gen_txt_clr)
    plt.figtext( 1.08,  0.49, f"{mag(thermo['dcape'])} J/kg", fontsize=15, color='cornflowerblue', alpha=0.7, weight='bold')
    #0-3KM LR
    plt.figtext( 0.86, 0.46, f"Γ₀₋₃:", weight='bold', fontsize=16, color=gen_txt_clr)
    plt.figtext( 0.92, 0.46, f"{mag(thermo['lr_03km'])} Δ°C/km", fontsize=15, color='saddlebrown', weight='bold')
    #3-6km LR
    plt.figtext( 1.01, 0.46, f"Γ₃₋₆:", weight='bold', fontsize=16, color=gen_txt_clr)
    plt.figtext( 1.08, 0.46, f"{mag(thermo['lr_36km'])} Δ°C/km", fontsize=15, color='saddlebrown', weight='bold')
    #################################################################
    
        
    # PLOT TITLES ----------------------------------------------------------------------------------
    plt.figtext( 0.23, 0.89, left_title, ha='left', fontsize=16, color=gen_txt_clr)
    plt.figtext( 1.20, 0.89, f'{right_title}  ', ha='right', fontsize=16, color=gen_txt_clr)
    plt.figtext( 0.23, 0.92, top_title, ha='left', weight='bold', fontsize=22, color=gen_txt_clr) 
    plt.figtext( 0.23, 0.94, ' ') 
    plt.figtext( 0.225, 0.09, f'SOUNDERPY VERTICAL PROFILE ANALYSIS TOOL | (C) KYLE J GILLETT 2024',
                ha='left', color='cornflowerblue', alpha=0.8, weight='bold', fontsize=10)
    
    img = Image.open(urlopen('https://user-images.githubusercontent.com/100786530/251580013-2e9477c9-e36a-4163-accb-fe46780058dd.png'))
    #                  side-side  up-down  size   size
    imgax = fig.add_axes([0.22, 0.12, 0.08, 0.08], anchor='SE')
    imgax.imshow(img)
    imgax.axis('off')
    
    plt.tight_layout()
    
    elapsed_time = time.time() - st
    print('> RUNTIME:', time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
    
    return plt










#########################################################################
######################## COMPOSITE SOUNDING #############################

def __composite_sounding(data_list, shade_between, cmap, colors_to_use,
                         ls_to_use, alphas_to_use, lw_to_use, dark_mode): 
    
    
    ###############################################
    ### DETERMINE WHAT COLORS TO USE
    ###############################################
    # may be a user defined colormap, a user defined list of colors, or the default colormap (viridis)
    
    def get_colors(data_list, cmap):
        '''
        convert a linear colormap to a list of colors
        '''
        # get the cmap
        cmap = plt.get_cmap(cmap)
        # get value slices based on data list size
        values = np.linspace(0, 1, len(data_list))
        # Get colors from the colormap
        colors = [cmap(value) for value in values]
        # return the hex values in a list
        return [mcolors.rgb2hex(color) for color in colors]
    
    
    if (colors_to_use != 'none'):
        # user has provided a list of colors
        colors_to_use = colors_to_use
    else:
        # user has provided a colormap
        # or the default cmap will be used
        colors_to_use = get_colors(data_list, cmap)        
    ########################################################
    
    
    #
    if alphas_to_use == 'none':
        alphas_to_use = []
        for i in range(0, len(data_list)):
            alphas_to_use.append(1)
    
    if ls_to_use == 'none':
        ls_to_use = []
        for i in range(0, len(data_list)):
            ls_to_use.append('-')   
  
    if lw_to_use == 'none':
        lw_to_use = []
        for i in range(0, len(data_list)):
            lw_to_use.append(3)
    
    # Create an empty list to store letters
    letter_list = []
    # Loop through the data and append letters to the list
    for i in range(0, len(data_list)):
        # Append each letter to the list
        letter_list.append(str.upper(chr(ord('a')+i)))
    
    
    # interpolate function to intrp u and v data for clean plotting w/ hgt
    def interpolate(var,hgts,step):
            levels=np.arange(0,np.max(hgts),step)
            varinterp=np.zeros(len(levels))
            for i in range(0,len(levels)):
                lower=np.where(hgts-levels[i]<=0,hgts-levels[i],-np.inf).argmax()
                varinterp[i]=(((var[lower+1]-var[lower])/(hgts[lower+1]-hgts[lower]))*(levels[i]-hgts[lower])+var[lower])
            return varinterp 
    resolution=100

    
    ########################################
    ### PROFILE INTERPOLATION
    ########################################
    list_of_dicts = []
    for i in range(0, len(data_list)):
            list_of_dicts.append({})  
    for empty_dict, profile in zip(list_of_dicts, data_list):           
        intrp_p  = profile['p'].m[np.isnan(profile['z'])==False]
        intrp_u  = profile['u'].m[np.isnan(profile['z'])==False]
        intrp_v  = profile['v'].m[np.isnan(profile['z'])==False]
        intrp_z  = profile['z'].m[np.isnan(profile['z'])==False]
        intrp_z  = intrp_z-intrp_z[0]
        intrp_p  = interpolate(intrp_p, intrp_z, resolution)
        intrp_u  = interpolate(intrp_u, intrp_z, resolution)
        intrp_v  = interpolate(intrp_v, intrp_z, resolution)
        intrp_z  = interpolate(intrp_z, intrp_z, resolution)
        empty_dict['intrp_p'] = intrp_p
        empty_dict['intrp_u'] = intrp_u
        empty_dict['intrp_v'] = intrp_v
        empty_dict['intrp_z'] = intrp_z
    ########################################
        
        
    # Define light-mode and dark-mode properties
    if dark_mode == True:
        gen_txt_clr = 'white'
        bckgrnd_clr = 'black'
        brdr_clr    = 'white'
        barb_clr    = 'white'
        shade_alpha = 0.06
    else: 
        gen_txt_clr = 'black'
        bckgrnd_clr = 'white'
        brdr_clr    = 'black'
        barb_clr    = 'black'
        shade_alpha = 0.02
    
    
    # record process time 
    st = time.time()  

        
    
    #########################################################################
    ################################ SKEW-T ################################# 
    #########################################################################
    
    
    #################################################################
    ### CREATE FIGURE ###
    #################################################################
    # Define figure and skew-t
    fig = plt.figure(figsize=(22,13), linewidth=10, edgecolor=brdr_clr)      
    skew = SkewT(fig, rotation=47, rect=(0.1124, 0.1005, 0.60, 0.85))  
    skew.ax.set_box_aspect(0.87)
    skew.ax.zorder = 5
    # Define axis bounds 
    skew.ax.set_adjustable('datalim')
    skew.ax.set_ylim(1050, 100)    
    skew.ax.set_xlim(-45, 52)   
    # Define axis labels 
    plt.xlabel("  ", fontsize=12)
    plt.ylabel("  ", fontsize=12) 
    plt.xticks(fontsize=13)  
    plt.yticks(fontsize=13, ha='left')
    plt.tick_params(axis="x",direction="in", pad=-12, colors=gen_txt_clr)
    plt.tick_params(axis="y",direction="in", pad=-7, colors=gen_txt_clr)
    skew.ax.set_yticks([1000, 900, 800, 700, 600, 500, 400, 300, 200])
    skew.ax.set_yticklabels([1000, 900, 800, 700, 600, 500, 400, 300, 200], color=gen_txt_clr)
    skew.ax.spines["top"].set_color(brdr_clr)
    skew.ax.spines["left"].set_color(brdr_clr)
    skew.ax.spines["right"].set_color(brdr_clr)
    skew.ax.spines["bottom"].set_color(brdr_clr)
    skew.ax.spines["bottom"].set_color(brdr_clr)
    # Define background colors
    fig.set_facecolor(bckgrnd_clr)         
    skew.ax.set_facecolor(bckgrnd_clr)    
    # Add shaded isotherms
    x1 = np.linspace(-100, 40, 8)                                                          
    x2 = np.linspace(-90, 50, 8)                                                         
    y = [1200, 50]                                                                      
    for i in range(0, 8):              
        skew.shade_area(y=y, x1=x1[i], x2=x2[i], color='gray', alpha=shade_alpha, zorder=1)   
    #################################################################
    
    
    
    #################################################################
    ### PLOT SKEW T LINES ###
    #################################################################
    # Plot relevent Skew-T lines
    skew.ax.axvline(0 * units.degC, linestyle='--', color='blue', alpha=0.3)
    skew.ax.axvline(-20 * units.degC, linestyle='--', color='blue', alpha=0.3)
    skew.plot_dry_adiabats(color='black', linewidth=0.5, alpha=0.4)   
    skew.plot_moist_adiabats(linewidth=0.5, alpha=0.4)
    skew.plot_mixing_lines(linewidth=0.2, alpha=0.4)
    # add basic temperature lines
    
    for profile, i in zip(data_list, range(0, len(data_list))):
        skew.plot(profile['p'], profile['Td'], colors_to_use[i],
                  alpha=alphas_to_use[i], ls=ls_to_use[i], lw=lw_to_use[i])
        skew.plot(profile['p'], profile['T'], colors_to_use[i],
                  alpha=alphas_to_use[i], ls=ls_to_use[i], lw=lw_to_use[i],
                  label=f'PROFILE {letter_list[i]}')                 
        if shade_between == True:
            skew.shade_area(profile['p'], x1=profile['Td'], x2=profile['T'],
                            color=colors_to_use[i], alpha=0.05, zorder=1) 

    
    
    
    #################################################################
    ### CREATE HODOGRAPH OBJECT ###
    #################################################################
    hod_ax = plt.axes((0.53, 0.1003, 0.85, 0.85))
    h = Hodograph(hod_ax, component_range=150.)
    h.ax.set_xlim(-80,80)
    h.ax.set_ylim(-80,80)
    h.add_grid(increment=20, color=gen_txt_clr, linestyle='-', linewidth=1.5, alpha=0.4) 
    h.add_grid(increment=10, color=gen_txt_clr, linewidth=1, linestyle='--', alpha=0.4)
    h.ax.set_facecolor(bckgrnd_clr)
    h.ax.spines["top"].set_color(brdr_clr)
    h.ax.spines["left"].set_color(brdr_clr)
    h.ax.spines["right"].set_color(brdr_clr)
    h.ax.spines["bottom"].set_color(brdr_clr)
    h.ax.spines["bottom"].set_color(brdr_clr)
    h.ax.set_box_aspect(1) 
    h.ax.set_yticklabels([])
    h.ax.set_xticklabels([])
    h.ax.set_xticks([])
    h.ax.set_yticks([])
    h.ax.set_xlabel(' ')
    h.ax.set_ylabel(' ')
    
    # ADD HODOGRAPH AXIS MARKERS
    plt.xticks(np.arange(0,0,1))
    plt.yticks(np.arange(0,0,1))
    for i in range(10,200,10):
        h.ax.annotate(str(i),(i,0),xytext=(0,2),textcoords='offset pixels',clip_on=True,
                      color=gen_txt_clr, fontsize=12,weight='bold',alpha=0.5,zorder=0)
    for i in range(10,200,10):
        h.ax.annotate(str(i),(0,i),xytext=(0,2),textcoords='offset pixels',clip_on=True,
                      color=gen_txt_clr, fontsize=12,weight='bold',alpha=0.5,zorder=0)
    for i in range(10,200,10):
        h.ax.annotate(str(i),(-i,0),xytext=(0,2),textcoords='offset pixels',clip_on=True,
                      color=gen_txt_clr, fontsize=12,weight='bold',alpha=0.5,zorder=0)
    for i in range(10,200,10):
        h.ax.annotate(str(i),(0,-i),xytext=(0,2),textcoords='offset pixels',clip_on=True,
                      color=gen_txt_clr, fontsize=12,weight='bold',alpha=0.5,zorder=0)

    
    # PLOT EACH HODOGRAPH LINE FROM EACH PROFILE, DECREASING ALPHA = INCREASING HGHT
    # LINEWIDTHS CAN'T BE CHANGED HERE AS LARGER LINEWIDTHS LOOK BEST ON THE HODOGRAPH
    for profile, i in zip(list_of_dicts, range(0, len(data_list))):
        h.ax.plot(profile['intrp_u'][0:11], profile['intrp_v'][0:11], 
                  color=colors_to_use[i], alpha=1, linewidth=5, clip_on=True)
        h.ax.plot(profile['intrp_u'][10:31], profile['intrp_v'][10:31], 
                  color=colors_to_use[i], alpha=0.8, linewidth=5, clip_on=True)
        h.ax.plot(profile['intrp_u'][30:61], profile['intrp_v'][30:61], 
                  color=colors_to_use[i], alpha=0.6, linewidth=5, clip_on=True)
        h.ax.plot(profile['intrp_u'][60:91], profile['intrp_v'][60:91], 
                  color=colors_to_use[i], alpha=0.4, linewidth=5, clip_on=True)

    # CREATE MANUAL ALPHA-HGT LEGEND
    plt.figtext( 0.72, 0.12, '0-1km', color=gen_txt_clr, weight='bold', fontsize=20, alpha=1)
    plt.figtext( 0.77, 0.12, '1-3km', color=gen_txt_clr, weight='bold', fontsize=20, alpha=0.8)
    plt.figtext( 0.82, 0.12, '3-6km', color=gen_txt_clr, weight='bold', fontsize=20, alpha=0.6)
    plt.figtext( 0.87, 0.12, '6-9km', color=gen_txt_clr, weight='bold', fontsize=20, alpha=0.4)
    
    
    #########################################################################
    ############################# PLOT EXTRAS ############################### 
    ######################################################################### 
    
    # LOOP THROUGH THE PROFILES AND DYNAMICALLY ADD THEM TO THE BOTTOM OF THE PLOT
    list_of_dicts = []
    for i in range(0, len(data_list)):
            list_of_dicts.append({}) 

    locs      = []
    types     = []
    val_times = []
    int_times = []

    for profile in data_list:

        if profile['site_info']['source'] == 'RAOB OBSERVED PROFILE':
            types.append(f"RAOB")
            locs.append(f"{profile['site_info']['site-id']}, {profile['site_info']['site-name']}")
            val_times.append(f"{profile['site_info']['valid-time'][3]}Z | {profile['site_info']['valid-time'][1]}-{profile['site_info']['valid-time'][2]}-{profile['site_info']['valid-time'][0]}")
            int_times.append(f"{profile['site_info']['valid-time'][3]}Z | {profile['site_info']['valid-time'][1]}-{profile['site_info']['valid-time'][2]}-{profile['site_info']['valid-time'][0]}")

        if profile['site_info']['source'] == 'MODEL REANALYSIS':
            types.append(f"{profile['site_info']['fcst-hour']} {profile['site_info']['model']}")
            locs.append(f"{profile['site_info']['site-latlon']}")
            val_times.append(f"{profile['site_info']['valid-time'][3]}Z | {profile['site_info']['valid-time'][1]}-{profile['site_info']['valid-time'][2]}-{profile['site_info']['valid-time'][0]}")
            int_times.append(f"{profile['site_info']['run-time'][3]}Z | {profile['site_info']['run-time'][1]}-{profile['site_info']['run-time'][2]}-{profile['site_info']['run-time'][0]}")

        if profile['site_info']['source'] == 'BUFKIT FORECAST PROFILE':
            types.append(f"{profile['site_info']['fcst-hour']} {profile['site_info']['model']}")
            locs.append(f"{profile['site_info']['site-id']}, {profile['site_info']['site-name']}")
            val_times.append(f"{profile['site_info']['valid-time'][3]}Z | {profile['site_info']['valid-time'][1]}-{profile['site_info']['valid-time'][2]}-{profile['site_info']['valid-time'][0]}")
            int_times.append(f"{profile['site_info']['run-time'][3]}Z | {profile['site_info']['run-time'][1]}-{profile['site_info']['run-time'][2]}-{profile['site_info']['run-time'][0]}")
        
        if profile['site_info']['source'] == 'ACARS OBSERVED AIRCRAFT PROFILE':
            types.append(f"ACARS OB")
            locs.append(f"{profile['site_info']['site-id']}, {profile['site_info']['site-name']}")
            val_times.append(f"{profile['site_info']['valid-time'][3]}Z | {profile['site_info']['valid-time'][1]}-{profile['site_info']['valid-time'][2]}-{profile['site_info']['valid-time'][0]}")
            int_times.append(f"{profile['site_info']['valid-time'][3]}Z | {profile['site_info']['valid-time'][1]}-{profile['site_info']['valid-time'][2]}-{profile['site_info']['valid-time'][0]}")
        
        
    top_title = 'COMPOSITE VERTICAL PROFILE'
    
    
    plt.figtext( 0.19, 0.975, top_title, weight='bold', ha='left', fontsize=44, color=gen_txt_clr)
    plt.figtext( 0.20, 1.02, ' ', ha='left', fontsize=12)
    skew.ax.legend(fontsize=20, loc='upper right', framealpha=0.3, labelcolor=gen_txt_clr)
    
    plt.figtext( 0.15, 0.06, 'DATA TYPE', 
                color=gen_txt_clr, ha='left', weight='bold', fontsize=25, alpha=1)
    
    for txt, i in zip(types, range(0, len(types))):
            plt.figtext( 0.15, 0.03-(i/34), f'{letter_list[i]} | {txt}', 
                        color=colors_to_use[i], ha='left', weight='bold', fontsize=20, alpha=1)
    
    plt.figtext( 0.32, 0.06, 'LOCATIONS', 
                color=gen_txt_clr, ha='left', weight='bold', fontsize=25, alpha=1)
    
    for txt, i in zip(locs, range(0, len(locs))):
            plt.figtext( 0.32, 0.03-(i/34), f'{letter_list[i]} | {txt}', 
                        color=colors_to_use[i], ha='left', weight='bold', fontsize=20, alpha=1)
                        
    plt.figtext( 0.77, 0.06, 'RUN/LAUNCH TIME', 
                color=gen_txt_clr, ha='left', weight='bold', fontsize=25, alpha=1)
    
    for txt, i in zip(int_times, range(0, len(int_times))):
            plt.figtext( 0.77, 0.03-(i/34), f'{letter_list[i]} | {txt}', 
                        color=colors_to_use[i], ha='left', weight='bold', fontsize=20, alpha=1)

    plt.figtext( 0.97, 0.06, 'VALID TIME', 
                color=gen_txt_clr, ha='left', weight='bold', fontsize=25, alpha=1)
    
    for txt, i in zip(val_times, range(0, len(val_times))):
            plt.figtext( 0.97, 0.03-(i/34), f'{letter_list[i]} | {txt}', 
                        color=colors_to_use[i], ha='left', weight='bold', fontsize=20, alpha=1)
    plt.figtext( 0.75, 0.03-((i+1)/34), f'  ')      
    
    
    # plot author credit information 
    plt.figtext( 1, 0.96, 'SOUNDERPY (C) KYLE J GILLETT 2023        ', 
                fontsize=15, ha='left', color='cornflowerblue', weight='bold', alpha=0.6)

    # sounderpy logo
    img = Image.open(urlopen('https://user-images.githubusercontent.com/100786530/251580013-2e9477c9-e36a-4163-accb-fe46780058dd.png'))
    #                  side-side  up-down  size   size
    imgax = fig.add_axes([0.115, 0.955, 0.07, 0.07], anchor='SE', zorder=3)
    imgax.imshow(img)
    imgax.axis('off')
    
    plt.tight_layout()
    
    print('> COMPLETE --------')
    elapsed_time = time.time() - st
    print('> RUNTIME:', time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
    
    return plt













#########################################################################
########################### VAD HODOGRAPH ###############################

def __vad_hodograph(vad_data, dark_mode, storm_motion, sr_hodo):


    if dark_mode == True:
        gen_txt_clr = 'white'
        bckgrnd_clr = 'black'
        brdr_clr    = 'white'
        barb_clr    = 'white'
        shade_alpha = 0.06
    else: 
        gen_txt_clr = 'black'
        bckgrnd_clr = 'white'
        brdr_clr    = 'black'
        barb_clr    = 'black'
        shade_alpha = 0.02
    
    
    # record process time 
    st = time.time()  

    
    def find_nearest(array, value):
        array = np.asarray(array)
        nearest_idx = (np.abs(array - value)).argmin()
        return nearest_idx
        
    hodo_color = ['purple','red','darkorange','gold','#fff09f'] 
    
    def mag(param):
        if ma.is_masked(param):
            fixed = '--'
        else:
            try:
                fixed = int(param.m)
            except:
                try: 
                    fixed = int(param)
                except: fixed = param
        return fixed
    
    def mag_round(param, dec):
        if ma.is_masked(param):
            fixed = '--'
        else:
            fixed = np.round(param, dec)
        return fixed

    #################################################################
    ### SET UP THE DATA ###
    #################################################################
    # declare easy variable names for reuse from `clean_data` 
    z  = vad_data['z']
    u  = vad_data['u']
    v  = vad_data['v']
    
    # calculate other sounding parameters using SounderPy Calc
    kinem, intrp = vad_params(vad_data, storm_motion).calc()
    #################################################################
    
    hodo_title = 'HODOGRAPH'
    
    if sr_hodo == True:
        if ma.is_masked(kinem['sm_u']) == False:
            if np.isnan(kinem['sm_u']) == False:
                u = u - (kinem['sm_u'])
                v = v - (kinem['sm_v'])
                intrp['uINTRP'] = intrp['uINTRP'] - kinem['sm_u']
                intrp['vINTRP'] = intrp['vINTRP'] - kinem['sm_v']
                hodo_title = 'STORM RELATIVE HODOGRAPH'
            else:
                sr_hodo = False
                warnings.warn("This profile can not be plotted as storm relative because a storm-motion does not exist for this data"+
                              "This plot will feature a ground relative instead.", Warning)
        else:
            sr_hodo = False
            warnings.warn("This profile can not be plotted as storm relative because a storm-motion does not exist for this data"+
                  "This plot will feature a ground relative instead.", Warning)
    
    #################################################################
    ### DETERMINE PLOT TITLE BASED ON THE DATA ###
    #################################################################
    top_title = f"{vad_data['site_info']['source']} {hodo_title}"
    left_title = f"VALID: {vad_data['site_info']['valid-time'][1]}-{vad_data['site_info']['valid-time'][2]}-{vad_data['site_info']['valid-time'][0]} {vad_data['site_info']['valid-time'][3]}Z"
    right_title = f"{vad_data['site_info']['site-id']} - {vad_data['site_info']['site-name']} | {vad_data['site_info']['site-latlon'][0]}, {vad_data['site_info']['site-latlon'][1]}    " 

    ################################################################    
    
    
    
    #################################################################
    ### DEFINE HODOGRAPH BOUNDS ###
    #################################################################
    # determine max height of wind data to plot on hodograph in km (if hodo_layer = 9, 0-9km u and v are plotted)
    # remove nan values from base wind u and v component arrays to find min & max values.
    u_clean = u[np.logical_not(np.isnan(u))]
    v_clean = v[np.logical_not(np.isnan(v))]
    # restructure u and v, p, ws and z data arrays based on corrected u and v arrays and hodo_layer depth
    # define x and y min/max values from 'cleaned' and restructured u and v arrays
    x_min = u.min()
    y_min = v.min()
    x_max = u.max()
    y_max = v.max()
    # if statements to determine approprate x axis and y axis limits (to change dynamically with the data)
    if y_max >= 0:
        y_Maxlimit = (y_max + 15)
    if y_max < 0:
        y_Maxlimit = (y_max + 15)

    if x_max >= 0:
        x_Maxlimit = (x_max + 15)
    if x_max < 0:
        x_Maxlimit = (x_max + 15)

    if y_min >= 0:
        y_Minlimit = (y_min - 40)
    if y_min < 0:
        y_Minlimit = (y_min - 40)

    if x_min >= 0:
        x_Minlimit = (x_min - 40)
    if x_min < 0:
        x_Minlimit = (x_min - 40)
    #################################################################
        
        
    #################################################################
    ### CREATE HODOGRAPH OBJECT ###
    #################################################################
    fig = plt.figure(figsize=(16, 12), linewidth=1.5, edgecolor=brdr_clr)
    fig.set_facecolor(bckgrnd_clr)  
    
    hod_ax = plt.axes((0.13, 0.11, 0.77, 0.77))
    h = Hodograph(hod_ax, component_range=150.)
    try:
        h.ax.set_xlim(x_Minlimit, x_Maxlimit)                                  
        h.ax.set_ylim(y_Minlimit, y_Maxlimit)                             
    except:
        h.ax.set_xlim(-65,65)
        h.ax.set_ylim(-65,65)
        pass                                                                         
    h.add_grid(increment=20, color=gen_txt_clr, linestyle='-', linewidth=1.5, alpha=0.4) 
    h.add_grid(increment=10, color=gen_txt_clr, linewidth=1, linestyle='--', alpha=0.4) 
    h.ax.set_facecolor(bckgrnd_clr)
    h.ax.spines["top"].set_color(brdr_clr)
    h.ax.spines["left"].set_color(brdr_clr)
    h.ax.spines["right"].set_color(brdr_clr)
    h.ax.spines["bottom"].set_color(brdr_clr)
    h.ax.spines["bottom"].set_color(brdr_clr)
    h.ax.set_box_aspect(1) 
    h.ax.set_yticklabels([])
    h.ax.set_xticklabels([])
    h.ax.set_xticks([])
    h.ax.set_yticks([])
    h.ax.set_xlabel(' ')
    h.ax.set_ylabel(' ')
    #################################################################
    
    
    
    #################################################################
    ### PLOT HEIGHT MARKERS ###
    #################################################################
    plt.xticks(np.arange(0,0,1))
    plt.yticks(np.arange(0,0,1))
    for i in range(10,130,20):
        h.ax.annotate(str(i),(i,0),xytext=(0,2),textcoords='offset pixels',clip_on=True,fontsize=12,weight='bold',alpha=0.2,zorder=0, color=gen_txt_clr)
    for i in range(10,130,20):
        h.ax.annotate(str(i),(0,i),xytext=(0,2),textcoords='offset pixels',clip_on=True,fontsize=12,weight='bold',alpha=0.2,zorder=0, color=gen_txt_clr)
    for i in range(10,130,20):
        h.ax.annotate(str(i),(-i,0),xytext=(0,2),textcoords='offset pixels',clip_on=True,fontsize=12,weight='bold',alpha=0.2,zorder=0, color=gen_txt_clr)
    for i in range(10,130,20):
        h.ax.annotate(str(i),(0,-i),xytext=(0,2),textcoords='offset pixels',clip_on=True,fontsize=12,weight='bold',alpha=0.2,zorder=0, color=gen_txt_clr)

    h.plot(intrp['uINTRP'][intrp['hgt_lvls']['h05']],intrp['vINTRP'][intrp['hgt_lvls']['h05']],marker='.', markeredgecolor='black',
           color='white', alpha=1, markersize=30, clip_on=True, zorder=5)
    h.ax.annotate(str('.5'),(intrp['uINTRP'][intrp['hgt_lvls']['h05']],intrp['vINTRP'][intrp['hgt_lvls']['h05']]),
                  weight='bold', fontsize=11, color='black',xytext=(0.02,-5),textcoords='offset pixels',horizontalalignment='center',clip_on=True, zorder=6) 

    hgt_lvls = [] 
    for key in intrp['hgt_lvls'].keys():
        hgt_lvls.append(intrp['hgt_lvls'][key])
    hgt_lvls.pop(0) 

    for lvl in hgt_lvls[1::2]:
        if lvl < 130:
            h.plot(intrp['uINTRP'][lvl],intrp['vINTRP'][lvl], marker='.', color='white', markeredgecolor='black', alpha=1, markersize=30, zorder=5)
            h.ax.annotate(str(int(round(intrp['zINTRP'][lvl]/1000,0))),(intrp['uINTRP'][lvl],intrp['vINTRP'][lvl]), 
                          weight='bold', fontsize=11, color='black',xytext=(0.02,-5),textcoords='offset pixels',horizontalalignment='center',clip_on=True, zorder=5.1) 
    #################################################################
    
    
    
    
    
    #################################################################
    ### PLOT HODOGRAPH LINE ###
    #################################################################
    hodo_color = ['purple','red','darkorange','gold','#fff09f']

    h.ax.plot(intrp['uINTRP'][0:10+1],   intrp['vINTRP'][0:10+1],   color=hodo_color[0], linewidth=6, clip_on=True)
    h.ax.plot(intrp['uINTRP'][10:30+1],  intrp['vINTRP'][10:30+1],  color=hodo_color[1], linewidth=6, clip_on=True)
    h.ax.plot(intrp['uINTRP'][30:60+1],  intrp['vINTRP'][30:60+1],  color=hodo_color[2], linewidth=6, clip_on=True)
    h.ax.plot(intrp['uINTRP'][60:90+1],  intrp['vINTRP'][60:90+1],  color=hodo_color[3], linewidth=6, clip_on=True)
    h.ax.plot(intrp['uINTRP'][90:120+1], intrp['vINTRP'][90:120+1], color=hodo_color[4], linewidth=6, clip_on=True) 

    
    #################################################################
    ### ADD HODOGRAPH ANNOTATION ###
    #################################################################
    if ma.is_masked(kinem['sm_rm']) == False:
        # BUNKERS STORM MOTION
        if sr_hodo == False:
            h.ax.text((kinem['sm_rm'][0]+0.5), (kinem['sm_rm'][1]-0.5), 'RM', weight='bold', ha='left', fontsize=14, alpha=0.9, color=gen_txt_clr)
            h.ax.text((kinem['sm_lm'][0]+0.5), (kinem['sm_lm'][1]-0.5), 'LM', weight='bold', ha='left', fontsize=14, alpha=0.9, color=gen_txt_clr)
            h.ax.text((kinem['sm_mw'][0]+0.5), (kinem['sm_mw'][1]-0.5), 'MW', weight='bold', ha='left', fontsize=14, alpha=0.9, color=gen_txt_clr)
        elif sr_hodo == True:
            h.ax.text((kinem['sm_lm'][0] - kinem['sm_u'] +0.5), (kinem['sm_lm'][1] - kinem['sm_v'] -0.5), 'LM', weight='bold', ha='left', fontsize=14, alpha=0.9, color=gen_txt_clr)
            h.ax.text((kinem['sm_mw'][0] - kinem['sm_u'] +0.5), (kinem['sm_mw'][1] - kinem['sm_v'] -0.5), 'MW', weight='bold', ha='left', fontsize=14, alpha=0.9, color=gen_txt_clr)
    

    if ma.is_masked(kinem['sm_u']) == False:    
        # ADD SM POINT IF ITS A CUSTOM STORM MOTION
        if str(type(storm_motion)) == "<class 'list'>":
            h.ax.text((kinem['sm_u']+0.5), (kinem['sm_v']-0.5), 'SM', weight='bold', ha='left', fontsize=14, alpha=0.9, color=gen_txt_clr)
    
        h.ax.arrow(0,0,kinem['sm_u']-0.3, kinem['sm_v']-0.3, linewidth=3, color=gen_txt_clr, alpha=0.2, 
                label='SM Vector', length_includes_head=True, head_width=0.5)
        
        if sr_hodo == False:
            # DEVIANT TORNADO MOTION
            h.ax.text(kinem['dtm'][0], (kinem['dtm'][1] + 2), 'DTM', weight='bold', fontsize=10, color='brown', ha='center')
            h.plot(kinem['dtm'][0], kinem['dtm'][1], marker='v', color='brown', markersize=8, alpha=0.8, ls='', label='DEVIANT TORNADO MOTION')
            
        elif sr_hodo == True:
            # DEVIANT TORNADO MOTION
            h.ax.text(kinem['dtm'][0] - kinem['sm_u'], (kinem['dtm'][1] - kinem['sm_v'] + 2), 'DTM', weight='bold', fontsize=10, color='brown', ha='center')
            h.plot(kinem['dtm'][0] - kinem['sm_u'], kinem['dtm'][1] - kinem['sm_v'], marker='v', color='brown', markersize=8, alpha=0.8, ls='', label='DEVIANT TORNADO MOTION')
          
        
        
    # EFFECTIVE INFLOW LAYER
    if sr_hodo == False:
        ebot = h.ax.plot((kinem['sm_u'], intrp['uINTRP'][0]), (kinem['sm_v'], intrp['vINTRP'][0]),  
                     linestyle='--', linewidth=2.3, alpha=0.5, color='lightblue', label='Effective Inflow Layer')
        etop = h.ax.plot((kinem['sm_u'], intrp['uINTRP'][find_nearest(z,3000)]), (kinem['sm_v'], intrp['vINTRP'][find_nearest(z, 3000)]),  
                     linestyle='--', linewidth=2.3, alpha=0.5, color='lightblue')
        # EFFECTIVE INFLOW LAYER SRH FILL
        fill_srh = h.ax.fill(np.append(intrp['uINTRP'][0:find_nearest(z, 3000)+1], kinem['sm_u']), 
                             np.append(intrp['vINTRP'][0:find_nearest(z, 3000)+1], kinem['sm_v']),
                             'lightblue',alpha=0.1, label='EIL SRH')
        
    elif sr_hodo == True: 
        ebot = h.ax.plot((0, intrp['uINTRP'][0]), (0, intrp['vINTRP'][0]),  
                     linestyle='--', linewidth=2.3, alpha=0.5, color='lightblue', label='Effective Inflow Layer')
        etop = h.ax.plot((0, intrp['uINTRP'][find_nearest(z, 3000)]), (0, intrp['vINTRP'][find_nearest(z, 3000)]),  
                     linestyle='--', linewidth=2.3, alpha=0.5, color='lightblue')
        # EFFECTIVE INFLOW LAYER SRH FILL
        fill_srh = h.ax.fill(np.append(intrp['uINTRP'][0:find_nearest(z, 3000)+1], 0), 
                             np.append(intrp['vINTRP'][0:find_nearest(z, 3000)+1], 0),
                             'lightblue',alpha=0.1, label='EIL SRH')
    
    if sr_hodo == False:
        h.ax.text(kinem['mcs'][0], kinem['mcs'][1], 'UP', weight='bold', fontsize=12, color='orange', ha='center', alpha=0.5, clip_on=True)
        h.ax.text(kinem['mcs'][2], kinem['mcs'][3], 'DN', weight='bold', fontsize=12, color='orange', ha='center', alpha=0.5, clip_on=True)
    
    elif sr_hodo == True:
        h.ax.text(kinem['mcs'][0]- kinem['sm_u'], kinem['mcs'][1]- kinem['sm_v'], 'UP', weight='bold', fontsize=12, color='orange', ha='center', alpha=0.5, clip_on=True)
        h.ax.text(kinem['mcs'][2]- kinem['sm_u'], kinem['mcs'][3]- kinem['sm_v'], 'DN', weight='bold', fontsize=12, color='orange', ha='center', alpha=0.5, clip_on=True)
    
    
    
    # STORM MOTION PRINTOUT
    if ma.is_masked(kinem['sm_u']) == False:
        try:
            keys = ['sm_rm', 'sm_lm', 'sm_mw', 'dtm'] 

            speeds = []
            directions = []
            
            speeds.append(mpcalc.wind_speed(kinem['sm_u']*units.kts,kinem['sm_v']*units.kts).m) 
            directions.append(mpcalc.angle_to_direction((np.degrees(np.arctan2(kinem['sm_u'],kinem['sm_v']))), full=False, level=3)) 

            for key in keys:
                speeds.append(mpcalc.wind_speed(kinem[key][0]*units.kts,kinem[key][1]*units.kts).m) 
                directions.append(mpcalc.angle_to_direction((np.degrees(np.arctan2(kinem[key][0],kinem[key][1]))), full=False, level=3))

            speeds.append(mpcalc.wind_speed(kinem['mcs'][0]*units.kts,kinem['mcs'][1]*units.kts).m) 
            directions.append(mpcalc.angle_to_direction((np.degrees(np.arctan2(kinem['mcs'][0],kinem['mcs'][1]))), full=False, level=3))   

            speeds.append(mpcalc.wind_speed(kinem['mcs'][2]*units.kts,kinem['mcs'][3]*units.kts).m) 
            directions.append(mpcalc.angle_to_direction((np.degrees(np.arctan2(kinem['mcs'][2],kinem['mcs'][3]))), full=False, level=3)) 
                    
                    
            #plot Bunkers Storm Motion & DTM Data in box on Hodograph 
            plt.figtext(0.228, 0.845, 
                        f' RM: {directions[1]} @ {mag(speeds[1])} kts\n LM: {directions[2]} @ {mag(speeds[2])} kts\n' + 
                        f' MW: {directions[3]} @ {mag(speeds[3])} kts\nDTM: {directions[4]} @ {mag(speeds[4])} kts\n US: {directions[5]} @ {mag(speeds[5])} kts\n' + 
                        f' DS: {directions[6]} @ {mag(speeds[6])} kts\n',
                        color=gen_txt_clr, fontsize=10, verticalalignment='top', linespacing=2.2, alpha=0.6)  
        except IndexError:
            pass
        
        def sm_str(storm_motion, speeds, directions):
            
            if storm_motion in ['right_moving',
                                'left_moving',
                                'mean_wind']:
                return f"{str.upper(storm_motion.replace('_', ' '))} | {directions[0]} @ {mag(speeds[0])}"

            else:
                return f"USER DEFINED | {directions[0]} @ {mag(speeds[0])}"
            
        plt.figtext(0.228, 0.87, 
                f' SM: {sm_str(storm_motion, speeds, directions)} kts',
                color=gen_txt_clr, weight='bold', fontsize=10, verticalalignment='top', linespacing=2.2, alpha=0.6) 
    ################################################################
    
    
    
    
    #########################################################################
    ############################ OTHER AXES ################################# 
    #########################################################################
    
    
    
    #################################################################
    ### STREAMWISENESS AND RH W/HGT ###
    #################################################################
    # PLOT AXIS/LOC
    strmws_ax = plt.axes((0.8015, 0.11, 0.095, 0.32))
    plt.figtext(0.85, 0.40, f'SW ζ (%)', color=gen_txt_clr, weight='bold', fontsize=12, ha='center', alpha=0.7)
    strmws_ax.spines["top"].set_color(brdr_clr)
    strmws_ax.spines["left"].set_color(brdr_clr)
    strmws_ax.spines["right"].set_color(brdr_clr)
    strmws_ax.spines["bottom"].set_color(brdr_clr)
    strmws_ax.spines["bottom"].set_color(brdr_clr)   
    strmws_ax.set_facecolor(bckgrnd_clr) 

    #YTICKS
    strmws_ax.set_ylim(0, 3000)
    strmws_ax.set_yticklabels([])
    strmws_ax.set_ylabel(' ')
    strmws_ax.tick_params(axis='y', length = 0)
    strmws_ax.grid(True, axis='y')
    strmws_ax.tick_params(axis="x",direction="in", pad=-12)

    #XTICKS
    strmws_ax.set_xlim(40, 102)
    strmws_ax.set_xticks([50, 90])
    strmws_ax.set_xticklabels([50, 90], weight='bold', alpha=0.5, fontstyle='italic', color=gen_txt_clr)
    strmws_ax.set_xlabel(' ')

    #HGT LABLES 
    strmws_ax.text(47, 502 , '0.5km', fontsize=8, alpha=0.6, color=gen_txt_clr)
    strmws_ax.text(47, 1002, '1.0km', fontsize=8, alpha=0.6, color=gen_txt_clr)
    strmws_ax.text(47, 1502, '1.5km', fontsize=8, alpha=0.6, color=gen_txt_clr)
    strmws_ax.text(47, 2002, '2.0km', fontsize=8, alpha=0.6, color=gen_txt_clr)
    strmws_ax.text(47, 2502, '2.5km', fontsize=8, alpha=0.6, color=gen_txt_clr)

    if ma.is_masked(kinem['sm_u']) == False:
        plt.plot(kinem['swv_perc'][0:11],  intrp['zINTRP'][0:11],  color=hodo_color[0], lw=3, clip_on=True)
        plt.plot(kinem['swv_perc'][10:30], intrp['zINTRP'][10:30], color=hodo_color[1], lw=3, clip_on=True)
    
    else:
        warnings.warn("Streamwiseness could not be plotted (no valid storm motion/not enough data)", Warning)
    #################################################################

    
    
    
    #################################################################
    ### VORTICITY W/HGT ###
    #################################################################
    vort_ax = plt.axes((0.8965, 0.11, 0.095, 0.32))
    plt.figtext(0.945, 0.39, f'ζₜₒₜ & ζSW\n(/sec)', color=gen_txt_clr, weight='bold', fontsize=12, ha='center', alpha=0.7)
    vort_ax.spines["top"].set_color(brdr_clr)
    vort_ax.spines["left"].set_color(brdr_clr)
    vort_ax.spines["right"].set_color(brdr_clr)
    vort_ax.spines["bottom"].set_color(brdr_clr)
    vort_ax.spines["bottom"].set_color(brdr_clr)   
    vort_ax.set_facecolor(bckgrnd_clr) 

    #YTICKS
    vort_ax.tick_params(axis='y', length = 0)
    vort_ax.grid(True, axis='y')
    vort_ax.set_ylim(0, 3000)
    vort_ax.set_yticklabels([])
    vort_ax.set_ylabel(' ')
    vort_ax.tick_params(axis="x",direction="in", pad=-12)
        
    if ma.is_masked(kinem['sm_u']) == False: 
        #XTICKS
        vort_ax.set_xlabel(' ')
        vort_max = kinem['vort'][0:30].max()+0.005
        vort_min = kinem['vort'][0:30].min()-0.005
    
        vort_ax.set_xlim(vort_min-0.002, vort_max+0.002)
        vort_ax.set_xticks([(vort_min+0.005),(vort_max-0.005)])
        vort_ax.set_xticklabels([(np.round(vort_min+0.002,2)),(np.round(vort_max-0.002,2))], weight='bold', alpha=0.5, fontstyle='italic', color=gen_txt_clr)
 
        vort_ax.plot(kinem['swv'][0:30],  intrp['zINTRP'][0:30], color='orange', linewidth=3, alpha=0.8, label='SW ζ')
        vort_ax.plot(kinem['vort'][0:30], intrp['zINTRP'][0:30], color=gen_txt_clr,  linewidth=4, alpha=0.4, label='Total ζ')
        
    else:
        warnings.warn("Total Vorticity could not be plotted (no valid storm motion/not enough data)", Warning)
    #################################################################
    
    
    
    
    
    #################################################################
    ### SRW W/HGT ###
    #################################################################
    wind_ax = plt.axes((0.9915, 0.11, 0.095, 0.32))
    plt.figtext(1.037, 0.39, f'SR Wind\n(kts)', weight='bold', color=gen_txt_clr, fontsize=12, ha='center', alpha=0.7)
    wind_ax.spines["top"].set_color(brdr_clr)
    wind_ax.spines["left"].set_color(brdr_clr)
    wind_ax.spines["right"].set_color(brdr_clr)
    wind_ax.spines["bottom"].set_color(brdr_clr)
    wind_ax.spines["bottom"].set_color(brdr_clr)   
    wind_ax.set_facecolor(bckgrnd_clr) 
    plt.ylabel(' ')
    plt.xlabel(' ')

    #YTICKS
    wind_ax.set_ylim(0, 3000)
    wind_ax.grid(True, axis='y')
    wind_ax.set_yticklabels([])
    wind_ax.tick_params(axis='y', length = 0)
    wind_ax.tick_params(axis="x",direction="in", pad=-12)

    if ma.is_masked(kinem['sm_u']) == False:
        #XTICKS
        wind_max = kinem['srw'][0:30].max()+1
        wind_min = kinem['srw'][0:30].min()-1
        wind_ax.set_xlim(wind_min-5, wind_max+5)
        wind_ax.set_xticks([(wind_min)+2, (wind_max)-2])
        wind_ax.set_xticklabels([(int(wind_min)+2), (int(wind_max)-2)], weight='bold', alpha=0.5, fontstyle='italic', color=gen_txt_clr)

        #PLOT SR WIND  
        wind_ax.plot(kinem['srw'][0:11],  intrp['zINTRP'][0:11],  color=hodo_color[0], clip_on=True, linewidth=3, alpha=0.8, label='0-1 SR Wind')
        wind_ax.plot(kinem['srw'][10:30], intrp['zINTRP'][10:30], color=hodo_color[1], clip_on=True, linewidth=3, alpha=0.8, label='1-3 SR Wind')
        
    else:
        warnings.warn("Storm Relative Wind could not be plotted (no valid storm motion/not enough data)", Warning)
        
    #################################################################
        
        
        
    
    #################################################################
    ### WIND COMPONENTS PLOT ###
    #################################################################
    #PLOT AXES/LOC
    comp_ax = plt.axes((1.0865, 0.11, 0.095, 0.32))
    plt.figtext(1.135, 0.39, f'U-Wind (green)\nV-Wind (blue)', weight='bold', color=gen_txt_clr, fontsize=12, ha='center', alpha=0.7)
    comp_ax.spines["top"].set_color(brdr_clr)
    comp_ax.spines["left"].set_color(brdr_clr)
    comp_ax.spines["right"].set_color(brdr_clr)
    comp_ax.spines["bottom"].set_color(brdr_clr)
    comp_ax.spines["bottom"].set_color(brdr_clr)   
    comp_ax.set_facecolor(bckgrnd_clr) 

    max_comp = np.array([intrp['uINTRP'][0:30].max(), intrp['vINTRP'][0:30].max()]).max()
    min_comp = np.array([intrp['uINTRP'][0:30].min(), intrp['vINTRP'][0:30].min()]).min()

    #YTICKS
    comp_ax.set_ylim(intrp['zINTRP'][0], 3000)
    comp_ax.set_yticklabels([])
    plt.ylabel(' ')
    comp_ax.tick_params(axis='y', length = 0)
    comp_ax.grid(True, axis='y')

    #XTICKS
    comp_ax.set_xlim(min_comp - 10, max_comp + 5)
    comp_ax.set_xticks([(min_comp), (max_comp)])
    comp_ax.set_xticklabels([int(min_comp), int(max_comp)], weight='bold', alpha=0.5, fontstyle='italic', color=gen_txt_clr)

    comp_ax.tick_params(axis="x", direction="in", pad=-12)
    comp_ax.set_xlabel(' ')

    #PLOT THETA VS HGT
    plt.plot(intrp['uINTRP'], intrp['zINTRP'], color='forestgreen', linewidth=3.5, alpha=0.5, clip_on=True)
    plt.plot(intrp['vINTRP'], intrp['zINTRP'], color='cornflowerblue', linewidth=3.5, alpha=0.8, clip_on=True)

    #################################################################
    
    
    
    #########################################################################
    ############################## TEXT PLOTS ###############################
    #########################################################################

    
    
    #################################################################
    ### BOXES FOR TEXT ###
    #################################################################
    #                                 xloc   yloc   xsize  ysize
    fig.patches.extend([plt.Rectangle((0.8015, 0.43), 0.38, 0.45,
                                      edgecolor=brdr_clr, facecolor=bckgrnd_clr, linewidth=1, alpha=1,
                                      transform=fig.transFigure, figure=fig)])
    
    #################################################################
    ### KINEMATICS ###
    #################################################################
    # now some kinematic parameters
    met_per_sec = (units.m*units.m)/(units.sec*units.sec)
    plt.figtext( 0.88, 0.84, 'BS          SRH         SRW    SWζ%    SWζ', color=gen_txt_clr, weight='bold', fontsize=15, alpha=0.8)
    
    plt.figtext( 0.815, 0.81, f"0-.5ₖₘ:", weight='bold', fontsize=15, color=gen_txt_clr, alpha=0.9)
    plt.figtext( 0.875, 0.81, f"{mag(kinem['shear_0_to_500'])} kt", fontsize=15, color='deepskyblue', weight='bold')
    plt.figtext( 0.929, 0.81,  f"{mag(kinem['srh_0_to_500'])* met_per_sec:~P}", fontsize=15, color='deepskyblue', weight='bold')
    plt.figtext( 1.018, 0.81, f"{mag(kinem['srw_0_to_500'])} kt", fontsize=15, color='deepskyblue', weight='bold')
    plt.figtext( 1.08, 0.81, f"{mag(kinem['swv_perc_0_to_500'])}", fontsize=15, color='deepskyblue', weight='bold')
    plt.figtext( 1.125, 0.81,  f"{mag_round(kinem['swv_0_to_500'], 3)}", fontsize=15, color='deepskyblue', weight='bold')

    plt.figtext( 0.815, 0.78, f"0-1ₖₘ:", weight='bold', fontsize=15, color=gen_txt_clr, alpha=0.9)
    plt.figtext( 0.875, 0.78, f"{mag(kinem['shear_0_to_1000'])} kt", fontsize=15, color='mediumslateblue', weight='bold')
    plt.figtext( 0.929, 0.78,  f"{mag(kinem['srh_0_to_1000'])* met_per_sec:~P}", fontsize=15, color='mediumslateblue', weight='bold')
    plt.figtext( 1.018, 0.78, f"{mag(kinem['srw_0_to_1000'])} kt", fontsize=15, color='mediumslateblue', weight='bold')
    plt.figtext( 1.08, 0.78, f"{mag(kinem['swv_perc_0_to_1000'])}", fontsize=15, color='mediumslateblue', weight='bold')
    plt.figtext( 1.125, 0.78,  f"{mag_round(kinem['swv_0_to_1000'], 3)}", fontsize=15, color='mediumslateblue', weight='bold')

    plt.figtext( 0.815, 0.75, f"1-3ₖₘ:", weight='bold', fontsize=15, color=gen_txt_clr, alpha=0.9)
    plt.figtext( 0.875, 0.75, f"{mag(kinem['shear_1_to_3000'])} kt", fontsize=15, color='slateblue', weight='bold')
    plt.figtext( 0.929, 0.75,  f"{mag(kinem['srh_1_to_3000'])* met_per_sec:~P}", fontsize=15, color='slateblue', weight='bold')
    plt.figtext( 1.018, 0.75, f"{mag(kinem['srw_1_to_3000'])} kt", fontsize=15, color='slateblue', weight='bold')
    plt.figtext( 1.08, 0.75, f"{mag(kinem['swv_perc_1_to_3000'])}", fontsize=15, color='slateblue', weight='bold')
    plt.figtext( 1.125, 0.75,  f"{mag_round(kinem['swv_1_to_3000'], 3)}", fontsize=15, color='slateblue', weight='bold')

    plt.figtext( 0.815, 0.72, f"3-6ₖₘ:", weight='bold', fontsize=15, color=gen_txt_clr, alpha=0.9)
    plt.figtext( 0.875, 0.72, f"{mag(kinem['shear_3_to_6000'])} kt", fontsize=15, color='darkslateblue', weight='bold')
    plt.figtext( 0.929, 0.72,  f"{mag(kinem['srh_3_to_6000'])* met_per_sec:~P}", fontsize=15, color='darkslateblue', weight='bold')
    plt.figtext( 1.018, 0.72, f"{mag(kinem['srw_3_to_6000'])} kt", fontsize=15, color='darkslateblue', weight='bold')
    plt.figtext( 1.08, 0.72, f"{mag(kinem['swv_perc_3_to_6000'])}", fontsize=15, color='darkslateblue', weight='bold')
    plt.figtext( 1.125, 0.72,  f"{mag_round(kinem['swv_3_to_6000'], 3)}", fontsize=15, color='darkslateblue', weight='bold')
    
    
    plt.figtext( 0.815, 0.69,  f"6-9ₖₘ:", weight='bold', fontsize=15, color=gen_txt_clr, alpha=0.9)
    plt.figtext( 0.875, 0.69,  f"{mag(kinem['shear_6_to_9000'])} kt", fontsize=15, color='darkslateblue', weight='bold')
    plt.figtext( 0.929, 0.69,  f"{mag(kinem['srh_6_to_9000'])* met_per_sec:~P}", fontsize=15, color='darkslateblue', weight='bold')
    plt.figtext( 1.018, 0.69,  f"{mag(kinem['srw_6_to_9000'])} kt", fontsize=15, color='darkslateblue', weight='bold')
    plt.figtext( 1.08, 0.69,  f"{mag(kinem['swv_perc_6_to_9000'])}", fontsize=15, color='darkslateblue', weight='bold')
    plt.figtext( 1.125,  0.69,  f"{mag_round(kinem['swv_6_to_9000'], 3)}", fontsize=15, color='darkslateblue', weight='bold')
    
    #################################################################
    

    
        
    # PLOT TITLES ----------------------------------------------------------------------------------
    plt.figtext( 0.23, 0.89, left_title, ha='left', fontsize=16, color=gen_txt_clr)
    plt.figtext( 1.20, 0.89, f'{right_title}  ', ha='right', fontsize=16, color=gen_txt_clr)
    plt.figtext( 0.23, 0.92, top_title, ha='left', weight='bold', fontsize=22, color=gen_txt_clr) 
    plt.figtext( 0.23, 0.94, ' ') 
    plt.figtext( 0.225, 0.09, f'SOUNDERPY VERTICAL PROFILE ANALYSIS TOOL | (C) KYLE J GILLETT 2024',
                ha='left', color='cornflowerblue', alpha=0.8, weight='bold', fontsize=10)
    
    img = Image.open(urlopen('https://user-images.githubusercontent.com/100786530/251580013-2e9477c9-e36a-4163-accb-fe46780058dd.png'))
    #                  side-side  up-down  size   size
    imgax = fig.add_axes([0.22, 0.12, 0.08, 0.08], anchor='SE')
    imgax.imshow(img)
    imgax.axis('off')
    
    
    elapsed_time = time.time() - st
    print('> RUNTIME:', time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
    
    plt.tight_layout()
    
    return plt
    ########################################################################


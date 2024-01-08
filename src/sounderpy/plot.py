import metpy.calc as mpcalc
from metpy.units import units
from metpy.plots import SkewT, Hodograph

import matplotlib.lines    as mlines
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms

import numpy as np
import numpy.ma as ma

from urllib.request import urlopen
from PIL import Image
import warnings
import time

from .calc import *



#########################################################
#            SOUNDERPY SPYPLOT FUNCTIONS                #
# (C) KYLE J GILLETT, CENTRAL MICHIGAN UNIVERSTIY, 2023 #
#########################################################



#########################################################################
########################## FULL SOUNDING ################################

def __full_sounding(clean_data, color_blind, dark_mode):
    
    if dark_mode == True:
        gen_txt_clr = 'white'
        bckgrnd_clr = 'black'
        brdr_clr    = 'white'
        barb_clr    = 'white'
        shade_alpha = 0.06
        skw_ln_clr = 'white'
    else: 
        gen_txt_clr = 'black'
        bckgrnd_clr = 'white'
        brdr_clr    = 'black'
        barb_clr    = 'black'
        shade_alpha = 0.02
        skw_ln_clr = 'black'
    
    # record process time 
    st = time.time()  

    def find_nearest(array, value):
        array = np.asarray(array)
        nearest_idx = (np.abs(array - value)).argmin()
        return nearest_idx
    
    if color_blind == True:
        td_color = 'cornflowerblue'
    else:
        td_color = 'green'
        
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
    T  = clean_data['T']
    Td = clean_data['Td']
    p  = clean_data['p']
    z  = clean_data['z']
    u  = clean_data['u']
    v  = clean_data['v']
    wd = mpcalc.wind_direction(u, v)
    ws = mpcalc.wind_speed(u, v) 
    
    # calculate other sounding parameters using SounderPy Calc
    general, thermo, kinem, intrp = sounding_params(clean_data).calc()
    #################################################################
    
    
    
    #################################################################
    ### DETERMINE PLOT TITLE BASED ON THE DATA ###
    #################################################################
    
    if 'ACARS' in clean_data['site_info']['source']:
        top_title = f"ACARS AIRCRAFT OBSERVATION VERTICAL PROFILE"
        left_title = f"VALID: {clean_data['site_info']['valid-time'][1]}-{clean_data['site_info']['valid-time'][2]}-{clean_data['site_info']['valid-time'][0]} {clean_data['site_info']['valid-time'][3]}Z"
        right_title = f"{clean_data['site_info']['site-id']} - {clean_data['site_info']['site-name']} | {clean_data['site_info']['site-latlon'][0]}, {clean_data['site_info']['site-latlon'][1]}    " 

    elif 'BUFKIT' in clean_data['site_info']['source']:
        top_title = f"BUFKIT MODEL FORECAST PROFILE | {clean_data['site_info']['run-time'][3]}Z {clean_data['site_info']['model']} {clean_data['site_info']['fcst-hour']}"
        left_title = f" VALID: {clean_data['site_info']['valid-time'][1]}/{clean_data['site_info']['valid-time'][2]}/{clean_data['site_info']['valid-time'][0]} {clean_data['site_info']['valid-time'][3]}Z"
        right_title = f"{clean_data['site_info']['site-id']} - {clean_data['site_info']['site-name']} | {clean_data['site_info']['site-latlon'][0]}, {clean_data['site_info']['site-latlon'][1]}    " 

    elif 'RAOB' in clean_data['site_info']['source']:
        top_title = "RAOB OBSERVED VERTICAL PROFILE"
        left_title = f"VALID: {clean_data['site_info']['valid-time'][1]}-{clean_data['site_info']['valid-time'][2]}-{clean_data['site_info']['valid-time'][0]} {clean_data['site_info']['valid-time'][3]}Z"
        right_title = f"{clean_data['site_info']['site-id']} - {clean_data['site_info']['site-name']} | {clean_data['site_info']['site-latlon'][0]}, {clean_data['site_info']['site-latlon'][1]}    " 

    elif 'REANALYSIS' in clean_data['site_info']['source']:
        top_title = f"MODEL REANALYSIS VERTICAL PROFILE | {clean_data['site_info']['valid-time'][3]}Z {clean_data['site_info']['model']} {clean_data['site_info']['fcst-hour']}"
        left_title = f"{clean_data['site_info']['run-time'][3]}Z {clean_data['site_info']['model']} {clean_data['site_info']['fcst-hour']} | VALID: {clean_data['site_info']['valid-time'][1]}/{clean_data['site_info']['valid-time'][2]}/{clean_data['site_info']['valid-time'][0]} {clean_data['site_info']['valid-time'][3]}Z"
        right_title = f"{clean_data['site_info']['site-latlon'][0]}, {clean_data['site_info']['site-latlon'][1]}    " 
        
    else:
        top_title = clean_data['site_info']['source']
        left_title = f"VALID: {clean_data['site_info']['valid-time'][1]}-{clean_data['site_info']['valid-time'][2]}-{clean_data['site_info']['valid-time'][0]} {clean_data['site_info']['valid-time'][3]}Z"
        right_title = f"{clean_data['site_info']['site-id']} - {clean_data['site_info']['site-name']} | {clean_data['site_info']['site-latlon'][0]}, {clean_data['site_info']['site-latlon'][1]}    " 
        
    ################################################################    
        
    
    #########################################################################
    ################################ SKEW-T ################################# 
    #########################################################################
    
    
    
    #################################################################
    ### CREATE FIGURE ###
    #################################################################
    #################################################################
    # Define figure and skew-t
    fig = plt.figure(figsize=(22,13), linewidth=10, edgecolor=brdr_clr)          # create figure                                   # create figure
    skew = SkewT(fig, rotation=47, rect=(0.1124, 0.1005, 0.60, 0.85))  
    skew.ax.set_box_aspect(0.87)
    skew.ax.zorder = 5
    # Define axis bounds 
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
    skew.plot_dry_adiabats(color=skw_ln_clr, linewidth=0.5, alpha=0.7)   
    skew.plot_moist_adiabats(color=skw_ln_clr, linewidth=0.5, alpha=0.7)
    skew.plot_mixing_lines(color=skw_ln_clr, linewidth=0.2, alpha=0.7)
    # add basic temperature lines
    twline = skew.plot(p, general['wet_bulb'], '#3d8aff', linewidth=1, label='WETBULB TEMP', alpha=0.3)
    tvline = skew.plot(p, general['virt_temp'], '#0b6431', linestyle='--', linewidth=1, label='VIRTUAL TEMP', alpha=0.3)  
    tdline = skew.plot(p, Td, td_color, linewidth=5, label='DEWPOINT')
    tline1 = skew.plot(p, T, 'red', linewidth=5, label='TEMPERATURE')  
    
    # add parcel lines and CAPE shading
    if thermo['sbcape'] > 0:
        sbparcelline = skew.plot(thermo['sbP_trace'], thermo['sbT_trace'], linestyle='--',
                                 linewidth=1, color='red', alpha=1, label='SB PARCEL')    
    if thermo['mlcape'] > 0:
        mlparcelline = skew.plot(thermo['mlP_trace'], thermo['mlT_trace'], color='orangered', linestyle='--',
                                 linewidth=1, alpha=1, label='ML PARCEL') 
    if thermo['mucape'] > 0:
        muparcelline = skew.plot(thermo['muP_trace'], thermo['muT_trace'], color='orange', linestyle='--',  
                                 linewidth=1, alpha=1, label='MU PARCEL')

    skew.plot(thermo['dparcel_p'], thermo['dparcel_T'], linestyle='--',linewidth=0.7, color='purple', 
              alpha=0.8, label='DWNDRFT PARCEL')
    #################################################################
    
    
    
    #################################################################
    ### PLOT SKEW T ANNOTATIONS ###
    #################################################################    
    hgt_lvls =[]
    for key in intrp['hgt_lvls'].keys():
        hgt_lvls.append(intrp['hgt_lvls'][key])
    hgt_lvls.pop(0) 

    for key in hgt_lvls[1::4]:
        trans, _, _ = skew.ax.get_yaxis_text1_transform(0)
        skew.ax.text(0.048, intrp['pINTRP'][key], f"{int(intrp['zINTRP'][key]/1000)}km", 
                     fontsize=13, transform=trans, alpha=0.6, weight='bold', color=gen_txt_clr)   

    sfc = mpcalc.height_to_pressure_std(general['elevation']*units.m)
    skew.ax.text(0.048, p[0], '-SFC-', fontsize=13, transform=trans, alpha=0.6, weight='bold', color=gen_txt_clr) # plot 'SFC' @ surface pressure
    
    
    # SFC TEMPERATURE AND DEWPOINT ANNOTATIONS---------------------------------------------
    T_degF = np.round(T.to(units.degF), 1)
    T_degF_label = '{}°F'.format(int(T_degF[0].magnitude))                             
    plt.annotate(T_degF_label, (T[0], p[0]), textcoords="offset points", xytext=(16,-15),
                     fontsize=12, color='red', weight='bold', alpha=0.7, ha='center')   
    Td_degF = np.round(Td.to(units.degF), 1) 
    Td_degF_label = '{}°F'.format(int(Td_degF[0].magnitude))                             
    plt.annotate(Td_degF_label,(Td[0], p[0]),textcoords="offset points",xytext=(-16,-15), 
                     fontsize=12, color=td_color, weight='bold', alpha=0.7, ha='center') 

    # PARCEL HEIGHT ANNOTATIONS------------------------------------------------------------- 
    plt.text((0.82), (thermo['sb_lcl_p']), "←SBLCL", weight='bold',color='gray',
             alpha=0.9, fontsize=15, transform=skew.ax.get_yaxis_transform(), clip_on=True)
    if ma.is_masked(thermo['mu_lfc_z']) == True:
        plt.text((0.82), (thermo['sb_lfc_p']), "←SBLFC", weight='bold',color='gray',
                 alpha=0.9, fontsize=15, transform=skew.ax.get_yaxis_transform(), clip_on=True)
        plt.text((0.82), (thermo['sb_el_p']), "←SBEL", weight='bold',color='gray', 
                 alpha=0.9, fontsize=15, transform=skew.ax.get_yaxis_transform(), clip_on=True)
    else: 
        plt.text((0.82), (thermo['mu_lfc_p']), "←MULFC", weight='bold',color='gray',
                 alpha=0.9, fontsize=15, transform=skew.ax.get_yaxis_transform(), clip_on=True)
        plt.text((0.82), (thermo['mu_el_p']), "←MUEL", weight='bold',color='gray',
                 alpha=0.9, fontsize=15, transform=skew.ax.get_yaxis_transform(), clip_on=True)

    #FRREZING POINT ANNOTATION--------------------------------------------------------------
    if ma.is_masked(general['frz_pt_z']) == False:
        if general['frz_pt_z'] >= 50*units.m:
            plt.text((0.82), (general['frz_pt_p']), "←FRZ", weight='bold',color='cornflowerblue',  
                     alpha=0.6, fontsize=13.5, transform=skew.ax.get_yaxis_transform(), clip_on=True)

    #PBL TOP POINT ANNOTATION---------------------------------------------------------------                   
    plt.text((0.82), (thermo['pbl_top']), "←PBL", weight='bold',color='gray', 
             alpha=0.9, fontsize=10, transform=skew.ax.get_yaxis_transform(), clip_on=True)
    

    
    # 0-3km & 0-6km CAPE ANNOTATIONS--------------------------------------------------------
    if thermo['mu3cape'] > 10:
        idx = find_nearest(thermo['muZ_trace'], 3000)
        cape03_label = " ←{}J/kg".format(mag(thermo['mu3cape']))     
        plt.annotate(cape03_label,  ((thermo['muT_trace'][idx] + 7), intrp['pINTRP'][intrp['hgt_lvls']['h3']]),  textcoords="offset points",  xytext=(15, 0), 
                     color='red', alpha=0.4, fontsize=13.5, weight='bold', ha='right')  
            
    if thermo['mu6cape'] > thermo['mu3cape']:
        idx = find_nearest(thermo['muZ_trace'], 6000)
        cape06_label = " ←{}J/kg".format(mag(thermo['mu6cape']))                                            
        plt.annotate(cape06_label, ((thermo['muT_trace'][idx] + 7), intrp['pINTRP'][intrp['hgt_lvls']['h6']]), textcoords="offset points",  xytext=(15, 0), 
                         color='red', alpha=0.4, fontsize=13.5, weight='bold', ha='right') 
            
    if thermo['mucape'] > thermo['mu6cape']:
            cape_label = "←{}J/kg".format(mag(thermo['mucape']))                                            
            plt.annotate(cape_label,((thermo['mu_el_T'] + 7), thermo['mu_el_p']), textcoords="offset points",  xytext=(10, 0), 
                         color='red', alpha=0.4, fontsize=13.5, weight='bold', ha='center') 
            
    # MAX LAPSE RATE ANNOTATION---------------------------------
    x_start, x_end = 0.80, 0.81
    x_mid = (x_start + x_end)/2
    Lapse_line = plt.Line2D([x_mid, x_mid], (thermo['lr_max'][1], thermo['lr_max'][2]), color='purple', alpha=0.3, transform=skew.ax.get_yaxis_transform())
    plt.text((x_start-0.01), (thermo['lr_max'][2]-14), f"{np.round(thermo['lr_max'][0],1)}", 
             color='purple', weight='bold', fontsize=13, alpha=0.3, transform=skew.ax.get_yaxis_transform())
    skew.ax.add_artist(Lapse_line)

    # EFFECTIVE INFLOW LAYER ANNOTATION--------------------------
    x_start, x_end = 0.2, 0.22
    x_mid = (x_start + x_end)/2
    plt.text((x_start+0.01), (kinem['eil'][1]-10), "EIL", weight='bold',color='lightblue', alpha=0.95, ha='center', fontsize=13, transform=skew.ax.get_yaxis_transform())
    EIL_line = plt.Line2D([x_mid, x_mid], (kinem['eil'][0], kinem['eil'][1]), color='lightblue', alpha=0.95, transform=skew.ax.get_yaxis_transform())
    skew.ax.add_artist(EIL_line)

    if T[0].m < 5:
        dgz_p, dgz_T, dgz_Td = mpcalc.get_layer(p, T, Td, 
                                                bottom=thermo['dgz'][0]*units.hPa, depth=(thermo['dgz'][0]-thermo['dgz'][1])*units.hPa)
        skew.shade_area(y=dgz_p, x1=dgz_Td, x2=dgz_T, color='blue', alpha=0.03, zorder=1) 
        skew.plot(dgz_p, dgz_T, 'blue', linewidth=4)                                     # plot temperature trace 
        skew.plot(dgz_p, dgz_Td, 'blue', linewidth=4) 

    else:
        if thermo['hgz'][1] > clean_data['p'][0].m:
            #HAIL GROWTH ZONE ANNOTATION------------------------------------------------------------------------------------
            hgz_p, hgz_T, hgz_Td = mpcalc.get_layer(p, T, Td, 
                                                    bottom=thermo['hgz'][0]*units.hPa, depth=(thermo['hgz'][0]-thermo['hgz'][1])*units.hPa)
            skew.shade_area(y=hgz_p, x1=hgz_Td, x2=hgz_T, color='green', alpha=0.03, zorder=1) 
            skew.plot(hgz_p, hgz_T, 'limegreen', linewidth=4)                                     # plot temperature trace 
            skew.plot(hgz_p, hgz_Td, 'limegreen', linewidth=4) 
    #################################################################
    
    
    
    
    #################################################################
    ### PLOT SKEW T WIND BARBS ###
    #################################################################
    interval = np.logspace(2.113, 3, 30) *units.hPa # Arrange wind barbs for best fit
    idx = mpcalc.resample_nn_1d(p, interval) # Resample wind barbs for best fit

    # create blank barbs for small dot at the start of each actual barb
    blank_len = len(u[idx])     
    blank = np.zeros(blank_len)  
    skew.plot_barbs(pressure=p[idx], u=blank, v=blank, xloc=0.955, fill_empty=True, color=barb_clr,
                    sizes=dict(emptybarb=0.075, width=0.18, height=0.4))

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
    if z.max().m > 9001: 
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
    hod_ax = plt.axes((0.53, 0.1003, 0.85, 0.85))
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

    h.plot(intrp['uINTRP'][intrp['hgt_lvls']['h05']],intrp['vINTRP'][intrp['hgt_lvls']['h05']],marker='.', 
           color='white', alpha=1, markersize=30, clip_on=True, zorder=5)
    h.ax.annotate(str('.5'),(intrp['uINTRP'][intrp['hgt_lvls']['h05']],intrp['vINTRP'][intrp['hgt_lvls']['h05']]),
                  weight='bold', fontsize=12, color='black',xytext=(0,-3.2),textcoords='offset pixels',horizontalalignment='center',clip_on=True, zorder=6) 

    hgt_lvls = [] 
    for key in intrp['hgt_lvls'].keys():
        hgt_lvls.append(intrp['hgt_lvls'][key])
    hgt_lvls.pop(0) 

    for i in hgt_lvls[1::2]:
        h.plot(intrp['uINTRP'][i],intrp['vINTRP'][i], marker='.', color='white', alpha=1, markersize=30, zorder=5)
        h.ax.annotate(str(int(round(intrp['zINTRP'][i]/1000,0))),(intrp['uINTRP'][i],intrp['vINTRP'][i]), 
                      weight='bold', fontsize=12, color='black',xytext=(0,-3.2),textcoords='offset pixels',horizontalalignment='center',clip_on=True, zorder=6) 
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
        h.ax.text((kinem['sm_rm'][0]+0.5), (kinem['sm_rm'][1]-0.5), 'RM', weight='bold', ha='left', fontsize=14, alpha=0.9, color=gen_txt_clr)
        h.ax.text((kinem['sm_lm'][0]+0.5), (kinem['sm_lm'][1]-0.5), 'LM', weight='bold', ha='left', fontsize=14, alpha=0.9, color=gen_txt_clr)
        h.ax.text((kinem['sm_mw'][0]+0.5), (kinem['sm_mw'][1]-0.5), 'MW', weight='bold', ha='left', fontsize=14, alpha=0.9, color=gen_txt_clr)
        h.ax.arrow(0,0,kinem['sm_u']-0.3, kinem['sm_v']-0.3, linewidth=3, color=gen_txt_clr, alpha=0.2, 
                   label='Bunkers RM Vector', length_includes_head=True, head_width=0.5)
        # DEVIANT TORNADO MOTION
        h.ax.text(kinem['dtm'][0], (kinem['dtm'][1] + 2), 'DTM', weight='bold', fontsize=10, color='brown', ha='center')
        h.plot(kinem['dtm'][0], kinem['dtm'][1], marker='v', color='brown', markersize=8, alpha=0.8, ls='', label='DEVIANT TORNADO MOTION')
    
    # EFFECTIVE INFLOW LAYER
    ebot = h.ax.plot((kinem['sm_u'], intrp['uINTRP'][find_nearest(p,kinem['eil'][0])]), (kinem['sm_v'], intrp['vINTRP'][find_nearest(p,kinem['eil'][0])]),  
                     linestyle='--', linewidth=2.3, alpha=0.5, color='lightblue', label='Effective Inflow Layer')
    etop = h.ax.plot((kinem['sm_u'], intrp['uINTRP'][find_nearest(p,kinem['eil'][1])]), (kinem['sm_v'], intrp['vINTRP'][find_nearest(p,kinem['eil'][1])]),  
                     linestyle='--', linewidth=2.3, alpha=0.5, color='lightblue')
    # EFFECTIVE INFLOW LAYER SRH FILL
    fill_srh = h.ax.fill(np.append(intrp['uINTRP'][find_nearest(p, kinem['eil'][0]):find_nearest(p, kinem['eil'][1])+1], kinem['sm_u']), 
                         np.append(intrp['vINTRP'][find_nearest(p, kinem['eil'][0]):find_nearest(p, kinem['eil'][1])+1], kinem['sm_v']),
                         'lightblue',alpha=0.1, label='EIL SRH')
    
    h.ax.text(kinem['mcs'][0], kinem['mcs'][1], 'UP', weight='bold', fontsize=12, color='orange', ha='center', alpha=0.5, clip_on=True)
    h.ax.text(kinem['mcs'][2], kinem['mcs'][3], 'DN', weight='bold', fontsize=12, color='orange', ha='center', alpha=0.5, clip_on=True)
    
    try:
        keys = ['sm_rm', 'sm_lm', 'sm_mw', 'dtm'] 

        speeds = []
        directions = []

        for key in keys:
            speeds.append(mpcalc.wind_speed(kinem[key][0]*units.kts,kinem[key][1]*units.kts).m) 
            directions.append(mpcalc.angle_to_direction((np.degrees(np.arctan2(kinem[key][0],kinem[key][1]))), full=False, level=3))

        speeds.append(mpcalc.wind_speed(kinem['mcs'][0]*units.kts,kinem['mcs'][1]*units.kts).m) 
        directions.append(mpcalc.angle_to_direction((np.degrees(np.arctan2(kinem['mcs'][0],kinem['mcs'][1]))), full=False, level=3))   

        speeds.append(mpcalc.wind_speed(kinem['mcs'][2]*units.kts,kinem['mcs'][3]*units.kts).m) 
        directions.append(mpcalc.angle_to_direction((np.degrees(np.arctan2(kinem['mcs'][2],kinem['mcs'][3]))), full=False, level=3)) 
    
        #plot Bunkers Storm Motion & DTM Data in box on Hodograph
        plt.figtext(0.748, 0.9451, 
                    f' RM: {directions[0]} @ {mag(speeds[0])} kts\n LM: {directions[1]} @ {mag(speeds[1])} kts\n MW: {directions[2]} @ {mag(speeds[2])} kts\n'+
                    f' DTM: {directions[3]} @ {mag(speeds[3])} kts\n US: {directions[4]} @ {mag(speeds[4])} kts\n DS: {directions[5]} @ {mag(speeds[5])} kts', 
                    color=gen_txt_clr, fontsize=15, verticalalignment='top', linespacing=2.2, alpha=0.6)  
    except IndexError:
        pass
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
                        temp_adv_ax.barh(top_arr[i], thermo['temp_adv'][0][i], align='center', height=bot_arr[i]-top_arr[i], edgecolor='black', alpha=0.3, color=temp_adv_bxclr)
                        if thermo['temp_adv'][0][i] > 0:
                            temp_adv_ax.annotate((np.round(thermo['temp_adv'][0][i],1)), xy=(0.3, top_arr[i]+10), color=gen_txt_clr, textcoords='data', ha='left', weight='bold')
                        if thermo['temp_adv'][0][i] < 0:
                            temp_adv_ax.annotate((np.round(thermo['temp_adv'][0][i],1)), xy=(-0.3, top_arr[i]+10), color=gen_txt_clr, textcoords='data', ha='right', weight='bold')
    temp_adv_ax.axvline(x=0, color=gen_txt_clr, linewidth=1, linestyle='--', clip_on=True)
    #################################################################
    
    
    
    #################################################################
    ### STREAMWISENESS AND RH W/HGT ###
    #################################################################
    # PLOT AXIS/LOC
    strmws_ax = plt.axes((0.945, -0.13, 0.065, 0.23))
    plt.figtext(0.978, 0.07, f'SW ζ (%)', color=gen_txt_clr, weight='bold', fontsize=12, ha='center', alpha=0.7)
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

    if ma.is_masked(kinem['sm_rm']) == False:
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
    vort_ax = plt.axes((1.01, -0.13, 0.066, 0.23))
    plt.figtext(1.045, 0.06, f'ζₜₒₜ & ζSW\n(/sec)', color=gen_txt_clr, weight='bold', fontsize=12, ha='center', alpha=0.7)
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
        
    if ma.is_masked(kinem['sm_rm']) == False: 
        #XTICKS
        vort_ax.set_xlabel(' ')
        vort_max = kinem['vort'][0:30].max()+0.005
        vort_min = kinem['vort'][0:30].min()-0.005
    
        vort_ax.set_xlim(vort_min-0.002, vort_max+0.002)
        vort_ax.set_xticks([(vort_min+0.005),(vort_max-0.005)])
        vort_ax.set_xticklabels([(np.round(vort_min+0.002,2)),(np.round(vort_max-0.002,2))], weight='bold', alpha=0.5, fontstyle='italic', color=gen_txt_clr)
 
        vort_ax.plot(kinem['swv'][0:30],  intrp['zINTRP'][0:30], color='orange', linewidth=3, alpha=0.8, label='SW ζ')
        vort_ax.plot(kinem['vort'][0:30], intrp['zINTRP'][0:30], color='black',  linewidth=4, alpha=0.4, label='Total ζ')
        
        if ma.is_masked(kinem['eil_z'][0]) == False:
            vort_ax.fill_between(x=(vort_min-0.002, vort_max+0.002), y1=kinem['eil_z'][0], y2=kinem['eil_z'][1], color='lightblue', alpha=0.2)
    else:
        warnings.warn("Total Vorticity could not be plotted (no valid storm motion/not enough data)", Warning)
    #################################################################
    
    
    
    
    
    #################################################################
    ### SRW W/HGT ###
    #################################################################
    wind_ax = plt.axes((1.0755, -0.13, 0.066, 0.23))
    plt.figtext(1.108, 0.06, f'SR Wind\n(kts)', weight='bold', color=gen_txt_clr, fontsize=12, ha='center', alpha=0.7)
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

    if ma.is_masked(kinem['sm_rm']) == False:
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
    theta_ax = plt.axes((1.141, -0.13, 0.0653, 0.23))
    plt.figtext(1.175, 0.06, f'Theta-e &\nTheta (K)', weight='bold', color=gen_txt_clr, fontsize=12, ha='center', alpha=0.7)
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
    #################################################################
    
    
    
    #################################################################
    ### THERMODYNAMICS ###
    #################################################################
    plt.figtext( 0.17, 0.07, '   CAPE         6CAPE         3CAPE        NCAPE          CIN            LCL', color=gen_txt_clr, weight='bold', fontsize=15)
    #SBCAPE
    plt.figtext( 0.13, 0.04,  f"SB:", weight='bold',   fontsize=15, color=gen_txt_clr)
    plt.figtext( 0.17, 0.04,  f"{mag(thermo['sbcape'])} J/kg",  fontsize=15, color='orangered', weight='bold')
    plt.figtext( 0.235, 0.04, f"{mag(thermo['sb6cape'])} J/kg", fontsize=15, color='orangered', weight='bold')
    plt.figtext( 0.302, 0.04, f"{mag(thermo['sb3cape'])} J/kg", fontsize=15, color='orangered', weight='bold')
    plt.figtext( 0.37, 0.04,  f"{mag(thermo['sb_ncape'])}",     fontsize=15, color='orangered', weight='bold')
    plt.figtext( 0.425, 0.04, f"{mag(thermo['sbcin'])} J/kg",   fontsize=15, color='orangered', weight='bold')
    plt.figtext( 0.485, 0.04, f"{mag(thermo['sb_lcl_z'])} m",   fontsize=15, color='orangered', weight='bold')
    #MUCAPE
    plt.figtext( 0.13, 0.01,  f"MU:", weight='bold',   fontsize=15, color=gen_txt_clr)
    plt.figtext( 0.17, 0.01,  f"{mag(thermo['mucape'])} J/kg",  fontsize=15, color='red', weight='bold')
    plt.figtext( 0.235, 0.01, f"{mag(thermo['mu6cape'])} J/kg", fontsize=15, color='red', weight='bold')
    plt.figtext( 0.302, 0.01, f"{mag(thermo['mu3cape'])} J/kg", fontsize=15, color='red', weight='bold')
    plt.figtext( 0.37, 0.01,  f"{mag(thermo['mu_ncape'])}",     fontsize=15, color='red', weight='bold')
    plt.figtext( 0.425, 0.01, f"{mag(thermo['mucin'])} J/kg",   fontsize=15, color='red', weight='bold')
    plt.figtext( 0.485, 0.01, f"{mag(thermo['mu_lcl_z'])} m",   fontsize=15, color='red', weight='bold')
    #MLCAPE
    plt.figtext( 0.13, -0.02,  f'ML:', weight='bold', fontsize=15, color=gen_txt_clr)
    plt.figtext( 0.17, -0.02,  f"{mag(thermo['mlcape'])} J/kg",  fontsize=15, color='darkred', weight='bold')
    plt.figtext( 0.235, -0.02, f"{mag(thermo['ml6cape'])} J/kg", fontsize=15, color='darkred', weight='bold')
    plt.figtext( 0.302, -0.02, f"{mag(thermo['ml3cape'])} J/kg", fontsize=15, color='darkred', weight='bold')
    plt.figtext( 0.37, -0.02,  f"{mag(thermo['ml_ncape'])}",     fontsize=15, color='darkred', weight='bold')
    plt.figtext( 0.425, -0.02, f"{mag(thermo['mlcin'])} J/kg",   fontsize=15, color='darkred', weight='bold')
    plt.figtext( 0.485, -0.02, f"{mag(thermo['ml_lcl_z'])} m",   fontsize=15, color='darkred', weight='bold')
    # ECAPE
    plt.figtext( 0.13, -0.061, f"ECAPE:", weight='bold', fontsize=15, color=gen_txt_clr)
    plt.figtext( 0.17, -0.061, f" {mag(thermo['ecape'])} J/kg", fontsize=15, color='brown', alpha=0.7, weight='bold')
    #DCAPE
    plt.figtext( 0.228, -0.061, f"DCAPE:", weight='bold', fontsize=15, color=gen_txt_clr)
    plt.figtext( 0.27,  -0.061, f"{mag(thermo['dcape'])} J/kg", fontsize=15, color='cornflowerblue', alpha=0.7, weight='bold')
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

    plt.figtext( 0.54,  -0.115, f"PWAT:", fontsize=15, weight='bold', color=gen_txt_clr)
    plt.figtext( 0.588, -0.115, f"{mag(general['pwat'])} in", fontsize=15, color='darkgreen', ha='center', weight='bold')
    plt.figtext( 0.625, -0.115, f"Tw:", fontsize=15, weight='bold', color=gen_txt_clr, ha='center')
    plt.figtext( 0.65,  -0.115, f"{mag(general['wet_bulb'][0])} °C", fontsize=15, color='darkgreen', ha='center', weight='bold')
    #################################################################
    
    
    
    #################################################################
    ### KINEMATICS ###
    #################################################################   
    met_per_sec = (units.m*units.m)/(units.sec*units.sec)

    plt.figtext( 0.735, 0.07, 'BS         SRH       SRW    SWζ%    SWζ', color=gen_txt_clr, weight='bold', fontsize=15, alpha=0.8)

    plt.figtext( 0.689, 0.04, f"0-.5ₖₘ:", weight='bold', fontsize=15, color=gen_txt_clr, alpha=0.9)
    plt.figtext( 0.732, 0.04, f"{mag(kinem['shear_0_to_500'])} kt", fontsize=15, color='deepskyblue', weight='bold')
    plt.figtext( 0.769, 0.04, f"{mag(kinem['srh_0_to_500'])* met_per_sec:~P}", fontsize=15, color='deepskyblue', weight='bold')
    plt.figtext( 0.828, 0.04, f"{mag(kinem['srw_0_to_500'])} kt", fontsize=15, color='deepskyblue', weight='bold')
    plt.figtext( 0.870, 0.04, f"{mag(kinem['swv_perc_0_to_500'])}", fontsize=15, color='deepskyblue', weight='bold')
    plt.figtext( 0.905, 0.04, f"{mag_round(kinem['swv_0_to_500'], 3)}", fontsize=15, color='deepskyblue', weight='bold')

    plt.figtext( 0.689, 0.01, f"0-1ₖₘ:", weight='bold', fontsize=15, color=gen_txt_clr, alpha=0.9)
    plt.figtext( 0.732, 0.01, f"{mag(kinem['shear_0_to_1000'])} kt", fontsize=15, color='mediumslateblue', weight='bold')
    plt.figtext( 0.769, 0.01, f"{mag(kinem['srh_0_to_1000'])* met_per_sec:~P}", fontsize=15, color='mediumslateblue', weight='bold')
    plt.figtext( 0.828, 0.01, f"{mag(kinem['srw_0_to_1000'])} kt", fontsize=15, color='mediumslateblue', weight='bold')
    plt.figtext( 0.870, 0.01, f"{mag(kinem['swv_perc_1_to_3000'])}", fontsize=15, color='mediumslateblue', weight='bold')
    plt.figtext( 0.905, 0.01, f"{mag_round(kinem['swv_3_to_6000'], 3)}", fontsize=15, color='mediumslateblue', weight='bold')

    plt.figtext( 0.689, -0.02, f"1-3ₖₘ:", weight='bold', fontsize=15, color=gen_txt_clr, alpha=0.9)
    plt.figtext( 0.732, -0.02, f"{mag(kinem['shear_1_to_3000'])} kt", fontsize=15, color='slateblue', weight='bold')
    plt.figtext( 0.769, -0.02, f"{mag(kinem['srh_1_to_3000'])* met_per_sec:~P}", fontsize=15, color='slateblue', weight='bold')
    plt.figtext( 0.828, -0.02, f"{mag(kinem['srw_1_to_3000'])} kt", fontsize=15, color='slateblue', weight='bold')
    plt.figtext( 0.870, -0.02, f"{mag(kinem['swv_perc_1_to_3000'])}", fontsize=15, color='slateblue', weight='bold')
    plt.figtext( 0.905, -0.02, f"{mag_round(kinem['swv_1_to_3000'], 3)}", fontsize=15, color='slateblue', weight='bold')

    plt.figtext( 0.689, -0.05, f"3-6ₖₘ:", weight='bold', fontsize=15, color=gen_txt_clr, alpha=0.9)
    plt.figtext( 0.732, -0.05, f"{mag(kinem['shear_3_to_6000'])} kt", fontsize=15, color='darkslateblue', weight='bold')
    plt.figtext( 0.769, -0.05, f"{mag(kinem['srh_3_to_6000'])* met_per_sec:~P}", fontsize=15, color='darkslateblue', weight='bold')
    plt.figtext( 0.828, -0.05, f"{mag(kinem['srw_3_to_6000'])} kt", fontsize=15, color='darkslateblue', weight='bold')
    plt.figtext( 0.870, -0.05, f"{mag(kinem['swv_perc_3_to_6000'])}", fontsize=15, color='darkslateblue', weight='bold')
    plt.figtext( 0.905, -0.05, f"{mag_round(kinem['swv_3_to_6000'], 3)}", fontsize=15, color='darkslateblue', weight='bold')
    #################################################################
    
    
    
    #########################################################################
    ############################# PLOT EXTRAS ############################### 
    #########################################################################   
    plt.figtext( 0.125, 0.985, top_title, weight='bold', ha='left', fontsize=30, color=gen_txt_clr)
    plt.figtext( 0.125, 0.959, left_title, ha='left', fontsize=23, color=gen_txt_clr)
    plt.figtext( 1.22, 0.959, right_title, ha='right', fontsize=23, color=gen_txt_clr)
    skewleg1 = skew.ax.legend(loc='upper left', framealpha=0.3, labelcolor=gen_txt_clr)
    
    # plot author credit information 
    plt.figtext( 0.34, -0.105, 'SOUNDERPY VERTICAL PROFILE ANALYSIS TOOL', 
                fontsize=20, ha='center', color='cornflowerblue', weight='bold', alpha=0.6)
    plt.figtext( 0.34, -0.125, '(C) KYLE J GILLETT 2023, CENTRAL MICHIGAN UNIVERSITY | AVAILABLE ON GITHUB', 
                fontsize=12, ha='center', color='cornflowerblue', weight='bold', alpha=0.6)

    # sounderpy logo
    img = Image.open(urlopen('https://user-images.githubusercontent.com/100786530/251580013-2e9477c9-e36a-4163-accb-fe46780058dd.png'))
    #                  side-side  up-down  size   size
    imgax = fig.add_axes([0.115, -0.135, 0.05, 0.05], anchor='SE', zorder=3)
    imgax.imshow(img)
    imgax.axis('off')
    
    plt.tight_layout()
    
    print('> COMPLETE --------')
    elapsed_time = time.time() - st
    print('> RUNTIME:', time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
    
    return plt










#########################################################################
########################## SIMPLE SOUNDING ##############################

def __simple_sounding(clean_data, color_blind, dark_mode):
    
    if dark_mode == True:
        gen_txt_clr = 'white'
        bckgrnd_clr = 'black'
        brdr_clr    = 'white'
        barb_clr    = 'white'
        shade_alpha = 0.06
        skw_ln_clr  = 'white'
    else: 
        gen_txt_clr = 'black'
        bckgrnd_clr = 'white'
        brdr_clr    = 'black'
        barb_clr    = 'black'
        shade_alpha = 0.02
        skw_ln_clr  = 'black'
    
    # record process time 
    st = time.time()   
    
    def find_nearest(array, value):
        array = np.asarray(array)
        nearest_idx = (np.abs(array - value)).argmin()
        return nearest_idx
    
    if color_blind == True:
        td_color = 'cornflowerblue'
    else:
        td_color = 'green'
        
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
    T  = clean_data['T']
    Td = clean_data['Td']
    p  = clean_data['p']
    z  = clean_data['z']
    u  = clean_data['u']
    v  = clean_data['v']
    wd = mpcalc.wind_direction(u, v)
    ws = mpcalc.wind_speed(u, v) 
    
    # calculate other sounding parameters using SounderPy Calc
    general, thermo, kinem, intrp = sounding_params(clean_data).calc()
    #################################################################
    
    
    
    #################################################################
    ### DETERMINE PLOT TITLE BASED ON THE DATA ###
    #################################################################
    if 'ACARS' in clean_data['site_info']['source']:
        top_title = f"ACARS AIRCRAFT OBSERVATION VERTICAL PROFILE"
        left_title = f"VALID: {clean_data['site_info']['valid-time'][1]}-{clean_data['site_info']['valid-time'][2]}-{clean_data['site_info']['valid-time'][0]} {clean_data['site_info']['valid-time'][3]}Z"
        right_title = f"{clean_data['site_info']['site-id']} - {clean_data['site_info']['site-name']} | {clean_data['site_info']['site-latlon'][0]}, {clean_data['site_info']['site-latlon'][1]}    " 

    elif 'BUFKIT' in clean_data['site_info']['source']:
        top_title = f"BUFKIT MODEL FORECAST PROFILE | {clean_data['site_info']['run-time'][3]}Z {clean_data['site_info']['model']} {clean_data['site_info']['fcst-hour']}"
        left_title = f" VALID: {clean_data['site_info']['valid-time'][1]}/{clean_data['site_info']['valid-time'][2]}/{clean_data['site_info']['valid-time'][0]} {clean_data['site_info']['valid-time'][3]}Z"
        right_title = f"{clean_data['site_info']['site-id']} - {clean_data['site_info']['site-name']} | {clean_data['site_info']['site-latlon'][0]}, {clean_data['site_info']['site-latlon'][1]}    " 

    elif 'RAOB' in clean_data['site_info']['source']:
        top_title = "RAOB OBSERVED VERTICAL PROFILE"
        left_title = f"VALID: {clean_data['site_info']['valid-time'][1]}-{clean_data['site_info']['valid-time'][2]}-{clean_data['site_info']['valid-time'][0]} {clean_data['site_info']['valid-time'][3]}Z"
        right_title = f"{clean_data['site_info']['site-id']} - {clean_data['site_info']['site-name']} | {clean_data['site_info']['site-latlon'][0]}, {clean_data['site_info']['site-latlon'][1]}    " 


    elif 'REANALYSIS' in clean_data['site_info']['source']:
        top_title = f"MODEL REANALYSIS VERTICAL PROFILE | {clean_data['site_info']['valid-time'][3]}Z {clean_data['site_info']['model']} {clean_data['site_info']['fcst-hour']}"
        left_title = f"{clean_data['site_info']['run-time'][3]}Z {clean_data['site_info']['model']} {clean_data['site_info']['fcst-hour']} | VALID: {clean_data['site_info']['valid-time'][1]}/{clean_data['site_info']['valid-time'][2]}/{clean_data['site_info']['valid-time'][0]} {clean_data['site_info']['valid-time'][3]}Z"
        right_title = f"{clean_data['site_info']['site-latlon'][0]}, {clean_data['site_info']['site-latlon'][1]}    "    
        
    else:
        top_title = clean_data['site_info']['source']
        left_title = f"VALID: {clean_data['site_info']['valid-time'][1]}-{clean_data['site_info']['valid-time'][2]}-{clean_data['site_info']['valid-time'][0]} {clean_data['site_info']['valid-time'][3]}Z"
        right_title = f"{clean_data['site_info']['site-id']} - {clean_data['site_info']['site-name']} | {clean_data['site_info']['site-latlon'][0]}, {clean_data['site_info']['site-latlon'][1]}    " 
        
    ################################################################    
        
    
    #########################################################################
    ################################ SKEW-T ################################# 
    #########################################################################

    

    #################################################################
    ### CREATE FIGURE ###
    #################################################################
    fig = plt.figure(figsize=(18,12), linewidth=10, edgecolor=brdr_clr)                             
    skew = SkewT(fig, rotation=45, rect=(0, 0, 0.6, 1)) 
    # Change to adjust data limits and give it a semblance of what we want
    skew.ax.set_box_aspect(1)
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
    # Set some better labels than the default to increase readability
    skew.ax.set_xlabel(str.upper(f'Temperature ({T.units:~P})'), weight='bold', color=gen_txt_clr)
    skew.ax.set_ylabel(str.upper(f'Pressure ({p.units:~P})'), weight='bold', color=gen_txt_clr)
    # Set the facecolor of the Skew Object and the Figure to white
    skew.ax.spines["top"].set_color(brdr_clr)
    skew.ax.spines["left"].set_color(brdr_clr)
    skew.ax.spines["right"].set_color(brdr_clr)
    skew.ax.spines["bottom"].set_color(brdr_clr)
    skew.ax.spines["bottom"].set_color(brdr_clr)
    plt.tick_params(colors=gen_txt_clr)
    plt.tick_params(colors=gen_txt_clr)
    # Define background colors
    fig.set_facecolor(bckgrnd_clr)         
    skew.ax.set_facecolor(bckgrnd_clr)  
    # Here we can use some basic math and python functionality to make a cool
    # shaded isotherm pattern. 
    x1 = np.linspace(-100, 40, 8)                                                          
    x2 = np.linspace(-90, 50, 8)                                                         
    y = [1100, 50]                                                                      
    for i in range(0, 8):              
        skew.shade_area(y=y, x1=x1[i], x2=x2[i], color='gray', alpha=shade_alpha, zorder=1)   
    #################################################################
        

        
    #################################################################
    ### PLOT SKEW T LINES ###
    #################################################################
    # Plot relevent Skew-T lines
    skew.ax.axvline(0 * units.degC, linestyle='--', color='blue', alpha=0.3)
    skew.ax.axvline(-20 * units.degC, linestyle='--', color='blue', alpha=0.3)
    skew.plot_dry_adiabats(color=skw_ln_clr, linewidth=0.5, alpha=0.7)   
    skew.plot_moist_adiabats(color=skw_ln_clr, linewidth=0.5, alpha=0.7)
    skew.plot_mixing_lines(color=skw_ln_clr, linewidth=0.2, alpha=0.7)
    # add basic temperature lines
    twline = skew.plot(p, general['wet_bulb'], '#3d8aff', linewidth=1, label='WETBULB TEMP', alpha=0.3)
    tvline = skew.plot(p, general['virt_temp'], '#0b6431', linestyle='--', linewidth=1, label='VIRTUAL TEMP', alpha=0.3)  
    tdline = skew.plot(p, Td, td_color, linewidth=4, label='DEWPOINT')
    tline1 = skew.plot(p, T, 'red', linewidth=4, label='TEMPERATURE')  
    
    # add parcel lines and CAPE shading
    if thermo['sbcape'] > 0:
        sbparcelline = skew.plot(thermo['sbP_trace'], thermo['sbT_trace'], linestyle='--',
                                 linewidth=1, color='red', alpha=1, label='SB PARCEL')    
    if thermo['mlcape'] > 0:
        mlparcelline = skew.plot(thermo['mlP_trace'], thermo['mlT_trace'], color='orangered', linestyle='--',
                                 linewidth=1, alpha=1, label='ML PARCEL') 
    if thermo['mucape'] > 0:
        muparcelline = skew.plot(thermo['muP_trace'], thermo['muT_trace'], color='orange', linestyle='--',  
                                 linewidth=1, alpha=1, label='MU PARCEL')

    skew.plot(thermo['dparcel_p'], thermo['dparcel_T'], linestyle='--',linewidth=0.7, color='purple', 
              alpha=0.8, label='DWNDRFT PARCEL')
    #################################################################
    
    
    
    
    #################################################################
    ### PLOT SKEW T WIND BARBS ###
    #################################################################
    interval = np.logspace(2.113, 3, 40) * units.hPa
    idx = mpcalc.resample_nn_1d(p, interval) 
    # create blank barbs for small dot at the start of each actual barb
    blank_len = len(u[idx])         
    blank = np.zeros(blank_len)
    skew.plot_barbs(pressure=p[idx],u=blank,v=blank,xloc=0.955,color=barb_clr,fill_empty=True,
                    sizes=dict(emptybarb=0.075, width=0.18, height=0.4))
    skew.plot_barbs(pressure=p[idx], u=u[idx], v=v[idx],xloc=0.955,color=barb_clr,fill_empty=True,
                    sizes=dict(emptybarb=0.075, width=0.18, height=0.4), length=7)
    # Draw line underneath wind barbs
    line = mlines.Line2D([0.955, 0.955], [0.01,0.95],color=barb_clr,linewidth=0.5,
                         transform=skew.ax.transAxes,clip_on=False,zorder=1)
    skew.ax.add_line(line) 
    #################################################################

    
    
    #################################################################
    ### PLOT SKEW T ANNOTATIONS ###
    #################################################################    
    hgt_lvls =[]
    for key in intrp['hgt_lvls'].keys():
        hgt_lvls.append(intrp['hgt_lvls'][key])
    hgt_lvls.pop(0) 

    for key in hgt_lvls[1::4]:
        trans, _, _ = skew.ax.get_yaxis_text1_transform(0)
        skew.ax.text(0.048, intrp['pINTRP'][key], f"{int(intrp['zINTRP'][key]/1000)}km", 
                     fontsize=13, transform=trans, color=gen_txt_clr, alpha=0.6, weight='bold')   

    sfc = mpcalc.height_to_pressure_std(general['elevation']*units.m)
    skew.ax.text(0.048, p[0], '-SFC-', fontsize=13, transform=trans, color=gen_txt_clr, alpha=0.6, weight='bold') # plot 'SFC' @ surface pressure
    
    
    # SFC TEMPERATURE AND DEWPOINT ANNOTATIONS---------------------------------------------
    T_degF = np.round(T.to(units.degF), 1)
    T_degF_label = '{}°F'.format(int(T_degF[0].magnitude))                             
    plt.annotate(T_degF_label, (T[0], p[0]), textcoords="offset points", xytext=(16,-15),
                     fontsize=12, color='red', weight='bold', alpha=0.7, ha='center')   
    Td_degF = np.round(Td.to(units.degF), 1) 
    Td_degF_label = '{}°F'.format(int(Td_degF[0].magnitude))                             
    plt.annotate(Td_degF_label,(Td[0], p[0]),textcoords="offset points",xytext=(-16,-15), 
                     fontsize=12, color=td_color, weight='bold', alpha=0.7, ha='center') 

    # PARCEL HEIGHT ANNOTATIONS------------------------------------------------------------- 
    plt.text((0.82), (thermo['sb_lcl_p']), "←SBLCL", weight='bold',color='gray',
             alpha=0.9, fontsize=15, transform=skew.ax.get_yaxis_transform(), clip_on=True)
    if ma.is_masked(thermo['mu_lfc_z']) == True:
        plt.text((0.82), (thermo['sb_lfc_p']), "←SBLFC", weight='bold',color='gray',
                 alpha=0.9, fontsize=15, transform=skew.ax.get_yaxis_transform(), clip_on=True)
        plt.text((0.82), (thermo['sb_el_p']), "←SBEL", weight='bold',color='gray', 
                 alpha=0.9, fontsize=15, transform=skew.ax.get_yaxis_transform(), clip_on=True)
    else: 
        plt.text((0.82), (thermo['mu_lfc_p']), "←MULFC", weight='bold',color='gray',
                 alpha=0.9, fontsize=15, transform=skew.ax.get_yaxis_transform(), clip_on=True)
        plt.text((0.82), (thermo['mu_el_p']), "←MUEL", weight='bold',color='gray',
                 alpha=0.9, fontsize=15, transform=skew.ax.get_yaxis_transform(), clip_on=True)

    #FRREZING POINT ANNOTATION--------------------------------------------------------------
    if ma.is_masked(general['frz_pt_z']) == False:
        if general['frz_pt_z'] >= 50*units.m:
            plt.text((0.82), (general['frz_pt_p']), "←FRZ", weight='bold',color='cornflowerblue',  
                     alpha=0.3, fontsize=13.5, transform=skew.ax.get_yaxis_transform(), clip_on=True)
    #################################################################
           
        
        
        
    #########################################################################
    ############################# HODOGRAPH ################################# 
    #########################################################################
    
            
    #################################################################
    ### DEFINE HODOGRAPH BOUNDS ###
    #################################################################
    # restructure u and v, p, ws and z data arrays based on corrected u and v arrays and hodo_layer depth
    # define x and y min/max values from 'cleaned' and restructured u and v arrays
    x_min = intrp['uINTRP'][0:90].min()
    y_min = intrp['vINTRP'][0:90].min()
    x_max = intrp['uINTRP'][0:90].max()
    y_max = intrp['vINTRP'][0:90].max()
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


    # then we can create the metpy Hodograph
    hodo_ax = plt.axes((0.52, 0.45, 0.5, 0.5))
    h = Hodograph(hodo_ax, component_range=160.)
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
        h.ax.annotate(str(i),(i,0),xytext=(0,2),textcoords='offset pixels',clip_on=True,fontsize=10,weight='bold',alpha=0.2,zorder=0, color=gen_txt_clr)
    for i in range(10,130,20):
        h.ax.annotate(str(i),(0,i),xytext=(0,2),textcoords='offset pixels',clip_on=True,fontsize=10,weight='bold',alpha=0.2,zorder=0, color=gen_txt_clr)
    for i in range(10,130,20):
        h.ax.annotate(str(i),(-i,0),xytext=(0,2),textcoords='offset pixels',clip_on=True,fontsize=10,weight='bold',alpha=0.2,zorder=0, color=gen_txt_clr)
    for i in range(10,130,20):
        h.ax.annotate(str(i),(0,-i),xytext=(0,2),textcoords='offset pixels',clip_on=True,fontsize=10,weight='bold',alpha=0.2,zorder=0, color=gen_txt_clr)

    h.plot(intrp['uINTRP'][intrp['hgt_lvls']['h05']],intrp['vINTRP'][intrp['hgt_lvls']['h05']],marker='.', 
           color='white', alpha=1, markersize=20, clip_on=True, zorder=5)
    h.ax.annotate(str('.5'),(intrp['uINTRP'][intrp['hgt_lvls']['h05']],intrp['vINTRP'][intrp['hgt_lvls']['h05']]),
                  weight='bold', fontsize=8, color='black',xytext=(0,-3.2),textcoords='offset pixels',horizontalalignment='center',clip_on=True, zorder=6) 

    hgt_lvls = [] 
    for key in intrp['hgt_lvls'].keys():
        hgt_lvls.append(intrp['hgt_lvls'][key])
    hgt_lvls.pop(0) 

    for i in hgt_lvls[1::2]:
        h.plot(intrp['uINTRP'][i],intrp['vINTRP'][i], marker='.', color='white', alpha=1, markersize=20, zorder=5)
        h.ax.annotate(str(int(round(intrp['zINTRP'][i]/1000,0))),(intrp['uINTRP'][i],intrp['vINTRP'][i]), 
                      weight='bold', fontsize=8, color='black',xytext=(0,-3.2),textcoords='offset pixels',horizontalalignment='center',clip_on=True, zorder=6) 
    #################################################################
    
    
    #################################################################
    ### PLOT HODOGRAPH LINE ###
    #################################################################
    hodo_color = ['purple','red','darkorange','gold','#fff09f']

    h.ax.plot(intrp['uINTRP'][0:10+1],   intrp['vINTRP'][0:10+1],   color=hodo_color[0], linewidth=5, clip_on=True)
    h.ax.plot(intrp['uINTRP'][10:30+1],  intrp['vINTRP'][10:30+1],  color=hodo_color[1], linewidth=5, clip_on=True)
    h.ax.plot(intrp['uINTRP'][30:60+1],  intrp['vINTRP'][30:60+1],  color=hodo_color[2], linewidth=5, clip_on=True)
    h.ax.plot(intrp['uINTRP'][60:90+1],  intrp['vINTRP'][60:90+1],  color=hodo_color[3], linewidth=5, clip_on=True)
    h.ax.plot(intrp['uINTRP'][90:120+1], intrp['vINTRP'][90:120+1], color=hodo_color[4], linewidth=5, clip_on=True) 
    #################################################################
    
    
    
    #################################################################
    ### ADD HODOGRAPH ANNOTATION ###
    #################################################################
    if ma.is_masked(kinem['sm_rm']) == False:
        # BUNKERS STORM MOTION
        h.ax.text((kinem['sm_rm'][0]+0.5), (kinem['sm_rm'][1]-0.5), 'RM', weight='bold', ha='left', fontsize=14, alpha=0.9, color=gen_txt_clr)
        h.ax.text((kinem['sm_lm'][0]+0.5), (kinem['sm_lm'][1]-0.5), 'LM', weight='bold', ha='left', fontsize=14, alpha=0.9, color=gen_txt_clr)
        h.ax.text((kinem['sm_mw'][0]+0.5), (kinem['sm_mw'][1]-0.5), 'MW', weight='bold', ha='left', fontsize=14, alpha=0.9, color=gen_txt_clr)
        h.ax.arrow(0,0,kinem['sm_u']-0.3, kinem['sm_v']-0.3, linewidth=3, color=gen_txt_clr, alpha=0.2, 
                   label='Bunkers RM Vector', length_includes_head=True, head_width=0.5)
        # DEVIANT TORNADO MOTION
        h.ax.text(kinem['dtm'][0], (kinem['dtm'][1] + 2), 'DTM', weight='bold', fontsize=10, color='brown', ha='center')
        h.plot(kinem['dtm'][0], kinem['dtm'][1], marker='v', color='brown', markersize=8, alpha=0.8, ls='', label='DEVIANT TORNADO MOTION')
    
    # EFFECTIVE INFLOW LAYER
    ebot = h.ax.plot((kinem['sm_u'], intrp['uINTRP'][find_nearest(p,kinem['eil'][0])]), (kinem['sm_v'], intrp['vINTRP'][find_nearest(p,kinem['eil'][0])]),  
                     linestyle='--', linewidth=2.3, alpha=0.5, color='lightblue', label='Effective Inflow Layer')
    etop = h.ax.plot((kinem['sm_u'], intrp['uINTRP'][find_nearest(p,kinem['eil'][1])]), (kinem['sm_v'], intrp['vINTRP'][find_nearest(p,kinem['eil'][1])]),  
                     linestyle='--', linewidth=2.3, alpha=0.5, color='lightblue')
    # EFFECTIVE INFLOW LAYER SRH FILL
    fill_srh = h.ax.fill(np.append(intrp['uINTRP'][find_nearest(p, kinem['eil'][0]):find_nearest(p, kinem['eil'][1])+1], kinem['sm_u']), 
                         np.append(intrp['vINTRP'][find_nearest(p, kinem['eil'][0]):find_nearest(p, kinem['eil'][1])+1], kinem['sm_v']),
                         'lightblue',alpha=0.1, label='EIL SRH')
    
    h.ax.text(kinem['mcs'][0], kinem['mcs'][1], 'UP', weight='bold', fontsize=12, color='orange', ha='center', alpha=0.5, clip_on=True)
    h.ax.text(kinem['mcs'][2], kinem['mcs'][3], 'DN', weight='bold', fontsize=12, color='orange', ha='center', alpha=0.5, clip_on=True)
    
    try:
        keys = ['sm_rm', 'sm_lm', 'sm_mw', 'dtm'] 

        speeds = []
        directions = []

        for key in keys:
            speeds.append(mpcalc.wind_speed(kinem[key][0]*units.kts,kinem[key][1]*units.kts).m) 
            directions.append(mpcalc.angle_to_direction((np.degrees(np.arctan2(kinem[key][0],kinem[key][1]))), full=False, level=3))

        speeds.append(mpcalc.wind_speed(kinem['mcs'][0]*units.kts,kinem['mcs'][1]*units.kts).m) 
        directions.append(mpcalc.angle_to_direction((np.degrees(np.arctan2(kinem['mcs'][0],kinem['mcs'][1]))), full=False, level=3))   

        speeds.append(mpcalc.wind_speed(kinem['mcs'][2]*units.kts,kinem['mcs'][3]*units.kts).m) 
        directions.append(mpcalc.angle_to_direction((np.degrees(np.arctan2(kinem['mcs'][2],kinem['mcs'][3]))), full=False, level=3)) 
    
        #plot Bunkers Storm Motion & DTM Data in box on Hodograph
        plt.figtext(0.61, 0.9455, 
                    f' RM: {directions[0]} @ {mag(speeds[0])} kts\n LM: {directions[1]} @ {mag(speeds[1])} kts\n MW: {directions[2]} @ {mag(speeds[2])} kts\n'+
                    f' DTM: {directions[3]} @ {mag(speeds[3])} kts\n US: {directions[4]} @ {mag(speeds[4])} kts\n DS: {directions[5]} @ {mag(speeds[5])} kts', 
                    color=gen_txt_clr,fontsize=10, verticalalignment='top', linespacing=2.2, alpha=0.6)  
    except IndexError:
        pass
    ################################################################
    
    
    
    
    
    #########################################################################
    ############################# PLOT TEXT ################################# 
    #########################################################################
    
    fig.patches.extend([plt.Rectangle((0.603, 0.05), 0.334, 0.37,
                                      edgecolor=brdr_clr, facecolor=bckgrnd_clr, linewidth=1, alpha=1,
                                      transform=fig.transFigure, figure=fig)])



    # THERMODYNAMICS 
    plt.figtext( 0.62, 0.37,  f"SBCAPE: ", weight='bold', fontsize=15, color=gen_txt_clr, ha='left')
    plt.figtext( 0.75, 0.37,  f"{mag(thermo['sbcape'])} J/kg", weight='bold', fontsize=15, color='orangered', ha='right')
    plt.figtext( 0.62, 0.34,  f"SBCIN: ", weight='bold', fontsize=15, color=gen_txt_clr, ha='left')
    plt.figtext( 0.75, 0.34,  f"{mag(thermo['sbcin'])} J/kg", weight='bold', fontsize=15, color='lightblue', ha='right')
    plt.figtext( 0.62, 0.29,  f"MLCAPE: ", weight='bold', fontsize=15, color=gen_txt_clr, ha='left')
    plt.figtext( 0.75, 0.29,  f"{mag(thermo['mlcape'])} J/kg", weight='bold', fontsize=15, color='orangered', ha='right')
    plt.figtext( 0.62, 0.26,  f"MLCIN: ", weight='bold', fontsize=15, color=gen_txt_clr, ha='left')
    plt.figtext( 0.75, 0.26,  f"{mag(thermo['mlcin'])} J/kg", weight='bold', fontsize=15, color='lightblue', ha='right')
    plt.figtext( 0.62, 0.21,  f"MUCAPE: ", weight='bold', fontsize=15, color=gen_txt_clr, ha='left')
    plt.figtext( 0.75, 0.21,  f"{mag(thermo['mucape'])} J/kg", weight='bold', fontsize=15, color='orangered', ha='right')
    plt.figtext( 0.62, 0.18,  f"MUCIN: ", weight='bold', fontsize=15, color=gen_txt_clr, ha='left')
    plt.figtext( 0.75, 0.18,  f"{mag(thermo['mucin'])} J/kg", weight='bold', fontsize=15, color='lightblue', ha='right')
    plt.figtext( 0.62, 0.13,  f"3km MUCAPE: ", weight='bold', fontsize=15, color=gen_txt_clr, ha='left')
    plt.figtext( 0.75, 0.13,  f"{mag(thermo['sb3cape'])}", weight='bold', fontsize=15, color='orangered', ha='right')
    plt.figtext( 0.62, 0.10,  f"DCAPE:", weight='bold', fontsize=15, color=gen_txt_clr, ha='left')
    plt.figtext( 0.75, 0.10,  f"{mag(thermo['dcape'])}", weight='bold', fontsize=15, color='orangered', ha='right')
    plt.figtext( 0.62, 0.07,  f"ECAPE: ", weight='bold', fontsize=15, color=gen_txt_clr, ha='left')
    plt.figtext( 0.75, 0.07,  f"{mag(thermo['ecape'])} J/kg", weight='bold', fontsize=15, color='orangered', ha='right')
                                 
    # KINEMATICS 
    met_per_sec = (units.m*units.m)/(units.sec*units.sec)
    plt.figtext( 0.77, 0.37,  f"0-500m SRH: ", weight='bold', fontsize=15, color=gen_txt_clr, ha='left')
    plt.figtext( 0.92, 0.37,  f"{mag(kinem['srh_0_to_500'])* met_per_sec:~P}", weight='bold', fontsize=15, color='deepskyblue', ha='right')
    plt.figtext( 0.77, 0.34,  f"0-500m SHEAR: ", weight='bold', fontsize=15, color=gen_txt_clr, ha='left')
    plt.figtext( 0.92, 0.34,  f"{mag(kinem['shear_0_to_500'])} kts", weight='bold', fontsize=15, color='mediumslateblue', ha='right')
    plt.figtext( 0.77, 0.29,  f"0-1km SRH: ", weight='bold', fontsize=15, color=gen_txt_clr, ha='left')
    plt.figtext( 0.92, 0.29,  f"{mag(kinem['srh_0_to_1000'])* met_per_sec:~P}", weight='bold', fontsize=15, color='deepskyblue', ha='right')
    plt.figtext( 0.77, 0.26,  f"0-1km SHEAR: ", weight='bold', fontsize=15, color=gen_txt_clr, ha='left')
    plt.figtext( 0.92, 0.26,  f"{mag(kinem['shear_0_to_1000'])} kts", weight='bold', fontsize=15, color='mediumslateblue', ha='right')
    plt.figtext( 0.77, 0.21,  f"1-3km SRH: ", weight='bold', fontsize=15, color=gen_txt_clr, ha='left')
    plt.figtext( 0.92, 0.21,  f"{mag(kinem['srh_1_to_3000'])* met_per_sec:~P}", weight='bold', fontsize=15, color='deepskyblue', ha='right')
    plt.figtext( 0.77, 0.18,  f"1-3km SHEAR: ", weight='bold', fontsize=15, color=gen_txt_clr, ha='left')
    plt.figtext( 0.92, 0.18,  f"{mag(kinem['shear_1_to_3000'])} kts", weight='bold', fontsize=15, color='mediumslateblue', ha='right')
    plt.figtext( 0.77, 0.13,  f"3-6km SRH: ", weight='bold', fontsize=15, color=gen_txt_clr, ha='left')
    plt.figtext( 0.92, 0.13,  f"{mag(kinem['srh_3_to_6000'])* met_per_sec:~P}", weight='bold', fontsize=15, color='deepskyblue', ha='right')
    plt.figtext( 0.77, 0.10,  f"3-6km SHEAR: ", weight='bold', fontsize=15, color=gen_txt_clr, ha='left')
    plt.figtext( 0.92, 0.10,  f"{mag(kinem['shear_3_to_6000'])} kts", weight='bold', fontsize=15, color='mediumslateblue', ha='right')
    plt.figtext( 0.77, 0.07,  f"0-500m SWV: ", weight='bold', fontsize=15, color=gen_txt_clr, ha='left')
    plt.figtext( 0.92, 0.07,  f"{mag(kinem['swv_perc_0_to_500'])} %", weight='bold', fontsize=15, color='cornflowerblue', ha='right')


    
    #########################################################################
    ############################# PLOT EXTRAS ###############################
    #########################################################################
    #legend
    skewleg = skew.ax.legend(loc='upper right', prop = { "size": 10 }, framealpha=0.3)

    # PLOT TITLES ----------------------------------------------------------------------------------
    plt.figtext( 0.00, 0.96, left_title, ha='left', fontsize=20, color=gen_txt_clr)
    plt.figtext( 0.96, 0.96, f'{right_title}   ', ha='right', fontsize=20, color=gen_txt_clr)
    plt.figtext( 0.00, 0.986, top_title, ha='left', weight='bold', fontsize=27, color=gen_txt_clr) 
    plt.figtext( 0.955, 0.03,  f'SOUNDERPY VERTICAL PROFILE ANALYSIS TOOL | (C) KYLE J GILLETT          ', ha='right', color='cornflowerblue', alpha=0.8, weight='bold', fontsize=10)

    # sounderpy logo
    img = Image.open(urlopen('https://user-images.githubusercontent.com/100786530/251580013-2e9477c9-e36a-4163-accb-fe46780058dd.png'))
    #                  side-side  up-down  size   size
    imgax = fig.add_axes([-0.01, 0.84, 0.1, 0.1], anchor='SE')
    imgax.imshow(img)
    imgax.axis('off')
    
    print('> COMPLETE --------')
    elapsed_time = time.time() - st
    print('> RUNTIME:', time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
    
    plt.tight_layout()
    
    return plt












#########################################################################
############################ HODOGRAPH ##################################

def __full_hodograph(clean_data, dark_mode, sr_hodo):


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
    T  = clean_data['T']
    Td = clean_data['Td']
    p  = clean_data['p']
    z  = clean_data['z']
    u  = clean_data['u']
    v  = clean_data['v']
    wd = mpcalc.wind_direction(u, v)
    ws = mpcalc.wind_speed(u, v) 
    
    # calculate other sounding parameters using SounderPy Calc
    general, thermo, kinem, intrp = sounding_params(clean_data).calc()
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
    ### DETERMINE PLOT TITLE BASED ON THE DATA ###
    #################################################################
    if 'ACARS' in clean_data['site_info']['source']:
        top_title = f"ACARS AIRCRAFT OBSERVATION {hodo_title}"
        left_title = f"VALID: {clean_data['site_info']['valid-time'][1]}-{clean_data['site_info']['valid-time'][2]}-{clean_data['site_info']['valid-time'][0]} {clean_data['site_info']['valid-time'][3]}Z"
        right_title = f"{clean_data['site_info']['site-id']} - {clean_data['site_info']['site-name']} | {clean_data['site_info']['site-latlon'][0]}, {clean_data['site_info']['site-latlon'][1]}    " 

    elif 'BUFKIT' in clean_data['site_info']['source']:
        top_title = f"BUFKIT MODEL FORECAST {hodo_title} | {clean_data['site_info']['run-time'][3]}Z {clean_data['site_info']['model']} {clean_data['site_info']['fcst-hour']}"
        left_title = f"VALID: {clean_data['site_info']['valid-time'][1]}/{clean_data['site_info']['valid-time'][2]}/{clean_data['site_info']['valid-time'][0]} {clean_data['site_info']['valid-time'][3]}Z"
        right_title = f"{clean_data['site_info']['site-id']} - {clean_data['site_info']['site-name']} | {clean_data['site_info']['site-latlon'][0]}, {clean_data['site_info']['site-latlon'][1]}    " 

    elif 'RAOB' in clean_data['site_info']['source']:
        top_title = f"RAOB OBSERVED {hodo_title}"
        left_title = f"VALID: {clean_data['site_info']['valid-time'][1]}-{clean_data['site_info']['valid-time'][2]}-{clean_data['site_info']['valid-time'][0]} {clean_data['site_info']['valid-time'][3]}Z"
        right_title = f"{clean_data['site_info']['site-id']} - {clean_data['site_info']['site-name']} | {clean_data['site_info']['site-latlon'][0]}, {clean_data['site_info']['site-latlon'][1]}    " 


    elif 'REANALYSIS' in clean_data['site_info']['source']:
        top_title = f"MODEL REANALYSIS {hodo_title} | {clean_data['site_info']['valid-time'][3]}Z {clean_data['site_info']['model']} {clean_data['site_info']['fcst-hour']}"
        left_title = f"{clean_data['site_info']['run-time'][3]}Z {clean_data['site_info']['model']} {clean_data['site_info']['fcst-hour']} | VALID: {clean_data['site_info']['valid-time'][1]}/{clean_data['site_info']['valid-time'][2]}/{clean_data['site_info']['valid-time'][0]} {clean_data['site_info']['valid-time'][3]}Z"
        right_title = f"{clean_data['site_info']['site-latlon'][0]}, {clean_data['site_info']['site-latlon'][1]}    " 
        
    else:
        top_title = clean_data['site_info']['source']
        left_title = f"VALID: {clean_data['site_info']['valid-time'][1]}-{clean_data['site_info']['valid-time'][2]}-{clean_data['site_info']['valid-time'][0]} {clean_data['site_info']['valid-time'][3]}Z"
        right_title = f"{clean_data['site_info']['site-id']} - {clean_data['site_info']['site-name']} | {clean_data['site_info']['site-latlon'][0]}, {clean_data['site_info']['site-latlon'][1]}    " 
        
    ################################################################    
    
    
    
    #################################################################
    ### DEFINE HODOGRAPH BOUNDS ###
    #################################################################
    if z.max().m > 10000: 
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

    h.plot(intrp['uINTRP'][intrp['hgt_lvls']['h05']],intrp['vINTRP'][intrp['hgt_lvls']['h05']],marker='.', 
           color='white', alpha=1, markersize=30, clip_on=True, zorder=5)
    h.ax.annotate(str('.5'),(intrp['uINTRP'][intrp['hgt_lvls']['h05']],intrp['vINTRP'][intrp['hgt_lvls']['h05']]),
                  weight='bold', fontsize=12, color='black',xytext=(0,-3.2),textcoords='offset pixels',horizontalalignment='center',clip_on=True, zorder=6) 

    hgt_lvls = [] 
    for key in intrp['hgt_lvls'].keys():
        hgt_lvls.append(intrp['hgt_lvls'][key])
    hgt_lvls.pop(0) 

    for i in hgt_lvls[1::2]:
        h.plot(intrp['uINTRP'][i],intrp['vINTRP'][i], marker='.', color='white', alpha=1, markersize=30, zorder=5)
        h.ax.annotate(str(int(round(intrp['zINTRP'][i]/1000,0))),(intrp['uINTRP'][i],intrp['vINTRP'][i]), 
                      weight='bold', fontsize=12, color='black',xytext=(0,-3.2),textcoords='offset pixels',horizontalalignment='center',clip_on=True, zorder=6) 
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
            h.ax.arrow(0,0,kinem['sm_u']-0.3, kinem['sm_v']-0.3, linewidth=3, color=gen_txt_clr, alpha=0.2,label='Bunkers RM Vector', length_includes_head=True, head_width=0.5)
            # DEVIANT TORNADO MOTION
            h.ax.text(kinem['dtm'][0], (kinem['dtm'][1] + 2), 'DTM', weight='bold', fontsize=10, color='brown', ha='center')
            h.plot(kinem['dtm'][0], kinem['dtm'][1], marker='v', color='brown', markersize=8, alpha=0.8, ls='', label='DEVIANT TORNADO MOTION')
        elif sr_hodo == True:
            h.ax.text((kinem['sm_lm'][0] - kinem['sm_u'] +0.5), (kinem['sm_lm'][1] - kinem['sm_v'] -0.5), 'LM', weight='bold', ha='left', fontsize=14, alpha=0.9, color=gen_txt_clr)
            h.ax.text((kinem['sm_mw'][0] - kinem['sm_u'] +0.5), (kinem['sm_mw'][1] - kinem['sm_v'] -0.5), 'MW', weight='bold', ha='left', fontsize=14, alpha=0.9, color=gen_txt_clr)
            # DEVIANT TORNADO MOTION
            h.ax.text(kinem['dtm'][0] - kinem['sm_u'], (kinem['dtm'][1] - kinem['sm_v'] + 2), 'DTM', weight='bold', fontsize=10, color='brown', ha='center')
            h.plot(kinem['dtm'][0] - kinem['sm_u'], kinem['dtm'][1] - kinem['sm_v'], marker='v', color='brown', markersize=8, alpha=0.8, ls='', label='DEVIANT TORNADO MOTION')
          
        
    # EFFECTIVE INFLOW LAYER
    if sr_hodo == False:
        ebot = h.ax.plot((kinem['sm_u'], intrp['uINTRP'][find_nearest(p,kinem['eil'][0])]), (kinem['sm_v'], intrp['vINTRP'][find_nearest(p,kinem['eil'][0])]),  
                     linestyle='--', linewidth=2.3, alpha=0.5, color='lightblue', label='Effective Inflow Layer')
        etop = h.ax.plot((kinem['sm_u'], intrp['uINTRP'][find_nearest(p,kinem['eil'][1])]), (kinem['sm_v'], intrp['vINTRP'][find_nearest(p,kinem['eil'][1])]),  
                     linestyle='--', linewidth=2.3, alpha=0.5, color='lightblue')
        # EFFECTIVE INFLOW LAYER SRH FILL
        fill_srh = h.ax.fill(np.append(intrp['uINTRP'][find_nearest(p, kinem['eil'][0]):find_nearest(p, kinem['eil'][1])+1], kinem['sm_u']), 
                             np.append(intrp['vINTRP'][find_nearest(p, kinem['eil'][0]):find_nearest(p, kinem['eil'][1])+1], kinem['sm_v']),
                             'lightblue',alpha=0.1, label='EIL SRH')
        
    elif sr_hodo == True: 
        ebot = h.ax.plot((0, intrp['uINTRP'][find_nearest(p,kinem['eil'][0])]), (0, intrp['vINTRP'][find_nearest(p,kinem['eil'][0])]),  
                     linestyle='--', linewidth=2.3, alpha=0.5, color='lightblue', label='Effective Inflow Layer')
        etop = h.ax.plot((0, intrp['uINTRP'][find_nearest(p,kinem['eil'][1])]), (0, intrp['vINTRP'][find_nearest(p,kinem['eil'][1])]),  
                     linestyle='--', linewidth=2.3, alpha=0.5, color='lightblue')
        # EFFECTIVE INFLOW LAYER SRH FILL
        fill_srh = h.ax.fill(np.append(intrp['uINTRP'][find_nearest(p, kinem['eil'][0]):find_nearest(p, kinem['eil'][1])+1], 0), 
                             np.append(intrp['vINTRP'][find_nearest(p, kinem['eil'][0]):find_nearest(p, kinem['eil'][1])+1], 0),
                             'lightblue',alpha=0.1, label='EIL SRH')
    
    if sr_hodo == False:
        h.ax.text(kinem['mcs'][0], kinem['mcs'][1], 'UP', weight='bold', fontsize=12, color='orange', ha='center', alpha=0.5, clip_on=True)
        h.ax.text(kinem['mcs'][2], kinem['mcs'][3], 'DN', weight='bold', fontsize=12, color='orange', ha='center', alpha=0.5, clip_on=True)
    
    elif sr_hodo == True:
        h.ax.text(kinem['mcs'][0]- kinem['sm_u'], kinem['mcs'][1]- kinem['sm_v'], 'UP', weight='bold', fontsize=12, color='orange', ha='center', alpha=0.5, clip_on=True)
        h.ax.text(kinem['mcs'][2]- kinem['sm_u'], kinem['mcs'][3]- kinem['sm_v'], 'DN', weight='bold', fontsize=12, color='orange', ha='center', alpha=0.5, clip_on=True)
    
    try:
        keys = ['sm_rm', 'sm_lm', 'sm_mw', 'dtm'] 

        speeds = []
        directions = []

        for key in keys:
            speeds.append(mpcalc.wind_speed(kinem[key][0]*units.kts,kinem[key][1]*units.kts).m) 
            directions.append(mpcalc.angle_to_direction((np.degrees(np.arctan2(kinem[key][0],kinem[key][1]))), full=False, level=3))

        speeds.append(mpcalc.wind_speed(kinem['mcs'][0]*units.kts,kinem['mcs'][1]*units.kts).m) 
        directions.append(mpcalc.angle_to_direction((np.degrees(np.arctan2(kinem['mcs'][0],kinem['mcs'][1]))), full=False, level=3))   

        speeds.append(mpcalc.wind_speed(kinem['mcs'][2]*units.kts,kinem['mcs'][3]*units.kts).m) 
        directions.append(mpcalc.angle_to_direction((np.degrees(np.arctan2(kinem['mcs'][2],kinem['mcs'][3]))), full=False, level=3)) 
    
        #plot Bunkers Storm Motion & DTM Data in box on Hodograph
        plt.figtext(0.228, 0.86, 
                    f' RM: {directions[0]} @ {mag(speeds[0])} kts\n LM: {directions[1]} @ {mag(speeds[1])} kts\n MW: {directions[2]} @ {mag(speeds[2])} kts\n'+
                    f' DTM: {directions[3]} @ {mag(speeds[3])} kts\n US: {directions[4]} @ {mag(speeds[4])} kts\n DS: {directions[5]} @ {mag(speeds[5])} kts', 
                    color=gen_txt_clr, fontsize=10, verticalalignment='top', linespacing=2.2, alpha=0.6)  
    except IndexError:
        pass
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

    if ma.is_masked(kinem['sm_rm']) == False:
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
        
    if ma.is_masked(kinem['sm_rm']) == False: 
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

    if ma.is_masked(kinem['sm_rm']) == False:
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
    plt.figtext( 1.070, 0.78, f"{mag(kinem['swv_perc_1_to_3000'])}", fontsize=15, color='mediumslateblue', weight='bold')
    plt.figtext( 1.12, 0.78,  f"{mag_round(kinem['swv_3_to_6000'], 3)}", fontsize=15, color='mediumslateblue', weight='bold')

    plt.figtext( 0.815, 0.75, f"1-3ₖₘ:", weight='bold', fontsize=15, color=gen_txt_clr, alpha=0.9)
    plt.figtext( 0.875, 0.75, f"{mag(kinem['shear_1_to_3000'])} kt", fontsize=15, color='slateblue', weight='bold')
    plt.figtext( 0.929, 0.75,  f"{mag(kinem['srh_1_to_3000'])* met_per_sec:~P}", fontsize=15, color='slateblue', weight='bold')
    plt.figtext( 1.006,   0.75, f"{mag(kinem['srw_1_to_3000'])} kt", fontsize=15, color='slateblue', weight='bold')
    plt.figtext( 1.070, 0.75, f"{mag(kinem['swv_perc_1_to_3000'])}", fontsize=15, color='slateblue', weight='bold')
    plt.figtext( 1.12, 0.75,  f"{mag_round(kinem['swv_1_to_3000'], 3)}", fontsize=15, color='slateblue', weight='bold')

    plt.figtext( 0.815, 0.72, f"3-6ₖₘ:", weight='bold', fontsize=15, color=gen_txt_clr, alpha=0.9)
    plt.figtext( 0.875, 0.72, f"{mag(kinem['shear_3_to_6000'])} kt", fontsize=15, color='darkslateblue', weight='bold')
    plt.figtext( 0.929, 0.72,  f"{mag(kinem['srh_3_to_6000'])* met_per_sec:~P}", fontsize=15, color='darkslateblue', weight='bold')
    plt.figtext( 1.006,   0.72, f"{mag(kinem['srw_3_to_6000'])} kt", fontsize=15, color='darkslateblue', weight='bold')
    plt.figtext( 1.070, 0.72, f"{mag(kinem['swv_perc_3_to_6000'])}", fontsize=15, color='darkslateblue', weight='bold')
    plt.figtext( 1.12, 0.72,  f"{mag_round(kinem['swv_3_to_6000'], 3)}", fontsize=15, color='darkslateblue', weight='bold')
    #################################################################
    
    
    #################################################################
    ### THERMODYNAMICS ###
    #################################################################
    plt.figtext( 0.87, 0.65, 'CAPE         6CAPE         3CAPE        CIN', color=gen_txt_clr, weight='bold', fontsize=15)
    #SBCAPE
    plt.figtext( 0.815, 0.62,  f"SB:", weight='bold',   fontsize=15, color=gen_txt_clr)
    plt.figtext( 0.86, 0.62,   f"{mag(thermo['sbcape'])} J/kg",  fontsize=15, color='orangered', weight='bold')
    plt.figtext( 0.94, 0.62,   f"{mag(thermo['sb6cape'])} J/kg", fontsize=15, color='orangered', weight='bold')
    plt.figtext( 1.035, 0.62,  f"{mag(thermo['sb3cape'])} J/kg", fontsize=15, color='orangered', weight='bold')
    plt.figtext( 1.107, 0.62,  f"{mag(thermo['sbcin'])} J/kg",   fontsize=15, color='orangered', weight='bold')
    #MUCAPE
    plt.figtext( 0.815, 0.59,  f"MU:", weight='bold',   fontsize=15, color=gen_txt_clr)
    plt.figtext( 0.86, 0.59,   f"{mag(thermo['mucape'])} J/kg",  fontsize=15, color='red', weight='bold')
    plt.figtext( 0.94, 0.59,   f"{mag(thermo['mu6cape'])} J/kg", fontsize=15, color='red', weight='bold')
    plt.figtext( 1.035, 0.59,  f"{mag(thermo['mu3cape'])} J/kg", fontsize=15, color='red', weight='bold')
    plt.figtext( 1.107, 0.59,  f"{mag(thermo['mucin'])} J/kg",   fontsize=15, color='red', weight='bold')
    #MLCAPE
    plt.figtext( 0.815, 0.56,  f'ML:', weight='bold', fontsize=15, color=gen_txt_clr)
    plt.figtext( 0.86, 0.56,   f"{mag(thermo['mlcape'])} J/kg",  fontsize=15, color='darkred', weight='bold')
    plt.figtext( 0.94, 0.56,   f"{mag(thermo['ml6cape'])} J/kg", fontsize=15, color='darkred', weight='bold')
    plt.figtext( 1.035, 0.56,  f"{mag(thermo['ml3cape'])} J/kg", fontsize=15, color='darkred', weight='bold')
    plt.figtext( 1.107, 0.56,  f"{mag(thermo['mlcin'])} J/kg",   fontsize=15, color='darkred', weight='bold')
    # ECAPE
    plt.figtext( 0.86, 0.49, f"ECAPE:", weight='bold', fontsize=15, color=gen_txt_clr)
    plt.figtext( 0.92, 0.49, f" {mag(thermo['ecape'])} J/kg", fontsize=15, color='brown', alpha=0.7, weight='bold')
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
    plt.figtext( 0.23, 0.89, left_title, ha='left', weight='bold', fontsize=16, color=gen_txt_clr)
    plt.figtext( 1.20, 0.89, f'{right_title}  ', ha='right', weight='bold', fontsize=16, color=gen_txt_clr)
    plt.figtext( 0.23, 0.92, top_title, ha='left', weight='bold', fontsize=22, color=gen_txt_clr) 
    plt.figtext( 0.23, 0.94, ' ') 
    plt.figtext( 0.225, 0.09, f'SOUNDERPY VERTICAL PROFILE ANALYSIS TOOL | (C) KYLE J GILLETT 2023, CENTRAL MICHIGAN UNIVERSITY | AVAILABLE ON GITHUB & PYPI',
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

def __composite_sounding(data_list, shade_between, ls_to_use,
                    alphas_to_use, colors_to_use, lw_to_use, dark_mode): 
    
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
    
    if colors_to_use == 'none':
        colors_to_use = [
        'cornflowerblue', 'orange', 'darkgreen', 'purple', 'brown', 'cyan', 'red', 'gray', 'teal', 'pink',
        'darkcyan', 'olive', 'blue', 'coral', 'darkred', 'navy', 'magenta', 'indigo', 'darkkhaki',
        'violet', 'turquoise', 'tomato', 'sienna', 'slategray', 'green', 'mediumvioletred', 'mediumseagreen',
        'cadetblue', 'darkolivegreen', 'firebrick']
    
    
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
            lw_to_use.append(1)  
    
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

    # record process time 
    st = time.time()  

    def find_nearest(array, value):
        array = np.asarray(array)
        nearest_idx = (np.abs(array - value)).argmin()
        return nearest_idx
    
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
    
    plt.xticks(np.arange(0,0,1))
    plt.yticks(np.arange(0,0,1))
    for i in range(10,120,10):
        h.ax.annotate(str(i),(i,0),xytext=(0,2),textcoords='offset pixels',clip_on=True,
                      color=gen_txt_clr, fontsize=12,weight='bold',alpha=0.5,zorder=0)
    for i in range(10,120,10):
        h.ax.annotate(str(i),(0,i),xytext=(0,2),textcoords='offset pixels',clip_on=True,
                      color=gen_txt_clr, fontsize=12,weight='bold',alpha=0.5,zorder=0)
    for i in range(10,120,10):
        h.ax.annotate(str(i),(-i,0),xytext=(0,2),textcoords='offset pixels',clip_on=True,
                      color=gen_txt_clr, fontsize=12,weight='bold',alpha=0.5,zorder=0)
    for i in range(10,120,10):
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

        if profile['site_info']['source'] == 'MODEL REANALYSIS PROFILE':
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

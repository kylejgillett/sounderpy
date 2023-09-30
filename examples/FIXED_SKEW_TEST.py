import sounderpy as spy 
import time
from PIL import Image
from urllib.request import urlopen
import numpy as np
from numpy import loadtxt
import numpy.ma as ma
from ecape.calc import calc_ecape, _get_parcel_profile, calc_mse, calc_integral_arg, calc_lfc_height, calc_el_height
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import metpy.calc as mpcalc
from metpy.units import units
from metpy.plots import SkewT, Hodograph
#########################################################################################################



raw_data = spy.get_model_data('rap', [35.31, -97.57], '2013', '05', '20', '20')

clean_data = spy.parse_data(raw_data)



### METPY_SOUNDING() FUNCTION ###                                                            
#########################################################################################################

def __metpy_sounding(clean_data):
    
    # record process time 
    st = time.time()   
    
    '''
        this function inputs cleaned profile data through a MetPy sounding plot script. 

        COPYRIGHT
        Created by Kyle J Gillett (@wxkylegillett) 2023
    '''

    print(f'> SOUNDING PLOTTER FUNCTION --\n---------------------------------')
    
    # declare basic variables
    p = clean_data['p']
    T = clean_data['T']
    Td = clean_data['Td']
    z = clean_data['z']
    u = clean_data['u']
    v = clean_data['v'] 

    # ADD A SIMPLE PLOT TITLE---------------------------------------------------
    top_title = f"{clean_data['site_info']['source']}"


    if 'ACARS' in clean_data['site_info']['source']:
        left_title = f"VALID: {clean_data['site_info']['valid-time'][1]}-{clean_data['site_info']['valid-time'][2]}-{clean_data['site_info']['valid-time'][0]} {clean_data['site_info']['valid-time'][3]}Z"
        right_title = f"{clean_data['site_info']['site-id']} - {clean_data['site_info']['site-name']} | {clean_data['site_info']['site-latlon'][0]}, {clean_data['site_info']['site-latlon'][1]}" 

    elif 'BUFKIT' in clean_data['site_info']['source']:
        left_title = f"{clean_data['site_info']['run-time'][3]}Z {clean_data['site_info']['model']} {clean_data['site_info']['fcst-hour']} | VALID: {clean_data['site_info']['valid-time'][1]}/{clean_data['site_info']['valid-time'][2]}/{clean_data['site_info']['valid-time'][0]} {clean_data['site_info']['valid-time'][3]}Z"
        right_title = f"{clean_data['site_info']['site-id']} - {clean_data['site_info']['site-name']} | {clean_data['site_info']['site-latlon'][0]}, {clean_data['site_info']['site-latlon'][1]}" 

    elif 'OBSERVED' in clean_data['site_info']['source']:
        left_title = f"VALID: {clean_data['site_info']['valid-time'][1]}-{clean_data['site_info']['valid-time'][2]}-{clean_data['site_info']['valid-time'][0]} {clean_data['site_info']['valid-time'][3]}Z"
        right_title = f"{clean_data['site_info']['site-id']} - {clean_data['site_info']['site-name']} | {clean_data['site_info']['site-latlon'][0]}, {clean_data['site_info']['site-latlon'][1]}" 


    elif 'REANALYSIS' in clean_data['site_info']['source']:
        left_title = f"{clean_data['site_info']['run-time'][3]}Z {clean_data['site_info']['model']} {clean_data['site_info']['fcst-hour']} | VALID: {clean_data['site_info']['valid-time'][1]}/{clean_data['site_info']['valid-time'][2]}/{clean_data['site_info']['valid-time'][0]} {clean_data['site_info']['valid-time'][3]}Z"
        right_title = f"{clean_data['site_info']['site-latlon'][0]}, {clean_data['site_info']['site-latlon'][1]}" 



    # DEFINE find_nearest FUNCTION ----------------------------------------------- 
    def find_nearest(array, value):
        array = np.asarray(array)
        nearest_idx = (np.abs(array - value)).argmin()
        return nearest_idx

    # CREATE THE METPY SKEWT FIGURE ----------------------------------------------- 
    fig = plt.figure(figsize=(18,12))                             
    skew = SkewT(fig, rotation=45, rect=(0, 0, 0.6, 1)) 
    # Change to adjust data limits and give it a semblance of what we want
    skew.ax.set_box_aspect(1)
    skew.ax.set_adjustable('datalim')
    skew.ax.set_ylim(1050, 100)    
    if T[0].m <= -10:
        bounds = [-42, 42]
    elif T[0].m >= 30:
        bounds = [-22, 62]
    else: 
        bounds = [-32, 52]
    skew.ax.set_xlim(bounds[0], bounds[1])    
    # Set some better labels than the default to increase readability
    skew.ax.set_xlabel(str.upper(f'Temperature ({T.units:~P})'), weight='bold')
    skew.ax.set_ylabel(str.upper(f'Pressure ({p.units:~P})'), weight='bold')
    # Set the facecolor of the Skew Object and the Figure to white
    fig.set_facecolor('#ffffff')         
    skew.ax.set_facecolor('#ffffff')     
    # Here we can use some basic math and python functionality to make a cool
    # shaded isotherm pattern. 
    x1 = np.linspace(-100, 40, 8)                                                          
    x2 = np.linspace(-90, 50, 8)                                                         
    y = [1100, 50]                                                                      
    for i in range(0, 8):              
        skew.shade_area(y=y, x1=x1[i], x2=x2[i], color='gray', alpha=0.02, zorder=1)   

    # TEMPERATURE & DEWPOINT TRACE --------------------------------------------------
    # Plot the data using normal plotting functions, in this case using
    # log scaling in Y, as dictated by the typical meteorological plot
    # set the linewidth to 4 for increased readability. 
    # We will also add the 'label' kew word argument for our legend. 
    skew.plot(p, T, 'r', lw=4, label='TEMPERATURE')
    skew.plot(p, Td, 'g', lw=4, label='DEWPOINT')


    # WIND BARBS --------------------------------------------------------------------
    # again we can use some simple python math functionality to 'resample'
    # the wind barbs for a cleaner output with increased readability. 
    # Something like this would work.
    interval = np.logspace(2.113, 3, 40) * units.hPa
    idx = mpcalc.resample_nn_1d(p, interval) 
    # create blank barbs for small dot at the start of each actual barb
    blank_len = len(u[idx])         
    blank = np.zeros(blank_len)
    skew.plot_barbs(pressure=p[idx],u=blank,v=blank,xloc=0.955,fill_empty=True,
                    sizes=dict(emptybarb=0.075, width=0.18, height=0.4))
    skew.plot_barbs(pressure=p[idx], u=u[idx], v=v[idx],xloc=0.955,fill_empty=True,
                    sizes=dict(emptybarb=0.075, width=0.18, height=0.4), length=7)
    # Draw line underneath wind barbs
    line = mpl.lines.Line2D([0.955, 0.955], [0.01,0.95],color='black',linewidth=0.5,
                         transform=skew.ax.transAxes,clip_on=False,zorder=1)
    skew.ax.add_line(line) 


    # OTHER RELEVANT LINES -----------------------------------------------------------
    # Add the relevant special lines native to the Skew-T Log-P diagram &
    # provide basic adjustments to linewidth and alpha to increase readability
    # first we add a matplotlib axvline to highlight the 0 degree isotherm
    skew.ax.axvline(0 * units.degC, linestyle='--', color='blue', alpha=0.3)
    skew.plot_dry_adiabats(lw=1, alpha=0.3)
    skew.plot_moist_adiabats(lw=1, alpha=0.3)
    skew.plot_mixing_lines(lw=1, alpha=0.3)


    # PARCEL PROPERTIES --------------------------------------------------------------
    # Calculate full parcel profile and add to plot as black line
    # mixed layer parcel properties!
    ml_t, ml_td = mpcalc.mixed_layer(p, T, Td, depth=100 * units.hPa)
    ml_p, _, _ = mpcalc.mixed_parcel(p, T, Td, depth=100 * units.hPa)
    # most unstable parcel properties!
    mu_p, mu_t, mu_td, _ = mpcalc.most_unstable_parcel(p, T, Td, depth=100 * units.hPa)
    # Compute parcel profiles
    sb_prof = mpcalc.parcel_profile(p, T[0], Td[0]).to('degC')
    mu_prof = mpcalc.parcel_profile(p, mu_t, mu_td).to('degC')
    ml_prof = mpcalc.parcel_profile(p, ml_t, ml_td).to('degC')
    # compute CAPE & CIN
    mucape, mucin = mpcalc.most_unstable_cape_cin(p, T, Td, depth=50 * units.hPa)
    mlcape, mlcin = mpcalc.mixed_layer_cape_cin(p, T, ml_prof, depth=100 * units.hPa)
    sbcape, sbcin = mpcalc.surface_based_cape_cin(p, T, Td)

    try: 
        q = mpcalc.specific_humidity_from_dewpoint(p, Td)
        ecape = int(calc_ecape(z, p, T, q, u, v, 'most_unstable').m)
        mse_star, mse_bar = calc_mse(p, z, T, q)
        int_arg = calc_integral_arg(mse_bar, mse_star, T)
        ecape_lfc = calc_lfc_height(p, z, T, Td, parcel_func = mpcalc.mixed_parcel)
        ecape_el = calc_el_height(p, z, T, Td, parcel_func = mpcalc.mixed_parcel)
    except:
        ecape = -9999*units.joule/units.kilogram
        pass

    lcl_pressure, lcl_temperature = mpcalc.lcl(p[0], T[0], Td[0])
    lfc_pressure, lfc_temperature = mpcalc.lfc(p, T, Td, which='most_cape')
    el_pressure, el_temperature = mpcalc.el(p, T, Td, which='most_cape')

    # make minor corrections to CAPE amount to match SHARPPY (SPC) output
    sbcape = (sbcape.m + (sbcape.m*(0.10)))*units('J/kg')
    mucape = (mucape.m + (mucape.m*(0.10)))*units('J/kg')
    mlcape = (mlcape.m + (mlcape.m*(0.10)))*units('J/kg')


    # PLOT PARCEL-RELATED 'LEVELS' & SHADE CAPE/CIN ---------------------------------
    # always plot LCL
    plt.text((0.84), (lcl_pressure), "←LCL", weight='bold',color='gray',             
             alpha=0.9, fontsize=15, transform=skew.ax.get_yaxis_transform())

    # only plot LFC, EL, SB trace, & CAPE shade when SBCAPE is > 10
    if sbcape.m > 10:
        plt.text((0.84), (lfc_pressure), "←LFC", weight='bold',color='gray',          
                 alpha=0.9, fontsize=15, transform=skew.ax.get_yaxis_transform())
        plt.text((0.84), (el_pressure), "←EL", weight='bold',color='gray',             
                 alpha=0.9, fontsize=15, transform=skew.ax.get_yaxis_transform())
        skew.plot(p, sb_prof, 'red', linewidth=2, ls='--', alpha=0.6, label='SB PARCEL PATH')
        # Shade areas of CAPE and CIN
        skew.shade_cin(p, T, sb_prof, Td, alpha=0.2, label='SBCIN')
        skew.shade_cape(p, T, sb_prof, alpha=0.2, label='SBCAPE')

    # plot ML-trace when MLCAPE > 10
    if mlcape.m > 10:
        skew.plot(p, ml_prof, 'orangered', linewidth=2, ls='--', alpha=0.2, label='ML PARCEL PATH')            
        # Shade areas of CAPE and CIN
        # skew.shade_cin(p, T, ml_prof, Td, alpha=0.2, label='MLCIN')
        # skew.shade_cape(p, T, ml_prof, alpha=0.2, label='MLCAPE') 

    # plot MU-trace when MUCAPE > 10
    if mucape.m > 10:
        skew.plot(p, mu_prof, 'orange', linewidth=2, ls='--', alpha=0.4, label='MU PARCEL PATH')
        # Shade areas of CAPE and CIN
        #skew.shade_cin(p, T, mu_prof, Td, alpha=0.2, label='MUCIN')
        #skew.shade_cape(p, T, mu_prof, alpha=0.2, label='MUCAPE')


    # SURFACE TEMPERATURE & DEWPOINT ANNOTATIONS ------------------------------------
    T_degF = np.round(T.to(units.degF), 1)
    T_degF_label = '{}°F'.format(int(T_degF[0].magnitude))                             
    plt.annotate(T_degF_label, (T[0], p[0]), textcoords="offset points", xytext=(16,-15),
                     fontsize=15, color='red', weight='bold', alpha=0.7, ha='center')   
    Td_degF = np.round(Td.to(units.degF), 1) 
    Td_degF_label = '{}°F'.format(int(Td_degF[0].magnitude))                             
    plt.annotate(Td_degF_label,(Td[0], p[0]),textcoords="offset points",xytext=(-16,-15), 
                     fontsize=15, color='green', weight='bold', alpha=0.7, ha='center') 


    # FREEZING POINT ANNOTATION -----------------------------------------------------
    T_list = T.m.tolist()    
    try:
        res = list(filter(lambda i: i <= 0, T_list))[0]
        frz_pt_index = T_list.index(res) 
        frz_pt_p = p[frz_pt_index]     
        frz_pt_z = z[frz_pt_index]  
        if frz_pt_z >= 50*units.m:                                                       
            plt.text((0.84), (frz_pt_p), "←FRZ", weight='bold',color='blue',           
                     alpha=0.3, fontsize=15, transform=skew.ax.get_yaxis_transform())
    except:
        pass

    # CREATE METPY HODOGRAPH --------------------------------------------------------
    # Create a hodograph object: first we need to add an axis
    # determine max height of wind data to plot on hodograph in km (if hodo_layer = 9, 0-9km u and v are plotted)

    # remove nan values from base wind u and v component arrays to find min & max values.
    u_clean = u.magnitude[np.logical_not(np.isnan(u.magnitude))]
    v_clean = v.magnitude[np.logical_not(np.isnan(v.magnitude))]


    if z.max().m > 12000: 
        hodo_hgt = 12000*units.m
        # restructure u and v, p, ws and z data arrays based on corrected u and v arrays and hodo_layer depth
        p_hodo, u_hodo, v_hodo, z_hodo = mpcalc.get_layer(p, u, v, z, depth=hodo_hgt)
    else:
        p_hodo = p
        u_hodo = u
        v_hodo = v
        z_hodo = z

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
       

    # then we can create the metpy Hodograph
    hodo_ax = plt.axes((0.52, 0.45, 0.5, 0.5))
    h = Hodograph(hodo_ax, component_range=160.)

    # Add two seperate grid increments for a cooler look. This also
    # helps to increase readability
    h.add_grid(increment=20, ls='-', lw=1.5, alpha=0.5)
    h.add_grid(increment=10, ls='--', lw=1, alpha=0.2)

    try:
        h.ax.set_xlim(x_Minlimit, x_Maxlimit)                                  
        h.ax.set_ylim(y_Minlimit, y_Maxlimit)                             
    except:
        h.ax.set_xlim(-65,65)
        h.ax.set_ylim(-65,65)
        pass
    
    # The next few steps makes for a clean hodograph inset, removing the 
    # tick marks, tick labels and axis labels
    h.ax.set_box_aspect(1) 
    h.ax.set_yticklabels([])
    h.ax.set_xticklabels([])
    h.ax.set_xticks([])
    h.ax.set_yticks([])
    h.ax.set_xlabel(' ')
    h.ax.set_ylabel(' ')
    h.ax.set_xlim(x_Minlimit, x_Maxlimit)
    h.ax.set_ylim(y_Minlimit, y_Maxlimit)  

    # Here we can add a simple python for loop that adds tick marks to the inside 
    # of the hodograph plot to increase readability! 
    plt.xticks(np.arange(0,0,1))
    plt.yticks(np.arange(0,0,1))
    for i in range(10,130,10):
        h.ax.annotate(str(i),(i,0),xytext=(0,2),textcoords='offset pixels',clip_on=True,fontsize=12,weight='bold',alpha=0.2,zorder=0)
    for i in range(10,130,10):
        h.ax.annotate(str(i),(0,i),xytext=(0,2),textcoords='offset pixels',clip_on=True,fontsize=12,weight='bold',alpha=0.2,zorder=0)
    for i in range(10,130,10):
        h.ax.annotate(str(i),(-i,0),xytext=(0,2),textcoords='offset pixels',clip_on=True,fontsize=12,weight='bold',alpha=0.2,zorder=0)
    for i in range(10,130,10):
        h.ax.annotate(str(i),(0,-i),xytext=(0,2),textcoords='offset pixels',clip_on=True,fontsize=12,weight='bold',alpha=0.2,zorder=0)


    every = int(len(z)/3)
    for zi, ui, vi in zip(z_hodo[2 :: every], u_hodo[2 :: every], v_hodo[2 :: every]): 
        h.plot(ui,vi, marker='.', color='white', alpha=0.8, markersize=30)
        h.ax.annotate(str(int(zi.m)),(ui,vi), weight='bold', fontsize=13, color='black',xytext=(0,-3.2),textcoords='offset pixels',horizontalalignment='center',clip_on=True) 


    # plot the hodograph itself, using plot_colormapped, colored
    # by height 
    h.plot_colormapped(u_hodo, v_hodo, c=z_hodo, linewidth=6, label="0-9km WIND")

    try:
        # compute Bunkers storm motion so we can plot it on the hodograph! 
        RM, LM, MW = mpcalc.bunkers_storm_motion(p, u, v, z)
        h.ax.text((RM[0].m), (RM[1].m), 'RM', weight='bold', ha='center', fontsize=15, alpha=0.5) 
        h.ax.text((LM[0].m), (LM[1].m), 'LM', weight='bold', ha='center', fontsize=15, alpha=0.5) 
        h.ax.text((MW[0].m), (MW[1].m), 'MW', weight='bold', ha='center', fontsize=15, alpha=0.5) 
        h.ax.arrow(0,0,RM[0].m, RM[1].m, linewidth=2, color='black', alpha=0.2, label='Bunkers RM Vector', 
                   length_includes_head=True, head_width=1)

        #SHARPPY MEAN WIND FROM 0-300m
        idx_300m = find_nearest(z, z[0].m+300)
        MW_300m_u = sum(u[0:idx_300m].m)/len(u[0:idx_300m].m)
        MW_300m_v = sum(v[0:idx_300m].m)/len(v[0:idx_300m].m)

        #DTM CALC
        DTM_u = (RM[0].m+MW_300m_u)/2
        DTM_v = (RM[1].m+MW_300m_v)/2
        h.ax.text((DTM_u), (DTM_v + 1.4), 'DTM', weight='bold', fontsize=13, color='brown', ha='center')
        h.plot(DTM_u, DTM_v, marker='v', color='brown', markersize=10, alpha=0.6)

    except:
        pass

    # First we want to actually add values of data to the plot for easy viewing
    # to do this, lets first add a simple rectangle using matplotlib's 'patches' 
    # fucntionality to add some simple layout for plotting calculated parameters 
    #                                  xloc   yloc   xsize  ysize
    fig.patches.extend([plt.Rectangle((0.603, 0.05), 0.334, 0.37,
                                      edgecolor='black', facecolor='white', linewidth=1, alpha=1,
                                      transform=fig.transFigure, figure=fig)])


    # COMPUTE METPY PARAMETERS -----------------------------------------------------------
    # now lets take a moment to calculate some simple severe-weather parameters using
    # metpy's calculations 
    # here are some classic severe parameters!
    kindex = mpcalc.k_index(p, T, Td)
    total_totals = mpcalc.total_totals_index(p, T, Td)
    # Compute SRH 

    try:
        (u_storm, v_storm), *_ = mpcalc.bunkers_storm_motion(p, u, v, z)
        *_, total_helicity1 = mpcalc.storm_relative_helicity(z, u, v, depth=1 * units.km,
                                                            storm_u=u_storm, storm_v=v_storm)
        *_, total_helicity3 = mpcalc.storm_relative_helicity(z, u, v, depth=3 * units.km,
                                                            storm_u=u_storm, storm_v=v_storm)
        *_, total_helicity6 = mpcalc.storm_relative_helicity(z, u, v, depth=6 * units.km,
                                                            storm_u=u_storm, storm_v=v_storm)
    except:
        total_helicity1 = float("Nan")
        total_helicity3 = float("Nan")
        total_helicity6 = float("Nan")
        pass

    # Copmute Bulk Shear components and then magnitude
    ubshr1, vbshr1 = mpcalc.bulk_shear(p, u, v, height=z, depth=1 * units.km)
    bshear1 = mpcalc.wind_speed(ubshr1, vbshr1)
    try:
        ubshr3, vbshr3 = mpcalc.bulk_shear(p, u, v, height=z, depth=3 * units.km)
        bshear3 = mpcalc.wind_speed(ubshr3, vbshr3)
    except:
        bshear3 = float("Nan")
    try:
        ubshr6, vbshr6 = mpcalc.bulk_shear(p, u, v, height=z, depth=6 * units.km)
        bshear6 = mpcalc.wind_speed(ubshr6, vbshr6)
    except:
        bshear6 = float("Nan")

    # Estimate height of LCL in meters from hydrostatic thickness (for sig_tor)
    new_p = np.append(p[p > lcl_pressure], lcl_pressure)
    new_t = np.append(T[p > lcl_pressure], lcl_temperature)
    lcl_height = mpcalc.thickness_hydrostatic(new_p, new_t)

    try:
        # Use all computed pieces to calculate the Significant Tornado parameter
        sig_tor = mpcalc.significant_tornado(sbcape, lcl_height,
                                             total_helicity3, bshear3).to_base_units()

        # Perform the calculation of supercell composite if an effective layer exists
        super_comp = mpcalc.supercell_composite(mucape, total_helicity3, bshear3)
    except:
        sig_tor = float("Nan")
        super_comp = float("Nan")
        pass


    def mag_int(param):
        try:
            param = int(param.m)
        except:
            try: 
                param = param.m
            except: param = param
        return param


    # PRINT VALUES OF PARAMETERS TO THE PLOT ------------------------------------------------
    # there is a lot we can do with this data operationally, so lets plot some of
    # these values right on the plot, in the box we made
    # first lets plo5 some thermodynamic parameters
    plt.figtext( 0.62, 0.37,  f'SBCAPE: ', weight='bold', fontsize=15, color='black', ha='left')
    plt.figtext( 0.75, 0.37,  f'{mag_int(sbcape)} J/kg', weight='bold', fontsize=15, color='orangered', ha='right')
    plt.figtext( 0.62, 0.34,  f'SBCIN: ', weight='bold', fontsize=15, color='black', ha='left')
    plt.figtext( 0.75, 0.34,  f'{mag_int(sbcin)} J/kg', weight='bold', fontsize=15, color='lightblue', ha='right')
    plt.figtext( 0.62, 0.29,  f'MLCAPE: ', weight='bold', fontsize=15, color='black', ha='left')
    plt.figtext( 0.75, 0.29,  f'{mag_int(mlcape)} J/kg', weight='bold', fontsize=15, color='orangered', ha='right')
    plt.figtext( 0.62, 0.26,  f'MLCIN: ', weight='bold', fontsize=15, color='black', ha='left')
    plt.figtext( 0.75, 0.26,  f'{mag_int(mlcin)} J/kg', weight='bold', fontsize=15, color='lightblue', ha='right')
    plt.figtext( 0.62, 0.21,  f'MUCAPE: ', weight='bold', fontsize=15, color='black', ha='left')
    plt.figtext( 0.75, 0.21,  f'{mag_int(mucape)} J/kg', weight='bold', fontsize=15, color='orangered', ha='right')
    plt.figtext( 0.62, 0.18,  f'MUCIN: ', weight='bold', fontsize=15, color='black', ha='left')
    plt.figtext( 0.75, 0.18,  f'{mag_int(mucin)} J/kg', weight='bold', fontsize=15, color='lightblue', ha='right')
    plt.figtext( 0.62, 0.13,  f'TT-INDEX: ', weight='bold', fontsize=15, color='black', ha='left')
    plt.figtext( 0.75, 0.13,  f'{mag_int(total_totals)} Δ°C', weight='bold', fontsize=15, color='orangered', ha='right')
    plt.figtext( 0.62, 0.10,  f'K-INDEX: ', weight='bold', fontsize=15, color='black', ha='left')
    plt.figtext( 0.75, 0.10,  f'{mag_int(kindex)} °C', weight='bold', fontsize=15, color='orangered', ha='right')
    # now some kinematic parameters
    met_per_sec = (units.m*units.m)/(units.sec*units.sec)
    plt.figtext( 0.77, 0.37,  f'0-1km SRH: ', weight='bold', fontsize=15, color='black', ha='left')
    plt.figtext( 0.92, 0.37,  f'{mag_int(total_helicity1)* met_per_sec:~P}', weight='bold', fontsize=15, color='navy', ha='right')
    plt.figtext( 0.77, 0.34,  f'0-1km SHEAR: ', weight='bold', fontsize=15, color='black', ha='left')
    plt.figtext( 0.92, 0.34,  f'{mag_int(bshear1)} kts', weight='bold', fontsize=15, color='blue', ha='right')
    plt.figtext( 0.77, 0.29,  f'0-3km SRH: ', weight='bold', fontsize=15, color='black', ha='left')
    plt.figtext( 0.92, 0.29,  f'{mag_int(total_helicity3)* met_per_sec:~P}', weight='bold', fontsize=15, color='navy', ha='right')
    plt.figtext( 0.77, 0.26,  f'0-3km SHEAR: ', weight='bold', fontsize=15, color='black', ha='left')
    plt.figtext( 0.92, 0.26,  f'{mag_int(bshear3)} kts', weight='bold', fontsize=15, color='blue', ha='right')
    plt.figtext( 0.77, 0.21,  f'0-6km SRH: ', weight='bold', fontsize=15, color='black', ha='left')
    plt.figtext( 0.92, 0.21,  f'{mag_int(total_helicity6)* met_per_sec:~P}', weight='bold', fontsize=15, color='navy', ha='right')
    plt.figtext( 0.77, 0.18,  f'0-6km SHEAR: ', weight='bold', fontsize=15, color='black', ha='left')
    plt.figtext( 0.92, 0.18,  f'{mag_int(bshear6)} kts', weight='bold', fontsize=15, color='blue', ha='right')
    plt.figtext( 0.77, 0.13,  f'SIG TORNADO: ', weight='bold', fontsize=15, color='black', ha='left')
    plt.figtext( 0.92, 0.13,  f'{mag_int(sig_tor)}', weight='bold', fontsize=15, color='orangered', ha='right')
    plt.figtext( 0.77, 0.10,  f'SUPERCELL COMP: ', weight='bold', fontsize=15, color='black', ha='left')
    plt.figtext( 0.92, 0.10,  f'{mag_int(super_comp)}', weight='bold', fontsize=15, color='orangered', ha='right')
    plt.figtext( 0.763, 0.065,  f'ECAPE: ', weight='bold', fontsize=15, color='black', ha='right')
    plt.figtext( 0.767, 0.065,  f'{mag_int(ecape)} J/kg', weight='bold', fontsize=15, color='orangered', ha='left')


    # LEGEND --------------------------------------------------------------------------------------
    skewleg = skew.ax.legend(loc='upper right', prop = { "size": 10 })
    hodoleg = h.ax.legend(loc='upper left')


    # PLOT TITLES ----------------------------------------------------------------------------------
    plt.figtext( 0.00, 0.96, left_title, ha='left', weight='bold', fontsize=20)
    plt.figtext( 0.96, 0.96, f'{right_title}     ', ha='right', weight='bold', fontsize=20)
    plt.figtext( 0.00, 0.986, top_title, ha='left', weight='bold', fontsize=27) 
    plt.figtext( 0.955, 0.03,  f'POWERED BY METPY & SOUNDERPY | (@WXKYLEGILLETT)     ', ha='right', color='blue', alpha=0.4, weight='bold', fontsize=12)

    img = Image.open(urlopen('https://user-images.githubusercontent.com/100786530/251580013-2e9477c9-e36a-4163-accb-fe46780058dd.png'))
    #                  side-side  up-down  size   size
    imgax = fig.add_axes([0.00, 0.06, 0.1, 0.1], anchor='SE')
    imgax.imshow(img)
    imgax.axis('off')
    
    print('> COMPLETE --------')
    elapsed_time = time.time() - st
    print('> RUNTIME:', time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
    
    plt.tight_layout()
    
    return plt




def metpy_sounding(clean_data, method='show', filename='sounderpy_sounding'):
    if method == 'show':
        __metpy_sounding(clean_data).show()
    elif method == 'save':
        __metpy_sounding(clean_data).savefig(filename, bbox_inches='tight')
    
#########################################################################################################    




metpy_sounding(clean_data, 'save')








### IMPORT SOFTWARE ###
#########################################################################################################
# BUILT IN
from datetime import datetime, timedelta
import time
import csv
import sys
import requests
from urllib.request import urlopen
from urllib.error import HTTPError   
import warnings
import bs4
# OTHER
import cdsapi
import pandas as pd
import xarray as xr
import numpy as np
from numpy import loadtxt
import numpy.ma as ma
from scipy import interpolate
from PIL import Image
from ecape.calc import calc_ecape, _get_parcel_profile, calc_mse, calc_integral_arg, calc_lfc_height, calc_el_height
# MATPLOTLIB
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
# METPY
import metpy.calc as mpcalc
from metpy.units import units
from metpy.plots import SkewT, Hodograph
# SIPHON 
from siphon.catalog import TDSCatalog
from siphon.ncss import NCSS
from siphon.simplewebservice.wyoming import WyomingUpperAir
from siphon.simplewebservice.iastate import IAStateUpperAir
from siphon.simplewebservice.igra2 import IGRAUpperAir
#########################################################################################################





### COPYRIGHT ###
#########################################################################################################
'''
        VERTICAL PROFILE DATA RETRIEVAL TOOL
        ------------------------------------
        This script it used to access vertical profile data for calculations or plotting of a vertical profile (sounding).
        
        RELEASE
        -------
        Version: 2.0.6 | Oct. 4, 2023
        
        DOCUMENTATION
        -------
        GiHub Wiki: https://github.com/kylejgillett/sounderpy/wiki

        COPYRIGHT
        ---------
        Created by Kyle J Gillett (@wxkylegillett) 2023

'''



citation_text = f"""
## ------------------ VERTICAL PROFILE DATA RETRIEVAL TOOL ------------------------ ##
##                    v2.0.6 | Oct 2023 | By Kyle J Gillett                         ##
##        RAOB, IGRA, RAP, RUC, NCEP, ERA5, RAP-ANALYSIS, BUFKIT & ACARS DATA       ##
## -------------------- THANK YOU FOR USING THIS PACKAGE! ------------------------- ##
"""
print(citation_text)

#########################################################################################################








### GET_MODEL_DATA() FUNCTION ###                                                               
#########################################################################################################

def get_model_data(method, latlon, year, month, day, hour, domain='point'):
    st = time.time()
    
    '''
    get_model_data(method, domain, latlons, year, month, day, hour)

    function used to access model reanalysis profile data from ERA5, RAP or RUC reanalysis archives
    or realtime RAP analysis data
    
    COPYRIGHT
    Created by Kyle J Gillett (@wxkylegillett) 2023

    '''
    
    # send error message if given method is invaild 
    if method.casefold() not in ['era', 'era5', 'rap', 'ruc', 'rap-ruc', 'rap-now', 'ncep-fnl', 'ncep']:
        #warnings.filterwarnings("ignore")
        sys.exit("Invalid 'method' kwarg. Valid methods inlcude ['era5', 'rap-ruc', 'rap-now', 'ncep-fnl']")
    
    # define global variables that will be needed in later 
    # functions (parse_data())
    global latlon1, source, year1, month1, day1, hour1
    latlon1 = latlon
    year1  = year
    month1 = month
    day1   = day
    hour1  = hour
    
    # create full lat-lon domain
    #               + lat         # - lat         # - lon           # + lon
    latlons = [latlon[0] + .5, latlon[0] - .5, latlon[1] - .5, latlon[1] + .5]

    
    
    
    
    ### ERA 5 REANALYSIS ###
    #########################################################################################################
    if method.casefold() in ['era', 'era5']:
        # define source 
        source = 'ERA5'
        
        print(f'> ERA5 REANALYSIS DATA ACCESS FUNCTION --\n------------------------------------------')
    
        # define ERA5 dataset names we want to acess data from
        dataset_presLvls = 'reanalysis-era5-pressure-levels'
        dataset_singleLvls = 'reanalysis-era5-single-levels'
        download_flag = 'false' 
        
        # define lat-lon box based on user defined domain (currently still in dev)
        if domain == 'point':
            latlon_list = [latlons[0], latlons[2], latlons[1], latlons[3]]
        elif domain == 'map':
            latlon_list = [latlons[0]+10, latlons[2]-10, latlons[1]-10, latlons[3]+10]
        else: 
            warnings.filterwarnings("ignore")
            sys.exit("'domain' kwarg must be 'point' or 'map'")

        # set up cds api call for pressure level data 
        c = cdsapi.Client()
        params = {
                'product_type':'reanalysis',
                'variable': ['temperature', 'geopotential', 'relative humidity', 'U WIND COMPONENT', 'V WIND COMPONENT', 'vertical_velocity'],
                'pressure_level': [
                    '100', '125',
                    '150', '175', '200',
                    '225', '250', '300',
                    '350', '400', '450',
                    '500', '550', '600',
                    '650', '700', '750',
                    '775', '800', '825',
                    '850', '875', '900',
                    '925', '950', '975',
                    '1000',
                ],
                'year'  : year,
                'month' : month,
                'day'   : day,
                'time'  : f'{hour}:00',
                'format': 'netcdf',
                'area'  : latlon_list
                }

        # set up cds api call for surface data 
        c2 = cdsapi.Client()
        params2 = {
                'product_type':'reanalysis',
                'variable': ['2m_temperature', '2m_dewpoint_temperature', 'surface_pressure', '10u', '10v', 'z', 'msl', 'cape'],
                'year'  : year,
                'month' : month,
                'day'   : day,
                'time'  : f'{hour}:00',
                'format': 'netcdf',
                'area'  : latlon_list
                }

        # retrieve data from CDS
        fl = c.retrieve(dataset_presLvls , params)
        print('> DATASET ACCESSED: '+dataset_presLvls )
        fl2 = c2.retrieve(dataset_singleLvls, params2)
        print('> DATASET ACCESSED: '+dataset_singleLvls )
        # load data to memory via output.nc files
        fl.download("./output.nc")
        fl2.download("./output.nc")
        
        # create xarray datasets from .nc files
        with urlopen(fl.location) as f:
            ds = xr.open_dataset(f.read())
        with urlopen(fl2.location) as f:
            ds2 = xr.open_dataset(f.read())
        
        # merge the two datasets together
        ds2 = ds2.rename({'z':'hgts','sp':'ps','t2m':'Ts','d2m':'tds','u10':'us','v10':'vs', })
        raw_data = xr.merge([ds,ds2])
        
        
        print('> COMPLETE --------')
        elapsed_time = time.time() - st
        print('> RUNTIME:', time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
        return raw_data
    #########################################################################################################

    
    
    
    
    ### RAP REANALYSIS ###
    #########################################################################################################
    if method in ['rap', 'ruc', 'rap-ruc']:
        
        print(f'> RAP REANALYSIS DATA ACCESS FUNCTION --\n-----------------------------------------')

        # define lat-lon box based on user defined domain (currently still in dev)
        if domain == 'point':
            latlon_list = [latlons[2],latlons[3],latlons[1],latlons[0]]
        elif domain == 'map':
            latlon_list = [latlons[2]-10 ,latlons[3]+10,latlons[1]-10,latlons[0]+10]
        else:
            warnings.filterwarnings("ignore")
            sys.exit("'domain' kwarg must be 'point' or 'map'")

        # create dict of RAP-data urls for the different versions of RAP and RUC from NCEI    
        urls = {
        'RAP_13km' : 'https://www.ncei.noaa.gov/thredds/ncss/model-rap130/'+str(year)+str(month)+'/'+str(year)+str(month)+str(day)+'/rap_130_'+str(year)+str(month)+str(day)+'_'+str(hour)+'00_000.grb2',
        'RAP_13km_old' : 'https://www.ncdc.noaa.gov/thredds/ncss/model-rap130-old/'+str(year)+str(month)+'/'+str(year)+str(month)+str(day)+'/rap_130_'+str(year)+str(month)+str(day)+'_'+str(hour)+'00_000.grb2',
            
        'RAP_13km_anl' : 'https://www.ncei.noaa.gov/thredds/ncss/model-rap130anl/'+str(year)+str(month)+'/'+str(year)+str(month)+str(day)+'/rap_130_'+str(year)+str(month)+str(day)+'_'+str(hour)+'00_000.grb2',
        'RAP_13km_anl_old' : 'https://www.ncdc.noaa.gov/thredds/ncss/model-rap130anl-old/'+str(year)+str(month)+'/'+str(year)+str(month)+str(day)+'/rap_130_'+str(year)+str(month)+str(day)+'_'+str(hour)+'00_000.grb2',
            
        'RAP_25km' : 'https://www.ncei.noaa.gov/thredds/ncss/model-rap252/'+str(year)+str(month)+'/'+str(year)+str(month)+str(day)+'/rap_252_'+str(year)+str(month)+str(day)+'_'+str(hour)+'00_000.grb2',
        'RAP_25km_old' : 'https://www.ncei.noaa.gov/thredds/ncss/model-rap252-old/'+str(year)+str(month)+'/'+str(year)+str(month)+str(day)+'/rap_252_'+str(year)+str(month)+str(day)+'_'+str(hour)+'00_000.grb2',
            
        'RAP_25km_anl' : 'https://www.ncei.noaa.gov/thredds/ncss/model-rap252anl/'+str(year)+str(month)+'/'+str(year)+str(month)+str(day)+'/rap_252_'+str(year)+str(month)+str(day)+'_'+str(hour)+'00_000.grb2',
        'RAP_25km_anl_old' : 'https://www.ncei.noaa.gov/thredds/ncss/model-rap252anl-old/'+str(year)+str(month)+'/'+str(year)+str(month)+str(day)+'/rap_252_'+str(year)+str(month)+str(day)+'_'+str(hour)+'00_000.grb2',
            
        'RUC_13km' : 'https://www.ncei.noaa.gov/thredds/ncss/model-ruc130anl/'+str(year)+str(month)+'/'+str(year)+str(month)+str(day)+'/ruc2anl_130_'+str(year)+str(month)+str(day)+'_'+str(hour)+'00_000.grb2',
        'RUC_13km_old' : 'https://www.ncei.noaa.gov/thredds/ncss/model-ruc130anl-old/'+str(year)+str(month)+'/'+str(year)+str(month)+str(day)+'/ruc2anl_130_'+str(year)+str(month)+str(day)+'_'+str(hour)+'00_000.grb2',
            
        'RUC_25km' : 'https://www.ncei.noaa.gov/thredds/ncss/model-ruc252anl/'+str(year)+str(month)+'/'+str(year)+str(month)+str(day)+'/ruc2anl_252_'+str(year)+str(month)+str(day)+'_'+str(hour)+'00_000.grb',
        'RUC_25km_old' : 'https://www.ncei.noaa.gov/thredds/ncss/model-ruc252anl/'+str(year)+str(month)+'/'+str(year)+str(month)+str(day)+'/ruc2anl_252_'+str(year)+str(month)+str(day)+'_'+str(hour)+'00_000.grb'
        }
        
        # create a simple test for each URL, use the first one that works 
        try:
            for url, key in zip(urls.values(), urls.keys()):
                try:
                    NCSS(url)
                    print(f'> DATASET USED: {key}')
                    url_to_use = url
                    source = str(key)[0:3]
                    break
                except:
                    pass
        except:
            pass
        try:
            data = NCSS(url_to_use)
        except:
            warnings.filterwarnings("ignore")
            sys.exit('NCSS URL FAILED -- THIS HAPPENS WHEN A BAD REQUEST IS MADE.\n> CHECK TO MAKE SURE YOU ENTERED THE CORRECT DATES.\n> NOTE: DATA IS NOT AVAILIABLE FOR EVERY DATE\n> THIS CATALOG OFTEN EXPERIENCES MISSING DATA/OUTAGES\n> MAKE SURE DATES ARE STRINGS -- MONTH, DAY AND HOUR MUST BE TWO DIGITS (EX: 18, 06, 00)')
            pass
        
        # set up TDS query 
        query = data.query()
        
        # subset data by variable names for RAP & RUC (of course they have to be different)
        if source == 'rap':
            query.variables('Pressure_surface',
                        'Geopotential_height_isobaric', 'Geopotential_height_surface',
                        'Temperature_isobaric', 'Temperature_height_above_ground',
                        'Relative_humidity_isobaric', 'Dewpoint_temperature_height_above_ground',
                        'Relative_humidity_height_above_ground',
                        'u-component_of_wind_height_above_ground', 'v-component_of_wind_height_above_ground', 
                        'u-component_of_wind_isobaric', 'v-component_of_wind_isobaric').add_lonlat()
        else:
            query.variables('Pressure_surface',
                        'Geopotential_height_isobaric', 'Geopotential_height_surface',
                        'Temperature_isobaric', 'Temperature_height_above_ground',
                        'Relative_humidity_isobaric','Dewpoint_temperature_height_above_ground',
                        'Relative_humidity_height_above_ground',
                        'u-component_of_wind_height_above_ground', 'v-component_of_wind_height_above_ground', 
                        'u-component_of_wind_isobaric', 'v-component_of_wind_isobaric').add_lonlat()
            
        # subset data by requested domain
        query.lonlat_box(latlon_list[0], latlon_list[1], latlon_list[2], latlon_list[3])
        
        # laod the data from TDS
        raw_data = data.get_data(query)

        print('> COMPLETE --------')
        
        elapsed_time = time.time() - st
        print('> RUNTIME:', time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
        return raw_data
    #########################################################################################################
    
    
    
    
    ### NCEP FNL REANALYSIS ###
    #########################################################################################################
    if method.casefold() in ['ncep-fnl', 'ncep']:
        # define source 
        source = 'NCEP-FNL'
        
        print(f'> NCEP-FNL REANALYSIS DATA ACCESS FUNCTION --\n------------------------------------------')
        
        # define lat-lon box based on user defined domain (currently still in dev)
        if domain == 'point':
            latlon_list = [latlons[0], latlons[2], latlons[1], latlons[3]]
        elif domain == 'map':
            latlon_list = [latlons[0]+10, latlons[2]-10, latlons[1]-10, latlons[3]+10]
        else: 
            warnings.filterwarnings("ignore")
            sys.exit("'domain' kwarg must be 'point' or 'map'")
            
        
        # set up NCEP FNL url
        url = f"https://thredds.rda.ucar.edu/thredds/ncss/grid/files/g/ds083.3/{year}/{year}{month}/gdas1.fnl0p25.{year}{month}{day}{hour}.f00.grib2"

        # access ncss thredds server 
        data = NCSS(url)

        # set up TDS query 
        query = data.query()

        query.variables(
                    'Geopotential_height_isobaric', 'Geopotential_height_surface',
                    'Temperature_isobaric', 'Temperature_height_above_ground',
                    'Relative_humidity_isobaric', 'Dewpoint_temperature_height_above_ground',
                    'Relative_humidity_height_above_ground', 'Pressure_surface',
                    'u-component_of_wind_height_above_ground', 'v-component_of_wind_height_above_ground', 
                    'u-component_of_wind_isobaric', 'v-component_of_wind_isobaric').add_lonlat()

        # subset data by requested domain

        # north=90.000&    west=-.125&    east=-.125&    south=-90.000
        query.lonlat_box(latlon_list[1], latlon_list[3], latlon_list[2], latlon_list[0])

        # laod the data from TDS
        raw_data = data.get_data(query)

        print('> COMPLETE --------')

        elapsed_time = time.time() - st
        print('> RUNTIME:', time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
        return raw_data
    #########################################################################################################
    
    
    
    
    
    ### RAP ANALYSIS ###
    #########################################################################################################
    
    if method in ['rap-now']:

        print(f'> RAP REANALYSIS DATA ACCESS FUNCTION --\n-----------------------------------------')
        
        # define lat-lon box based on user defined domain (currently still in dev)
        if domain == 'point':
            latlon_list = [latlons[2],latlons[3],latlons[1],latlons[0]]
        elif domain == 'map':
            latlon_list = [latlons[2]-10 ,latlons[3]+10,latlons[1]-10,latlons[0]+10]
        else:
            warnings.filterwarnings("ignore")
            sys.exit("'domain' kwarg must be 'point' or 'map'")
        
        # define dataset URL & try to access it to make sure it works 
        url = 'http://thredds.ucar.edu/thredds/catalog/grib/NCEP/RAP/CONUS_13km/latest.xml'
        try:
            cat = TDSCatalog('http://thredds.ucar.edu/thredds/catalog/grib/NCEP/RAP/CONUS_13km/latest.xml')
            source = 'RAP ANALYSIS'
        except:
            sys.exit("NCSS URL FAILED -- THIS HAPPENS WHEN A BAD REQUEST IS MADE. RAP Analysis data may not be available at this time.")
            pass
        
        # set up TDS query    
        latest_ds = list(cat.datasets.values())[0]
        ncss = NCSS(latest_ds.access_urls['NetcdfSubset'])
        query = ncss.query()
        # Find start time
        start_time = ncss.metadata.time_span['begin']
        fcst_date = datetime.strptime(start_time, '%Y-%m-%dT%H:%M:%SZ')
        year1  = fcst_date.strftime('%Y')
        month1 = fcst_date.strftime('%m')
        day1   = fcst_date.strftime('%d')
        hour1  = fcst_date.strftime('%H')
        # Subset data by time
        query.time(fcst_date).accept('netcdf4')
        # Subsets data by variables 
        query.variables('MSLP_MAPS_System_Reduction_msl',
                'Pressure_surface',
                'Geopotential_height_isobaric',
                'Temperature_isobaric',
                'Relative_humidity_isobaric',
                'Temperature_height_above_ground',
                'Relative_humidity_height_above_ground',
                'u-component_of_wind_height_above_ground', 'v-component_of_wind_height_above_ground', 
                'u-component_of_wind_isobaric', 'v-component_of_wind_isobaric',
                'Convective_available_potential_energy_surface',
                'Storm_relative_helicity_height_above_ground_layer',
                'U-Component_Storm_Motion_height_above_ground_layer',
                'V-Component_Storm_Motion_height_above_ground_layer').add_lonlat()

        # Subset data by lat-lon domain 
        query.lonlat_box(latlon_list[0], latlon_list[1], latlon_list[2], latlon_list[3])
        
        # Gets data
        raw_data = ncss.get_data(query)

        print('> COMPLETE --------')
        elapsed_time = time.time() - st
        print('> RUNTIME:', time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
        
        return raw_data
    #########################################################################################################
    
#########################################################################################################    
    
    
    
    
    
    
    
### GET_OBS_DATA() FUNCTION ###
#########################################################################################################


def get_obs_data(station, year, month, day, hour):
    
    # record process time
    st = time.time()
    
    
    """
    get_obs_data(station, year, month, day, hour)

    function used to access and parse RAOB profile data
    this function will search UW, ISU & IRGA2 for desired 
    station and date
        
    COPYRIGHT
    Created by Kyle J Gillett (@wxkylegillett) 2023
    
    """
    
    
    print(f'> OBSERVED DATA ACCESS FUNCTION --\n-----------------------------------')
    
    
    
    # set up global variables 
    global station1, source, year1, month1, day1, hour1
    source = 'obs'
    station1 = station
    year1  = year
    month1 = month
    day1   = day
    hour1  = hour
    
    # get station lists from SounderPy GitHub Repo
    RAOB_STATIONS = pd.read_csv(f'https://raw.githubusercontent.com/kylejgillett/sounderpy/main/src/RAOB-STATIONS.txt', 
                                skiprows=7, skipinitialspace = True)
    IGRA_STATIONS = pd.read_csv(f'https://raw.githubusercontent.com/kylejgillett/sounderpy/main/src/IGRA-STATIONS.txt', 
                                skiprows=7, skipinitialspace = True)
    
    # create dt object from user date for siphon data access 
    dt = datetime(int(year), int(month), int(day), int(hour))
    
    got_data = False
    
    # set up siphon API call for raob data -- if station ID is found in RAOB_STATIONS, it is
    # a RAOB ID and siphon UW or ISU must be used to get data
    if len(station) == 11:
        search_for = 'igra'
    else:
        search_for = 'raob'
        
    
    ### RAOB OBSERVATIONS ###
    #########################################################################################################
    if search_for == 'raob':
        # try this process 10 times, sometimes requests fail due to temporary 404 errors
        for i in range(1, 11):
            try: 
                # try UW data request
                df = WyomingUpperAir.request_data(dt, station)
                got_data = True
                if got_data == True:
                    print(f'> FOUND RAOB: {station} on {month}/{day}/{year} at {hour}z | From UW')
                    break
            except:
                if i == 10:
                    print('! NO DATA FOR THIS STATION & DATE IN THE UW ARCHIVE, TRYING IEM...')
                    got_data = False          
            pass
        
        # if UW fails, try ISU data request 
        if got_data == False: 
            for i in range(1, 11):
                try: 
                    df = IAStateUpperAir.request_data(dt, station)
                    got_data = True
                    if got_data == True:
                        print(f'> FOUND RAOB: {station} on {month}/{day}/{year} at {hour}z | From IEM')
                        break
                except:
                    if i == 10:
                        warnings.filterwarnings("ignore")
                        sys.exit(f"There is likely no available data for station {station} on {month}/{day}/{year} at {hour}z\nPlease make sure you entered a valid station ID and date\nor try a different time, date, or station")
                        got_data = False          
                pass
        
        
        # if data is found, parse data and create a dict of clean data
        if got_data == True:
            station = RAOB_STATIONS['ICAO'][np.where(RAOB_STATIONS['ICAO'].str.contains(station, na=False, case=True))[0]].values[0].strip()

            # create dict of data
            new_keys = ['p', 'z', 'T', 'Td', 'u', 'v']
            old_keys = ['pressure', 'height', 'temperature', 'dewpoint', 'u_wind', 'v_wind'] # 'latitude', 'longitude']
            units_list = ['hPa', 'meter', 'degC', 'degC', 'kt', 'kt']
            clean_data = {}
            non_dups = np.concatenate(([True], np.diff(df.to_dict('list')['pressure']) != 0))
            for old_key, new_key, unit in zip (old_keys, new_keys, units_list):
                clean_data[new_key] = np.array(df.to_dict('list')[old_key])[non_dups]*units(unit)
            clean_data['site_info'] = {
                'site-id'   : RAOB_STATIONS[RAOB_STATIONS['ICAO']==station]['ICAO'].values[0],
                'site-name' : RAOB_STATIONS[RAOB_STATIONS['ICAO']==station]['NAME'].values[0],
                'site-lctn' : RAOB_STATIONS[RAOB_STATIONS['ICAO']==station]['LOC'].values[0],
                'site-latlon' : get_latlon('raob', station),
                'site-elv'  : RAOB_STATIONS[RAOB_STATIONS['ICAO']==station]['EL(m)'].values[0],
                'source'    : 'RAOB OBSERVED PROFILE',
                'model'     : 'no-model',
                'fcst-hour' : 'no-fcst-hour',
                'run-time'  : ['no-run-time'],
                'valid-time': [year, month, day, hour]}
            
            try:
                # trim data to 98hPa and below for less process time 
                slc = (len(clean_data['p']) - np.where(clean_data['p']<=98.*units('hPa'))[0][0])
                for key in new_keys:
                    clean_data[key] = clean_data[key][:-slc]
            except:
                pass
            
            
            elapsed_time = time.time() - st
            print('> RUNTIME:', time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
            return clean_data 
    #########################################################################################################
    
    
    

    ### IGRAv2 OBSERVATIONS ###    
    #########################################################################################################
    # if user station ID is not a RAOB ID, it may be a IGRA2 ID, so try searching for it in IGRA_STATIONS csv        
    elif search_for =='igra': #IGRA_STATIONS['ID'].str.contains(station, na=False, case=True).sum() >= 1:
        for i in range(1, 3):
            try: 
                # try siphon IGRA request 
                df = IGRAUpperAir.request_data(dt, station)
                got_data = True
                if got_data == True:
                    print(f'> FOUND DATA: {station} on {month}/{day}/{year} at {hour}z | From IGRAv2')
                    break
            except:
                if i == 3:
                    warnings.filterwarnings("ignore")
                    sys.exit(f"There is likely no available data for station {station} on {month}/{day}/{year} at {hour}z\nPlease make sure you entered a valid station ID and date\nor try a different time, date, or station")
                    got_data = False        
        
        
        # if data is found, parse data and create a dict of clean data
        if got_data == True:
            station = IGRA_STATIONS['ID'][np.where(IGRA_STATIONS['ID'].str.contains(station, na=False, case=True))[0]].values[0].strip()
            
            # create dict of data
            head = df[1]
            df = df[0]
            new_keys = ['p', 'z', 'T', 'Td', 'u', 'v']
            old_keys = ['pressure', 'height', 'temperature', 'dewpoint', 'u_wind', 'v_wind'] # 'latitude', 'longitude']
            units_list = ['hPa', 'meter', 'degC', 'degC', 'kt', 'kt']
            clean_data = {}
            zflag=np.array(df['zflag'])
            pflag=np.array(df['pflag'])
            tflag=np.array(df['tflag'])
            for old_key, new_key, unit in zip (old_keys, new_keys, units_list):
                clean_data[new_key] = np.array(df.to_dict('list')[old_key])[zflag+pflag+tflag>=4]*units(unit)
            clean_data['site_info'] = {
                    'site-id'   : IGRA_STATIONS[IGRA_STATIONS['ID']==station]['ID'].str.strip().values[0],
                    'site-name' : IGRA_STATIONS[IGRA_STATIONS['ID']==station]['NAME'].str.strip().values[0],
                    'site-lctn' : 'no-site-location',
                    'site-latlon' : get_latlon('igra', station),
                    'site-elv'  : IGRA_STATIONS[IGRA_STATIONS['ID']==station]['EL(m)'].values[0],
                    'source'    : 'IGRA OBSERVED PROFILE',
                    'model'     : 'no-model',
                    'fcst-hour' : 'no-fcst-hour',
                    'run-time'  : ['no-run-time'],
                    'valid-time': [year, month, day, hour]}
            # correct u & v units
            clean_data['u'] = clean_data['u']*1.94384
            clean_data['v'] = clean_data['v']*1.94384
           
        
            print('> COMPLETE --------')
            elapsed_time = time.time() - st
            print('> RUNTIME:', time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
            
            return clean_data     
    #########################################################################################################
        
        

        else: 
            warnings.filterwarnings("ignore")
            sys.exit(f"There is likely no available data for station {station} on {month}/{day}/{year} at {hour}z\nPlease make sure you entered a valid station ID and date\nor try a different time, date, or station")
        
    else:
        warnings.filterwarnings("ignore")
        sys.exit(f"There is likely no available data for station {station} on {month}/{day}/{year} at {hour}z\nPlease make sure you entered a valid station ID and date\nor try a different time, date, or station")
        
        
    
#########################################################################################################

    

    
    
    
    
    
    

### GET_BUFKIT_DATA() FUNCTION ###                                                             
#########################################################################################################


def get_bufkit_data(model, station, fcst_hour, run_year=None, run_month=None, run_day=None, run_hour=None):
    
    # record process time
    st = time.time()
    
    
    '''
    function used to load and parse BUFKIT sounding data from PSU & ISU

    COPYRIGHT
    Created by Kyle J Gillett (@wxkylegillett) 2023

    '''
    
    
    print(f'> BUFKIT DATA ACCESS FUNCTION --\n---------------------------------')
    
    # make sure model is lower-case
    model = str.lower(model)
    
    # GET MOST-RECENT RUNS FROM PSU SERVERS 
    # if date variables (year, month, day) are not given, the user has 'selected' a most
    # recent forecast run, get that from PSU
    if run_year == None:
        if model not in ['gfs', 'nam', 'namnest', 'rap', 'hrrr', 'sref', 'hiresw']:
            sys.exit("NOT A VALID MODEL -- VALID MODELS ARE ['GFS', 'NAM', 'NAMNEST', 'RAP', 'HRRR', 'SREF', 'HIRESW']")
        if model == 'gfs':
            model3 = 'gfs3' 
        else:
            model3 = model
        data_conn = f'http://www.meteo.psu.edu/bufkit/data/{model.upper()}/{model3}_{station.lower()}.buf'
     
    
    # GET ARCHIVE DATA FROM THE IEM SERVERS. CORRECT GFS & NAM MODEL NAMES
    # if date variables (year, month, day) are given, the user has 'selected' a 
    # archived forecast for the given date 
    else:
        if model not in ['gfs', 'nam', 'namnest', 'rap', 'hrrr']:
            sys.exit("NOT A VALID MODEL -- VALID MODELS ARE ['GFS', 'NAM', 'NAMNEST', 'RAP', 'HRRR']")
        if model == 'namnest':
            model == 'nam4km'
        if model == 'gfs':
            model3 = 'gfs3' 
        else:
            model3 = model
        data_conn = f'https://mtarchive.geol.iastate.edu/{run_year}/{run_month}/{run_day}/bufkit/{run_hour}/{model}/{model3}_{station.lower()}.buf'  

        
    # GET BUFKIT FILE 
    # CONVERT LINES OF BYTES TO STRINGS 
    buf_file = urlopen(data_conn)
    buf_file = [str(line).replace("b'", "").replace("\\r\\n'", "") for line in buf_file]
    
    # GET BUFKIT STATIONS LISTING FROM SOUNDERPY GITHUB REPO
    BUFKIT_STATIONS = pd.read_csv(f'https://raw.githubusercontent.com/kylejgillett/sounderpy/main/src/BUFKIT-STATIONS-MASTER.txt', 
                                  skiprows=7, skipinitialspace = True)
    
    # FIND THE STATION DATA 
    station = BUFKIT_STATIONS['ID'][np.where(BUFKIT_STATIONS['ID'].str.contains(station, na=False, case=True))[0]].values[0]
    
    # CREATE TEMP DATA
    tmp_data, sounding_headers, derived_headers = [], '', ''
    recordSounding = False
    
    # SET UP DATE / TIME OBJECTS FROM THE BUFKIT FILE
    run_time = buf_file[4][buf_file[4].index('TIME') + 7:(buf_file[4].index('TIME')+9)+9]
    run_dt = datetime(int(f'20{run_time[0:2]}'), int(run_time[2:4]), int(run_time[4:6]), int(run_time[7:9])) 
    fct_dt = run_dt + timedelta(hours = fcst_hour)
    stp_dt = fct_dt + timedelta(hours = 1)
    
    # Loop over each line in data file
    for line in buf_file:
        # Find start of sounding data
        if f'TIME = {fct_dt.strftime("%Y")[2:4]}{fct_dt.strftime("%m")}{fct_dt.strftime("%d")}/{fct_dt.strftime("%H")}00' in line:
            recordSounding=True  
        if 'SNPARM' in line:
            sounding_headers=line[line.index('=')+2:].replace(' ', '').split(';')
        if 'STNPRM' in line:
            derived_headers=line[line.index('=')+2:].replace(' ', '').split(';')
        # Append data line to temp data list
        if recordSounding:
            tmp_data.append(line)
        # Break out of loop when end key reached
        if f'TIME = {stp_dt.strftime("%Y")[2:4]}{stp_dt.strftime("%m")}{stp_dt.strftime("%d")}/{stp_dt.strftime("%H")}00' in line:
            tmp_data.pop(-1)
            break
    
    # SET UP UTILS
    station_headers=['STID', 'STNM', 'TIME', 'SLAT', 'SLON', 'SELV', 'STIM']
    tmp_str=''
    recordStationInfo, recordDerivedQty, recordSoundingQty = False, False, True
    station_metadata, derived_data, sounding_data = [], [], []

    # PARSE THROUGH FILE, SPLIT LINES AND RECORD DATA WE WANT TO KEEP
    for line in tmp_data:
        # Check for station infromation
        if recordStationInfo and line=='':
            # Break values up to only be seperated by one whitespace
            station_info=(tmp_str.replace(' = ', ' '))
            # Split values into list
            station_info=station_info.split(' ')
            # Remove label values
            station_info=[x for x in station_info if x not in station_headers]
            while '' in station_info:
                station_info.remove('')
            # Add to main list
            station_metadata.append(station_info)
            # Reset temp vars
            tmp_str=''
            recordStationInfo=False
        if any(var in line for var in station_headers):
            recordStationInfo=True
            tmp_str+=(' ' + line)
        # Check for derived sounding quantities
        if recordDerivedQty==True and line=='':
            # Break values up to only be seperated by one whitespace
            derived_qty=(tmp_str.replace(' = ', ' '))
            # Split values into list
            derived_qty=derived_qty.split(' ')
            # Remove non-numeric values
            derived_qty=[x for x in derived_qty if x not in derived_headers]
            while '' in derived_qty:
                derived_qty.remove('')
            # Add to main list
            derived_data.append(derived_qty)
            # Reset temp vars
            tmp_str=''
            recordDerivedQty=False
        if any(var in line for var in derived_headers):
            recordDerivedQty=True
            tmp_str+=(' ' + line)
        # Check for sounding quantities
        if any(var in line for var in sounding_headers):
            recordSoundingQty=True
        if recordSoundingQty and line=='':
            level_list=[]
            # Split data string into values
            data_list=tmp_str.split(' ')
            # Remove empty indices
            while '' in data_list:
                data_list.remove('')
            # Break data up into pressure levels
            for i in range(0, len(data_list), len(sounding_headers)):
                level_list.append(data_list[i:len(sounding_headers)+i])
        elif recordSoundingQty:
            if any(var in line for var in sounding_headers)==False:
                tmp_str+=(' ' + line)

                
    # CREATE BLANK LISTS 
    p = []
    z = []
    T = []
    Td = []
    ws = []
    wd = []

    
    # APPEND LISTS WITH DATA FROM BUFKIT FILES
    if model in ['gfs']:
        for i in range(0, len(level_list)):
            p.append(float(level_list[i][0]))
            z.append(float(level_list[i][8]))
            T.append(float(level_list[i][1]))
            Td.append(float(level_list[i][3]))
            ws.append(float(level_list[i][6]))
            wd.append(float(level_list[i][5]))  
    else:
        for i in range(0, len(level_list)):
            p.append(float(level_list[i][0]))
            z.append(float(level_list[i][9]))
            T.append(float(level_list[i][1]))
            Td.append(float(level_list[i][3]))
            ws.append(float(level_list[i][6]))
            wd.append(float(level_list[i][5])) 
     
    
    # CALCULATE U AND V COMPONENTS 
    u = list(mpcalc.wind_components(ws*units.kts, wd*units.degrees)[0].m)
    v = list(mpcalc.wind_components(ws*units.kts, wd*units.degrees)[1].m)
    
    
    # DEFINE find_nearest() FUNCTION 
    def find_nearest(array, value):
        array = np.asarray(array)
        nearest_idx = (np.abs(array - value)).argmin()
        return nearest_idx
    
    
    # FIND P LEVEL AT 50hPa
    hPa50 = find_nearest(p, 50)

    
    # ARRANGE DATA IN CLEAN_DATA DICT
    clean_data = {}
    lists = [p[0:hPa50], z[0:hPa50], T[0:hPa50], Td[0:hPa50], u[0:hPa50], v[0:hPa50]]
    keys = ['p', 'z', 'T', 'Td', 'u', 'v']
    units_list = ['hPa', 'meter', 'degC', 'degC', 'kt', 'kt']
    for key, lst, unit in zip(keys, lists, units_list):
        clean_data[key] = lst*units(unit) 
    
    # create dict of data
    clean_data['site_info'] = {
                'site-id'   : BUFKIT_STATIONS[BUFKIT_STATIONS['ID']==station]['ID'].str.strip().values[0],
                'site-name' : BUFKIT_STATIONS[BUFKIT_STATIONS['ID']==station]['NAME'].str.strip().values[0],
                'site-lctn' : BUFKIT_STATIONS[BUFKIT_STATIONS['ID']==station]['LOC'].str.strip().values[0],
                'site-latlon' : [BUFKIT_STATIONS[BUFKIT_STATIONS['ID']==station]['LAT'].values[0],
                                 BUFKIT_STATIONS[BUFKIT_STATIONS['ID']==station]['LON'].values[0]],
                'site-elv'  : BUFKIT_STATIONS[BUFKIT_STATIONS['ID']==station]['EL(m)'].values[0],
                'source'    : 'BUFKIT FORECAST PROFILE',
                'model'     : str.upper(model),
                'fcst-hour' : f'F0{fcst_hour}',
                'run-time'  : [run_dt.strftime("%Y"), run_dt.strftime("%m"), run_dt.strftime("%d"), run_dt.strftime("%H")],
                'valid-time': [fct_dt.strftime("%Y"), fct_dt.strftime("%m"), fct_dt.strftime("%d"), fct_dt.strftime("%H")]} 
    
    
    print('> COMPLETE --------')
    elapsed_time = time.time() - st
    print('> RUNTIME:', time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
    return clean_data
#########################################################################################################  
    
    
    

    
    
    
    
    
### PARSE_DATA() FUNCTION ###
#########################################################################################################
    
def parse_data(raw_data):
    
    # record process time
    st = time.time()
    
    
    '''
    function used to parse sounding data after using a get_model_data() function

    COPYRIGHT
    Created by Kyle J Gillett (@wxkylegillett) 2023

    '''
    
    # if dataset is a xarray.core dataset, it came from the ERA5
    # and specific processing of this data is needed
    if str(type(raw_data)) == "<class 'xarray.core.dataset.Dataset'>":
        
        print(f'> ERA5 REANALYSIS DATA PARSE FUNCTION --\n------------------------------------------')
        vert_data = {
        'vert_T' : (np.array(raw_data['t'][0,:,:,:]-273.15)),
        'vert_p' : (np.array(raw_data['level'])),
        'vert_z' : (np.array(raw_data['z'][0])/9.80665),
        'vert_rh': (np.array(raw_data['r'][0])),
        'vert_u' : (np.array(raw_data['u'][0])*1.94384),
        'vert_v' : (np.array(raw_data['v'][0])*1.94384),
        'vert_Td': (mpcalc.dewpoint_from_relative_humidity(
            (np.array(raw_data['t'][0,:,:,:]-273.15))*units.degC, 
            (np.array(raw_data['r'][0]))*units.percent)).m,
        }  
        sfc_data = {
        'sfc_T' : (np.array(raw_data['Ts'][0,:])-273.15),
        'sfc_p' : (np.array(raw_data['ps'][0,:])/100),
        'sfc_z' : (np.array(raw_data['hgts'][0,:])/9.80665),
        'sfc_Td': (np.array(raw_data['tds'][0,:])-273.15),
        'sfc_rh': (mpcalc.relative_humidity_from_dewpoint(
            (np.array(raw_data['Ts'][0,:])-273.15)*units.degC, 
            (np.array(raw_data['tds'][0,:])-273.15)*units.degC)*100),
        'sfc_u' : (np.array(raw_data['us'][0,:])*1.94384),
        'sfc_v' : (np.array(raw_data['vs'][0,:])*1.94384)
        }
        latlon_data = {
        'data_lat'    : (raw_data['latitude'][:]),
        'data_lon'    : (raw_data['latitude'][:]),
        'data_latnum' : (raw_data['latitude'][:]).shape[0],
        'data_lonnum' : (raw_data['latitude'][:]).shape[0]
        }


    
    
    # if data is a netCDF4 dataset, it is RAP, RUC or NCEP data
    if str(type(raw_data)) == "<class 'netCDF4._netCDF4.Dataset'>":
            
            print(f'> REANALYSIS DATA PARSE FUNCTION --\n------------------------------------------')
            
            # if Geopotential_height_surface exists within the dataset, it is reanalysis data
            if "Geopotential_height_surface" in raw_data.variables.keys():
                
                # try to determine how many isobaric levels exist in the dataset
                # this is mainly for the NCEP-FNL dataset which annoyingly 
                # comes with varying isobaric level intervals 
                
                try:
                    pressures = (ma.getdata(raw_data.variables['isobaric'])).data/100
                except:
                    try:
                        pressures = (ma.getdata(raw_data.variables['isobaric1'])).data/100
                    except:
                        try:
                            pressures = (ma.getdata(raw_data.variables['isobaric2'])).data/100
                        except:
                            try:
                                pressures = (ma.getdata(raw_data.variables['isobaric3'])).data/100
                            except:
                                pass
                            pass
                        pass
                    pass
                       
                # pressures = (np.array([100, 125, 150, 175, 200, 225, 250, 275, 
                #                                300, 325, 350, 375, 400, 425, 450, 475, 
                #                                500, 525, 550, 575, 600, 625, 650, 675, 
                #                                700, 725, 750, 775, 800, 825, 850, 875, 900, 925, 950, 975, 1000]))

                # create a dict of vertical data 
                vert_data = {
                    'vert_T' : (ma.getdata(raw_data.variables['Temperature_isobaric'][0,:,:,:]-273.15)),
                    'vert_p' : pressures,
                    'vert_z' : (ma.getdata(raw_data.variables['Geopotential_height_isobaric'][0,:,:,:])),
                    'vert_rh': (ma.getdata(raw_data.variables['Relative_humidity_isobaric'][0,:,:,:])),
                    'vert_u' : (ma.getdata(raw_data.variables['u-component_of_wind_isobaric'][0,:,:,:]*1.94384)),
                    'vert_v' : (ma.getdata(raw_data.variables['v-component_of_wind_isobaric'][0,:,:,:]*1.94384)),
                    'vert_Td': (mpcalc.dewpoint_from_relative_humidity(
                        (ma.getdata(raw_data.variables['Temperature_isobaric'][0,:,:,:]-273.15))*units.degC, 
                        (ma.getdata(raw_data.variables['Relative_humidity_isobaric'][0,:,:,:])*units.percent))).m,
                } 

                # create a dict of sfc data 
                sfc_data = {
                    'sfc_T' : (ma.getdata(raw_data.variables['Temperature_height_above_ground'][0,0,:,:]-273.15)),
                    'sfc_p' : (ma.getdata(raw_data.variables['Pressure_surface'][0,:,:]/100)),
                    'sfc_z' : (ma.getdata(raw_data.variables['Geopotential_height_surface'][0,:,:])),
                    'sfc_rh': (ma.getdata(raw_data.variables['Relative_humidity_height_above_ground'][0,0,:,:])),
                    'sfc_u' : (ma.getdata(raw_data.variables['u-component_of_wind_height_above_ground'][0,0,:,:]*1.94384)),
                    'sfc_v' : (ma.getdata(raw_data.variables['v-component_of_wind_height_above_ground'][0,0,:,:]*1.94384)),
                    'sfc_Td': (mpcalc.dewpoint_from_relative_humidity(
                        (ma.getdata(raw_data.variables['Temperature_height_above_ground'][0,0,:,:]-273.15))*units.degC, 
                        (ma.getdata(raw_data.variables['Relative_humidity_height_above_ground'][0,0,:,:])*units.percent))).m,
                } 

                # create a dict of lat/lon information 
                latlon_data = {
                'data_lat'    : (ma.getdata(1000*(raw_data.variables[raw_data.variables['Pressure_surface'].dimensions[1]][:]))),
                'data_lon'    : (ma.getdata(1000*(raw_data.variables[raw_data.variables['Pressure_surface'].dimensions[2]][:]))),
                'data_latnum' : (ma.getdata(1000*(raw_data.variables[raw_data.variables['Pressure_surface'].dimensions[1]][:]))).shape[0],
                'data_lonnum' : (ma.getdata(1000*(raw_data.variables[raw_data.variables['Pressure_surface'].dimensions[2]][:]))).shape[0]
                }

            else: 
                
                # if Geopotential_height_surface is not in the dataset, then it is RAP analysis data 
                # declare pressures 
                pressures = (np.array([100, 125, 150, 175, 200, 225, 250, 275, 
                               300, 325, 350, 375, 400, 425, 450, 475, 
                               500, 525, 550, 575, 600, 625, 650, 675, 
                               700, 725, 750, 775, 800, 825, 850, 875, 900, 925, 950, 975, 1000]))
                
                # create dict of vertical data 
                vert_data = {
                'vert_T' : (ma.getdata(raw_data.variables['Temperature_isobaric'][0,:,:,:]-273.15)),
                'vert_p' : pressures,
                'vert_z' : (ma.getdata(raw_data.variables['Geopotential_height_isobaric'][0,:,:,:])),
                'vert_rh': (ma.getdata(raw_data.variables['Relative_humidity_isobaric'][0,:,:,:])),
                'vert_u' : (ma.getdata(raw_data.variables['u-component_of_wind_isobaric'][0,:,:,:]*1.94384)),
                'vert_v' : (ma.getdata(raw_data.variables['v-component_of_wind_isobaric'][0,:,:,:]*1.94384)),
                'vert_Td': (mpcalc.dewpoint_from_relative_humidity(
                    (ma.getdata(raw_data.variables['Temperature_isobaric'][0,:,:,:]-273.15))*units.degC, 
                    (ma.getdata(raw_data.variables['Relative_humidity_isobaric'][0,:,:,:])*units.percent))).m,
            } 
                
                # create dict of surface data 
                sfc_data = {
                'sfc_T' : (ma.getdata(raw_data.variables['Temperature_height_above_ground'][0,0,:,:]-273.15)),
                'sfc_p' : (ma.getdata(raw_data.variables['Pressure_surface'][0,:,:]/100)),
                #'sfc_z' : (ma.getdata(raw_data.variables['Geopotential_height_surface'][0,:,:])),
                'sfc_rh': (ma.getdata(raw_data.variables['Relative_humidity_height_above_ground'][0,0,:,:])),
                'sfc_u' : (ma.getdata(raw_data.variables['u-component_of_wind_height_above_ground'][0,0,:,:]*1.94384)),
                'sfc_v' : (ma.getdata(raw_data.variables['v-component_of_wind_height_above_ground'][0,0,:,:]*1.94384)),
                'sfc_Td': (mpcalc.dewpoint_from_relative_humidity(
                    (ma.getdata(raw_data.variables['Temperature_height_above_ground'][0,0,:,:]-273.15))*units.degC, 
                    (ma.getdata(raw_data.variables['Relative_humidity_height_above_ground'][0,0,:,:])*units.percent))).m,
                } 
                
                # calculate surface heights 
                sfc_data['sfc_z'] = np.zeros(sfc_data['sfc_p'].shape)
                for i in range(0,sfc_data['sfc_p'].shape[0]):
                    for j in range(0,sfc_data['sfc_p'].shape[1]):
                        sfc_data['sfc_z'][i,j]=np.interp([sfc_data['sfc_p'][i,j]],vert_data['vert_p'],vert_data['vert_z'][:,i,j])[0]
                        vert_data['vert_z'][:,i,j]=vert_data['vert_z'][:,i,j]-sfc_data['sfc_z'][i,j]
                
                # create dict of latlon information 
                latlon_data = {
                'data_lat'    : (ma.getdata(1000*(raw_data.variables[raw_data.variables['Pressure_surface'].dimensions[1]][:]))),
                'data_lon'    : (ma.getdata(1000*(raw_data.variables[raw_data.variables['Pressure_surface'].dimensions[2]][:]))),
                'data_latnum' : (ma.getdata(1000*(raw_data.variables[raw_data.variables['Pressure_surface'].dimensions[1]][:]))).shape[0],
                'data_lonnum' : (ma.getdata(1000*(raw_data.variables[raw_data.variables['Pressure_surface'].dimensions[2]][:]))).shape[0]
                }
   

    # NOW THAT DATA IS PARSED AND IN THE SAME FORMAT, CLEAN IT UP:    
    # define a few functions needed here
    def get_sfc_index(height_arr):
        i = 0
        while i < len(height_arr):
            if height_arr[i] >= 0:
                return i
            else:
                i += 1
        return -1
    
    def make_sfc_based(arr, sfc_val, sfc_index):
        mod_arr = np.empty(len(arr) - sfc_index)
        i = 0
        while i < len(mod_arr):
            mod_arr[i] = arr[i + sfc_index]
            i = i + 1
        mod_arr = np.insert(mod_arr, 0, sfc_val)
        return mod_arr
    
    # create a mesh grid to find a single point in the data 
    clons, clats = np.meshgrid(latlon_data['data_lon'], latlon_data['data_lat'])
    tlons = clons 
    tlats = clats
    
    # find a single point using the haversine function 
    def haversine(lon_arr, lat_arr, lon_point, lat_point): 
        distance = 6371*(2*np.arcsin(np.sqrt(np.sin((np.radians(lat_arr)-np.radians(lat_point))/2)**2 + np.cos(np.radians(lat_point)) * np.cos(np.radians(lat_arr)) * np.sin((np.radians(lon_arr)-np.radians(lon_point))/2)**2)))
        return distance
    location_index = np.where(haversine(tlons,tlats,latlon1[1],latlon1[0]) == np.min(haversine(tlons,tlats,latlon1[1],latlon1[0])))
    
    
    # declare empty dicts 
    sfc_loc_idx  = {}
    vert_loc_idx = {}
    
    # fill dicts with data after being 'location-indexed'
    for key in sfc_data.keys():
        sfc_loc_idx[key] = sfc_data[key][location_index[0],location_index[1]][0]
    for key in ['vert_T', 'vert_z', 'vert_rh', 'vert_u', 'vert_v', 'vert_Td']:
        vert_loc_idx[key] = vert_data[key][::-1,location_index[0],location_index[1]]
        vert_loc_idx['vert_p'] = vert_data['vert_p'][::-1]
    
        
    # find the surface of the vertical data and combine them 
    sfc_index  = get_sfc_index((vert_loc_idx['vert_z']-sfc_loc_idx['sfc_z']))
    surface_height = sfc_loc_idx['sfc_z']
    sfc_loc_idx['sfc_z'] = 0
    sb_dict = {}
    new_keys = ['T', 'Td', 'rh', 'u', 'v', 'z', 'p']
    sfc_keys = ['sfc_T', 'sfc_Td', 'sfc_rh', 'sfc_u', 'sfc_v', 'sfc_z', 'sfc_p']
    vert_keys = ['vert_T', 'vert_Td', 'vert_rh', 'vert_u', 'vert_v', 'vert_z', 'vert_p']
    
    
    # create a dict of surface-based data 
    for vert_key, sfc_key, new_key in zip(vert_keys, sfc_keys, new_keys):
        sb_dict[new_key] = make_sfc_based(vert_loc_idx[vert_key],  sfc_loc_idx[sfc_key],  sfc_index)
    

    # Interpolates data
    dz = 250 
    soundingtop_hght = sb_dict['z'][-1]
    toplvl      = int(soundingtop_hght/dz)*dz
    numlvls     = int(toplvl/dz)
    interp_lvls = np.linspace(0,toplvl,numlvls+1)
    
    # prepare new dicts 
    keys = ['T', 'Td', 'rh', 'u', 'v', 'z', 'p']
    units_list = ['degC', 'degC', 'percent', 'kt', 'kt', 'm', 'hPa']
    interp_dict = {}
    zeros_dict  = {}
    clean_data  = {}
    
    # create dict of clean data 
    for key in keys:
        interp_dict[key] = (interpolate.interp1d(sb_dict['z'], sb_dict[key]))
        zeros_dict[key]  = np.zeros((len(interp_lvls)))
        for zeros_arr in zeros_dict.values():
            zeros_dict[key][0] = sb_dict[key][0]
        for i in range(1,len(zeros_dict[key]),1):
            zeros_dict[key][i] = interp_dict[key](dz*i)
    for i, unit, key in zip(range(0, len(units_list)), units_list, keys):
        clean_data[key] = zeros_dict[key]*units(unit)
    clean_data['zAGL'] = clean_data['z'] + surface_height*units.m
    clean_data['site_info'] = {
                'site-id'   : 'no-site-id',
                'site-name' : 'no-site-name',
                'site-lctn' : 'no-site-location',
                'site-latlon' : [latlon1[0], latlon1[1]],
                'site-elv'  : surface_height,
                'source'    : 'MODEL REANALYSIS PROFILE',
                'model'     : source,
                'fcst-hour' : 'F00',
                'run-time'  : [year1, month1, day1, hour1],
                'valid-time': [year1, month1, day1, hour1]}
    
    
    print('> COMPLETE --------')
    elapsed_time = time.time() - st
    print('> RUNTIME:', time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
    return clean_data
#########################################################################################################







### ACARS_DATA CLASS ###   
#########################################################################################################

# define class
class acars_data:
    
    # init 
    def __init__(self, year, month, day, hour):
        self.year = year
        self.hour = hour
        self.month = month
        self.day = day
        self.hour = hour
        
    ### LIST PROFILES ###
    #########################################################################################################
    def list_profiles(self):
        st = time.time()
        print(f'> LIST ACARS PROFILES FUNCTION --\n---------------------------------')

        """
        function used to search through available ACARS profiles
        for a given date and hour. The user must use this function
        to determine what ACARS profiles can be used for getting 
        data.

        COPYRIGHT
        Created by Kyle J Gillett (@wxkylegillett) 2023

        """

        # SET UP OU DIRECTORY REF
        # SEARCH FOR DIRECTORY FOR THE USER-SPECIFIED DATE
        data_dir = f'https://sharp.weather.ou.edu//soundings//acars//{self.year}//{self.month}//{self.day}//{self.hour}'

        # ACCESS THE RAW WEBSITE HTML & FIND THE ID_TIME KEYS 
        print(f"> AVAILABLE ACARS PROFILES FOR {self.year}-{self.month}-{self.day} {self.hour}Z...")
        # SET UP BEAUTIFUL SOUP TO PARSE HTML 
        # THIS WORKS AS A SORT OF JERRY-RIGGED WAY
        # TO REVEAL ALL AVAILABLE ACARS PROFILES
        # FOR A GIVEN DATE/TIME
        body = urlopen(data_dir).read().decode("utf-8")
        soup = bs4.BeautifulSoup(body, features="html.parser")

        # ADD PROFILES TO A LIST
        profiles_list = []
        for link in soup.select('a[href$=".txt"]'):
            profiles_list.append(link.get("href")[0:8]) 
        
        
        print('> COMPLETE --------')
        elapsed_time = time.time() - st
        print('> RUNTIME:', time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))

        return profiles_list
    #########################################################################################################


    ### GET PROFILE DATA ###
    #########################################################################################################
    def get_profile(self, acars_profile):
        
        """
        function used to access the raw data of a given
        ACARS profile.

        COPYRIGHT
        Created by Kyle J Gillett (@wxkylegillett) 2023
        """
        st = time.time()
        print(f'> ACARS DATA ACCESS FUNCTION --\n---------------------------------')
        
        # SET UP OU DIR REF TO THE SPECIFIC PROFILE
        profile_url    = f'https://sharp.weather.ou.edu//soundings//acars//{self.year}//{self.month}//{self.day}//{self.hour}//{acars_profile}.txt'
        
        # SEPERATE DATA BETWEEN HEADER AND ACTUAL DATA 
        try:
            data   = loadtxt(urlopen(profile_url).readlines()[6:-1], dtype='str', comments="%", unpack=True)
            header = loadtxt(urlopen(profile_url).readlines()[0:3], dtype='str', comments="%", unpack=True)
        except HTTPError as err:
            if err.code == 404:
                sys.exit('! ERROR ! -- Invalid profile, try again with a valid profile (ex: BNA_2320)')
            else:
                raise
        
        # PARSE DATE INFO FROM OU FILE
        year  = f'20{header[1][0:2]}'
        month = header[1][2:4]
        day   = header[1][4:6]
        hour  = f'{header[1][7:9]}:{header[1][9:11]}'

        # PARSE PROFILE DATA FROM OU FILE IN DICT
        new_keys = ['p', 'z', 'T', 'Td', 'u', 'v']
        units_list = ['hPa', 'meter', 'degC', 'degC']
        clean_data = {}
        for new_key, idx, unit in zip (new_keys, range(0,4), units_list):
            clean_data[new_key] = np.array([float(ele) for ele in [ele[0:-1] for ele in data[idx]]])*units(unit)
            
        clean_data['u'], clean_data['v'] = mpcalc.wind_components([float(ele) for ele in [ele[0:-1] for ele in data[5]]]*units.kts,
                                                             [float(ele) for ele in [ele[0:-1] for ele in data[4]]]*units.deg)
        
        # GET AIRPORT INFO FROM GITHUB AIRPORTS.CSV
        airports_csv = pd.read_csv(f'https://raw.githubusercontent.com/kylejgillett/sounderpy/main/src/AIRPORTS.csv',
                skiprows=7, skipinitialspace = True)
        where = [np.where(airports_csv['IATA'].str.contains(header[0], na=False, case=True))[0]][0][0]
        
        # ADD AIRPORT DATA INTO DICT
        keys = ['Name', 'City', 'Country', 'Latitude', 'Longitude', 'Altitude']
        airport_info = []
        for key in keys:
            airport_info.append(airports_csv[key][where])
        clean_data['site_info'] = {
                        'site-id'   : header[0],
                        'site-name' : airport_info[0],
                        'site-lctn' : airport_info[2],
                        'site-latlon' : [np.round(airport_info[3],2), np.round(airport_info[4],2)],
                        'site-elv'  : str(int(airport_info[5])),
                        'source'    : 'ACARS OBSERVED AIRCRAFT PROFILE',
                        'model'     : 'no-model',
                        'fcst-hour' : 'no-fcst-hour',
                        'run-time'  : ['no-run-time'],
                        'valid-time': [f'20{header[1][0:2]}', header[1][2:4], header[1][4:6], f'{header[1][7:9]}:{header[1][9:11]}']}

        print('> COMPLETE --------')
        elapsed_time = time.time() - st
        print('> RUNTIME:', time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
        return clean_data
    #########################################################################################################
    
#########################################################################################################











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
    
    # DEFINE find_nearest FUNCTION ----------------------------------------------- 
    def find_nearest(array, value):
        array = np.asarray(array)
        nearest_idx = (np.abs(array - value)).argmin()
        return nearest_idx
    
    
    # declare basic variables
    p = clean_data['p']
    T = clean_data['T']
    Td = clean_data['Td']
    z = clean_data['z']
    u = clean_data['u']
    v = clean_data['v'] 
    
    Td[Td < -120*units.degC] = np.nan
    
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


    # CREATE THE METPY SKEWT FIGURE ----------------------------------------------- 
    fig = plt.figure(figsize=(18,12))                             
    skew = SkewT(fig, rotation=43, rect=(0, 0, 0.6, 1)) 
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
    ml_t, ml_td = mpcalc.mixed_layer(p, T, Td)
    ml_p, _, _ = mpcalc.mixed_parcel(p, T, Td)
    # most unstable parcel properties!
    mu_p, mu_t, mu_td, _ = mpcalc.most_unstable_parcel(p, T, Td)
    # Compute parcel profiles
    sb_prof = mpcalc.parcel_profile(p, T[0], Td[0]).to('degC')
    mu_prof = mpcalc.parcel_profile(p, mu_t, mu_td).to('degC')
    ml_prof = mpcalc.parcel_profile(p, ml_t, ml_td).to('degC')
    # compute CAPE & CIN
    mlcape, mlcin = mpcalc.cape_cin(p, T, Td, ml_prof, which_lfc='bottom', which_el='top')
    mucape, mucin = mpcalc.cape_cin(p, T, Td, mu_prof, which_lfc='bottom', which_el='top')
    sbcape, sbcin = mpcalc.cape_cin(p, T, Td, sb_prof, which_lfc='bottom', which_el='top')
    
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
    imgax = fig.add_axes([-0.01, 0.84, 0.1, 0.1], anchor='SE')
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
    
    
    
    
    
    

    
    
    
    
### METPY_HODOGRAPH() FUNCTION ###   
#########################################################################################################

def __metpy_hodograph(clean_data):
    
    # record process time 
    st = time.time()   
    
    '''
        metpy_hodograph(clean_data)

        this function inputs cleaned profile data through a slightly modified MetPy sounding plot script. 

        COPYRIGHT
        Created by Kyle J Gillett (@wxkylegillett) 2023
    '''
    
    print(f'> HODOGRAPH PLOTTER FUNCTION --\n---------------------------------')
    
    # DEFINE find_nearest FUNCTION ----------------------------------------------- 
    def find_nearest(array, value):
        array = np.asarray(array)
        nearest_idx = (np.abs(array - value)).argmin()
        return nearest_idx

    # declare basic variables
    p = clean_data['p']
    T = clean_data['T']
    Td = clean_data['Td']
    z = clean_data['z']
    u = clean_data['u']
    v = clean_data['v'] 
    
    Td[Td < -120*units.degC] = np.nan
    
    # Create a hodograph
    fig = plt.figure(figsize=(16, 12), edgecolor="#04253a")
    ax = plt.subplot()
    h = Hodograph(ax, component_range=100)
    
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
   
    try:
        h.ax.set_xlim(x_Minlimit, x_Maxlimit)                                  
        h.ax.set_ylim(y_Minlimit, y_Maxlimit)                             
    except:
        h.ax.set_xlim(-65,65)
        h.ax.set_ylim(-65,65)
        pass
    
    h.add_grid(increment=20, color='black', linestyle='-', linewidth=1.5, alpha=0.3)   
    h.add_grid(increment=10, color='silver', linewidth=1, linestyle='--', alpha=0.5) 
    fig.set_facecolor('#ffffff')        
    h.ax.set_facecolor('#ffffff')
    h.ax.set_facecolor('#ffffff')                                      
    h.ax.set_box_aspect(1)                          
    h.ax.set_yticklabels([])
    h.ax.set_xticklabels([])
    h.ax.set_xticks([])
    h.ax.set_yticks([])
    h.ax.set_xlabel(' ')
    h.ax.set_ylabel(' ')
    # plot the hodograph itself, using plot_colormapped, colored
    # by height 
    h.plot_colormapped(u_hodo, v_hodo, c=z_hodo, linewidth=6, label="0-9km WIND")

    try:
        # compute Bunkers storm motion so we can plot it on the hodograph! 
        RM, LM, MW = mpcalc.bunkers_storm_motion(p, u, v, z)
        h.ax.text((RM[0].m), (RM[1].m), 'RM', weight='bold', ha='center', fontsize=13, alpha=0.6) 
        h.ax.text((LM[0].m), (LM[1].m), 'LM', weight='bold', ha='center', fontsize=13, alpha=0.6) 
        h.ax.text((MW[0].m), (MW[1].m), 'MW', weight='bold', ha='center', fontsize=13, alpha=0.6) 
        h.ax.arrow(0,0,RM[0].m, RM[1].m, linewidth=2, color='black', alpha=0.2, label='Bunkers RM Vector', 
                   length_includes_head=True, head_width=1)
        
    except:
        pass
    
    try:
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

    # Here we can add a simple python for loop that adds tick marks to the inside 
    # of the hodograph plot to increase readability! 
    plt.xticks(np.arange(0,0,1))
    plt.yticks(np.arange(0,0,1))
    for i in range(10,120,10):
        h.ax.annotate(str(i),(i,0),xytext=(0,2),textcoords='offset pixels',clip_on=True,fontsize=10,weight='bold',alpha=0.3,zorder=0)
    for i in range(10,120,10):
        h.ax.annotate(str(i),(0,i),xytext=(0,2),textcoords='offset pixels',clip_on=True,fontsize=10,weight='bold',alpha=0.3,zorder=0)
    for i in range(10,120,10):
        h.ax.annotate(str(i),(-i,0),xytext=(0,2),textcoords='offset pixels',clip_on=True,fontsize=10,weight='bold',alpha=0.3,zorder=0)
    for i in range(10,120,10):
        h.ax.annotate(str(i),(0,-i),xytext=(0,2),textcoords='offset pixels',clip_on=True,fontsize=10,weight='bold',alpha=0.3,zorder=0)
    
    every = int(len(z)/8)
    for zi, ui, vi in zip(z_hodo[2 :: every], u_hodo[2 :: every], v_hodo[2 :: every]): 
        h.plot(ui,vi, marker='.', color='white', alpha=0.8, markersize=30)
        h.ax.annotate(str(int(zi.m)),(ui,vi), weight='bold', fontsize=11, color='black',xytext=(0,-3.2),textcoords='offset pixels',horizontalalignment='center',clip_on=True) 
    
    # First we want to actually add values of data to the plot for easy viewing
    # to do this, lets first add a simple rectangle using matplotlib's 'patches' 
    # fucntionality to add some simple layout for plotting calculated parameters 
    #                                  xloc   yloc   xsize  ysize
    fig.patches.extend([plt.Rectangle((0.8, 0.125), 0.25, 0.755,
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
    # compute CAPE & CIN
    mlcape, mlcin = mpcalc.cape_cin(p, T, Td, ml_prof, which_lfc='bottom', which_el='top')
    mucape, mucin = mpcalc.cape_cin(p, T, Td, mu_prof, which_lfc='bottom', which_el='top')
    sbcape, sbcin = mpcalc.cape_cin(p, T, Td, sb_prof, which_lfc='bottom', which_el='top')

    lcl_pressure, lcl_temperature = mpcalc.lcl(p[0], T[0], Td[0])
    lfc_pressure, lfc_temperature = mpcalc.lfc(p, T, Td, which='most_cape')
    el_pressure, el_temperature = mpcalc.el(p, T, Td, which='most_cape')

    # make minor corrections to CAPE amount to match SHARPPY (SPC) output
    sbcape = (sbcape.m + (sbcape.m*(0.10)))*units('J/kg')
    mucape = (mucape.m + (mucape.m*(0.10)))*units('J/kg')
    mlcape = (mlcape.m + (mlcape.m*(0.10)))*units('J/kg')
    
    new_p = np.append(p[p > lcl_pressure], lcl_pressure)
    new_t = np.append(T[p > lcl_pressure], lcl_temperature)
    lcl_height = mpcalc.thickness_hydrostatic(new_p, new_t)
    
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
    
    # now some kinematic parameters
    met_per_sec = (units.m*units.m)/(units.sec*units.sec)
    plt.figtext( 0.83, 0.83,  f'0-1km SRH: ', weight='bold', fontsize=15, color='black', ha='left')
    plt.figtext( 1.01, 0.83,  f'{mag_int(total_helicity1)* met_per_sec:~P}', weight='bold', fontsize=15, color='navy', ha='right')
    plt.figtext( 0.83, 0.80,  f'0-1km SHEAR: ', weight='bold', fontsize=15, color='black', ha='left')
    plt.figtext( 1.01, 0.80,  f'{mag_int(bshear1)} kts', weight='bold', fontsize=15, color='blue', ha='right')
    plt.figtext( 0.83, 0.75,  f'0-3km SRH: ', weight='bold', fontsize=15, color='black', ha='left')
    plt.figtext( 1.01, 0.75,  f'{mag_int(total_helicity3)* met_per_sec:~P}', weight='bold', fontsize=15, color='navy', ha='right')
    plt.figtext( 0.83, 0.72,  f'0-3km SHEAR: ', weight='bold', fontsize=15, color='black', ha='left')
    plt.figtext( 1.01, 0.72,  f'{mag_int(bshear3)} kts', weight='bold', fontsize=15, color='blue', ha='right')
    plt.figtext( 0.83, 0.67,  f'0-6km SRH: ', weight='bold', fontsize=15, color='black', ha='left')
    plt.figtext( 1.01, 0.67,  f'{mag_int(total_helicity6)* met_per_sec:~P}', weight='bold', fontsize=15, color='navy', ha='right')
    plt.figtext( 0.83, 0.64,  f'0-6km SHEAR: ', weight='bold', fontsize=15, color='black', ha='left')
    plt.figtext( 1.01, 0.64,  f'{mag_int(bshear6)} kts', weight='bold', fontsize=15, color='blue', ha='right')
    plt.figtext( 0.83, 0.59,  f'SIG TORNADO: ', weight='bold', fontsize=15, color='black', ha='left')
    plt.figtext( 1.01, 0.59,  f'{mag_int(sig_tor)}', weight='bold', fontsize=15, color='orangered', ha='right')
    plt.figtext( 0.83, 0.56,  f'SUPERCELL COMP: ', weight='bold', fontsize=15, color='black', ha='left')
    plt.figtext( 1.01, 0.56,  f'{mag_int(super_comp)}', weight='bold', fontsize=15, color='orangered', ha='right')
    
    plt.figtext( 0.83, 0.51,  f'SBCAPE: ', weight='bold', fontsize=15, color='black', ha='left')
    plt.figtext( 1.01, 0.51,  f'{mag_int(sbcape)} J/kg', weight='bold', fontsize=15, color='orangered', ha='right')
    plt.figtext( 0.83, 0.48,  f'SBCIN: ', weight='bold', fontsize=15, color='black', ha='left')
    plt.figtext( 1.01, 0.48,  f'{mag_int(sbcin)} J/kg', weight='bold', fontsize=15, color='lightblue', ha='right')
    plt.figtext( 0.83, 0.43,  f'MLCAPE: ', weight='bold', fontsize=15, color='black', ha='left')
    plt.figtext( 1.01, 0.43,  f'{mag_int(mlcape)} J/kg', weight='bold', fontsize=15, color='orangered', ha='right')
    plt.figtext( 0.83, 0.40,  f'MLCIN: ', weight='bold', fontsize=15, color='black', ha='left')
    plt.figtext( 1.01, 0.40,  f'{mag_int(mlcin)} J/kg', weight='bold', fontsize=15, color='lightblue', ha='right')
    plt.figtext( 0.83, 0.35,  f'MUCAPE: ', weight='bold', fontsize=15, color='black', ha='left')
    plt.figtext( 1.01, 0.35,  f'{mag_int(mucape)} J/kg', weight='bold', fontsize=15, color='orangered', ha='right')
    plt.figtext( 0.83, 0.32,  f'MUCIN: ', weight='bold', fontsize=15, color='black', ha='left')
    plt.figtext( 1.01, 0.32,  f'{mag_int(mucin)} J/kg', weight='bold', fontsize=15, color='lightblue', ha='right')
    plt.figtext( 0.83, 0.27,  f'TT-INDEX: ', weight='bold', fontsize=15, color='black', ha='left')
    plt.figtext( 1.01, 0.27,  f'{mag_int(total_totals)} Δ°C', weight='bold', fontsize=15, color='orangered', ha='right')
    plt.figtext( 0.83, 0.24,  f'K-INDEX: ', weight='bold', fontsize=15, color='black', ha='left')
    plt.figtext( 1.01, 0.24,  f'{mag_int(kindex)} °C', weight='bold', fontsize=15, color='orangered', ha='right')
    
    plt.figtext( 0.83, 0.19,  f'ECAPE: ', weight='bold', fontsize=15, color='black', ha='left')
    plt.figtext( 1.01, 0.19,  f'{mag_int(ecape)} J/kg', weight='bold', fontsize=15, color='orangered', ha='right')
        
    # PLOT TITLES ----------------------------------------------------------------------------------
    plt.figtext( 0.23, 0.89, left_title, ha='left', weight='bold', fontsize=16)
    plt.figtext( 1.07, 0.89, f'{right_title}     ', ha='right', weight='bold', fontsize=16)
    plt.figtext( 0.23, 0.91, top_title, ha='left', weight='bold', fontsize=22) 
    plt.figtext( 1.07, 0.1,  f'POWERED BY METPY & SOUNDERPY | (@WXKYLEGILLETT)     ', ha='right', color='blue', alpha=0.4, weight='bold', fontsize=12)
    
    img = Image.open(urlopen('https://user-images.githubusercontent.com/100786530/251580013-2e9477c9-e36a-4163-accb-fe46780058dd.png'))
    #                  side-side  up-down  size   size
    imgax = fig.add_axes([0.23, 0.13, 0.1, 0.1], anchor='SE')
    imgax.imshow(img)
    imgax.axis('off')
    
    hodoleg = h.ax.legend(loc='upper left')
    
    elapsed_time = time.time() - st
    print('> RUNTIME:', time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
    
    return plt


def metpy_hodograph(clean_data, method='show', filename='sounderpy_hodograph'):
    if method == 'show':
        __metpy_hodograph(clean_data).show()
    elif method == 'save':
        __metpy_hodograph(clean_data).savefig(filename, bbox_inches='tight')
    
#########################################################################################################    
    
    





    

    
    
### SOUNDERPY HELPER FUCNTIONS  ###                                                            
#########################################################################################################
'''
     
     A collection of helper fucntions that users may find useful for processing vertical profile data but are not nessacary to use the basic 
     functions of SounderPy.
     
'''
    
####################################################### DATA INTERPOLATION ################################################################# 
def interp_data(variable, heights, step=100):
    
    '''
        interpolate vertical profile over vertical height array
    '''
    
    try:
        variable.units
        variable = variable.m
    except:
        variable = variable
    try:
        heights.units
        heights = heights.m
    except:
        heights = heights
        
    levels=np.arange(0,np.max(heights),step)
    varinterp=np.zeros(len(levels))
    for i in range(0,len(levels)):
        lower=np.where(heights-levels[i]<=0,heights-levels[i],-np.inf).argmax()
        varinterp[i]=(((variable[lower+1]-variable[lower])/(heights[lower+1]-heights[lower]))*(levels[i]-heights[lower])+variable[lower])
    return varinterp 

    
########################################################## FIND NEAREST ####################################################################### 
def find_nearest(array, value):
    
    '''
        search through an array to find the index of the value nearest to a given value
        
    '''
    
    array = np.asarray(array)
    nearest_idx = (np.abs(array - value)).argmin()
    return nearest_idx

    

    
##################################################### MAKE SURFACED BASED #################################################################### 

def get_sfc_index(height_arr):
    
    """
        This function takes an array of height values and finds the index before the first positive value

    """
    
    i = 0
    # Search the array for the sfc index
    while i < len(height_arr):
        if height_arr[i] >= 0:
            return i
        else:
            i += 1
    # Did not find a positive index
    return -1




def make_sfc_based(arr, sfc_val, sfc_index):
    
    """
        takes an array and a valid index in that array, then returns a copy of the
        array beginning at the provided index i.e., chops off below-ground values

    """

    # Initialize an empty numpy array as the modified array
    mod_arr = np.empty(len(arr) - sfc_index)
    # Insert the sfc and higher values into the modified array
    i = 0
    while i < len(mod_arr):
        mod_arr[i] = arr[i + sfc_index]
        i = i + 1
    # Inserts surface values
    mod_arr = np.insert(mod_arr, 0, sfc_val)
    return mod_arr
    

    
    
def make_sfc_based_3D(arr, sfc_arr):
    
    '''
        takes a 3D array of mandatory level data and a 2D array of surface data, 
        appends the surface data onto the mandatory level array, and returns a single array
        of both surface and mandatory level data
        
    '''
    
    mod_arr = np.zeros((np.shape(arr)[0]+1,np.shape(arr)[1],np.shape(arr)[2]))
    for j in range(np.shape(arr)[1]):
        for k in range(np.shape(arr)[2]):
            mod_arr[0,j,k] = sfc_arr[j,k]
            mod_arr[1:,j,k] = arr[:,j,k]
    return mod_arr


    
    
    
    
########################################################### TO CM1 FILE ####################################################################     

def to_file(file_type, clean_data, filename=None):
    
    
    # set file name 
    if filename == None:
        filename = f'sounderpy_data'
    else:
        filename = filename
        
        
        
    ####################################### CM1 #######################################
    if file_type == 'cm1':
        '''
        function, creates CM1 input sounding file for CM1 integration
        COPYRIGHT
        Created by Kyle J Gillett (@wxkylegillett) 2023 & Kelton Halbert 2021
        
        Derived from Kelton Halbert / Leigh Orf via github: 
        https://github.com/leighorf/LOFS-read/blob/master/bin/sndmod
        '''
            
        # create file    
        outfile = open(filename, 'w')
        num_lines = len(list(clean_data.items())[0][1])
        delimiter=''

        # use metpy to find parameters that CM1 likes
        clean_data['theta'] = mpcalc.potential_temperature(clean_data['p'], clean_data['T'])
        clean_data['relhm'] = mpcalc.relative_humidity_from_dewpoint(clean_data['T'], clean_data['Td'])
        clean_data['mixrt'] = mpcalc.mixing_ratio_from_relative_humidity(clean_data['p'], clean_data['T'], clean_data['relhm'])

        # add data to lines 
        for idx in range(num_lines):
                line_str = ""
                line_str += "%12s" % str(format(np.around(clean_data["z"][idx].m, 6), "0.6f")) + delimiter + str("\t")
                line_str += "%12s" % str(format(np.around(clean_data["theta"][idx].m, 6), "0.6f")) + delimiter + str("\t")
                line_str += "%12s" % str(format(np.around(clean_data["mixrt"][idx].m, 6), "0.6f")) + delimiter + str("\t")
                line_str += "%12s" % str(format(np.around(clean_data["u"][idx].m, 6), "0.6f")) + delimiter + str("\t")
                line_str += "%12s" % str(format(np.around(clean_data["v"][idx].m, 6), "0.6f")) + str("\n")
                outfile.write(line_str)

        outfile.close()
        
        
        
    ####################################### CSV #######################################
    elif file_type == 'csv':
        '''
        function, creates CSV file of sounding data

        COPYRIGHT
        Created by Kyle J Gillett (@wxkylegillett) 2023
        '''

        # remove units from data
        no_units = {}
        for key in ['p', 'z', 'T', 'Td', 'u', 'v']:
                no_units[key] = clean_data[key].m
        # open and write to CSV 
        with open(filename, "w") as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(no_units.keys())
            writer.writerows(zip(*no_units.values()))   
            
            
            
    ####################################### SHARPPY #######################################    
    elif file_type == 'sharppy':
        '''
        function, creates NSHARP input sounding file for SharpPy integration
        COPYRIGHT
        Created by Kyle J Gillett (@wxkylegillett) 2023 & Kelton Halbert 2021

        Derived from Kelton Halbert / Leigh Orf via github: 
        https://github.com/leighorf/LOFS-read/blob/master/bin/sndmod
        '''


        outfile_file = open(filename, 'w')

        outfile_loc = ("****")

        dt = datetime(int(clean_data['site_info']['valid-time'][0]), int(clean_data['site_info']['valid-time'][1]), 
                       int(clean_data['site_info']['valid-time'][2]), int(clean_data['site_info']['valid-time'][3][0:2]))

        outfile_file.write("%TITLE%\n")
        outfile_file.write("%s   %s\n" % (clean_data['site_info']['site-id'], dt.strftime("%y%m%d/%H%M")))
        outfile_file.write("   LEVEL       HGHT       TEMP       DWPT       WDIR       WSPD\n")
        outfile_file.write("-------------------------------------------------------------------\n")
        outfile_file.write("%RAW%\n")

        ws = mpcalc.wind_speed(clean_data['u'], clean_data['v'])
        wd = mpcalc.wind_direction(clean_data['u'], clean_data['v'])

        new_data = {
            'p' : clean_data['p'],
            'z' : clean_data['z'],
            'T' : clean_data['T'],
            'Td': clean_data['Td'],
            'wd': wd,
            'ws': ws,
        }

        for idx in range(new_data['p'].shape[0]):
            string = ""
            for col in ['p', 'z', 'T', 'Td', 'wd', 'ws']:
                string += "%12.6f,  " % new_data[col][idx].m

            outfile_file.write(string[:-3] + "\n")
        outfile_file.write("%END%\n")
        outfile_file.close()    

##########################################################################################################################################    
    

    
    
    
    
    
    
  
    
    
########################################################### get_latlon ####################################################################    
 
def get_latlon(station_type, id):
    
    # valid methods include 'METAR', 'RAOB', 'IGRA', 'BUOY', 'BUFKIT'
    
    # DMS to decimal degrees
    def dms2dd_min(degrees, minutes, direction):
        dd = float(degrees) + float(minutes) / 60
        if direction == "S" or direction == "W":
            dd *= -1
        return dd
    
    def dms2dd(degrees, direction):
        dd = float(degrees)
        if direction == "S" or direction == "W":
            dd *= -1
        return dd
    
    

    ############################################### METAR #####################################################
    if station_type.casefold() == 'metar':
        '''
        takes a METAR site id such as 'KMBS', searches over 9000 station IDs  and returns
        a list including the lat/lon for the METAR site
        '''

        request = requests.get("http://aviationweather.gov/docs/metar/stations.txt", stream=True)
        stations = {}
        for line in request.iter_lines():
            data = line.decode("ascii")
            if data:
                if data[0] == "!" or len(data) != 83:
                    continue
                province = data[0:2]
                station = data[3:19].strip()
                icao = data[20:24].strip()
                lat = dms2dd_min(data[39:41], data[42:44], data[44:45])
                lon = dms2dd_min(data[47:50], data[51:53], data[53:54])
                altitude = int(data[55:59])
                country = data[81:83]

                if icao:
                    stations[icao] = {"name": station, "lat": lat, "lon": lon, "altitude": altitude, "country": country}
        try:
            latlon = [np.round(stations.get(id)['lat'],2), np.round(stations.get(id)['lon'],2)]
            return latlon
        except:
            sys.exit("NO METAR SITE FOUND FOR THAT ID -- Make sure you entered a valid METAR ID -- Make sure you entered the METAR ID as a string, ex: 'KMOP'")
            pass
        
        
        
    ############################################### BUFKIT #####################################################    
    elif station_type.casefold() == 'bufkit': 
        '''
        takes a BUFKIT site id such as 'KMOP', searches over 1200 station IDs and returns
        a list including the lat/lon for the BUFKIT site
        '''
        
        BUFKIT_STATIONS = pd.read_csv(f'https://raw.githubusercontent.com/kylejgillett/sounderpy/main/src/BUFKIT-STATIONS-MASTER.txt', skiprows=7, skipinitialspace = True)
        
        try:
            station = BUFKIT_STATIONS['ID'][np.where(BUFKIT_STATIONS['ID'].str.contains(id, na=False, case=True))[0]].values[0]

            lat = (BUFKIT_STATIONS[BUFKIT_STATIONS['ID']==station]['LAT'].values[0])
            lon = (BUFKIT_STATIONS[BUFKIT_STATIONS['ID']==station]['LON'].values[0])
            return [lat, lon]

        except:
            sys.exit("NO BUFKIT SITE FOUND FOR THAT ID -- Make sure you entered a valid BUFKIT ID -- Make sure you entered the BUFKIT ID as a string, ex: 'KMOP'")
            pass
        
        
        
    ############################################### RAOB #####################################################    
    elif station_type.casefold() == 'raob':
        
        '''
        takes a RAOB site id such as 'DTX', searches over 9000 station IDs  and returns
        a list including the lat/lon for the RAOB site
        '''
        # add raob stations list 
        RAOB_STATIONS = pd.read_csv(f'https://raw.githubusercontent.com/kylejgillett/sounderpy/main/src/RAOB-STATIONS.txt', skiprows=7, skipinitialspace = True)

        # find lat-lon from stations list if it exists 
        try:
            station = RAOB_STATIONS['ICAO'][np.where(RAOB_STATIONS['ICAO'].str.contains(id, na=False, case=True))[0]].values[0]
            lat = dms2dd(RAOB_STATIONS[RAOB_STATIONS['ICAO']==station]['LAT'].values[0], RAOB_STATIONS[RAOB_STATIONS['ICAO']==station]['A'].values[0])
            lon = dms2dd(RAOB_STATIONS[RAOB_STATIONS['ICAO']==station]['LON'].values[0], RAOB_STATIONS[RAOB_STATIONS['ICAO']==station]['B'].values[0])
            return [lat, lon]
        except:
            try: 
                station = RAOB_STATIONS['WMO'][RAOB_STATIONS[RAOB_STATIONS['WMO']==int(id)].index[0]]
                lat = dms2dd(RAOB_STATIONS[RAOB_STATIONS['WMO']==station]['LAT'].values[0], RAOB_STATIONS[RAOB_STATIONS['WMO']==station]['A'].values[0])
                lon = dms2dd(RAOB_STATIONS[RAOB_STATIONS['WMO']==station]['LON'].values[0], RAOB_STATIONS[RAOB_STATIONS['WMO']==station]['B'].values[0])
                return [lat, lon]
            except:
                sys.exit("NO RAOB SITE FOUND FOR THAT ID -- Make sure you entered a valid RAOB ID -- Make sure you entered the RAOB ID as a string, ex: 'DTX' or '72632'")
                pass
    
    
    
    ############################################### IGRA #####################################################
    elif station_type.casefold() == 'igra':
        '''
        takes a IGRA2 site id such as 'GMM00010393', searches nearly 3000 station IDs  and returns
        a list including the lat/lon for the IGRA2 site
        '''

        # get IGRA stations list
        IGRA_STATIONS = pd.read_csv(f'https://raw.githubusercontent.com/kylejgillett/sounderpy/main/src/IGRA-STATIONS.txt', skiprows=7, skipinitialspace = True)

        # find lat-lon from stations list if it exists 
        try:
            station = IGRA_STATIONS['ID'][np.where(IGRA_STATIONS['ID'].str.contains(id, na=False, case=True))[0]].values[0]

            lat = np.round(IGRA_STATIONS[IGRA_STATIONS['ID']==station]['LAT'].values[0], 2)
            lon = np.round(IGRA_STATIONS[IGRA_STATIONS['ID']==station]['LON'].values[0], 2)
            return [lat, lon]

        except:
            sys.exit("NO IGRAv2 SITE FOUND FOR THAT ID -- Make sure you entered a valid IGRAv2 ID -- Make sure you entered the RAOB ID as a string, ex: 'GMM00010393'")
            pass
    
    
    
    ############################################### BUOY #####################################################
    elif station_type.casefold() == 'buoy':
        '''
        takes a BUOY/CMAN site id such as '41001', searches through a number of station IDs and returns
        a list including the lat/lon for the BUOY/CMAN  site
        '''

        # get buoy stations list
        BUOY_STATIONS = pd.read_csv('https://raw.githubusercontent.com/kylejgillett/sounderpy/main/src/BUOY-STATIONS.txt', skiprows=7, skipinitialspace = True)

        # find lat-lon from stations list if it exists
        try:
            station = BUOY_STATIONS['ID'][np.where(BUOY_STATIONS['ID'].str.contains(id, na=False, case=True))[0]].values[0]

            lat = dms2dd(BUOY_STATIONS[BUOY_STATIONS['ID']==station]['LAT'].values[0], BUOY_STATIONS[BUOY_STATIONS['ID']==station]['A'].values[0])
            lon = dms2dd(BUOY_STATIONS[BUOY_STATIONS['ID']==station]['LON'].values[0], BUOY_STATIONS[BUOY_STATIONS['ID']==station]['B'].values[0])
            return [lat, lon]

        except:
            sys.exit("NO BUOY SITE FOUND FOR THAT ID -- Make sure you entered a valid BUOY ID -- Make sure you entered the BUOY ID as a string, ex: '41001'")
            pass
        
    else:
        sys.exit("ERROR: you may have passed an incorrect station_type. Valid station_types are 'metar', 'raob', 'igra', 'bufkit', 'buoy'")  
        
        
##########################################################################################################################################    

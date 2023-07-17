import cdsapi
from sys import exit
import requests
from   metpy.units  import units
import metpy.calc   as mpcalc
import xarray as xr
from urllib.request import urlopen
import numpy as np
import numpy.ma       as ma
from siphon.catalog   import TDSCatalog
from siphon.ncss      import NCSS
from scipy import interpolate
import cartopy.crs as ccrs
from siphon.simplewebservice.wyoming import WyomingUpperAir
from siphon.simplewebservice.iastate import IAStateUpperAir
from datetime         import datetime
import time
import pandas as pd

version = 'v1.0.0'
version_date = 'July 17, 2023'


'''
        VERTICAL PROFILE DATA RETRIEVAL TOOL
        ------------------------------------
        
        This script it used to access vertical profile data for calculations or plotting of a vertical profile (sounding)
        This script is capable of accessing data from the ECMWF CDS (ERA5 reanalysis), the unidata THREDDS TDS (RAP reanalysis & real time),
        the University of Wyoming RAOB profile archive (observed), the Iowa State University's RAOB archive (observed) or a manually uploaded
        .txt file of sounding data.
        
        RELEASE
        -------
        Version: 1.0.0 | July, 17, 2023

        COPYRIGHT
        ---------
        Created by Kyle J Gillett (@wxkylegillett) 2023

'''



citation_text = f"""
## ------------------ VERTICAL PROFILE DATA RETRIEVAL TOOL ------------------------ ##
##                           {version} | {version_date}                                  ##
##                             (C) KYLE J GILLETT                                   ##
##         THIS TOOL CAN LOAD RAOB, RAP, RUC, ERA5 AND VERTICAL PROFILE DATA        ##
##  CALL THE FUNCTION '.get_docs()' TO PRINT SIMPLE DOCUMENTATION FOR THIS PACKAGE  ##
## -------------------- THANK YOU FOR USING THIS PACKAGE -------------------------- ##
"""
print(citation_text)










#################################################################################################################################################
#                                                      GET_MODEL_DATA() FUNCTION                                                                 
#################################################################################################################################################

def get_model_data(method, latlon, year, month, day, hour, domain='point'):
    st = time.time()
    
    '''
    get_model_data(method, domain, latlons, year, month, day, hour)

    function used to access model profile data from the ERA5, RAP or RUC

    PARAMETERS
    ----------
      name      dtype  required  example         explanation 
      method:    str,    yes,    'rap-ruc'       model data source: 'rap-ruc', 'era5' 
      latlon:    list,   yes,    [42.3, -97.3]   ordered list of lat, lon points 
      year:      str,    yes,    '2014'          4 digit year (UTC)
      month:     str,    yes,    '06'            2 digit month (UTC)
      day:       str,    yes,    '16'            2 digit day (UTC)
      hour:      str,    yes,    '20'            2 digit hour (UTC)
      domain:    str,    no,     'point'         default='point', single lat/lon 'point' or 'map' for lat/lon box *under dev.

    RETURNS
    -------
      raw_data: for ERA5, an xarray.Dataset of ERA5 reanalysis data or
      raw_data: for RAP/RUC, a netCDF4._netCDF4.Dataset of RAP/RUC reanalysis data

    EXAMPLES
    --------
      raw_data = spy.get_model_data('rap', 'point', [42.3, -97.3], '2014', '06', '16', '18')
      raw_data = spy.get_model_data('era5', 'point', [42.3, -97.3], '2014', '06', '16', '18')

    note: method 'rap-ruc' will try for RAP data first, if the data is not availible, then it will try for RUC data
    
    COPYRIGHT
    ---------
       Created by Kyle J Gillett (@wxkylegillett) 2023

    '''
    
    global latlon1, source, year1, month1, day1, hour1
    latlon1 = latlon
    source = str.upper(method)
    year1  = year
    month1 = month
    day1   = day
    hour1  = hour

                 #lat            #lat           #lon               #lon
    latlons = [latlon[0] + .5, latlon[0] - .5, latlon[1] - .5, latlon[1] + .5]

    
    ######################################################## ERA5 FUNCTION ###################################################################
    if method in ['era', 'era5', 'ERA', 'ERA5']:

        print(f'-- ERA5 REANALYSIS DATA ACCESS FUNCTION --\n------------------------------------------')

        dataset_presLvls = 'reanalysis-era5-pressure-levels'
        dataset_singleLvls = 'reanalysis-era5-single-levels'
        download_flag = 'false' 

        if domain == 'point':
            latlon_list = [latlons[0], latlons[2], latlons[1], latlons[3]]
        elif domain == 'map':
            latlon_list = [latlons[0]+10, latlons[2]-10, latlons[1]-10, latlons[3]+10]
        else:
            latlon_list = [9999, 9999, 9999, 9999] 
            print('[!] FUNCTION ERROR [!]\n[!] method str must be point, or map, this code will now end [!]')
            exit()

        # SET UP CALL FOR PRESSURE LEVEL DATA 
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

        # SET UP CALL FOR SINGLE LEVEL DATA 
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

        # retrieves the path to the file
        fl = c.retrieve(dataset_presLvls , params)# download the file
        print('DATASET ACCESSED: '+dataset_presLvls )
        fl2 = c2.retrieve(dataset_singleLvls, params2)
        print('DATASET ACCESSED: '+dataset_singleLvls )
        #if download_flag:
        fl.download("./output.nc")
        fl2.download("./output.nc")# load into memory

        with urlopen(fl.location) as f:
            ds = xr.open_dataset(f.read())
        with urlopen(fl2.location) as f:
            ds2 = xr.open_dataset(f.read())

        ds2 = ds2.rename({'z':'hgts','sp':'ps','t2m':'Ts','d2m':'tds','u10':'us','v10':'vs', })

        raw_data = xr.merge([ds,ds2])
        print('- COMPLETE -')
        elapsed_time = time.time() - st
        print('RUNTIME:', time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
        return raw_data



    ######################################################## RAP FUNCTION ###################################################################
    if method in ['rap', 'ruc', 'rap-ruc', 'RAP', 'RUC', 'RAP-RUC']:

        print(f'-- RAP REANALYSIS DATA ACCESS FUNCTION --\n-----------------------------------------')

        if domain == 'point':
            latlon_list = [latlons[2],latlons[3],latlons[1],latlons[0]]
        elif domain == 'map':
            latlon_list = [latlons[2]-10 ,latlons[3]+10,latlons[1]-10,latlons[0]+10]
        else:
            latlon_list = [9999, 9999, 9999, 9999] 
            print('[!] FUNCTION ERROR [!]\n[!] method str must be point, or map, this code will now end [!]')
            exit()

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
        try:
            for url, key in zip(urls.values(), urls.keys()):
                try:
                    NCSS(url)
                    print(f'DATASET USED: {key}')
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
            print('[!] ERROR [!]\nNCSS URL FAILED -- THIS HAPPENS WHEN A BAD REQUEST IS MADE. POSSIBLE ERRORS: \n> CHECK TO MAKE SURE YOU ENTERED THE CORRECT DATES.\n> NOTE: DATA IS NOT AVAILIABLE FOR EVERY DATE\n> THIS CATALOG OFTEN EXPERIENCES MISSING DATA/OUTAGES\n> MAKE SURE DATES ARE STRINGS -- MONTH, DAY AND HOUR MUST BE TWO DIGITS (EX: 18, 06, 00)')
            pass
        query = data.query()
        # Subsets by data
        if source == 'rap':
            query.variables('Pressure_surface',
                        'Geopotential_height_isobaric', 'Geopotential_height_surface',
                        'Temperature_isobaric', 'Temperature_height_above_ground',
                        'Relative_humidity_isobaric', 'Dewpoint_temperature_height_above_ground',
                        'Relative_humidity_height_above_ground','Vertical_velocity_pressure_isobaric',
                        'u-component_of_wind_height_above_ground', 'v-component_of_wind_height_above_ground', 
                        'u-component_of_wind_isobaric', 'v-component_of_wind_isobaric').add_lonlat()
        else:
            query.variables('Pressure_surface',
                        'Geopotential_height_isobaric', 'Geopotential_height_surface',
                        'Temperature_isobaric', 'Temperature_height_above_ground',
                        'Relative_humidity_isobaric','Dewpoint_temperature_height_above_ground',
                        'Relative_humidity_height_above_ground', 'Vertical_velocity_pressure_isobaric',
                        'u-component_of_wind_height_above_ground', 'v-component_of_wind_height_above_ground', 
                        'u-component_of_wind_isobaric', 'v-component_of_wind_isobaric').add_lonlat()
        # Subsets by area
        query.lonlat_box(latlon_list[0], latlon_list[1], latlon_list[2], latlon_list[3])
        # Gets data
        raw_data = data.get_data(query)

        print('- COMPLETE -')
        
        elapsed_time = time.time() - st
        print('RUNTIME:', time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
        return raw_data

    else:
        print("[!] ERROR [!]\n[!] You may have entered an invalid method.\nValid methods include 'RAP' or 'ERA5'")
    
    
    
    
    
    
    
##################################################################################################################################################
#                                                     GET_OBS_DATA() FUNCTION                                                                    #
##################################################################################################################################################

def get_obs_data(station, year, month, day, hour):
    st = time.time()
    """
    get_obs_data(station, year, month, day, hour)

    function used to access and parse RAOB profile data
    this function will search both UW & ISU for desired station and date

    PARAMETERS
    ----------
      name      dtype  required  example         explanation 
      station:   str,    yes,    'OAX'           3 digit RAOB station identifier
      year:      str,    yes,    '2014'          4 digit year (UTC)
      month:     str,    yes,    '06'            2 digit month (UTC)
      day:       str,    yes,    '16'            2 digit day (UTC)
      hour:      str,    yes,    '18'            2 digit hour (UTC)

    RETURNS
    -------
      clean_data: a dict of cleaned profile data ready to use for plotting

    EXAMPLES
    --------
      clean_data = spy.get_obs_data('OAX', '2014', '06', '16', '18')
        
    COPYRIGHT
    ---------
      Created by Kyle J Gillett (@wxkylegillett) 2023
    
    """
    
    global station1, source, year1, month1, day1, hour1
    source = 'obs'
    station1 = station
    year1  = year
    month1 = month
    day1   = day
    hour1  = hour
    RAOB_STATIONS = pd.read_csv(f'https://raw.githubusercontent.com/kylejgillett/sounderpy/main/RAOB-STATIONS.txt', skiprows=7, skipinitialspace = True)
    
    
    dt = datetime(int(year), int(month), int(day), int(hour)) 
    got_data = False
    for i in range(1, 11):
        try: 
            df = WyomingUpperAir.request_data(dt, station)
            got_data = True
            if got_data == True:
                print(f'FOUND RAOB: {station} on {month}/{day}/{year} at {hour}z')
                break
        except:
            if i == 10:
                print('NO DATA FOR THIS STATION & DATE IN THE UW ARCHIVE, TRYING IEM...')
                got_data = False          
        pass
    if got_data == False: 
        for i in range(1, 11):
            try: 
                df = IAStateUpperAir.request_data(dt, station)
                got_data = True
                if got_data == True:
                    print(f'FOUND RAOB: {station} on {month}/{day}/{year} at {hour}z')
                    break
            except:
                if i == 10:
                    print(f'[!] OPERATION FAILED [!]\n There is likely no available data for station {station} on {month}/{day}/{year} at {hour}z\nPlease make sure you entered a valid station ID and date\nor try a different time, date, or station')
                    got_data = False          
            pass
 
    if got_data == True:
        station = RAOB_STATIONS['ICAO'][np.where(RAOB_STATIONS['ICAO'].str.contains(station, na=False, case=True))[0]].values[0]
        
        new_keys = ['p', 'z', 'T', 'Td', 'u', 'v']
        old_keys = ['pressure', 'height', 'temperature', 'dewpoint', 'u_wind', 'v_wind'] # 'latitude', 'longitude']
        units_list = ['hPa', 'meter', 'degC', 'degC', 'kt', 'kt']
        clean_data = {}
        non_dups = np.concatenate(([True], np.diff(df.to_dict('list')['pressure']) != 0))
        for old_key, new_key, unit in zip (old_keys, new_keys, units_list):
            clean_data[new_key] = np.array(df.to_dict('list')[old_key])[non_dups]*units(unit)
        clean_data['station_info'] = [RAOB_STATIONS[RAOB_STATIONS['ICAO']==station]['ICAO'].values[0],
                                      RAOB_STATIONS[RAOB_STATIONS['ICAO']==station]['NAME'].values[0],
                                      RAOB_STATIONS[RAOB_STATIONS['ICAO']==station]['LOC'].values[0],
                                      raob_latlon(station),
                                      RAOB_STATIONS[RAOB_STATIONS['ICAO']==station]['EL(m)'].values[0]]
        try:
            slc = (len(clean_data['p']) - np.where(clean_data['p']<=98.*units('hPa'))[0][0])
            for key in new_keys:
                clean_data[key] = clean_data[key][:-slc]
        except:
            pass
        elapsed_time = time.time() - st
        print('RUNTIME:', time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
        return clean_data    

    
    
    
    
    
    
    
    
    
    

    
##################################################################################################################################################
#                                                        PARSE_DATA() FUNCTION                                                                   #
##################################################################################################################################################

    
def parse_data(raw_data):
    st = time.time()
    '''
    parse_data() [function, used to parse sounding data after using a get_data() function]

    PARAMETERS
    ----------
     name          dtype      required   explanation 
     raw_data:     dataset,   required,  xr.DataSet containing ERA5 renalaysis data
    
    RETURNS
    -------
     name         dtype    explanation 
     clean_data   dict,    cleaned profile data with units in a dict format. 
    
    COPYRIGHT
    ---------
    Created by Kyle J Gillett (@wxkylegillett) 2023

    '''
    
    if str(type(raw_data)) == "<class 'xarray.core.dataset.Dataset'>":
        
        print(f'-- ERA5 REANALYSIS DATA PARSE FUNCTION --\n------------------------------------------')
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
        'vert_w': (np.array(raw_data['w'][0])),
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


    
    
    
    if str(type(raw_data)) == "<class 'netCDF4._netCDF4.Dataset'>":

        print(f'-- RAP REANALYSIS DATA PARSE FUNCTION --\n------------------------------------------')

    
        vert_data = {
            'vert_T' : (ma.getdata(raw_data.variables['Temperature_isobaric'][0,:,:,:]-273.15)),
            'vert_p' : (np.array([100, 125, 150, 175, 200, 225, 250, 275, 300, 325, 350, 375, 400, 425, 450, 475, 500, 525, 550, 575, 600, 625, 650, 675, 700, 725, 750, 775, 800, 825, 850, 875, 900, 925, 950, 975, 1000])),
            'vert_z' : (ma.getdata(raw_data.variables['Geopotential_height_isobaric'][0,:,:,:])),
            'vert_rh': (ma.getdata(raw_data.variables['Relative_humidity_isobaric'][0,:,:,:])),
            'vert_u' : (ma.getdata(raw_data.variables['u-component_of_wind_isobaric'][0,:,:,:]*1.94384)),
            'vert_v' : (ma.getdata(raw_data.variables['v-component_of_wind_isobaric'][0,:,:,:]*1.94384)),
            'vert_w' : (ma.getdata(raw_data.variables['Vertical_velocity_pressure_isobaric'][0,:,:,:])),
            'vert_Td': (mpcalc.dewpoint_from_relative_humidity(
                (ma.getdata(raw_data.variables['Temperature_isobaric'][0,:,:,:]-273.15))*units.degC, 
                (ma.getdata(raw_data.variables['Relative_humidity_isobaric'][0,:,:,:])*units.percent))).m,
        } 

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

        latlon_data = {
        'data_lat'    : (ma.getdata(1000*(raw_data.variables[raw_data.variables['Pressure_surface'].dimensions[1]][:]))),
        'data_lon'    : (ma.getdata(1000*(raw_data.variables[raw_data.variables['Pressure_surface'].dimensions[2]][:]))),
        'data_latnum' : (ma.getdata(1000*(raw_data.variables[raw_data.variables['Pressure_surface'].dimensions[1]][:]))).shape[0],
        'data_lonnum' : (ma.getdata(1000*(raw_data.variables[raw_data.variables['Pressure_surface'].dimensions[2]][:]))).shape[0]
        }
    
    
    
        
        
    # NOW THAT DATA IS PARSED AND IN THE SAME FORMAT, CLEAN IT UP:    
    
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
    
    crs = ccrs.PlateCarree()
    clons, clats = np.meshgrid(latlon_data['data_lon'], latlon_data['data_lat'])
    tlons = clons 
    tlats = clats
    def haversine(lon_arr, lat_arr, lon_point, lat_point): 
        distance = 6371*(2*np.arcsin(np.sqrt(np.sin((np.radians(lat_arr)-np.radians(lat_point))/2)**2 + np.cos(np.radians(lat_point)) * np.cos(np.radians(lat_arr)) * np.sin((np.radians(lon_arr)-np.radians(lon_point))/2)**2)))
        return distance
    location_index = np.where(haversine(tlons,tlats,latlon1[1],latlon1[0]) == np.min(haversine(tlons,tlats,latlon1[1],latlon1[0])))
    
    
    # surface data 
    sfc_loc_idx  = {}
    vert_loc_idx = {}
    
    for key in sfc_data.keys():
        sfc_loc_idx[key] = sfc_data[key][location_index[0],location_index[1]][0]
    for key in ['vert_T', 'vert_z', 'vert_rh', 'vert_u', 'vert_v', 'vert_w', 'vert_Td']:
        vert_loc_idx[key] = vert_data[key][::-1,location_index[0],location_index[1]]
        vert_loc_idx['vert_p'] = vert_data['vert_p'][::-1]

    # Creates Arrays for Vertical Profile
    sfc_index  = get_sfc_index((vert_loc_idx['vert_z']-sfc_loc_idx['sfc_z']))
    surface_height = sfc_loc_idx['sfc_z']
    sfc_loc_idx['sfc_z'] = 0
    sb_dict = {}
    new_keys = ['T', 'Td', 'rh', 'u', 'v', 'z', 'p']
    sfc_keys = ['sfc_T', 'sfc_Td', 'sfc_rh', 'sfc_u', 'sfc_v', 'sfc_z', 'sfc_p']
    vert_keys = ['vert_T', 'vert_Td', 'vert_rh', 'vert_u', 'vert_v', 'vert_z', 'vert_p']
    
    for vert_key, sfc_key, new_key in zip(vert_keys, sfc_keys, new_keys):
        sb_dict[new_key] = make_sfc_based(vert_loc_idx[vert_key],  sfc_loc_idx[sfc_key],  sfc_index)
    sb_dict['w'] = make_sfc_based(vert_loc_idx['vert_w'],  vert_loc_idx['vert_w'][0],  sfc_index)
    # if sb_dict['p'][0] < sb_dict['p'][1]:
    #     sb_dict['p'][0] = sb_dict['p'][0] + 25
    
    # Interpolates data
    dz = 250 
    soundingtop_hght = sb_dict['z'][-1]
    toplvl      = int(soundingtop_hght/dz)*dz
    numlvls     = int(toplvl/dz)
    interp_lvls = np.linspace(0,toplvl,numlvls+1)
    
    keys = ['T', 'Td', 'rh', 'u', 'v', 'z', 'p', 'w']
    units_list = ['degC', 'degC', 'percent', 'kt', 'kt', 'm', 'hPa']
    interp_dict = {}
    zeros_dict  = {}
    final_profile  = {}
    
    for key in keys:
        interp_dict[key] = (interpolate.interp1d(sb_dict['z'], sb_dict[key]))
        zeros_dict[key]  = np.zeros((len(interp_lvls)))
        for zeros_arr in zeros_dict.values():
            zeros_dict[key][0] = sb_dict[key][0]
        for i in range(1,len(zeros_dict[key]),1):
            zeros_dict[key][i] = interp_dict[key](dz*i)
    for i, unit, key in zip(range(0, len(units_list)), units_list, keys):
        final_profile[key] = zeros_dict[key]*units(unit)
    final_profile['w'] = zeros_dict['w']/units('sec') 
    final_profile['zAGL'] = final_profile['z'] + surface_height*units.m
    clean_data = final_profile 
    
    print('- COMPLETE -')
    elapsed_time = time.time() - st
    print('RUNTIME:', time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
    return clean_data







##################################################################################################################################################
#                                                       METPY_SOUNDING() FUNCTION                                                                #
##################################################################################################################################################

def metpy_sounding(clean_data):
    st = time.time()
    import matplotlib.pyplot as plt
    import metpy.calc as mpcalc
    from metpy.plots import add_metpy_logo, SkewT, Hodograph
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    from metpy.units import units
    
    
    '''
        metpy_sounding(clean_data)

        this function inputs cleaned profile data through a slightly modified MetPy sounding plot script. 

        PARAMETERS
        ----------
        name         dtype    required   example      explanation 
        clean_data:  dict,      yes,    clean_data    a dict of cleaned profile data acquired from the above functions. 

        RETURNS
        -------
        img: MetPy Sounding Plot

        EXAMPLES
        --------
        spy.metpy_sounding(clean_sounding)

        COPYRIGHT
        ---------
        Created by Kyle J Gillett (@wxkylegillett) 2023

    '''
    
    fig = plt.figure(figsize=(16, 12), edgecolor="#04253a")
    add_metpy_logo(fig, 80, 80)
    skew = SkewT(fig, rotation=45)
    fig.set_facecolor('#ffffff')        
    skew.ax.set_facecolor('#ffffff')                                                    # set facecolor to white                         
    skew.plot(clean_data['p'], clean_data['T'], 'r', lw=3)                              # plot T
    skew.plot(clean_data['p'], clean_data['Td'], 'g', lw=3)                             # plot Td
    interval = np.logspace(2.113, 3, 30) * units.hPa                                    # Arrange wind barbs for best fit
    idx = mpcalc.resample_nn_1d(clean_data['p'], interval)                              # Resample wind barbs for bes
    skew.plot_barbs(clean_data['p'][idx], clean_data['u'][idx], clean_data['v'][idx])   # plot wind barbs         
    skew.ax.set_ylim(1070, 95)                                                          # set y limits
    skew.ax.set_xlim(-30, 70)                                                           # set x limits
    plt.xticks(fontsize=13)                                                             # define x axis tick mark size
    plt.yticks(fontsize=13, ha='left')                                                  # define y axis tick mark size
    plt.tick_params(axis="x",direction="in", pad=-12)                                   # move x ticks in
    plt.tick_params(axis="y",direction="in", pad=-7)                                    # move y ticks in
    plt.xlabel("  ", fontsize=12)                                                       # remove x axis label
    plt.ylabel("  ", fontsize=12)                                                       # remove y axis label
    # ADD A SIMPLE PLOT TITLE
    if source == 'obs':
        plt.title(f"{clean_data['station_info'][0]} |  VALID: {month1}/{day1}/{year1} {hour1}Z  |  {clean_data['station_info'][3][0]}, {clean_data['station_info'][3][1]}",loc='left', weight='bold', fontsize=15)
    else:
        plt.title(f'{source} F00  |  VALID: {month1}/{day1}/{year1} {hour1}Z  |  {latlon1}',loc='left', weight='bold', fontsize=15)
    plt.title("POWERED BY METPY & SOUNDERPY", fontsize=12, alpha=0.5, weight='bold', color='blue', loc='right')
    # Set some better labels than the default
    skew.ax.set_xlabel(f"Temperature ({clean_data['T'].units:~P})")
    skew.ax.set_ylabel(f"Pressure ({clean_data['p'].units:~P})")
    # Calculate LCL height and plot as black dot. Because `p`'s first value is
    # ~1000 mb and its last value is ~250 mb, the `0` index is selected for
    # `p`, `T`, and `Td` to lift the parcel from the surface. If `p` was inverted,
    # i.e. start from low value, 250 mb, to a high value, 1000 mb, the `-1` index
    # should be selected.
    lcl_pressure, lcl_temperature = mpcalc.lcl(clean_data['p'][0], clean_data['T'][0], clean_data['Td'][0])
    skew.plot(lcl_pressure, lcl_temperature, 'ko', markerfacecolor='black')
    # Calculate full parcel profile and add to plot as black line
    prof = mpcalc.parcel_profile(clean_data['p'], clean_data['T'][0], clean_data['Td'][0]).to('degC')
    skew.plot(clean_data['p'], prof, 'k', linestyle='--', linewidth=2)
    try:
        # Shade areas of CAPE and CIN
        skew.shade_cin(clean_data['p'], clean_data['T'], prof, clean_data['Td'])
        skew.shade_cape(clean_data['p'], clean_data['T'], prof)
    except:
        pass
    #highlight 0degC
    skew.ax.axvline(0, color='c', linestyle='--', linewidth=2)
    # Add the relevant special lines
    skew.plot_dry_adiabats(color='black', linewidth=0.5, alpha=0.3)            
    skew.plot_moist_adiabats(linewidth=0.5, alpha=0.3)         
    skew.plot_mixing_lines(linewidth=0.2, alpha=0.3)      
    x1 = np.linspace(-100, 40, 8)                                                          # The starting x values for the shaded regions
    x2 = np.linspace(-90, 50, 8)                                                           # The ending x values for the shaded regions
    y = [1100, 50]                                                                        # The range of y values that the shades regions should cover
    for i in range(0, 8):                                                                  # shade every 10 degC isotherms
        skew.shade_area(y=y, x1=x1[i], x2=x2[i], color='gray', alpha=0.02, zorder=1)       
    # Create a hodograph
    # Create an inset axes object that is 40% width and height of the
    # figure and put it in the upper right hand corner.
    # determine max height of wind data to plot on hodograph in km (if hodo_layer = 9, 0-9km u and v are plotted)
    hodo_layer = int(9)
    # remove nan values from base wind u and v component arrays to find min & max values.
    u_clean = clean_data['u'].magnitude[np.logical_not(np.isnan(clean_data['u'].magnitude))]
    v_clean = clean_data['v'].magnitude[np.logical_not(np.isnan(clean_data['v'].magnitude))]
    # restructure u and v, p, ws and z data arrays based on corrected u and v arrays and hodo_layer depth
    p_hodo, u_hodo, v_hodo, z_hodo = mpcalc.get_layer(clean_data['p'], clean_data['u'], clean_data['v'], clean_data['z'], depth=hodo_layer*units.km)
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
    
    ax_hod = inset_axes(skew.ax, '40%', '40%', loc=1)
    h = Hodograph(ax_hod, component_range=80.)
    try:
        h.ax.set_xlim(x_Minlimit, x_Maxlimit)                                                 # set x axis bounds
        h.ax.set_ylim(y_Minlimit, y_Maxlimit)                                                 # set y axis bounds    
    except:
        pass
    h.add_grid(increment=20, linestyle='-', linewidth=1.5, alpha=0.2)                     # define 1st hodo grid
    h.add_grid(increment=10, color='black', linewidth=1, linestyle='--', alpha=0.2)       # define 2nd hodo grid
    h.ax.set_facecolor('#ffffff')                                                         # hodo background color
    h.ax.set_box_aspect(1)                                                                # set hodo aspect ratio
    h.ax.set_yticklabels([])
    h.ax.set_xticklabels([])
    h.ax.set_xticks([])
    h.ax.set_yticks([])
    h.ax.set_xlabel(' ')
    h.ax.set_ylabel(' ')
    hodo_plot = h.plot_colormapped(u_hodo, v_hodo, z_hodo, lw=4)  # Plot a line colored by wind speed
    plt.xticks(np.arange(0,0,1))
    plt.yticks(np.arange(0,0,1))
    for i in range(10,120,10):
        h.ax.annotate(str(i),(i,0),xytext=(0,2),textcoords='offset pixels',clip_on=True,fontsize=10,weight='bold',alpha=0.3,zorder=0)
    for i in range(10,120,10):
        h.ax.annotate(str(i),(0,i),xytext=(0,2),textcoords='offset pixels',clip_on=True,fontsize=10,weight='bold',alpha=0.3,zorder=0)
    # Show the plot
    plt.show()
    elapsed_time = time.time() - st
    print('RUNTIME:', time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
    
    
    
    
    
    
    

##################################################################################################################################################
#                                                           GET_DOCS() FUNCTION                                                                  #
##################################################################################################################################################
 
def get_docs():
    docs_text = f"""
## ------------ VERTICAL PROFILE DATA RETRIEVAL TOOL DOCUMENTATION ----------------
## check out the full documentation here: https://github.com/kylejgillett/sounderpy/blob/main/DOCUMENTATION.md
## 
## SounderPy can access and parse data from:
## > ECMWF CDS ERA5 reanalysis [1940-present] *note: you must set up an account through the CDS to unlock ERA5 data. 
##    (see: https://cds.climate.copernicus.eu/api-how-to)
## > UNIDATA THREDDS TDS RAP reanalysis [2005-present]
## > UNIDATA THREDDS TDS RUC reanalysis [2005-2020]
## > The University of Wyoming RAOB archive [1973-present, depending on station]
## > Iowa State University's RAOB archive [1945-present, depending on station]
##
##
## COPYRIGHT
## ---------
## Created by Kyle J Gillett (@wxkylegillett) 2023
##
## Version: {version} | {version_date}
## --------------------------------------------------------------------------------
##
## TOOLS INCLUDED IN THIS PACKAGE:
## -------------------------------
## get_docs(): [function, prints quick documentation for this package]
##
## get_model_data(method, domain, latlons, year, month, day, hour): [function, used to access model sounding data] - returns raw_data 
##
## get_obs_data(station, year, month, day, hour): [function, used to access RAOB sounding data] - returns clean_data 
##
## parse_data(raw_data) [function, used to parse sounding data after using a get_model_data() function] - returns clean_data
##
## metpy_sounding(clean_data)  [function, used to metpy's skewt plotting functions to create a quick sounding 
##                                of the parsed data the original script from MetPy is found here: 
##                                https://unidata.github.io/MetPy/latest/tutorials/upperair_soundings.html
## metar_latlon(metar_site) [function, get lat-lon pair for a metar site] - returns list of lat/lon pair
##
## raob_latlon(raob_site) [function, get lat-lon pair for a raob site] - returns list of lat/lon pair
##
##
## EXAMPLES:
##     >>> import sounderpy as spy
##     >>> raw_data   = spy.get_model_data('rap', [43.45, -84.57], '2014', '06', '16', '21')   
##     >>> clean_data = spy.parse_data(raw_data)
##     >>> spy.metpy_sounding(clean_data)
"""
    print(docs_text)

    

    
    

    
    
##################################################################################################################################################
#                                                      SOUNDERPY HELPER FUCNTIONS                                                                #
##################################################################################################################################################
'''
     
     A collection of helper fucntions that users may find useful for processing vertical profile data but are not nessacary to use the basic 
     functions of SounderPy.
     
'''
    
####################################################### DATA INTERPOLATION ################################################################# 
def interp_data(variable, heights, step=100):
    
    '''
        interpolate vertical profile over vertical height array
        
        PARAMETERS
        ----------
        name        dtype    required    explanation 
        variable:   array    yes         the data to be interpolated. Must be same length as height array
        heights:    array    yes         heights corresponding to the vertical profile used to interpolate
        step:       int      no          resolution of interpolation. Default is 100 (recommended value is 100)
        
        RETURNS
        -------
        varinterp:  array,  array of interpolated data 
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
        
        PARAMETERS
        ----------
        name         dtype                        required    explanation 
        array        array of ints or floats      yes         the array of data search through to 'find the nearest'
        value        int or float (same as array) yes         the value used to campare against the array of data
        
        RETURNS
        -------
        nearest_idx:  int,  index of the data array that corresponds with the nearest value to the given value. 
    '''
    
    array = np.asarray(array)
    nearest_idx = (np.abs(array - value)).argmin()
    return nearest_idx

    

    
##################################################### MAKE SURFACED BASED #################################################################### 

def get_sfc_index(height_arr):
    
    """
        This function takes an array of height values and finds the index before the first positive value

        PARAMETERS
        ----------
        name         dtype         required     explanation 
        height_arr   array         yes          array of height values
        RETURNS
        -------
        index:       int     indicating the index just before the first positive height value. If -1 is returned,
               a positive index could not be found.
    
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

        PARAMETERS
        ----------
        name       dtype         required     explanation 
        arr        array         yes          the array to manipulate
        sfc_val    int or float  yes          the value that corresponds to the true surface value
        sfc_index  int           yes          index indicating where the true sfc is (where copying occurs) (found via get_sfc_index())

        RETURNS
        -------
        mod_arr:   array   array of values from arr[index:len(arr)]
    
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


        PARAMETERS
        ----------
        name       dtype         required     explanation 
        arr        3D array      yes          the 3D vertical profile(s) array to manipulate
        sfc_arr    2D array      yes          the true surface values to be added to the vertical profile(s)

        RETURNS
        -------
        mod_arr:   array   array of values from arr[index:len(arr)]
    
    '''
    mod_arr = np.zeros((np.shape(arr)[0]+1,np.shape(arr)[1],np.shape(arr)[2]))
    for j in range(np.shape(arr)[1]):
        for k in range(np.shape(arr)[2]):
            mod_arr[0,j,k] = sfc_arr[j,k]
            mod_arr[1:,j,k] = arr[:,j,k]
    return mod_arr





#################################################### FIND METAR SITE LATLON ################################################################## 

def metar_latlon(metar_site):
    
    '''
        takes a METAR site id such as 'KMBS', searches over 9000 station IDs  and returns
        a list including the lat/lon for the METAR site

        PARAMETERS
        ----------
        name        dtype    required   example   explanation 
        metar_site  str      yes        'KMBS'    4 digit METAR site id

        RETURNS
        -------
        latlon:     list      lat lon pair for requested METAR site in a list
    '''
    
    def dms2dd(degrees, minutes, direction):
        dd = float(degrees) + float(minutes) / 60
        if direction == "S" or direction == "W":
            dd *= -1
        return dd
    
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
            lat = dms2dd(data[39:41], data[42:44], data[44:45])
            lon = dms2dd(data[47:50], data[51:53], data[53:54])
            altitude = int(data[55:59])
            country = data[81:83]

            if icao:
                stations[icao] = {"name": station, "lat": lat, "lon": lon, "altitude": altitude, "country": country}
    try:
        latlon = [np.round(stations.get(metar_site)['lat'],2), np.round(stations.get(metar_site)['lon'],2)]
        return latlon
    except:
        print('[!] ERROR [!]\nMETAR site lat/lon not available')
        pass
    
    
    
    
    
#################################################### FIND RAOB SITE LATLON ##################################################################    
    
def raob_latlon(raob_site):

    def dms2dd(degrees, direction):
        dd = float(degrees)
        if direction == "S" or direction == "W":
            dd *= -1
        return dd
    RAOB_STATIONS = pd.read_csv(f'https://raw.githubusercontent.com/kylejgillett/sounderpy/main/RAOB-STATIONS.txt', skiprows=7, skipinitialspace = True)
    try:
        station = RAOB_STATIONS['ICAO'][np.where(RAOB_STATIONS['ICAO'].str.contains(raob_site, na=False, case=True))[0]].values[0]

        lat = dms2dd(RAOB_STATIONS[RAOB_STATIONS['ICAO']==station]['LAT'].values[0], RAOB_STATIONS[RAOB_STATIONS['ICAO']==station]['A'].values[0])
        lon = dms2dd(RAOB_STATIONS[RAOB_STATIONS['ICAO']==station]['LON'].values[0], RAOB_STATIONS[RAOB_STATIONS['ICAO']==station]['B'].values[0])
        return [lat, lon]
    
    except:
        print("[!] ERROR [!]\n[!] NO RAOB SITE FOUND FOR THAT ID[!]\n > Make sure you entered a valid RAOB ID\n > Make sure you entered the RAOB ID as a string, ex: 'DTX'")
        pass

from datetime import datetime
import warnings
import time

import xarray as xr
import numpy as np
import numpy.ma as ma

from urllib.request import urlopen
import cdsapi
from siphon.ncss import NCSS
import netCDF4
from siphon.catalog import TDSCatalog

from scipy import interpolate
import metpy.calc as mpcalc
from metpy.units import units

from .calc import *





"""
    SOUNDERPY MODEL REANALYSIS 'GET-DATA' FUNCTIONS  

    Purpose of module: 

    House function for loading and parsing ERA5, RAP, RUC and NCEP-FNL
    model reanalysis data. Functions here are referenced by sounderpy.py


    (C) KYLE J GILLETT, UNIVERSITY OF NORTH DAKOTA, 2024
"""






#######################
# MODEL REANALYSIS DATA
#########################################################################

def fetch_model(model, latlon, year, month, day, hour, dataset, box_avg_size, hush, clean_it):

    # record process time
    st = time.time()

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


    # send error message if given model is invalid
    if model.casefold() not in ['era', 'era5', 'rap', 'ruc', 'rap-ruc', 'rap-now', 'ncep-fnl', 'ncep']:
        raise ValueError(
            f"The model you requested, '{model}', is not a valid model. Valid models for this function include ['rap-ruc', 'era5', 'ncep']")

    # create list of lat-lon points for box-average domain
    #                    + lat                   # - lat                   # - lon                   # + lon
    latlons = [latlon[0] + box_avg_size, latlon[0] - box_avg_size, latlon[1] - box_avg_size, latlon[1] + box_avg_size]






    ### ERA 5 REANALYSIS ###
    #########################################################################################################
    '''
    Get ERA-5 reanalysis data the CDS API, return an xarray dataset of the data.
    '''
    if model.casefold() in ['era', 'era5']:
        # define source
        source = 'ERA5'
        dtype = 'reanalysis'

        print(f'> ERA5 REANALYSIS DATA ACCESS FUNCTION\n  ------------------------------------------')
        print(f'    > some messages from ECMWF CDS...')

        # define ERA5 dataset names we want to access data from
        dataset_presLvls = 'reanalysis-era5-pressure-levels'
        dataset_singleLvls = 'reanalysis-era5-single-levels'

        # rearrange lat-lon list for my sanity
        latlon_list = [latlons[0], latlons[2], latlons[1], latlons[3]]

        # set up cds api call for pressure level data
        c = cdsapi.Client()
        params = {
            'product_type': 'reanalysis',
            'variable': ['temperature', 'geopotential', 'relative humidity', 'U WIND COMPONENT', 'V WIND COMPONENT',
                         'vertical_velocity'],
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
            'year': year,
            'month': month,
            'day': day,
            'time': f'{hour}:00',
            'format': 'netcdf',
            'area': latlon_list
        }

        # set up cds api call for surface data
        c2 = cdsapi.Client()
        params2 = {
            'product_type': 'reanalysis',
            'variable': ['2m_temperature', '2m_dewpoint_temperature', 'surface_pressure', '10u', '10v', 'z', 'msl'],
            'year': year,
            'month': month,
            'day': day,
            'time': f'{hour}:00',
            'format': 'netcdf',
            'area': latlon_list
        }

        # retrieve data from CDS
        fl = c.retrieve(dataset_presLvls, params)
        print('    > DATASET ACCESSED: ' + dataset_presLvls)
        fl2 = c2.retrieve(dataset_singleLvls, params2)
        print('    > DATASET ACCESSED: ' + dataset_singleLvls)
        # load data to memory via -output.nc files
        fl.download("./presLvls-output.nc")
        fl2.download("./singleLvls-output.nc")

        # create xarray datasets from .nc files
        ds = xr.open_dataset("./presLvls-output.nc")
        ds2 = xr.open_dataset("./singleLvls-output.nc")

        # merge the two datasets together
        ds2 = ds2.rename({'z': 'hgts', 'sp': 'ps', 't2m': 'Ts', 'd2m': 'tds', 'u10': 'us', 'v10': 'vs', })
        raw_data = xr.merge([ds, ds2])

        ds.close()
        ds2.close()

        #########################################################################################################








    ### RAP REANALYSIS ###
    #########################################################################################################
    '''
    Get RAP reanalysis data from NCEI THREDDS Server, return a netcdf4 dataset
    '''
    if model in ['rap', 'ruc', 'rap-ruc']:

        dtype = 'reanalysis'

        print(f'> RAP REANALYSIS DATA ACCESS FUNCTION\n  -----------------------------------------')

        # rearrange latlon list for my sanity
        latlon_list = [latlons[2], latlons[3], latlons[1], latlons[0]]

        # create dict of RAP-data urls for the different datasets of RAP and RUC from NCEI
        urls = {
            'RAP_25km': 'https://www.ncei.noaa.gov/thredds/ncss/model-rap252/' + str(year) + str(month) + '/' + str(
                year) + str(month) + str(day) + '/rap_252_' + str(year) + str(month) + str(day) + '_' + str(
                hour) + '00_000.grb2',
            'RAP_25km_old': 'https://www.ncei.noaa.gov/thredds/ncss/model-rap252-old/' + str(year) + str(
                month) + '/' + str(year) + str(month) + str(day) + '/rap_252_' + str(year) + str(month) + str(
                day) + '_' + str(hour) + '00_000.grb2',

            'RAP_25km_anl': 'https://www.ncei.noaa.gov/thredds/ncss/model-rap252anl/' + str(year) + str(
                month) + '/' + str(year) + str(month) + str(day) + '/rap_252_' + str(year) + str(month) + str(
                day) + '_' + str(hour) + '00_000.grb2',
            'RAP_25km_anl_old': 'https://www.ncei.noaa.gov/thredds/ncss/model-rap252anl-old/' + str(year) + str(
                month) + '/' + str(year) + str(month) + str(day) + '/rap_252_' + str(year) + str(month) + str(
                day) + '_' + str(hour) + '00_000.grb2',

            'RAP_13km': 'https://www.ncei.noaa.gov/thredds/ncss/model-rap130/' + str(year) + str(month) + '/' + str(
                year) + str(month) + str(day) + '/rap_130_' + str(year) + str(month) + str(day) + '_' + str(
                hour) + '00_000.grb2',
            'RAP_13km_old': 'https://www.ncdc.noaa.gov/thredds/ncss/model-rap130-old/' + str(year) + str(
                month) + '/' + str(year) + str(month) + str(day) + '/rap_130_' + str(year) + str(month) + str(
                day) + '_' + str(hour) + '00_000.grb2',

            'RAP_13km_anl': 'https://www.ncei.noaa.gov/thredds/ncss/model-rap130anl/' + str(year) + str(
                month) + '/' + str(year) + str(month) + str(day) + '/rap_130_' + str(year) + str(month) + str(
                day) + '_' + str(hour) + '00_000.grb2',
            'RAP_13km_anl_old': 'https://www.ncdc.noaa.gov/thredds/ncss/model-rap130anl-old/' + str(year) + str(
                month) + '/' + str(year) + str(month) + str(day) + '/rap_130_' + str(year) + str(month) + str(
                day) + '_' + str(hour) + '00_000.grb2',

            'RUC_13km': 'https://www.ncei.noaa.gov/thredds/ncss/model-ruc130anl/' + str(year) + str(month) + '/' + str(
                year) + str(month) + str(day) + '/ruc2anl_130_' + str(year) + str(month) + str(day) + '_' + str(
                hour) + '00_000.grb2',
            'RUC_13km_old': 'https://www.ncei.noaa.gov/thredds/ncss/model-ruc130anl-old/' + str(year) + str(
                month) + '/' + str(year) + str(month) + str(day) + '/ruc2anl_130_' + str(year) + str(month) + str(
                day) + '_' + str(hour) + '00_000.grb2',

            'RUC_25km': 'https://www.ncei.noaa.gov/thredds/ncss/model-ruc252anl/' + str(year) + str(month) + '/' + str(
                year) + str(month) + str(day) + '/ruc2anl_252_' + str(year) + str(month) + str(day) + '_' + str(
                hour) + '00_000.grb',
            'RUC_25km_old': 'https://www.ncei.noaa.gov/thredds/ncss/model-ruc252anl-old/' + str(year) + str(
                month) + '/' + str(year) + str(month) + str(day) + '/ruc2anl_252_' + str(year) + str(month) + str(
                day) + '_' + str(hour) + '00_000.grb'
        }

        # if a user defined a target dataset, try it.
        tries = 0
        if dataset is not None:
            try:
                print(f'    > SEARCHING FOR {dataset}...')
                url = NCSS(urls[dataset])
                print(f'    > DATASET FOUND: {dataset}')
                url_to_use = urls[dataset]
                source = str(dataset)[0:3]
                data = NCSS(url_to_use)
            except:
                raise ValueError(
                    f'NCSS Connection failed -- the date and dataset you requested is invalid. Ensure you have the correct dates, model, & dataset name\n' +
                    f'Note: data may not be available for every date/time. This catalog experiences periodic outages and may host missing data.\n' +
                    f'The date you entered is: {year}-{month}-{day}-{hour}z. An example of a valid date is: 2014-06-16-18z')

        # else, create a simple test for each URL, use the first one that works
        # and return its data
        else:

            tries = 0
            for url, key in zip(urls.values(), urls.keys()):
                try:
                    NCSS(url)
                    print(f'    > DATASET USED: {key}')
                    url_to_use = url
                    source = str(key)[0:3]
                    break
                except:
                    tries += 1
                    pass

            if tries == 12:
                raise ValueError(
                    f'NCSS Connection failed -- ensure you have the correct dates and corresponding model\n' +
                    f'Note: data may not be available for every date/time. This catalog experiences periodic outages and may host missing data.\n' +
                    f'The date you entered is: {year}-{month}-{day}-{hour}z. An example of a valid date is: 2014-06-16-18z')
            else:
                data = NCSS(url_to_use)

        # set up TDS query
        query = data.query()

        # subset data by variable names for RAP & RUC (of course they have to be different)
        if source in ['rap', 'RAP']:
            query.variables('Pressure_surface',
                            'Geopotential_height_isobaric', 'Geopotential_height_surface',
                            'Temperature_isobaric', 'Temperature_height_above_ground',
                            'Relative_humidity_isobaric', 'Dewpoint_temperature_height_above_ground',
                            'Relative_humidity_height_above_ground', 'Vertical_velocity_pressure_isobaric',
                            'u-component_of_wind_height_above_ground', 'v-component_of_wind_height_above_ground',
                            'u-component_of_wind_isobaric', 'v-component_of_wind_isobaric').add_lonlat()
        else:
            query.variables('Pressure_surface',
                            'Geopotential_height_isobaric', 'Geopotential_height_surface',
                            'Temperature_isobaric', 'Temperature_height_above_ground',
                            'Relative_humidity_isobaric', 'Dewpoint_temperature_height_above_ground',
                            'Relative_humidity_height_above_ground', 'Vertical_velocity_pressure_isobaric',
                            'u-component_of_wind_height_above_ground', 'v-component_of_wind_height_above_ground',
                            'u-component_of_wind_isobaric', 'v-component_of_wind_isobaric').add_lonlat()

        # subset data by requested domain
        query.lonlat_box(latlon_list[0], latlon_list[1], latlon_list[2], latlon_list[3])

        # load the data from TDS
        raw_data = data.get_data(query)

    #########################################################################################################








    ### NCEP FNL REANALYSIS ###
    #########################################################################################################
    '''
    Get NCEP-FNL reanalysis data from NCEI THREDDS server, return a netcdf4 dataset
    '''

    if model.casefold() in ['ncep-fnl', 'ncep', 'fnl']:
        # define source
        source = 'NCEP-FNL'
        dtype = 'reanalysis'

        print(f'> NCEP-FNL REANALYSIS DATA ACCESS FUNCTION\n  ------------------------------------------')

        latlon_list = [latlons[0], latlons[2], latlons[1], latlons[3]]

        # access ncss thredds server, if access failed, exit function and return an error
        try:
            data = NCSS(
                f"https://thredds.rda.ucar.edu/thredds/ncss/grid/files/g/ds083.3/{year}/{year}{month}/gdas1.fnl0p25.{year}{month}{day}{hour}.f00.grib2")
            worked = True
        except:
            worked = False
            pass

        if worked:
            # set up TDS query
            query = data.query()

            query.variables('Geopotential_height_isobaric', 'Geopotential_height_surface',
                            'Temperature_isobaric', 'Temperature_height_above_ground',
                            'Relative_humidity_isobaric', 'Dewpoint_temperature_height_above_ground',
                            'Relative_humidity_height_above_ground', 'Pressure_surface', 'Vertical_velocity_pressure_isobaric',
                            'u-component_of_wind_height_above_ground', 'v-component_of_wind_height_above_ground',
                            'u-component_of_wind_isobaric', 'v-component_of_wind_isobaric').add_lonlat()

            # subset data by requested domain
            # north=90.000&    west=-.125&    east=-.125&    south=-90.000
            query.lonlat_box(latlon_list[1], latlon_list[3], latlon_list[2], latlon_list[0])

            # load the data from TDS
            raw_data = data.get_data(query)

        else:
            raise ValueError(f'NCSS Connection failed -- ensure you have the correct dates and corresponding model\n' +
                             f'Note: data may not be available for every date/time. This catalog experiences periodic outages and may host missing data.\n' +
                             f'The date you entered is: {year}-{month}-{day}-{hour}z. An example of a valid date is: 2014-06-16-18z')
    #########################################################################################################







    ### RAP ANALYSIS ###
    #########################################################################################################
    '''
    Get latest RAP analysis from NCEI THREDDS Server 
    '''
    if model in ['rap-now']:

        print(f'> RAP ANALYSIS DATA ACCESS FUNCTION\n  -----------------------------------------')

        latlon_list = [latlons[2], latlons[3], latlons[1], latlons[0]]

        # define dataset URL & try to access it to make sure it works
        url = 'https://thredds.ucar.edu/thredds/catalog/grib/NCEP/RAP/CONUS_13km/latest.xml'
        try:
            cat = TDSCatalog('https://thredds.ucar.edu/thredds/catalog/grib/NCEP/RAP/CONUS_13km/latest.xml')
            source = 'RAP'
            dtype = 'analysis'
            worked = True
        except:
            worked = False
            pass
        if worked:
            # set up TDS query
            latest_ds = list(cat.datasets.values())[0]
            ncss = NCSS(latest_ds.access_urls['NetcdfSubset'])
            query = ncss.query()
            # Find start time
            start_time = ncss.metadata.time_span['begin']
            fcst_date = datetime.strptime(start_time, '%Y-%m-%dT%H:%M:%SZ')
            year1 = fcst_date.strftime('%Y')
            month1 = fcst_date.strftime('%m')
            day1 = fcst_date.strftime('%d')
            hour1 = fcst_date.strftime('%H')
            # Subset data by time
            query.time(fcst_date).accept('netcdf4')
            # Subsets data by variables
            query.variables('MSLP_MAPS_System_Reduction_msl', 'Pressure_surface', 'Geopotential_height_isobaric',
                            'Temperature_isobaric', 'Relative_humidity_isobaric', 'Temperature_height_above_ground',
                            'Relative_humidity_height_above_ground', 'u-component_of_wind_height_above_ground',
                            'v-component_of_wind_height_above_ground', 'u-component_of_wind_isobaric',
                            'Vertical_velocity_pressure_isobaric', 'v-component_of_wind_isobaric').add_lonlat()

            # Subset data by lat-lon domain
            query.lonlat_box(latlon_list[0], latlon_list[1], latlon_list[2], latlon_list[3])

            # Gets data
            raw_data = ncss.get_data(query)

        else:
            raise ValueError(f'NCSS Connection failed -- RAP data may not be available at this time')
    #########################################################################################################










    #########################################################################################################
    ### PARSE RAW MODEL DATA ################################################################################
    #########################################################################################################
    '''
    Convert raw datasets into SounderPy 'clean_data' dicts
    '''

    def parse_data(raw_data: dataset, latlon: list, box_avg_size: float, hush: bool):

        r"""Get model reanalysis vertical profile data

           :param raw_data: raw datasets from data retrieval methods above
           :type raw_data: dataset, required
           :param latlon: list of lat\lon points
           :type latlon: list, required
           :param box_avg_size: optional, determine an area-averaged box size in degrees, default is 0.10 degrees
           :type raw_data: float, required

           :return: clean_data, a dict of ready-to-use vertical profile data including pressure, height, temperature, dewpoint, u-wind, v-wind, & model information
           :rtype: dict
        """

        # if dataset is a xarray.core dataset, it came from the ERA5
        # and specific processing of this data is needed
        if str(type(raw_data)) == "<class 'xarray.core.dataset.Dataset'>":
            vert_data = {
                'vert_T': np.flip(np.mean(np.array(raw_data['t'][0, :, :, :] - 273.15), axis=(1, 2))),
                'vert_p': np.flip(np.array(raw_data['pressure_level'])),
                'vert_z': np.flip(np.mean(np.array(raw_data['z'][0]) / 9.80665, axis=(1, 2))),
                'vert_rh': np.flip(np.mean(np.array(raw_data['r'][0]), axis=(1, 2))),
                'vert_u': np.flip(np.mean(np.array(raw_data['u'][0]) * 1.94384, axis=(1, 2))),
                'vert_v': np.flip(np.mean(np.array(raw_data['v'][0]) * 1.94384, axis=(1, 2))),
                'vert_omega': np.flip(np.mean(np.array(raw_data['w'][0]), axis=(1, 2))),
                'vert_Td': np.flip((mpcalc.dewpoint_from_relative_humidity(
                    np.mean(np.array(raw_data['t'][0, :, :, :] - 273.15), axis=(1, 2)) * units.degC,
                    np.mean(np.array(raw_data['r'][0]), axis=(1, 2)) * units.percent)).m),
            }
            sfc_data = {
                'sfc_T': np.mean(np.array(raw_data['Ts'][0, :]) - 273.15),
                'sfc_p': np.mean(np.array(raw_data['ps'][0, :]) / 100),
                'sfc_z': np.mean(np.array(raw_data['hgts'][0, :]) / 9.80665),
                'sfc_Td': np.mean(np.array(raw_data['tds'][0, :]) - 273.15),
                'sfc_rh': (mpcalc.relative_humidity_from_dewpoint(
                    np.mean(np.array(raw_data['Ts'][0, :]) - 273.15) * units.degC,
                    np.mean(np.array(raw_data['tds'][0, :]) - 273.15) * units.degC) * 100),
                'sfc_u': (np.array(raw_data['us'][0, :]) * 1.94384),
                'sfc_v': (np.array(raw_data['vs'][0, :]) * 1.94384),
                'sfc_omega': (0)
            }

            # parse out raw dataset date and time
            vtime = np.datetime_as_string(raw_data.valid_time.values[0])
            strftime = [vtime[0:4], vtime[5:7], vtime[8:10], vtime[11:13]]

            latlon_data = {
                'data_lat': (raw_data['latitude'][:]),
                'data_lon': (raw_data['latitude'][:]),
                'data_latnum': (raw_data['latitude'][:]).shape[0],
                'data_lonnum': (raw_data['latitude'][:]).shape[0],
                'data_time': strftime
            }

            # sfc_data['sfc_w'] = (np.zeros((latlon_data['data_lonnum'], latlon_data['data_lonnum'])))



        # if data is a netCDF4 dataset, it is RAP, RUC or NCEP data
        if str(type(raw_data)) == "<class 'netCDF4._netCDF4.Dataset'>":
            # if Geopotential_height_surface exists within the dataset, it is reanalysis data
            if "Geopotential_height_surface" in raw_data.variables.keys():

                # try to determine how many isobaric levels exist in the dataset
                # this is mainly for the NCEP-FNL dataset which annoyingly
                # comes with varying isobaric level intervals

                try:
                    pressures = (ma.getdata(raw_data.variables['isobaric'])).data / 100
                except:
                    try:
                        pressures = (ma.getdata(raw_data.variables['isobaric1'])).data / 100
                    except:
                        try:
                            pressures = (ma.getdata(raw_data.variables['isobaric2'])).data / 100
                        except:
                            try:
                                pressures = (ma.getdata(raw_data.variables['isobaric3'])).data / 100
                            except:
                                pass
                            pass
                        pass
                    pass

                # RUC pressure comes in units of hPa instead of Pa like RAP data
                # so check if the values are in hPa or
                if (pressures[0] / 100) < 1:
                    pressures = pressures*100


                # create a dict of vertical data
                vert_data = {
                    'vert_T': np.mean(ma.getdata(raw_data.variables['Temperature_isobaric'][0, :, :, :] - 273.15),
                                      axis=(1, 2)),
                    'vert_p': pressures,
                    'vert_z': np.mean(ma.getdata(raw_data.variables['Geopotential_height_isobaric'][0, :, :, :]),
                                      axis=(1, 2)),
                    'vert_rh': np.mean(ma.getdata(raw_data.variables['Relative_humidity_isobaric'][0, :, :, :]),
                                       axis=(1, 2)),
                    'vert_u': np.mean(
                        ma.getdata(raw_data.variables['u-component_of_wind_isobaric'][0, :, :, :] * 1.94384),
                        axis=(1, 2)),
                    'vert_v': np.mean(
                        ma.getdata(raw_data.variables['v-component_of_wind_isobaric'][0, :, :, :] * 1.94384),
                        axis=(1, 2)),
                    'vert_omega': np.mean(
                        ma.getdata(raw_data.variables['Vertical_velocity_pressure_isobaric'][0, :, :, :]),
                        axis=(1, 2)),
                    'vert_Td': (mpcalc.dewpoint_from_relative_humidity(
                        np.mean(ma.getdata(raw_data.variables['Temperature_isobaric'][0, :, :, :] - 273.15),
                                axis=(1, 2)) * units.degC,
                        np.mean(ma.getdata(raw_data.variables['Relative_humidity_isobaric'][0, :, :, :]),
                                axis=(1, 2)) * units.percent)).m,
                }

                # create a dict of sfc data
                sfc_data = {
                    'sfc_T': np.mean(
                        ma.getdata(raw_data.variables['Temperature_height_above_ground'][0, 0, :, :] - 273.15)),
                    'sfc_p': np.mean(ma.getdata(raw_data.variables['Pressure_surface'][0, :, :] / 100)),
                    'sfc_z': np.mean(ma.getdata(raw_data.variables['Geopotential_height_surface'][0, :, :])),
                    'sfc_rh': np.mean(
                        ma.getdata(raw_data.variables['Relative_humidity_height_above_ground'][0, 0, :, :])),
                    'sfc_u': np.mean(ma.getdata(
                        raw_data.variables['u-component_of_wind_height_above_ground'][0, 0, :, :] * 1.94384)),
                    'sfc_v': np.mean(ma.getdata(
                        raw_data.variables['v-component_of_wind_height_above_ground'][0, 0, :, :] * 1.94384)),
                    'sfc_omega': (0),
                    'sfc_Td': (mpcalc.dewpoint_from_relative_humidity(
                        np.mean(ma.getdata(
                            raw_data.variables['Temperature_height_above_ground'][0, 0, :, :] - 273.15)) * units.degC,
                        np.mean(ma.getdata(raw_data.variables['Relative_humidity_height_above_ground'][0, 0, :,
                                           :]) * units.percent))).m,
                }

                # parse out raw data date and time
                dtime = raw_data.variables['Pressure_surface'].dimensions[0]
                vtime = netCDF4.num2date(raw_data.variables[dtime][:], raw_data.variables[dtime].units)[0]
                strftime = [vtime.strftime('%Y'), vtime.strftime('%m'), vtime.strftime('%d'), vtime.strftime('%H')]

                # create a dict of lat/lon information
                latlon_data = {
                    'data_lat': (ma.getdata(
                        1000 * (raw_data.variables[raw_data.variables['Pressure_surface'].dimensions[1]][:]))),
                    'data_lon': (ma.getdata(
                        1000 * (raw_data.variables[raw_data.variables['Pressure_surface'].dimensions[2]][:]))),
                    'data_latnum': (ma.getdata(
                        1000 * (raw_data.variables[raw_data.variables['Pressure_surface'].dimensions[1]][:]))).shape[0],
                    'data_lonnum': (ma.getdata(
                        1000 * (raw_data.variables[raw_data.variables['Pressure_surface'].dimensions[2]][:]))).shape[0],
                    'data_time': strftime
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
                    'vert_T': np.mean(ma.getdata(raw_data.variables['Temperature_isobaric'][0, :, :, :] - 273.15),
                                      axis=(1, 2)),
                    'vert_p': pressures,
                    'vert_z': np.mean(ma.getdata(raw_data.variables['Geopotential_height_isobaric'][0, :, :, :]),
                                      axis=(1, 2)),
                    'vert_rh': np.mean(ma.getdata(raw_data.variables['Relative_humidity_isobaric'][0, :, :, :]),
                                       axis=(1, 2)),
                    'vert_u': np.mean(
                        ma.getdata(raw_data.variables['u-component_of_wind_isobaric'][0, :, :, :] * 1.94384),
                        axis=(1, 2)),
                    'vert_v': np.mean(
                        ma.getdata(raw_data.variables['v-component_of_wind_isobaric'][0, :, :, :] * 1.94384),
                        axis=(1, 2)),
                    'vert_omega': np.mean(
                        ma.getdata(raw_data.variables['Vertical_velocity_pressure_isobaric'][0, :, :, :]),
                        axis=(1, 2)),
                    'vert_Td': (mpcalc.dewpoint_from_relative_humidity(
                        np.mean(ma.getdata(raw_data.variables['Temperature_isobaric'][0, :, :, :] - 273.15),
                                axis=(1, 2)) * units.degC,
                        np.mean(ma.getdata(raw_data.variables['Relative_humidity_isobaric'][0, :, :, :]),
                                axis=(1, 2)) * units.percent)).m,
                }

                # create dict of surface data
                sfc_data = {
                    'sfc_T': np.mean(
                        ma.getdata(raw_data.variables['Temperature_height_above_ground'][0, 0, :, :] - 273.15)),
                    'sfc_p': np.mean(ma.getdata(raw_data.variables['Pressure_surface'][0, :, :] / 100)),
                    'sfc_rh': np.mean(
                        ma.getdata(raw_data.variables['Relative_humidity_height_above_ground'][0, 0, :, :])),
                    'sfc_u': np.mean(ma.getdata(
                        raw_data.variables['u-component_of_wind_height_above_ground'][0, 0, :, :] * 1.94384)),
                    'sfc_v': np.mean(ma.getdata(
                        raw_data.variables['v-component_of_wind_height_above_ground'][0, 0, :, :] * 1.94384)),
                    'sfc_omega': (0),
                    'sfc_Td': (mpcalc.dewpoint_from_relative_humidity(
                        np.mean(ma.getdata(
                            raw_data.variables['Temperature_height_above_ground'][0, 0, :, :] - 273.15)) * units.degC,
                        np.mean(ma.getdata(raw_data.variables['Relative_humidity_height_above_ground'][0, 0, :,
                                           :]) * units.percent))).m,
                }

                # calculate surface heights
                sfc_data['sfc_z'] = np.interp([sfc_data['sfc_p']], vert_data['vert_p'], vert_data['vert_z'])[0]

                # parse out raw data date and time
                dtime = raw_data.variables['Pressure_surface'].dimensions[0]
                vtime = netCDF4.num2date(raw_data.variables[dtime][:], raw_data.variables[dtime].units)[0]
                strftime = [vtime.strftime('%Y'), vtime.strftime('%m'), vtime.strftime('%d'), vtime.strftime('%H')]

                # create dict of latlon information
                latlon_data = {
                    'data_lat': (ma.getdata(
                        1000 * (raw_data.variables[raw_data.variables['Pressure_surface'].dimensions[1]][:]))),
                    'data_lon': (ma.getdata(
                        1000 * (raw_data.variables[raw_data.variables['Pressure_surface'].dimensions[2]][:]))),
                    'data_latnum': (ma.getdata(
                        1000 * (raw_data.variables[raw_data.variables['Pressure_surface'].dimensions[1]][:]))).shape[0],
                    'data_lonnum': (ma.getdata(
                        1000 * (raw_data.variables[raw_data.variables['Pressure_surface'].dimensions[2]][:]))).shape[0],
                    'data_time': strftime
                }

        sb_dict = {}
        new_keys = ['T', 'Td', 'rh', 'u', 'v', 'z', 'p', 'omega']
        sfc_keys = ['sfc_T', 'sfc_Td', 'sfc_rh', 'sfc_u', 'sfc_v', 'sfc_z', 'sfc_p', 'sfc_omega']
        vert_keys = ['vert_T', 'vert_Td', 'vert_rh', 'vert_u', 'vert_v', 'vert_z', 'vert_p', 'vert_omega']

        # create a dict of surface-based data
        for vert_key, sfc_key, new_key in zip(vert_keys, sfc_keys, new_keys):
            sb_dict[new_key] = np.insert(
                np.flip(vert_data[vert_key])[np.flip(vert_data['vert_z']) >= sfc_data['sfc_z']], 0, sfc_data[sfc_key])
        sb_dict['z'] = sb_dict['z'] - sb_dict['z'][0]

        # Interpolates data
        dz = 250
        soundingtop_hght = sb_dict['z'][-1]
        toplvl = int(soundingtop_hght / dz) * dz
        numlvls = int(toplvl / dz)
        interp_lvls = np.linspace(0, toplvl, numlvls + 1)

        # prepare new dicts
        keys = ['T', 'Td', 'rh', 'u', 'v', 'z', 'p', 'omega']
        units_list = ['degC', 'degC', 'percent', 'kt', 'kt', 'm', 'hPa', 'Pa/sec']
        interp_dict = {}
        zeros_dict = {}
        clean_data = {}

        surface_height = sfc_data['sfc_z']

        # create dict of clean data
        for key in keys:
            interp_dict[key] = (interpolate.interp1d(sb_dict['z'], sb_dict[key]))
            zeros_dict[key] = np.zeros((len(interp_lvls)))
            for zeros_arr in zeros_dict.values():
                zeros_dict[key][0] = sb_dict[key][0]
            for i in range(1, len(zeros_dict[key]), 1):
                zeros_dict[key][i] = interp_dict[key](dz * i)

        for i, unit, key in zip(range(0, len(units_list)), units_list, keys):
            clean_data[key] = zeros_dict[key] * units(unit)

        # add site & profile info to clean_data dict
        clean_data['site_info'] = {
            'site-id': 'no-site-id',
            'site-name': 'no-site-name',
            'site-lctn': 'no-site-location',
            'site-latlon': [latlon[0], latlon[1]],
            'site-elv': surface_height,
            'source': f'MODEL {str.upper(dtype)}',
            'model': source,
            'fcst-hour': 'F00',
            'run-time': latlon_data['data_time'],
            'valid-time': latlon_data['data_time'],
            'box_area': f'{box_avg_size}° BOX AVG'}

        # add plot title to clean_data dict
        clean_data['titles'] = {
            'top_title': f"MODEL REANALYSIS VERTICAL PROFILE | {clean_data['site_info']['valid-time'][3]}Z {clean_data['site_info']['model']} {clean_data['site_info']['fcst-hour']}",
            'left_title': f"{clean_data['site_info']['run-time'][3]}Z {clean_data['site_info']['model']} {clean_data['site_info']['fcst-hour']} | VALID: {clean_data['site_info']['valid-time'][1]}/{clean_data['site_info']['valid-time'][2]}/{clean_data['site_info']['valid-time'][0]} {clean_data['site_info']['valid-time'][3]}Z",
            'right_title': f"{clean_data['site_info']['site-latlon'][0]}, {clean_data['site_info']['site-latlon'][1]} | {clean_data['site_info']['box_area']}    "
        }


        print('    > COMPLETE --------')
        elapsed_time = time.time() - st
        print('    > RUNTIME:', time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))

        if not hush:
            print(
                f"    > SUMMARY: {clean_data['site_info']['run-time'][3]}Z {clean_data['site_info']['model']} {clean_data['site_info']['fcst-hour']} for"
                f"{clean_data['site_info']['site-latlon']} at {clean_data['site_info']['valid-time'][1]}-{clean_data['site_info']['valid-time'][2]}-{clean_data['site_info']['valid-time'][0]}-{clean_data['site_info']['valid-time'][3]}Z")

            warnings.filterwarnings("ignore")

            sounding_params(clean_data).print_vals()

        return clean_data

    if clean_it:

        return parse_data(raw_data, latlon, box_avg_size, hush)

    else:
        # if user sets `clean_it=False` return raw data ('dirty') format
        return raw_data
from pyproj import Geod
import numpy as np
from matplotlib.colors import ListedColormap
from datetime import datetime, timedelta
import fsspec

import os
from contextlib import contextmanager, redirect_stderr, redirect_stdout
@contextmanager
def suppress_stdout_stderr():
    """A context manager that redirects stdout and stderr to devnull"""
    with open(os.devnull, 'w') as fnull:
        with redirect_stderr(fnull) as err, redirect_stdout(fnull) as out:
            yield (err, out)

with suppress_stdout_stderr():
    import pyart

"""
    SOUNDERPY RADAR UTILITY FUNCTIONS 

    Purpose of module: 

     A collection of radar-related helper functions for plotting radar data
     on SounderPY figures


    (C) KYLE J GILLETT, UNIVERSITY OF NORTH DAKOTA, 2025
"""


#########################
# RADAR DATA DATE HANDLER
#########################################################################
def define_radar_time(clean_data, radar_time):
    data_dt = datetime(int(clean_data['site_info']['valid-time'][0]), int(clean_data['site_info']['valid-time'][1]),
                       int(clean_data['site_info']['valid-time'][2]),
                       int(clean_data['site_info']['valid-time'][3][0:2]))

    # get current time
    curr_dt = datetime.utcnow()

    # if data time is after current time, default to now
    if data_dt > curr_dt:
        time = curr_dt.strftime('%H%M')
        datestr = curr_dt.strftime('%Y%m%d')
        dt_obj = curr_dt

    # if radar time is given as now, time is now
    elif radar_time == 'now':
        time = curr_dt.strftime('%H%M')
        datestr = curr_dt.strftime('%Y%m%d')
        dt_obj = curr_dt

    # if radar time is given as sounding, use the valid-time info from clean_data
    elif radar_time == 'sounding':
        if len(clean_data['site_info']['valid-time'][3]) == 2:
            time = f"{clean_data['site_info']['valid-time'][3]}00"
        else:
            time = clean_data['site_info']['valid-time'][3]
        datestr = data_dt.strftime('%Y%m%d')
        dt_obj = data_dt

    else:
        time = radar_time.strftime('%H%M')
        datestr = radar_time.strftime('%Y%m%d')
        dt_obj = radar_time

    return time, datestr, dt_obj
#########################################################################



#########################
# CREATE RADARSCOPE COLORMAP
#########################################################################
script_dir = os.path.dirname(os.path.abspath(__file__))
txt_file_path = os.path.join(script_dir, "rs_reflectivity.txt")
rs_data = np.loadtxt(txt_file_path, skiprows=3, usecols=(1, 2, 3, 4))
rgb = rs_data[:, 1:] / 255.0
levels = rs_data[:, 0]
rs_expertreflect_cmap = ListedColormap(rgb)
#########################################################################



#########################
# DMS to decimal coords
#########################################################################
def dms_to_decimal(deg, minutes, seconds):
    return deg + minutes / 60 + seconds / 3600
#########################################################################



#########################
# PARSE DMS STRING FROM FILE
#########################################################################
def parse_dms_string(dms_str, is_lat=True):
    if is_lat:
        deg = int(dms_str[0:2])
        minutes = int(dms_str[2:4])
        seconds = int(dms_str[4:6])
    else:
        deg = int(dms_str[0:3])
        minutes = int(dms_str[3:5])
        seconds = int(dms_str[5:7])
    return dms_to_decimal(deg, minutes, seconds)
#########################################################################


#########################
# OPEN AND PARSE NEXRAD LOCS FILE
#########################################################################
def parse_nexrad_locs(filename):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    txt_file_path = os.path.join(script_dir, filename)
    nexrad_locs = np.genfromtxt(txt_file_path, dtype='str', delimiter=',', skip_header=1)
    site_ids = nexrad_locs[:, 0]
    latitudes = []
    longitudes = []
    for site in nexrad_locs[:, 1]:
        lat = parse_dms_string(site[0:6], is_lat=True)
        lon = -parse_dms_string(site[9:16], is_lat=False)
        latitudes.append(lat)
        longitudes.append(lon)

    return site_ids, latitudes, longitudes
#########################################################################



#########################
# FIND NEAREST NEXRAD SITE
#########################################################################
def find_nearest_station(site_ids, lats, lons, target_lat, target_lon, max_distance_km=400):
    g = Geod(ellps='sphere')
    distances_km = np.array([g.inv(target_lon, target_lat, lon, lat)[2] / 1000.0 for lat, lon in zip(lats, lons)])

    if np.all(np.isnan(distances_km)):
        return None

    min_dist = np.nanmin(distances_km)
    min_idx = np.nanargmin(distances_km)

    if min_dist < max_distance_km:
        site = site_ids[min_idx]
        site_lat = lats[min_idx]
        site_lon = lons[min_idx]
        # Handle special site replacement rule
        if site == 'KJAN':
            site = 'KDGX'
        return site[:], site_lat, site_lon
    return None, None, None
#########################################################################



#########################
# LOAD NEXRAD DATA FROM AWS
#########################################################################
def get_radar_data(nexrad_site, clean_data, radar_time):

    time, datestr, dt_obj = define_radar_time(clean_data, radar_time)

    fs = fsspec.filesystem('s3', anon=True)
    # try current time, then fallback to previous hour if needed
    for hour_offset in [0, -1]:
        radar_dt = dt_obj + timedelta(hours=hour_offset)
        time_candidates = [radar_dt + timedelta(minutes=delta) for delta in range(-4, 5)]

        files = []
        for t in time_candidates:
            prefix = t.strftime(f"s3://noaa-nexrad-level2/%Y/%m/%d/{nexrad_site}/{nexrad_site}%Y%m%d_%H%M")
            matches = fs.glob(prefix + "*")
            matches = [f for f in matches if not f.endswith('_MDM')]
            files.extend(matches)

        if files:
            files = sorted(files)
            file = files[0]
            break
    else:
        raise FileNotFoundError("    > RADAR ERROR: No single-site radar scan found")

    radar = pyart.io.read_nexrad_archive("s3://" + file)

    # extract timestamp from filename
    scan_time = str(files[0][48:50] + ':' + files[0][50:52])
    scan_date = str(files[0][39:43] + '-' + files[0][43:45] + '-' + files[-1][45:47])
    scan_timestamp = str(scan_date + ' | ' + scan_time)
    print(f"    > RADAR SCAN FOUND: Site: {nexrad_site} Scan: {scan_timestamp}")

    return radar, scan_timestamp, nexrad_site
#########################################################################



#########################
# LOAD RADAR MOSAIC DATA
#########################################################################
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

    # GET TIME INFO -----------------------------------------------------------------------------------
    time, datestr, dt_obj = define_radar_time(clean_data, radar_time)

    # PULL SOUNDING DATA LAT/LON ----------------------------------------------------------------------
    data_lat = clean_data['site_info']['site-latlon'][0]
    data_lon = clean_data['site_info']['site-latlon'][1]

    # GET RADAR DATA ----------------------------------------------------------------------------------
    composite_url = 'https://thredds.ucar.edu/thredds/catalog/nexrad/composite/gini/dhr/1km/' + datestr + '/catalog.xml'
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
        composite_query.lonlat_box(north=data_lat + (map_zoom + 1),
                                   south=data_lat - (map_zoom + 1),
                                   east=data_lon + (map_zoom + 1),
                                   west=data_lon - (map_zoom + 1))
        composite_query.add_lonlat(value=True)
        composite_query.accept('netcdf4')
        composite_query.variables('Reflectivity')
        radar_data = filename.get_data(composite_query)
        radar_data = open_dataset(NetCDF4DataStore(radar_data))

        return radar_data
#########################################################################
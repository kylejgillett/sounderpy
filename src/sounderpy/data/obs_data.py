import re
import struct
import time
import warnings
import zlib
from datetime import datetime, timezone

import pandas as pd
import requests
import numpy as np

from metpy.units import units
import metpy.calc as mpcalc

from siphon.simplewebservice.wyoming import WyomingUpperAir
from siphon.simplewebservice.igra2 import IGRAUpperAir

from ..utils import find_raob_station, raob_station_has_no_known_data, split_display_name
from ..calc import sounding_params







"""
    SOUNDERPY RAOB OBSERVATIONS 'GET-DATA' FUNCTIONS  

    Purpose of module: 

    House function for loading and parsing RAOB & IGRAv2 observation
    data. Functions here are referenced by sounderpy.py


    (C) KYLE J GILLETT, UNIVERSITY OF NORTH DAKOTA, 2024, 2026
"""




def _choose_profile_stride(values, target_len=300):
    """Return a slicing stride that keeps an array length near ``target_len``."""
    length = len(values)
    if length <= target_len:
        return 1

    interval = max(1, int(np.round(length / target_len)))
    best_interval = interval
    best_diff = abs(int(np.ceil(length / interval)) - target_len)

    for candidate in [max(1, interval - 1), interval, interval + 1]:
        sampled_len = int(np.ceil(length / candidate))
        diff = abs(sampled_len - target_len)
        if diff < best_diff:
            best_diff = diff
            best_interval = candidate

    return best_interval




#######################
# IGRAv2 STATION FILE READER
#########################################################################
# siphon's IGRAUpperAir downloads a station's whole period-of-record zip (82 MB for 72357), unzips all of it
# into a list of Python strings (243 MB of text, 4.5 million lines, ~1.1 GB of memory) and then keeps one
# sounding. Here the zip is decompressed while it downloads and only the requested sounding is kept, so memory
# stays at a few MB; the download stops once the sounding is read or the file has passed its date (the file
# is in time order); and recent dates are read from IGRA's much smaller year-to-date file.

IGRA_URL = 'https://www.ncei.noaa.gov/pub/data/igra/data/'


class IGRANoData(ValueError):
    '''IGRAv2 answered, but has no sounding for this station and time'''


class _IGRAFormatError(Exception):
    '''the station file isn't a zip this reader can stream; siphon's full-file reader is used instead'''


class _IGRASounding(IGRAUpperAir):
    '''siphon's IGRA parser, fed the lines of one sounding instead of the whole station file'''

    def __init__(self, lines):
        super().__init__()
        self._lines = lines

    def _get_data_raw(self):
        return self._select_date_range(self._lines)


def _unzip_stream(chunks):
    '''Yield the decompressed bytes of the first file in a zip archive, from an iterator of its raw bytes.'''
    head = b''
    for part in chunks:
        head += part
        if len(head) >= 30:
            break
    if len(head) < 30:
        raise _IGRAFormatError('station file is empty')
    sig, _, _, method, _, _, _, _, _, name_len, extra_len = struct.unpack('<IHHHHHIIIHH', head[:30])
    if sig != 0x04034b50 or method != 8:   # local file header, deflate
        raise _IGRAFormatError('station file is not a deflated zip')
    start = 30 + name_len + extra_len
    while len(head) < start:
        part = next(chunks, None)
        if part is None:
            raise _IGRAFormatError('station file is truncated')
        head += part
    unzip = zlib.decompressobj(-zlib.MAX_WBITS)
    out = unzip.decompress(head[start:])
    if out:
        yield out
    for part in chunks:
        if unzip.eof:
            return
        out = unzip.decompress(part)
        if out:
            yield out


def _find_igra_sounding(chunks, station, when):
    '''
    Return the lines (header first) of the sounding at ``when`` from an IGRAv2 station file, given as an
    iterator of its decompressed bytes, or ``None`` if the file doesn't have it.

    Like siphon, a sounding with no launch hour (``99``) counts as 00z.
    '''
    key = f'#{station} {when:%Y %m %d} '.encode()
    hours = {f'{when.hour:02d}'.encode()} | ({b'99'} if when.hour == 0 else set())
    day = f'{when:%Y %m %d}'.encode()
    buf = b''
    for chunk in chunks:
        buf += chunk
        end = buf.rfind(b'\n') + 1          # search whole lines only
        pos = buf.find(key, 0, end)
        while pos >= 0:
            if (pos == 0 or buf[pos - 1] == 10) and buf[pos + len(key):pos + len(key) + 2] in hours:
                block = buf[pos:]
                n_lines = int(block[32:36]) + 1  # header + its levels
                while block.count(b'\n') < n_lines:
                    more = next(chunks, None)
                    if more is None:
                        break
                    block += more
                return [line.decode('utf-8') + '\n' for line in block.split(b'\n')[:n_lines]]
            pos = buf.find(key, pos + 1, end)
        # soundings are in time order: once a header is past the requested day, the file doesn't have it
        last = buf.rfind(b'#', 0, end)
        if last >= 0 and (last == 0 or buf[last - 1] == 10) and buf[last + 13:last + 23] > day:
            return None
        buf = buf[end:]
    return None


def _igra_station_files(station, when, now=None):
    '''
    URLs of the IGRAv2 files that could hold ``when``, smallest first. The first one that exists is
    the only one read: a year-to-date file is updated as often as the period-of-record file.
    '''
    now = now or datetime.now(timezone.utc)
    # data-y2d holds everything since January 1 of the current or previous year, named for that year
    files = [f'{IGRA_URL}data-y2d/{station}-data-beg{year}.txt.zip'
             for year in (now.year, now.year - 1) if year <= when.year]
    return files + [f'{IGRA_URL}data-por/{station}-data.txt.zip']


def request_igra(station, when, chunk_size=1 << 20, timeout=(15, 60)):
    '''
    Return siphon's ``(body, header)`` DataFrames for one IGRAv2 sounding, reading only as much of the
    station file as needed. Raises ``IGRANoData`` when IGRAv2 doesn't have the sounding.
    '''
    for url in _igra_station_files(station, when):
        with requests.get(url, stream=True, timeout=timeout) as resp:
            if resp.status_code == 404:
                continue   # no year-to-date file under this year's name (or no station file at all)
            resp.raise_for_status()
            chunks = _unzip_stream(resp.iter_content(chunk_size))
            lines = _find_igra_sounding(chunks, station, when)
        if lines is None:
            break
        igra = _IGRASounding(lines)
        igra.begin_date = igra.end_date = when
        igra.site_id = station
        return igra._get_data()
    raise IGRANoData(f'No IGRAv2 sounding for station {station} at {when:%Y-%m-%d %HZ}')


#######################
# RAOB AND IRGAv2 OBSERVATIONS DATA
#########################################################################

def _data_source(row, archive):
    '''Footer text naming the archive a profile came from, plus the station's IDs and IGRAv2 record years'''
    ids = ' / '.join(f'{name} {row[key]}' for name, key in (('WMO', 'WMO'), ('ICAO', 'ICAO'), ('IGRA2', 'IGRA_ID')) if row[key])
    record = ''
    if row['IGRA_FIRST_YEAR'] and row['IGRA_LAST_YEAR']:
        record = f"RECORD {row['IGRA_FIRST_YEAR']}–{row['IGRA_LAST_YEAR']}" + (' (ACTIVE)' if row['ACTIVE'] == 'Y' else '')
    return ' · '.join(part for part in (archive, ids, record) if part)


def _requested_icao(station, row):
    '''
    The ICAO code the user asked for, as the station list spells it: ``'NKX'`` -> ``'KNKX'`` even when the row's
    main code is ``'KMYF'``. Falls back to the row's main ICAO when the request matched some other way.
    '''
    codes = [row['ICAO']] + row['ICAO_ALT'].split()
    for code in (station, 'K' + station):
        if code in codes:
            return code
    return row['ICAO']


def _site_info(row, site_id, year, month, day, hour, archive):
    '''site_info for an observed profile, from a station list row (see utils.find_raob_station)'''
    name, region = split_display_name(row['DISPLAY_NAME'] or row['NAME'])
    has_pos = bool(row['LAT'] and row['LON'])
    return {
        'site-id': site_id,
        'site-name': name,
        'site-lctn': region,
        'site-latlon': [round(float(row['LAT']), 2), round(float(row['LON']), 2)] if has_pos else [np.nan, np.nan],
        'site-elv': float(row['ELEV_M']) if row['ELEV_M'] else np.nan,
        'source': 'RAOB OBSERVED PROFILE',
        'model': 'no-model',
        'fcst-hour': 'no-fcst-hour',
        'run-time': ['none', 'none', 'none', 'none'],
        'valid-time': [year, month, day, hour],
        'data-source': _data_source(row, archive)}


def fetch_obs(station, year, month, day, hour, hush, clean_it):

    # record process time
    st = time.time()

    r"""
       Return a ``dict`` of 'cleaned up' observed profile data
    
       :param station: a 3-4 digit RAOB identifier (such as: 'KDTX') or 11 digit IGRAv2 identifier (such as: 'GMM00010393')
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

    print(f'> OBSERVED DATA ACCESS FUNCTION\n  -----------------------------------')

    station = str.upper(station).strip()

    # look the station up in SounderPy's packaged station list (WMO, ICAO & IGRAv2 IDs)
    try:
        row = find_raob_station(station)
    except ValueError as e:
        raise ValueError(
            f"{e}\nIf you think this is an error, contact the author: "
            "https://kylejgillett.github.io/sounderpy/user_guide/about.html#about-the-author"
        ) from None

    if raob_station_has_no_known_data(row):
        warnings.warn(
            f"Station '{station}' ({row['DISPLAY_NAME'] or row['NAME']}) has no IGRAv2 record and hasn't reported to "
            "the Wyoming archive recently, so data may not be available. Older soundings may still be in "
            "Wyoming's archive.", stacklevel=2)

    got_data = False

    # an 11-character IGRAv2 ID ('GMM00010393') is fetched from IGRAv2; anything else from UW
    search_for = 'igra' if len(station) == 11 else 'raob'



    ### RAOB OBSERVATIONS ###
    #########################################################################################################
    if search_for == 'raob':
        station_icao = row['ICAO']
        station_wmo  = row['WMO']      # text, so leading zeros are kept ('01001')
        if not (station_wmo or station_icao):
            raise ValueError(
                f"Station '{station}' has no WMO or ICAO ID, so it can't be requested from the Wyoming archive. "
                f"Try its IGRAv2 ID instead: '{row['IGRA_ID']}'")

        # UW sometimes files a station under its ICAO ID only (e.g. 42410 / VEGT), so a WMO
        # request can be rejected even when the data exists. Try the WMO ID first, then the ICAO ID.
        request_ids = [i for i in (station_wmo, station_icao) if i]

        request_time = datetime(int(year), int(month), int(day), int(hour))
        # why each ID failed: 'no data' when UW answered without a profile, otherwise the server problem
        failures = []
        for request_id in request_ids:
            # try up to 10 times, UW often returns temporary 'server busy' errors
            for i in range(1, 11):
                try:
                    df = WyomingUpperAir.request_data(request_time, request_id)
                    got_data = True
                    break
                except Exception as e:
                    got_data = False
                    # a 400/404 means UW has no data under this ID -- retrying won't help.
                    # siphon only puts the status code in the message: 'Server Error ( 400: ...)'
                    cause = e.__cause__ or e
                    status = re.search(r'Server Error \(\s*(\d{3})', str(cause))
                    if status and status.group(1) in ('400', '404'):
                        failure = 'no data'
                        break
                    # other HTTP errors (503 'server busy', 5xx) and no response at all mean UW couldn't answer;
                    # anything else (e.g. a reply with no profile in it) is treated as no data
                    if status:
                        failure = f'HTTP {status.group(1)}'
                    elif isinstance(cause, requests.RequestException):
                        failure = 'no response'
                    else:
                        failure = 'no data'
                    # UW's 'server busy' (503) spells last a few seconds; back-to-back retries all fail
                    if i < 10:
                        time.sleep(2)
            if not got_data:
                failures.append(failure)
            if got_data:
                print(f'    > PROFILE FOUND: {station} on {month}/{day}/{year} at {hour}z | From UW (ID: {request_id})')
                break

        # search through RAOB sites list with provided RAOB ID, first try ICAO ID, then WMO ID
        if got_data:
            if clean_it:
                # begin loading data
                # create dict of data
                new_keys = ['p', 'z', 'T', 'Td', 'u', 'v']
                old_keys = ['pressure', 'height', 'temperature', 'dewpoint', 'u_wind', 'v_wind']  # 'latitude', 'longitude']
                units_list = ['hPa', 'meter', 'degC', 'degC', 'kt', 'kt']
                clean_data = {}
                non_dups = np.concatenate(([True], np.diff(df.to_dict('list')['pressure']) != 0))
                profile_values = np.array(df.to_dict('list')['pressure'])[non_dups]
                stride = _choose_profile_stride(profile_values)

                for old_key, new_key, unit in zip(old_keys, new_keys, units_list):
                    if old_key in ['u_wind', 'v_wind']:
                        clean_data[new_key] = (np.array(df.to_dict('list')[old_key])[non_dups][::stride] * 1.94384 )* units(unit)
                    else:   
                        clean_data[new_key] = np.array(df.to_dict('list')[old_key])[non_dups][::stride] * units(unit)

                # label the site the way the user asked for it: by WMO number or by the ICAO code they typed
                site_id = station_wmo if station.isdigit() or not station_icao else _requested_icao(station, row)
                clean_data['site_info'] = _site_info(row, site_id, year, month, day, hour, 'UNIV. OF WYOMING')

                try:
                    # trim data to 98hPa and below for less process time
                    slc = (len(clean_data['p']) - np.where(clean_data['p'] <= 98. * units('hPa'))[0][0])
                    for key in new_keys:
                        clean_data[key] = clean_data[key][:-slc]
                except:
                    pass

        else:
            # only call it missing data when UW answered for every ID; if it couldn't answer, the data may exist
            problems = sorted({failure for failure in failures if failure != 'no data'})
            if problems:
                raise ValueError(
                    f"The Wyoming Upper Air Archive is unavailable right now ({', '.join(problems)}), so station "
                    f"{station} on {month}/{day}/{year} at {hour}z could not be checked -- try again later")
            raise ValueError(
                f'The Wyoming Upper Air Archive has no data for station {station} on {month}/{day}/{year} at {hour}z '
                f'-- ensure you have the correct date, time, and station identifier')
    #########################################################################################################







    ### IGRAv2 OBSERVATIONS ###
    #########################################################################################################
    elif search_for == 'igra':
        request_time = datetime(int(year), int(month), int(day), int(hour))
        for i in range(1, 3):
            try:
                try:
                    df = request_igra(station, request_time)
                except _IGRAFormatError:
                    # a station file this reader can't stream: let siphon read the whole thing
                    df = IGRAUpperAir.request_data(request_time, station)
                got_data = True
                print(f'    > PROFILE FOUND: {station} on {month}/{day}/{year} at {hour}z | From IGRAv2')
                break
            except IGRANoData:
                # IGRAv2 answered without the sounding; asking again would download the file again for nothing
                got_data = False
                break
            except Exception:
                got_data = False

        # if data is found, parse data and create a dict of clean data
        if got_data:
            if clean_it:
                station = row['IGRA_ID'] or station

                # create dict of data
                head = df[1]
                df = df[0]
                new_keys = ['p', 'z', 'T', 'Td', 'u', 'v']
                old_keys = ['pressure', 'height', 'temperature', 'dewpoint', 'u_wind', 'v_wind']  # 'latitude', 'longitude']
                units_list = ['hPa', 'meter', 'degC', 'degC', 'kt', 'kt']
                clean_data = {}
                zflag = np.array(df['zflag'])
                pflag = np.array(df['pflag'])
                tflag = np.array(df['tflag'])
                mask = zflag + pflag + tflag >= 4
                profile_values = np.array(df.to_dict('list')['pressure'])[mask]
                stride = _choose_profile_stride(profile_values)
                # older soundings often report relative humidity but no dewpoint depression, and siphon
                # builds the dewpoint from the depression only: fill those dewpoints from the humidity
                df = df.copy()
                rh = df['relative_humidity'].to_numpy(dtype=float)
                fill = df['dewpoint'].isna().to_numpy() & np.isfinite(rh) & (rh > 0) & df['temperature'].notna().to_numpy()
                if fill.any():
                    df.loc[fill, 'dewpoint'] = mpcalc.dewpoint_from_relative_humidity(
                        df.loc[fill, 'temperature'].to_numpy(dtype=float) * units.degC,
                        np.minimum(rh[fill], 100.) * units.percent).m
                for old_key, new_key, unit in zip(old_keys, new_keys, units_list):
                    clean_data[new_key] = np.array(df.to_dict('list')[old_key])[mask][::stride] * units(unit)
                clean_data['site_info'] = _site_info(row, station, year, month, day, hour, 'NOAA NCEI IGRA2')
                # correct u & v units
                clean_data['u'] = clean_data['u'] * 1.94384
                clean_data['v'] = clean_data['v'] * 1.94384

        else:
            raise ValueError(
                f'IGRAv2 Dataset connection failed -- ensure you have the correct dates and corresponding station identifier\n' +
                f'There is likely no available data for station {station} on {month}/{day}/{year} at {hour}z')
    #########################################################################################################




    if clean_it:

        ### CHECK AND FIX BAD HEIGHT DATA ###
        ######################################
        def are_heights_bad(arr):
            # check if heights are bad

            # inc flag
            increasing = False

            for i in range(1, len(arr)):
                if arr[i] > arr[i - 1]:
                    increasing = True
                elif increasing:  # are_heights_bad == True
                    return True

            return False


        def fix_bad_heights(heights):
            # fix bad heights

            increasing = False

            last_increasing_index = 0
            # search through height array to find the last
            # value in the sequence to 'increase' correctly
            for i in range(1, len(heights)):
                if heights[i] > heights[i - 1]:  # do so by checking if i is > than i-1
                    # record its index
                    last_increasing_index = i
                else:
                    # the last increasing value is found
                    # return its index
                    break
            last_increasing_value = heights[last_increasing_index]

            # for every 'bad' value afterward, add the last_increasing_index value
            # to it, to correct the height data
            for i in range(last_increasing_index + 1, len(heights)):
                heights[i] += last_increasing_value

                fixed_heights = heights

            return fixed_heights


        bad_heights = are_heights_bad(clean_data['z'])

        if bad_heights:
            clean_data['z'] = fix_bad_heights(clean_data['z'])



        # add plot title to clean_data dict
        clean_data['titles'] = {
            # title specific to RAOB Obs
            'top_title': "RAOB OBSERVED VERTICAL PROFILE",
            'left_title': f"VALID: {clean_data['site_info']['valid-time'][1]}-{clean_data['site_info']['valid-time'][2]}-{clean_data['site_info']['valid-time'][0]} {clean_data['site_info']['valid-time'][3]}Z",
            'right_title': f"{clean_data['site_info']['site-id']} - {clean_data['site_info']['site-name']}"
                           + (f", {clean_data['site_info']['site-lctn']}" if clean_data['site_info']['site-lctn'] else '')
                           + f" | {clean_data['site_info']['site-latlon'][0]}, {clean_data['site_info']['site-latlon'][1]}    "
        }




        print('    > COMPLETE --------')
        elapsed_time = time.time() - st
        print('    > RUNTIME:', time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))

        if not hush:
            print(
                f"    > SUMMARY: {clean_data['site_info']['valid-time'][3]}Z Launch for {clean_data['site_info']['site-id']}, {clean_data['site_info']['site-name']} at {clean_data['site_info']['valid-time'][1]}"
                f"-{clean_data['site_info']['valid-time'][2]}-{clean_data['site_info']['valid-time'][0]}-{clean_data['site_info']['valid-time'][3]}Z")

            warnings.filterwarnings("ignore")

            sounding_params(clean_data).print_vals()

        return clean_data


    else:
        # if user sets `clean_it=False` return raw data ('dirty') format
        return df

    #########################################################################################################
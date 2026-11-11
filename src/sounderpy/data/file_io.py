import os
import re
import warnings
from datetime import datetime
import csv
import numpy as np
import metpy.calc as mpcalc
from metpy.units import units



"""
    SOUNDERPY FILE READING AND WRITING FUNCTIONS 

    Purpose of module: 

     Functions for writing files CM1/SHARPPY/CSV and reading SHARPPY, GRAW, iMet, and Windsond files


    (C) KYLE J GILLETT, UNIVERSITY OF NORTH DAKOTA, 2026
"""



#########################
# FILE READING FUNCTIONS
#########################################################################
def from_sharppy(filepath):
    """
    Read a SHARPpy sounding file and return a SounderPy-style clean_data dict.
    """

    with open(filepath, "r", encoding="utf-8") as f:
        lines = [line.rstrip("\n") for line in f]

    if not lines:
        raise ValueError("The SHARPpy file is empty.")

    title_idx = None
    raw_idx = None
    for idx, line in enumerate(lines):
        if line.startswith("%TITLE%"):
            title_idx = idx
        elif line.startswith("%RAW%"):
            raw_idx = idx
            break

    if raw_idx is None:
        raise ValueError("This does not appear to be a SHARPpy sounding file.")

    site_info = {
        "site-id": None,
        "valid-time": None,
    }

    if title_idx is not None and title_idx + 1 < len(lines):
        header_line = lines[title_idx + 1].strip()
        if header_line:
            parts = re.split(r"\s+", header_line)
            if parts:
                site_info["site-id"] = parts[0]
            if len(parts) >= 2:
                stamp = parts[1]
                if "/" in stamp and len(stamp) >= 9:
                    try:
                        yy = int(stamp[0:2])
                        mm = int(stamp[2:4])
                        dd = int(stamp[4:6])
                        hh = int(stamp[7:9])
                        year = 2000 + yy if yy < 100 else yy
                        site_info["valid-time"] = [year, mm, dd, hh]
                    except ValueError:
                        pass

    if site_info["valid-time"] is None:
        try:
            dt = datetime.strptime(lines[title_idx + 1].split()[1], "%y%m%d/%H%M")
            site_info["valid-time"] = [dt.year, dt.month, dt.day, dt.hour]
        except Exception:
            pass

    rows = []
    for line in lines[raw_idx + 1:]:
        if not line or line.startswith("%"):
            if line == "%END%":
                break
            continue

        parts = [part.strip() for part in line.split(",")]
        if len(parts) < 6:
            continue

        try:
            rows.append([float(part) for part in parts[:6]])
        except ValueError:
            continue

    if not rows:
        raise ValueError("No sounding rows were found in the SHARPpy file.")

    # SHARPpy marks missing values with -9999; make them NaN, and drop
    # levels missing a pressure or height since they can't be placed
    rows = np.array(rows)
    rows[rows <= -9998] = np.nan
    rows = rows[~np.isnan(rows[:, 0]) & ~np.isnan(rows[:, 1])]

    p = rows[:, 0] * units.hPa
    z = rows[:, 1] * units.m
    T = rows[:, 2] * units.degC
    Td = rows[:, 3] * units.degC
    wd = rows[:, 4] * units.degrees
    ws = rows[:, 5] * units.kts

    u, v = mpcalc.wind_components(ws, wd)

    print("**NOTE**: to plot your SHARPPY file data with SounderPy, you'll have to manually provide additional site meta-data.\nThis process is fairly simple and is discussed in detail at https://kylejgillett.github.io/sounderpy/reference/customdatasources.html#ingesting-your-own-data-into-sounderpy.\nBe sure to include a site ID, name, location, lat/lon, and elevation. You must also create plot titles for the sounding plot.")


    return {
        "p": p,
        "z": z,
        "T": T,
        "Td": Td,
        "u": u,
        "v": v,
        "wd": wd,
        "ws": ws,
        "site_info": site_info,
    }




#########################
# RADIOSONDE FILE READERS
#########################################################################
# GRAW, iMet, and Windsond readers, built on the conversion scripts that Alec Sczepanski (University of
# North Dakota) contributed for the radiosondes UND launches. Each returns the same dict as from_sharppy(),
# with any site meta-data the file holds already in site_info.

# value the iMet software writes for a missing field
IMET_MISSING = 99999
# the iMet ascent rate (m/s) that marks the launch; earlier rows are the sonde sitting before launch
IMET_LAUNCH_ASCENT = 1.0


def _column_index(header, name, file_type):
    """Index of the first header column whose label starts with `name` (case-insensitive)."""
    for idx, label in enumerate(header):
        if label.strip().lower().startswith(name.lower()):
            return idx
    raise ValueError(f"This does not appear to be a {file_type} file (no '{name}' column was found).")


def _to_float(value, missing=None):
    """A field as a float, or NaN when it is blank, not a number, or the file's missing value."""
    try:
        number = float(str(value).strip())
    except ValueError:
        return np.nan
    return np.nan if missing is not None and number == missing else number


def _ascent_only(rows, p_col, z_col):
    """The rows of the balloon's ascent, from the first row to the burst (lowest pressure).

    Rows missing a pressure or height are dropped, as are rows whose pressure doesn't fall below every
    level under them (the balloon bobbing, or repeated pressures), so pressure always decreases upward.
    """
    rows = rows[np.isfinite(rows[:, p_col]) & np.isfinite(rows[:, z_col])]
    if len(rows) == 0:
        raise ValueError("No sounding rows with a pressure and height were found in the file.")
    rows = rows[: int(np.argmin(rows[:, p_col])) + 1]
    lowest_below = np.minimum.accumulate(np.concatenate([[np.inf], rows[:-1, p_col]]))
    return rows[rows[:, p_col] < lowest_below]


def _radiosonde_clean_data(p, z, T, Td, ws_ms, wd, site_info, file_type):
    """Assemble the clean_data dict (as from_sharppy() returns it) from radiosonde arrays."""
    p = p * units.hPa
    z = z * units.m
    T = T * units.degC
    Td = Td * units.degC
    ws = (ws_ms * units('m/s')).to(units.kts)
    wd = wd * units.degrees

    u, v = mpcalc.wind_components(ws, wd)

    missing = [name for name, key in (("site ID", "site-id"), ("name", "site-name"), ("location", "site-lctn"),
                                      ("lat/lon", "site-latlon"), ("elevation", "site-elv"))
               if site_info.get(key) is None]
    print(f"**NOTE**: to plot your {file_type} file data with SounderPy, you'll have to manually provide additional site meta-data.\n"
          f"This process is fairly simple and is discussed in detail at https://kylejgillett.github.io/sounderpy/reference/customdatasources.html#ingesting-your-own-data-into-sounderpy.\n"
          f"Still needed: {', '.join(missing) if missing else 'nothing'} in site_info, and plot titles for the sounding plot.")

    return {
        "p": p,
        "z": z,
        "T": T,
        "Td": Td,
        "u": u,
        "v": v,
        "wd": wd,
        "ws": ws,
        "site_info": site_info,
    }


def _dewpoint_from_rh(T, rh):
    """Dewpoint (degC) from temperature (degC) and relative humidity (%); NaN where RH is missing or <= 0."""
    rh = np.where(rh > 0, rh, np.nan)
    return mpcalc.dewpoint_from_relative_humidity(T * units.degC, rh * units.percent).to('degC').m


def from_graw(filepath, date=None):
    """
    Read a GRAW radiosonde file and return a SounderPy-style clean_data dict.

    Based on the GRAW conversion script contributed by Alec Sczepanski (University of North Dakota).

    The file is the tab-separated table GRAW's software exports (Time, Time (UTC), P, T, Hu, Ws, Wd, Alt,
    ..., Lon, Lat, Dew, ...). Its first row is the ground station's surface observation at launch; the
    profile runs from there to the burst. Heights are above sea level.

    GRAW files hold the time of day but not the date, so give `date` ('YYYY-MM-DD', 'YYYYMMDD', or a
    datetime.date). Without it, the date is read from a file name that starts with YYYYMMDD
    (e.g. '20260125_2100.txt').
    """

    with open(filepath, "r", encoding="latin-1") as f:
        lines = [line.rstrip("\r\n") for line in f if line.strip()]

    if not lines:
        raise ValueError("The GRAW file is empty.")

    header = lines[0].split("\t")
    col = {name: _column_index(header, label, "GRAW") for name, label in
           (("time", "Time (UTC)"), ("p", "P "), ("T", "T "), ("rh", "Hu"), ("ws", "Ws"), ("wd", "Wd"),
            ("z", "Alt"), ("lon", "Lon"), ("lat", "Lat"), ("Td", "Dew"))}

    rows, times = [], []
    for line in lines[1:]:
        parts = line.split("\t")
        if len(parts) < len(header) - 1:
            continue
        rows.append([_to_float(parts[col[name]]) for name in ("p", "z", "T", "Td", "rh", "ws", "wd", "lat", "lon")])
        times.append(parts[col["time"]].strip())

    if not rows:
        raise ValueError("No sounding rows were found in the GRAW file.")

    # remember each row's time, then keep the ascent
    rows = np.column_stack([np.array(rows), np.arange(len(rows))])
    rows = _ascent_only(rows, 0, 1)
    p, z, T, Td, rh, ws, wd, lat, lon, idx = rows.T
    Td = np.where(np.isfinite(Td), Td, _dewpoint_from_rh(T, rh))

    # the launch date, from `date` or the file name
    if date is None:
        match = re.match(r"(\d{8})", os.path.basename(str(filepath)))
        if match is None:
            raise ValueError("GRAW files don't hold the launch date: give it with date='YYYY-MM-DD', "
                             "or name the file starting with YYYYMMDD (e.g. '20260125_2100.txt').")
        date = match.group(1)
    if isinstance(date, str):
        date = datetime.strptime(date.replace("-", ""), "%Y%m%d").date()

    # the launch time: the first time after the surface observation's placeholder (00:00:00)
    launch = next((times[int(i)] for i in idx if times[int(i)] not in ("", "00:00:00")), None)
    hour = int(launch.split(":")[0]) if launch else None

    site_info = {
        "site-id": None,
        "site-name": None,
        "site-lctn": None,
        "site-latlon": [round(float(lat[0]), 2), round(float(lon[0]), 2)] if np.isfinite([lat[0], lon[0]]).all() else None,
        "site-elv": int(round(float(z[0]))),
        "source": "GRAW RADIOSONDE",
        "valid-time": [date.year, date.month, date.day, hour] if hour is not None else None,
    }

    return _radiosonde_clean_data(p, z, T, Td, ws, wd, site_info, "GRAW")


def from_imet(filepath):
    """
    Read an iMet radiosonde file (from InterMet's SkySonde software) and return a SounderPy-style clean_data dict.

    Based on the iMet conversion script contributed by Alec Sczepanski (University of North Dakota).

    The file is SkySonde's flight data .csv: a block of station meta-data, then one row per second with
    the date, time, iMet PTU altitude (km above sea level), pressure, corrected temperature, humidity,
    ascent rate, and GPS winds. Rows before launch (ascent rate under 1 m/s) and after the burst are
    dropped. Dewpoint is calculated from the temperature and relative humidity.
    """

    with open(filepath, "r", encoding="utf-8", errors="replace", newline="") as f:
        lines = [row for row in csv.reader(f)]

    if not lines:
        raise ValueError("The iMet file is empty.")

    # station meta-data ('Station Name:, Grand Forks, ND') above the column header
    meta, header_idx = {}, None
    for idx, row in enumerate(lines):
        if row and row[0].strip().lower().startswith("date"):
            header_idx = idx
            break
        if len(row) >= 2 and row[0].strip().endswith(":"):
            meta[row[0].strip().rstrip(":").lower()] = ",".join(row[1:]).strip()

    if header_idx is None:
        raise ValueError("This does not appear to be an iMet (SkySonde) file.")

    header = [label.strip() for label in lines[header_idx]]
    names = (("z", "altitude (from iMet PTU)"), ("p", "iMet pressure"), ("T", "iMet air temperature (corrected)"),
             ("rh", "iMet humidity"), ("ascent", "iMet ascent rate"), ("ws", "GPS wind speed"),
             ("wd", "GPS wind direction"))
    col = {name: _column_index(header, label, "iMet (SkySonde)") for name, label in names}
    date_col, time_col = 0, _column_index(header, "time", "iMet (SkySonde)")

    rows, stamps = [], []
    for row in lines[header_idx + 1:]:
        if len(row) < len(header) - 1:
            continue
        rows.append([_to_float(row[col[name]], IMET_MISSING) for name, _ in names])
        stamps.append((row[date_col].strip(), row[time_col].strip()))

    if not rows:
        raise ValueError("No sounding rows were found in the iMet file.")

    rows = np.array(rows)
    rows[:, 0] = rows[:, 0] * 1000   # km to m

    # start at the launch: the first row rising at IMET_LAUNCH_ASCENT or faster
    launched = np.where(rows[:, 4] >= IMET_LAUNCH_ASCENT)[0]
    if launched.size == 0:
        raise ValueError("No launch was found in the iMet file (the ascent rate never reaches 1 m/s).")
    start = int(launched[0])
    rows, stamps = rows[start:], stamps[start:]

    rows = _ascent_only(rows, 1, 0)
    z, p, T, rh, _, ws, wd = rows.T
    Td = _dewpoint_from_rh(T, rh)

    try:
        launch = datetime.strptime(" ".join(stamps[0]), "%Y-%m-%d %H:%M:%S")
        valid_time = [launch.year, launch.month, launch.day, launch.hour]
    except ValueError:
        valid_time = None

    lat = _to_float(meta.get("station latitude"), IMET_MISSING)
    lon = _to_float(meta.get("station longitude"), IMET_MISSING)
    elevation = _to_float(str(meta.get("station altitude", "")).split("[")[0], IMET_MISSING)

    site_info = {
        "site-id": None,
        "site-name": meta.get("station name") or None,
        "site-lctn": None,
        "site-latlon": [round(lat, 2), round(lon, 2)] if np.isfinite([lat, lon]).all() else None,
        "site-elv": int(round(elevation)) if np.isfinite(elevation) else int(round(float(z[0]))),
        "source": "IMET RADIOSONDE",
        "valid-time": valid_time,
    }

    return _radiosonde_clean_data(p, z, T, Td, ws, wd, site_info, "iMet")


def from_windsond(filepath, elevation=None):
    """
    Read a Windsond radiosonde file and return a SounderPy-style clean_data dict.

    Based on the Windsond conversion script contributed by Alec Sczepanski (University of North Dakota).

    The file is the Windsond software's .sounding.csv: a comment line with the launch lat/lon and time,
    then height (m AGL), pressure, temperature, relative humidity, wind speed, and wind direction.
    Dewpoint is calculated from the temperature and relative humidity.

    Windsond heights are above ground level and the file holds no station elevation, so give the
    launch site's `elevation` (m) to put the heights above sea level like other SounderPy data. Without
    it, the heights stay above ground and site_info has no elevation.
    """

    with open(filepath, "r", encoding="utf-8", errors="replace", newline="") as f:
        lines = [row for row in csv.reader(f)]

    if not lines:
        raise ValueError("The Windsond file is empty.")

    # launch lat/lon/time from the comment line ('# Radiation correction v2.4. Params: lat=..., lon=..., utc_time=...')
    params = {}
    for row in lines:
        if row and row[0].lstrip().startswith("#"):
            for part in row:
                key, sep, value = part.partition("=")
                if sep:
                    params[key.strip().split()[-1].lower()] = value.strip()

    header_idx = next((idx for idx, row in enumerate(lines) if row and row[0].strip().lower().startswith("height")), None)
    if header_idx is None:
        raise ValueError("This does not appear to be a Windsond file.")

    header = [label.strip() for label in lines[header_idx]]
    names = (("z", "Height"), ("p", "Pressure"), ("T", "Temperature"), ("rh", "Relative humidity"),
             ("ws", "Wind speed"), ("wd", "Wind direction"))
    col = {name: _column_index(header, label, "Windsond") for name, label in names}

    rows = [[_to_float(row[col[name]]) for name, _ in names] for row in lines[header_idx + 1:]
            if len(row) >= len(header) and not row[0].lstrip().startswith("#")]

    if not rows:
        raise ValueError("No sounding rows were found in the Windsond file.")

    rows = _ascent_only(np.array(rows), 1, 0)
    z, p, T, rh, ws, wd = rows.T
    Td = _dewpoint_from_rh(T, rh)
    if elevation is not None:
        z = z + float(elevation)

    try:
        launch = datetime.strptime(params.get("utc_time", ""), "%Y-%m-%d %H:%M")
        valid_time = [launch.year, launch.month, launch.day, launch.hour]
    except ValueError:
        valid_time = None
    lat, lon = _to_float(params.get("lat")), _to_float(params.get("lon"))

    site_info = {
        "site-id": None,
        "site-name": None,
        "site-lctn": None,
        "site-latlon": [round(lat, 2), round(lon, 2)] if np.isfinite([lat, lon]).all() else None,
        "site-elv": int(round(float(elevation))) if elevation is not None else None,
        "source": "WINDSOND RADIOSONDE",
        "valid-time": valid_time,
    }

    return _radiosonde_clean_data(p, z, T, Td, ws, wd, site_info, "Windsond")



#########################################################################






#########################
# FILE WRITING FUNCTIONS
#########################################################################

def to_file(file_type, clean_data, filename=None, convert_to_AGL=True):
    '''
    Create a file of 'cleaned' SounderPy data

   :param file_type: a `str` representing the file type you'd like to export data to: 'csv', 'sharppy', or 'cm1'.
   :type file_type: str, required
   :param clean_data: 'cleaned' SounderPy data `dict`
   :type clean_data: dict, required
   :param filename: the name you'd like to give the file
   :type filename: str, required
   :return: a file of SounderPy data.
    '''

    file_type = str(file_type).casefold()
    if file_type not in ('cm1', 'csv', 'sharppy'):
        raise ValueError(f"'{file_type}' is not a valid file type. Valid file types are 'csv', 'sharppy', and 'cm1'")

    # set file name
    if filename is None:
        filename = f'sounderpy_data'
    else:
        filename = filename

    ####################################### CM1 #######################################
    if file_type == 'cm1':
        '''
        creates CM1 input sounding file for CM1 integration

        Derived from Kelton Halbert / Leigh Orf via github: 
        https://github.com/leighorf/LOFS-read/blob/master/bin/sndmod
        '''

        num_lines = len(clean_data['p'])
        delimiter = ''

        # use metpy to find parameters that CM1 likes, in a copy so the caller's data is unchanged
        clean_data = dict(clean_data)
        clean_data['theta'] = mpcalc.potential_temperature(clean_data['p'], clean_data['T'])
        clean_data['relhm'] = mpcalc.relative_humidity_from_dewpoint(clean_data['T'], clean_data['Td'])
        clean_data['mixrt'] = mpcalc.mixing_ratio_from_relative_humidity(clean_data['p'], clean_data['T'],
                                                                         clean_data['relhm']) * 1000

        # CM1 can't read missing values: the surface line needs all of its values,
        # and levels above it missing any value are left out
        values = {key: np.ma.filled(np.ma.asarray(clean_data[key].m, dtype=float), np.nan)
                  for key in ['p', 'z', 'theta', 'mixrt', 'u', 'v']}
        if not np.all(np.isfinite([values['p'][0], values['theta'][0], values['mixrt'][0]])):
            raise ValueError("The lowest level is missing a pressure, temperature, or dewpoint, so a CM1 surface line can't be written")
        valid = np.all([np.isfinite(values[key]) for key in ['z', 'theta', 'mixrt', 'u', 'v']], axis=0)
        skipped = [idx for idx in range(1, num_lines) if not valid[idx]]
        if skipped:
            skipped_z = [int(round(values['z'][idx] - values['z'][0])) if np.isfinite(values['z'][idx]) else 'unknown' for idx in skipped]
            warnings.warn(f"Left {len(skipped)} level(s) with missing values out of the CM1 file, at heights (m AGL): {skipped_z}", Warning)

        # create file
        outfile = open(filename, 'w')

        # create the sfc values line
        top_line = (
                "%12s" % str(format(np.around(clean_data["p"][0].m, 6), "0.6f")) + delimiter + "\t" +
                "%12s" % str(format(np.around(clean_data["theta"][0].m, 6), "0.6f")) + delimiter + "\t" +
                "%12s" % str(format(np.around(clean_data["mixrt"][0].m, 6), "0.6f")) + "\n"
        )

        # write the sfc values line to the file
        outfile.write(top_line)

        # add data to lines
        for idx in range(1, num_lines):
            if not valid[idx]:
                continue
            line_str = ""
            if convert_to_AGL:
                line_str += "%12s" % str(format(np.around(clean_data["z"][idx].m - clean_data["z"][0].m, 6), "0.6f")) + delimiter + str("\t")
            else:
                line_str += "%12s" % str(
                    format(np.around(clean_data["z"][idx].m, 6), "0.6f")) + delimiter + str("\t")
            line_str += "%12s" % str(format(np.around(clean_data["theta"][idx].m, 6), "0.6f")) + delimiter + str("\t")
            line_str += "%12s" % str(format(np.around(clean_data["mixrt"][idx].m, 6), "0.6f")) + delimiter + str("\t")
            line_str += "%12s" % str(format(np.around(clean_data["u"][idx].m / 1.94384, 6), "0.6f")) + delimiter + str(
                "\t")
            line_str += "%12s" % str(format(np.around(clean_data["v"][idx].m / 1.94384, 6), "0.6f")) + str("\n")
            outfile.write(line_str)

        outfile.close()



    ####################################### CSV #######################################
    elif file_type == 'csv':
        '''
        creates CSV file of sounding data
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
        creates NSHARP input sounding file for SharpPy integration

        Derived from Kelton Halbert / Leigh Orf via github: 
        https://github.com/leighorf/LOFS-read/blob/master/bin/sndmod
        '''

        outfile_file = open(filename, 'w')

        outfile_loc = ("****")

        # valid-time parts may be ints or strings, and the hour may be 'HH' or 'HH:MM'
        try:
            year, month, day, hour = [str(part) for part in clean_data['site_info']['valid-time'][0:4]]
            hour, _, minute = hour.partition(':')
            dt = datetime(int(year), int(month), int(day), int(hour), int(minute or 0))
        except (KeyError, TypeError, ValueError):
            # SHARPpy needs a date, so data without one (e.g. CM1) gets a placeholder
            warnings.warn("No usable valid time in clean_data['site_info']['valid-time']; writing 2000-01-01 00Z as a placeholder", Warning)
            dt = datetime(2000, 1, 1, 0)

        outfile_file.write("%TITLE%\n")
        outfile_file.write("%s   %s\n" % (clean_data['site_info']['site-id'], dt.strftime("%y%m%d/%H%M")))
        outfile_file.write("   LEVEL       HGHT       TEMP       DWPT       WDIR       WSPD\n")
        outfile_file.write("-------------------------------------------------------------------\n")
        outfile_file.write("%RAW%\n")

        ws = mpcalc.wind_speed(clean_data['u'], clean_data['v'])
        wd = mpcalc.wind_direction(clean_data['u'], clean_data['v'])

        # SHARPpy marks missing values with -9999, so fill masked and NaN values with it
        new_data = {}
        for col, values in zip(['p', 'z', 'T', 'Td', 'wd', 'ws'],
                               [clean_data['p'], clean_data['z'], clean_data['T'], clean_data['Td'], wd, ws]):
            new_data[col] = np.ma.masked_invalid(np.ma.asarray(values.m, dtype=float)).filled(-9999.)

        for idx in range(new_data['p'].shape[0]):
            string = ""
            for col in ['p', 'z', 'T', 'Td', 'wd', 'ws']:
                string += "%12.6f,  " % new_data[col][idx]

            outfile_file.write(string[:-3] + "\n")
        outfile_file.write("%END%\n")
        outfile_file.close()

    ##########################################################################################################################################






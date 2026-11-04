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

     Functions for writing files CM1/SHARPPY/CSV and reading SHARPPY files


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






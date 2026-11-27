"""
    SOUNDERPY VAD DATA

    Purpose of module:

    Get NEXRAD VAD Wind Profile (VWP) data and parse it into a SounderPy `vad_data` dict for
    build_vad_hodograph() and vad_params(). The VWP is the radar's own VAD product (NEXRAD Level III
    product 48, 'NVW'), issued every volume scan. Files come from Unidata's NEXRAD Level III archive
    on AWS (March 2020 on, about 10 KB each) and are read with MetPy.

        vad_data = fetch_vad('KDMX')                                  # latest
        vad_data = fetch_vad('KDMX', '2024', '05', '21', '20', '30')  # closest to 20:30 UTC

    The VWP starts a few hundred feet above the radar, so by default the nearest METAR wind is added as
    the surface (0 m AGL) point, or the lowest level of a BUFKIT model profile (RAP by default) when there
    is no METAR. When the VAD doesn't reach 6 km, the model's winds are added as 'sm_profile' for Bunkers
    storm motion.

    Docs:
        VAD hodographs:  https://sounderpy.app/docs/reference/plottingdata.html#building-vad-hodographs

    References:
        NOAA NEXRAD Level III on AWS (Unidata): https://registry.opendata.aws/noaa-nexrad/ (VWP files)
        Browning, K. A., and R. Wexler, 1968: The determination of kinematic properties of a wind field
            using Doppler radar. J. Appl. Meteor., 7, 105-113,
            https://doi.org/10.1175/1520-0450(1968)007<0105:TDOKPO>2.0.CO;2 (VAD wind profiles)
        May, R. M., S. C. Arms, P. Marsh, E. Bruning, J. R. Leeman, K. Goebbert, J. E. Thielen, Z. S.
            Bruick, and M. D. Camron, 2023: MetPy: A Python package for meteorological data. Unidata,
            https://doi.org/10.5065/D6WW7G29 (Level III reader, station names)
        Iowa Environmental Mesonet ASOS archive: https://mesonet.agron.iastate.edu/request/download.phtml
            (surface winds)
        Bunkers, M. J., B. A. Klimowski, J. W. Zeitler, R. L. Thompson, and M. L. Weisman, 2000: Predicting
            supercell motion using a new hodograph technique. Wea. Forecasting, 15, 61-79,
            https://doi.org/10.1175/1520-0434(2000)015<0061:PSMUAN>2.0.CO;2 (storm motion)

    (C) KYLE J GILLETT, UNIVERSITY OF NORTH DAKOTA, 2026
"""

### IMPORT SOFTWARE ###
#########################################################################################################
import csv
import io
import math
import re
import time
import urllib.request
import warnings
from datetime import datetime, timedelta, timezone
from xml.etree import ElementTree

import numpy as np

LEVEL3_BUCKET_URL = "https://unidata-nexrad-level3.s3.amazonaws.com"
FIRST_ARCHIVE_DAY = datetime(2020, 3, 31, tzinfo=timezone.utc)
MS_TO_KTS = 1.94384
FT_TO_M = 0.3048


#########################
# FIND VWP FILES
#########################################################################
def _level3_id(site):
    """The 3-letter Level III ID: 'KDMX' -> 'DMX', 'PHKI' -> 'HKI', 'TJUA' -> 'JUA'; 'DMX' stays 'DMX'."""
    site = str(site).strip().upper()
    if len(site) == 4:
        return site[1:]
    if len(site) == 3:
        return site
    raise ValueError(f"'{site}' is not a NEXRAD site ID; use the 4-letter ID, such as 'KDMX'")


def _list_keys(prefix, timeout=20):
    """Every key in the bucket that starts with prefix."""
    keys, token = [], None
    while True:
        url = f"{LEVEL3_BUCKET_URL}/?list-type=2&prefix={prefix}"
        if token:
            url += "&continuation-token=" + urllib.request.quote(token, safe="")
        with urllib.request.urlopen(url, timeout=timeout) as response:
            root = ElementTree.fromstring(response.read())
        ns = {"s3": root.tag.split("}")[0].strip("{")}
        keys += [k.text for k in root.findall("s3:Contents/s3:Key", ns)]
        token_el = root.find("s3:NextContinuationToken", ns)
        if token_el is None:
            return keys
        token = token_el.text


def _key_time(key):
    """'DMX_NVW_2026_11_26_00_04_04' -> its UTC datetime."""
    return datetime.strptime(key.split("_NVW_")[1], "%Y_%m_%d_%H_%M_%S").replace(tzinfo=timezone.utc)


def _find_keys(radar_id, target, max_offset_min):
    """VWP files closest to target (UTC datetime) first, or newest first when target is None, all within
    max_offset_min of the first one, so a scan without winds can fall back to the next."""
    window = timedelta(minutes=max_offset_min)
    if target is None:
        now = datetime.now(timezone.utc)
        keys = []
        for day in (now - timedelta(days=1), now):
            keys += _list_keys(f"{radar_id}_NVW_{day:%Y_%m_%d}")
        if not keys:
            raise ValueError(f"No VAD Wind Profile files from {radar_id} in the last day; "
                             "the radar may be down or the ID may be wrong")
        newest = max(_key_time(k) for k in keys)
        return sorted((k for k in keys if newest - _key_time(k) <= window), key=_key_time, reverse=True)

    if target < FIRST_ARCHIVE_DAY:
        raise ValueError(f"VAD Wind Profile files are available from {FIRST_ARCHIVE_DAY:%Y-%m-%d} on")

    # the day before or after can hold the closest file near midnight
    days = {target.date()}
    days |= {(target - window).date(), (target + window).date()}
    keys = []
    for day in sorted(days):
        keys += _list_keys(f"{radar_id}_NVW_{day:%Y_%m_%d}")
    if not keys:
        raise ValueError(f"No VAD Wind Profile files from {radar_id} on {target:%Y-%m-%d}; "
                         "the radar may have been down or the ID may be wrong")

    keys = sorted(keys, key=lambda k: abs(_key_time(k) - target))
    offset = abs(_key_time(keys[0]) - target)
    if offset > window:
        raise ValueError(f"The closest VAD Wind Profile from {radar_id} is {offset.total_seconds() / 60:.0f} minutes "
                         f"from {target:%Y-%m-%d %H:%M}Z ({_key_time(keys[0]):%H:%M}Z); the radar may have been down")
    return [k for k in keys if abs(_key_time(k) - target) <= window]


#########################
# PARSE THE VWP TABLE
#########################################################################
def _parse_table(tab_pages):
    """Rows of the 'VAD Algorithm Output' table as a dict of arrays: alt (100s of ft MSL), u, v (m/s), rms (kts)."""
    rows = []
    for page in tab_pages:
        lines = page.splitlines()
        header = next((line.split() for line in lines if line.split()[:3] == ["ALT", "U", "V"]), None)
        if header is None:
            continue
        cols = {name: header.index(name) for name in ("ALT", "U", "V", "RMS")}
        for line in lines:
            parts = line.split()
            if len(parts) != len(header) or not re.fullmatch(r"\d+", parts[0]):
                continue
            try:
                rows.append(tuple(float(parts[cols[c]]) for c in ("ALT", "U", "V", "RMS")))
            except ValueError:  # 'NA'
                continue
    if not rows:
        raise ValueError("This VAD Wind Profile has no winds (the radar may not have detected enough returns)")

    # the table mixes fixed heights with the slant-range VAD, which can repeat a height: keep the lower RMS
    best = {}
    for alt, u, v, rms in rows:
        if alt not in best or rms < best[alt][3]:
            best[alt] = (alt, u, v, rms)
    alt, u, v, rms = (np.array(col) for col in zip(*sorted(best.values())))
    return {"alt": alt, "u": u, "v": v, "rms": rms}


#########################
# SURFACE WIND
#########################################################################
def _surface_wind(radar_site, lat, lon, valid, max_km, window_min):
    """The nearest METAR wind (u, v in kts, plus the obs dict), or None."""
    from .sfc_obs import nearest_stations, IEM_ASOS_URL

    if valid > datetime.now(timezone.utc) + timedelta(minutes=10):
        return None
    # the station table also has the radar itself and sites without METARs (such as 'Q' IDs and offices), so
    # keep IDs that start like the radar's ('K' in the CONUS, 'P' in Alaska and Hawaii, 'T' in Puerto Rico)
    # and ask for the nearest several in one request, which also keeps clear of IEM's rate limit
    candidates = [c for c in nearest_stations(lat, lon, n=40, max_km=max_km)
                  if c[0] != radar_site and c[0][0] == radar_site[0]][:8]
    if not candidates:
        return None
    sts = (valid - timedelta(minutes=window_min)).strftime("%Y-%m-%dT%H:%MZ")
    ets = (valid + timedelta(minutes=window_min)).strftime("%Y-%m-%dT%H:%MZ")
    query = "&".join([f"station={c[0]}" for c in candidates] +
                     ["data=drct", "data=sknt", f"sts={sts}", f"ets={ets}", "tz=Etc/UTC", "format=onlycomma",
                      "latlon=no", "missing=M", "direct=no", "report_type=3", "report_type=4"])
    request = urllib.request.Request(f"{IEM_ASOS_URL}?{query}", headers={"User-Agent": "sounderpy"})
    try:
        with urllib.request.urlopen(request, timeout=15) as response:
            text = response.read().decode("utf-8", "replace")
    except Exception:  # network trouble or IEM rate limiting: no surface wind
        return None
    if not text.startswith("station"):
        return None

    # each station's report closest to the valid time
    reports = {}
    for row in csv.DictReader(io.StringIO(text)):
        obs_time = datetime.strptime(row["valid"], "%Y-%m-%d %H:%M").replace(tzinfo=timezone.utc)
        offset = abs((obs_time - valid).total_seconds())
        if row["station"] not in reports or offset < reports[row["station"]][0]:
            reports[row["station"]] = (offset, obs_time, row)

    # IEM drops the 'K' from CONUS IDs
    for station, name, km in candidates:
        found = reports.get(station) or reports.get(station[1:])
        if found is None:
            continue
        _, obs_time, row = found
        try:
            speed = float(row["sknt"])
        except ValueError:
            continue
        obs = {"time": obs_time, "sknt": speed, "drct": row["drct"]}
        try:
            obs["drct"] = float(row["drct"])
        except ValueError:
            pass
        if speed == 0:
            u, v = 0.0, 0.0
        elif isinstance(obs.get("drct"), float):
            rad = math.radians(obs["drct"])
            u, v = -speed * math.sin(rad), -speed * math.cos(rad)
        else:  # a variable wind has no direction
            continue
        obs.update(station=station, name=name, distance_km=km)
        return u, v, obs
    return None


#########################
# MODEL PROFILE (BUFKIT)
#########################################################################
_BUFKIT_SITES = None  # [(id, lat, lon, elevation m)], loaded once


def _bufkit_sites():
    """The BUFKIT station list (the same one get_bufkit_data() checks IDs against), loaded once."""
    global _BUFKIT_SITES
    if _BUFKIT_SITES is None:
        import pandas as pd
        table = pd.read_csv('https://raw.githubusercontent.com/kylejgillett/sounderpy/main/src/BUFKIT-STATIONS-MASTER.txt',
                            skiprows=7, skipinitialspace=True)
        sites = {}
        for sid, lat, lon, elev in zip(table['ID'], table['LAT'], table['LON'], table['EL(m)']):
            sid = str(sid).strip()
            # IDs with '#' (such as 'P#A') are model grid points without a place name for the figure
            if sid and '#' not in sid and sid not in sites:
                sites[sid] = (sid, float(lat), float(lon), float(elev))
        _BUFKIT_SITES = list(sites.values())
    return _BUFKIT_SITES


def _model_profile(model, lat, lon, valid, max_km=100.0, tries=2):
    """
    The wind profile from the BUFKIT site nearest (lat, lon), from the model run valid closest to `valid`:
    a dict of 'z' (m AGL, the lowest model level as the surface), 'u', 'v' (kts), and 'source' (such as
    'RAP OUN 23Z F00'), or None. Recent runs may not be archived yet, so up to 3 earlier runs are tried.
    """
    import contextlib
    from .bufkit_data import fetch_bufkit
    from .sfc_obs import _distance_km

    try:
        sites = _bufkit_sites()
    except Exception:
        return None
    ranked = sorted((_distance_km(lat, lon, s_lat, s_lon), sid, elev) for sid, s_lat, s_lon, elev in sites)
    valid_hour = (valid + timedelta(minutes=30)).replace(minute=0, second=0, microsecond=0)

    for km, sid, elev in ranked[:tries]:
        if km > max_km:
            break
        for fcst_hour in range(4):
            run = valid_hour - timedelta(hours=fcst_hour)
            try:
                with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    data = fetch_bufkit(model, sid, fcst_hour, f"{run:%Y}", f"{run:%m}", f"{run:%d}", f"{run:%H}",
                                        True, True)
            except Exception:
                continue
            z = np.asarray(data['z'].m, dtype=float)
            z = z - z[0]
            return {"z": z, "u": np.asarray(data['u'].m, dtype=float), "v": np.asarray(data['v'].m, dtype=float),
                    "source": f"{model.upper()} {sid} {run:%H}Z F{fcst_hour:02d}", "station": sid, "distance_km": km,
                    "model": model.upper(), "run": run, "fcst_hour": fcst_hour, "valid": valid_hour}
    return None


#########################
# FETCH VAD DATA
#########################################################################
def fetch_vad(site, year=None, month=None, day=None, hour=None, minute=0, sfc_wind=True, model='rap', max_rms=None,
              max_offset_min=30, hush=False):
    """
    Get the NEXRAD VAD Wind Profile closest to a time (or the latest, when no date is given) as a
    `vad_data` dict: 'z' (m AGL), 'u', 'v' (kts), 'rms' (kts, NaN for the surface wind), 'site_info', and,
    when the VAD doesn't reach 6 km, 'sm_profile' (model winds that vad_params() uses for Bunkers storm motion).

    sfc_wind: True for the nearest METAR wind, then the model's lowest level; a [direction, speed (kts)] list
    for your own surface wind; or False for none (the profile then starts at the lowest VAD height).
    model: the BUFKIT model ('rap', 'hrrr', 'gfs', ...) for the surface wind and storm motion fallbacks, or None.
    """
    st = time.time()
    if not hush:
        print(f'> VAD DATA ACCESS FUNCTION\n  -----------------------------------')

    radar_id = _level3_id(site)
    site_id = str(site).strip().upper() if len(str(site).strip()) == 4 else f"K{radar_id}"

    if year is None:
        target = None
    else:
        if None in (month, day, hour):
            raise ValueError("Give year, month, day, and hour (or none of them for the latest VAD)")
        target = datetime(int(year), int(month), int(day), int(hour), int(minute or 0), tzinfo=timezone.utc)

    custom_sfc = isinstance(sfc_wind, (list, tuple))
    if custom_sfc and len(sfc_wind) != 2:
        raise ValueError("sfc_wind as a list is [direction (degrees), speed (kts)], such as [180, 15]")

    from metpy.io import Level3File, station_info

    # a scan without winds (clear air, or the radar starting up) falls back to the next closest, up to 3
    keys = _find_keys(radar_id, target, max_offset_min)
    for n, key in enumerate(keys[:3]):
        with urllib.request.urlopen(f"{LEVEL3_BUCKET_URL}/{key}", timeout=30) as response:
            vwp = Level3File(io.BytesIO(response.read()))
        try:
            table = _parse_table(vwp.tab_pages)
            break
        except ValueError:
            if n == len(keys[:3]) - 1:
                raise ValueError(f"The {len(keys[:3])} VAD Wind Profiles from {radar_id} closest to the requested time have no "
                                 "winds; the radar may not be detecting enough returns (common in clear, dry air)")

    if max_rms is not None:
        keep = table["rms"] <= max_rms
        if not keep.any():
            raise ValueError(f"No VAD winds have an RMS error at or below {max_rms} kts")
        table = {k: v[keep] for k, v in table.items()}

    if len(table["alt"]) < 5:
        warnings.warn(f"This VAD Wind Profile has only {len(table['alt'])} height(s); the radar may not be detecting "
                      "enough returns (common in clear, dry air)", Warning)

    # heights are hundreds of feet above sea level; the file gives the radar's elevation in feet
    radar_elev_m = vwp.height * FT_TO_M
    z = table["alt"] * 100 * FT_TO_M - radar_elev_m
    u = table["u"] * MS_TO_KTS
    v = table["v"] * MS_TO_KTS
    rms = table["rms"]

    valid = vwp.metadata["vol_time"].replace(tzinfo=timezone.utc)
    lat, lon = round(float(vwp.lat), 3), round(float(vwp.lon), 3)
    station = station_info.get(site_id)
    site_name = station.name.upper() if station is not None else "NEXRAD"

    # the model profile is downloaded at most once, for the surface wind and/or the storm motion
    model_profile = None
    model_tried = False

    def get_model_profile():
        nonlocal model_profile, model_tried
        if not model_tried and model:
            model_tried = True
            model_profile = _model_profile(model, lat, lon, valid)
        return model_profile

    # SURFACE WIND: METAR, then the model's lowest level, or the user's own
    sfc = None
    if custom_sfc:
        direction, speed = float(sfc_wind[0]), float(sfc_wind[1])
        rad = math.radians(direction)
        sfc = {"source": "CUSTOM", "u": -speed * math.sin(rad), "v": -speed * math.cos(rad), "drct": direction, "sknt": speed}
    elif sfc_wind:
        found = _surface_wind(site_id, lat, lon, valid, max_km=50.0, window_min=60)
        if found is not None:
            u_sfc, v_sfc, obs = found
            sfc = {"source": "METAR", "u": u_sfc, "v": v_sfc, "station": obs["station"], "name": obs["name"],
                   "distance_km": obs["distance_km"], "time": obs["time"], "drct": obs["drct"], "sknt": obs["sknt"]}
        elif get_model_profile() is not None:
            mp = model_profile
            speed = float(np.hypot(mp["u"][0], mp["v"][0]))
            sfc = {"source": mp["source"], "u": float(mp["u"][0]), "v": float(mp["v"][0]), "station": mp["station"],
                   "distance_km": mp["distance_km"], "drct": round(float(np.degrees(np.arctan2(-mp["u"][0], -mp["v"][0])) % 360)),
                   "sknt": round(speed, 1), "model": mp["model"], "run": mp["run"], "fcst_hour": mp["fcst_hour"],
                   "time": mp["valid"]}
    if sfc is not None:
        z, u, v = np.r_[0.0, z], np.r_[sfc["u"], u], np.r_[sfc["v"], v]
        rms = np.r_[np.nan, rms]
    elif not hush:
        print('    > no surface wind; the profile starts at the lowest VAD height')

    # STORM MOTION: Bunkers needs winds to 6 km; below that, the model's winds stand in for it
    sm_profile = None
    if z[-1] < 6000 and get_model_profile() is not None and model_profile["z"][-1] >= 6000:
        sm_profile = {k: model_profile[k] for k in ("z", "u", "v", "source", "station", "model", "run", "fcst_hour", "valid")}

    data_source = "NOAA NEXRAD LEVEL III VWP"
    if sfc is not None:
        if sfc["source"] == "METAR":
            data_source += f" · SFC WIND {sfc['station']} {sfc['time']:%H%M}Z"
        else:
            data_source += f" · SFC WIND {'USER DEFINED' if sfc['source'] == 'CUSTOM' else sfc['source']}"
    if sm_profile is not None:
        data_source += f" · BUNKERS {sm_profile['source']}"

    vad_data = {
        "z": np.round(z, 1),
        "u": np.round(u, 2),
        "v": np.round(v, 2),
        "rms": rms,
        "site_info": {
            "source": "NEXRAD VAD WIND PROFILE",
            "data-source": data_source,
            "site-id": site_id,
            "site-name": site_name,
            "site-lctn": "",
            "site-latlon": [lat, lon],
            "site-elv": round(radar_elev_m),
            "valid-time": [f"{valid:%Y}", f"{valid:%m}", f"{valid:%d}", f"{valid:%H:%M}"],
            "sfc-wind": None if sfc is None else {k: val for k, val in sfc.items() if k not in ("u", "v")},
        },
    }
    if sm_profile is not None:
        vad_data["sm_profile"] = sm_profile

    if not hush:
        sfc_note = "" if sfc is None else f" + surface wind from {sfc.get('station', 'CUSTOM') if sfc['source'] in ('METAR', 'CUSTOM') else sfc['source']}"
        print(f"    > PROFILE FOUND: {site_id} VWP at {valid:%H:%M}Z on {valid:%m/%d/%Y}{sfc_note}")
        print(f"    > {len(table['alt'])} VAD heights, {z[-1] if len(z) else 0:.0f} m AGL top")
        if sm_profile is not None:
            print(f"    > the VAD doesn't reach 6 km, so Bunkers storm motion will come from {sm_profile['source']}")
        elif z[-1] < 6000:
            print("    > the VAD doesn't reach 6 km and no model profile was found, so Bunkers storm motion can't be computed; "
                  "give a storm motion, such as storm_motion=[240, 30]")
        print('    > COMPLETE --------')
        print('    > RUNTIME:', time.strftime("%H:%M:%S", time.gmtime(time.time() - st)))
    return vad_data

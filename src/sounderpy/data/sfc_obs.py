"""
    SOUNDERPY SURFACE OBSERVATION

    Purpose of module:

    Finds the ASOS/AWOS stations nearest a sounding and fetches the METAR closest to its valid time
    from the Iowa Environmental Mesonet archive (no API key). The winter sounding draws it as a
    station model in the hodograph's upper-right corner (build_sounding(..., mode='winter', show_obs=True)).

        obs = sounding_obs(clean_data)          # None if nothing within reach
        obs["station"], obs["time"], obs["tmpf"], obs["wxcodes"], obs["timing"] ...

    Station coordinates come from MetPy's built-in station table (downloaded once to ~/.cache/metpy,
    then read locally). IEM rate-limits rapid requests, so callers should treat a failure as "no obs".

    (C) KYLE J GILLETT, 2026
"""
import csv
import io
import math
import time
import urllib.request
from datetime import datetime, timedelta, timezone

IEM_ASOS_URL = "https://mesonet.agron.iastate.edu/cgi-bin/request/asos.py"
FIELDS = ("tmpf", "dwpf", "drct", "sknt", "gust", "vsby", "wxcodes", "mslp", "skyc1", "skyc2", "skyc3")

_STATIONS = None  # [(id, name, lat, lon)], loaded once


def _station_table():
    """4-letter ICAO stations (ASOS/AWOS sites) from MetPy's station table, loaded once."""
    global _STATIONS
    if _STATIONS is None:
        from metpy.io import station_info
        seen, stations = set(), []
        for sid, st in station_info.items():
            if len(sid) == 4 and sid.isalpha() and sid not in seen and st.latitude == st.latitude:
                seen.add(sid)
                stations.append((sid, st.name, float(st.latitude), float(st.longitude)))
        _STATIONS = stations
    return _STATIONS


def _distance_km(lat1, lon1, lat2, lon2):
    p1, p2 = math.radians(lat1), math.radians(lat2)
    dlat, dlon = p2 - p1, math.radians(lon2 - lon1)
    a = math.sin(dlat / 2) ** 2 + math.cos(p1) * math.cos(p2) * math.sin(dlon / 2) ** 2
    return 6371.0 * 2 * math.asin(math.sqrt(a))


def nearest_stations(lat, lon, n=5, max_km=75.0):
    """The n stations closest to (lat, lon) within max_km, nearest first: [(id, name, km)]."""
    ranked = sorted((_distance_km(lat, lon, s_lat, s_lon), sid, name) for sid, name, s_lat, s_lon in _station_table())
    return [(sid, name, km) for km, sid, name in ranked[:n] if km <= max_km]


def _value(text):
    if text in ("M", "", None):
        return None
    try:
        return float(text)
    except ValueError:
        return text


def fetch_obs(station, valid, window_min=60, timeout=10):
    """The report from `station` closest to `valid` (UTC datetime) within +/- window_min, or None."""
    sts = (valid - timedelta(minutes=window_min)).strftime("%Y-%m-%dT%H:%MZ")
    ets = (valid + timedelta(minutes=window_min)).strftime("%Y-%m-%dT%H:%MZ")
    query = "&".join([f"station={station}"] + [f"data={f}" for f in FIELDS] +
                     [f"sts={sts}", f"ets={ets}", "tz=Etc/UTC", "format=onlycomma", "latlon=no", "missing=M",
                      "trace=T", "direct=no", "report_type=3", "report_type=4"])
    request = urllib.request.Request(f"{IEM_ASOS_URL}?{query}", headers={"User-Agent": "sounderpy"})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        text = response.read().decode("utf-8", "replace")
    if not text.startswith("station"):  # e.g. "Too many requests from your IP address"
        raise RuntimeError(text.strip()[:120])

    best = None
    for row in csv.DictReader(io.StringIO(text)):
        obs_time = datetime.strptime(row["valid"], "%Y-%m-%d %H:%M").replace(tzinfo=timezone.utc)
        offset = abs((obs_time - valid).total_seconds())
        if best is None or offset < best[0]:
            best = (offset, obs_time, row)
    if best is None:
        return None
    offset, obs_time, row = best
    obs = {f: _value(row.get(f)) for f in FIELDS}
    obs.update(time=obs_time, offset_min=offset / 60)
    return obs


def sounding_valid_time(clean_data):
    """The valid time as a UTC datetime. BUFKIT and RAOB data give the hour ("14"); ACARS gives hour and
    minute ("14:18")."""
    year, month, day, hour = clean_data["site_info"]["valid-time"][:4]
    hour, _, minute = str(hour).partition(":")
    return datetime(int(year), int(month), int(day), int(hour), int(minute or 0), tzinfo=timezone.utc)


def sounding_obs(clean_data, max_km=75.0, tries=3, window_min=60):
    """
    The observation nearest the sounding in space and time: tries the closest `tries` stations within
    max_km until one has a report within +/- window_min of the valid time. Returns None when nothing is
    found (no nearby station, a future valid time, or the archive refusing the request).
    obs["timing"] has the seconds spent finding the station and fetching the report.
    """
    timing = {}
    t0 = time.perf_counter()
    lat, lon = (float(v) for v in clean_data["site_info"]["site-latlon"])
    candidates = nearest_stations(lat, lon, n=tries, max_km=max_km)
    timing["find_station"] = time.perf_counter() - t0

    valid = sounding_valid_time(clean_data)
    if valid > datetime.now(timezone.utc) + timedelta(minutes=10):
        return None  # a forecast valid in the future has no observation yet

    t0 = time.perf_counter()
    for station, name, km in candidates:
        obs = None
        for attempt in range(2):
            try:
                obs = fetch_obs(station, valid, window_min=window_min)
                break
            except Exception as ex:  # network trouble or rate limiting: treat as no obs
                timing["error"] = str(ex)[:120]
                # IEM answers rapid requests with "Too many requests"; wait and retry once
                if attempt == 0 and "too many" in str(ex).lower():
                    time.sleep(2.0)
                    continue
                break
        if obs is not None:
            timing["fetch"] = time.perf_counter() - t0
            obs.update(station=station, name=name, distance_km=km, timing=timing)
            return obs
    timing["fetch"] = time.perf_counter() - t0
    return None

from __future__ import annotations

import numpy as np
import pandas as pd
import xarray as xr
import metpy.calc as mpcalc

from metpy.units import units
from scipy.interpolate import interp1d



"""
    SOUNDERPY RAP/RUC GRIB PARSER

    The purpose of this module is to parse the *annoying* grib
    files for RAP/RUC data -- now the only option from NOAA through
    AWS and NCEI. 

    * the development of this script was aided by artificial intelligence 

    (C) KYLE J GILLETT, UNIVERSITY OF NORTH DAKOTA, 2026
"""




DZ = 250.0
MS_TO_KT = 1.94384


###############
# OPEN FILE
#########################################################################

def _open_grib(path, filter_by_keys):
    return xr.open_dataset(path, engine="cfgrib", backend_kwargs={"filter_by_keys": filter_by_keys, "errors": "ignore", "indexpath": "",})


def open_rapruc_grib(pressure_file, native_file=None):
    
    native_file = native_file or pressure_file

    ds_vertical = _open_grib(pressure_file, {"typeOfLevel": "isobaricInhPa"})

    ds_2m = _open_grib(native_file, {"typeOfLevel": "heightAboveGround", "level": 2},)

    ds_10m = _open_grib(native_file, {"typeOfLevel": "heightAboveGround", "level": 10,})

    try:
        ds_sfc = _open_grib(native_file, {"typeOfLevel": "surface","stepType": "instant",})
    except Exception:
        ds_sfc = _open_grib(native_file, {"typeOfLevel": "surface"})

    return ds_vertical, ds_2m, ds_10m, ds_sfc







###############
# GRIB UTIL FUNCTIONS
#########################################################################

def _get_var(ds, *names, required=True):
    for name in names:
        if name in ds.variables:
            return ds[name]

    if required:
        raise KeyError(
            f"None of {names} found. "
            f"Available variables: {list(ds.data_vars)}")

    return None





def _normalize_lon(lon):
    return (np.asarray(lon, dtype=float) + 180.0) % 360.0 - 180.0





def _latlon_arrays(ds):
    lat = np.asarray(ds["latitude"].values, dtype=float)
    lon = _normalize_lon(ds["longitude"].values)

    if lat.ndim == 1 and lon.ndim == 1:
        lon, lat = np.meshgrid(lon, lat)

    return lat, lon





def _box_mask(ds, lat0, lon0, box_avg_size):
    # define domain for box avg profile
    lat, lon = _latlon_arrays(ds)
    half = box_avg_size / 2.0

    mask = (
        (lat >= lat0 - half)
        & (lat <= lat0 + half)
        & (lon >= lon0 - half)
        & (lon <= lon0 + half))

    if not np.any(mask):
        dlat = lat - lat0
        dlon = (lon - lon0) * np.cos(np.deg2rad(lat0))

        iy, ix = np.unravel_index(np.nanargmin(dlat**2 + dlon**2), lat.shape)

        mask = np.zeros_like(lat, dtype=bool)
        mask[iy, ix] = True

    return mask





def _mean_2d(da, ds, lat, lon, box_avg_size):
    # perform mean on box-sounding domain
    arr = np.squeeze(np.asarray(da.values, dtype=float))

    if arr.ndim != 2:
        raise ValueError(
            f"{da.name}: expected 2-D field, got shape {arr.shape}")

    mask = _box_mask(ds, lat, lon, box_avg_size)

    return float(np.nanmean(arr[mask]))





def _mean_profile(da, ds, lat, lon, box_avg_size):
    p = np.asarray(ds["isobaricInhPa"].values, dtype=float)
    arr = np.asarray(da.values, dtype=float)

    axis = da.get_axis_num("isobaricInhPa")
    arr = np.moveaxis(arr, axis, 0)
    arr = np.squeeze(arr)

    if arr.ndim != 3:
        raise ValueError(
            f"{da.name}: expected (pressure,y,x), got shape {arr.shape}")

    mask = _box_mask(ds, lat, lon, box_avg_size)
    vals = np.array([np.nanmean(level[mask]) for level in arr], dtype=float,)

    return p, vals





def _dataset_time(ds):
    # extract time info from grib
    for name in ("valid_time", "time"):
        if name in ds.coords:
            try:
                return pd.Timestamp(np.asarray(ds[name].values).squeeze())
            except Exception:
                pass

    return pd.NaT




###############
# PROFILE CONSTRUCTION
#########################################################################

def build_clean_dict(ds_vertical, ds_2m, ds_10m, ds_sfc, latlon, box_avg_size=0.25, model="RAP", valid_time=None):
    # build sounderpy clean_data dict
    
    lat = float(latlon[0])
    lon = float(_normalize_lon(latlon[1]))

    # PRESSURE LEVEL FIELDS
    t = _get_var(ds_vertical, "t")
    gh = _get_var(ds_vertical, "gh", "z")
    rh = _get_var(ds_vertical, "r", "rh")
    u = _get_var(ds_vertical, "u")
    v = _get_var(ds_vertical, "v")
    omega = _get_var(ds_vertical, "w", "vvel", required=False)

    p, T = _mean_profile(t, ds_vertical, lat, lon, box_avg_size)
    p2, Z = _mean_profile(gh, ds_vertical, lat, lon, box_avg_size)
    p3, RH = _mean_profile(rh, ds_vertical, lat, lon, box_avg_size)
    p4, U = _mean_profile(u, ds_vertical, lat, lon, box_avg_size)
    p5, V = _mean_profile(v, ds_vertical, lat, lon, box_avg_size)

    if omega is not None:
        p6, OMG = _mean_profile(omega, ds_vertical, lat, lon, box_avg_size)
    else:
        p6 = p.copy()
        OMG = np.full_like(p, np.nan, dtype=float)

    for other in (p2, p3, p4, p5, p6):
        if len(other) != len(p) or not np.allclose(other, p):
            raise RuntimeError(
                "RAP/RUC pressure coordinates do not match.")

    T = T - 273.15
    U = U * MS_TO_KT
    V = V * MS_TO_KT

    Td = mpcalc.dewpoint_from_relative_humidity(
        T * units.degC,
        RH * units.percent,
    ).to("degC").m


    # SURFACE, 2m, 10m
    sp = _get_var(ds_sfc, "sp", "pres")
    terrain = _get_var(ds_sfc,"orog","gh","hgt",required=False)

    sfc_p = (_mean_2d(sp, ds_sfc, lat, lon, box_avg_size)/100.0)

    if terrain is not None:
        sfc_z = _mean_2d(terrain, ds_sfc, lat, lon, box_avg_size)
    else:
        order = np.argsort(p)
        sfc_z = float(np.interp(sfc_p, p[order], Z[order]))

    t2 = _get_var(ds_2m,"t2m","t",)
    sfc_T = (_mean_2d(t2, ds_2m, lat, lon, box_avg_size)- 273.15)


    q2 = _get_var( ds_2m, "sh2", "q", "q2", required=False)
    d2 = _get_var( ds_2m, "d2m", "dpt", "td", required=False)
    r2 = _get_var(ds_2m, "r2", "r", "rh", required=False)

    if q2 is not None:
        sfc_q = _mean_2d(q2, ds_2m, lat, lon, box_avg_size)

        sfc_Td = (mpcalc.dewpoint_from_specific_humidity(sfc_p * units.hPa, sfc_T * units.degC, sfc_q).to("degC").m)

    elif d2 is not None:
        sfc_Td = (_mean_2d(d2, ds_2m, lat, lon, box_avg_size)- 273.15)

    elif r2 is not None:
        sfc_rh_native = _mean_2d(r2, ds_2m, lat, lon, box_avg_size)

        sfc_Td = (mpcalc.dewpoint_from_relative_humidity(sfc_T * units.degC, sfc_rh_native * units.percent).to("degC").m)

    else:
        raise KeyError("No usable 2-m moisture variable found.")

    sfc_rh = (mpcalc.relative_humidity_from_dewpoint(sfc_T * units.degC, sfc_Td * units.degC).to("percent").m)

    u10 = _get_var( ds_10m, "u10", "u")
    v10 = _get_var( ds_10m, "v10", "v")

    sfc_u = (_mean_2d(u10, ds_10m, lat, lon, box_avg_size)* MS_TO_KT)

    sfc_v = (_mean_2d( v10, ds_10m, lat, lon, box_avg_size)* MS_TO_KT)


    # COMBINE FIELDS AND APPLY SFC MASKING
    keep = (
        np.isfinite(p)
        & np.isfinite(Z)
        & np.isfinite(T)
        & np.isfinite(Td)
        & np.isfinite(RH)
        & np.isfinite(U)
        & np.isfinite(V)
        & (p < sfc_p)
        & (Z > sfc_z))

    z_msl = np.concatenate(([sfc_z], Z[keep]))

    native = {
        "T": np.concatenate(
            ([sfc_T], T[keep])),
        "Td": np.concatenate(
            ([sfc_Td], Td[keep])),
        "rh": np.concatenate(
            ([sfc_rh], RH[keep])),
        "u": np.concatenate(
            ([sfc_u], U[keep])),
        "v": np.concatenate(
            ([sfc_v], V[keep])),
        "p": np.concatenate(
            ([sfc_p], p[keep])),
        "omega": np.concatenate(
            ([0.0], OMG[keep]))}

    order = np.argsort(z_msl)
    z_msl = z_msl[order]

    for key in native:
        native[key] = native[key][order]

    _, unique = np.unique(z_msl, return_index=True)
    unique = np.sort(unique)

    z_msl = z_msl[unique]

    for key in native:
        native[key] = native[key][unique]

    z_agl = z_msl - sfc_z
    z_agl[0] = 0.0


    # Z-LEVEL INTERPOLATION
    top = (np.floor(z_agl[-1] / DZ)* DZ)

    levels = np.arange(0.0, top + DZ, DZ)

    out = {"z": levels}

    for key in ("T","Td","rh","u","v","p","omega"):
        arr = native[key]

        valid = (np.isfinite(z_agl) & np.isfinite(arr))

        if np.count_nonzero(valid) < 2:
            out[key] = np.full(levels.shape, np.nan)
            continue
 
        f = interp1d(z_agl[valid], arr[valid], kind="linear", bounds_error=False, fill_value=np.nan)

        out[key] = np.asarray(f(levels), dtype=float)

        out[key][0] = arr[0]


    # CLEAN DATA DICT PROFILE METADATA 
    # Requested retrieval time is authoritative. Fall back to GRIB metadata if parser is called independently.
    if valid_time is not None:
        dt = pd.Timestamp(valid_time)
    else:
        dt = _dataset_time(ds_vertical)
        if pd.isna(dt):
            dt = _dataset_time(ds_sfc)

    if pd.isna(dt):
        raise ValueError(
            "Could not determine RAP/RUC valid time. "
            "Pass valid_time= when calling parse_rapruc_grib()."
        )

    time_list = [f"{dt.year:04d}", f"{dt.month:02d}", f"{dt.day:02d}", f"{dt.hour:02d}",]

    clean_dict = {
        "T": out["T"] * units.degC,
        "Td": out["Td"] * units.degC,
        "rh": out["rh"] * units.percent,
        "u": out["u"] * units.kt,
        "v": out["v"] * units.kt,
        "z": out["z"] * units.m,
        "p": out["p"] * units.hPa,
        "omega": out["omega"] * units("Pa/sec"),
    }

    clean_dict["site_info"] = {
        "site-id": "no-site-id",
        "site-name": "no-site-name",
        "site-lctn": "no-site-location",
        "site-latlon":[lat, lon],
        "site-elv": float(sfc_z),
        "source": "MODEL REANALYSIS",
        "model": model.upper(),
        "fcst-hour": "F00",
        "run-time": time_list,
        "valid-time": time_list,
        "box_area": (f"{box_avg_size}° x {box_avg_size}° BOX AVG")}

    clean_dict["titles"] = {"top_title": f"MODEL REANALYSIS VERTICAL PROFILE | {time_list[3]}Z {model.upper()} F00",
                            "left_title": f"{time_list[3]}Z {model.upper()} F00 | VALID: {time_list[1]}/{time_list[2]}/{time_list[0]} {time_list[3]}Z",
                            "right_title": f"{lat}, {lon} | {box_avg_size}° x {box_avg_size}° BOX AVG    "}

    return clean_dict


def parse_rapruc_grib(pressure_file, latlon, *, native_file=None, box_avg_size=0.25, model="RAP", valid_time=None,):
    """
    Parse downloaded RAP/RUC GRIB file(s) into a SounderPy clean_dict.
    """
    datasets = open_rapruc_grib(pressure_file, native_file=native_file)

    try:
        return build_clean_dict(*datasets, latlon=latlon, box_avg_size=box_avg_size, model=model, valid_time=valid_time)
    finally:
        for ds in datasets:
            ds.close()
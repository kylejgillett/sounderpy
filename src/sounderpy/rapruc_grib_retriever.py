from __future__ import annotations

import tempfile
from pathlib import Path

import fsspec
import pandas as pd
import requests
import time
import warnings

from .rapruc_grib_parser import parse_rapruc_grib
from .calc import *


"""
    SOUNDERPY NEW RAP-RUC DATA RETRIEVAL - from NCEI and AWS   

    Purpose of module: 

    A new way to download and plot RAP/RUC data as previous option using .nc
    has been essentially discontinued by NCEI. This version searches AWS and 
    existsing NCEI *grib* datasets for the requested date.
    
    Part of this script is based on the AWS method used by Chase Archive.


    (C) KYLE J GILLETT, UNIVERSITY OF NORTH DAKOTA, 2026
"""



AWS_BUCKET = "noaa-rap-pds"
NCEI_FILESERVER = "https://www.ncei.noaa.gov/thredds/fileServer"
NCEI_HEALTH_URL = "https://www.ncei.noaa.gov/thredds/catalog/model/model.html"



###############
# DATE CONSTRAINTS FOR DATA ACCESS
#########################################################################
'''
- start search for data in NCEI RAP before 2021, AWS after 2021, RUC before 2012
- can be overrid by directly searching in a specific dataset using `dataset` kwarg.
'''
RAP_START = pd.Timestamp("2012-05-01")
AWS_PREFERRED_START = pd.Timestamp("2021-05-01")
HEALTH_TIMEOUT = (2, 4)      # NCEI server check timeout
PROBE_TIMEOUT = (2, 4)       # check dataset url for grib file timeout
DOWNLOAD_TIMEOUT = (10, 120) # download file from server timeout

class NCEIServerUnavailable(ConnectionError):
    # raise error when NCEI is down and timeout is met.
    pass





###############
# AWS GRIB RAP DATA ACCESS
#########################################################################
def _aws_uri(dt, product):
    return (f"s3://{AWS_BUCKET}/rap.{dt:%Y%m%d}/"
            f"rap.t{dt:%H}z.{product}.grib2")



def _download_aws_file(uri, destination):
    # download aws grib file to local
    fs = fsspec.filesystem("s3",anon=True,)

    if not fs.exists(uri):
        raise FileNotFoundError(uri)

    fs.get(uri, str(destination))

    return str(destination)



def try_aws(dt, tempdir, resolution=13):
    """
    Try RAP AWS.

    AWS uses separate pressure-level and native-level files.
    """
    if resolution == 13:
        pressure_product = "awp130pgrbf00"
        native_product = "awp130bgrbf00"

    elif resolution == 20:
        pressure_product = "awp252pgrbf00"
        native_product = "awp252bgrbf00"

    pressure_uri = _aws_uri(dt, pressure_product) # get pressure level data
    native_uri = _aws_uri(dt,native_product) # get sfc level data

    print(f"    > searching AWS RAP {resolution}-km dataset")

    pressure_local = (
        Path(tempdir)
        / f"rap_{resolution}km_pressure.grib2"
    )
    native_local = (
        Path(tempdir)
        / f"rap_{resolution}km_native.grib2"
    )

    try:
        _download_aws_file(pressure_uri, pressure_local)
        _download_aws_file(native_uri, native_local)

    except FileNotFoundError:
        # remove a partial file if anything downloaded
        pressure_local.unlink(missing_ok=True)
        native_local.unlink(missing_ok=True)
        return None

    return {"source": f"AWS RAP {resolution}-km",
            "model": "RAP",
            "pressure_file": str(pressure_local),
            "native_file": str(native_local)}




###############
# NCEI GRIB RAP DATA ACCESS
#########################################################################

def _ncei_candidates(dt, model=None):
    """
    based on date info provided by user, search for data in dirs that 
    make sense -- e.g.: skip RAP datasets if date is 2009 and will thus be in RUC set
    """

    ymd = dt.strftime("%Y%m%d")
    hh = dt.strftime("%H")


    # construct possible url components for searching
    rap130 = f"rap_130_{ymd}_{hh}00_000.grb2"
    rap252 = f"rap_252_{ymd}_{hh}00_000.grb2"
    ruc130 = f"ruc2anl_130_{ymd}_{hh}00_000.grb2"
    ruc252 = f"ruc2anl_252_{ymd}_{hh}00_000.grb"

    rap = [
        # dataset key        model  ncei url key         date url component
        ("RAP_13km_anl",     "RAP", "model-rap130anl",     rap130),
        ("RAP_13km_anl_old", "RAP", "model-rap130anl-old", rap130),
        ("RAP_13km",         "RAP", "model-rap130",        rap130),
        ("RAP_13km_old",     "RAP", "model-rap130-old",    rap130),
        ("RAP_25km_anl",     "RAP", "model-rap252anl",     rap252),
        ("RAP_25km_anl_old", "RAP", "model-rap252anl-old", rap252),
        ("RAP_25km",         "RAP", "model-rap252",        rap252),
        ("RAP_25km_old",     "RAP", "model-rap252-old",    rap252),
    ]

    ruc = [
        ("RUC_13km",         "RUC", "model-ruc130anl",     ruc130),
        ("RUC_13km_old",     "RUC", "model-ruc130anl-old", ruc130),
        ("RUC_25km",         "RUC", "model-ruc252anl",     ruc252),
        ("RUC_25km_old",     "RUC", "model-ruc252anl-old", ruc252),
    ]

    if model is not None:
        model = str(model).upper()
        if model == "RAP":
            return rap
        if model == "RUC":
            return ruc
        raise ValueError("model must be 'RAP' or 'RUC'")

    return rap if dt >= RAP_START else ruc




def _ncei_url(dt, collection, filename):
    # construct ncei url 
    return (
        f"{NCEI_FILESERVER}/"
        f"{collection}/"
        f"{dt:%Y%m}/"
        f"{dt:%Y%m%d}/"
        f"{filename}"
    )




def _check_ncei_server():
    """
    Check NCEI THREDDS once before searching individual archive files.

    This distinguishes a server outage/network failure from an ordinary
    "file not found" result and avoids waiting through every candidate URL
    when NCEI itself is unavailable.
    """
    print("    > checking NCEI server availability first...")

    try:
        with requests.get(
            NCEI_HEALTH_URL,
            stream=True,
            timeout=HEALTH_TIMEOUT,
        ) as r:
            # Any normal HTTP response below 500 means the host responded.
            # 4xx responses can reflect a path/access issue, but the server
            # itself is still reachable.
            if r.status_code >= 500:
                raise NCEIServerUnavailable(
                    f"NCEI THREDDS is currently returning HTTP "
                    f"{r.status_code}. Please try again later."
                )

    except requests.Timeout as exc:
        raise NCEIServerUnavailable(
            "NCEI THREDDS could not be reached within the timeout window. "
            "The server may be temporarily unavailable. Please try again later."
        ) from exc

    except requests.ConnectionError as exc:
        raise NCEIServerUnavailable(
            "NCEI THREDDS could not be reached. The server may be temporarily "
            "unavailable, or the network connection to NCEI failed. "
            "Please try again later."
        ) from exc

    except requests.RequestException as exc:
        raise NCEIServerUnavailable(
            f"NCEI THREDDS availability check failed: {exc}"
        ) from exc

    print("    > NCEI server test was successful")




def _url_exists(url):
    # test ncei urls to see if the requested date exists within the given dataset
    try:
        with requests.get(url,headers={"Range": "bytes=0-0"}, stream=True, timeout=PROBE_TIMEOUT) as r: 
            return r.status_code in (200, 206)
    except requests.RequestException:
        return False




def _download_http(url, destination):
    # begin downloading data if search was successful
    print(f"    > downloading: {url}")

    with requests.get(url, stream=True, timeout=DOWNLOAD_TIMEOUT) as r:
        r.raise_for_status()

        with open(destination, "wb") as f:
            for chunk in r.iter_content(
                chunk_size=1024 * 1024):
                if chunk:
                    f.write(chunk)
    return str(destination)




def try_ncei(dt,tempdir,dataset=None):
    # perform server check
    _check_ncei_server()

    # search dataset if provided by user, else search for data
    if dataset is not None:
        all_candidates = (_ncei_candidates(dt, model="RAP") + _ncei_candidates(dt, model="RUC"))
        candidates = [item for item in all_candidates if item[0].lower() == dataset.lower()]

        if not candidates:
            valid = [item[0] for item in all_candidates]

            raise ValueError(
                # raise error if dataset not found
                f"Unknown dataset '{dataset}'. "
                f"Valid choices: {valid}")
        
    else:
        candidates = _ncei_candidates(dt)

    # search for data by date constraints 
    era = "RAP" if dt >= RAP_START else "RUC"
    print(
        f"    > searching NCEI {era} archives "
        f"({len(candidates)} likely collections)")

    for (name, model, collection, filename) in candidates:

        url = _ncei_url(dt, collection, filename)

        print(f"      trying {name}")

        if not _url_exists(url):
            continue

        # if the data is found, grab it
        suffix = (".grb2" if filename.endswith(".grb2") else ".grb")
        local = (Path(tempdir) / f"{name}_{dt:%Y%m%d_%H}{suffix}")

        print(f"    > NCEI DATASET FOUND: {name}")

        # download the data
        _download_http(url, local)

        return {
            "source": f"NCEI {name}",
            "model": model,
            "pressure_file": str(local),
            "native_file": None,
        }

    return None







###############
# USER-FACING FUNCTION 
#########################################################################
def get_rapruc_data(latlon, year, month, day, hour, dataset, box_avg_size, hush):
    st = time.time()
    """
    Retrieve RAP/RUC analysis data and return a SounderPy clean_dict.

    Search order
    ------------
    * >= 2021-05-01: AWS RAP 13-km, then NCEI RAP archives
    * 2012-05-01 through 2021-04-30: NCEI RAP archives directly
    * before 2012-05-01: NCEI RUC archives directly

    ``dataset=`` bypasses automatic era selection and tests only that named
    NCEI archive family.

    All downloaded GRIB files are stored only inside a TemporaryDirectory
    and are automatically deleted before this function returns.

    Example
    -------
    clean_dict = get_rapruc_data("2024", "08", "28", "18", [44.58, -100.82])
    """

    dt = pd.Timestamp(year=int(year),month=int(month),day=int(day),hour=int(hour),)

    print("> RAP/RUC REANALYSIS DATA ACCESS FUNCTION")
    print("  ------------------------------------------")

    print(
        f"    > request: "
        f"{dt:%Y-%m-%d %HZ} | "
        f"{latlon[0]}, {latlon[1]}"
    )

    # Everything downloaded in this block is deleted automatically
    # when the block exits, including if parsing raises an exception.
    with tempfile.TemporaryDirectory(prefix="sounderpy_rapruc_") as tempdir:

        result = None

        if dataset is not None:
            print(f"    > requested NCEI dataset: {dataset}")
            result = try_ncei(dt, tempdir, dataset=dataset)

        elif dt >= AWS_PREFERRED_START:
            print("    > requested RAP date is > 05-2021 - trying AWS")
            result = try_aws(dt, tempdir, resolution=13)

            if result is None:
                print("    > AWS unavailable; falling back to NCEI RAP")
                result = try_ncei(dt, tempdir)

        elif dt >= RAP_START:
            print("    > requested RAP data is < 05-2021: trying NCEI")
            result = try_ncei(dt, tempdir)

        else:
            print("    > requested data is < 05-2012: going directly to NCEI RUC")
            result = try_ncei(dt,tempdir)

        if result is None:
            raise FileNotFoundError(
                "The data sources were reachable, but no matching RAP/RUC "
                f"file was found for {dt:%Y-%m-%d %HZ}.")

        print(f"    > DATA FOUND AT: {result['source']}")



        # PARSE GRIB FILES USING PARSER MODULE
        clean_data = parse_rapruc_grib(
            pressure_file=result["pressure_file"],
            native_file=result["native_file"],
            latlon=latlon,
            box_avg_size=box_avg_size,
            model=result["model"],
            valid_time=dt)

        clean_data["site_info"]["archive"] = result["source"]

        print("    > parsing complete")
        print("    > temporary GRIB files will now be deleted")


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



if __name__ == "__main__":

    data = get_rapruc_data("2024","08","28","18",[44.58, -100.82])

    print(data["site_info"])
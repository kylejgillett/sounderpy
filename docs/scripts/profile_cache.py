"""
Download each profile only once for the docs figure scripts, and keep it for later runs.

The figure scripts use the same soundings in many figures (e.g. OAX 2014-06-16 18Z). These stand-ins for
spy.get_obs_data(), spy.get_bufkit_data(), spy.get_model_data(), and an ACARS profile download save every profile
they fetch as a pickle in docs/scripts/cached_profiles/ (git-ignored), and read it from there the next time, in
this run or a later one, by any of the scripts. Each caller gets its own copy, so a figure that changes its
profile (smooth_profile, modify_sfc, a new title) doesn't change it for the next figure.

Requests count as the same when their arguments match, ignoring case, `hush`, and int vs. str ('08' and 8 differ,
as they do to SounderPy). A failed download isn't saved, so it is tried again the next time. Requests for the
latest data (a BUFKIT run with no date, or 'rap-now') are kept for one run only, since their answer changes.

Delete docs/scripts/cached_profiles/ to download everything again, e.g. after a change to how SounderPy reads or
cleans data (the pickles hold the cleaned profiles).

    from profile_cache import get_obs_data, get_bufkit_data, get_model_data, get_acars_profile
"""

from __future__ import annotations

import copy
import hashlib
import pickle
import re
from pathlib import Path

import sounderpy as spy

CACHE_DIR = Path(__file__).resolve().parent / "cached_profiles"

_PROFILES: dict[tuple, object] = {}
_STATS = {"downloaded": 0, "read from disk": 0, "reused": 0}


def _same(value):
    """A value as it's compared: strings and numbers ignoring case and type, anything else by its repr."""
    if isinstance(value, (str, int, float)) and not isinstance(value, bool):
        return str(value).strip().upper()
    return repr(value)


def _file_for(key: tuple) -> Path:
    """The pickle for a request: a readable name (the function and its arguments) and a short hash of the request."""
    name, args, kwargs = key
    words = [name, *args, *(f"{k}-{v}" for k, v in kwargs)]
    readable = re.sub(r"[^A-Za-z0-9.-]+", "_", "_".join(words)).strip("_")[:120]
    digest = hashlib.sha1(repr(key).encode()).hexdigest()[:10]
    return CACHE_DIR / f"{readable}_{digest}.pkl"


def _read(path: Path):
    try:
        with open(path, "rb") as file:
            return pickle.load(file)
    except Exception as exc:      # a partial or outdated pickle: download it again
        print(f"[profile cache] couldn't read {path.name} ({type(exc).__name__}); downloading it again")
        return None


def _write(path: Path, profile) -> None:
    CACHE_DIR.mkdir(exist_ok=True)
    partial = path.with_suffix(".tmp")
    with open(partial, "wb") as file:
        pickle.dump(profile, file, protocol=pickle.HIGHEST_PROTOCOL)
    partial.replace(path)         # whole files only, even if the script is stopped while writing


def _cached(name: str, fetch, *args, keep_on_disk: bool = True, **kwargs):
    key = (name, tuple(_same(arg) for arg in args),
           tuple(sorted((k, _same(v)) for k, v in kwargs.items() if k != "hush")))
    label = f"{name}{tuple(args)}"
    if key in _PROFILES:
        _STATS["reused"] += 1
        print(f"[profile cache] reusing {label}")
        return copy.deepcopy(_PROFILES[key])

    path = _file_for(key)
    profile = _read(path) if keep_on_disk and path.exists() else None
    if profile is not None:
        _STATS["read from disk"] += 1
        print(f"[profile cache] read {label} from {CACHE_DIR.name}/")
    else:
        profile = fetch(*args, **kwargs)
        _STATS["downloaded"] += 1
        if keep_on_disk:
            _write(path, profile)
    _PROFILES[key] = profile
    return copy.deepcopy(profile)


def get_obs_data(*args, **kwargs):
    return _cached("get_obs_data", spy.get_obs_data, *args, **kwargs)


def get_bufkit_data(*args, **kwargs):
    # with no run date, BUFKIT returns the latest run, which changes
    dated = len(args) >= 7 and args[3] is not None or kwargs.get("run_year") is not None
    return _cached("get_bufkit_data", spy.get_bufkit_data, *args, keep_on_disk=dated, **kwargs)


def get_model_data(*args, **kwargs):
    model = str(args[0] if args else kwargs.get("model", "")).casefold()
    return _cached("get_model_data", spy.get_model_data, *args, keep_on_disk=model != "rap-now", **kwargs)


def get_acars_profile(year, month, day, hour, profile_id, **kwargs):
    """One ACARS profile, like spy.acars_data(year, month, day, hour).get_profile(profile_id)."""
    def fetch(year, month, day, hour, profile_id, **kwargs):
        return spy.acars_data(year, month, day, hour).get_profile(profile_id, **kwargs)
    return _cached("acars", fetch, year, month, day, hour, profile_id, **kwargs)


def summary() -> str:
    """Where this run's profiles came from."""
    return (f"[profile cache] {_STATS['downloaded']} downloaded, {_STATS['read from disk']} read from "
            f"{CACHE_DIR.name}/, {_STATS['reused']} reused in this run")

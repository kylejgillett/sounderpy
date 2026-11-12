# FIXES NEEDED BEFORE RELEASE

- [ ] `_get_zoomed_background_lines`: on a thin or inverted zoomed view (e.g. fzra at `ylim=(1000, 900)`, sleet at (1000, 700)) the x-axis hits `ZOOM_MAX_SPAN`, the hot top-right corner gives MetPy `wet_bulb_potential_temperature` values in the thousands of °C, and up to ~1200 moist adiabats are drawn instead of ~12 (strict xfail in `tests/test_plot_components.py`)
- [ ] `utils.mag()` truncates instead of rounding (`int()`), so readouts drop toward zero (999.9 J/kg shows 999, -0.6 shows 0)
- [ ] Severe and winter docs/gallery images predate the PBL changes (one PBL bar, PBL row in the severe SRH/BWD box, winter PBL gust/mean/depth/LR rows and hodograph PBL arrow)
- [ ] Figures compute parameters on the `max_levels`-thinned profile, so level-snapped values change with thinning: GRAW 1 s sounding (UND 2026-01-25 20:26Z) PBL depth/LR is 520 m / 8.7 °C/km at full resolution, 589 m / 7.5 at the default 250 levels
- [ ] Skew-t omega panel: the "+10" / "-10" reference labels at 150 hPa overlap the "13km" height label (e.g. KGFK HRRR 2026-01-18 16Z)
- [ ] Skew-t: when the DGZ sits at the surface (surface T near -12 to -17 °C) its "DGZ" label overlaps the "-SFC (257m)-" label (KGFK HRRR 2026-01-18 16Z)
- [ ] `calc.py:824-845` ECAPE: on a zero-CAPE profile (no LFC, e.g. KABR 2024-08-28 12Z) `ecape_parcel` crashes on `lfc_z=None`, the bare `except` sets MU/ML/SB-ECAPE to masked (readout shows `--` instead of 0) and warns "could not be computed (calculation error)" three times; should be 0 with no warning when the parcel has no CAPE
- [x] Radiosonde readers (`from_graw`/`from_imet`/`from_windsond`): data gaps, launch time, GRAW midnight date, missing surface wind, print() notes, NaN temperatures, too-short profiles, Windsond elevation, Windsond BOM (uncommitted, 2026-11-11)
- [ ] `from_imet`: the profile starts at the sonde's PTU altitude (268 m in UND 2026-01-16 18:37Z) but the plot's surface label uses the station altitude (256 m)



# FUTURE FIXES
- [ ] `set_watermark`: settings are shared across a Python process, so a shared web app worker would affect all users
- [ ] MetPy's skew-T uses a Matplotlib function deprecated in 3.11 and due for removal in 3.13; confirm a MetPy fix before relaxing or raising version limits; MetPy 1.7.1's `StationPlot` text already fails to draw under Matplotlib 3.11 (its `TextCollection.draw` unpacks 4 layout values, 3.11 returns 3), so any station plot needs plain text + MetPy's symbol font until MetPy fixes it
- [ ] Station lists: `BUFKIT-STATIONS-MASTER.txt` has 8,527 rows for 2,093 IDs, and 23 IDs have conflicting entries (typos such as `B44` at lon -100.29 and -110.29, `B48` at -80.51 and -81.51; different stations sharing `CYQX`); the BUFKIT, METAR, buoy and airport lists are still downloaded from GitHub `main` on every call (RAOB/IGRA now ship in the package as `RAOB-STATIONS-MASTER.csv`)
- [ ] Full sounding hodograph: the storm-motion text block (top left) can overlap a ring speed label, e.g. "SSE @ 26 kt" over the "30kt" label on the 9ZT9MRK ship sounding, 2026-11-10 00Z



### Timeseries, Publication Figures, VAD hodos
- [ ] `build_publication_sounding`: reworked 2026-11-11 (panels, layouts, widths, map, parameter table, `(fig, axes)`); tests in `sandbox/plot_api/test_publication.py` until the module is tracked; no docs yet
- [ ] `build_editorial_sounding`: experimental; choose whether to keep it, and it has no tests or docs yet
- [ ] `build_timeseries`: parameters with the same unit share an axis (e.g. 0-3 km MUCAPE flattens next to MU-ECAPE); only the built-in SCP and CIN axes can be split, users can't choose
- [ ] Timeseries docs: no docs page yet for `timeseries_params()` / `build_timeseries()`
- [ ] `build_timeseries` kinematics panel: six lines on two axes get busy on long runs (e.g. 84-hour RRFS). Options to revisit: split into `srh` and `shear` panels; note solid = SRH, dashed = BWD on the axis labels; align right-axis ticks to the left gridlines; group the legend by SRH/BWD; drop a line
- [ ] VAD `calc`: 0-1 km SRH is really 1-3 km and 1-3 km is really 1-4 km (`calc.py:1086-1091`); Corfidi jet search never increments `i`, so it uses the surface wind (waiting until there is a VAD data fetcher)
- [ ] `plot/timeseries.py`: surface wind barbs are clipped to the shaft; right-axis parameters without a style reuse the left axis's colors
- [ ] `build_publication_sounding`, `build_editorial_sounding`: Type 42 font settings are lost when the returned figure is saved later
- [ ] Tests: nothing covers `timeseries_params`, `build_timeseries`, `build_sounding(sr_hodo=True)`, or `all_fcst_hours`
- [ ] Timeseries (`timeseries.py`, `plot/timeseries.py`) and publication/editorial figures (`plot/publication.py`) are on hold: unhooked from `import sounderpy` and git-ignored (use locally with `from sounderpy.timeseries import timeseries_params, build_timeseries`). `timeseries.py` and `publication.py` are still tracked until `git rm --cached`; setuptools packages any `.py` on disk, so build the pip release from a clean clone. Their open items above wait until they come back
- [ ] Old module paths (`sounderpy.bufkit_data`, `sounderpy.file_io`, etc.) are stand-ins for `sounderpy.data.*` that show a `FutureWarning`; remove the block at the end of `sounderpy/__init__.py` in a later release
- [ ] Old reader names `make_cm1_profile()` / `make_wrf_profile()` are stand-ins for `from_cm1()` / `from_wrf()` that show a `FutureWarning`; remove the wrappers at the end of `data/cm1_utils.py` and `data/wrf_utils.py` (and their exports) in a later release

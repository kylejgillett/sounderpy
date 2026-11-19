# FIXES NEEDED BEFORE RELEASE

- [x] Deviant tornado motion (`calc.py`, sounding and VAD): uses the 0-300 m mean wind, but Nixon & Allen (2021) define DTM with the 0-500 m ("lowest half-kilometer") mean wind; decide whether to change it (2026-11-19: user keeps 0-300 m)
- [x] `_get_zoomed_background_lines`: on a thin or inverted zoomed view (e.g. fzra at `ylim=(1000, 900)`, sleet at (1000, 700)) the x-axis hits `ZOOM_MAX_SPAN`, the hot top-right corner gives MetPy `wet_bulb_potential_temperature` values in the thousands of °C, and up to ~1200 moist adiabats are drawn instead of ~12 (strict xfail in `tests/test_plot_components.py`) (fixed 2026-11-18: moist adiabats follow `moist_lapse` to 1000 hPa)
- [x] `utils.mag()` truncates instead of rounding (`int()`), so readouts drop toward zero (999.9 J/kg shows 999, -0.6 shows 0) (fixed 2026-11-18: rounds)
- [ ] Severe and winter docs/gallery images predate the PBL changes (one PBL bar, PBL row in the severe SRH/BWD box, winter PBL gust/mean/depth/LR rows and hodograph PBL arrow)
- [x] Skew-t: when the DGZ sits at the surface (surface T near -12 to -17 °C) its "DGZ" label overlaps the "-SFC (257m)-" label (KGFK HRRR 2026-01-18 16Z) (fixed 2026-11-18 in the winter figure: label lifted above the SFC label)
- [x] Radiosonde readers (`from_graw`/`from_imet`/`from_windsond`): data gaps, launch time, GRAW midnight date, missing surface wind, print() notes, NaN temperatures, too-short profiles, Windsond elevation, Windsond BOM (uncommitted, 2026-11-11)
- [ ] `from_imet`: the profile starts at the sonde's PTU altitude (268 m in UND 2026-01-16 18:37Z) but the plot's surface label uses the station altitude (256 m)



# PUBLISHING PROCESS TO-DOS
- [ ] What's New (`docs/source/index.rst`): add GRAW/iMET/Windsond readers, `make_cm1_profile`/`make_wrf_profile` -> `from_cm1`/`from_wrf`, RRFS BUFKIT, packaged RAOB station list (ICAO fallback, leading-zero WMO IDs), IGRA RH -> Td, data source statements on figures, new composite params, SHARPpy HGZ fix, darker EIL with top/bottom heights, `smooth_profile()` Savitzky-Golay default, HRRR & RAP NOMADS BUFR fallback, `map_layers` / `map_label`, 0-10 km BWD & Δθe (Cohen et al. 2007), EIL depth, Parameters & References docs page
- [ ] README.md: replace the "SounderPy 3.2" section with 3.3
- [x] Clean up code comments and the module header blocks at the top of each file (purpose, (C) years) (done 2026-11-19: every module opens with a SOUNDERPY header docstring: purpose, docs links, verified references, (C) line; Sczepanski et al. cited as the 2024 AMS abstract until the journal paper is out)
- [x] Audit every comment in the package for ones written to the developer instead of the user (e.g. "didn't change the param, left as is", "kept the old behavior", notes about edits made): remove them, or rewrite real changes as user-facing notes like "changed in v3.3.0: updated axis ticks" (done 2026-11-19: [prototype] tags -> [SounderPy], ECAPE/kwarg/sys.exit/deprecation notes rewritten, RAP-RUC module docstring rewritten; user kept the first-person comments in model_data.py & bufkit_data.py)
- [ ] Build the docs and check them (winter docs never built)
- [ ] Build the 3.3.0 wheel/sdist from a clean clone (git-ignored `timeseries.py` / `publication.py` are on disk)



# NEW FEATURES
- [x] Map inset: let users draw their own 2D x/y data on it (e.g. `.contourf()`) (done 2026-11-18: `map_layers=` functions + `map_label=` on build_sounding (severe & winter), build_hodograph, build_vad_hodograph; Python only, not in the CLI)
- [x] Map inset tutorial for `map_layers` / `map_label` (done 2026-11-18: `tutorials/custom_map_data.rst`, IEM METARs + AWS RAP on the ABR 2024-08-29 00Z RAOB; images from `make_gallery_figures.py --sections custom-map`)
- [x] BUFKIT: add the HRRR to the NOMADS fallback (only RRFS falls back now; PSU HRRR runs behind), if NOMADS serves HRRR BUFR soundings (done 2026-11-18: RRFS, HRRR & RAP share `ncep_bufr.NOMADS_MODELS`)
- [x] Skew-t: EIL top and bottom heights next to the EIL annotation (done 2026-11-18: heights beside the end ticks, EIL moved to x=0.22, color #4a9fd0)
- [x] User request: 0-10 km shear and the min-max theta-e difference (Cohen et al. 2007 derecho discriminators, https://journals.ametsoc.org/view/journals/wefo/22/5/waf1040_1.xml) (done 2026-11-19: `kinem['shear_0_to_10000']`, `thermo['thetae_diff']` = min theta-e from the 0-3 km max level to 7 km minus that max; new "0-10km BWD, Δθe" row in the SRH/BWD box beside the hodograph)
- [x] Severe sounding: show the EIL depth in the thermo parameter table (done 2026-11-19: `kinem['eil_depth']`, "EIL:" in the thermo table's top row)
- [x] `smooth_profile()` (`profile_utils.py`): add a kwarg to choose a Savitzky-Golay filter (new default) or the current Gaussian (`sigma`); update `helpertools.rst` (done 2026-11-18: `method=`, `window_length=5`, `polyorder=2`; `sigma` alone still means Gaussian)
- [x] Docs: Parameters & References page (`reference/references.rst`), generated by `docs/scripts/make_references_page.py` from `docs/scripts/sounderpy_references.py`, the same list `docs/scripts/module_headers.py` writes into the module headers (both scripts have `--check`)


# POSSIBLE FUTURE ADDITIONS
- [ ] Built-in METARs on the map: `map_obs=None | 'dewpoint' | 'temperature'` (+ `--map-obs` CLI flag) on all map figures; fetch stations inside the map extent in one IEM request (reuse `data/sfc_obs.py` station table), draw/thin like the custom map tutorial; one value per station (MetPy `StationPlot` fails under Matplotlib 3.11); source: IEM archive at the sounding's valid time (works for any date), with an optional fast path from the aviationweather.gov bulk METAR cache when the valid time is recent (whole-world file, latest reports only; user has a working AWC reader in ~/und-atsc-weather-wall); forecast profiles have no obs at their valid time, so offer a 'now' option like radar_time
- [ ] User request: plot EPV, e.g. to show instability above a frontal zone/inversion; true EPV needs horizontal theta gradients a single sounding doesn't have, so a sounding-only version would estimate them from the thermal wind (shear vector)



# FUTURE FIXES
- [ ] Calc module refactor: move `calc.py` and `calc_winter.py` into a `calc/` dir (like `data/`) with a `calc_utils.py` of shared helpers (e.g. `layer_mean`, `layer_max`, `value_at` from `calc_winter.py`) and `FutureWarning` stand-ins for the old paths; revisit once new calc work (Cohen 2007 params, EPV) shares helpers
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

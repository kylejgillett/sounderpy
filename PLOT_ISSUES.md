# FIXES NEEDED BEFORE RELEASE

- [ ] `_get_zoomed_background_lines`: on a thin or inverted zoomed view (e.g. fzra at `ylim=(1000, 900)`, sleet at (1000, 700)) the x-axis hits `ZOOM_MAX_SPAN`, the hot top-right corner gives MetPy `wet_bulb_potential_temperature` values in the thousands of °C, and up to ~1200 moist adiabats are drawn instead of ~12 (strict xfail in `tests/test_plot_components.py`)
- [ ] `utils.mag()` truncates instead of rounding (`int()`), so readouts drop toward zero (999.9 J/kg shows 999, -0.6 shows 0)
- [ ] `plot/components.py`: no public API for building custom figures; every drawing function is private; will need to determine the best option for this, what users should have and shouldn't have
- [ ] Check the CLI. Test for bugs, errors, and missing features.


# FUTURE FIXES
- [ ] `set_watermark`: settings are shared across a Python process, so a shared web app worker would affect all users
- [ ] MetPy's skew-T uses a Matplotlib function deprecated in 3.11 and due for removal in 3.13; confirm a MetPy fix before relaxing or raising version limits; MetPy 1.7.1's `StationPlot` text already fails to draw under Matplotlib 3.11 (its `TextCollection.draw` unpacks 4 layout values, 3.11 returns 3), so any station plot needs plain text + MetPy's symbol font until MetPy fixes it
- [ ] Station lists: `BUFKIT-STATIONS-MASTER.txt` has 8,527 rows for 2,093 IDs, and 23 IDs have conflicting entries (typos such as `B44` at lon -100.29 and -110.29, `B48` at -80.51 and -81.51; different stations sharing `CYQX`); all station lists are downloaded from GitHub `main` on every call rather than read from the package



### Timeseries, Publication Figures, VAD hodos
- [ ] `build_publication_sounding`: experimental; hodograph height labels can overlap with several profiles or tight traces, and no tests or docs yet
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

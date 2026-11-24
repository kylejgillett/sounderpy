# SounderPy tests

Everything runs offline on synthetic profiles unless you ask for live data.

```bash
.venv/bin/python -m pytest -q                  # everything offline (~1.5 min)
.venv/bin/python -m pytest -q -m "not slow"    # skip full-figure renders (~20 s)
.venv/bin/python -m pytest -q --network        # also contact live data services
```

| File | Covers |
| --- | --- |
| `helpers.py` | synthetic profiles (`make_clean_profile`, `make_winter_profile`, `WINTER_CASES`, ...) and shared assertions |
| `conftest.py` | headless Matplotlib, the `--network` switch, profile fixtures, and session-cached `severe_params` / `winter_params_for` |
| `test_calc.py` | `modify_surface` (mix, replace, blend) and `sounding_params(parcel_depth=...)` |
| `test_calc_winter.py` | winter calc: precip type, Bourgouin energies, mixed layer, DGZ, snow ratio |
| `test_plot_components.py` | `plot/components.py` helpers: thinning, `ylim` zoom, labels, watermark, missing values, panels and tables |
| `test_plot_components_winter.py` | `plot/components_winter.py` helpers, drawn one at a time from the winter cases |
| `test_figures.py` | end-to-end `build_sounding` (severe and winter), `build_hodograph`, `build_composite`, `build_vad_hodograph`, `return_fig`, and `map_layers` / `map_label` (drawn offline, without the satellite image) |
| `test_publication.py` | `build_publication_sounding`: panels, layouts, widths, options and their errors, the parameter table, map and radar, and the saved size |
| `test_retrieval.py` | routing of the `get_*` functions, the RAOB/IGRA station list and lookup, surface observations, and live retrieval (`network`) |
| `test_file_io.py` | SHARPpy and CSV files |
| `test_readout.py` | web readout files (checked against the saved PNG's pixels), `apply_profile_edits`, and the EDITED PROFILE note |
| `test_readout_js.py` + `js/test_readout.mjs`, `js/test_readout_winter.mjs` | the browser script, `sounderpy-readout.js`, run with Node's test runner on a readout file SounderPy just made: pixel mappings, interpolation, shear, SRH, storm-relative wind, storm motions, SB/ML/MU parcels, layer stats, snapping and thickness, arrow keys, edits, and units; and, for each idealized winter case, precipitation type, melting/refreezing layers, and cloud layers against `calc_winter` (`js/readout_view.mjs` is shared setup; skipped without Node.js 18+) |
| `test_cli.py` | the command-line interface |

Markers: `slow` renders a full figure; `network` contacts a live service and is skipped without `--network`.
Only `test_build_sounding` draws the default MUECAPE parcel (about 15 s in the ecape_parcel library); other severe
renders pass `special_parcels="simple"`.

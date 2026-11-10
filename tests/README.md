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
| `test_figures.py` | end-to-end `build_sounding` (severe and winter), `build_hodograph`, `build_composite`, `build_vad_hodograph` |
| `test_retrieval.py` | routing of the `get_*` functions, surface observations, and live retrieval (`network`) |
| `test_file_io.py` | SHARPpy and CSV files |
| `test_cli.py` | the command-line interface |

Markers: `slow` renders a full figure; `network` contacts a live service and is skipped without `--network`.
Only `test_build_sounding` draws the default MUECAPE parcel (about 15 s in the ecape_parcel library); other severe
renders pass `special_parcels="simple"`.

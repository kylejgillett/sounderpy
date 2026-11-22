<div align="center">

# SounderPy

### Atmospheric vertical-profile retrieval, analysis, and visualization for Python

Retrieve observations, forecasts, and reanalyses.  
Analyze thermodynamic and kinematic environments.  
Create publication-ready soundings and hodographs.

[![PyPI](https://img.shields.io/pypi/v/sounderpy.svg)](https://pypi.org/project/sounderpy/)
[![Python](https://img.shields.io/pypi/pyversions/sounderpy.svg)](https://pypi.org/project/sounderpy/)
[![Conda Version](https://img.shields.io/conda/vn/conda-forge/sounderpy.svg)](https://anaconda.org/conda-forge/sounderpy)
[![License](https://img.shields.io/pypi/l/sounderpy.svg)](https://github.com/kylejgillett/sounderpy/blob/main/LICENSES)
[![JOSS](https://joss.theoj.org/papers/10.21105/joss.08087/status.svg)](https://doi.org/10.21105/joss.08087)

**[Documentation](https://sounderpy.app/docs/)** ·
**[Tutorials](https://sounderpy.app/docs/tutorials/index.html)** ·
**[Plot Gallery](https://sounderpy.app/docs/gallery/examplegallery.html)** ·
**[Operational Site](https://sounderpy.app/)** ·
**[PyPI](https://pypi.org/project/sounderpy/)**

</div>

---

![SounderPy example sounding](https://sounderpy.app/docs/_static/gallery/themes/sounding_light.png)

## What is SounderPy?

**ABSTRACT:** SounderPy is a simple, open-source Python package for retrieving and plotting 
vertical profile (sounding) data. Built for simplicity and reliability for all uses and users, this 
project’s goal is to provide a uniform method for sounding analysis across multiple data types. 
Severe weather analysis and forecasting requires a sound comprehension of thermodynamic and 
kinematic properties of the environment. SounderPy makes this possible with robust access to 
data and custom visualizations. The tool creates complex yet effective sounding and hodograph 
plots with high readability which are designed specifically for severe weather analysis and 
forecasting. SounderPy is capable of retrieving and plotting model forecast data, observed 
radiosonde data, Aircraft Communications Addressing and Reporting System (ACARS) 
observation data, and model reanalysis data. All of this functionality can be completed in three
simple lines of code or less, making SounderPy an accessible tool for both Python experts and 
novices. A number of scientific Python libraries build the base of SounderPy’s efficient and 
durable functionality, such as NumPy, Matplotlib, xarray, Metpy, and SHARPpy. SounderPy is 
available through GitHub and PyPi and is distributed under an MIT license. 

SounderPy has been used by several institutions. For example, this tool has been implemented by the Des Moines, Columbia, and Grand Rapids National Weather Service Offices, the State University of New York at Albany, Mississippi State University, the University of North Dakota, and others. Many students at various universities have used SounderPy in projects, posters, and papers, such as students at The University of Oklahoma, Ohio State University, Central Michigan University, Iowa State University, & Rizal Technological University.

## Why SounderPy?

- **Retrieve data without writing retrieval/parsing code yourself.**
- **Use one workflow across observations, forecasts, reanalyses, and model data.**
- **Calculate severe-weather thermodynamic and kinematic parameters.**
- **Create sounding, hodograph, VAD, and composite-sounding figures.**
- **Customize parcels, storm motion, radar, maps, themes, and plot accessibility.**
- **Export profiles to CSV, SHARPpy, and CM1 formats.**
- **Use SounderPy from Python or directly from the command line.**

---

## Installation

### pip

```bash
pip install sounderpy
```

### conda-forge

```bash
conda install conda-forge::sounderpy
```

or:

```bash
mamba install conda-forge::sounderpy
```

Then import SounderPy:

```python
import sounderpy as spy
```

<br>
<br>

## Quickstart

Retrieve the 18 UTC OAX observed sounding from 16 June 2014:

```python
import sounderpy as spy

data = spy.get_obs_data(
    "OAX",
    "2014", "06", "16", "18",
    hush=True,
)

spy.build_sounding(data)
```

That's it.

`data` is now a SounderPy `clean_data` dictionary containing the vertical profile
and associated metadata.

Typical fields include:

```text
p          pressure
z          height
T          temperature
Td         dewpoint
u          u-component wind
v          v-component wind
site_info  profile metadata
```

Explore more examples in the
**[SounderPy Plot Gallery](https://sounderpy.app/docs/gallery/examplegallery.html)**.

<br>
<br>

## Command Line Interface

SounderPy can also be used without writing a Python script.

### Retrieve an observation

```bash
sounderpy obs OAX 2014-06-16 18
```

### Retrieve RAP/RUC

```bash
sounderpy model rap-ruc 44.58 -100.82 2024-08-28 18 \
    --box-size 0.25
```

### Create a plot

```bash
sounderpy obs OAX 2014-06-16 18 \
    --plot sounding \
    --map-zoom 0 \
    --plot-file oax_sounding.png
```

### Export data

```bash
sounderpy obs OAX 2014-06-16 18 \
    --output oax.csv
```

---

<br>
<br>

## Case Studies and Tutorials

The documentation includes task-oriented tutorials that build complete
SounderPy workflows rather than simply listing function arguments.

Topics include:

- plot appearance and accessibility;
- maps and radar;
- parcel settings;
- hodograph configuration;
- observations, forecasts, and reanalyses;
- composite soundings;
- data export;
- command-line workflows;
- real severe-weather case studies.

Start here: **[SounderPy Tutorials →](https://sounderpy.app/docs/tutorials/index.html)**

---

<br>
<br>

## What's New in SounderPy 3.3

Version 3.3 rebuilds the plotting module, adds winter and publication figures, and makes every
figure easier to customize.

### 📈 Figure Updates

- **Rebuilt plot module**: most figures are about 50% faster, with custom colors (`theme`), a custom
  watermark (`spy.set_watermark()`), and a `ylim` to zoom the skew-T.
- **New `build_sounding()` kwargs**: `max_levels` (reduce the number of vertical levels of a profile), `ylim` (zoom the skew-T into a particular layer), `theme` (custom figure colors), `sr_hodo` (storm-relative hodograph), `mode` (`'severe'`, `'winter'`, or `'auto'`), `show_obs` (the nearest surface observation for `mode='winter'` figures), `map_layers` and `map_label` (adding your own data to the map inset), and `return_fig` (return the figure and its panels for customization). Many of these also work on `build_hodograph()`, `build_composite()`, and `build_vad_hodograph()`.
- **Publication figures (experimental)**: `spy.build_publication_sounding()` draws journal-sized figures
  for AMS and AGU journals, and every `build_*` function takes `return_fig=True` to hand back its figure
  ([custom figures tutorial](https://sounderpy.app/docs/tutorials/custom_figures.html)).
- **Your own data on the map inset** with `map_layers` and `map_label`
  ([custom map data tutorial](https://sounderpy.app/docs/tutorials/custom_map_data.html)).
- **Rebuilt composite soundings** with figure improvements and parameter read-outs, and rebuilt VAD hodographs with map, radar, and boundary options.
- **New parameter readouts** on the full sounding: composite indices in the thermodynamics box, EIL depth and
  heights, 0-10 km bulk shear and theta-e difference, PBL height and transport wind, improved layout for ML/MU parcel heights, updated watermark, new data source attribution line, and numerous other small improvements to the plot layout.

### ❄️ NEW Winter Soundings

- `build_sounding(mode='winter')` draws a winter-focused sounding built around precipitation type, snow growth, and
  winter hazards. Setting `mode='auto'` will generate a winter figure when the surface wet-bulb is below 35°F, and a new `show_obs=True` kwarg adds the latest (nearest) surface observation to the figure ([See the winter soundings tutorial!](https://sounderpy.app/docs/tutorials/winter_soundings.html)).
- `spy.winter_params(clean_data).calc()` returns the dictionary of winter parameters.

### 🧮 Calculations

- NEW severe parameters: 0-10 km bulk shear and the low-to-midlevel theta-e difference (for MCSs; Cohen et al. 2007), PBL depth and transport wind (for air quality), EIL bottom height, top height, and depth, SHIP, MLBRN, TT, KI, and 700-500 hPa lapse rate.
- UPDATED/CORRECTED parameters: Streamwise vorticity bug counted antistreamwise layers as streamwise; STP/SCP used EIL shear in kts (now m/s).
- Every parameter, its units, and its source are listed on the
  [Parameters & References](https://sounderpy.app/docs/reference/references.html) documentation page.
- Settings: `modify_sfc` adds `'mix'`, `'blend'`, and `'replace'` methods, and the NEW `spy.modify_surface()` returns the modified profile itself, to compare or save; `spy.smooth_profile()` adds a Savitzky-Golay filter.
 
### 🌦️ Data Sources

- BUFKIT now serves the RRFS; NAM, NAMNEST, SREF, and HIRESW are no longer available as "latest runs" (discontinued; but NAM and NAMNEST are available as archived runs)
- RRFS, HRRR, and RAP may now fall back to NOMADS NCEP's BUFR dataset when Penn State is behind or an archived run is missing.
- NEW Station list: observed soundings use a packaged WMO/ICAO/IGRAv2 station list, with an ICAO fallback for the
  Wyoming archive. This list is a best-guess combination of 5 data sources to limit incorrect site names, IDs, and locations.
- Dewpoints are now calculated (when possible) for older IGRA v2 soundings that don't ship dewpoint.
- Thread-safe RAP/RUC GRIB parsing.


### 💻 Command Line Interface Updates

- NEW `file` command: read a local SHARPpy-format sounding (`--latlon LAT LON` adds its location for the map).
- NEW options:
  - `--params`: print key calculated parameters (CAPE, SRH, shear, and precip type with `--mode winter`);
  - `--mode {severe,winter,auto}` and `--obs`: winter soundings and the nearest surface observation;
  - `--storm-motion`: `right_moving`, `left_moving`, `mean_wind`, or a direction and speed such as `270 25`;
  - `--modify-sfc`: change the surface, e.g. `T=30 Td=21` or `T=-6 method=blend depth=300`;
  - `--radar {mosaic,single,none}` and `--radar-time {sounding,now}`: the map's radar;
  - `--parcels` and `--background-parcels`: custom parcel codes.
- Options now match the Python functions, with clearer warnings.

### 🐞 Bug Fixes

- ACARS airport elevations corrected, CIN doesn't accumulate above the LFC, HGZ identification, extra moist adiabats on thin zoomed skew-Ts, and figures now close after `save=True` to preserve memory.

### 📚 NEW HOME for Documentation!

- Freshly redesigned docs now live at [sounderpy.app/docs](https://sounderpy.app/docs/), with new Plot Settings, Publication Figures, and Parameters & References pages.
- Also, each source module now includes a text header with its docs and references.

### ⚠️ Changes to Note

- `sounding_params(...).calc()` and `vad_params(...).calc()` return a dictionary instead of a tuple.
- Renamed, with the old names still working (with a `FutureWarning`) until a later release:
  `build_composite()`'s `colors_to_use`, `ls_to_use`, `alphas_to_use`, and `lw_to_use` (now `colors`,
  `linestyles`, `alphas`, and `linewidths`), `make_cm1_profile()` / `make_wrf_profile()` (now
  `spy.from_cm1()` / `spy.from_wrf()`), and the data modules (now in `sounderpy.data`).
- `spy.smooth_profile()` defaults to the Savitzky-Golay filter (passing `sigma` alone still uses the
  Gaussian filter).

See the full [What's New](https://sounderpy.app/docs/) on the docs landing page and the repository's
**[release history](https://github.com/kylejgillett/sounderpy/releases)**.

---

<br>
<br>

## Documentation and Project Links

- **Documentation:** https://sounderpy.app/docs/
- **Tutorials:** https://sounderpy.app/docs/tutorials/index.html
- **Plot Gallery:** https://sounderpy.app/docs/gallery/examplegallery.html
- **GitHub:** https://github.com/kylejgillett/sounderpy
- **PyPI:** https://pypi.org/project/sounderpy/
- **conda-forge:** https://anaconda.org/conda-forge/sounderpy
- **Operational SounderPy site:** https://sounderpy.app/

Questions, bugs, and feature ideas are welcome through
**[GitHub Issues](https://github.com/kylejgillett/sounderpy/issues)**.

---

<br>
<br>

## Citing SounderPy

If SounderPy contributes to research or a publication, please cite the
SounderPy JOSS article:

> Gillett, K. J., 2025: SounderPy: An atmospheric sounding visualization and
> analysis tool for Python. *J. Open Source Software*, **10** (112), 8087.
> https://doi.org/10.21105/joss.08087

BibTeX:

```bibtex
@article{Gillett2025SounderPy,
  author  = {Gillett, Kyle J.},
  title   = {SounderPy: An atmospheric sounding visualization and analysis tool for Python},
  journal = {Journal of Open Source Software},
  year    = {2025},
  volume  = {10},
  number  = {112},
  pages   = {8087},
  doi     = {10.21105/joss.08087}
}
```

---

<br>
<br>

## Authors and Contributors

**Author**

- **Kyle J. Gillett** — University of North Dakota

**Contributors**

- **Scott Thomas**, NWS Grand Rapids — VWP hodograph and buoy-site listing
- **Amelia R. H. Urquhart**, University of Oklahoma — `ecape-parcels`
- **Daryl Herzmann**, Iowa State University — SounderPy conda-forge feedstock
- **Ryan Vandersmith** — stepwise CAPE/CIN plot
- **David E. Levin**, NOAA/NWS — Alaskan RAOB site identifiers

Additional contributors are recognized through the repository history.

---

<br>
<br>


## License

SounderPy is distributed under the
**[MIT License](https://github.com/kylejgillett/sounderpy/blob/main/LICENSES)**.

---

<div align="center">

[Documentation](https://sounderpy.app/docs/) ·
[GitHub](https://github.com/kylejgillett/sounderpy) ·
[PyPI](https://pypi.org/project/sounderpy/) ·
[JOSS](https://doi.org/10.21105/joss.08087)

</div>

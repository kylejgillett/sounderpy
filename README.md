<div align="center">

# SounderPy

### Atmospheric vertical-profile retrieval, analysis, and visualization for Python

Retrieve observations, forecasts, and reanalyses.  
Analyze thermodynamic and kinematic environments.  
Create publication-ready soundings and hodographs.

[![PyPI](https://img.shields.io/pypi/v/sounderpy.svg)](https://pypi.org/project/sounderpy/)
[![Python](https://img.shields.io/pypi/pyversions/sounderpy.svg)](https://pypi.org/project/sounderpy/)
[![Conda Version](https://img.shields.io/conda/vn/conda-forge/sounderpy.svg)](https://anaconda.org/conda-forge/sounderpy)
[![License](https://img.shields.io/pypi/l/sounderpy.svg)](https://github.com/kylejgillett/sounderpy/blob/main/LICENSE)
[![JOSS](https://joss.theoj.org/papers/10.21105/joss.08087/status.svg)](https://doi.org/10.21105/joss.08087)

**[Documentation](https://kylejgillett.github.io/sounderpy/)** ·
**[Tutorials](https://kylejgillett.github.io/sounderpy/tutorials/index.html)** ·
**[Plot Gallery](https://kylejgillett.github.io/sounderpy/examplegallery.html)** ·
**[Operational Site](https://sounderpysoundings.anvil.app/)** ·
**[PyPI](https://pypi.org/project/sounderpy/)**

</div>

---

![SounderPy example sounding](https://kylejgillett.github.io/sounderpy/_static/gallery/themes/sounding_light.png)

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
hodograph visualizations.

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
**[SounderPy Plot Gallery](https://kylejgillett.github.io/sounderpy/examplegallery.html)**.

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

Start here: **[SounderPy Tutorials →](https://kylejgillett.github.io/sounderpy/tutorials/index.html)**

---

<br>
<br>

## SounderPy 3.2

Version 3.2 expands both SounderPy's data-access workflow and the ways it can be
used.

Major changes include:

- **new command line interface** for retrieval, plotting, JSON, and file export;
- **modernized RAP/RUC retrieval**, including GRIB2-based archive access;
- expanded automated testing and packaging infrastructure;
- substantially expanded documentation, tutorials, and visual examples.

See the repository's
**[release history](https://github.com/kylejgillett/sounderpy/releases)**
for version-specific details.

---

<br>
<br>

## Documentation and Project Links

- **Documentation:** https://kylejgillett.github.io/sounderpy/
- **Tutorials:** https://kylejgillett.github.io/sounderpy/tutorials/index.html
- **Plot Gallery:** https://kylejgillett.github.io/sounderpy/examplegallery.html
- **GitHub:** https://github.com/kylejgillett/sounderpy
- **PyPI:** https://pypi.org/project/sounderpy/
- **conda-forge:** https://anaconda.org/conda-forge/sounderpy
- **Operational SounderPy site:** https://sounderpysoundings.anvil.app/

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
**[MIT License](https://github.com/kylejgillett/sounderpy/blob/main/LICENSE)**.

---

<div align="center">

[Documentation](https://kylejgillett.github.io/sounderpy/) ·
[GitHub](https://github.com/kylejgillett/sounderpy) ·
[PyPI](https://pypi.org/project/sounderpy/) ·
[JOSS](https://doi.org/10.21105/joss.08087)

</div>

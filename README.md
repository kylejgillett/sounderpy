<div align="center">
<img src="https://github.com/kylejgillett/sounderpy/assets/100786530/2e9477c9-e36a-4163-accb-fe46780058dd" width="250">

# SOUNDERPY | Vertical Profile Data Retrieval and Analysis Tool For Python
LATEST VERSION: v2.0.6 |  RELEASED: October 4th, 2023  |  COPYRIGHT Kyle J Gillett, 2023
### [VISIT SOUNDERPY DOCUMENTATION HERE](https://github.com/kylejgillett/sounderpy/wiki)
#### [CHECK OUT AN EXAMPLE NOTEBOOK](https://github.com/kylejgillett/sounderpy/blob/main/examples/sounderpy_tutorial.ipynb)
A Python package that helps you to access and plot vertical profile data for meteorological analysis 

[![PyPI Package](https://img.shields.io/pypi/v/sounderpy.svg)](https://pypi.python.org/pypi/sounderpy/)
[![PyPI Downloads](https://img.shields.io/pypi/dm/sounderpy.svg)](https://pypi.python.org/pypi/sounderpy/)
[![PyPI license](https://img.shields.io/pypi/l/ansicolortags.svg)](https://github.com/kylejgillett/sounderpy/blob/main/LICENSE.txt)
[![PyPI pyversions](https://img.shields.io/pypi/pyversions/sounderpy.svg)](https://pypi.python.org/pypi/sounderpy/)
[![GitHub commits](https://badgen.net/github/commits/kylejgillett/sounderpy)](https://GitHub.com/kylejgillett/sounderpy/commit/)
[![Maintainer](https://img.shields.io/badge/maintainer-kylejgillett-blue)](https://github.com/kylejgillett)
[![made-with-python](https://img.shields.io/badge/Made%20with-Python-1f425f.svg)](https://www.python.org/)

</div>



-----
## What is SounderPy:

- SounderPy is a Python package used to access vertical profile data for calculations or plotting of a vertical profile (sounding). SounderPy's main use is for getting the data, but some basic plotting tools are included.

## SounderPy is currently capable of accessing and processing data from:

- ECMWF CDS ERA5 reanalysis [1940-present] *note: you must set up an account through the CDS to unlock ERA5 data. (see: https://cds.climate.copernicus.eu/api-how-to)
- UNIDATA THREDDS TDS RAP reanalysis [2005-present]
- UNIDATA THREDDS TDS RUC reanalysis [2005-2020]
- The University of Wyoming RAOB archive [1973-present, depending on station]
- Iowa State University's RAOB archive [1945-present, depending on station]
- The IGRAv2 Observed profile archive [1905-present, depending on station]
- Iowa State University's BUFKIT archive [2011-present, depending on station & model]
- Penn State University's BUFKIT archive [most recent run, depending on station & model]
- UNIDATA THREDDS TDS RAP-most-recent-analysis [now, most recent analysis only]
- OU Aircraft Communications, Addressing and Reporting System (ACARS) [2019-present]
- NCEP FNL 0.25deg Gridded Reanalysis Dataset [2021-present]


## Why SounderPy?
- Sometimes data is tough to find, and often times is even tougher to get it in the format you like. SounderPy gets you this data!
- The code needed for loading and parsing vertical data (especially from models) can be large and messy. SounderPy keeps it hidden away in a PyPi package -- just import and call sounderPy functions to keep your code clean!

-------

## How to use SounderPy:
1. Make sure your environment has the required dependencies:
   - cdsapi>=0.6.1
   - ecape>0.0.0
   - matplotlib>=3.3.0, <=3.7.1
   - metpy>=1.5.1
   - netcdf4>=1.6.4
   - numpy>=1.20.0
   - pandas>=1.2.0
   - siphon>=0.9
   - scipy>= 1.10.1
   - xarray>=0.18.0

2. ```
   pip install sounderpy
   ```
   Find it at https://pypi.org/project/sounderpy/
3. ```
   import sounderpy as spy
   ```
4. ```
   year  = '2011' 
   month = '04'
   day   = '27'
   hour  = '22'
   latlon = [33.19, -87.46]
   method = 'rap' 
   ```
5. ```
   raw_data = spy.get_model_data(method, latlon, year, month, day, hour)
   ```
6. ```
   clean_data = spy.parse_data(raw_data)
   ```
------
  and boom! Now you have a callable dictionary of vertical profile reanalysis data including... 
  
1. Temperature
2. Dewpoint
3. Relative Humidity
4. Pressure
5. Height 
6. U-component Wind 
7. V-component Wind


You can make a quick sounding plot of the data using built-in MetPy plotting functions! Just call...
`spy.metpy_sounding(clean_data)`
<div align="center">
<img src="https://raw.githubusercontent.com/kylejgillett/sounderpy/main/images/sounderpy_v2-0-6_example-sounding.png" width="600">
</div>


or for a hodograph-only plot...

`spy.metpy_hodograph(clean_data)`


<div align="center">
<img src="https://raw.githubusercontent.com/kylejgillett/sounderpy/main/images/sounderpy_v2-0-6_example-hodograph.png" width="600">
</div>

------

## AUTHORS AND CONTRIBUTORS
### **AUTHOR: Kyle J Gillett, Central Michigan University** 
#### CONTRIBUTOR: Scott Thomas, NWS Grand Rapids 

------


## CITING SOUNDERPY
in AMS format:

- Gillett, K., 2023: SounderPy: Vertical Profile Data Retrieval & Analysis Tool for Python (Version 2.0.6). Py-Pi, https://pypi.org/project/sounderpy/

------



## REFERENCES
- Harris, C.R., Millman, K.J., van der Walt, S.J. et al. Array programming with NumPy. Nature 585, 357–362 (2020). DOI: 10.1038/s41586-020-2649-2.
  
- Hoyer, S. & Hamman, J., (2017). xarray: N-D labeled Arrays and Datasets in Python. Journal of Open Research Software. 5(1), p.10. DOI: https://doi.org/10.5334/jors.148
  
- J. D. Hunter, "Matplotlib: A 2D Graphics Environment", Computing in Science & Engineering, vol. 9, no. 3, pp. 90-95, 2007.

- Ryan M. May, Sean C. Arms, Patrick Marsh, Eric Bruning, John R. Leeman, Kevin Goebbert, Jonathan E. Thielen, Zachary S Bruick, and M. Drew. Camron. Metpy: a Python package for meteorological data. 2023. URL: Unidata/MetPy, doi:10.5065/D6WW7G29.

- Ryan M. May, Sean C. Arms, John R. Leeman, and Chastang, J. Siphon: A collection of Python Utilities for Accessing Remote Atmospheric and Oceanic Datasets. Unidata. 2017. [Available online at https://github.com/Unidata/siphon.] doi:10.5065/D6CN72NW.

- Pauli Virtanen, Ralf Gommers, Travis E. Oliphant, Matt Haberland, Tyler Reddy, David Cournapeau, Evgeni Burovski, Pearu Peterson, Warren Weckesser, Jonathan Bright, Stéfan J. van der Walt, Matthew Brett, Joshua Wilson, K. Jarrod Millman, Nikolay Mayorov, Andrew R. J. Nelson, Eric Jones, Robert Kern, Eric Larson, CJ Carey, İlhan Polat, Yu Feng, Eric W. Moore, Jake VanderPlas, Denis Laxalde, Josef Perktold, Robert Cimrman, Ian Henriksen, E.A. Quintero, Charles R Harris, Anne M. Archibald, Antônio H. Ribeiro, Fabian Pedregosa, Paul van Mulbregt, and SciPy 1.0 Contributors. (2020) SciPy 1.0: Fundamental Algorithms for Scientific Computing in Python. Nature Methods, 17(3), 261-272.

<div align="center">
<img src="https://github.com/kylejgillett/sounderpy/assets/100786530/2e9477c9-e36a-4163-accb-fe46780058dd" width="250">

</div>

# SounderPy, the vertical profile data retrieval and analysis tool for Python
LATEST VERSION: v3.0.0 |  RELEASED: January, 2024  |  COPYRIGHT Kyle J Gillett, 2023, 2024

A Python package that helps you to access and plot vertical profile data for meteorological analysis 

[![PyPI Package](https://img.shields.io/pypi/v/sounderpy.svg)](https://pypi.python.org/pypi/sounderpy/)
[![PyPI Downloads](https://img.shields.io/pypi/dm/sounderpy.svg)](https://pypi.python.org/pypi/sounderpy/)
[![PyPI license](https://img.shields.io/pypi/l/ansicolortags.svg)](https://github.com/kylejgillett/sounderpy/blob/main/LICENSE.txt)
[![PyPI pyversions](https://img.shields.io/pypi/pyversions/sounderpy.svg)](https://pypi.python.org/pypi/sounderpy/)
[![GitHub commits](https://badgen.net/github/commits/kylejgillett/sounderpy)](https://GitHub.com/kylejgillett/sounderpy/commit/)
[![Maintainer](https://img.shields.io/badge/maintainer-kylejgillett-blue)](https://github.com/kylejgillett)
[![made-with-python](https://img.shields.io/badge/Made%20with-Python-1f425f.svg)](https://www.python.org/)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10011851.svg)](https://doi.org/10.5281/zenodo.10011851)

-----

<div align="center">
   
### [VISIT SOUNDERPY DOCUMENTATION](https://kylejgillett.github.io/sounderpy/) | [CHECK OUT AN EXAMPLES & TUTORIALS](https://kylejgillett.github.io/sounderpy/examplescripts.html)
</div>


## What is SounderPy:

SounderPy is an open-source atmospheric science Python package for vertical profile analysis. This tool is designed to get data, ‘clean it up’ for simple use, and plot the data on advanced-sounding plots. SounderPy was developed with the goal in mind to keep the code simple and efficient for users of all experience levels and for reliability in all use cases.

SounderPy has been used by several institutions. For example, this tool has been implemented by the Des Moines National Weather Service Office, the State University of New York at Albany, Mississippi State University, and others. Many have used SounderPy in projects and papers, such as students at Ohio State University, Central Michigan University & Rizal Technological University.



## Why SounderPy?
- Sometimes data is tough to find, and often times it’s even tougher to get it in the format you like. SounderPy gets you this data!
- The code needed for loading and parsing meteorological data, especially from models, can be large and messy. SounderPy keeps it hidden away in a PyPi package – just import and call sounderPy functions to keep your code clean!
- SounderPy functions are designed to be simple and quick making for reliable use in research, forecast/analysis operations, and simply for fun!



## What kind of data?:

| **DATA**                        | **FUNCTION**       | **TYPE**          | **TIME RANGE**   |
|---------------------------------|--------------------|-------------------|-------------------|
| ECMWF CDS ERA5 reanalysis*      | get_model_data()   | Reanalysis        | 1940-present      |
| UNIDATA THREDDS TDS RAP         | get_model_data()   | Reanalysis        | 2005-present      |
| UNIDATA THREDDS TDS RUC         | get_model_data()   | Reanalysis        | 2005-2020         |
| UNIDATA THREDDS NCEP-FNL        | get_model_data()   | Reanalysis        | 2005-2020         |
| ISU's BUFKIT archive             | get_bufkit_data()  | Model Forecast    | 2011-present      |
| PSU's BUFKIT feed               | get_bufkit_data()  | Model Forecast    | Most recent runs  |
| UNIDATA THREDDS TDS RAP         | get_model_data()   | Model Analysis    | Most recent run   |
| OU ACARS Archive                | acars_data()       | Observations      | 2019-present      |
| The Unv. of WY RAOB Archive      | get_obs_data()     | Observations      | 1973-present      |
| IGRAv2 Observation archive       | get_obs_data()     | Observations      | 1905-present      |


-------

## Installation

1. #### Install the SounderPy software via pip:
   ```
   pip install sounderpy
   ```
   Find it at https://pypi.org/project/sounderpy/
   
3. #### Import SounderPy into your Python project:
   ```py
   import sounderpy as spy
   ```
   
5. #### Lets declare a few simple variables we can use to get data:
   ```py
   year  = '2014' 
   month = '06'
   day   = '16'
   hour  = '18'
   station = 'OAX'
   ```
7. #### Get some data!
   ```py
   # this will get us 18z observations on June 16th, 2014 from OAX (Omaha, Neb)
   clean_data = spy.get_obs_data(station, year, month, day, hour)
   ```
   and boom! Now you have a callable dictionary of vertical profile reanalysis data including... 

   + Temperature
   + Dewpoint
   + Pressure
   + Height 
   + U-component Wind 
   + V-component Wind
   
### SounderPy can also plot profile data on unique sounding and hodograph figures!
   ```py 
   spy.build_sounding(clean_data, color_blind=True)
   ```
   
   ```py 
   spy.build_hodograph(clean_data, dark_mode=True)
   ```
   <div align="center">
   <img src="https://kylejgillett.github.io/sounderpy/_images/example-sounding_light.png" width="600">
   
   <img src="https://kylejgillett.github.io/sounderpy/_images/example-hodograph_dark.png" width="600">
   </div>

## To learn more about what you can do with SounderPy, [check out the documentation](https://kylejgillett.github.io/sounderpy/)
------

## AUTHORS AND CONTRIBUTORS
### **AUTHOR: Kyle J Gillett, Central Michigan University** 
##### CONTRIBUTOR: Scott Thomas, NWS Grand Rapids 

------


## CITING SOUNDERPY
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10011851.svg)](https://doi.org/10.5281/zenodo.10011851)


in AMS format:

- Gillett, K., 2023: SounderPy: Vertical Profile Data Retrieval & Analysis Tool for Python (Version 3.0.0). Py-Pi, https://pypi.org/project/sounderpy/

------



## REFERENCES
- Harris, C.R., Millman, K.J., van der Walt, S.J. et al. Array programming with NumPy. Nature 585, 357–362 (2020). DOI: 10.1038/s41586-020-2649-2.
  
- Hoyer, S. & Hamman, J., (2017). xarray: N-D labeled Arrays and Datasets in Python. Journal of Open Research Software. 5(1), p.10. DOI: https://doi.org/10.5334/jors.148
  
- J. D. Hunter, "Matplotlib: A 2D Graphics Environment", Computing in Science & Engineering, vol. 9, no. 3, pp. 90-95, 2007.

- Ryan M. May, Sean C. Arms, Patrick Marsh, Eric Bruning, John R. Leeman, Kevin Goebbert, Jonathan E. Thielen, Zachary S Bruick, and M. Drew. Camron. Metpy: a Python package for meteorological data. 2023. URL: Unidata/MetPy, doi:10.5065/D6WW7G29.

- Ryan M. May, Sean C. Arms, John R. Leeman, and Chastang, J. Siphon: A collection of Python Utilities for Accessing Remote Atmospheric and Oceanic Datasets. Unidata. 2017. [Available online at https://github.com/Unidata/siphon.] doi:10.5065/D6CN72NW.

- Pauli Virtanen, Ralf Gommers, Travis E. Oliphant, Matt Haberland, Tyler Reddy, David Cournapeau, Evgeni Burovski, Pearu Peterson, Warren Weckesser, Jonathan Bright, Stéfan J. van der Walt, Matthew Brett, Joshua Wilson, K. Jarrod Millman, Nikolay Mayorov, Andrew R. J. Nelson, Eric Jones, Robert Kern, Eric Larson, CJ Carey, İlhan Polat, Yu Feng, Eric W. Moore, Jake VanderPlas, Denis Laxalde, Josef Perktold, Robert Cimrman, Ian Henriksen, E.A. Quintero, Charles R Harris, Anne M. Archibald, Antônio H. Ribeiro, Fabian Pedregosa, Paul van Mulbregt, and SciPy 1.0 Contributors. (2020) SciPy 1.0: Fundamental Algorithms for Scientific Computing in Python. Nature Methods, 17(3), 261-272.

- Marsh, P., Halbert, K., Blumberg, G., Supinie, T., Esmaili, R., Szkodzinski, J., "SHARPpy: Sounding/Hodograph Analysis and Research Program in Python." GitHub. Available at: https://github.com/sharppy/SHARPpy.
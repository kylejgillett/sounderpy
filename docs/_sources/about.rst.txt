About
======


Installation
------------

SounderPy is available on PyPi and can be installed via ``pip``:

.. code-block:: console

   pip install sounderpy 

In your Python document, its fun to import SounderPy as ``spy``!:

.. code-block:: python

	import sounderpy as spy



***************************************************************

Sample Basic Use 
-----------------

SounderPy is designed for simple and efficient use. Below is a basic example for plotting the latest HRRR forecast profile at hour 4 for BUFKIT site 'KMOP':

.. code-block:: python

   import sounderpy as spy
   clean_data = spy.get_bufkit_data('hrrr', 'kmop', 4)
   spy.build_sounding(clean_data, style='simple', dark_mode=True)

Those three lines will make this!:

.. image:: _static/images/example-sounding_dark.png
   :alt: Example SounderPy Sounding





***************************************************************



Authors and Contributors 
-------------------------
	**AUTHOR: Kyle J Gillett, Central Michigan University**

	*CONTRIBUTOR: Scott Thomas, NWS Grand Rapids*



Citing SounderPy
-----------------
	.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.10011851.svg
	   :target: https://doi.org/10.5281/zenodo.10011851
	   :alt: DOI

	in AMS format:
	    * Gillett, K., 2023: SounderPy: Vertical Profile Data Retrieval & Analysis Tool for Python (Version 2.0.6). Py-Pi, https://pypi.org/project/sounderpy/


Refrences 
----------

	* Harris, C.R., Millman, K.J., van der Walt, S.J. et al. Array programming with NumPy. Nature 585, 357–362 (2020). DOI: 10.1038/s41586-020-2649-2.


	* Hoyer, S. & Hamman, J., (2017). xarray: N-D labeled Arrays and Datasets in Python. Journal of Open Research Software. 5(1), p.10. DOI: https://doi.org/10.5334/jors.148


	* J. D. Hunter, "Matplotlib: A 2D Graphics Environment", Computing in Science & Engineering, vol. 9, no. 3, pp. 90-95, 2007.


	* Ryan M. May, Sean C. Arms, Patrick Marsh, Eric Bruning, John R. Leeman, Kevin Goebbert, Jonathan E. Thielen, Zachary S Bruick, and M. Drew. Camron. Metpy: a Python package for meteorological data. 2023. URL: Unidata/MetPy, doi:10.5065/D6WW7G29.


	* Ryan M. May, Sean C. Arms, John R. Leeman, and Chastang, J. Siphon: A collection of Python Utilities for Accessing Remote Atmospheric and Oceanic Datasets. Unidata. 2017. [Available online at https://github.com/Unidata/siphon.] doi:10.5065/D6CN72NW.


	* Pauli Virtanen, Ralf Gommers, Travis E. Oliphant, Matt Haberland, Tyler Reddy, David Cournapeau, Evgeni Burovski, Pearu Peterson, Warren Weckesser, Jonathan Bright, Stéfan J. van der Walt, Matthew Brett, Joshua Wilson, K. Jarrod Millman, Nikolay Mayorov, Andrew R. J. Nelson, Eric Jones, Robert Kern, Eric Larson, CJ Carey, İlhan Polat, Yu Feng, Eric W. Moore, Jake VanderPlas, Denis Laxalde, Josef Perktold, Robert Cimrman, Ian Henriksen, E.A. Quintero, Charles R Harris, Anne M. Archibald, Antônio H. Ribeiro, Fabian Pedregosa, Paul van Mulbregt, and SciPy 1.0 Contributors. (2020) SciPy 1.0: Fundamental Algorithms for Scientific Computing in Python. Nature Methods, 17(3), 261-272.



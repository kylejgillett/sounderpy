📖 About
=========

☕ SounderPy is an open-source package developed on my own time. Would you like to support continued SounderPy development? Consider "`Buying me a coffee <https://www.buymeacoffee.com/kylejgillett>`_"! ☕

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

SounderPy is designed for simple and efficient use. Below is a basic example plotting an 03/31/2023 12z HRRR forecast profile at forecast hour 8 for BUFKIT site 'KMLI':

.. code-block:: python

   import sounderpy as spy
   clean_data = spy.get_bufkit_data('hrrr', 'kmli', 8, '2023', '03', '31', '12')
   spy.build_sounding(clean_data, style='simple', dark_mode=True)

Those three lines will make this!:

.. image:: _static/images/example-sounding_simple_dark.png
   :alt: Example SounderPy Sounding


More examples of these plots are available :ref:`here<gallery>`


***************************************************************



Authors and Contributors 
-------------------------
	**AUTHOR: Kyle J Gillett, Central Michigan University**

	*CONTRIBUTOR: Scott Thomas, NWS Grand Rapids*



Citing SounderPy
-----------------
	.. image:: https://zenodo.org/badge/DOI/10.5281/zenodo.10443609.svg
	   :target: https://doi.org/10.5281/zenodo.10443609
	   :alt: DOI


	in AMS format:
	     Gillett, K., 2024: SounderPy: Vertical Profile Data Retrieval & Analysis Tool for Python (Version 3.0.1). Py-Pi, https://pypi.org/project/sounderpy/


***************************************************************


References 
----------

	* Harris, C.R., Millman, K.J., van der Walt, S.J. et al. Array programming with NumPy. Nature 585, 357–362 (2020). DOI: 10.1038/s41586-020-2649-2.
      

	* Hoyer, S. & Hamman, J., (2017). xarray: N-D labeled Arrays and Datasets in Python. Journal of Open Research Software. 5(1), p.10. DOI: https://doi.org/10.5334/jors.148

       
	* J. D. Hunter, "Matplotlib: A 2D Graphics Environment", Computing in Science & Engineering, vol. 9, no. 3, pp. 90-95, 2007.

      
	* Ryan M. May, Sean C. Arms, Patrick Marsh, Eric Bruning, John R. Leeman, Kevin Goebbert, Jonathan E. Thielen, Zachary S Bruick, and M. Drew. Camron. Metpy: a Python package for meteorological data. 2023. URL: Unidata/MetPy, doi:10.5065/D6WW7G29.
      

	* Ryan M. May, Sean C. Arms, John R. Leeman, and Chastang, J. Siphon: A collection of Python Utilities for Accessing Remote Atmospheric and Oceanic Datasets. Unidata. 2017. [Available online at https://github.com/Unidata/siphon.] doi:10.5065/D6CN72NW.
      

	* Pauli Virtanen, Ralf Gommers, Travis E. Oliphant, Matt Haberland, Tyler Reddy, David Cournapeau, Evgeni Burovski, Pearu Peterson, Warren Weckesser, Jonathan Bright, Stéfan J. van der Walt, Matthew Brett, Joshua Wilson, K. Jarrod Millman, Nikolay Mayorov, Andrew R. J. Nelson, Eric Jones, Robert Kern, Eric Larson, CJ Carey, İlhan Polat, Yu Feng, Eric W. Moore, Jake VanderPlas, Denis Laxalde, Josef Perktold, Robert Cimrman, Ian Henriksen, E.A. Quintero, Charles R Harris, Anne M. Archibald, Antônio H. Ribeiro, Fabian Pedregosa, Paul van Mulbregt, and SciPy 1.0 Contributors. (2020) SciPy 1.0: Fundamental Algorithms for Scientific Computing in Python. Nature Methods, 17(3), 261-272.

       
	* Marsh, P., Halbert, K., Blumberg, G., Supinie, T., Esmaili, R., Szkodzinski, J., "SHARPpy: Sounding/Hodograph Analysis and Research Program in Python." GitHub. Available at: https://github.com/sharppy/SHARPpy.


****************************************

About the Author
-----------------

Hey! 

Thank you so much for checking out and using SounderPy. My name is Kyle Gillett and I am the developer of SounderPy. This tool started out as a way for me to internally house all of my data retrieval functions for plotting soundings. As you can see, it has since blossomed into a full-scale Python package. I am currently an undergraduate student of Meteorology at Central Michigan University, the President of the Central Michigan University Student Chapter of the American Meteorological Society and an on-air meteorologist with WNEM TV-5 in Saginaw, MI.

SounderPy is published on PyPi and the source code is available on GitHub -- this tool is an open source project. If you have found SounderPy useful in your work, I'd love to hear about it! The coolest part of this project has been hearing how many folks have been using this software. If you'd like to support continued SounderPy development, consider "`Buying me a coffee <https://www.buymeacoffee.com/kylejgillett>`_"! ☕. 

*Have an issue?* You can `open a GitHub Issue <https://github.com/kylejgillett/sounderpy/issues>`_ or just shoot me a DM on `Twitter <https://twitter.com/wxkylegillett>`_!

**Useful Links**

+ Check out SounderPy `on GitHub <https://github.com/kylejgillett/sounderpy>`_
+ Check out SounderPy `on PyPi <https://pypi.org/project/sounderpy/>`_
+ Check out my `website <https://kylegillettwx.wordpress.com/>`_
+ Get updates on SounderPy development on `Twitter <https://twitter.com/wxkylegillett>`_
+ Support SounderPy by "`Buying me a coffee <https://www.buymeacoffee.com/kylejgillett>`_"


Thanks for using SounderPy!

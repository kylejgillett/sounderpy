Tools for Plotting Data
========================

.. image:: _static/images/example-sounding_light.png
   :width: 300 px
   :align: right


Initially, SounderPy started as a tool meant for getting and parsing data with little focus on plotting said data. However, recent releases have put more emphasis on plotting capabilites. Version 3.0.0 features a number of significant & exciting upgrades to SounderPy's plotting abilities. 

SounderPy can create general :ref:`sounding<soundings>` and :ref:`hodograph<hodographs>` plots as well as :ref:`composite sounding<compsoundings>` plots! 

The `full` sounding plots that SounderPy creates are complex figures with unique design geared towards severe convective storm enviroment analysis. 

These sounding plots are 'my baby' and I hope you find them useful! :)

***************************************************************




.. _soundings:

Building Soundings
----------------------------------

We can use the simple ``spy.build_sounding()`` function:

.. py:function:: spy.build_sounding(clean_data, style='full', save=False, filename='sounderpy_sounding', color_blind=True, dark_mode=False)


   Return a full sounding plot of SounderPy data, ``plt`` 

   :param clean_data: the dictionary of data to be plotted (see :doc:`gettingdata`)
   :type clean_data: dict, required
   :param style: may be `simple` or `full`. Default is `full`.
   :type style: str, optional
   :param save: whether to show the plot inline or save to a file. Default is ``False`` which displays the file inline.
   :type save: bool, optional
   :param filename: the filename by which a file should be saved to if ``save = True``. Default is `sounderpy_sounding`.
   :type filename: str, optional
   :param color_blind: whether or not to change the dewpoint trace line from green to blue for improved readability for color deficient users/readers. Default is ``False``
   :type color_blind: bool, optional
   :param dark_mode: ``True`` will invert the color scheme for a 'dark-mode' sounding. Default is ``False``.
   :type dark_mode: bool, optional
   :return: plt, a SounderPy sounding built with Matplotlib, MetPy, SharpPy, & SounderPy.
   :rtype: plt

Examples
^^^^^^^^^^

.. code-block:: python

   import sounderpy as spy
     
   # get data | Note: any sounderpy data will work!
   clean_data = spy.get_obs_data('OAX', '2014', '06', '16', '18')

   # build the sounding! 
   spy.build_sounding(clean_data)

.. image:: _static/images/example-sounding_light.png
   :width: 800 px

**************************************************







.. _hodographs:

Building Hodographs
----------------------------------

Very similarly to soundings, we can use the simple ``spy.build_hodograph()`` function:

.. py:function:: spy.build_hodograph(clean_data, save=False, filename='sounderpy_sounding', dark_mode=False)


   Return a full sounding plot of SounderPy data, ``plt`` 

   :param clean_data: the dictionary of data to be plotted (see :doc:`gettingdata`)
   :type clean_data: dict, required
   :param save: whether to show the plot inline or save to a file. Default is ``False`` which displays the file inline.
   :type save: bool, optional
   :param filename: the filename by which a file should be saved to if ``save = True``. Default is `sounderpy_sounding`.
   :type filename: str, optional
   :param dark_mode: ``True`` will invert the color scheme for a 'dark-mode' sounding. Default is ``False``.
   :type dark_mode: bool, optional
   :return: plt, a SounderPy sounding built with Matplotlib, MetPy, SharpPy, & SounderPy.
   :rtype: plt



Examples
^^^^^^^^^^

.. code-block:: python

   import sounderpy as spy
     
   # get data | Note: any sounderpy data will work!
   clean_data = spy.get_obs_data('OAX', '2014', '06', '16', '18')

   # build the hodograph! 
   spy.build_hodograph(clean_data)

.. image:: _static/images/example-hodograph_light.png
   :width: 800 px

**************************************************





.. _compsoundings:

Building Composite Soundings
----------------------------------

Sometimes we want to compare two or more profiles against each other. Perhaps at different locations or times, or we may want to compare different models or model run-times. SounderPy allows you to do this!

To do so, a list of `clean_data` dicts is needed. If you want to customize the look of each profile, you can create equal length lists with alphas, linestyles, linewidths, & colors. See below:

.. py:function:: spy.build_composite(data_list, dark_mode=False, shade_between=False,alphas_to_use='none', colors_to_use='none', ls_to_use='none', lw_to_use='none', save=False, filename='sounderpy_sounding')


   Return a composite sounding plot of multiple profiles, ``plt`` 

   :param data_list: a list of data dictionaries for each profile to be plotted
   :type data_list: list of dicts, required
   :param dark_mode: ``True`` will invert the color scheme for a 'dark-mode' sounding. Default is ``False``.
   :type dark_mode: bool, optional
   :param shade_between: Lightly shade between the dewpoint & temperature trace. In many cases, this improves readability. Default is ``True``.
   :type shade_between: bool, optional
   :param alphas_to_use: A list of custom alphas (0.0-1.0). List length must match the number of profiles listed in ``data_list``. Default is 'none'. Default alpha is 1.
   :type alphas_to_use: list of floats, optional
   :param colors_to_use: A list of custom matplotlib color name stings. List length must match the number of profiles listed in ``data_list``. Default is 'none'.
   :type colors_to_use: list of strings, optional
   :param ls_to_use: A list of custom matplotlib linestyles. List length must match the number of profiles listed in ``data_list``. Default is 'none'. Default linestyle is '-'.
   :type ls_to_use: list of stings, optional
   :param lw_to_use: A list of custom linewidths. List length must match the number of profiles listed in ``data_list``. Default is 'none'. Default linewidth is 3.
   :type lw_to_use: list of floats, optional
   :param save: whether to show the plot inline or save to a file. Default is ``False`` which displays the file inline.
   :type save: bool, optional
   :param filename: the filename by which a file should be saved to if ``save = True``. Default is `sounderpy_sounding`.
   :type filename: str, optional
   :return: plt, a SounderPy composite sounding built with Matplotlib, MetPy, SharpPy, & SounderPy.
   :rtype: plt


Examples
^^^^^^^^^^

.. code-block:: python

   import sounderpy as spy
     
   # get data | Note: any sounderpy data will work!
   # this example looks at 3 profiles from OAX on Pilger-day.
   clean_data1 = spy.get_obs_data('oax', '2014', '06', '16', '12')
   clean_data2 = spy.get_obs_data('oax', '2014', '06', '16', '18')
   clean_data3 = spy.get_obs_data('oax', '2014', '06', '17', '00')
     
   # add each dict of data to a list
   data_list = [clean_data1, clean_data2, clean_data3]

   # build the composite! 
   spy.build_composite(data_list)

.. image:: _static/images/example-composite_light.png
   :width: 800 px

.. code-block:: python

   import sounderpy as spy
     
   # get data | Note: any sounderpy data will work!
   # this example looks at 3 profiles from OAX on Pilger-day.
   clean_data1 = spy.get_obs_data('oax', '2014', '06', '16', '12')
   clean_data2 = spy.get_obs_data('oax', '2014', '06', '16', '18')
   clean_data3 = spy.get_obs_data('oax', '2014', '06', '17', '00')
     
   # add each dict of data to a list
   data_list = [clean_data1, clean_data2, clean_data3]

   # lets add some custom parameters:
   ls_to_use     = ['-', '--', ':']
   lw_to_use     = [1, 2, 3]
   alphas_to_use = [0.6, 0.8, 1]
   colors_to_use = ['cornflowerblue', 'orange', 'lime'] 

   # build the composite with our custom parameters!
   # and make it dark-mode for fun! 
   spy.build_composite(data_list, ls_to_use=ls_to_use, lw_to_use=lw_to_use
   alphas_to_use=alphas_to_use, colors_to_use=colors_to_use, dark_mode=True)

.. image:: _static/images/example-composite_dark.png
   :width: 800 px

**************************************************








Printing data to the console
-----------------------------

.. py:function:: spy.print_variables(clean_data):

   :param clean_data: the dictionary of data to be plotted (see :doc:`gettingdata`)
   :type clean_data: dict, required
	:return: prints a number of thermodynamic and kinematic variables to the console.


.. code-block:: python

   > THERMODYNAMICS --------------------------------------------- 
   --- SBCAPE: 5771.5 | MUCAPE: 5771.5 | MLCAPE: 4322.8 | ECAPE: 5136.3
   --- MU 0-3: 111.7 | MU 0-6: 1174.7 | SB 0-3: 111.7 | SB 0-6: 1174.7
    
   > KINEMATICS ------------------------------------------------- 
   --- 0-500 SRW: 34.8 knot | 0-500 SWV: 0.014 | 0-500 SHEAR: 16.6 | 0-500 SRH: 153.8
   --- 1-3km SRW: 30.3 knot | 1-3km SWV: 0.008 | 1-3km SHEAR: 35.6 | | 1-3km SRH: 234.2

**************************************************





About These Plots
-----------------

This plot style has been developed in a way that acts to provide as much information to the user as possible with attributes designed specifically for the analysis of severe convective environments, and supercells/tornadoes in particular. You will also find that this particular plot style does not host many of the common and popular severe weather composite indices – that was intentional. Most, if not all, of the data provided on this plot, are considered, for the lack of a better word, ‘true’ observations of the atmosphere though most are still subject to heavy assumptions. 

The data on these plots are considered, by most, to be useful in determining critical characteristics of the atmosphere related to mesoscylonegenesis and tornadogenesis.
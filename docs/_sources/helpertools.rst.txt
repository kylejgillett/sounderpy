🛠️ Helper Tools
================

A collection of helper tools included within SounderPy.

*************************

Printing a dictionary of profile parameters
--------------------------------------------

Return a massive dictionary of common calculated paramters and special variables for a vertical profile.

Data returned in this dictionary include...
      - an interpolated version of the profile
      - basic parameters such as mixing ratio, theta-e, wet-bulb, etc
      - Thermodynamic parameters such as SB, MU, ML parcel properties (LCL, LFC, EL, CAPE, CIN, etc)
      - Kinematic parameters such as storm motion, storm-relative wind, streamwise vorticity, SRH, etc.

.. class:: sounding_params()

      :param clean_data: a dictionary of SounderPy data (see :doc:`gettingdata`)
      :type clean_data: dict, required

   .. py:function:: .calc()

      :return: special sounding parameters
      :rtype: dict


.. note::
   These calculated parameters are the same used on SounderPy plots. This function simply returns the master dictionary of all of these values.



****************************************************************

Finding site lat/lon pairs
---------------------------

.. py:function:: spy.get_latlon(station_type, station_id)


   Return a latitude-longitude float pair in a ``list``

   :param station_type: the station 'type' that corresponds with the given station ID
   :type station_type: str, required
   :param station_id: the station ID for the given station type
   :type station_id: str, required
   :return: lat/lon float pair
   :rtype: list

Example:

.. code-block:: python

   spy.get_latlon('metar', 'kmop')
   spy.get_latlon('bufkit', 'apx')
   spy.get_latlon('raob', 'oun') 
   spy.get_latlon('buoy', '45210')


* note: you can use this lat/lon pair list when calling the function :ref:`get_model_data<modeldata>`


***************************************************************

Saving data to a file
----------------------

.. py:function:: spy.to_file(file_type, clean_data, filename=None)


   Create a file of 'cleaned' SounderPy data

   :param file_type: a `str` representing the file type you'd like to export data to.
   :type file_type: str, required
   :param clean_data: 'cleaned' SounderPy data `dict`
   :type clean_data: dict, required
   :param filename: the name you'd like to give the file
   :type filename: str, required
   :return: a file of SounderPy data.

Example:
	* File options include `csv`, `cm1`, & `sharppy`

.. code-block:: python

   spy.to_file('csv', clean_data)
   spy.to_file('cm1', clean_data)
   spy.to_file('sharppy', clean_data)  


***************************************************************

Interpolating a vertical profile
---------------------------------

.. py:function:: spy.interp_data(variable, heights, step=100)


   Interpolate a 1D array of data (such as a temperature profile) over a given interval (step) based on a corresponding array of height values. 

   :param variable: an array of data to be interpolated. Must be same length as height array.
   :type variable: arr, required
   :param heights: heights corresponding to the vertical profile used to interpolate. Must be same length as variable array.
   :type heights: arr, required
   :param step: the resolution of interpolation. Default is 100 (recommended value is 100)
   :type step: int, optional
   :return: interp_var, an array of interpolated data.
   :rtype: arr

Example:

.. code-block:: python

   spy.interp_data(temperature_array, height_array, step=100)  


***************************************************************

Finding a 'nearest' value
--------------------------

.. py:function:: spy.find_nearest(array, value)


	Return a value of an index of an array who's value is closest to a define value.

   :param array: an array of data to be searched through
   :type array: arr, required
   :param heights: the value used to compare against the array of data
   :type heights: int or float, required
   :return: nearest_idx, index of the data array that corresponds with the nearest value to the given value
   :rtype: int

Example:

.. code-block:: python

   z_equals_500m = spy.interp_data(z, 500)
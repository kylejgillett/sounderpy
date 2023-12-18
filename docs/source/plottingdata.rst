Potting Data with SounderPy
===========================



***************************************************************


Building Soundings
----------------------------------

We can use the simple ``spy.build_sounding()`` function:

.. py:function:: spy.build_sounding(clean_data, style='full', save=False, filename='sounderpy_sounding', color_blind=True, dark_mode=False)


   Return a full sounding plot of SounderPy data, ``plt`` 

   :param model: required, the requested model to use (rap-ruc, era5, ncep)
   :type model: str
   :param latlon: required, the latitude & longitude pair for sounding ([44.92, -84.72])
   :type latlon: list
   :param year: required, valid year
   :type year: str
   :param month: required, valid month
   :type month: str
   :param day: required, valid day
   :type day: str
   :param hour: required, valid hour
   :type hour: str
   :return: clean_data, a dict of ready-to-use vertical profile data including pressure, height, temperature, dewpoint, u-wind, v-wind, & model information
   :rtype: dict


Bulding Hodographs
----------------------------------

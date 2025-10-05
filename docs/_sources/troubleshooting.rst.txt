❓ Troubleshooting
===================


Like any software, SounderPy is prone to errors now and then -- some of them are known with fixes! If you are running into issues using SounderPy, check out this guide to troubleshooting below! 

If you have an error that you don't see on this page, you can `open a GitHub Issue <https://github.com/kylejgillett/sounderpy/issues>`_ or feel free to shoot me a DM on `Twitter <https://twitter.com/wxkylegillett>`_!


Testing SounderPy
-------------------

You can test your installation of SounderPy against this figure. To do so, run the following block of code. The figure produced should be the same.

.. image:: _static/images/test-figure_abr-2024-08-29-00z.png
   :width: 300 px
   :align: center


.. code-block:: python

    import sounderpy as spy

    data = spy.get_obs_data("kabr", "2024", "08", "29", "00")

    spy.build_sounding(data, color_blind=True, radar='single')




***************************************



Funky Plot Scaling 
-------------------

.. image:: _static/images/funky_layout_example.jpg
   :width: 300 px
   :align: left

When using SounderPy's plot function, sometimes the matplotlib formatting operations the plot functions use (``build_sounding()``, ``build_hodograph()``, ``build_composite()``) don't agree with some screen sizes. This throws off the scaling of the plots when plotting them 'inline'. Thankfully there is a fix! 

Each of SounderPy's plotting functions have a ``save`` key-word-argument (kwarg). By default these are set to *False*, but when setting the kwarg to *True* you can save the plot to a .png file. 

When the plot is saved to a file it will preserve the correct scaling and dimensions. 

The Fix: saving the plot to a file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	Add the ``save=True`` kwarg to your plot function (``build_sounding()``, ``build_hodograph()``, ``build_composite()``). 

	You can optionally add the ``filename`` kwarg that allows you to set a specific filename. By default, this filename is *'sounderpy_sounding.png'*.







***************************************

Goofy ACARS Profiles
---------------------

.. image:: _static/images/example-goofy-acars.png
   :width: 300 px
   :align: left

Sometimes ACARS data just comes with errors that can not be fixed. Because these observations are taken by aircraft during arrivals and departures, they can be prone to errors such as too-dry dewpoints, too-strong wind velocity, and strange height data. Unfortunately these errors can not be fixed.

The best way to check if ACARS observations are correct is by comparing them to recent forecast model profiles (at forecast hour 0) for that location.






***************************************

ERA5 data access not working
-----------------------------

To access ERA5 data, some extra steps are involved. Make sure you create an account with the ECMWF Climate Data Store (CDS) to gain access to their API.

Also note that CDS has recently updated their API (September 2024). Thus, even if you have done this before, you have to migrate to the new API.

Basic steps include:
    - creating a CDS API account
    - Setting up a CDS API personal access token
    - Creating a $HOME/.cdsapirc file

   Follow the instructions on the CDSAPI "how to" documentation -- See: https://cds.climate.copernicus.eu/how-to-api






******************************************

RAP data won't load
-------------------

Sometimes, this just happens. The NCEP server (THREDDS) is prone to high latency or even simply being down. Therefore, data access may take several minutes or data access is not possible.

Sometimes these issues are resolved within hours, other times possibly a few days. Trying again at a later time is generally the best fix.






***************************************

SounderPy can't find a RAOB sounding you know exists
----------------------------------------------------

Ah, this one is fun. Different databases store RAOB data differently. As such, you may find that a sounding that you believe was launched at '18z' is actually stored as '19z' or '17z'. Its annoying.

This will cause the `get_obs_data()` to "not find" the sounding you want because the database it's searching through has it listed under a different time.

Try increasing or decreasing the launch_hour by 1 hour.


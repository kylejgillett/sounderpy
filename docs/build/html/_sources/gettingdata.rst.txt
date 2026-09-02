🌐 Tools for Getting Data
===========================

SounderPy provides a common interface for retrieving atmospheric vertical
profiles from observations, forecasts, analyses, and reanalyses.

Regardless of source, supported retrieval functions return SounderPy's
``clean_data`` structure, allowing the same plotting and analysis workflow
to be used with each dataset.

See :ref:`clean_data schema <clean_data_scheme>` for details on the returned
data structure.


Available Data Sources
-----------------------

.. list-table::
   :header-rows: 1
   :widths: 28 25 22 25

   * - Data source
     - Function
     - Type
     - Availability
   * - ERA5
     - ``get_model_data()``
     - Reanalysis
     - 1940–present
   * - RAP / RUC
     - ``get_model_data()``
     - Analysis / reanalysis
     - Historical archive
   * - NCEP-FNL
     - ``get_model_data()``
     - Analysis
     - Historical archive
   * - BUFKIT
     - ``get_bufkit_data()``
     - Model forecast
     - Recent + archived runs
   * - RAOB / IGRAv2
     - ``get_obs_data()``
     - Observation
     - Historical archive
   * - ACARS
     - ``acars_data()``
     - Observation
     - Archive-dependent
                    

**********************************************************************

.. _modeldata:

Model Reanalysis Data | RAP/RUC, ERA5, NCEP 
---------------------------------------------

SounderPy provides ``get_model_data()`` for retrieving vertical profiles from
RAP/RUC analyses, ERA5 reanalysis, and NCEP-FNL analysis data.

.. important::
   **RAP/RUC access was modernized in SounderPy 3.2.0.**

   Modern RAP data are retrieved from NOAA's AWS GRIB2 archive, while older
   RAP/RUC data are retrieved through the NCEI historical GRIB archive.

   Availability of historical NCEI data may occasionally be affected by
   upstream archive outages.

- The retrieval workflow accesses pressure-level and near-surface model fields, samples them over a small geographic box surrounding the requested latitude/longitude, and converts the resulting vertical profile into SounderPy's ``clean_data`` format.

- The size of the averaging area is controlled with ``box_avg_size``. See :ref:`clean_data schema <clean_data_scheme>` for details on the returned profile structure.

- For RAP/RUC requests, SounderPy automatically determines the appropriate archive and dataset for the requested date and time. The optional ``dataset``  argument may be used to target a specific historical RAP/RUC dataset rather than allowing SounderPy to search automatically. This option applies to archive datasets that expose dataset identifiers; NOAA AWS RAP retrieval does not use dataset keys.

We can use the simple ``spy.get_model_data()`` function:

.. py:function:: spy.get_model_data(model, latlon, year, month, day, hour, dataset=None, box_avg_size=0.10, hush=False, clean_it=True)


   Return a ``dict`` of 'cleaned up' model reanalysis data from a given model, for a given location, date, and time

   :param model: the requested model to use ("rap-ruc", "era5", "ncep")
   :type model: str, required
   :param latlon: the latitude & longitude pair for sounding ([44.92, -84.72])
   :type latlon: list, required
   :param year: valid year
   :type year: str, required
   :param month: valid month
   :type month: str, required
   :param day: valid day
   :type day: str, required
   :param hour: required, valid hour
   :type hour: str, required
   :param dataset: target a specific dataset instead of searching for the first one with data ("rap-ruc" only).
   :type dataset: str, optional, default is `None`
   :param box_avg_size: Width, in degrees, of the geographic box used to spatially average gridded model data around the requested location.
   :type box_avg_size: float, optional, default is ``0.10``
   :param hush: whether to 'hush' a read-out of thermodynamic and kinematic parameters when getting data.
   :type hush: bool, optional, default is `False`
   :param clean_it: whether to return the raw_data object or a clean_data dict.
   :type clean_it: bool, optional, default is `True`
   :return: A SounderPy ``clean_data`` dictionary containing the retrieved vertical profile and associated metadata. See :ref:`clean_data schema <clean_data_scheme>`.
   :rtype: dict



.. _reanalysiskeys:

Model key names 
^^^^^^^^^^^^^^^
  - ``'era5'``: ECMWF renalysis v5 (ERA5), reanalysis

  - ``'rap'``, or ``'rap-ruc'``: NCEP Rapid Refresh model (RAP) / Rapid Update Cycle model (RUC), reanalysis

  - ``'ncep'``: NCEP Global Data Assimilation System/Final 0.25 degree (ncep-fnl), reanalysis

  - ``'rap-now'``: NCEP Rapid Refresh model, latest analysis


.. _datasetkeys:



Dataset key names 
^^^^^^^^^^^^^^^^^^

Most users should leave ``dataset=None`` and allow SounderPy to select the
appropriate archive automatically.

For advanced troubleshooting or reproducibility, a specific historical
RAP/RUC dataset may be selected with ``dataset=``:

  - ``'RAP_25km'``
  - ``'RAP_25km_old'``
      
  - ``'RAP_25km_anl'``
  - ``'RAP_25km_anl_old'``
          
  - ``'RAP_13km'``
  - ``'RAP_13km_old'``
      
  - ``'RAP_13km_anl'``
  - ``'RAP_13km_anl_old'``
      
  - ``'RUC_13km'``
  - ``'RUC_13km_old'``
      
  - ``'RUC_25km'``
  - ``'RUC_25km_old'``


.. _latlonpairs:

Latitude-Longitude pairs
^^^^^^^^^^^^^^^^^^^^^^^^^
Locations are supplied as a two-element list in the form:

.. code-block:: python

   [latitude, longitude]

For example:

.. code-block:: python

   [44.92, -84.72]



.. note::
   **Scientific consideration:** Model analyses and reanalyses are gridded
   representations of the atmosphere and should not be interpreted as direct
   observations. Their spatial resolution, assimilation system, model physics,
   and the ``box_avg_size`` used by SounderPy should be considered when
   interpreting a retrieved profile.

.. tip::
   **Using ERA5 requires ECMWF CDS API access.**

   Before retrieving ERA5 data, you must:

   - create a Climate Data Store account;
   - configure your CDS API personal access token; and
   - create the required ``$HOME/.cdsapirc`` file.

   Follow the official CDS API setup instructions:
   https://cds.climate.copernicus.eu/how-to-api


.. tip::
   **Is data access taking forever?** Sometimes the NCEP (RAP-RUC, NCEP-FNL) & ECMWF CDS (ERA5) servers are down and not able to be accessed. Sometimes these issues are resolved within hours, other times possibly a few days. 


Example
^^^^^^^

Retrieve a RAP analysis near central South Dakota:

.. code-block:: python

   import sounderpy as spy

   data = spy.get_model_data("rap-ruc", [44.58, -100.82], "2024", "08", "28", "18",
                            box_avg_size=0.25, hush=True)

The returned ``data`` object can be passed directly to SounderPy's plotting
and analysis tools:

.. code-block:: python

   spy.build_sounding(data)



***************************************************************




.. _bufkitdata:

Model Forecast Data | BUFKIT Profiles
---------------------------------------

SounderPy provides ``get_bufkit_data()`` for retrieving model forecast
vertical profiles from BUFKIT sites.

.. py:function:: spy.get_bufkit_data(model, station, fcst_hour, run_year=None, run_month=None, run_day=None, run_hour=None, hush=False, clean_it=True)

   Retrieve a BUFKIT model forecast profile for a requested model, site,
   forecast hour, and optionally a specific model run.

   :param model: the model :ref:`'key' <forecastkeys>` name to request data from 
   :type model: str, required
   :param station: a 3-4 digit BUFKIT site identifier
   :type station: str, required
   :param fcst_hour: valid forecast hour
   :type fcst_hour: int, required
   :param run_year: valid year
   :type run_year: str, optional, Default=None
   :param run_month: valid month
   :type run_month: str, optional, Default=None
   :param run_day: valid day
   :type run_day: str, optional, Default=None
   :param run_hour: valid hour
   :type run_hour: str, optional, Default=None
   :param hush: whether to 'hush' a read-out of thermodynamic and kinematic parameters when getting data.
   :type hush: bool, optional, default is `False`
   :param clean_it: whether to return the raw_data object or a clean_data dict.
   :type clean_it: bool, optional, default is `True`
   :return: :ref:`clean_data<datadescription>`, a dict of ready-to-use vertical profile data including pressure, height, temperature, dewpoint, u-wind, v-wind, omega, & model information
   :rtype: dict

.. _forecastsites:

Available BUFKIT Sites:
^^^^^^^^^^^^^^^^^^^^^^^^
.. raw:: html

    <embed>
        <iframe src="https://kylejgillett.github.io/sounderpy/bufkit_map" width="100%" height="500"></iframe>
    </embed>

\ 

.. _forecastmodels:

Available Models:
^^^^^^^^^^^^^^^^^^

**Recent model runs**: Recent runs are retrieved through the Penn State BUFKIT feed.

- GFS
- NAM
- NAMNEST
- RAP
- HRRR
- SREF
- HIRESW


**Archived model runs**: Archived runs are retrieved through the Iowa State BUFKIT archive.

- GFS
- NAM
- NAMNEST
- RAP
- HRRR

\ 

.. _forecastkeys:

Model key names 
^^^^^^^^^^^^^^^
  - ``hrrr``: High Resolution Rapid Refresh, analysis (F00) & forecast; out to forecast hour 48
  - ``rap``: Rapid Refresh Model, analysis (F00) & forecast; out to forecast hour 51
  - ``nam``: North American Mesoscale Model, analysis (F00) & forecast; out to forecast hour 48
  - ``namnest``: Nested North American Mesoscale model, analysis (F00) & forecast; out to forecast hour 60
  - ``gfs``: Global Forecast System, analysis (F00) & forecast; out to forecast hour 180
  - ``sref``: Short Range Ensemble Forecast, analysis (F00) & forecast; out to forecast hour 84
  - ``hiresw``: High Resolution Window Forecast System, analysis (F00) & forecast; out to forecast hour 48


.. tip::
   Running the ``get_bufkit_data()`` function without date kwargs will return the latest available forecast. Example:

   .. code-block:: python
       :linenos:

       # RAP model for site KGFK at forecast hour 5
       spy.get_bufkit_data('rap', 'kgfk', 5)



.. note::
   BUFKIT profiles are model forecast data,
   not observations. Profiles are available only at designated BUFKIT
   locations and represent the model atmosphere at those sites.

   Interpretation should therefore consider the model, initialization time,
   forecast hour, and model resolution.



Example
^^^^^^^

Retrieve a GFS BUFKIT forecast from the 12 UTC 5 August 2023 model run at
``KMOP``, valid at forecast hour 6:

.. code-block:: python

   import sounderpy as spy

   data = spy.get_bufkit_data(
       "gfs",
       "KMOP",
       6,
       "2023", "08", "05", "12",
       hush=True,
   )

The returned profile can be used directly with SounderPy plotting and
analysis tools:

.. code-block:: python

   spy.build_sounding(data)




Latest Available Forecast
^^^^^^^^^^^^^^^^^^^^^^^^^

If the model-run date arguments are omitted, SounderPy retrieves the latest
available BUFKIT forecast. For example:

.. code-block:: python

   data = spy.get_bufkit_data(
       "rap",
       "KGFK",
       5,
   )

This requests the latest available RAP forecast for ``KGFK`` at forecast
hour 5.

***************************************************************



.. _obsdata:

Observed Data | RAOB & IGRAv2 Profiles
---------------------------------------

SounderPy provides ``get_obs_data()`` for retrieving observed radiosonde
profiles from supported RAOB and IGRAv2 archives.

The appropriate archive is selected automatically based on the supplied
station identifier.


.. py:function:: spy.get_obs_data(station, year, month, day, hour, hush=False, clean_it=True)

   Retrieve an observed atmospheric profile for a requested station,
   date, and launch hour.

   :param station: Station identifier. Supported formats include ICAO, WMO, and IGRAv2 identifiers. See :ref:`'site ids' <siteids>`for details.
   :type station: str, required
   :param year: launch year
   :type year: str, required
   :param month: launch month
   :type month: str, required
   :param day: launch day
   :type day: str, required
   :param hour: launch hour
   :type hour: str, required
   :param hush: whether to 'hush' a read-out of thermodynamic and kinematic parameters when getting data.
   :type hush: bool, optional, default is `False`
   :param clean_it: whether to return the raw_data object or a clean_data dict.
   :type clean_it: bool, optional, default is `True`
   :return: :ref:`clean_data<datadescription>`, a dict of ready-to-use vertical profile data including pressure, height, temperature, dewpoint, u-wind, v-wind, & profile information
   :rtype: dict

.. note::
   Archived observations may occasionally be missing, incomplete, or assigned
   to a nearby launch hour. If a known sounding cannot be found, try the
   preceding or following UTC hour.


\ 

.. _siteids:

Station Identifiers
^^^^^^^^^^^^^^^^^^^

``get_obs_data()`` accepts several station-identifier formats:

- ICAO identifier: ``"DTX"`` or ``"KDTX"`` -- **note that using 4-character ICAOs is recommended**.
- WMO identifier: ``"72317"``
- IGRAv2 identifier: ``"GMM00010393"``

SounderPy uses the identifier format to determine the appropriate
observational archive automatically.


.. tip::
   Some stations share the same three-character suffix. For example,
   ``"PABR"`` and ``"KABR"`` both end in ``"ABR"``.

   When ambiguity is possible, use the full four-character station identifier.



\ 

.. _raobsites:

Available RAOB Sites:
^^^^^^^^^^^^^^^^^^^^^^^^
.. raw:: html

    <embed>
        <iframe src="https://kylejgillett.github.io/sounderpy/raob_map" width="100%" height="500"></iframe>
    </embed>
    



Example
^^^^^^^

Retrieve the 18 UTC Omaha, Nebraska sounding from 16 June 2014:

.. code-block:: python

   import sounderpy as spy

   data = spy.get_obs_data("OAX", "2014", "06", "16", "18", hush=True,)

The returned profile can be passed directly to SounderPy's plotting and
analysis tools:

.. code-block:: python

   spy.build_sounding(data)


***************************************************************



.. _acarsdata:


Observed Data | ACARS Aircraft Profiles
-----------------------------------------

SounderPy provides the ``acars_data`` class for retrieving observed aircraft
vertical profiles from the ACARS archive.

Unlike the retrieval functions above, ACARS access uses a two-step workflow:

#. Create an ``acars_data`` object for a requested date and hour.
#. Use ``.list_profiles()`` to find available profiles, then retrieve one with ``.get_profile()``.

- To learn more about ACARS, check out the 'AIRCRAFT' section of this webpage: `NOAA Observing Systems <https://www.weather.gov/about/observation-equipment>`_

.. class:: acars_data(year, month, day, hour)

   :param year: observation year
   :type year: str, required
   :param month: observation month
   :type month: str, required
   :param day: observation day
   :type day: str, required
   :param hour: observation hour
   :type hour: str, required


   .. method:: list_profiles()

      Return a list of profile identifiers available for the selected
      date and hour.

      :return: Available ACARS profile identifiers.
      :rtype: list


   .. method:: get_profile(profile, hush=False, clean_it=True)

      Retrieve and process a selected ACARS profile.

      :param profile: Profile identifier returned by ``list_profiles()``.
      :type profile: str, required

      :param hush: Suppress the default parameter-summary output.
      :type hush: bool, optional, default is ``False``

      :param clean_it: Return SounderPy ``clean_data`` when ``True``; return
                       the raw source data when ``False``.
      :type clean_it: bool, optional, default is ``True``

      :return: A SounderPy ``clean_data`` dictionary containing the aircraft
               profile and associated metadata. See
               :ref:`clean_data schema <clean_data_scheme>`.
      :rtype: dict





.. _acarslists:

ACARS Data Retrieval Example
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Create an ACARS archive connection:


.. code-block:: python

   import sounderpy as spy

   acars = spy.acars_data("2024", "05", "21", "18")


List the profiles available during that hour:


.. code-block:: python

   profiles = acars.list_profiles()

   print(profiles)


Retrieve one of the returned profile IDs:

.. code-block:: python

   data = acars.get_profile(profiles[0], hush=True)


The returned profile can then be used with SounderPy's plotting and analysis
tools:

.. code-block:: python

   spy.build_sounding(data)


Profile identifiers typically contain an airport/site identifier and observation
time, for example:

.. code-block:: text

   ATL_1450
   AUS_1430
   BNA_1420

.. note::
   **Scientific consideration:** ACARS profiles are aircraft observations and
   are often shallower than radiosonde profiles. Many profiles extend only
   through the lower or middle troposphere rather than reaching traditional
   sounding-top levels near 100 hPa.

   Aircraft profiles may also contain occasional quality-control issues,
   including unrealistic dewpoints or wind values. Users should inspect the
   profile before interpreting derived thermodynamic or kinematic parameters.

.. tip::
   The contents of ``list_profiles()`` vary by date and hour. If no profiles
   are returned, try a nearby hour or a different date.



****************************************



.. _datadescription:

What Does SounderPy Return?
------------------------------

SounderPy retrieval functions return a ``clean_data`` dictionary containing
the atmospheric profile and associated metadata.

The core fields are:

.. code-block:: text

   p          pressure
   z          height
   T          temperature
   Td         dewpoint
   u          u-component wind
   v          v-component wind
   site_info  profile metadata

Some sources may also provide fields such as ``rh``, ``omega``, and ``titles``.

For the complete structure, units, metadata fields, and custom-data
requirements, see :ref:`clean_data schema <clean_data_scheme>`.
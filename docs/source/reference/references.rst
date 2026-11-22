.. _references:

==============================
📖 Parameters & References
==============================

The parameters SounderPy computes, where each one comes from, and the full list of references behind
SounderPy's calculations, data sources, and plots. The same references are listed in the header at the
top of each SounderPy module.

***************************************************************

Citing SounderPy
================

If SounderPy contributes to your work, please cite its Journal of Open Source Software article:

    Gillett, K. J., 2025: SounderPy: An atmospheric sounding visualization and analysis tool for Python. J. Open Source Software, 10 (112), 8087, https://doi.org/10.21105/joss.08087

and the references below for the parameters and data you use.

***************************************************************

Parameters
==========

``spy.sounding_params(clean_data).calc()`` returns these values in four dictionaries: ``general``,
``thermo``, ``kinem``, and ``intrp`` (the profile interpolated to 100 m levels). See
the :doc:`helper tools </reference/helpertools>` page.

``spy.winter_params(clean_data).calc()`` adds a ``winter`` dictionary; its parameters and their sources
are described in the :doc:`winter soundings tutorial </tutorials/winter_soundings>`, and listed under
:ref:`ref-winter` below.

General values (``general``)
----------------------------

Surface values and values at every level of the profile.

.. list-table::
   :header-rows: 1
   :widths: 30 44 12 14

   * - Key(s)
     - Parameter
     - Units
     - Source
   * - ``elevation, sfc_pressure``
     - Surface elevation and pressure
     - m, hPa
     - SounderPy
   * - ``wet_bulb, rel_humidity, spec_humidity, mix_ratio``
     - Wet-bulb temperature, relative humidity, specific humidity, and mixing ratio at each level
     - °C, %, g/kg, g/kg
     - :ref:`May et al. 2023 (MetPy) <ref-metpy>`
   * - ``virt_temp, theta, theta_e``
     - Virtual, potential, and equivalent potential temperature at each level
     - °C, K, K
     - :ref:`May et al. 2023 (MetPy) <ref-metpy>`, :ref:`Bolton 1980 <ref-bolton1980>`
   * - ``pwat``
     - Precipitable water
     - in
     - :ref:`May et al. 2023 (MetPy) <ref-metpy>`
   * - ``frz_pt_p, frz_pt_z, wb_frz_pt_p, wb_frz_pt_z``
     - Freezing level and wet-bulb zero height
     - hPa, m AGL
     - SounderPy
   * - ``rh_0_500 ... rh_6_9000, w_0_500 ... w_6_9000``
     - Layer-mean relative humidity and mixing ratio (0-0.5, 0-1, 1-3, 3-6, and 6-9 km)
     - %, g/kg
     - SounderPy

Thermodynamic values (``thermo``)
---------------------------------

``sb_``, ``ml_``, and ``mu_`` are the surface-based, mixed-layer (lowest 100 hPa), and most-unstable parcels.

.. list-table::
   :header-rows: 1
   :widths: 30 44 12 14

   * - Key(s)
     - Parameter
     - Units
     - Source
   * - ``sbcape, mlcape, mucape, sbcin, mlcin, mucin``
     - CAPE and CIN, with the virtual temperature correction
     - J/kg
     - :ref:`Blumberg et al. 2017 (SHARPpy) <ref-sharppy>`, :ref:`Doswell and Rasmussen 1994 <ref-doswell_rasmussen1994>`
   * - ``sb3cape, sb6cape (and ml, mu)``
     - CAPE in the lowest 3 and 6 km
     - J/kg
     - :ref:`Blumberg et al. 2017 (SHARPpy) <ref-sharppy>`
   * - ``sb_lcl_p/_z, sb_lfc_p/_z, sb_el_p/_z, sb_el_T, sb_mpl_p (and ml, mu)``
     - Lifted condensation level, level of free convection, equilibrium level, and maximum parcel level
     - hPa, m AGL, °C
     - :ref:`Blumberg et al. 2017 (SHARPpy) <ref-sharppy>`
   * - ``sbT_trace, sbP_trace (and ml, mu); sbZ_trace, muZ_trace``
     - Parcel temperature trace
     - °C, hPa, m
     - :ref:`Blumberg et al. 2017 (SHARPpy) <ref-sharppy>`
   * - ``sb_ncape, ml_ncape, mu_ncape``
     - Normalized CAPE (CAPE divided by the LFC-to-EL depth)
     - m/s²
     - :ref:`Blanchard 1998 <ref-blanchard1998>`
   * - ``sb_ecape, ml_ecape, mu_ecape``
     - Entraining CAPE
     - J/kg
     - :ref:`Peters et al. 2023 <ref-peters2023>`
   * - ``dcape, dcin, dparcel_p, dparcel_T``
     - Downdraft CAPE and downdraft CIN, and the downdraft parcel trace
     - J/kg, hPa, °C
     - :ref:`Emanuel 1994 <ref-emanuel1994>`, :ref:`Gilmore and Wicker 1998 <ref-gilmore_wicker1998>`, :ref:`Market et al. 2017 <ref-market_dcape>`
   * - ``pbl_top, pbl_depth, pbl_lapse_rate``
     - Planetary boundary layer top, depth, and lapse rate (θv parcel + 0.5 K)
     - hPa, m, °C/km
     - :ref:`Blumberg et al. 2017 (SHARPpy) <ref-sharppy>`
   * - ``dgz, hgz``
     - Dendritic growth zone (-12 to -17 °C) and hail growth zone (-10 to -30 °C), as (bottom, top)
     - hPa
     - :ref:`Blumberg et al. 2017 (SHARPpy) <ref-sharppy>`
   * - ``lr_03km, lr_36km, lr_75``
     - 0-3 km, 3-6 km, and 700-500 hPa lapse rates
     - °C/km
     - :ref:`Blumberg et al. 2017 (SHARPpy) <ref-sharppy>`
   * - ``lr_max``
     - Steepest 2 km lapse rate below 6 km, as (lapse rate, bottom, top)
     - °C/km, hPa
     - :ref:`Blumberg et al. 2017 (SHARPpy) <ref-sharppy>`
   * - ``thetae_diff``
     - Theta-e difference: the minimum theta-e from the 0-3 km theta-e maximum up to 7 km, minus that maximum
     - K
     - :ref:`Cohen et al. 2007 <ref-cohen2007>`
   * - ``k_index, t_totals``
     - K-index and Total Totals
     - 
     - :ref:`George 1960 <ref-george1960>`, :ref:`Miller 1972 <ref-miller1972>`
   * - ``sb_brn, ml_brn, mu_brn``
     - Bulk Richardson number
     - 
     - :ref:`Weisman and Klemp 1982 <ref-weisman_klemp1982>`
   * - ``ship``
     - Significant hail parameter
     - 
     - :ref:`SPC <ref-spc_ship>`, :ref:`Johnson and Sugden 2014 <ref-johnson_sugden2014>`
   * - ``temp_adv``
     - Inferred temperature advection from the thermal wind, as (values, layer bounds)
     - °C/hr, hPa
     - :ref:`Blumberg et al. 2017 (SHARPpy) <ref-sharppy>`

Kinematic values (``kinem``)
----------------------------

Layers are in m AGL (``0_to_500`` is 0-500 m); ``eil`` is the effective inflow layer. Wind vectors are (u, v) in kt.

.. list-table::
   :header-rows: 1
   :widths: 30 44 12 14

   * - Key(s)
     - Parameter
     - Units
     - Source
   * - ``eil, eil_z, eil_depth``
     - Effective inflow layer (bottom, top) and its depth
     - hPa, m AGL, m
     - :ref:`Thompson et al. 2007 <ref-thompson2007>`
   * - ``sm_rm, sm_lm, sm_mw, sm_u, sm_v``
     - Bunkers right- and left-mover motions, the 850-250 hPa mean wind, and the storm motion in use (u, v)
     - kt
     - :ref:`Bunkers et al. 2000 <ref-bunkers2000>`, :ref:`Bunkers et al. 2014 <ref-bunkers2014>`
   * - ``dtm``
     - Deviant tornado motion (storm motion averaged with the 0-300 m mean wind)
     - kt
     - :ref:`Nixon and Allen 2021 <ref-nixon_allen2021>`
   * - ``mcs``
     - Corfidi MCS motion, as (upshear u, v, downshear u, v)
     - kt
     - :ref:`Corfidi 2003 <ref-corfidi2003>`
   * - ``mean_wind_pbl``
     - PBL mean wind
     - kt
     - SounderPy
   * - ``shear_0_to_500 ... shear_6_to_9000, shear_0_to_3000, shear_0_to_6000, shear_eil``
     - Bulk wind difference over each layer
     - kt
     - :ref:`Blumberg et al. 2017 (SHARPpy) <ref-sharppy>`
   * - ``shear_0_to_10000``
     - 0-10 km bulk wind difference
     - kt
     - :ref:`Cohen et al. 2007 <ref-cohen2007>`
   * - ``ebwd``
     - Effective bulk wind difference (EIL base to half the MU EL height)
     - kt
     - :ref:`Thompson et al. 2007 <ref-thompson2007>`
   * - ``srh_0_to_500 ... srh_6_to_9000, srh_0_to_3000, srh_0_to_6000, srh_eil``
     - Storm-relative helicity over each layer
     - m²/s²
     - :ref:`Davies-Jones et al. 1990 <ref-davies_jones1990>`
   * - ``srw, srw_0_to_500 ... srw_eil``
     - Storm-relative wind at each level and its layer means
     - kt
     - SounderPy
   * - ``vort, swv, swv_perc, swv_0_to_500 ... swv_eil, swv_perc_0_to_500 ... swv_perc_eil``
     - Horizontal and streamwise vorticity, and streamwiseness, at each level and as layer means
     - 1/s, %
     - :ref:`Davies-Jones 1984 <ref-davies_jones1984>`, :ref:`Brandt 2022 <ref-brandt2022>`
   * - ``eil_stp``
     - Significant tornado parameter (effective layer, with CIN)
     - 
     - :ref:`Thompson et al. 2012 <ref-thompson2012>`
   * - ``eil_scp``
     - Supercell composite parameter (effective layer)
     - 
     - :ref:`Thompson et al. 2004 <ref-thompson2004>`
   * - ``ehi_0_to_1000, ehi_0_to_3000``
     - Energy-helicity index (MUCAPE with 0-1 and 0-3 km SRH)
     - 
     - :ref:`Hart and Korotky 1991 <ref-hart_korotky1991>`, :ref:`Davies 1993 <ref-davies1993>`

***************************************************************

References
==========

.. _ref-thermodynamics:

Thermodynamics
--------------

.. _ref-doswell_rasmussen1994:

Doswell, C. A., III, and E. N. Rasmussen, 1994: The effect of neglecting the virtual temperature correction on CAPE calculations. Wea. Forecasting, 9, 625-629, `https://doi.org/10.1175/1520-0434(1994)009\<0625:TEONTV\>2.0.CO;2 <https://doi.org/10.1175/1520-0434(1994)009%3C0625:TEONTV%3E2.0.CO;2>`__ *(virtual temperature correction in CAPE)*

.. _ref-blanchard1998:

Blanchard, D. O., 1998: Assessing the vertical distribution of convective available potential energy. Wea. Forecasting, 13, 870-877, `https://doi.org/10.1175/1520-0434(1998)013\<0870:ATVDOC\>2.0.CO;2 <https://doi.org/10.1175/1520-0434(1998)013%3C0870:ATVDOC%3E2.0.CO;2>`__ *(normalized CAPE)*

.. _ref-peters2023:

Peters, J. M., D. R. Chavas, C.-Y. Su, H. Morrison, and B. E. Coffer, 2023: An analytic formula for entraining CAPE in midlatitude storm environments. J. Atmos. Sci., 80, 2165-2186, https://doi.org/10.1175/JAS-D-23-0003.1; computed with the ecape-parcel package, https://github.com/a-urq/ecape-parcel-py *(entraining CAPE)*

.. _ref-emanuel1994:

Emanuel, K. A., 1994: Atmospheric Convection. Oxford University Press, 592 pp. *(DCAPE)*

.. _ref-gilmore_wicker1998:

Gilmore, M. S., and L. J. Wicker, 1998: The influence of midtropospheric dryness on supercell morphology and evolution. Mon. Wea. Rev., 126, 943-958, `https://doi.org/10.1175/1520-0493(1998)126\<0943:TIOMDO\>2.0.CO;2 <https://doi.org/10.1175/1520-0493(1998)126%3C0943:TIOMDO%3E2.0.CO;2>`__ *(DCAPE)*

.. _ref-market_dcape:

Market, P. S., S. M. Rochette, J. Shewchuk, R. Difani, J. S. Kastman, C. B. Henson, and N. I. Fox, 2017: Evaluating elevated convection with the downdraft convective inhibition. Atmos. Sci. Lett., 18, 76-81, https://doi.org/10.1002/asl.727 *(DCIN)*

.. _ref-bolton1980:

Bolton, D., 1980: The computation of equivalent potential temperature. Mon. Wea. Rev., 108, 1046-1053, `https://doi.org/10.1175/1520-0493(1980)108\<1046:TCOEPT\>2.0.CO;2 <https://doi.org/10.1175/1520-0493(1980)108%3C1046:TCOEPT%3E2.0.CO;2>`__ *(dewpoint and equivalent potential temperature)*

.. _ref-weisman_klemp1982:

Weisman, M. L., and J. B. Klemp, 1982: The dependence of numerically simulated convective storms on vertical wind shear and buoyancy. Mon. Wea. Rev., 110, 504-520, `https://doi.org/10.1175/1520-0493(1982)110\<0504:TDONSC\>2.0.CO;2 <https://doi.org/10.1175/1520-0493(1982)110%3C0504:TDONSC%3E2.0.CO;2>`__ *(bulk Richardson number)*

.. _ref-george1960:

George, J. J., 1960: Weather Forecasting for Aeronautics. Academic Press, 673 pp. *(K-index)*

.. _ref-miller1972:

Miller, R. C., 1972: Notes on analysis and severe-storm forecasting procedures of the Air Force Global Weather Central. AWS Tech. Rep. 200 (Rev.), Air Weather Service, Scott AFB, IL, 190 pp. *(Total Totals)*

.. _ref-holzworth1964:

Holzworth, G. C., 1964: Estimates of mean maximum mixing depths in the contiguous United States. Mon. Wea. Rev., 92, 235-242 *(surface modification mixing depth)*

.. _ref-betts1973:

Betts, A. K., 1973: Non-precipitating cumulus convection and its parameterization. Quart. J. Roy. Meteor. Soc., 99, 178-196 *(mixed-layer surface modification)*

.. _ref-carson1973:

Carson, D. J., 1973: The development of a dry inversion-capped convectively unstable boundary layer. Quart. J. Roy. Meteor. Soc., 99, 450-467 *(mixed-layer surface modification)*

.. _ref-tennekes1973:

Tennekes, H., 1973: A model for the dynamics of the inversion above a convective boundary layer. J. Atmos. Sci., 30, 558-567 *(mixed-layer surface modification)*

.. _ref-kinematics-and-storm-motion:

Kinematics and storm motion
---------------------------

.. _ref-bunkers2000:

Bunkers, M. J., B. A. Klimowski, J. W. Zeitler, R. L. Thompson, and M. L. Weisman, 2000: Predicting supercell motion using a new hodograph technique. Wea. Forecasting, 15, 61-79, `https://doi.org/10.1175/1520-0434(2000)015\<0061:PSMUAN\>2.0.CO;2 <https://doi.org/10.1175/1520-0434(2000)015%3C0061:PSMUAN%3E2.0.CO;2>`__ *(Bunkers storm motion)*

.. _ref-bunkers2014:

Bunkers, M. J., D. A. Barber, R. L. Thompson, R. Edwards, and J. Garner, 2014: Choosing a universal mean wind for supercell motion prediction. J. Operational Meteor., 2 (11), 115-129, https://doi.org/10.15191/nwajom.2014.0211 *(Bunkers storm motion from the effective inflow layer)*

.. _ref-corfidi2003:

Corfidi, S. F., 2003: Cold pools and MCS propagation: Forecasting the motion of downwind-developing MCSs. Wea. Forecasting, 18, 997-1017, `https://doi.org/10.1175/1520-0434(2003)018\<0997:CPAMPF\>2.0.CO;2 <https://doi.org/10.1175/1520-0434(2003)018%3C0997:CPAMPF%3E2.0.CO;2>`__ *(Corfidi MCS motion)*

.. _ref-nixon_allen2021:

Nixon, C. J., and J. T. Allen, 2021: Anticipating deviant tornado motion using a simple hodograph technique. Wea. Forecasting, 36, 219-235, https://doi.org/10.1175/WAF-D-20-0056.1 *(deviant tornado motion)*

.. _ref-davies_jones1984:

Davies-Jones, R., 1984: Streamwise vorticity: The origin of updraft rotation in supercell storms. J. Atmos. Sci., 41, 2991-3006, `https://doi.org/10.1175/1520-0469(1984)041\<2991:SVTOOU\>2.0.CO;2 <https://doi.org/10.1175/1520-0469(1984)041%3C2991:SVTOOU%3E2.0.CO;2>`__ *(streamwise vorticity)*

.. _ref-davies_jones1990:

Davies-Jones, R., D. W. Burgess, and M. P. Foster, 1990: Test of helicity as a tornado forecast parameter. Preprints, 16th Conf. on Severe Local Storms, Kananaskis Park, AB, Canada, Amer. Meteor. Soc., 588-592 *(SRH)*

.. _ref-brandt2022:

Brandt, S., 2022: Storm-Relative Hodograph Plotter (Streamwise_Hodograph_Operational.py). GitHub, https://github.com/SamBrandtMeteo/Storm-Relative-Hodograph-Plotter *(streamwise vorticity)*

.. _ref-thompson2007:

Thompson, R. L., C. M. Mead, and R. Edwards, 2007: Effective storm-relative helicity and bulk shear in supercell thunderstorm environments. Wea. Forecasting, 22, 102-115, https://doi.org/10.1175/WAF969.1 *(effective inflow layer, effective bulk shear, and effective SRH)*

.. _ref-browning_wexler1968:

Browning, K. A., and R. Wexler, 1968: The determination of kinematic properties of a wind field using Doppler radar. J. Appl. Meteor., 7, 105-113, `https://doi.org/10.1175/1520-0450(1968)007\<0105:TDOKPO\>2.0.CO;2 <https://doi.org/10.1175/1520-0450(1968)007%3C0105:TDOKPO%3E2.0.CO;2>`__ *(VAD wind profiles)*

.. _ref-composite-parameters:

Composite parameters
--------------------

.. _ref-thompson2004:

Thompson, R. L., R. Edwards, and C. M. Mead, 2004: An update to the supercell composite and significant tornado parameters. Preprints, 22nd Conf. on Severe Local Storms, Hyannis, MA, Amer. Meteor. Soc., P8.1, https://www.spc.noaa.gov/publications/thompson/stp_scp.pdf *(SCP)*

.. _ref-thompson2012:

Thompson, R. L., B. T. Smith, J. S. Grams, A. R. Dean, and C. Broyles, 2012: Convective modes for significant severe thunderstorms in the contiguous United States. Part II: Supercell and QLCS tornado environments. Wea. Forecasting, 27, 1136-1154, https://doi.org/10.1175/WAF-D-11-00116.1 *(STP)*

.. _ref-hart_korotky1991:

Hart, J. A., and W. Korotky, 1991: The SHARP Workstation v1.50 Users Guide. NOAA/National Weather Service, 30 pp. *(EHI)*

.. _ref-davies1993:

Davies, J. M., 1993: Hourly helicity, instability, and EHI in forecasting supercell tornadoes. Preprints, 17th Conf. on Severe Local Storms, St. Louis, MO, Amer. Meteor. Soc., 107-111 *(EHI)*

.. _ref-spc_ship:

Storm Prediction Center: Significant Hail Parameter (SHIP). https://www.spc.noaa.gov/exper/mesoanalysis/help/help_sigh.html *(SHIP)*

.. _ref-johnson_sugden2014:

Johnson, A. W., and K. E. Sugden, 2014: Evaluation of sounding-derived thermodynamic and wind-related parameters associated with large hail events. Electronic J. Severe Storms Meteor., 9 (5), 1-42, https://doi.org/10.55599/ejssm.v9i5.57 *(SHIP)*

.. _ref-cohen2007:

Cohen, A. E., M. C. Coniglio, S. F. Corfidi, and S. J. Corfidi, 2007: Discrimination of mesoscale convective system environments using sounding observations. Wea. Forecasting, 22, 1045-1062, https://doi.org/10.1175/WAF1040.1 *(0-10 km bulk shear and theta-e difference)*

.. _ref-winter:

Winter parameters
-----------------

.. _ref-alduchov1996:

Alduchov, O. A., and R. E. Eskridge, 1996: Improved Magnus form approximation of saturation vapor pressure. J. Appl. Meteor., 35, 601-609 *(RH over ice)*

.. _ref-banacos2014:

Banacos, P. C., A. N. Loconto, and G. A. DeVoir, 2014: Snow squalls: Forecasting and hazard mitigation. J. Operational Meteor., 2 (12), 130-151, https://doi.org/10.15191/nwajom.2014.0212 *(snow squall parameter)*

.. _ref-bernstein2005:

Bernstein, B. C., and Coauthors, 2005: Current icing potential: Algorithm description and comparison with aircraft observations. J. Appl. Meteor., 44, 969-986, https://doi.org/10.1175/JAM2246.1 *(icing layers)*

.. _ref-birk2021:

Birk, K., E. Lenning, K. Donofrio, and M. T. Friedlein, 2021: A revised Bourgouin precipitation-type algorithm. Wea. Forecasting, 36, 425-438, https://doi.org/10.1175/WAF-D-20-0118.1 *(wet-bulb energies; thresholds for warm layers aloft)*

.. _ref-bourgouin2000:

Bourgouin, P., 2000: A method to determine precipitation types. Wea. Forecasting, 15, 583-592 *(precipitation type: warm surface layer)*

.. _ref-gultepe2007:

Gultepe, I., and Coauthors, 2007: Fog research: A review of past achievements and future perspectives. Pure Appl. Geophys., 164, 1121-1159, https://doi.org/10.1007/s00024-007-0211-x *(freezing fog)*

.. _ref-huffman1988:

Huffman, G. J., and G. A. Norman Jr., 1988: The supercooled warm rain process and the specification of freezing precipitation. Mon. Wea. Rev., 116, 2172-2182 *(precipitation source phase, freezing drizzle)*

.. _ref-libbrecht2005:

Libbrecht, K. G., 2005: The physics of snow crystals. Rep. Prog. Phys., 68, 855-895, https://doi.org/10.1088/0034-4885/68/4/R03 *(DGZ)*

.. _ref-mcdonough2004:

McDonough, F., B. C. Bernstein, M. K. Politovich, and C. A. Wolff, 2004: The Forecast Icing Potential (FIP) algorithm. Preprints, 11th Conf. on Aviation, Range, and Aerospace Meteorology, Hyannis, MA, Amer. Meteor. Soc., 3.2 *(icing layers, freezing drizzle)*

.. _ref-kuchera:

NOAA VLab: NWP snow accumulation products. Forecast Guide *(snow-to-liquid ratio)*

.. _ref-nws_windchill:

NWS, 2001: Wind chill chart. https://www.weather.gov/safety/cold-wind-chill-chart *(wind chill)*

.. _ref-osczevski2005:

Osczevski, R., and M. Bluestein, 2005: The new wind chill equivalent temperature chart. Bull. Amer. Meteor. Soc., 86, 1453-1458, https://doi.org/10.1175/BAMS-86-10-1453 *(wind chill)*

.. _ref-rauber2000:

Rauber, R. M., L. S. Olthoff, M. K. Ramamurthy, and K. E. Kunkel, 2000: The relative importance of warm rain and melting processes in freezing precipitation events. J. Appl. Meteor., 39, 1185-1195 *(precipitation source phase, freezing drizzle)*

.. _ref-takahashi1991:

Takahashi, T., T. Endoh, G. Wakahama, and N. Fukuta, 1991: Vapor diffusional growth of free-falling snow crystals between -3 and -23 C. J. Meteor. Soc. Japan, 69, 15-30, https://doi.org/10.2151/jmsj1965.69.1_15 *(DGZ)*

.. _ref-sczepanski2026:

Sczepanski, A., A. D. Kennedy, N. Wood, D. H. Bromwich, and S. D. Shuvo, 2024: Synthesis of in situ and remotely sensed observations in North Dakota to understand Northern Great Plains ground blizzards. 24th Symp. on Meteorological Observation and Instrumentation, Baltimore, MD, Amer. Meteor. Soc., 645, https://ams.confex.com/ams/104ANNUAL/webprogram/Paper433440.html *(PBL depth, lapse rate, and mixing for blowing snow)*

.. _ref-smoothing:

Smoothing
---------

.. _ref-savitzky_golay1964:

Savitzky, A., and M. J. E. Golay, 1964: Smoothing and differentiation of data by simplified least squares procedures. Anal. Chem., 36, 1627-1639, https://doi.org/10.1021/ac60214a047 *(Savitzky-Golay smoothing)*

.. _ref-models-and-data-sources:

Models and data sources
-----------------------

.. _ref-benjamin2016:

Benjamin, S. G., and Coauthors, 2016: A North American hourly assimilation and model forecast cycle: The Rapid Refresh. Mon. Wea. Rev., 144, 1669-1694, https://doi.org/10.1175/MWR-D-15-0242.1 *(RAP)*

.. _ref-dowell2022:

Dowell, D. C., and Coauthors, 2022: The High-Resolution Rapid Refresh (HRRR): An hourly updating convection-allowing forecast model. Part I: Motivation and system description. Wea. Forecasting, 37, 1371-1395, https://doi.org/10.1175/WAF-D-21-0151.1 *(HRRR)*

.. _ref-hersbach2020:

Hersbach, H., and Coauthors, 2020: The ERA5 global reanalysis. Quart. J. Roy. Meteor. Soc., 146, 1999-2049, https://doi.org/10.1002/qj.3803 *(ERA5)*

.. _ref-ncep_fnl:

National Centers for Environmental Prediction/National Weather Service/NOAA/U.S. Department of Commerce, 2000: NCEP FNL Operational Model Global Tropospheric Analyses, continuing from July 1999. NSF National Center for Atmospheric Research, https://doi.org/10.5065/D6M043C6 *(NCEP-FNL)*

.. _ref-durre2018:

Durre, I., X. Yin, R. S. Vose, S. Applequist, and J. Arnfield, 2018: Enhancing the data coverage in the Integrated Global Radiosonde Archive. J. Atmos. Oceanic Technol., 35, 1753-1770, https://doi.org/10.1175/JTECH-D-17-0223.1 *(IGRA v2)*

.. _ref-moninger2003:

Moninger, W. R., R. D. Mamrosh, and P. M. Pauley, 2003: Automated meteorological reports from commercial aircraft. Bull. Amer. Meteor. Soc., 84, 203-216, https://doi.org/10.1175/BAMS-84-2-203 *(ACARS aircraft observations)*

.. _ref-bryan_fritsch2002:

Bryan, G. H., and J. M. Fritsch, 2002: A benchmark simulation for moist nonhydrostatic numerical models. Mon. Wea. Rev., 130, 2917-2928, `https://doi.org/10.1175/1520-0493(2002)130\<2917:ABSFMN\>2.0.CO;2 <https://doi.org/10.1175/1520-0493(2002)130%3C2917:ABSFMN%3E2.0.CO;2>`__ *(CM1)*

.. _ref-skamarock2019:

Skamarock, W. C., and Coauthors, 2019: A Description of the Advanced Research WRF Model Version 4. NCAR Tech. Note NCAR/TN-556+STR, https://doi.org/10.5065/1dfh-6p97 *(WRF-ARW)*

.. _ref-wmo_bufr:

WMO: Manual on Codes, Volume I.2: International Codes, Part B (Binary Codes), FM 94 BUFR. WMO-No. 306, World Meteorological Organization, Geneva, https://library.wmo.int/idurl/4/35625 *(BUFR format)*

.. _ref-lofs_sndmod:

Halbert, K., and L. Orf: LOFS-read sndmod (CM1 input_sounding writer). GitHub, https://github.com/leighorf/LOFS-read/blob/master/bin/sndmod *(CM1 input_sounding writer)*

.. _ref-psu_bufkit:

Penn State BUFKIT data: https://www.meteo.psu.edu/bufkit/ *(latest BUFKIT runs)*

.. _ref-iem_bufkit:

Iowa Environmental Mesonet BUFKIT archive: https://mtarchive.geol.iastate.edu/ *(archived BUFKIT runs)*

.. _ref-nomads:

NCEP NOMADS: https://nomads.ncep.noaa.gov/ *(RRFS, HRRR, and RAP BUFR soundings)*

.. _ref-uwyo:

University of Wyoming upper air soundings: https://weather.uwyo.edu/upperair/sounding.shtml *(observed soundings)*

.. _ref-iem_asos:

Iowa Environmental Mesonet ASOS archive: https://mesonet.agron.iastate.edu/request/download.phtml *(METARs)*

.. _ref-ou_acars:

University of Oklahoma ACARS sounding archive: https://sharp.weather.ou.edu/soundings/acars/ *(ACARS profiles)*

.. _ref-nexrad_aws:

NOAA NEXRAD on AWS (Unidata Level II archive): https://registry.opendata.aws/noaa-nexrad/ *(single-site radar)*

.. _ref-noaa_rap_aws:

NOAA Rapid Refresh (RAP) on AWS: https://registry.opendata.aws/noaa-rap/ *(RAP GRIB files, 2021 on)*

.. _ref-ncei_rap:

NCEI Rapid Refresh / Rapid Update Cycle archive: https://www.ncei.noaa.gov/products/weather-climate-models/rapid-refresh-update *(RAP/RUC GRIB archive)*

.. _ref-python-packages:

Python packages
---------------

.. _ref-sharppy:

Blumberg, W. G., K. T. Halbert, T. A. Supinie, P. T. Marsh, R. L. Thompson, and J. A. Hart, 2017: SHARPpy: An open-source sounding analysis toolkit for the atmospheric sciences. Bull. Amer. Meteor. Soc., 98, 1625-1636, https://doi.org/10.1175/BAMS-D-15-00309.1 *(parcels, indices, and storm motions)*

.. _ref-metpy:

May, R. M., S. C. Arms, P. Marsh, E. Bruning, J. R. Leeman, K. Goebbert, J. E. Thielen, Z. S. Bruick, and M. D. Camron, 2023: MetPy: A Python package for meteorological data. Unidata, https://doi.org/10.5065/D6WW7G29 *(thermodynamic calculations)*

.. _ref-siphon:

May, R. M., S. C. Arms, J. R. Leeman, and J. Chastang, 2017: Siphon: A collection of Python utilities for accessing remote atmospheric and oceanic datasets. Unidata, https://doi.org/10.5065/D6CN72NW *(remote data access)*

.. _ref-xarray:

Hoyer, S., and J. Hamman, 2017: xarray: N-D labeled arrays and datasets in Python. J. Open Res. Software, 5 (1), 10, https://doi.org/10.5334/jors.148 *(gridded data)*

.. _ref-matplotlib:

Hunter, J. D., 2007: Matplotlib: A 2D graphics environment. Computing in Science & Engineering, 9 (3), 90-95 *(figures)*

.. _ref-scipy:

Virtanen, P., R. Gommers, T. E. Oliphant, and Coauthors, 2020: SciPy 1.0: Fundamental algorithms for scientific computing in Python. Nature Methods, 17 (3), 261-272 *(filters and interpolation)*

.. _ref-cartopy:

Met Office, 2010-2015: Cartopy: a cartographic Python library with a Matplotlib interface. Exeter, Devon, https://cartopy.readthedocs.io *(map insets)*

.. _ref-pyart:

Helmus, J. J., and S. M. Collis, 2016: The Python ARM Radar Toolkit (Py-ART), a library for working with weather radar data in the Python programming language. J. Open Res. Software, 4 (1), e25, https://doi.org/10.5334/jors.119 *(NEXRAD Level II radar data)*

.. _ref-cfgrib:

ECMWF: cfgrib, a Python interface to map GRIB files to xarray. https://github.com/ecmwf/cfgrib *(GRIB decoding)*

.. _ref-cdsapi:

ECMWF: cdsapi, the Copernicus Climate Data Store API client. https://github.com/ecmwf/cdsapi *(ERA5 downloads)*

.. _ref-fsspec:

fsspec / s3fs: filesystem interfaces for Python. https://filesystem-spec.readthedocs.io/ *(AWS downloads)*

.. _ref-bs4:

Richardson, L.: Beautiful Soup. https://www.crummy.com/software/BeautifulSoup/ *(ACARS archive listings)*

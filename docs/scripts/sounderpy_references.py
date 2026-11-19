#!/usr/bin/env python3
"""
SounderPy's reference list, shared by the module headers and the docs.

    REFS        every verified reference, keyed by a short name
    TAGS        what each reference is cited for, added after it as "(...)"
    SHORT       short names for the docs' "Source" column
    GROUPS      how the docs page groups the reference list
    PARAMETERS  the parameters sounding_params(...).calc() returns, with their sources

Used by module_headers.py (the header at the top of each module in src/sounderpy) and
make_references_page.py (docs/source/reference/references.rst). Journal references were checked
against Crossref or DataCite; add new ones the same way, then rerun both scripts.
"""

REFS = {'sharppy': 'Blumberg, W. G., K. T. Halbert, T. A. Supinie, P. T. Marsh, R. L. Thompson, and J. A. Hart, 2017: '
            'SHARPpy: An open-source sounding analysis toolkit for the atmospheric sciences. Bull. Amer. Meteor. '
            'Soc., 98, 1625-1636, https://doi.org/10.1175/BAMS-D-15-00309.1',
 'metpy': 'May, R. M., S. C. Arms, P. Marsh, E. Bruning, J. R. Leeman, K. Goebbert, J. E. Thielen, Z. S. Bruick, and '
          'M. D. Camron, 2023: MetPy: A Python package for meteorological data. Unidata, '
          'https://doi.org/10.5065/D6WW7G29',
 'siphon': 'May, R. M., S. C. Arms, J. R. Leeman, and J. Chastang, 2017: Siphon: A collection of Python utilities '
           'for accessing remote atmospheric and oceanic datasets. Unidata, https://doi.org/10.5065/D6CN72NW',
 'xarray': 'Hoyer, S., and J. Hamman, 2017: xarray: N-D labeled arrays and datasets in Python. J. Open Res. '
           'Software, 5 (1), 10, https://doi.org/10.5334/jors.148',
 'matplotlib': 'Hunter, J. D., 2007: Matplotlib: A 2D graphics environment. Computing in Science & Engineering, 9 '
               '(3), 90-95',
 'scipy': 'Virtanen, P., R. Gommers, T. E. Oliphant, and Coauthors, 2020: SciPy 1.0: Fundamental algorithms for '
          'scientific computing in Python. Nature Methods, 17 (3), 261-272',
 'cartopy': 'Met Office, 2010-2015: Cartopy: a cartographic Python library with a Matplotlib interface. Exeter, '
            'Devon, https://cartopy.readthedocs.io',
 'pyart': 'Helmus, J. J., and S. M. Collis, 2016: The Python ARM Radar Toolkit (Py-ART), a library for working with '
          'weather radar data in the Python programming language. J. Open Res. Software, 4 (1), e25, '
          'https://doi.org/10.5334/jors.119',
 'cfgrib': 'ECMWF: cfgrib, a Python interface to map GRIB files to xarray. https://github.com/ecmwf/cfgrib',
 'cdsapi': 'ECMWF: cdsapi, the Copernicus Climate Data Store API client. https://github.com/ecmwf/cdsapi',
 'fsspec': 'fsspec / s3fs: filesystem interfaces for Python. https://filesystem-spec.readthedocs.io/',
 'bs4': 'Richardson, L.: Beautiful Soup. https://www.crummy.com/software/BeautifulSoup/',
 'sounderpy_joss': 'Gillett, K. J., 2025: SounderPy: An atmospheric sounding visualization and analysis tool for '
                   'Python. J. Open Source Software, 10 (112), 8087, https://doi.org/10.21105/joss.08087',
 'peters2023': 'Peters, J. M., D. R. Chavas, C.-Y. Su, H. Morrison, and B. E. Coffer, 2023: An analytic formula for '
               'entraining CAPE in midlatitude storm environments. J. Atmos. Sci., 80, 2165-2186, '
               'https://doi.org/10.1175/JAS-D-23-0003.1; computed with the ecape-parcel package, '
               'https://github.com/a-urq/ecape-parcel-py',
 'emanuel1994': 'Emanuel, K. A., 1994: Atmospheric Convection. Oxford University Press, 592 pp.',
 'market_dcape': 'Market, P. S., S. M. Rochette, J. Shewchuk, R. Difani, J. S. Kastman, C. B. Henson, and N. I. Fox, '
                 '2017: Evaluating elevated convection with the downdraft convective inhibition. Atmos. Sci. Lett., '
                 '18, 76-81, https://doi.org/10.1002/asl.727',
 'gilmore_wicker1998': 'Gilmore, M. S., and L. J. Wicker, 1998: The influence of midtropospheric dryness on '
                       'supercell morphology and evolution. Mon. Wea. Rev., 126, 943-958, '
                       'https://doi.org/10.1175/1520-0493(1998)126<0943:TIOMDO>2.0.CO;2',
 'bunkers2000': 'Bunkers, M. J., B. A. Klimowski, J. W. Zeitler, R. L. Thompson, and M. L. Weisman, 2000: Predicting '
                'supercell motion using a new hodograph technique. Wea. Forecasting, 15, 61-79, '
                'https://doi.org/10.1175/1520-0434(2000)015<0061:PSMUAN>2.0.CO;2',
 'corfidi2003': 'Corfidi, S. F., 2003: Cold pools and MCS propagation: Forecasting the motion of downwind-developing '
                'MCSs. Wea. Forecasting, 18, 997-1017, '
                'https://doi.org/10.1175/1520-0434(2003)018<0997:CPAMPF>2.0.CO;2',
 'davies_jones1984': 'Davies-Jones, R., 1984: Streamwise vorticity: The origin of updraft rotation in supercell '
                     'storms. J. Atmos. Sci., 41, 2991-3006, '
                     'https://doi.org/10.1175/1520-0469(1984)041<2991:SVTOOU>2.0.CO;2',
 'davies_jones1990': 'Davies-Jones, R., D. W. Burgess, and M. P. Foster, 1990: Test of helicity as a tornado '
                     'forecast parameter. Preprints, 16th Conf. on Severe Local Storms, Kananaskis Park, AB, Canada, '
                     'Amer. Meteor. Soc., 588-592',
 'thompson2004': 'Thompson, R. L., R. Edwards, and C. M. Mead, 2004: An update to the supercell composite and '
                 'significant tornado parameters. Preprints, 22nd Conf. on Severe Local Storms, Hyannis, MA, Amer. '
                 'Meteor. Soc., P8.1, https://www.spc.noaa.gov/publications/thompson/stp_scp.pdf',
 'thompson2007': 'Thompson, R. L., C. M. Mead, and R. Edwards, 2007: Effective storm-relative helicity and bulk '
                 'shear in supercell thunderstorm environments. Wea. Forecasting, 22, 102-115, '
                 'https://doi.org/10.1175/WAF969.1',
 'thompson2012': 'Thompson, R. L., B. T. Smith, J. S. Grams, A. R. Dean, and C. Broyles, 2012: Convective modes for '
                 'significant severe thunderstorms in the contiguous United States. Part II: Supercell and QLCS '
                 'tornado environments. Wea. Forecasting, 27, 1136-1154, https://doi.org/10.1175/WAF-D-11-00116.1',
 'weisman_klemp1982': 'Weisman, M. L., and J. B. Klemp, 1982: The dependence of numerically simulated convective '
                      'storms on vertical wind shear and buoyancy. Mon. Wea. Rev., 110, 504-520, '
                      'https://doi.org/10.1175/1520-0493(1982)110<0504:TDONSC>2.0.CO;2',
 'george1960': 'George, J. J., 1960: Weather Forecasting for Aeronautics. Academic Press, 673 pp.',
 'miller1972': 'Miller, R. C., 1972: Notes on analysis and severe-storm forecasting procedures of the Air Force '
               'Global Weather Central. AWS Tech. Rep. 200 (Rev.), Air Weather Service, Scott AFB, IL, 190 pp.',
 'spc_ship': 'Storm Prediction Center: Significant Hail Parameter (SHIP). '
             'https://www.spc.noaa.gov/exper/mesoanalysis/help/help_sigh.html',
 'cohen2007': 'Cohen, A. E., M. C. Coniglio, S. F. Corfidi, and S. J. Corfidi, 2007: Discrimination of mesoscale '
              'convective system environments using sounding observations. Wea. Forecasting, 22, 1045-1062, '
              'https://doi.org/10.1175/WAF1040.1',
 'holzworth1964': 'Holzworth, G. C., 1964: Estimates of mean maximum mixing depths in the contiguous United States. '
                  'Mon. Wea. Rev., 92, 235-242',
 'betts1973': 'Betts, A. K., 1973: Non-precipitating cumulus convection and its parameterization. Quart. J. Roy. '
              'Meteor. Soc., 99, 178-196',
 'carson1973': 'Carson, D. J., 1973: The development of a dry inversion-capped convectively unstable boundary layer. '
               'Quart. J. Roy. Meteor. Soc., 99, 450-467',
 'tennekes1973': 'Tennekes, H., 1973: A model for the dynamics of the inversion above a convective boundary layer. '
                 'J. Atmos. Sci., 30, 558-567',
 'savitzky_golay1964': 'Savitzky, A., and M. J. E. Golay, 1964: Smoothing and differentiation of data by simplified '
                       'least squares procedures. Anal. Chem., 36, 1627-1639, https://doi.org/10.1021/ac60214a047',
 'alduchov1996': 'Alduchov, O. A., and R. E. Eskridge, 1996: Improved Magnus form approximation of saturation vapor '
                 'pressure. J. Appl. Meteor., 35, 601-609',
 'banacos2014': 'Banacos, P. C., A. N. Loconto, and G. A. DeVoir, 2014: Snow squalls: Forecasting and hazard '
                'mitigation. J. Operational Meteor., 2 (12), 130-151, https://doi.org/10.15191/nwajom.2014.0212',
 'bernstein2005': 'Bernstein, B. C., and Coauthors, 2005: Current icing potential: Algorithm description and '
                  'comparison with aircraft observations. J. Appl. Meteor., 44, 969-986, '
                  'https://doi.org/10.1175/JAM2246.1',
 'birk2021': 'Birk, K., E. Lenning, K. Donofrio, and M. T. Friedlein, 2021: A revised Bourgouin precipitation-type '
             'algorithm. Wea. Forecasting, 36, 425-438, https://doi.org/10.1175/WAF-D-20-0118.1',
 'bourgouin2000': 'Bourgouin, P., 2000: A method to determine precipitation types. Wea. Forecasting, 15, 583-592',
 'gultepe2007': 'Gultepe, I., and Coauthors, 2007: Fog research: A review of past achievements and future '
                'perspectives. Pure Appl. Geophys., 164, 1121-1159, https://doi.org/10.1007/s00024-007-0211-x',
 'huffman1988': 'Huffman, G. J., and G. A. Norman Jr., 1988: The supercooled warm rain process and the specification '
                'of freezing precipitation. Mon. Wea. Rev., 116, 2172-2182',
 'libbrecht2005': 'Libbrecht, K. G., 2005: The physics of snow crystals. Rep. Prog. Phys., 68, 855-895, '
                  'https://doi.org/10.1088/0034-4885/68/4/R03',
 'mcdonough2004': 'McDonough, F., B. C. Bernstein, M. K. Politovich, and C. A. Wolff, 2004: The Forecast Icing '
                  'Potential (FIP) algorithm. Preprints, 11th Conf. on Aviation, Range, and Aerospace Meteorology, '
                  'Hyannis, MA, Amer. Meteor. Soc., 3.2',
 'kuchera': 'NOAA VLab: NWP snow accumulation products. Forecast Guide',
 'nws_windchill': 'NWS, 2001: Wind chill chart. https://www.weather.gov/safety/cold-wind-chill-chart',
 'osczevski2005': 'Osczevski, R., and M. Bluestein, 2005: The new wind chill equivalent temperature chart. Bull. '
                  'Amer. Meteor. Soc., 86, 1453-1458, https://doi.org/10.1175/BAMS-86-10-1453',
 'rauber2000': 'Rauber, R. M., L. S. Olthoff, M. K. Ramamurthy, and K. E. Kunkel, 2000: The relative importance of '
               'warm rain and melting processes in freezing precipitation events. J. Appl. Meteor., 39, 1185-1195',
 'takahashi1991': 'Takahashi, T., T. Endoh, G. Wakahama, and N. Fukuta, 1991: Vapor diffusional growth of '
                  'free-falling snow crystals between -3 and -23 C. J. Meteor. Soc. Japan, 69, 15-30, '
                  'https://doi.org/10.2151/jmsj1965.69.1_15',
 'benjamin2016': 'Benjamin, S. G., and Coauthors, 2016: A North American hourly assimilation and model forecast '
                 'cycle: The Rapid Refresh. Mon. Wea. Rev., 144, 1669-1694, https://doi.org/10.1175/MWR-D-15-0242.1',
 'dowell2022': 'Dowell, D. C., and Coauthors, 2022: The High-Resolution Rapid Refresh (HRRR): An hourly updating '
               'convection-allowing forecast model. Part I: Motivation and system description. Wea. Forecasting, 37, '
               '1371-1395, https://doi.org/10.1175/WAF-D-21-0151.1',
 'hersbach2020': 'Hersbach, H., and Coauthors, 2020: The ERA5 global reanalysis. Quart. J. Roy. Meteor. Soc., 146, '
                 '1999-2049, https://doi.org/10.1002/qj.3803',
 'ncep_fnl': 'National Centers for Environmental Prediction/National Weather Service/NOAA/U.S. Department of '
             'Commerce, 2000: NCEP FNL Operational Model Global Tropospheric Analyses, continuing from July 1999. '
             'NSF National Center for Atmospheric Research, https://doi.org/10.5065/D6M043C6',
 'durre2018': 'Durre, I., X. Yin, R. S. Vose, S. Applequist, and J. Arnfield, 2018: Enhancing the data coverage in '
              'the Integrated Global Radiosonde Archive. J. Atmos. Oceanic Technol., 35, 1753-1770, '
              'https://doi.org/10.1175/JTECH-D-17-0223.1',
 'moninger2003': 'Moninger, W. R., R. D. Mamrosh, and P. M. Pauley, 2003: Automated meteorological reports from '
                 'commercial aircraft. Bull. Amer. Meteor. Soc., 84, 203-216, https://doi.org/10.1175/BAMS-84-2-203',
 'bryan_fritsch2002': 'Bryan, G. H., and J. M. Fritsch, 2002: A benchmark simulation for moist nonhydrostatic '
                      'numerical models. Mon. Wea. Rev., 130, 2917-2928, '
                      'https://doi.org/10.1175/1520-0493(2002)130<2917:ABSFMN>2.0.CO;2',
 'skamarock2019': 'Skamarock, W. C., and Coauthors, 2019: A Description of the Advanced Research WRF Model Version '
                  '4. NCAR Tech. Note NCAR/TN-556+STR, https://doi.org/10.5065/1dfh-6p97',
 'bolton1980': 'Bolton, D., 1980: The computation of equivalent potential temperature. Mon. Wea. Rev., 108, '
               '1046-1053, https://doi.org/10.1175/1520-0493(1980)108<1046:TCOEPT>2.0.CO;2',
 'wmo_bufr': 'WMO: Manual on Codes, Volume I.2: International Codes, Part B (Binary Codes), FM 94 BUFR. WMO-No. 306, '
             'World Meteorological Organization, Geneva, https://library.wmo.int/idurl/4/35625',
 'lofs_sndmod': 'Halbert, K., and L. Orf: LOFS-read sndmod (CM1 input_sounding writer). GitHub, '
                'https://github.com/leighorf/LOFS-read/blob/master/bin/sndmod',
 'psu_bufkit': 'Penn State BUFKIT data: https://www.meteo.psu.edu/bufkit/',
 'iem_bufkit': 'Iowa Environmental Mesonet BUFKIT archive: https://mtarchive.geol.iastate.edu/',
 'nomads': 'NCEP NOMADS: https://nomads.ncep.noaa.gov/',
 'uwyo': 'University of Wyoming upper air soundings: https://weather.uwyo.edu/upperair/sounding.shtml',
 'iem_asos': 'Iowa Environmental Mesonet ASOS archive: https://mesonet.agron.iastate.edu/request/download.phtml',
 'ou_acars': 'University of Oklahoma ACARS sounding archive: https://sharp.weather.ou.edu/soundings/acars/',
 'nexrad_aws': 'NOAA NEXRAD on AWS (Unidata Level II archive): https://registry.opendata.aws/noaa-nexrad/',
 'noaa_rap_aws': 'NOAA Rapid Refresh (RAP) on AWS: https://registry.opendata.aws/noaa-rap/',
 'ncei_rap': 'NCEI Rapid Refresh / Rapid Update Cycle archive: '
             'https://www.ncei.noaa.gov/products/weather-climate-models/rapid-refresh-update',
 'brandt2022': 'Brandt, S., 2022: Storm-Relative Hodograph Plotter (Streamwise_Hodograph_Operational.py). GitHub, '
               'https://github.com/SamBrandtMeteo/Storm-Relative-Hodograph-Plotter',
 'sczepanski2026': 'Sczepanski, A., A. D. Kennedy, N. Wood, D. H. Bromwich, and S. D. Shuvo, 2024: Synthesis of in '
                   'situ and remotely sensed observations in North Dakota to understand Northern Great Plains ground '
                   'blizzards. 24th Symp. on Meteorological Observation and Instrumentation, Baltimore, MD, Amer. '
                   'Meteor. Soc., 645, https://ams.confex.com/ams/104ANNUAL/webprogram/Paper433440.html',
 'blanchard1998': 'Blanchard, D. O., 1998: Assessing the vertical distribution of convective available potential '
                  'energy. Wea. Forecasting, 13, 870-877, '
                  'https://doi.org/10.1175/1520-0434(1998)013<0870:ATVDOC>2.0.CO;2',
 'hart_korotky1991': 'Hart, J. A., and W. Korotky, 1991: The SHARP Workstation v1.50 Users Guide. NOAA/National '
                     'Weather Service, 30 pp.',
 'davies1993': 'Davies, J. M., 1993: Hourly helicity, instability, and EHI in forecasting supercell tornadoes. '
               'Preprints, 17th Conf. on Severe Local Storms, St. Louis, MO, Amer. Meteor. Soc., 107-111',
 'nixon_allen2021': 'Nixon, C. J., and J. T. Allen, 2021: Anticipating deviant tornado motion using a simple '
                    'hodograph technique. Wea. Forecasting, 36, 219-235, https://doi.org/10.1175/WAF-D-20-0056.1',
 'bunkers2014': 'Bunkers, M. J., D. A. Barber, R. L. Thompson, R. Edwards, and J. Garner, 2014: Choosing a universal '
                'mean wind for supercell motion prediction. J. Operational Meteor., 2 (11), 115-129, '
                'https://doi.org/10.15191/nwajom.2014.0211',
 'johnson_sugden2014': 'Johnson, A. W., and K. E. Sugden, 2014: Evaluation of sounding-derived thermodynamic and '
                       'wind-related parameters associated with large hail events. Electronic J. Severe Storms '
                       'Meteor., 9 (5), 1-42, https://doi.org/10.55599/ejssm.v9i5.57',
 'browning_wexler1968': 'Browning, K. A., and R. Wexler, 1968: The determination of kinematic properties of a wind '
                        'field using Doppler radar. J. Appl. Meteor., 7, 105-113, '
                        'https://doi.org/10.1175/1520-0450(1968)007<0105:TDOKPO>2.0.CO;2',
 'doswell_rasmussen1994': 'Doswell, C. A., III, and E. N. Rasmussen, 1994: The effect of neglecting the virtual '
                          'temperature correction on CAPE calculations. Wea. Forecasting, 9, 625-629, '
                          'https://doi.org/10.1175/1520-0434(1994)009<0625:TEONTV>2.0.CO;2'}

TAGS = {'sharppy': 'parcels, indices, and storm motions',
 'metpy': 'thermodynamic calculations',
 'siphon': 'remote data access',
 'xarray': 'gridded data',
 'matplotlib': 'figures',
 'scipy': 'filters and interpolation',
 'cartopy': 'map insets',
 'pyart': 'NEXRAD Level II radar data',
 'cfgrib': 'GRIB decoding',
 'cdsapi': 'ERA5 downloads',
 'fsspec': 'AWS downloads',
 'bs4': 'ACARS archive listings',
 'sounderpy_joss': 'SounderPy',
 'peters2023': 'entraining CAPE',
 'emanuel1994': 'DCAPE',
 'market_dcape': 'DCIN',
 'gilmore_wicker1998': 'DCAPE',
 'bunkers2000': 'Bunkers storm motion',
 'bunkers2014': 'Bunkers storm motion from the effective inflow layer',
 'corfidi2003': 'Corfidi MCS motion',
 'davies_jones1984': 'streamwise vorticity',
 'davies_jones1990': 'SRH',
 'brandt2022': 'streamwise vorticity',
 'thompson2004': 'SCP',
 'thompson2007': 'effective inflow layer, effective bulk shear, and effective SRH',
 'thompson2012': 'STP',
 'weisman_klemp1982': 'bulk Richardson number',
 'george1960': 'K-index',
 'miller1972': 'Total Totals',
 'spc_ship': 'SHIP',
 'johnson_sugden2014': 'SHIP',
 'cohen2007': '0-10 km bulk shear and theta-e difference',
 'blanchard1998': 'normalized CAPE',
 'hart_korotky1991': 'EHI',
 'davies1993': 'EHI',
 'nixon_allen2021': 'deviant tornado motion',
 'browning_wexler1968': 'VAD wind profiles',
 'doswell_rasmussen1994': 'virtual temperature correction in CAPE',
 'holzworth1964': 'surface modification mixing depth',
 'betts1973': 'mixed-layer surface modification',
 'carson1973': 'mixed-layer surface modification',
 'tennekes1973': 'mixed-layer surface modification',
 'savitzky_golay1964': 'Savitzky-Golay smoothing',
 'alduchov1996': 'RH over ice',
 'banacos2014': 'snow squall parameter',
 'bernstein2005': 'icing layers',
 'birk2021': 'wet-bulb energies; thresholds for warm layers aloft',
 'bourgouin2000': 'precipitation type: warm surface layer',
 'gultepe2007': 'freezing fog',
 'huffman1988': 'precipitation source phase, freezing drizzle',
 'libbrecht2005': 'DGZ',
 'mcdonough2004': 'icing layers, freezing drizzle',
 'kuchera': 'snow-to-liquid ratio',
 'nws_windchill': 'wind chill',
 'osczevski2005': 'wind chill',
 'rauber2000': 'precipitation source phase, freezing drizzle',
 'takahashi1991': 'DGZ',
 'sczepanski2026': 'PBL depth, lapse rate, and mixing for blowing snow',
 'benjamin2016': 'RAP',
 'dowell2022': 'HRRR',
 'hersbach2020': 'ERA5',
 'ncep_fnl': 'NCEP-FNL',
 'durre2018': 'IGRA v2',
 'moninger2003': 'ACARS aircraft observations',
 'bryan_fritsch2002': 'CM1',
 'skamarock2019': 'WRF-ARW',
 'bolton1980': 'dewpoint and equivalent potential temperature',
 'wmo_bufr': 'BUFR format',
 'lofs_sndmod': 'CM1 input_sounding writer',
 'psu_bufkit': 'latest BUFKIT runs',
 'iem_bufkit': 'archived BUFKIT runs',
 'nomads': 'RRFS, HRRR, and RAP BUFR soundings',
 'uwyo': 'observed soundings',
 'iem_asos': 'METARs',
 'ou_acars': 'ACARS profiles',
 'nexrad_aws': 'single-site radar',
 'noaa_rap_aws': 'RAP GRIB files, 2021 on',
 'ncei_rap': 'RAP/RUC GRIB archive'}

# Short names for the docs page's "Source" column (one per reference)
SHORT = {
    "sharppy": "Blumberg et al. 2017 (SHARPpy)", "metpy": "May et al. 2023 (MetPy)", "siphon": "May et al. 2017 (Siphon)",
    "xarray": "Hoyer and Hamman 2017 (xarray)", "matplotlib": "Hunter 2007 (Matplotlib)", "scipy": "Virtanen et al. 2020 (SciPy)",
    "cartopy": "Met Office (Cartopy)", "pyart": "Helmus and Collis 2016 (Py-ART)", "cfgrib": "ECMWF (cfgrib)",
    "cdsapi": "ECMWF (cdsapi)", "fsspec": "fsspec / s3fs", "bs4": "Beautiful Soup", "sounderpy_joss": "Gillett 2025",
    "peters2023": "Peters et al. 2023", "emanuel1994": "Emanuel 1994", "market_dcape": "Market et al. 2017",
    "gilmore_wicker1998": "Gilmore and Wicker 1998", "bunkers2000": "Bunkers et al. 2000", "bunkers2014": "Bunkers et al. 2014",
    "corfidi2003": "Corfidi 2003", "davies_jones1984": "Davies-Jones 1984", "davies_jones1990": "Davies-Jones et al. 1990",
    "brandt2022": "Brandt 2022", "thompson2004": "Thompson et al. 2004", "thompson2007": "Thompson et al. 2007",
    "thompson2012": "Thompson et al. 2012", "weisman_klemp1982": "Weisman and Klemp 1982", "george1960": "George 1960",
    "miller1972": "Miller 1972", "spc_ship": "SPC", "johnson_sugden2014": "Johnson and Sugden 2014",
    "cohen2007": "Cohen et al. 2007", "blanchard1998": "Blanchard 1998", "hart_korotky1991": "Hart and Korotky 1991",
    "davies1993": "Davies 1993", "nixon_allen2021": "Nixon and Allen 2021", "browning_wexler1968": "Browning and Wexler 1968",
    "doswell_rasmussen1994": "Doswell and Rasmussen 1994", "holzworth1964": "Holzworth 1964", "betts1973": "Betts 1973",
    "carson1973": "Carson 1973", "tennekes1973": "Tennekes 1973", "savitzky_golay1964": "Savitzky and Golay 1964",
    "alduchov1996": "Alduchov and Eskridge 1996", "banacos2014": "Banacos et al. 2014", "bernstein2005": "Bernstein et al. 2005",
    "birk2021": "Birk et al. 2021", "bourgouin2000": "Bourgouin 2000", "gultepe2007": "Gultepe et al. 2007",
    "huffman1988": "Huffman and Norman 1988", "libbrecht2005": "Libbrecht 2005", "mcdonough2004": "McDonough et al. 2004",
    "kuchera": "NOAA VLab (Kuchera)", "nws_windchill": "NWS 2001", "osczevski2005": "Osczevski and Bluestein 2005",
    "rauber2000": "Rauber et al. 2000", "takahashi1991": "Takahashi et al. 1991", "sczepanski2026": "Sczepanski et al. 2024",
    "benjamin2016": "Benjamin et al. 2016", "dowell2022": "Dowell et al. 2022", "hersbach2020": "Hersbach et al. 2020",
    "ncep_fnl": "NCEP 2000", "durre2018": "Durre et al. 2018", "moninger2003": "Moninger et al. 2003",
    "bryan_fritsch2002": "Bryan and Fritsch 2002", "skamarock2019": "Skamarock et al. 2019", "bolton1980": "Bolton 1980",
    "wmo_bufr": "WMO (BUFR)", "lofs_sndmod": "Halbert and Orf (sndmod)", "psu_bufkit": "Penn State BUFKIT",
    "iem_bufkit": "IEM BUFKIT archive", "nomads": "NCEP NOMADS", "uwyo": "University of Wyoming", "iem_asos": "IEM ASOS archive",
    "ou_acars": "OU ACARS archive", "nexrad_aws": "NEXRAD on AWS", "noaa_rap_aws": "RAP on AWS", "ncei_rap": "NCEI RAP/RUC archive",
}

# How the docs page groups the reference list
GROUPS = [
    ("Thermodynamics", ["doswell_rasmussen1994", "blanchard1998", "peters2023", "emanuel1994", "gilmore_wicker1998",
                        "market_dcape", "bolton1980", "weisman_klemp1982", "george1960", "miller1972",
                        "holzworth1964", "betts1973", "carson1973", "tennekes1973"]),
    ("Kinematics and storm motion", ["bunkers2000", "bunkers2014", "corfidi2003", "nixon_allen2021", "davies_jones1984",
                                     "davies_jones1990", "brandt2022", "thompson2007", "browning_wexler1968"]),
    ("Composite parameters", ["thompson2004", "thompson2012", "hart_korotky1991", "davies1993", "spc_ship",
                              "johnson_sugden2014", "cohen2007"]),
    ("Winter parameters", ["alduchov1996", "banacos2014", "bernstein2005", "birk2021", "bourgouin2000", "gultepe2007",
                           "huffman1988", "libbrecht2005", "mcdonough2004", "kuchera", "nws_windchill", "osczevski2005",
                           "rauber2000", "takahashi1991", "sczepanski2026"]),
    ("Smoothing", ["savitzky_golay1964"]),
    ("Models and data sources", ["benjamin2016", "dowell2022", "hersbach2020", "ncep_fnl", "durre2018", "moninger2003",
                                 "bryan_fritsch2002", "skamarock2019", "wmo_bufr", "lofs_sndmod", "psu_bufkit", "iem_bufkit",
                                 "nomads", "uwyo", "iem_asos", "ou_acars", "nexrad_aws", "noaa_rap_aws", "ncei_rap"]),
    ("Python packages", ["sharppy", "metpy", "siphon", "xarray", "matplotlib", "scipy", "cartopy", "pyart", "cfgrib",
                         "cdsapi", "fsspec", "bs4"]),
]

# The parameters sounding_params(...).calc() returns: (dict, key(s), parameter, units, source reference keys)
# sb_/ml_/mu_ prefixes are the surface-based, mixed-layer (lowest 100 hPa), and most-unstable parcels.
PARAMETERS = [
    ("general", "elevation, sfc_pressure", "Surface elevation and pressure", "m, hPa", []),
    ("general", "wet_bulb, rel_humidity, spec_humidity, mix_ratio", "Wet-bulb temperature, relative humidity, specific humidity, and mixing ratio at each level", "°C, %, g/kg, g/kg", ["metpy"]),
    ("general", "virt_temp, theta, theta_e", "Virtual, potential, and equivalent potential temperature at each level", "°C, K, K", ["metpy", "bolton1980"]),
    ("general", "pwat", "Precipitable water", "in", ["metpy"]),
    ("general", "frz_pt_p, frz_pt_z, wb_frz_pt_p, wb_frz_pt_z", "Freezing level and wet-bulb zero height", "hPa, m AGL", []),
    ("general", "rh_0_500 ... rh_6_9000, w_0_500 ... w_6_9000", "Layer-mean relative humidity and mixing ratio (0-0.5, 0-1, 1-3, 3-6, and 6-9 km)", "%, g/kg", []),
    ("thermo", "sbcape, mlcape, mucape, sbcin, mlcin, mucin", "CAPE and CIN, with the virtual temperature correction", "J/kg", ["sharppy", "doswell_rasmussen1994"]),
    ("thermo", "sb3cape, sb6cape (and ml, mu)", "CAPE in the lowest 3 and 6 km", "J/kg", ["sharppy"]),
    ("thermo", "sb_lcl_p/_z, sb_lfc_p/_z, sb_el_p/_z, sb_el_T, sb_mpl_p (and ml, mu)", "Lifted condensation level, level of free convection, equilibrium level, and maximum parcel level", "hPa, m AGL, °C", ["sharppy"]),
    ("thermo", "sbT_trace, sbP_trace, sbZ_trace (and ml, mu)", "Parcel temperature trace", "°C, hPa, m", ["sharppy"]),
    ("thermo", "sb_ncape, ml_ncape, mu_ncape", "Normalized CAPE (CAPE divided by the LFC-to-EL depth)", "m/s²", ["blanchard1998"]),
    ("thermo", "sb_ecape, ml_ecape, mu_ecape", "Entraining CAPE", "J/kg", ["peters2023"]),
    ("thermo", "dcape, dcin, dparcel_p, dparcel_T", "Downdraft CAPE and downdraft CIN, and the downdraft parcel trace", "J/kg, hPa, °C", ["emanuel1994", "gilmore_wicker1998", "market_dcape"]),
    ("thermo", "pbl_top, pbl_depth, pbl_lapse_rate", "Planetary boundary layer top, depth, and lapse rate (θv parcel + 0.5 K)", "hPa, m, °C/km", ["sharppy"]),
    ("thermo", "dgz, hgz", "Dendritic growth zone (-12 to -17 °C) and hail growth zone (-10 to -30 °C), as (bottom, top)", "hPa", ["sharppy"]),
    ("thermo", "lr_03km, lr_36km, lr_75", "0-3 km, 3-6 km, and 700-500 hPa lapse rates", "°C/km", ["sharppy"]),
    ("thermo", "lr_max", "Steepest 2 km lapse rate below 6 km, as (lapse rate, bottom, top)", "°C/km, hPa", ["sharppy"]),
    ("thermo", "thetae_diff", "Theta-e difference: the minimum theta-e from the 0-3 km theta-e maximum up to 7 km, minus that maximum", "K", ["cohen2007"]),
    ("thermo", "k_index, t_totals", "K-index and Total Totals", "", ["george1960", "miller1972"]),
    ("thermo", "sb_brn, ml_brn, mu_brn", "Bulk Richardson number", "", ["weisman_klemp1982"]),
    ("thermo", "ship", "Significant hail parameter", "", ["spc_ship", "johnson_sugden2014"]),
    ("thermo", "temp_adv", "Inferred temperature advection from the thermal wind, as (values, layer bounds)", "°C/hr, hPa", ["sharppy"]),
    ("kinem", "eil, eil_z, eil_depth", "Effective inflow layer (bottom, top) and its depth", "hPa, m AGL, m", ["thompson2007"]),
    ("kinem", "sm_rm, sm_lm, sm_mw, sm_u, sm_v", "Bunkers right- and left-mover motions, the 850-250 hPa mean wind, and the storm motion in use (u, v)", "kt", ["bunkers2000", "bunkers2014"]),
    ("kinem", "dtm", "Deviant tornado motion (storm motion averaged with the 0-300 m mean wind)", "kt", ["nixon_allen2021"]),
    ("kinem", "mcs", "Corfidi MCS motion, as (upshear u, v, downshear u, v)", "kt", ["corfidi2003"]),
    ("kinem", "mean_wind_pbl", "PBL mean wind", "kt", []),
    ("kinem", "shear_0_to_500 ... shear_6_to_9000, shear_0_to_3000, shear_0_to_6000, shear_eil", "Bulk wind difference over each layer", "kt", ["sharppy"]),
    ("kinem", "shear_0_to_10000", "0-10 km bulk wind difference", "kt", ["cohen2007"]),
    ("kinem", "ebwd", "Effective bulk wind difference (EIL base to half the MU EL height)", "kt", ["thompson2007"]),
    ("kinem", "srh_0_to_500 ... srh_6_to_9000, srh_0_to_3000, srh_0_to_6000, srh_eil", "Storm-relative helicity over each layer", "m²/s²", ["davies_jones1990"]),
    ("kinem", "srw, srw_0_to_500 ... srw_eil", "Storm-relative wind at each level and its layer means", "kt", []),
    ("kinem", "vort, swv, swv_perc, swv_0_to_500 ... swv_eil, swv_perc_0_to_500 ... swv_perc_eil", "Horizontal and streamwise vorticity, and streamwiseness, at each level and as layer means", "1/s, %", ["davies_jones1984", "brandt2022"]),
    ("kinem", "eil_stp", "Significant tornado parameter (effective layer, with CIN)", "", ["thompson2012"]),
    ("kinem", "eil_scp", "Supercell composite parameter (effective layer)", "", ["thompson2004"]),
    ("kinem", "ehi_0_to_1000, ehi_0_to_3000", "Energy-helicity index (MUCAPE with 0-1 and 0-3 km SRH)", "", ["hart_korotky1991", "davies1993"]),
]

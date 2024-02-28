import sounderpy as spy

def __test_remove_in_prod__():
    # This file uses real-world meteorological data as a test for the ECAPE parcel code. 
    # May be removed from repository later on if any circular dependency issues come up
    # year  = '2023' 
    # month = '04'
    # day   = '19'
    # hour  = '23'

    # # # latlon = [32.79, -96.80]
    # # # latlon = [33.01, -96.52]
    # # # latlon = [32.91, -96.57]
    # # # latlon = [34.24, -99.66]
    # latlon = [35.18, -97.44]
    # method = 'rap' 

    # clean_data = spy.get_model_data(method, latlon, year, month, day, hour)
    
    year  = '2023' 
    month = '04'
    day   = '19'
    hour  = '18'
    station = 'OUN'
    clean_data = spy.get_obs_data(station, year, month, day, hour)

    spy.build_sounding(clean_data, dark_mode=False, save=True, filename="RAOB_Norman-OK_20230419-19_ml-ia.png", parcel_highlight = ["ml_ia_ecape"], parcel_background = ["ml_ia_cape"])

__test_remove_in_prod__()
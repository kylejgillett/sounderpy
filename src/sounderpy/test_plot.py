import sys
# This ensures Python looks in your local 'src' folder first for the sounderpy package
sys.path.insert(0, './src')

import sounderpy as spy

print("Step 1: Fetching observation data...")
# We are using the same UW data parameters as before
clean_data = spy.get_obs_data(
    station='PAYA', 
    year='2026', 
    month='7', 
    day='29', 
    hour='12'
)

print("\nStep 2: Building sounding plot...")
# Passing the data into the plotting function and saving it as an image file
spy.build_sounding(
    clean_data=clean_data, 
    save=True, 
    filename='uw_test_sounding.png'
)

print("\nSuccess! The end-to-end test is complete. Check your file explorer for 'uw_test_sounding.png'.")
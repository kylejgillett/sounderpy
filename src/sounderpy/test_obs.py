# Import the newly updated function from your obs_data.py file
from .obs_data import fetch_obs

# 1. Define the test parameters
# We are using the exact date and station ID from your provided CSV file
station_id = 'PAYA'  # WMO ID for the station
year = '2026'
month = '7'
day = '29'
hour = '12'

print("Starting test for fetch_obs using the new UW API...")

try:
    # 2. Call the function
    # We set hush=False to see the summary output and clean_it=True to test our new parsing logic
    clean_data = fetch_obs(
        station=station_id, 
        year=year, 
        month=month, 
        day=day, 
        hour=hour, 
        hush=False, 
        clean_it=True
    )
    
    # 3. Verify the output
    print("\n--- Test Results ---")
    
    # Check if we got the expected dictionary back
    if isinstance(clean_data, dict) and 'p' in clean_data:
        print("Success! The function returned a dictionary of clean data.")
        
        # Print a few data points to verify our new column mapping and wind calculations
        print(f"\nNumber of data levels parsed: {len(clean_data['p'])}")
        print(f"First 3 Pressure levels (hPa): {clean_data['p'][:3]}")
        print(f"First 3 Heights (m): {clean_data['z'][:3]}")
        
        # Verify the u and v winds were calculated properly
        print(f"First 3 U-wind components (kt): {clean_data['u'][:3]}")
        print(f"First 3 V-wind components (kt): {clean_data['v'][:3]}")
        
        print("\nTest completed successfully!")
    else:
        print("Test failed: The function did not return the expected clean data dictionary.")

except Exception as e:
    # If anything goes wrong, this will catch the error and print it so we can debug
    print(f"\nTest encountered an error:\n{e}")
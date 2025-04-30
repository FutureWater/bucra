import os
import glob
import calendar
from pathlib import Path
import rasterio
import numpy as np
import pandas as pd
from IPython.core.debugger import set_trace


"""
This script produces land suitability maps by:
1) Loading suitability layers from different parameters (NDVI, Temperature, Water, etc.)
2) Applying weightings to each parameter based on its importance
3) Combining all weighted layers to create a final land suitability map
4) Creating variations for different temperature/precipitation conditions
"""

####################################################################################################
###################### Define directories, constants and file paths ################################
####################################################################################################
# Define current and parent working directories
current_wd = os.getcwd()
parent_wd = os.path.dirname(current_wd)
angola_wd = "/Users/thomasfuturewater/FutureWater Dropbox/Team/Projects/Completed/2019/2019019_G4AW_MavoDiami_Angola/Data/2019019_MavoDiami_LV/2019019_MavoDiami"

# Gets province from subprocess in 000_Run_All.py
PROVINCE_NAME = os.environ.get("PROVINCE")
PROVINCE_NAME = "Zaire"                         # dummy variable for testing.

# Define other folders
DATA_DIR = os.path.join(angola_wd, "01_Data")
GIS_DIR = os.path.join(angola_wd, "02_GIS")     # Directory with shapefiles
RESULTS_DIR = os.path.join(parent_wd, "04_Results", PROVINCE_NAME)
TEMP_DIR = os.path.join(parent_wd, "05_Temp")
HHS_DATA_DIR = os.path.join(DATA_DIR, "Soil_Hydraulic_Properties")
LS_RESULTS = os.path.join(RESULTS_DIR, '_LS_Results')

# Define constants
RES = 250  # Resolution in meters
NO_DATA_VALUE = 0  # No data value
LOCAL_PROJECTION = "EPSG:32733"

T_PERC = [0.75, 0.95]  # Temperature percentiles
T_PERC_NAMES = ["Warmer", "Much Warmer"]  # Names for temperature scenarios
P_PERC = [0.05, 0.25]  # Precipitation percentiles
P_PERC_NAMES = ["Much Drier", "Drier"]  # Names for precipitation scenarios

# Create a mock cropping calendar (replace with actual data loading)
CROPPING_CAL = pd.read_csv(os.path.join(current_wd, "Cropping_calendar_A.csv"))
PARAMS = pd.read_csv(os.path.join(current_wd, "Parameters.csv"))

####################################################################################################
############################## Produce Suitability maps ############################################
####################################################################################################
# Create output directory
new_dir = os.path.join(RESULTS_DIR, "_LS_Results", "Weighted")
os.makedirs(new_dir, exist_ok=True)

# Get input files for each parameter
input_files_ndvi = glob.glob(os.path.join(RESULTS_DIR,
                                          "_LS_Results", "NDVI", "*.tif"))
input_files_temperature = glob.glob(os.path.join(RESULTS_DIR,
                                                 "_LS_Results", "Temperature", "**", "*.tif"), recursive=True)
input_files_water = glob.glob(os.path.join(RESULTS_DIR,
                                           "_LS_Results", "Water", "**", "*.tif"), recursive=True)
# input_files_hhs = glob.glob(os.path.join(RESULTS_DIR,
#                                          "_LS_Results", "Soil_Hydraulic_Properties", "*.tif"))
input_files_snc = glob.glob(os.path.join(RESULTS_DIR,
                                         "_LS_Results", "Soil_Nutrient_Content", "*.tif"))
input_files_slope = glob.glob(os.path.join(RESULTS_DIR,
                                           "_LS_Results", "Elevation", "*.tif"))
# set_trace()
# Function to find file matching a pattern


def find_file(file_list, pattern):
    matching_files = [f for f in file_list if pattern in f]
    return matching_files[0] if matching_files else None


# Process each crop
print("Processing land suitability maps...")
for _, crop_row in CROPPING_CAL.iterrows():
    crop = crop_row['Crop']
    start_month = crop_row['Start_growing_season']
    end_month = crop_row['End_growing_season']

    months = list(range(1, 13))
    month_abbrs = [calendar.month_abbr[m] for m in months]
    start_month = month_abbrs[start_month-1]
    end_month = month_abbrs[end_month-1]

    # Process each temperature scenario
    for t_idx, t_perc in enumerate(T_PERC + [None]):
        if t_perc is None:
            warm = "Average"
            per_t = "Monthly_Mean"
        else:
            warm = T_PERC_NAMES[t_idx]
            per_t = f"Tmax_{str(t_perc).replace('.', '')}perc"

        # Process each precipitation scenario
        for p_idx, p_perc in enumerate(P_PERC + [None]):
            if p_perc is None:
                dry = "Average"
                per_p = "Mean_Monthly"
            else:
                dry = P_PERC_NAMES[p_idx]
                per_p = f"{str(p_perc).replace('.', '')}perc"

            print(f"    Processing {crop} for {dry}/{warm} conditions...")

            # Create output filename
            file_name = f"Land_Suitability_{dry}_{warm}_{crop}_{start_month}-{end_month}.tif"

            # Find and load each parameter file with appropriate weighting
            # NDVI
            ndvi_pattern = f"{start_month}-{end_month}"
            ndvi_file = find_file(input_files_ndvi, ndvi_pattern)
            if not ndvi_file:
                print(
                    f"  NDVI file not found for {crop}, {start_month}-{end_month}")
                continue

            # Temperature
            temp_pattern = f"{crop}_{start_month}-{end_month}"
            temp_files = [
                f for f in input_files_temperature if per_t in f and temp_pattern in f]
            if not temp_files:
                print(f"    Temperature file not found for {crop}, {per_t}")
                continue

            # Water
            water_pattern = f"{crop}_{start_month}-{end_month}_higher_"
            water_files = [
                f for f in input_files_water if per_p in f and water_pattern in f]
            if not water_files:
                print(f"    Water file not found for {crop}, {per_p}")
                continue

            # Find soil and slope files
            # ksat_file = find_file(input_files_hhs, "Ksat")
            # wcavail_file = find_file(input_files_hhs, "WCavail")
            potassium_file = find_file(input_files_snc, "_K_")
            phosphorus_file = find_file(input_files_snc, "_P_")
            slope_file = find_file(input_files_slope, "lower")

            if not all([potassium_file, phosphorus_file, slope_file]):  # ksat_file, wcavail_file,
                print(f"    One or more soil/slope files not found")
                continue

            # Get weights from parameters table
            ndvi_weight = PARAMS.loc[PARAMS['Parameter']
                                     == 'NDVI', 'Weight'].values[0]
            temp_weight = PARAMS.loc[PARAMS['Parameter']
                                     == 'Temperature', 'Weight'].values[0]
            water_weight = PARAMS.loc[PARAMS['Parameter']
                                      == 'Water', 'Weight'].values[0]
            ksat_weight = PARAMS.loc[PARAMS['Parameter']
                                     == 'Ksat', 'Weight'].values[0]
            wcavail_weight = PARAMS.loc[PARAMS['Parameter']
                                        == 'WCavail', 'Weight'].values[0]
            potassium_weight = PARAMS.loc[PARAMS['Parameter']
                                          == 'Potasium', 'Weight'].values[0]
            phosphorus_weight = PARAMS.loc[PARAMS['Parameter']
                                           == 'Phosphorus', 'Weight'].values[0]
            slope_weight = PARAMS.loc[PARAMS['Parameter']
                                      == 'Slope', 'Weight'].values[0]

            # Read and apply weights to each layer
            with rasterio.open(ndvi_file) as src:
                ndvi_data = src.read(1) * ndvi_weight
                output_meta = src.meta.copy()

            with rasterio.open(temp_files[0]) as src:
                temp_data = src.read(1) * temp_weight

            with rasterio.open(water_files[0]) as src:
                water_data = src.read(1) * water_weight

            # with rasterio.open(ksat_file) as src:
            #     ksat_data = src.read(1) * ksat_weight

            # with rasterio.open(wcavail_file) as src:
            #     wcavail_data = src.read(1) * wcavail_weight

            with rasterio.open(potassium_file) as src:
                potassium_data = src.read(1) * potassium_weight

            with rasterio.open(phosphorus_file) as src:
                phosphorus_data = src.read(1) * phosphorus_weight

            with rasterio.open(slope_file) as src:
                slope_data = src.read(1) * slope_weight

            # Sum all weighted layers
            ls_data = (ndvi_data + temp_data + water_data +
                       potassium_data + phosphorus_data + slope_data)
            # ksat_data + wcavail_data +

            # Save the result
            output_path = os.path.join(new_dir, file_name)
            output_meta.update({
                'dtype': 'float32',
                'count': 1,
                'nodata': NO_DATA_VALUE
            })
            with rasterio.open(output_path, 'w', **output_meta) as dst:
                dst.write(ls_data, 1)

print("Land suitability map production complete!")

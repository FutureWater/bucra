import os
import glob
import calendar
import rasterio
import numpy as np
from rasterio.merge import merge
import pandas as pd
import matplotlib.pyplot as plt

"""
This script calculates temperature limits for different crops by:
1) Processing minimum and maximum temperature files.
2) Applying crop-specific base and upper temperature thresholds.
    a. For tmin, the monthly mean min temp is used.
    b. For tmax, the 75 and 95 percentile are used.
3) Creating temperature suitability maps based on weighted averages
"""

####################################################################################################
# Functions
def fuzzy_membership(raster, function_type, parameters):
    result = raster.copy().astype(float)

    if function_type == "Increasing":
        # Two parameters: [a, b]
        # below a = 0, between a-b = 0 to 1, above b = 1
        a = parameters[0]
        b = parameters[1]
        
        result[raster <= a] = 0
        result[raster >= b] = 1
        mask = (a < raster) & (raster < b)
        result[mask] = (raster[mask] - a) / (b-a)
        
    elif function_type == "Decreasing":
        # Two parameters: [a, b]
        # below a = 1, between a-b = 1 to 0, above b = 0
        a = parameters[0]
        b = parameters[1]
        
        result[raster <= a] = 1
        result[raster >= b] = 0
        mask = (a < raster) & (raster < b)
        result[mask] = 1 - (raster[mask] - a) / (b - a)
    
    elif function_type == "Triangular":
        # Three parameters: [a, b, c]
        # below a = 0, a-b = 0 to 1, b-c = 1 to 0, above c = 0
        a = parameters[0]
        b = parameters[1]  # peak point
        c = parameters[2]

        result[raster <= a] = 0
        result[raster >= c] = 0

        mask1 = (a < raster) & (raster < b)
        mask2 = (b < raster) & (raster < c)
        result[mask1] = (raster[mask1] - a) / (b - a)
        result[mask2] = 1 - (raster[mask2] - b) / (c - b)
    
    elif function_type == "Trapezoidal":
        # Four parameters: [a, b, c, d]
        # below a = 0, a-b = 0 to 1, b-c = 1, c-d = 1 to 0, above d = 0
        a = parameters[0]
        b = parameters[1]  # start of plateau
        c = parameters[2]  # end of plateau
        d = parameters[3]

        result[raster <= a] = 0
        result[raster >= d] = 0

        mask1 = (a < raster) & (raster < b)
        mask2 = (b < raster) & (raster < c)
        mask3 = (c < raster) & (raster <d)
        result[mask1] = (raster[mask1] - a) / (b - a)
        result[mask2] = 1
        result[mask3] = 1 - (raster[mask3] -c) / (d - c)
    
    if np.all(result == raster):
        raise ValueError("Result raster is the same as input raster. Fuzzy logic not applied")
    
    return result

def plot_raster(raster):
    plt.imshow(np.ma.masked_where(raster == NO_DATA_VALUE, raster), cmap='viridis')
    plt.colorbar()
####################################################################################################
###################### Define directories, constants and file paths ################################
####################################################################################################
# Define current and parent working directories
current_wd = os.getcwd()
parent_wd = os.path.dirname(current_wd)
base_wd = os.path.dirname(os.path.dirname(parent_wd))
angola_wd = "/Users/thomasfuturewater/FutureWater Dropbox/Team/Projects/Completed/2019/2019019_G4AW_MavoDiami_Angola/Data/2019019_MavoDiami_LV/2019019_MavoDiami"

# Gets province from subprocess in 000_Run_All.py
PROVINCE_NAME = os.environ.get("PROVINCE")
if not os.environ.get("PROVINCE"):
    PROVINCE_NAME = "Sharkia"                           # dummy variable for testing.

# Define other folders
DATA_DIR = os.path.join(parent_wd, "01_Data")
GIS_DIR = os.path.join(base_wd, "GIS")     # Directory with shapefiles
RESULTS_DIR = os.path.join(parent_wd, "04_Results", PROVINCE_NAME)
os.makedirs(RESULTS_DIR, exist_ok=True)
TEMP_DIR = os.path.join(parent_wd, "05_Temp")

# Define constants
RES = 250  # Resolution in meters
NO_DATA_VALUE = -9999.0  # No data value
LOCAL_PROJ = "EPSG:32636"  # Local projection

# Define constants
RES = 250  # Resolution in meters
NO_DATA_VALUE = -9999.0  # No data value
VARIABLES = ['P75', 'P95', 'Mean']  # Percentiles to process

# Load cropping calendar and parameters
CROPPING_CAL = pd.read_csv(os.path.join(current_wd, "Cropping_calendar.csv"), sep = ",")
PARAMS = pd.read_csv(os.path.join(current_wd, "Parameters.csv"))


# Get input temperature files
input_files_tmin = glob.glob(os.path.join(
    RESULTS_DIR, "Temperature", "Tmin*", "**", "*.tif"), recursive=True)
input_files_tmax = glob.glob(os.path.join(
    RESULTS_DIR, "Temperature", "Tmax*", "**", "*.tif"), recursive=True)

# Create output directory
NEW_RESULTS_SUBDIR = os.path.join(RESULTS_DIR,
                                  "_LS_Results", "Temperature")

####################################################################################################
###################### Process parameters and create limit maps ####################################
####################################################################################################
# Process each percentile (and Monthly_Mean)
for p_idx, per in enumerate(VARIABLES):
    print(f"Processing temperature limits for {per}...")

    # Process each crop
    for idx in CROPPING_CAL.index:
        crop_data = CROPPING_CAL.loc[idx]
        crop_name = crop_data["Crop"]
        start_month = int(crop_data["Start_growing_season"])
        end_month = int(crop_data["End_growing_season"])

        # Determine months in growing season
        if start_month > end_month:
            months = list(range(start_month, 13)) + \
                list(range(1, end_month + 1))
        else:
            months = list(range(start_month, end_month + 1))

        # Get weights for each month
        weights = [float(crop_data[f"w_{i+1}"]) for i in range(len(months))]

        # Get temperature thresholds
        t_base = crop_data["T_Base"]
        t_optimal_start = crop_data["T_Optimal_Start"]
        t_optimal_end = crop_data["T_Optimal_End"]
        t_upper = crop_data["T_Upper"]

        # Get month abbreviations for filtering files
        month_abbrs = [calendar.month_abbr[m] for m in months]

        # Filter tmin files for these months with Monthly_Mean
        tmin_files = []
        for month_abbr in month_abbrs:
            tmin_files.extend([f for f in input_files_tmin
                               if month_abbr in os.path.basename(f) in f])

        # Filter tmax files for these months with specified percentiles
        tmax_files = []
        for month_abbr in month_abbrs:
            tmax_files.extend([f for f in input_files_tmax
                               if month_abbr in os.path.basename(f)
                               and per in f])

        # Make sure files are in the same order as the months
        tmin_files_ordered = []
        tmax_files_ordered = []
        for month_abbr in month_abbrs:
            for file in tmin_files:
                if month_abbr in os.path.basename(file):
                    tmin_files_ordered.append(file)
                    break

            for file in tmax_files:
                if month_abbr in os.path.basename(file):
                    tmax_files_ordered.append(file)
                    break

        # Initialize accumulation arrays
        with rasterio.open(tmin_files_ordered[0]) as src:
            shape = src.shape
            tmin_w_data = np.zeros(shape, dtype=np.float32)
            tmax_w_data = np.zeros(shape, dtype=np.float32)
            tmin_valid = np.zeros(shape, dtype=bool)
            tmax_valid = np.zeros(shape, dtype=bool)
            meta = src.meta.copy()

        # Get weighted average of Tmin for whole growing period. 
        for i, tmin_file in enumerate(tmin_files_ordered):
            with rasterio.open(tmin_file) as src:
                tmin_data = src.read(1)
                mask = tmin_data != NO_DATA_VALUE
                tmin_data[mask] *= weights[i]   # Only multiply valid cells in the month.
                tmin_data[~mask] = 0            # Set invalid cells to zero.
                tmin_w_data += tmin_data        # Sum all weighted values to get the average.
                tmin_valid |= mask              # Update total validity mask.
                
        # Set cells that were invalid in all rasters to no data value.
        tmin_w_data[~tmin_valid] = NO_DATA_VALUE

        # Get weighted average of Tmax for whole growing period. Add zero for invalid cells.
        for i, tmax_file in enumerate(tmax_files_ordered):
            with rasterio.open(tmax_file) as src:
                tmax_data = src.read(1)
                mask = tmax_data != NO_DATA_VALUE
                tmax_data[mask] *= weights[i]   # Only multiply valid cells.
                tmax_data[~mask] = 0            # Set invalid cells to zero.
                tmax_w_data += tmax_data        # Sum all weighted values to get the average.
                tmax_valid |= mask              # Update total validity mask.

        # Set cells that were invalid in all rasters to no data value.
        tmax_w_data[~tmax_valid] = NO_DATA_VALUE

        # Create temperature suitability maps
        if tmin_w_data is not None and tmax_w_data is not None:
            # Apply fuzzy logic
            temperature_lower_limit = fuzzy_membership(tmin_w_data, 'Increasing', [t_base, t_optimal_start])
            temperature_upper_limit = fuzzy_membership(tmax_w_data, 'Decreasing', [t_optimal_end, t_upper])
            temperature_limit = (temperature_lower_limit + temperature_upper_limit) * 0.5
            temperature_limit = temperature_limit.astype(np.float32)
            temperature_limit[~tmax_valid] = NO_DATA_VALUE

            # Save output
            season_label = f"{month_abbrs[0]}-{month_abbrs[-1]}"
            folder_name = f"Tmax_{per}"
            output_path = os.path.join(
                NEW_RESULTS_SUBDIR, folder_name,
                f"Temperature_suitability_{crop_name}_{season_label}.tif")
            os.makedirs(os.path.join(NEW_RESULTS_SUBDIR,
                        folder_name), exist_ok=True)

            # Update metadata for output raster. Set all nodata values to 0.
            meta.update({
                'dtype': 'float32',
                'count': 1,
                'nodata': -9999
            })

            # Write output raster
            with rasterio.open(output_path, 'w', **meta) as dst:
                dst.write(temperature_limit, 1)

            print(
                f"  Created temperature suitability map for {crop_name}, {season_label}, {per}")
        else:
            print(
                f"  Missing temperature data for {crop_name}, {season_label}, {per}")

print("Temperature limits calculation complete!")

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
LOCAL_PROJECTION = "EPSG:32733"

# Define constants
RES = 250  # Resolution in meters
NO_DATA_VALUE = -9999.0  # No data value
T_PERC = [0.75, 0.95]  # Percentiles to process

# Load cropping calendar and parameters
CROPPING_CAL = pd.read_csv(os.path.join(current_wd, "Cropping_calendar_A.csv"))
PARAMS = pd.read_csv(os.path.join(current_wd, "Parameters.csv"))


# Get input temperature files
input_files_tmin = glob.glob(os.path.join(
    RESULTS_DIR, "Temperature", "tmin", "**", "*.tif"), recursive=True)
input_files_tmax = glob.glob(os.path.join(
    RESULTS_DIR, "Temperature", "tmax", "**", "*.tif"), recursive=True)

# Create output directory
NEW_RESULTS_SUBDIR = os.path.join(RESULTS_DIR,
                                  "_LS_Results", "Temperature")

####################################################################################################
###################### Process parameters and create limit maps ####################################
####################################################################################################
# Process each percentile (and Monthly_Mean)
for p_idx, p_val in enumerate(T_PERC + ["Monthly_Mean"]):
    print(f"Processing temperature limits for {p_val}...")
    # Handle percentile naming and folder setup
    if p_val == "Monthly_Mean":
        per = "Monthly_Mean"
        folder_name = "Tmax_Monthly_Mean"
        start_name = "Tmax_MM_between"
    else:
        per = f"{str(p_val).replace('.', '')}perc"
        folder_name = f"Tmax_{per}"
        start_name = f"Tmax_{per}_between"

    # Create output directory
    os.makedirs(os.path.join(NEW_RESULTS_SUBDIR, folder_name), exist_ok=True)

    # Process each crop
    for idx in CROPPING_CAL.index:
        crop_data = CROPPING_CAL.loc[idx]
        crop_name = crop_data["Crop"]
        start_month = crop_data["Start_growing_season"]
        end_month = crop_data["End_growing_season"]

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
        t_upper = crop_data["T_Upper"]

        # Get month abbreviations for filtering files
        month_abbrs = [calendar.month_abbr[m] for m in months]

        # Filter tmin files for these months with Monthly_Mean
        tmin_files = []
        for month_abbr in month_abbrs:
            tmin_files.extend([f for f in input_files_tmin
                               if month_abbr in os.path.basename(f)
                               and "Monthly_Mean" in f])

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

        # Sum all tmin data. Add zero for invalid cells.
        for i, tmin_file in enumerate(tmin_files_ordered):
            with rasterio.open(tmin_file) as src:
                tmin_data = src.read(1)
                mask = tmin_data != NO_DATA_VALUE
                tmin_data[mask] *= weights[i]  # only multiply valid cells.
                tmin_data[~mask] = 0
                tmin_w_data += tmin_data
                tmin_valid |= mask  # Update total validity mask.
        # Set cells that were invalid in all rasters to no data value.
        tmin_w_data[~tmin_valid] = NO_DATA_VALUE

        # Sum all tmax data. Add zero for invalid cells.
        for i, tmax_file in enumerate(tmax_files_ordered):
            with rasterio.open(tmax_file) as src:
                tmax_data = src.read(1)
                mask = tmax_data != NO_DATA_VALUE
                tmax_data[mask] *= weights[i]  # only multiply valid cells.
                tmax_data[~mask] = 0
                tmax_w_data += tmax_data
                tmax_valid |= mask

        # Set cells that were invalid in all rasters to no data value.
        tmax_w_data[~tmax_valid] = NO_DATA_VALUE

        # Create temperature suitability maps
        if tmin_w_data is not None and tmax_w_data is not None:
            # Check if tmin is higher than base temperature
            tmin_w_higher_tbase = (tmin_w_data > t_base).astype(np.uint8)

            # Check if tmax is lower than upper temperature
            tmax_w_lower_tupper = (tmax_w_data < t_upper).astype(np.uint8)

            # Combined suitability (both conditions must be met)
            two_limits_w = (tmin_w_higher_tbase) & (
                tmax_w_lower_tupper) & tmax_valid & tmin_valid
            two_limits_w.astype(np.uint8)

            # Save output
            season_label = f"{month_abbrs[0]}-{month_abbrs[-1]}"
            output_path = os.path.join(
                NEW_RESULTS_SUBDIR, folder_name,
                f"Temp_between_{t_base}_and_{t_upper}°C_{crop_name}_{season_label}.tif")

            # Update metadata for output raster. Set all nodata values to 0.
            meta.update({
                'dtype': 'uint8',
                'count': 1,
                'nodata': 0
            })

            # Write output raster
            with rasterio.open(output_path, 'w', **meta) as dst:
                dst.write(two_limits_w, 1)

            print(
                f"  Created temperature suitability map for {crop_name}, {season_label}, {per}")
        else:
            print(
                f"  Missing temperature data for {crop_name}, {season_label}, {per}")

print("Temperature limits calculation complete!")

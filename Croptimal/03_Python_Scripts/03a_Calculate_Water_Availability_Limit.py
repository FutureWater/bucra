import os
import glob
from pathlib import Path
import calendar
import rasterio
import numpy as np
import pandas as pd
from rasterio.merge import merge
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable


"""
This script calculates water availability indices for different crops by:
1) Calculating total ETc (crop evapotranspiration) for specific growth periods
2) Calculating total precipitation for those same periods
3) Creating a water availability index (P/ETc) to identify suitable areas

Kc values come from FAO: "Guidelines for computing crop water requirements - FAO Paper 56"
"""

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
P_PERC = [0.05, 0.25]  # Percentiles to process

# Load cropping calendar and parameters
CROPPING_CAL = pd.read_csv(os.path.join(current_wd, "Cropping_calendar_A.csv"))
PARAMS = pd.read_csv(os.path.join(current_wd, "Parameters.csv"))

# Get input filepaths
INPUT_FILES_ETC = glob.glob(os.path.join(
    RESULTS_DIR, "Ref_ET", "Mean_Monthly", "*.tif"))
INPUT_FILES_P = glob.glob(os.path.join(
    RESULTS_DIR, "Rainfall", "**", "*.tif"), recursive=True)

# Create output directories
NEW_RESULTS_SUBDIR = os.path.join(RESULTS_DIR, "_LS_Results")
for folder in ["Water", "ETc", "Rainfall"]:
    os.makedirs(os.path.join(NEW_RESULTS_SUBDIR, folder), exist_ok=True)

###################### Process each percentile and mean ######################
for p_idx, p_val in enumerate(P_PERC + ["Mean_Monthly"]):
    print(f"Processing water availability with {p_val} Precipitation...")
    # Handle percentile naming
    if p_val == "Mean_Monthly":
        per = "Mean_Monthly"
    else:
        per = f"{str(p_val).replace('.', '')}perc"

    # Create subdirectories
    os.makedirs(os.path.join(NEW_RESULTS_SUBDIR,
                "Rainfall", per), exist_ok=True)
    os.makedirs(os.path.join(NEW_RESULTS_SUBDIR, "Water", per), exist_ok=True)

    ###################### Process each crop ######################
    # (iterate over rows of cropping calendar)
    for _, crop_data in CROPPING_CAL.iterrows():
        crop_name = crop_data["Crop"]
        start_month = crop_data["Start_growing_season"]
        end_month = crop_data["End_growing_season"]

        # Create list of month numbers in the season, and use these to extract kc values from cropping calendar.
        if start_month > end_month:
            # Handle seasons that cross year boundary.
            months = list(range(start_month, 13)) + \
                list(range(1, end_month + 1))
        else:
            months = list(range(start_month, end_month + 1))

        # Set seasonal values.
        kc = [float(crop_data[f"kc_{m}"]) for m in months]
        month_abbrs = [calendar.month_abbr[m] for m in months]
        season_label = f"{calendar.month_abbr[start_month]}-{calendar.month_abbr[end_month]}"

        # Find ET and P filepaths for these months
        et_files = [f for f in INPUT_FILES_ETC if any(
            m in os.path.basename(f) for m in month_abbrs)]
        p_files = [f for f in INPUT_FILES_P if per in f and any(
            m in os.path.basename(f) for m in month_abbrs)]

        ###################### Read and sum data ######################
        # Read ET rasters, multiply with monthly kc value and sum all rasters into a single seasonal ET raster.
        # Initialize accumulation arrays. Set all cells to value of zero.
        with rasterio.open(et_files[0]) as src:
            shape = src.read(1).shape
            meta = src.meta.copy()
            meta.update(dtype='float32', nodata=NO_DATA_VALUE)

        et_sum = np.zeros(shape, dtype=np.float32)
        p_sum = np.zeros(shape, dtype=np.float32)
        et_valid = np.zeros(shape, dtype=bool)
        p_valid = np.zeros(shape, dtype=bool)
        # --- ET SUM ---
        for et_file in et_files:
            with rasterio.open(et_file) as src:
                # Read et raster data.
                et_data = src.read(1)
                mask = et_data != NO_DATA_VALUE

                month_abbr = os.path.basename(et_file)[4:7]
                kc_factor = kc[month_abbrs.index(month_abbr)]

                # Only add data summed raster for valid cells. All other cells get 0  added.
                et_data[mask] *= kc_factor
                et_data[~mask] = 0  # exclude nodata from sum
                et_sum += et_data
                et_valid |= mask

        # Set all cells that never were valid to no data value.
        et_sum[~et_valid] = NO_DATA_VALUE

        # --- P SUM ---
        for p_file in p_files:
            with rasterio.open(p_file) as src:
                p_data = src.read(1)
                mask = p_data != NO_DATA_VALUE
                p_data[~mask] = 0
                p_sum += p_data
                p_valid |= mask

        p_sum[~p_valid] = NO_DATA_VALUE

        # --- WATER AVAILABILITY ---
        water = np.full(shape, NO_DATA_VALUE, dtype=np.float32)
        valid_mask = et_valid & p_valid
        water[valid_mask] = p_sum[valid_mask] / et_sum[valid_mask]

        # --- WATER LIMIT SUITABILITY ---
        limit = float(PARAMS["Limit"][2])
        # water_limit = np.full(shape, NO_DATA_VALUE, dtype=np.float32)
        water_limit = (water > limit) & valid_mask
        water_limit.astype(np.uint8)

        # unsuitable_mask = (water < limit) & valid_mask
        # water_limit[suitable_mask] = water[suitable_mask]
        # water_limit[unsuitable_mask] = 0

        # ############################# Plot the rasters ##################################
        # fig, ((ax3, ax2), (ax4, ax1)) = plt.subplots(2, 2, figsize=(10, 10))

        # # Plot your data on each subplot
        # im1 = ax1.imshow(np.ma.masked_where(
        #     water_limit == NO_DATA_VALUE, water_limit))
        # im2 = ax2.imshow(np.ma.masked_where(
        #     et_sum == NO_DATA_VALUE, et_sum))
        # im3 = ax3.imshow(np.ma.masked_where(
        #     p_sum == NO_DATA_VALUE, p_sum))
        # im4 = ax4.imshow(np.ma.masked_where(
        #     water == NO_DATA_VALUE, water))

        # # Create colorbars with the same height as the plots
        # for ax, im, title in zip([ax1, ax2, ax3, ax4],
        #                          [im1, im2, im3, im4],
        #                          ["Water Limit Data", "ET Sum Data", "P Sum Data", "Water Data"]):
        #     divider = make_axes_locatable(ax)
        #     cax = divider.append_axes("right", size="5%", pad=0.05)
        #     fig.colorbar(im, cax=cax)
        #     ax.set_title(title)
        # plt.tight_layout()
        # plt.show()

        # # Plot difference in masks
        # plt.figure()
        # plt.imshow(valid_mask != suitable_mask)
        # plt.title("Mask mismatch between ET and Suitability")
        # plt.colorbar()
        # plt.show()

        ######################## Save outputs ########################
        # Create season lables
        season_label = f"{calendar.month_abbr[start_month]}-{calendar.month_abbr[end_month]}"

        # Save rasters
        out_paths = {
            "ETc": os.path.join(NEW_RESULTS_SUBDIR, "ETc", f"ETc_{crop_name}_{season_label}.tif"),
            "P": os.path.join(NEW_RESULTS_SUBDIR, "Rainfall", per, f"P_{per}_{crop_name}_{season_label}.tif"),
            "Water": os.path.join(NEW_RESULTS_SUBDIR, "Water", per,
                                  f"Water_availability_{per}_{crop_name}_{season_label}_higher_than_{limit}.tif")
        }

        for label, arr in zip(out_paths, [et_sum, p_sum, water_limit]):
            if label == 'Water':
                meta.update({'dtype': 'uint8',
                             'nodata': 0})
            with rasterio.open(out_paths[label], "w", **meta) as dst:
                dst.write(arr, 1)

        print(f"    Processed {crop_name} in {season_label} for {per}")


print("Water availability analysis complete!")

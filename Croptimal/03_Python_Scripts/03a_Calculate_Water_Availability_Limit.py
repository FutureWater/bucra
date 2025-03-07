import os
import glob
import calendar
import rasterio
import numpy as np
from rasterio.merge import merge
from pathlib import Path

"""
This script calculates water availability indices for different crops by:
1) Calculating total ETc (crop evapotranspiration) for specific growth periods
2) Calculating total precipitation for those same periods
3) Creating a water availability index (P/ETc) to identify suitable areas

Kc values come from FAO: "Guidelines for computing crop water requirements - FAO Paper 56"
"""

# Define constants
DATA_DIR = "your_data_directory"  # Replace with actual path
RESULTS_DIR = "your_results_directory"  # Replace with actual path
PROVINCE_NAME = "your_province_name"  # Replace with actual province name
P_PERC = [0.05, 0.25]  # Percentiles to process

# Create a mock cropping calendar (replace with actual data loading)
# Structure needs: Crop, Start_growing_season, End_growing_season, Kc values for each month
cropping_cal = [
    {"Crop": "Wheat", "Start_growing_season": 11, "End_growing_season": 4,
     "Kc_1": 0.3, "Kc_2": 0.3, "Kc_3": 0.3, "Kc_4": 0.4, "Kc_5": 0.8, "Kc_6": 0.3,
     "Kc_7": 0.3, "Kc_8": 0.3, "Kc_9": 0.3, "Kc_10": 0.3, "Kc_11": 0.3, "Kc_12": 0.3},
    # Add more crops as needed
]

# Parameters dictionary (replace with actual data loading)
params = {"Parameter": ["Water"], "Limit": [0.5]}

# Get input files
input_files_etc = glob.glob(os.path.join(
    RESULTS_DIR, PROVINCE_NAME, "Ref_ET", "Mean_Monthly", "*.tif"))
input_files_p = glob.glob(os.path.join(
    RESULTS_DIR, PROVINCE_NAME, "Rainfall", "**", "*.tif"), recursive=True)

# Create output directories
new_dir = os.path.join(RESULTS_DIR, PROVINCE_NAME, "_LS_Results")
for folder in ["Water", "ETc", "Rainfall"]:
    os.makedirs(os.path.join(new_dir, folder), exist_ok=True)

# Process each percentile (and mean)
for p_idx, p_val in enumerate(P_PERC + ["Mean_Monthly"]):
    # Handle percentile naming
    if p_val == "Mean_Monthly":
        per = "Mean_Monthly"
    else:
        per = f"{str(p_val).replace('.', '')}perc"

    # Create subdirectories
    os.makedirs(os.path.join(new_dir, "Rainfall", per), exist_ok=True)
    os.makedirs(os.path.join(new_dir, "Water", per), exist_ok=True)

    # Process each crop
    for crop_data in cropping_cal:
        crop_name = crop_data["Crop"]
        start_month = crop_data["Start_growing_season"]
        end_month = crop_data["End_growing_season"]

        # Handle season that spans year boundary
        if start_month > end_month:
            months = list(range(start_month, 13)) + \
                list(range(1, end_month + 1))
            kc = [crop_data[f"Kc_{m}"] for m in (
                list(range(start_month, 13)) + list(range(1, end_month + 1)))]
        else:
            months = list(range(start_month, end_month + 1))
            kc = [crop_data[f"Kc_{m}"]
                  for m in range(start_month, end_month + 1)]

        # Get month abbreviations for filtering
        month_abbrs = [calendar.month_abbr[m] for m in months]

        # Find ET files for these months
        et_files = []
        for month_abbr in month_abbrs:
            et_files.extend(
                [f for f in input_files_etc if month_abbr in os.path.basename(f)])

        # Find P files for these months with correct percentile
        p_files = []
        for month_abbr in month_abbrs:
            p_files.extend(
                [f for f in input_files_p if month_abbr in os.path.basename(f) and per in f])

        # Read and sum ET rasters
        et_sum_data = None
        for et_file in et_files:
            with rasterio.open(et_file) as src:
                et_data = src.read(
                    1) * kc[month_abbrs.index(os.path.basename(et_file)[:3])]
                if et_sum_data is None:
                    et_sum_data = et_data
                    et_meta = src.meta.copy()
                else:
                    et_sum_data += et_data

        # Read and sum P rasters
        p_sum_data = None
        for p_file in p_files:
            with rasterio.open(p_file) as src:
                p_data = src.read(1)
                if p_sum_data is None:
                    p_sum_data = p_data
                    p_meta = src.meta.copy()
                else:
                    p_sum_data += p_data

        # Calculate water availability index
        water_data = p_sum_data / et_sum_data

        # Apply threshold to create binary suitability map
        limit = float(params["Limit"][0])  # Get water availability threshold
        water_limit_data = (water_data > limit).astype('uint8')

        # Save outputs
        season_label = f"{calendar.month_abbr[start_month]}-{calendar.month_abbr[end_month]}"

        # Save ET sum
        et_path = os.path.join(
            new_dir, "ETc", f"ETc_{crop_name}_{season_label}.tif")
        with rasterio.open(et_path, 'w', **et_meta) as dst:
            dst.write(et_sum_data, 1)

        # Save P sum
        p_path = os.path.join(new_dir, "Rainfall", per,
                              f"P_{crop_name}_{season_label}.tif")
        with rasterio.open(p_path, 'w', **p_meta) as dst:
            dst.write(p_sum_data, 1)

        # Save water limit
        water_path = os.path.join(
            new_dir, "Water", per,
            f"Water_rain_{per}_{crop_name}_{season_label}_higher_{limit}_{PROVINCE_NAME}.tif")
        with rasterio.open(water_path, 'w', **p_meta) as dst:
            dst.write(water_limit_data, 1)

        print(f"Processed {crop_name} for {per}")

print("Water availability analysis complete!")

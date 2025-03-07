import os
import glob
import calendar
import rasterio
import numpy as np
from rasterio.merge import merge

"""
This script calculates temperature limits for different crops by:
1) Processing minimum and maximum temperature files
2) Applying crop-specific base and upper temperature thresholds
3) Creating temperature suitability maps based on weighted averages
"""

# Define constants
DATA_DIR = "your_data_directory"  # Replace with actual path
RESULTS_DIR = "your_results_directory"  # Replace with actual path
PROVINCE_NAME = "your_province_name"  # Replace with actual province name
T_PERC = [0.05, 0.25]  # Percentiles to process

# Create a mock cropping calendar (replace with actual data loading)
# Structure needs: Crop, Start_growing_season, End_growing_season, T_Base, T_Upper, and weights
cropping_cal = [
    {
        "Crop": "Wheat",
        "Start_growing_season": 11,
        "End_growing_season": 4,
        "T_Base": 5.0,
        "T_Upper": 35.0,
        "w_1": 0.2, "w_2": 0.2, "w_3": 0.2, "w_4": 0.2, "w_5": 0.2, "w_6": 0.0
    },
    # Add more crops as needed
]

# Get input temperature files
input_files_tmin = glob.glob(os.path.join(
    RESULTS_DIR, PROVINCE_NAME, "Temperature", "tmin", "**", "*.tif"), recursive=True)
input_files_tmax = glob.glob(os.path.join(
    RESULTS_DIR, PROVINCE_NAME, "Temperature", "tmax", "**", "*.tif"), recursive=True)

# Create output directory
new_dir = os.path.join(RESULTS_DIR, PROVINCE_NAME,
                       "_LS_Results", "Temperature")

# Process each percentile (and Monthly_Mean)
for p_idx, p_val in enumerate(T_PERC + ["Monthly_Mean"]):
    # Handle percentile naming and folder setup
    if p_val == "Monthly_Mean":
        per = "Monthly_Mean"
        folder_name = "Tbase_Tupper"
        start_name = "T_two_limits_W_"
    else:
        per = f"{str(p_val).replace('.', '')}perc"
        folder_name = f"Tmax_{per}"
        start_name = f"T_two_limits_Tmax_{per}_W_"

    # Create output directory
    os.makedirs(os.path.join(new_dir, folder_name), exist_ok=True)

    # Process each crop
    for crop_data in cropping_cal:
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
        weights = [crop_data[f"w_{i+1}"] for i in range(len(months))]

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

        # Filter tmax files for these months with specified percentile
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
                if month_abbr == os.path.basename(file)[:3]:
                    tmin_files_ordered.append(file)
                    break

            for file in tmax_files:
                if month_abbr == os.path.basename(file)[:3]:
                    tmax_files_ordered.append(file)
                    break

        # Load and calculate weighted tmin
        tmin_w_data = None
        for i, tmin_file in enumerate(tmin_files_ordered):
            with rasterio.open(tmin_file) as src:
                tmin_data = src.read(1) * weights[i]
                if tmin_w_data is None:
                    tmin_w_data = tmin_data
                    meta = src.meta.copy()
                else:
                    tmin_w_data += tmin_data

        # Load and calculate weighted tmax
        tmax_w_data = None
        for i, tmax_file in enumerate(tmax_files_ordered):
            with rasterio.open(tmax_file) as src:
                tmax_data = src.read(1) * weights[i]
                if tmax_w_data is None:
                    tmax_w_data = tmax_data
                    if tmin_w_data is None:  # If we didn't get tmin data
                        meta = src.meta.copy()
                else:
                    tmax_w_data += tmax_data

        # Create temperature suitability maps
        if tmin_w_data is not None and tmax_w_data is not None:
            # Check if tmin is higher than base temperature
            tmin_w_higher_tbase = (tmin_w_data > t_base).astype(np.uint8)

            # Check if tmax is lower than upper temperature
            tmax_w_lower_tupper = (tmax_w_data < t_upper).astype(np.uint8)

            # Combined suitability (both conditions must be met)
            two_limits_w = tmin_w_higher_tbase * tmax_w_lower_tupper

            # Save output
            season_label = f"{calendar.month_abbr[start_month]}-{calendar.month_abbr[end_month]}"
            output_path = os.path.join(
                new_dir, folder_name,
                f"{start_name}{crop_name}_{season_label}.tif")

            # Update metadata for output raster
            meta.update({
                'dtype': 'uint8',
                'count': 1
            })

            # Write output raster
            with rasterio.open(output_path, 'w', **meta) as dst:
                dst.write(two_limits_w, 1)

            print(
                f"Created temperature suitability map for {crop_name}, {season_label}, {per}")
        else:
            print(
                f"Missing temperature data for {crop_name}, {season_label}, {per}")

print("Temperature limits calculation complete!")

import os
import glob
import rasterio
import numpy as np
from pathlib import Path
import pandas as pd

"""
This script produces land suitability maps by:
1) Loading suitability layers from different parameters (NDVI, Temperature, Water, etc.)
2) Applying weightings to each parameter based on its importance
3) Combining all weighted layers to create a final land suitability map
4) Creating variations for different temperature/precipitation conditions
"""

# Define constants
RESULTS_DIR = "your_results_directory"  # Replace with actual path
PROVINCE_NAME = "your_province_name"  # Replace with actual province name
T_PERC = [0.05, 0.25]  # Temperature percentiles
T_PERC_NAMES = ["Cold", "Cool"]  # Names for temperature scenarios
P_PERC = [0.05, 0.25]  # Precipitation percentiles
P_PERC_NAMES = ["Dry", "Normal"]  # Names for precipitation scenarios

# Create a mock cropping calendar (replace with actual data loading)
cropping_cal = pd.DataFrame({
    'Crop': ['Wheat', 'Maize'],
    'Start_growing_season': [11, 5],
    'End_growing_season': [4, 9]
})

# Create a mock parameters table (replace with actual data loading)
params = pd.DataFrame({
    'Parameter': ['NDVI', 'Temperature', 'Water', 'Ksat', 'WCavail', 'Potasium', 'Phosphorus', 'Slope'],
    'Weight': [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
})

# Create output directory
new_dir = os.path.join(RESULTS_DIR, PROVINCE_NAME, "_LS_Results", "_Weighted")
os.makedirs(new_dir, exist_ok=True)

# Get input files for each parameter
input_files_ndvi = glob.glob(os.path.join(RESULTS_DIR, PROVINCE_NAME,
                                          "_LS_Results", "NDVI", "*.tif"))
input_files_temperature = glob.glob(os.path.join(RESULTS_DIR, PROVINCE_NAME,
                                                 "_LS_Results", "Temperature", "**", "*.tif"), recursive=True)
input_files_water = glob.glob(os.path.join(RESULTS_DIR, PROVINCE_NAME,
                                           "_LS_Results", "Water", "**", "*.tif"), recursive=True)
input_files_hhs = glob.glob(os.path.join(RESULTS_DIR, PROVINCE_NAME,
                                         "_LS_Results", "Soil_Hydraulic_Properties", "*.tif"))
input_files_snc = glob.glob(os.path.join(RESULTS_DIR, PROVINCE_NAME,
                                         "_LS_Results", "Soil_Nutrient_Content", "*.tif"))
input_files_slope = glob.glob(os.path.join(RESULTS_DIR, PROVINCE_NAME,
                                           "_LS_Results", "Elevation", "*.tif"))

# Function to find file matching a pattern


def find_file(file_list, pattern):
    matching_files = [f for f in file_list if pattern in f]
    return matching_files[0] if matching_files else None


# Process each crop
for _, crop_row in cropping_cal.iterrows():
    crop = crop_row['Crop']
    start_month_num = crop_row['Start_growing_season']
    end_month_num = crop_row['End_growing_season']

    # Convert month numbers to abbreviations
    month_abbrs = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
                   'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
    start_month = month_abbrs[start_month_num - 1]
    end_month = month_abbrs[end_month_num - 1]

    # Process each temperature scenario
    for t_idx, t_perc in enumerate(T_PERC + [None]):
        if t_perc is None:
            warm = "Average"
            per_t = "Tbase_Tupper"
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

            print(f"Processing {crop} for {dry}/{warm} conditions...")

            # Create output filename
            file_name = f"Land_Suitability_{dry}_{warm}_{crop}_{start_month}-{end_month}.tif"

            # Find and load each parameter file with appropriate weighting
            # NDVI
            ndvi_pattern = f"{start_month}-{end_month}"
            ndvi_file = find_file(input_files_ndvi, ndvi_pattern)
            if not ndvi_file:
                print(
                    f"NDVI file not found for {crop}, {start_month}-{end_month}")
                continue

            # Temperature
            temp_pattern = f"{crop}_{start_month}-{end_month}"
            temp_files = [
                f for f in input_files_temperature if per_t in f and temp_pattern in f]
            if not temp_files:
                print(f"Temperature file not found for {crop}, {per_t}")
                continue

            # Water
            water_pattern = f"{crop}_{start_month}-{end_month}_higher_"
            water_files = [
                f for f in input_files_water if per_p in f and water_pattern in f]
            if not water_files:
                print(f"Water file not found for {crop}, {per_p}")
                continue

            # Find soil and slope files
            ksat_file = find_file(input_files_hhs, "Ksat")
            wcavail_file = find_file(input_files_hhs, "WCavail")
            potassium_file = find_file(input_files_snc, "_K_")
            phosphorus_file = find_file(input_files_snc, "_P_")
            slope_file = find_file(input_files_slope, "lower")

            if not all([ksat_file, wcavail_file, potassium_file, phosphorus_file, slope_file]):
                print(f"One or more soil/slope files not found")
                continue

            # Get weights from parameters table
            ndvi_weight = params.loc[params['Parameter']
                                     == 'NDVI', 'Weight'].values[0]
            temp_weight = params.loc[params['Parameter']
                                     == 'Temperature', 'Weight'].values[0]
            water_weight = params.loc[params['Parameter']
                                      == 'Water', 'Weight'].values[0]
            ksat_weight = params.loc[params['Parameter']
                                     == 'Ksat', 'Weight'].values[0]
            wcavail_weight = params.loc[params['Parameter']
                                        == 'WCavail', 'Weight'].values[0]
            potassium_weight = params.loc[params['Parameter']
                                          == 'Potasium', 'Weight'].values[0]
            phosphorus_weight = params.loc[params['Parameter']
                                           == 'Phosphorus', 'Weight'].values[0]
            slope_weight = params.loc[params['Parameter']
                                      == 'Slope', 'Weight'].values[0]

            # Read and apply weights to each layer
            with rasterio.open(ndvi_file) as src:
                ndvi_data = src.read(1) * ndvi_weight
                output_meta = src.meta.copy()

            with rasterio.open(temp_files[0]) as src:
                temp_data = src.read(1) * temp_weight

            with rasterio.open(water_files[0]) as src:
                water_data = src.read(1) * water_weight

            with rasterio.open(ksat_file) as src:
                ksat_data = src.read(1) * ksat_weight

            with rasterio.open(wcavail_file) as src:
                wcavail_data = src.read(1) * wcavail_weight

            with rasterio.open(potassium_file) as src:
                potassium_data = src.read(1) * potassium_weight

            with rasterio.open(phosphorus_file) as src:
                phosphorus_data = src.read(1) * phosphorus_weight

            with rasterio.open(slope_file) as src:
                slope_data = src.read(1) * slope_weight

            # Sum all weighted layers
            ls_data = (ndvi_data + temp_data + water_data + ksat_data +
                       wcavail_data + potassium_data + phosphorus_data + slope_data)

            # Save the result
            output_path = os.path.join(new_dir, file_name)

            with rasterio.open(output_path, 'w', **output_meta) as dst:
                dst.write(ls_data, 1)

            print(f"Created land suitability map: {file_name}")

print("Land suitability map production complete!")

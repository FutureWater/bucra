import os
import glob
import rasterio
import numpy as np
import pandas as pd
import geopandas as gpd
from rasterstats import zonal_stats

"""
This script calculates average temperature and rainfall values per commune by:
1. Reading temperature and rainfall raster files
2. Extracting mean values for each commune using zonal statistics
3. Saving results to CSV files organized by variable type
"""

# Define directories and paths (constants)
DATA_DIR = "your_data_directory"  # Replace with actual path
RESULTS_DIR = "your_results_directory"  # Replace with actual path
PROVINCE_NAME = "your_province_name"  # Replace with actual province name

# Create output directories
new_dir_t = os.path.join(RESULTS_DIR, "_LS_Results", "Temperature")
new_dir_p = os.path.join(RESULTS_DIR, "_LS_Results", "Rainfall")
os.makedirs(new_dir_t, exist_ok=True)
os.makedirs(new_dir_p, exist_ok=True)

# Get temperature files
t_avg_files = glob.glob(os.path.join(
    RESULTS_DIR, PROVINCE_NAME, "Temperature", "tavg", "*.tif"))
t_max_files = glob.glob(os.path.join(
    RESULTS_DIR, PROVINCE_NAME, "Temperature", "tmax", "*.tif"))
t_min_files = glob.glob(os.path.join(
    RESULTS_DIR, PROVINCE_NAME, "Temperature", "tmin", "*.tif"))
t_files_all = t_avg_files + t_max_files + t_min_files

# Get rainfall files
p_avg_files = glob.glob(os.path.join(
    RESULTS_DIR, PROVINCE_NAME, "Rainfall", "005perc", "*.tif"))
p_max_files = glob.glob(os.path.join(
    RESULTS_DIR, PROVINCE_NAME, "Rainfall", "025perc", "*.tif"))
p_min_files = glob.glob(os.path.join(
    RESULTS_DIR, PROVINCE_NAME, "Rainfall", "Mean_Monthly", "*.tif"))
p_files_all = p_avg_files + p_max_files + p_min_files

# Load communes shapefile
# Note: In the R script, there's a filtering on communes_df$Province == province_name
# Here we assume communes shapefile already contains only the relevant communes
communes_shp_path = os.path.join(
    RESULTS_DIR, PROVINCE_NAME, "Shapefiles", f"{PROVINCE_NAME}.shp")
communes_shp = gpd.read_file(communes_shp_path)

# Process temperature files


def process_files(file_list, output_path, variable_type):
    """Process raster files and extract zonal statistics for each commune"""
    results = {}
    # Assuming NAME_3 field exists
    commune_names = communes_shp['NAME_3'].tolist()

    # Initialize empty dataframe with communes as rows
    for commune in commune_names:
        results[commune] = {}

    # Process each file
    for i, file_path in enumerate(file_list):
        print(f"Processing file {i+1} of {len(file_list)}")

        # Extract file information for naming
        if variable_type == "temperature":
            var_type = os.path.basename(
                os.path.dirname(os.path.dirname(file_path)))
            percentile = os.path.basename(os.path.dirname(file_path))
        else:  # rainfall
            var_type = "P"
            percentile = os.path.basename(os.path.dirname(file_path))

        # Get month from filename (first 3 characters)
        month = os.path.basename(file_path)[:3]
        column_name = f"{var_type}_{percentile}_{month}"

        # Calculate zonal statistics for each commune
        stats = zonal_stats(
            communes_shp,
            file_path,
            stats="mean",
            geojson_out=True
        )

        # Extract results
        for idx, stat in enumerate(stats):
            commune_name = communes_shp.iloc[idx]['NAME_3']
            mean_value = stat['properties']['mean']
            results[commune_name][column_name] = mean_value

    # Convert to dataframe and save
    df_results = pd.DataFrame.from_dict(results, orient='index')
    df_results.to_csv(output_path, index=True)


# Process temperature data
process_files(
    t_files_all,
    os.path.join(
        new_dir_t, f"Average_Temperature_per_Commune_{PROVINCE_NAME}.csv"),
    "temperature"
)

# Process rainfall data
process_files(
    p_files_all,
    os.path.join(
        new_dir_p, f"Average_Rainfall_per_Commune_{PROVINCE_NAME}.csv"),
    "rainfall"
)

print(f"Processing complete for {PROVINCE_NAME}")

import os
import glob
import rasterio
import numpy as np
import pandas as pd
import geopandas as gpd
from rasterio.mask import mask
import rioxarray

"""
This script calculates average temperature and rainfall values per commune by:
1. Reading temperature and rainfall raster files
2. Extracting mean values for each commune using zonal statistics
3. Saving results to CSV files organized by variable type
"""
print("Processing Average Rainfall and Temperature for communes...")
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


# Get temperature files
t_avg_files = glob.glob(os.path.join(
    RESULTS_DIR, "Temperature", "tavg", "Monthly_Mean", "*.tif"))
t_max_files = glob.glob(os.path.join(
    RESULTS_DIR, "Temperature", "tmax", "Monthly_Mean", "*.tif"))
t_min_files = glob.glob(os.path.join(
    RESULTS_DIR, "Temperature", "tmin", "Monthly_Mean", "*.tif"))
t_files_all = t_avg_files + t_max_files + t_min_files

# Get rainfall files
p_avg_files = glob.glob(os.path.join(
    RESULTS_DIR, "Rainfall", "005perc", "*.tif"))
p_max_files = glob.glob(os.path.join(
    RESULTS_DIR, "Rainfall", "025perc", "*.tif"))
p_min_files = glob.glob(os.path.join(
    RESULTS_DIR, "Rainfall", "Mean_Monthly", "*.tif"))
p_files_all = p_avg_files + p_max_files + p_min_files

var_dict = {"Temperature": t_files_all, "Rainfall": p_files_all}


# Load commune shapefile
commune_filepath = os.path.join(GIS_DIR, "Shapefiles", "AGO_adm3.shp")
commune_shp = gpd.read_file(commune_filepath)
commune_shp_sel = commune_shp[commune_shp["NAME_1"] == PROVINCE_NAME]
commune_shp_sel_reproj = commune_shp_sel.to_crs(LOCAL_PROJECTION)

# Process each variable
for var, file_list in var_dict.items():
    print(f"Processing {var} files...")
    # Create empty dictionary for results, to be converted to pandas df
    results = {}

    # Create subdirectory for new results
    new_results_subdir = os.path.join(RESULTS_DIR, "_LS_Results", var)
    os.makedirs(new_results_subdir, exist_ok=True)

    # Get list of commune_names in province.
    commune_names = commune_shp_sel['NAME_3'].tolist()

    # Initialize empty dictionary with communes as keys and a dictionary as value.
    for commune in commune_names:
        results[commune] = {}

    # Process each min, max and avg monthly file
    for i, file_path in enumerate(file_list):
        print(f"    Processing file {i+1} of {len(file_list)}")

        # Get variable type, percentile and month abbreviatons and create a column name out  of it.
        var_type = os.path.basename(
            os.path.dirname(os.path.dirname(file_path)))
        # gets percentile or monthly mean as string
        percentile = os.path.basename(os.path.dirname(file_path))
        month = os.path.basename(file_path)[:3]
        column_name = f"{var_type}_{percentile}_{month}"

        # Calculate mean value of raster
        with rioxarray.open_rasterio(file_path) as raster:
            for idx, row in commune_shp_sel_reproj.iterrows():
                commune_name = row['NAME_3']

                # Clip raster to commune geometry and calculate mean
                clipped = raster.rio.clip([row.geometry], drop=False)
                mean_value = float(clipped.mean())
                results[commune_name][column_name] = mean_value

            # Convert to dataframe and save
            output_path = os.path.join(
                new_results_subdir, f"Average_{var}_per_Commune_of_{PROVINCE_NAME}.csv")
            df_results = pd.DataFrame.from_dict(results, orient='index')
            df_results.to_csv(output_path, index=True)


print(f"Processing complete for {PROVINCE_NAME}")

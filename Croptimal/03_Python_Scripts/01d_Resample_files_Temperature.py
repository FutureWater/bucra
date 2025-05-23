import os
import glob
import calendar
import rasterio
import numpy as np
import xarray as xr
import geopandas as gpd
from rasterio.warp import reproject, Resampling
from rasterio.mask import mask
from rasterio.transform import from_origin
import matplotlib.pyplot as plt


"""
This script processes temperature data files by:
1. Reading monthly temperature data from NetCDF files
2. Calculating monthly means across years
3. Cropping to a buffered province boundary
4. Resampling to match a reference DEM
5. Applying lapse rate correction based on elevation
6. Saving the processed data
"""

####################################################################################################
########################## Define directories and file paths #######################################
####################################################################################################
# Define current and parent working directories
current_wd = os.getcwd()
parent_wd = os.path.dirname(current_wd)
base_wd = os.path.dirname(os.path.dirname(parent_wd))
angola_wd = "/Users/thomasfuturewater/FutureWater Dropbox/Team/Projects/Completed/2019/2019019_G4AW_MavoDiami_Angola/Data/2019019_MavoDiami_LV/2019019_MavoDiami"

# Gets province from subprocess in 000_Run_All.py
PROVINCE_NAME = os.environ.get("PROVINCE")
if not os.environ.get("PROVINCE"):
    PROVINCE_NAME = "Sharkia"                              # dummy variable for testing.

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

# Define temperature variables
T_VARS = ["Tavg", "Tmin", "Tmax_P75", "Tmax_P95",
          "Tmax_Mean"]  # Temperature variables to process
T_LAPSE_RATE = -0.0065  # Temperature lapse rate (°C/m)

# Month abbreviations for output file naming
month_abbrs = [calendar.month_abbr[i] for i in range(1, 13)]

# Define input and output paths
INPUT_FILES = glob.glob(os.path.join(DATA_DIR, "Temperature", "*.tif"))
T_VARS = [os.path.basename(file)[:-12] for file in INPUT_FILES]
RESULTS_TEMPERATURE_DIR = os.path.join(RESULTS_DIR, "Temperature")
os.makedirs(RESULTS_TEMPERATURE_DIR, exist_ok=True)


####################################################################################################
####################### Loading DEM and province shapefile #########################################
# Import DEM for reference extent and resolution
DEM_PATH = os.path.join(RESULTS_DIR, "DEM", f"DEM_{PROVINCE_NAME}_{RES}m.tif")

with rasterio.open(DEM_PATH) as dem_src:
    DEM_DATA = dem_src.read(1)
    DEM_PROFILE = dem_src.profile.copy()
    DEM_BOUNDS = dem_src.bounds
    DEM_TRANSFORM = dem_src.transform
    DEM_CRS = dem_src.crs
    DEM_HEIGHT = dem_src.height
    DEM_WIDTH = dem_src.width

# Load province shapefile and set to crs of DEM
provinces_filepath = os.path.join(GIS_DIR, "Nile_delta_bnd_adm1.shp")
provinces_shp = gpd.read_file(provinces_filepath)
province_shp_sel = provinces_shp[provinces_shp["ADM1_EN"] == PROVINCE_NAME]
province_shp_reproj = province_shp_sel.to_crs(DEM_CRS)


####################################################################################################
########################## Check Temperature files and add lapse reate #############################
####################################################################################################
for i, (file_name, var) in enumerate(zip(INPUT_FILES, T_VARS)):
    # Load temperature file
    print(f"Processing {var} temperature data...")

    # Create output directory
    var_results_subdir = os.path.join(RESULTS_TEMPERATURE_DIR, var)
    os.makedirs(var_results_subdir, exist_ok=True)

    # Process each month (band) in the rainfall geotiff
    with rasterio.open(file_name) as src:
        num_bands = src.count

        for i in range(1, num_bands + 1):
            # Read the data and profile
            temperature_data = src.read(i)
            temperature_profile = src.profile

            # Get the month abbreviation
            month = month_abbrs[i-1]
            print(f"  Processing {var} for {month}...")

            # Create empty destination array
            destination_array = np.full(
                (DEM_HEIGHT, DEM_WIDTH), NO_DATA_VALUE, dtype=temperature_data.dtype)

            # Reproject and write
            reproject(
                source=temperature_data,
                destination=destination_array,
                src_transform=src.transform,
                src_crs=src.crs,
                dst_transform=DEM_TRANSFORM,
                dst_crs=DEM_CRS,
                resampling=Resampling.bilinear,
                src_nodata=src.nodata,
                dst_nodata=NO_DATA_VALUE
            )

            # Create a memory file with the reprojected data for masking
            with rasterio.MemoryFile() as memfile:
                with memfile.open(**DEM_PROFILE) as temp_dst:
                    temp_dst.write(destination_array, 1)

                    # Perform the masking operation
                    masked_data, masked_transform = mask(
                        temp_dst,
                        province_shp_reproj.geometry,
                        crop=True,
                        nodata=NO_DATA_VALUE
                    )

                    # Update profile for the masked result
                    masked_profile = temp_dst.profile.copy()
                    masked_profile.update({
                        'height': masked_data.shape[1],
                        'width': masked_data.shape[2],
                        'transform': masked_transform,
                    })

            # Apply lapse rate correction based on elevation difference
            no_data_mask = masked_data != NO_DATA_VALUE
            final_data = np.where(
                no_data_mask, destination_array + DEM_DATA * T_LAPSE_RATE, NO_DATA_VALUE)

            # Save final temperature raster
            output_path = os.path.join(
                var_results_subdir, f"{var}_{month}.tif")
            with rasterio.open(output_path, 'w', **masked_profile) as dst:
                dst.write(final_data)

    print(f"Completed processing {var} temperature data")

print("Temperature data processing complete!")

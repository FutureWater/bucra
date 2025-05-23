import os
import glob
import calendar
import rasterio
import numpy as np
from rasterio.warp import reproject, Resampling, calculate_default_transform
from rasterio.mask import mask
from rasterio.transform import from_origin
from rasterio.plot import show
# import rioxarray
import geopandas as gpd
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt

"""
This script processes reference evapotranspiration (ET) data by:
1. Reading monthly ET data from NetCDF files
2. Calculating monthly means
3. Cropping to a province boundary
4. Resampling to match a reference DEM
5. Applying monthly day multipliers
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
# PROVINCE_NAME = "Sharkia"                         # dummy variable for testing.

# Define other folders
DATA_DIR = os.path.join(parent_wd, "01_Data")
GIS_DIR = os.path.join(base_wd, "GIS")     # Directory with shapefiles
RESULTS_DIR = os.path.join(parent_wd, "04_Results", PROVINCE_NAME)
os.makedirs(RESULTS_DIR, exist_ok=True)
TEMP_DIR = os.path.join(parent_wd, "05_Temp")

# Define constants
RES = 250  # Set resolution in meters
NO_DATA_VALUE = -9999.0
DEM_PATH = os.path.join(
    RESULTS_DIR, "DEM", f"DEM_{PROVINCE_NAME}_{RES}m.tif")  # Path to DEM file

# Define input and output paths
INPUT_FILE = glob.glob(os.path.join(DATA_DIR, "RET", "*.tif"))
REF_ET_RESULTS_DIR = os.path.join(RESULTS_DIR, "Ref_ET")
os.makedirs(REF_ET_RESULTS_DIR, exist_ok=True)

####################################################################################################
####################### Loading DEM and province shapefile #########################################
####################################################################################################
# Import DEM for reference extent and resolution
DEM_PATH = os.path.join(RESULTS_DIR, "DEM", f"DEM_{PROVINCE_NAME}_{RES}m.tif")

with rasterio.open(DEM_PATH) as dem_src:
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
################################# Resample ET raster to match DEM ##################################
####################################################################################################
# Multiplier for days in each month
MULTIPLIER = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])

MONTH_ABBRS = [calendar.month_abbr[i] for i in range(1, 13)]

# Process each month (band) in the rainfall geotiff
with rasterio.open(INPUT_FILE[0]) as src:
    print("Opening RET files and calculating monthly mean RET")

    num_bands = src.count

    for i in range(1, num_bands + 1):
        # Read the data and profile
        RET_data = src.read(i)
        RET_profile = src.profile

        # Get the month abbreviation
        month = MONTH_ABBRS[i-1]
        print(f"  Processing RET for {month}...")

        # Create empty destination array
        destination_array = np.full(
            (DEM_HEIGHT, DEM_WIDTH), NO_DATA_VALUE, dtype=RET_data.dtype)

        # Reproject and write
        reproject(
            source=RET_data,
            destination=destination_array,
            src_transform=src.transform,
            src_crs=src.crs,
            dst_transform=DEM_TRANSFORM,
            dst_crs=DEM_CRS,
            resampling=Resampling.bilinear,
            src_nodata=NO_DATA_VALUE,
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

        # Write the final masked result to file
        OUTPUT_PATH = os.path.join(REF_ET_RESULTS_DIR, f"RET_{month}.tif")
        with rasterio.open(OUTPUT_PATH, 'w', **masked_profile) as dst:
            dst.write(masked_data)


# Show final ET raster
# with rasterio.open(output_path) as src:
#     show(src, title="Final ET raster with rasterio")

print(f"Processing complete for {PROVINCE_NAME}")

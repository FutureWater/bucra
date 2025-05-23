import os
import glob
import rasterio
import numpy as np
from rasterio.warp import reproject, Resampling, calculate_default_transform
from rasterio.mask import mask
import geopandas as gpd
from rasterio.merge import merge
import matplotlib.pyplot as plt

"""
This script processes Soil Nutrient Content (SNC) data by:
1. Reading soil nutrient content TIF files
2. Cropping to a buffered province boundary
3. Resampling to match a reference DEM
4. Saving the processed data
"""

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
    PROVINCE_NAME = "Sharkia"                       # dummy variable for testing.

# Define other folders
DATA_DIR = os.path.join(parent_wd, "01_Data")
GIS_DIR = os.path.join(base_wd, "GIS")     # Directory with shapefiles
RESULTS_DIR = os.path.join(parent_wd, "04_Results", PROVINCE_NAME)
os.makedirs(RESULTS_DIR, exist_ok=True)
TEMP_DIR = os.path.join(parent_wd, "05_Temp")
PROVINCES_FILEPATH = os.path.join(GIS_DIR, "Shapefiles", "AGO_adm1.shp")

# Define constants
RES = 250  # Set resolution in meters
LOCAL_PROJ = "EPSG:32636"
NO_DATA_VALUE = -9999.0

# Create output directory for processed SNC data
SNC_RESULTS_DIR = os.path.join(RESULTS_DIR, "Soil_Nutrient_Content")
os.makedirs(SNC_RESULTS_DIR, exist_ok=True)

# Get input SNC files
INPUT_FILES = glob.glob(os.path.join(
    DATA_DIR, "Soil_Nutrient_Content", "*.tif"))

if not INPUT_FILES:
    raise FileNotFoundError(
        f"No SNC files found in {os.path.join(DATA_DIR, 'Soil_Nutrient_Content')}")

print(f"Found {len(INPUT_FILES)} SNC files")

# Get file names without extension and replace "af" with province_name
VAR_NAMES = [os.path.basename(file) for file in INPUT_FILES]
VAR_NAMES_2 = [name.replace(
    "0to20cm", f"{PROVINCE_NAME}") for name in VAR_NAMES]

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
################################# Process SNC files ################################################
####################################################################################################
# Process each SNC file
for i, (input_file, output_name) in enumerate(zip(INPUT_FILES, VAR_NAMES_2)):
    print(
        f"  Processing {os.path.basename(input_file)} ({i+1}/{len(INPUT_FILES)})")

    # Read input file
    with rasterio.open(input_file) as src:
        # Calculate window/subset that covers the mask extent
        province_shp_window_crs = province_shp_sel.to_crs(src.crs)
        window = src.window(*province_shp_window_crs.total_bounds)
        temp_profile = src.profile

        # Read only that subset of the raster
        subset_data = src.read(1, window=window)
        subset_transform = src.window_transform(window)

        # Resample to match DEM resolution and extent
        print("  Resampling to match DEM...")
        destination_array = np.full(
            (DEM_HEIGHT, DEM_WIDTH), NO_DATA_VALUE, dtype=np.float32)

        reproject(
            source=subset_data,
            destination=destination_array,
            src_transform=subset_transform,
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

        # Write the final masked result to file
        OUTPUT_PATH = os.path.join(SNC_RESULTS_DIR, 'Extractable_' + output_name)
        with rasterio.open(OUTPUT_PATH, 'w', **masked_profile) as dst:
            dst.write(masked_data.astype(rasterio.float32))

print("SNC resampling complete!")

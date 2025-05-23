import os
import glob
import rasterio
import numpy as np
from rasterio.warp import reproject, Resampling, calculate_default_transform
from rasterio.mask import mask
from rasterio.windows import Window
from rasterio.transform import array_bounds
import geopandas as gpd
import matplotlib.pyplot as plt

"""
This script processes Soil Hydraulic Propertie data by:
1. Reading soil hydraulic properites ontent TIF files
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

# Gets province from subprocess in 000_Run_All.py
PROVINCE_NAME = os.environ.get("PROVINCE")
if not os.environ.get("PROVINCE"):
    PROVINCE_NAME = "Sharkia"                         # dummy variable for testing.

# Define other folders
DATA_DIR = os.path.join(parent_wd, "01_Data")
GIS_DIR = os.path.join(base_wd, "GIS")     # Directory with shapefiles
RESULTS_DIR = os.path.join(parent_wd, "04_Results", PROVINCE_NAME)
os.makedirs(RESULTS_DIR, exist_ok=True)
TEMP_DIR = os.path.join(parent_wd, "05_Temp")

# Define constants
RES = 250  # Set resolution in meters
LOCAL_PROJ = "EPSG:32636"
NO_DATA_VALUE = -9999.0
VAR_NAMES = ['WCavail', 'Ksat']

# Create output directory for processed SHP data
SHP_RESULTS_DIR = os.path.join(RESULTS_DIR, "Soil_Hydraulic_Properties")
os.makedirs(SHP_RESULTS_DIR, exist_ok=True)

# Get input SHP files
# INPUT_FILES = [file for var in VAR_NAMES for file in glob.glob(os.path.join(
#     '/Users/thomasfuturewater/FutureWater Dropbox/Team/Data/Global/Soil/HiHydroSoil_250m/Top_Subsoil', var, '*.tif'))]
# if not INPUT_FILES:
#     raise FileNotFoundError(
#         f"No SHP files found in 'FutureWater Dropbox/Team/Data/Global/Soil/HiHydroSoil_250m'")

INPUT_FILES = glob.glob(os.path.join(
    DATA_DIR, "Soil_Nutrient_Content", "*.tif"))
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
for i, (input_file, var) in enumerate(zip(INPUT_FILES, VAR_NAMES)):
    print(
        f"  Processing {os.path.basename(input_file)} ({i+1}/{len(INPUT_FILES)})")

    # Read input file
    with rasterio.open(input_file) as src:
        # Calculate window/subset that covers the mask extent
        province_shp_window_crs = province_shp_sel.to_crs(src.crs).buffer(10000)
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
            resampling=Resampling.nearest,
            src_nodata=src.nodata,
            dst_nodata=NO_DATA_VALUE
        )
        temp_output = os.path.join(TEMP_DIR, f"temp_{var}_reprojected_nearest.tif")
        with rasterio.open(temp_output, 'w', **DEM_PROFILE) as temp_dst:
            temp_dst.write(destination_array, 1)

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
        top_sub = os.path.basename(input_file)[-11:-4]
        output_name = var + '_' + top_sub + '_' + PROVINCE_NAME + ".tif"
        OUTPUT_PATH = os.path.join(SHP_RESULTS_DIR, output_name)
        with rasterio.open(OUTPUT_PATH, 'w', **masked_profile) as dst:
            dst.write(masked_data.astype(rasterio.float32))

print("SHP resampling complete!")

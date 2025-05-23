import os
import glob
import calendar
import rasterio
import numpy as np
from rasterio.warp import reproject, Resampling
from rasterio.windows import from_bounds
from rasterio.mask import mask
import geopandas as gpd
import matplotlib.pyplot as plt
from rasterio.plot import show


"""
This script processes NDVI files by:
1. Defining directories and file paths
2. Creating output directories
3. Loading input NDVI files
4. Importing DEM for reference extent and resolution
5. Processing each NDVI file individually to avoid memory issues:
   a. Reprojecting/resampling data to match DEM
   b. Cropping raster to extent of the DEM
   c. Saving the resampled and cropped raster
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
    PROVINCE_NAME = "Sharkia"                           # dummy variable for testing.

# Define other folders
DATA_DIR = os.path.join(parent_wd, "01_Data")
GIS_DIR = os.path.join(base_wd, "GIS")     # Directory with shapefiles
RESULTS_DIR = os.path.join(parent_wd, "04_Results", PROVINCE_NAME)
os.makedirs(RESULTS_DIR, exist_ok=True)

NDVI_MM_DIR = os.path.join(RESULTS_DIR, "NDVI", "Mean_Monthly")
os.makedirs(NDVI_MM_DIR, exist_ok=True)

# Set constants
NO_DATA_VALUE = -9999.0
RES = 250  # Set resolution
LOCAL_PROJ = "EPSG:32636"  # Local projection

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
################################# Resample NDVI files ##############################################
####################################################################################################
# Get input files
INPUT_FILE = os.path.join(DATA_DIR, "NDVI", "NDVI_mean_monthly_stack.tif")

# Process each NDVI file individually to avoid memory issues
month_abbrs = [calendar.month_abbr[i] for i in range(1, 13)]

# Iterate over al the months (bands in the NDVI geotiff)
print("Resampling NDVI files...")
with rasterio.open(INPUT_FILE) as src:
    num_bands = src.count
    for i in range(1, num_bands + 1):
        temp_data = src.read(i)
        temp_profile = src.profile

        month = month_abbrs[i-1]
        print(f"  Processing NDVI for {month}...")

        # Create empty destination array
        destination_array = np.full(
            (DEM_HEIGHT, DEM_WIDTH), NO_DATA_VALUE, dtype=temp_data.dtype)

        # Reproject and write resampled raster
        reproject(
            source=temp_data,
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

        output_path = os.path.join(
            NDVI_MM_DIR, f"NDVI_{month}.tif")
        with rasterio.open(output_path, 'w', **masked_profile) as dst:
            dst.write(masked_data.astype(rasterio.float32))

print("NDVI resampling complete!")

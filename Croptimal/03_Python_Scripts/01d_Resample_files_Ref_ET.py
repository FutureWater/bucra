import os
import glob
import rasterio
import numpy as np
from rasterio.warp import reproject, Resampling
from rasterio.mask import mask
import geopandas as gpd
import netCDF4 as nc
import pandas as pd

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
# Define directories and file paths
current_wd = os.getcwd()
parent_wd = os.path.dirname(current_wd)
angola_wd = "/Users/thomasfuturewater/FutureWater Dropbox/Team/Projects/Completed/2019/2019019_G4AW_MavoDiami_Angola/Data/2019019_MavoDiami_LV/2019019_MavoDiami"

# Gets province from subprocess in 000_Run_All.py
PROVINCE_NAME = os.environ.get("PROVINCE")
PROVINCE_NAME = "Zaire"                         # dummy variable for testing.

# Set directories
DATA_DIR = os.path.join(angola_wd, "01_Data")
GIS_DIR = os.path.join(angola_wd, "02_GIS")     # Directory with shapefiles
RESULTS_DIR = os.path.join(parent_wd, "04_Results", PROVINCE_NAME)
TEMP_DIR = os.path.join(parent_wd, "05_Temp")
RES = 250  # Set resolution in meters
LOCAL_PROJ = "EPSG:32733"

# Define input and output paths
INPUT_FILES = glob.glob(os.path.join(DATA_DIR, "Ref_ET/*.ymonmean.nc"))
REF_ET_RESULTS_DIR = os.path.join(RESULTS_DIR, "Ref_ET/Mean_Monthly/")
os.makedirs(REF_ET_RESULTS_DIR, exist_ok=True)

####################################################################################################
################################# Read and process NetCDF data ######################################
####################################################################################################
# Multiplier for days in each month
MULTIPLIER = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])

# Open NetCDF file to get time information
with nc.Dataset(INPUT_FILES[0], "r") as dataset:
    time_nc = dataset.variables["time"][:]
    st_date = nc.num2date(time_nc, dataset.variables["time"].units)
    index_names = [date.strftime("%m") for date in st_date]

# Read and process NetCDF data
brick_data = []
for file in INPUT_FILES:
    with nc.Dataset(file, "r") as dataset:
        data = dataset.variables["ET_ref"][:]
        brick_data.append(data)

data_mean = np.mean(np.stack(brick_data), axis=0)

# Crop raster to extent of the province
print(f"Cropping ET raster for {PROVINCE_NAME}")
PROVINCES_FILEPATH = os.path.join(GIS_DIR, "Shapefiles", "AGO_adm1.shp")
provinces_shp = gpd.read_file(PROVINCES_FILEPATH)
PROVINCE_SHP_SEL = provinces_shp[provinces_shp["NAME_1"] == PROVINCE_NAME]
province_shp_proj = PROVINCE_SHP_SEL.to_crs(LOCAL_PROJ)

# Crop raster to province boundary and get transform, extent, and metadata.
with rasterio.open(INPUT_FILES[0]) as src:
    out_image, out_transform = mask(src, province_shp_proj.geometry, crop=True)
    out_meta = src.meta.copy()
    out_meta.update({
        "driver": "GTiff",
        "height": out_image.shape[1],
        "width": out_image.shape[2],
        "transform": out_transform
    })


####################################################################################################
################################# Resample ET raster to match DEM ##################################
####################################################################################################
# Read resampled DEM raster and get profile
dem_diff_path = os.path.join(
    RESULTS_DIR, "DEM", f"DEM_{PROVINCE_NAME}_{RES}m.tif")  # Used to be: "DEM_{PROVINCE_NAME}_{RES}m_diff.tif"
with rasterio.open(dem_diff_path) as dem_src:
    dem_profile = dem_src.profile

# Resample ET raster to match DEM
print(f"Resampling ET raster to match DEM for {PROVINCE_NAME}")
resampled_raster = np.empty_like(data_mean)
reproject(
    source=data_mean,
    destination=resampled_raster,
    src_transform=out_transform,
    src_crs=LOCAL_PROJ,
    dst_transform=dem_src.transform,
    dst_crs=LOCAL_PROJ,
    resampling=Resampling.bilinear
)

# Apply multiplier and write final ET rasters
print(f"Writing final ET raster: {REF_ET_RESULTS_DIR, PROVINCE_NAME}")
MONTH_ABBRS = [calendar.month_abbr[i] for i in range(1, 13)]

for i, month in enumerate(MONTH_ABBRS):
    output_path = os.path.join(REF_ET_RESULTS_DIR, f"{month}.tif")
    with rasterio.open(output_path, "w", **dem_profile) as dst:
        dst.write(resampled_raster * MULTIPLIER[i], 1)

print(f"Processing complete for {PROVINCE_NAME}")

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
# Define directories and file paths (constants in UPPERCASE)
DATA_DIR = "your_data_directory"
RESULTS_DIR = "your_results_directory"
PROVINCE_NAME = "your_province_name"
RES = 250  # Set resolution in meters
LOCAL_PROJ = "EPSG:32733"

# Define input and output paths
INPUT_FILES = glob.glob(os.path.join(DATA_DIR, "Ref_ET/*.ymonmean.nc"))
NEW_DIR = os.path.join(RESULTS_DIR, PROVINCE_NAME, "Ref_ET/Mean_Monthly/")
os.makedirs(NEW_DIR, exist_ok=True)

####################################################################################################
################################# Read and process NetCDF data ######################################
####################################################################################################
# Multiplier for days in each month
MULTIPLIER = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])

# Open NetCDF file
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
province_shp_path = os.path.join(
    RESULTS_DIR, PROVINCE_NAME, "Shapefiles", f"{PROVINCE_NAME}.shp")
province_shp = gpd.read_file(province_shp_path)
province_shp_proj = province_shp.to_crs(LOCAL_PROJ)

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
# Read resampled DEM raster
dem_diff_path = os.path.join(
    RESULTS_DIR, PROVINCE_NAME, "DEM", f"DEM_{PROVINCE_NAME}_{RES}m_diff.tif")
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
print(f"Writing final ET rasters for {PROVINCE_NAME}")
MONTHS = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
          "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]

for i, month in enumerate(MONTHS):
    output_path = os.path.join(NEW_DIR, f"{month}.tif")
    with rasterio.open(output_path, "w", **dem_profile) as dst:
        dst.write(resampled_raster * MULTIPLIER[i], 1)

print(f"Processing complete for {PROVINCE_NAME}")

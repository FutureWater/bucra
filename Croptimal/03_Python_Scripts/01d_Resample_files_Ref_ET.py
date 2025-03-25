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
import netCDF4 as nc
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
SRC_CRS = "EPSG:4326"  # CRS from Netcdf files
NO_DATA_VALUE = -9999.0

# Define input and output paths
INPUT_FILES = glob.glob(os.path.join(DATA_DIR, "Ref_ET/*ymonmean.nc"))
REF_ET_RESULTS_DIR = os.path.join(RESULTS_DIR, "Ref_ET/Mean_Monthly/")
os.makedirs(REF_ET_RESULTS_DIR, exist_ok=True)

####################################################################################################
################################# Read and process NetCDF data ######################################
####################################################################################################
# Multiplier for days in each month
MULTIPLIER = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])

# Get province boundary
provinces_filepath = os.path.join(GIS_DIR, "Shapefiles", "AGO_adm1.shp")
provinces_shp = gpd.read_file(provinces_filepath)
province_shp_sel = provinces_shp[provinces_shp["NAME_1"] == PROVINCE_NAME]
province_shp_proj = province_shp_sel.to_crs(SRC_CRS)

# Open NetCDF file and calculate monthly means
ds = xr.open_dataset(INPUT_FILES[0])       # Open Dataset
time_coord = ds['time']
ref_et = ds['PMMT']
lon = ds['lon']
lat = ds['lat']
monthly_mean = ref_et.groupby(time_coord.dt.month).mean(dim='time')
ds.close()                                          # Close datasets


# Read resampled DEM raster and get profile
DEM_PATH = os.path.join(
    RESULTS_DIR, "DEM", f"DEM_{PROVINCE_NAME}_{RES}m.tif")  # Used to be: "DEM_{PROVINCE_NAME}_{RES}m_diff.tif"

with rasterio.open(DEM_PATH) as dem_src:
    DEM_PROFILE = dem_src.profile
    DEM_CRS = dem_src.crs
    DEM_TRANSFORM = dem_src.transform
    DEM_SHAPE = (dem_src.height, dem_src.width)

####################################################################################################
################################# Resample ET raster to match DEM ##################################
####################################################################################################
print(f"Writing final ET raster: {REF_ET_RESULTS_DIR, PROVINCE_NAME}")
MONTH_ABBRS = [calendar.month_abbr[i] for i in range(1, 13)]

for i, month in enumerate(MONTH_ABBRS[:1]):
    # Select data for month from multi-month dataset
    month_data = monthly_mean.isel(month=i).values * MULTIPLIER[i]

    """
    Transform = (most western point, most high point, distance of moving 1 pixel to the right, distance of moving one pixel down the raster)
    Geographic coordinate systems are descending in latitude (from north to south), 
    as you move down the raster. Therefore, in the transform it needs to be positive, 
    as latitude increases as you move down the raster. Most netcdf files are descending in order, so positive lat value.
    """
    # Save mean monthly RET as raster in temp folder, with WGS84 CRS and transform from dataset.
    src_transform = from_origin(lon.min(), lat.max(),
                                abs(lon[1]-lon[0]), abs(lat[1]-lat[0]))

    temp_path = os.path.join(TEMP_DIR, f"temp_{month}.tif")

    with rasterio.open(
        temp_path, 'w',
        driver='GTiff',
        height=monthly_mean.lat.size,
        width=monthly_mean.lon.size,
        count=1,
        dtype=month_data.dtype,
        crs=SRC_CRS,
        transform=src_transform
    ) as tmp:
        tmp.write(month_data, 1)

    # Open temporary raster and reproject and resample it to DEM profile.
    with rasterio.open(temp_path) as src:
        # Clip to province boundary
        out_image, out_transform = mask(
            src, province_shp_proj.geometry, crop=True, nodata=NO_DATA_VALUE)

        # Create destination transform
        dst_transform, dst_width, dst_height = calculate_default_transform(
            src.crs, DEM_CRS, src.width, src.height,
            *rasterio.transform.array_bounds(out_image.shape[1], out_image.shape[2], out_transform),
            resolution=RES
        )

        # Set profile for output raster
        dst_profile = src.profile.copy()
        dst_profile.update({
            'crs': DEM_CRS,
            'transform': dst_transform,
            'width': dst_width,
            'height': dst_height,
            'nodata': NO_DATA_VALUE  # Set a nodata value
        })

        # Set output path
        output_path = os.path.join(REF_ET_RESULTS_DIR, f"RET_{month}.tif")

        # Resample temp RET raster to DEM resolution.
        destination_array = np.full(
            (dst_height, dst_width), NO_DATA_VALUE, dtype=out_image.dtype)

        with rasterio.open(output_path, 'w', **dst_profile) as dst:
            # Reproject and write
            reproject(
                source=out_image[0],
                destination=destination_array,
                src_transform=out_transform,
                src_crs=SRC_CRS,
                dst_transform=dst_transform,
                dst_crs=DEM_CRS,
                resampling=Resampling.bilinear,
                src_nodata=NO_DATA_VALUE,
                dst_nodata=NO_DATA_VALUE
            )
            # Write the reprojected data to the output file
            dst.write(destination_array.astype(rasterio.float32), 1)
    os.remove(temp_path)

# Show final ET raster
# with rasterio.open(output_path) as src:
#     show(src, title="Final ET raster with rasterio")

print(f"Processing complete for {PROVINCE_NAME}")

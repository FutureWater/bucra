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
angola_wd = "/Users/thomasfuturewater/FutureWater Dropbox/Team/Projects/Completed/2019/2019019_G4AW_MavoDiami_Angola/Data/2019019_MavoDiami_LV/2019019_MavoDiami"


# Gets province from subprocess in 000_Run_All.py
PROVINCE_NAME = os.environ.get("PROVINCE")
# PROVINCE_NAME = "Zaire"                         # dummy variable for testing.

# Define other folders
DATA_DIR = os.path.join(angola_wd, "01_Data")
GIS_DIR = os.path.join(angola_wd, "02_GIS")     # Directory with shapefiles
RESULTS_DIR = os.path.join(parent_wd, "04_Results", PROVINCE_NAME)
TEMP_DIR = os.path.join(parent_wd, "05_Temp")

# Temperature variables and parameters
RES = 250  # Resolution in meters
NO_DATA_VALUE = -9999.0  # No data value
LOCAL_PROJ = "EPSG:32733"  # Local projection
SRC_CRS = "EPSG:4326"  # Projection from NetCDF files: WGS84 projection
T_VARS = ["tavg", "tmin", "tmax"]  # Temperature variables to process
T_LAPSE_RATE = -0.0065  # Temperature lapse rate (°C/m)

# Month abbreviations for output file naming
month_abbrs = [calendar.month_abbr[i] for i in range(1, 13)]


####################################################################################################
####################### Loading DEM and province shapefile #########################################
# print("Loading reference DEM...")
# Load reference DEM
dem_path = os.path.join(RESULTS_DIR, "DEM",
                        # Used to be "DEM_{PROVINCE_NAME}_{RES}m_diff.tif"
                        f"DEM_{PROVINCE_NAME}_{RES}m.tif")

with rasterio.open(dem_path) as dem_src:
    DEM_DATA = dem_src.read(1)
    DEM_META = dem_src.meta.copy()
    DEM_TRANSFORM = dem_src.transform
    DEM_CRS = dem_src.crs
    DEM_HEIGHT = dem_src.height
    DEM_WIDTH = dem_src.width
    DEM_PROFILE = dem_src.profile
# Load province shapefile and set to crs of DEM
provinces_filepath = os.path.join(GIS_DIR, "Shapefiles", "AGO_adm1.shp")
provinces_shp = gpd.read_file(provinces_filepath)
province_shp_sel = provinces_shp[provinces_shp["NAME_1"] == PROVINCE_NAME]
province_shp_proj = province_shp_sel.to_crs(DEM_CRS)

####################################################################################################
################################# Process temperature variables ###################################
####################################################################################################
# Loop through each temperature variable (tavg, tmin, tmax)
for var in T_VARS:
    print(f"Processing {var} temperature data...")

    # Get input files for this temperature variable. Files are per year, with monthly temp data.
    input_files = glob.glob(os.path.join(DATA_DIR, "Temperature", var, "*.nc"))

    if not input_files:
        print(f"No {var} files found. Skipping.")
        continue

    # Create output directory
    var_results_subdir = os.path.join(RESULTS_DIR,
                                      "Temperature", var, "Monthly_Mean")
    os.makedirs(var_results_subdir, exist_ok=True)

    # Read and process NetCDF data.
    datasets = []
    for file in input_files:
        with xr.open_dataset(file) as ds:
            datasets.append(ds)

    # Combine all datasets
    combined_data = xr.concat(datasets, dim="time")

    # Group by month and calculate mean for each month
    monthly_means = combined_data.groupby('time.month').mean(dim='time')

    # Process each month
    for month_idx in range(1, 13):
        month_name = month_abbrs[month_idx - 1]
        month_data = monthly_means['t2m'].sel(month=month_idx).values

        # Create a temporary raster file for cropping and reprojection
        temp_path = os.path.join(TEMP_DIR, f"temp_{var}_{month_name}.tif")

        with xr.open_dataset(input_files[0]) as ds:
            lon = ds['longitude']
            lat = ds['latitude']
            ds_transform = from_origin(lon.min(), lat.max(),
                                       abs(lon[1]-lon[0]), abs(lat[1]-lat[0]))

        # Create a profile for the raster
        temp_profile = {
            'driver': 'GTiff',
            'height': month_data.shape[0],
            'width': month_data.shape[1],
            'count': 1,
            'dtype': month_data.dtype,
            'crs': SRC_CRS,
            # 'transform': ds.rio.transform()
            'transform': ds_transform
        }

        # Write each month to temporary file
        with rasterio.open(temp_path, 'w', **temp_profile) as dst:
            dst.write(month_data.astype(rasterio.float32), 1)

        # Crop month raster to the buffered province boundary
        print(f"    Cropping {month_name} to {PROVINCE_NAME} boundary...")
        with rasterio.open(temp_path) as src:
            temp_data = src.read(1)
            temp_profile = src.profile

            # Create empty destination array
            destination_array = np.full(
                (DEM_HEIGHT, DEM_WIDTH), NO_DATA_VALUE, dtype=temp_data.dtype)

            # Resample to match DEM resolution.
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
                        province_shp_proj.geometry,
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
                no_data_mask, masked_data + DEM_DATA * T_LAPSE_RATE, NO_DATA_VALUE)

            # Save final temperature raster
            output_path = os.path.join(
                var_results_subdir, f"{var}_{month_name}.tif")
            with rasterio.open(output_path, 'w', **masked_profile) as dst:
                dst.write(final_data)

        os.remove(temp_path)

    print(f"Completed processing {var} temperature data")

print("Temperature data processing complete!")

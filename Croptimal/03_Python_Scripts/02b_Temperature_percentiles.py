import os
import glob
import calendar
import rasterio
import numpy as np
import xarray as xr
import geopandas as gpd
from rasterio.warp import reproject, Resampling
from rasterio.mask import mask
from collections import defaultdict
import matplotlib.pyplot as plt

"""
This script calculates temperature percentiles for each month by:
1. Processing each temperature variable (tmin, tmax, tavg)
2. For each variable, calculating specific percentiles (5th, 25th, etc.)
3. Grouping data by month and calculating percentiles across years
4. Cropping results to a buffered province boundary
5. Resampling to match a reference DEM
6. Saving the processed percentile data
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
PROVINCE_NAME = "Zaire"                         # dummy variable for testing.

# Define other folders
DATA_DIR = os.path.join(angola_wd, "01_Data")
GIS_DIR = os.path.join(angola_wd, "02_GIS")     # Directory with shapefiles
RESULTS_DIR = os.path.join(parent_wd, "04_Results", PROVINCE_NAME)
TEMP_DIR = os.path.join(parent_wd, "05_Temp")

# Define constants
RES = 250  # Resolution in meters
T_VARS = ["tavg", "tmin", "tmax"]  # Temperature variables
T_PERC = [0.05, 0.25]  # Percentiles to calculate (5% and 25%)
NO_DATA_VALUE = -9999.0  # No data value


# Load reference DEM
print("Loading reference DEM...")
dem_path = os.path.join(
    RESULTS_DIR, "DEM", f"DEM_{PROVINCE_NAME}_{RES}m.tif")  # used to be: "DEM_{PROVINCE_NAME}_{RES}m_diff.tif"
with rasterio.open(dem_path) as dem_src:
    DEM_PROFILE = dem_src.profile.copy()
    DEM_TRANSFORM = dem_src.transform
    DEM_CRS = dem_src.crs
    DEM_SHAPE = (dem_src.height, dem_src.width)
    DEM_HEIGHT = dem_src.height
    DEM_WIDTH = dem_src.width

# Load province shapefile
print("Loading province shapefile...")
provinces_filepath = os.path.join(GIS_DIR, "Shapefiles", "AGO_adm1.shp")
provinces_shp = gpd.read_file(provinces_filepath)
province_shp_sel = provinces_shp[provinces_shp["NAME_1"] == PROVINCE_NAME]
province_shp_sel_reproj = province_shp_sel.to_crs(DEM_CRS)


####################################################################################################
########################## Process each temperature variable and percentile #######################
####################################################################################################
# Process each temperature variable (tavg, tmin, tmax)
for var in T_VARS:
    print(f"Processing {var} temperature data...")

    # Process each percentile
    for p_t in T_PERC:
        p_str = str(p_t).replace(".", "")
        print(f"  Calculating {p_t*100}% percentile for {var}...")

        # Create output directory
        new_results_subdirults_subdir = os.path.join(RESULTS_DIR,
                                                     "Temperature", var, f"{p_str}perc")
        os.makedirs(new_results_subdirults_subdir, exist_ok=True)

        # Get input files for this temperature variable
        input_files = glob.glob(os.path.join(
            DATA_DIR, "Temperature", var, "*.nc"))

        if not input_files:
            print(f"  No {var} files found. Skipping.")
            continue

        # Create a dictionary to store data for each month
        month_data_dict = defaultdict(list)

        # Process each NetCDF file. Each file contains one year of data, and data is stored per month.
        for file in input_files:  # Go over all the year.
            try:
                # Open the NetCDF file
                with xr.open_dataset(file) as ds:
                    data = ds["t2m"]
                    lon = ds.longitude
                    lat = ds.latitude

                    # Extract data for each month from the dataset(1-12)
                    for month in range(1, 13):
                        # Filter data for this month
                        month_data = ds['t2m'].sel(
                            time=ds.time.dt.month == month)

                        if not month_data.size:
                            continue

                        # Convert to numpy array and add to the list for this month
                        month_data_dict[month].append(month_data.values)
            except Exception as e:
                print(f"  Error processing {file}: {e}")
                continue

        # Dictionary to store percentile results for each month
        out_layers = {}

        # Calculate percentiles for each month
        for month, data_list in month_data_dict.items():
            if not data_list:
                print(f"  No data found for month {month}")
                continue

            try:
                # Stack all data for this month across all years/files
                stacked_data = np.concatenate(
                    [d.reshape(d.shape[0], -1) for d in data_list], axis=0)

                # Calculate the specified percentile
                percentile_data = np.nanpercentile(
                    stacked_data, p_t * 100, axis=0)

                # Reshape back to the original spatial dimensions
                original_shape = data_list[0].shape[1:]  # Spatial dimensions
                percentile_data_reshaped = percentile_data.reshape(
                    original_shape)

                # Store the result for this month
                out_layers[month] = percentile_data_reshaped
            except Exception as e:
                print(f"  Error calculating percentile for month {month}: {e}")
                continue

        # If we don't have any processed data, skip to next percentile
        if not out_layers:
            print(
                f"  No percentile data could be calculated for {var} at {p_t*100}%")
            continue

        ####################################################################################################
        ########################## Crop, resample and save results #########################################
        ####################################################################################################
        print(f"  Resampling and cropping results for {var} at {p_t*100}%...")

        # Get metadata from one of the input files to use for creating GeoTIFFs
        # This is needed because NetCDF files may not directly provide the necessary GeoTIFF metadata
        sample_file = None
        for file in glob.glob(os.path.join(RESULTS_DIR, "Temperature", var, "Monthly_Mean", "*.tif*")):
            sample_file = file
            break

        # Get metadata from sample file
        with rasterio.open(sample_file) as src:
            sample_meta = src.meta.copy()
            sample_transform = src.transform
            sample_crs = src.crs

        # Process each month
        for month, percentile_data in out_layers.items():
            month_name = calendar.month_abbr[month]
            print(f"    Processing {month_name}...")

            # Create a temporary file to hold the percentile data
            temp_path = os.path.join(TEMP_DIR, f"temp_{month_name}.tif")

            # Create a profile for the temporary file
            temp_profile = sample_meta.copy()
            temp_profile.update({
                'height': percentile_data.shape[0],
                'width': percentile_data.shape[1],
                'count': 1,
            })

            # Write the percentile data to the temporary file
            with rasterio.open(temp_path, 'w', **temp_profile) as dst:
                dst.write(percentile_data.astype(rasterio.float32), 1)

            # Create an output array for the resampled data
            resampled_data = np.full(
                DEM_SHAPE, NO_DATA_VALUE, dtype=np.float32)

            # Reproject and resample the data to match the DEM
            with rasterio.open(temp_path) as src:
                temp_data = src.read(1)
                reproject(
                    source=temp_data,
                    destination=resampled_data,
                    src_transform=src.transform,
                    src_crs=src.crs,
                    dst_transform=DEM_TRANSFORM,
                    dst_crs=DEM_CRS,
                    resampling=Resampling.bilinear
                )

            # Create a memory file with the reprojected data for masking
            with rasterio.MemoryFile() as memfile:
                with memfile.open(**DEM_PROFILE) as temp_dst:
                    temp_dst.write(resampled_data, 1)

                    # Mask the data to the province boundary
                    masked_data, masked_transform = mask(
                        temp_dst,
                        province_shp_sel_reproj.geometry,
                        crop=True,
                        nodata=NO_DATA_VALUE)

                    # Save mask profile
                    masked_profile = temp_dst.profile.copy()
                    masked_profile.update({
                        "height": masked_data.shape[1],
                        "width": masked_data.shape[2],
                        "transform": masked_transform,
                    })

            # Save the masked data to anew_results_subdirle
            output_path = os.path.join(new_results_subdir, f"{month_name}.tif")
            with rasterio.open(output_path, 'w', **masked_profile) as dst:
                dst.write(masked_data.astype(rasterio.float32))
            os.remove(temp_path)

    print(f"Saved all results vor {var}")
print("Temperature percentiles processing complete!")

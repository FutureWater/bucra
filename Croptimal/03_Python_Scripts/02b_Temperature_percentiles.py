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
# Define true constants (in UPPERCASE)
DATA_DIR = "your_data_directory"  # Replace with actual path
RESULTS_DIR = "your_results_directory"  # Replace with actual path
PROVINCE_NAME = "your_province_name"  # Replace with actual province name
RES = 250  # Resolution in meters
T_VARS = ["tavg", "tmin", "tmax"]  # Temperature variables
T_PERC = [0.05, 0.25]  # Percentiles to calculate (5% and 25%)

# Load province shapefile
print("Loading province shapefile...")
province_shp_path = os.path.join(DATA_DIR, f"{PROVINCE_NAME}.shp")
province_shp = gpd.read_file(province_shp_path)

# Create a buffer around province shapefile
print("Creating buffer around province...")
province_buffer = province_shp.buffer(0.25)  # Buffer of 0.25 degrees

# Load reference DEM
print("Loading reference DEM...")
dem_path = os.path.join(
    RESULTS_DIR, PROVINCE_NAME, "DEM", f"DEM_{PROVINCE_NAME}_{RES}m_diff.tif")
with rasterio.open(dem_path) as dem_src:
    dem_meta = dem_src.meta.copy()
    dem_transform = dem_src.transform
    dem_crs = dem_src.crs
    dem_shape = (dem_src.height, dem_src.width)

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
        new_dir = os.path.join(RESULTS_DIR, PROVINCE_NAME,
                               "Temperature", var, f"{p_str}perc")
        os.makedirs(new_dir, exist_ok=True)

        # Get input files for this temperature variable
        input_files = glob.glob(os.path.join(
            DATA_DIR, "Temperature", var, "*.nc"))

        if not input_files:
            print(f"  No {var} files found. Skipping.")
            continue

        # Create a dictionary to store data for each month
        month_data_dict = defaultdict(list)

        # Process each NetCDF file
        print(f"  Reading NetCDF files for {var}...")

        # We'll use xarray which is well-suited for NetCDF data
        for file in input_files:
            try:
                # Open the NetCDF file
                with xr.open_dataset(file) as ds:
                    # Assuming the variable name is the same as var (tavg, tmin, tmax)
                    # Extract data for each month (1-12)
                    for month in range(1, 13):
                        # Filter data for this month
                        # Note: This assumes there's a 'time' dimension that can be used to filter by month
                        # The exact filtering depends on the structure of your NetCDF files
                        month_data = ds[var].sel(
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
            print(f"  Processing month {month}...")

            if not data_list:
                print(f"  No data found for month {month}")
                continue

            try:
                # Stack all data for this month across all years/files
                stacked_data = np.concatenate(
                    [d.reshape(d.shape[0], -1) for d in data_list], axis=0)

                # Calculate the specified percentile
                # Equivalent to calc(st, fun=function(x) raster::quantile(x, probs=p_T, na.rm=T))
                percentile_data = np.nanpercentile(
                    stacked_data, p_t * 100, axis=0)

                # Reshape back to the original spatial dimensions
                # Note: This assumes all files have the same spatial dimensions
                # You might need to adjust this based on your data structure
                original_shape = data_list[0].shape[1:]  # Spatial dimensions
                percentile_data = percentile_data.reshape(original_shape)

                # Store the result for this month
                out_layers[month] = percentile_data
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
        print(f"  Processing and saving results for {var} at {p_t*100}%...")

        # Get metadata from one of the input files to use for creating GeoTIFFs
        # This is needed because NetCDF files may not directly provide the necessary GeoTIFF metadata
        sample_file = None
        for file in glob.glob(os.path.join(RESULTS_DIR, PROVINCE_NAME, "Temperature", var, "Monthly_Mean", "*.tif*")):
            sample_file = file
            break

        if not sample_file:
            print(
                f"  Could not find sample GeoTIFF for {var}. Make sure to run the monthly mean processing first.")
            continue

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
            temp_path = os.path.join(new_dir, f"temp_{month_name}.tif")

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

            # Crop to province boundary
            with rasterio.open(temp_path) as src:
                # Ensure CRS compatibility
                if src.crs != province_shp.crs:
                    province_buffer_projected = province_buffer.to_crs(src.crs)
                else:
                    province_buffer_projected = province_buffer

                out_image, out_transform = mask(
                    src, province_buffer_projected.geometry, crop=True)
                out_meta = src.meta.copy()
                out_meta.update({
                    "driver": "GTiff",
                    "height": out_image.shape[1],
                    "width": out_image.shape[2],
                    "transform": out_transform,
                })

            # Save cropped data to another temporary file
            cropped_path = os.path.join(new_dir, f"cropped_{month_name}.tif")
            with rasterio.open(cropped_path, 'w', **out_meta) as dst:
                dst.write(out_image[0], 1)

            # Resample to match DEM resolution and extent
            # Equivalent to projectRaster(crop_var_raster, resampled_dem_diff_ras, res=res, method="bilinear")
            resampled_data = np.zeros(dem_shape, dtype=np.float32)

            with rasterio.open(cropped_path) as src:
                reproject(
                    source=rasterio.band(src, 1),
                    destination=resampled_data,
                    src_transform=src.transform,
                    src_crs=src.crs,
                    dst_transform=dem_transform,
                    dst_crs=dem_crs,
                    resampling=Resampling.bilinear
                )

            # Save the final resampled file
            output_path = os.path.join(new_dir, f"{month_name}.tif")
            with rasterio.open(output_path, 'w', **dem_meta) as dst:
                dst.write(resampled_data.astype(rasterio.float32), 1)

            print(f"    Saved {month_name}.tif")

            # Clean up temporary files
            try:
                os.remove(temp_path)
                os.remove(cropped_path)
            except Exception as e:
                print(f"    Warning: Could not remove temporary files: {e}")

print("Temperature percentiles processing complete!")

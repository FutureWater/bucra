import os
import glob
import calendar
import rasterio
import numpy as np
import netCDF4 as nc
import xarray as xr
import geopandas as gpd
from rasterio.warp import reproject, Resampling
from rasterio.mask import mask

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
# Define directories and file paths (only true constants in UPPERCASE)
DATA_DIR = "your_data_directory"  # Replace with actual path
RESULTS_DIR = "your_results_directory"  # Replace with actual path
PROVINCE_NAME = "your_province_name"  # Replace with actual province name
RES = 250  # Resolution in meters
LOCAL_PROJ = "EPSG:32733"  # Local projection

# Temperature variables and parameters
T_VARS = ["tavg", "tmin", "tmax"]  # Temperature variables to process
T_LAPSE_RATE = -0.0065  # Temperature lapse rate (°C/m)

# Month abbreviations for output file naming
month_abbrs = [calendar.month_abbr[i] for i in range(1, 13)]

# Load province shapefile
print("Loading province shapefile...")
province_shp_path = os.path.join(DATA_DIR, f"{PROVINCE_NAME}.shp")
province_shp = gpd.read_file(province_shp_path)

# Create a buffer around province shapefile (equivalent to gBuffer in R)
print("Creating buffer around province...")
province_buffer = province_shp.buffer(0.25)  # Buffer of 0.25 degrees

# Load reference DEM
print("Loading reference DEM...")
dem_path = os.path.join(RESULTS_DIR, PROVINCE_NAME, "DEM",
                        f"DEM_{PROVINCE_NAME}_{RES}m_diff.tif")
with rasterio.open(dem_path) as dem_src:
    dem_data = dem_src.read(1)
    dem_meta = dem_src.meta.copy()
    dem_transform = dem_src.transform
    dem_crs = dem_src.crs
    dem_shape = (dem_src.height, dem_src.width)

####################################################################################################
################################# Process temperature variables ###################################
####################################################################################################
# Loop through each temperature variable (tavg, tmin, tmax)
for var in T_VARS:
    print(f"Processing {var} temperature data...")

    # Get input files for this temperature variable
    input_files = glob.glob(os.path.join(DATA_DIR, "Temperature", var, "*.nc"))

    if not input_files:
        print(f"No {var} files found. Skipping.")
        continue

    # Create output directory
    new_dir = os.path.join(RESULTS_DIR, PROVINCE_NAME,
                           "Temperature", var, "Monthly_Mean")
    os.makedirs(new_dir, exist_ok=True)

    # Read and process NetCDF data
    # Note: This is a simplified approach to replicate R's stack() and stackApply()
    # functions. The exact implementation depends on the structure of your NetCDF files.
    print(f"  Reading NetCDF files for {var}...")

    # We'll use xarray which is well-suited for NetCDF data
    datasets = []
    for file in input_files:
        with xr.open_dataset(file) as ds:
            datasets.append(ds)

    # Combine all datasets
    combined_data = xr.concat(datasets, dim="time")

    # Group by month and calculate mean for each month
    # This replicates the stackApply function in R
    print(f"  Calculating monthly means for {var}...")
    monthly_means = combined_data.groupby('time.month').mean(dim='time')

    # Process each month
    for month_idx in range(1, 13):
        month_name = month_abbrs[month_idx - 1]
        print(f"Processing {month_name}...")

        # Extract data for this month
        try:
            # Get the temperature variable (assumes the main variable has the same name as var)
            month_data = monthly_means[var].sel(month=month_idx).values
        except:
            print(
                f"    Warning: Could not extract {month_name} data. Check variable names in NetCDF.")
            continue

        # Create a temporary raster file for cropping and reprojection
        temp_path = os.path.join(new_dir, f"temp_{month_name}.tif")

        # Get the transform and CRS information from the first file
        # Note: This assumes all files have the same geospatial reference
        with xr.open_dataset(input_files[0]) as ds:
            # Create a profile for the raster
            # Note: This part may need adjustment based on your NetCDF structure
            temp_profile = {
                'driver': 'GTiff',
                'height': month_data.shape[0],
                'width': month_data.shape[1],
                'count': 1,
                'dtype': rasterio.float32,
                'crs': ds.rio.crs,
                'transform': ds.rio.transform()
            }

        # Write to temporary file
        with rasterio.open(temp_path, 'w', **temp_profile) as dst:
            dst.write(month_data.astype(rasterio.float32), 1)

        # Crop raster to the buffered province boundary
        print(f"    Cropping {month_name} to province boundary...")
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
        print(f"    Resampling {month_name} to match DEM...")
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

        # Apply lapse rate correction based on elevation difference
        # This is the equivalent of: final_ras <- resampled_var_ras + resampled_dem_diff_ras*T_lapse_rate
        print(f"    Applying lapse rate correction for {month_name}...")
        final_data = resampled_data + dem_data * T_LAPSE_RATE

        # Save final temperature raster
        output_path = os.path.join(new_dir, f"{month_name}.tiff")
        with rasterio.open(output_path, 'w', **dem_meta) as dst:
            dst.write(final_data.astype(rasterio.float32), 1)

        # Clean up temporary files
        try:
            os.remove(temp_path)
            os.remove(cropped_path)
        except Exception as e:
            print(f"    Warning: Could not remove temporary files: {e}")

    print(f"Completed processing {var} temperature data")

print("Temperature data processing complete!")

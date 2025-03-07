import os
import glob
import rasterio
import numpy as np
import geopandas as gpd
import calendar
from rasterio.warp import reproject, Resampling
from rasterio.mask import mask

"""
This script calculates rainfall percentiles by:
1. Loading monthly rainfall data
2. Calculating specified percentiles (5th and 25th) for each month
3. Cropping to a buffered province boundary
4. Resampling to match a reference DEM
5. Saving the processed data
"""

####################################################################################################
########################## Define directories and file paths #######################################
####################################################################################################
# Define constants (in UPPERCASE)
DATA_DIR = "your_data_directory"  # Replace with actual path
RESULTS_DIR = "your_results_directory"  # Replace with actual path
PROVINCE_NAME = "your_province_name"  # Replace with actual province name
RES = 250  # Resolution in meters
P_PERC = [0.05, 0.25]  # Percentiles to calculate (5% and 25%)

# Get input rainfall files
input_files = glob.glob(os.path.join(
    DATA_DIR, "Rainfall", "Monthly_Sum", "*.tif"))
if not input_files:
    raise FileNotFoundError(
        f"No rainfall files found in {os.path.join(DATA_DIR, 'Rainfall', 'Monthly_Sum')}")

print(f"Found {len(input_files)} rainfall files")

# Extract month numbers from filenames
months_stack = [int(os.path.basename(file)[5:7]) for file in input_files]

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
########################## Process each percentile ################################################
####################################################################################################
# Process each percentile
for p_p in P_PERC:
    # Generate directory name by removing decimal point from percentile value
    p_str = str(p_p).replace(".", "")
    new_dir = os.path.join(RESULTS_DIR, PROVINCE_NAME,
                           "Rainfall", f"{p_str}perc")
    os.makedirs(new_dir, exist_ok=True)

    print(f"Processing {p_p*100}% percentile...")

    # Dictionary to store percentile results for each month
    out_layers = {}

    # Process each month (1-12)
    for month in range(1, 13):
        print(f"  Processing month {month}...")

        # Get files for this month
        month_files = [f for i, f in enumerate(
            input_files) if months_stack[i] == month]

        if not month_files:
            print(f"  No data found for month {month}")
            continue

        # Read all rasters for this month
        month_data = []
        for file in month_files:
            with rasterio.open(file) as src:
                # If this is the first file for this month, save the metadata
                if not month_data:
                    month_meta = src.meta.copy()
                    month_profile = src.profile.copy()
                    month_transform = src.transform
                    month_crs = src.crs

                # Read the data and append to our list
                month_data.append(src.read(1))

        # Calculate the percentile for this month across all years
        # This is equivalent to R's calc function with quantile
        if month_data:
            # Stack all data for this month
            stacked_data = np.stack(month_data)

            # Calculate the specified percentile (equivalent to quantile in R)
            percentile_data = np.nanpercentile(stacked_data, p_p * 100, axis=0)

            # Store the result for this month
            out_layers[month] = percentile_data

    # Crop each monthly percentile raster to the province boundary
    print("  Cropping raster to province boundary...")

    # Create a profile for temporary files
    temp_profile = month_profile.copy()

    # List to store cropped data
    cropped_layers = {}

    # Process each month
    for month, data in out_layers.items():
        # Create a temporary raster file for cropping
        temp_path = os.path.join(new_dir, f"temp_month_{month}.tif")

        with rasterio.open(temp_path, 'w', **temp_profile) as dst:
            dst.write(data.astype(rasterio.float32), 1)

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

            # Store the cropped data
            cropped_layers[month] = {
                'data': out_image[0],
                'transform': out_transform,
                'meta': out_meta
            }

        # Clean up temporary file
        try:
            os.remove(temp_path)
        except:
            print(f"  Warning: Could not remove temporary file {temp_path}")

    # Resample each cropped percentile raster to match DEM
    print("  Resampling to match DEM...")

    # Month abbreviations for naming output files
    month_abbrs = [calendar.month_abbr[i] for i in range(1, 13)]

    # Process each month
    for month, layer in cropped_layers.items():
        # Create an output array for the resampled data
        resampled_data = np.zeros(dem_shape, dtype=np.float32)

        # Reproject the data to match the DEM
        reproject(
            source=layer['data'],
            destination=resampled_data,
            src_transform=layer['transform'],
            src_crs=month_crs,  # Assuming all monthly rasters have the same CRS
            dst_transform=dem_transform,
            dst_crs=dem_crs,
            resampling=Resampling.bilinear
        )

        # Write the result to a file
        # Use month abbreviation for the filename
        output_path = os.path.join(new_dir, f"{month_abbrs[month-1]}.tif")

        with rasterio.open(output_path, 'w', **dem_meta) as dst:
            dst.write(resampled_data.astype(rasterio.float32), 1)

        print(f"  Saved {month_abbrs[month-1]} for {p_p*100}% percentile")

print("Rainfall percentiles processing complete!")

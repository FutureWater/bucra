import os
import glob
import rasterio
import numpy as np
from rasterio.warp import reproject, Resampling
from rasterio.windows import Window
from rasterio.windows import from_bounds


"""
This script processes NDVI files by:
1. Defining directories and file paths
2. Creating output directories
3. Loading input NDVI files
4. Importing DEM for reference extent and resolution
5. Processing each NDVI file individually to avoid memory issues:
   a. Cropping raster to extent of the DEM
   b. Reprojecting/resampling data to match DEM
   c. Saving the resampled and cropped raster
"""
####################################################################################################
########################## Define directories and file paths #######################################
####################################################################################################
# Define directories and file paths
DATA_DIR = "your_data_directory"
RESULTS_DIR = "your_results_directory"
PROVINCE_NAME = "your_province_name"
RES = 250  # Set resolution

# Create output directory
NEW_DIR = os.path.join(RESULTS_DIR, PROVINCE_NAME, "NDVI", "Mean_Monthly")
os.makedirs(NEW_DIR, exist_ok=True)

# Get input files
INPUT_FILES = glob.glob(os.path.join(
    DATA_DIR, "NDVI", PROVINCE_NAME, "Mean_Monthly", "*.tif"))
NAMES_RASTER = [os.path.basename(f) for f in INPUT_FILES]

# Import DEM for reference extent and resolution
DEM_DIFF_PATH = os.path.join(
    RESULTS_DIR, PROVINCE_NAME, "DEM", f"DEM_{PROVINCE_NAME}_{RES}m_diff.tif")

print(f"Loading reference DEM: {DEM_DIFF_PATH}")
with rasterio.open(DEM_DIFF_PATH) as dem_src:
    DEM_META = dem_src.meta.copy()
    DEM_BOUNDS = dem_src.bounds
    DEM_TRANSFORM = dem_src.transform
    DEM_CRS = dem_src.crs
    DEM_HEIGHT = dem_src.height
    DEM_WIDTH = dem_src.width

####################################################################################################
################################# Resample NDVI files ##############################################
####################################################################################################
# Process each NDVI file individually to avoid memory issues
for i, input_file in enumerate(INPUT_FILES):
    print(
        f"Processing {os.path.basename(input_file)} ({i+1}/{len(INPUT_FILES)})")

    with rasterio.open(input_file) as src:
        # Step 1: Crop raster to extent of the DEM
        # Create a window that represents the DEM bounds in the source raster's coordinate system
        src_bounds = src.bounds

        # Check for overlap
        if (src_bounds.left > DEM_BOUNDS.right or src_bounds.right < DEM_BOUNDS.left or
                src_bounds.bottom > DEM_BOUNDS.top or src_bounds.top < DEM_BOUNDS.bottom):
            print(
                f"Warning: {os.path.basename(input_file)} does not overlap with the DEM. Skipping.")
            continue

        # Create a destination array of the same shape as the DEM
        dst_data = np.zeros((DEM_HEIGHT, DEM_WIDTH), dtype=rasterio.float32)

        # Step 2: Reproject/resample data to match DEM
        print(f"  Resampling to match DEM...")
        reproject(
            source=rasterio.band(src, 1),
            destination=dst_data,
            src_transform=src.transform,
            src_crs=src.crs,
            dst_transform=DEM_TRANSFORM,
            dst_crs=DEM_CRS,
            resampling=Resampling.bilinear
        )

        # Step 3: Save the resampled and cropped raster
        output_path = os.path.join(
            NEW_DIR, f"{os.path.splitext(NAMES_RASTER[i])[0]}.tiff")
        output_meta = DEM_META.copy()
        output_meta.update({
            "driver": "GTiff",
            "height": DEM_HEIGHT,
            "width": DEM_WIDTH,
            "transform": DEM_TRANSFORM,
            "crs": DEM_CRS
        })

        print(f"  Saving to {output_path}")
        with rasterio.open(output_path, 'w', **output_meta) as dst:
            dst.write(dst_data, 1)

print("NDVI resampling complete!")

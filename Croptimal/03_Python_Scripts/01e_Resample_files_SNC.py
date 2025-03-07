import os
import glob
import rasterio
import numpy as np
from rasterio.warp import reproject, Resampling
from rasterio.mask import mask
import geopandas as gpd
from rasterio.merge import merge

"""
This script processes Soil Nutrient Content (SNC) data by:
1. Reading soil nutrient content TIF files
2. Cropping to a buffered province boundary
3. Resampling to match a reference DEM
4. Saving the processed data
"""

####################################################################################################
########################## Define directories and file paths #######################################
####################################################################################################
# Define directories and file paths (constants in UPPERCASE)
DATA_DIR = "your_data_directory"  # Replace with actual path
RESULTS_DIR = "your_results_directory"  # Replace with actual path
PROVINCE_NAME = "your_province_name"  # Replace with actual province name
RES = 250  # Resolution in meters
LOCAL_PROJ = "EPSG:32733"  # Local projection

# Create output directory for processed SNC data
NEW_DIR = os.path.join(RESULTS_DIR, PROVINCE_NAME, "Soil_Nutrient_Content")
os.makedirs(NEW_DIR, exist_ok=True)

# Get input SNC files
# Note: In R, paste0(data_dir, "Soil_Nutrient_Content\\") with pattern = "\\.tif$" is used
# In Python, we use os.path.join and glob for similar functionality
INPUT_FILES = glob.glob(os.path.join(
    DATA_DIR, "Soil_Nutrient_Content", "*.tif"))

if not INPUT_FILES:
    raise FileNotFoundError(
        f"No SNC files found in {os.path.join(DATA_DIR, 'Soil_Nutrient_Content')}")

print(f"Found {len(INPUT_FILES)} SNC files")

# Get file names without extension and replace "af" with province_name
VAR_NAMES = [os.path.basename(file) for file in INPUT_FILES]
VAR_NAMES_2 = [name.replace("af", f"{PROVINCE_NAME}_") for name in VAR_NAMES]

####################################################################################################
################################# Process SNC files ################################################
####################################################################################################
# Load province shapefile
print("Loading province shapefile...")
PROVINCE_SHP_PATH = os.path.join(
    DATA_DIR, f"{PROVINCE_NAME}.shp")  # Adjust path as needed
province_shp = gpd.read_file(PROVINCE_SHP_PATH)

# Create a buffer around province shapefile (equivalent to gBuffer in R)
print("Creating buffer around province...")
province_buffer = province_shp.buffer(0.25)  # Buffer of 0.25 degrees

# Load reference DEM
print("Loading reference DEM...")
DEM_PATH = os.path.join(RESULTS_DIR, PROVINCE_NAME, "DEM",
                        f"DEM_{PROVINCE_NAME}_{RES}m_diff.tif")
with rasterio.open(DEM_PATH) as dem_src:
    dem_meta = dem_src.meta.copy()
    dem_transform = dem_src.transform
    dem_crs = dem_src.crs
    dem_shape = (dem_src.height, dem_src.width)

# Process each SNC file
for i, (input_file, output_name) in enumerate(zip(INPUT_FILES, VAR_NAMES_2)):
    print(
        f"Processing {os.path.basename(input_file)} ({i+1}/{len(INPUT_FILES)})")

    # Read input file
    with rasterio.open(input_file) as src:
        # Ensure CRS compatibility
        if src.crs != province_shp.crs:
            province_buffer_projected = province_buffer.to_crs(src.crs)
        else:
            province_buffer_projected = province_buffer

        # Crop raster to the buffered province boundary
        print(f"  Cropping to province boundary...")
        out_image, out_transform = mask(
            src, province_buffer_projected.geometry, crop=True)
        out_meta = src.meta.copy()
        out_meta.update({
            "driver": "GTiff",
            "height": out_image.shape[1],
            "width": out_image.shape[2],
            "transform": out_transform,
        })

    # Need a temporary file for the cropped data
    temp_path = os.path.join(NEW_DIR, f"temp_{os.path.basename(input_file)}")
    with rasterio.open(temp_path, 'w', **out_meta) as temp:
        temp.write(out_image[0], 1)

    # Resample to match DEM resolution and extent
    print(f"  Resampling to match DEM...")
    resampled_data = np.zeros(DEM_SHAPE, dtype=np.float32)

    # Now reproject from the cropped temp file to match DEM
    with rasterio.open(temp_path) as src:
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
    OUTPUT_PATH = os.path.join(NEW_DIR, output_name)
    out_meta = dem_meta.copy()
    with rasterio.open(OUTPUT_PATH, 'w', **out_meta) as dst:
        dst.write(resampled_data.astype(rasterio.float32), 1)

    # Clean up temporary file
    try:
        os.remove(temp_path)
    except Exception as e:
        print(f"  Warning: Could not remove temporary file: {e}")

print("SNC resampling complete!")

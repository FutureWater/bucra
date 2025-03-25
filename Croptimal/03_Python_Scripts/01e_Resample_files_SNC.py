import os
import glob
import rasterio
import numpy as np
from rasterio.warp import reproject, Resampling, calculate_default_transform
from rasterio.mask import mask
import geopandas as gpd
from rasterio.merge import merge
import matplotlib.pyplot as plt

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
NO_DATA_VALUE = -9999.0

# Create output directory for processed SNC data
SNC_RESULTS_DIR = os.path.join(RESULTS_DIR, "Soil_Nutrient_Content")
os.makedirs(SNC_RESULTS_DIR, exist_ok=True)

# Get input SNC files
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
provinces_filepath = os.path.join(GIS_DIR, "Shapefiles", "AGO_adm1.shp")
provinces_shp = gpd.read_file(provinces_filepath)
province_shp_sel = provinces_shp[provinces_shp["NAME_1"] == PROVINCE_NAME]
province_shp_proj = province_shp_sel.to_crs(LOCAL_PROJ)

# Create a buffer around province shapefile
print("Creating buffer around province...")
province_buffer = province_shp_proj.buffer(500)  # Buffer of 500 meters

# Load reference DEM
print("Loading reference DEM...")
DEM_PATH = os.path.join(RESULTS_DIR, "DEM",
                        f"DEM_{PROVINCE_NAME}_{RES}m.tif")  # used to be: "DEM_{PROVINCE_NAME}_{RES}m_diff.tif"
with rasterio.open(DEM_PATH) as dem_src:
    DEM_META = dem_src.meta.copy()
    DEM_TRANSFORM = dem_src.transform
    DEM_CRS = dem_src.crs
    DEM_SHAPE = (dem_src.height, dem_src.width)

# Process each SNC file
for i, (input_file, output_name) in enumerate(zip(INPUT_FILES, VAR_NAMES_2)):
    print(
        f"Processing {os.path.basename(input_file)} ({i+1}/{len(INPUT_FILES)})")

    # Read input file
    with rasterio.open(input_file) as src:
        # Ensure CRS compatibility
        if src.crs != province_shp_proj.crs:
            province_buffer_projected = province_buffer.to_crs(src.crs)
        else:
            province_buffer_projected = province_buffer

        # Crop raster to the buffered province boundary
        print("Cropping to province boundary...")
        out_image, out_transform = mask(
            src, province_buffer_projected.geometry, crop=True, nodata=NO_DATA_VALUE)

        # Save metadata of raster
        out_meta = src.meta.copy()

        # Get transform and shape new raster converted to local projection
        dst_transform, dst_width, dst_height = calculate_default_transform(
            src.crs, DEM_CRS, src.width, src.height,
            *rasterio.transform.array_bounds(out_image.shape[1], out_image.shape[2], out_transform),
            resolution=RES
        )

        # Resample to match DEM resolution and extent
        print("  Resampling to match DEM...")
        destination_array = np.full(DEM_SHAPE, NO_DATA_VALUE, dtype=np.float32)

        out_meta.update({
            "crs": DEM_CRS,
            "transform": dst_transform,
            "width": dst_width,
            "height": dst_height,
            "nodata": NO_DATA_VALUE
        })

        OUTPUT_PATH = os.path.join(SNC_RESULTS_DIR, output_name)

        with rasterio.open(OUTPUT_PATH, 'w', **out_meta) as dst:
            reproject(
                source=out_image[0],
                destination=destination_array,
                src_transform=out_transform,
                src_crs=src.crs,
                dst_transform=DEM_TRANSFORM,
                dst_crs=DEM_CRS,
                resampling=Resampling.bilinear,
                src_nodata=NO_DATA_VALUE,
                dst_nodata=NO_DATA_VALUE
            )
            dst.write(destination_array.astype(rasterio.float32), 1)

print("SNC resampling complete!")

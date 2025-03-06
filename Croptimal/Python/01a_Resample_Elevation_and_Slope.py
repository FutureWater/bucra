import os
import math
import numpy as np
import rasterio
from rasterio.warp import reproject, Resampling, calculate_default_transform
from rasterio.mask import mask
import geopandas as gpd
from scipy.ndimage import sobel


##############################################################################################
################################### START OF DATA INPUT ######################################
##############################################################################################
# Define directories and file paths
DATA_DIR = "your_data_directory"
RESULTS_DIR = "your_results_directory"
PROVINCE_NAME = "your_province_name"
RES = 250  # Set resolution in meters
LOCAL_PROJ = "EPSG:32733"

# Function to ensure directory exists


def ensure_dir(directory):
    os.makedirs(directory, exist_ok=True)


# Create output directory for DEM
NEW_DIR = os.path.join(RESULTS_DIR, PROVINCE_NAME, "DEM")
ensure_dir(NEW_DIR)

# Create output directory for slope
INDIR_ELEV = os.path.join(DATA_DIR, "_LS_Results", "Elevation")
ensure_dir(INDIR_ELEV)

# Load original DEM (assuming it's already available)
DEM_PATH = os.path.join(DATA_DIR, "DEM.tif")  # Adjust this path as needed
print(f"Loading DEM: {DEM_PATH}")

# Load province shapefile
PROVINCE_SHP_PATH = os.path.join(
    DATA_DIR, f"{PROVINCE_NAME}.shp")  # Adjust as needed
PROVINCE_SHP = gpd.read_file(PROVINCE_SHP_PATH)

##############################################################################################
############################# CALCULATE ELEVATION AND SLOPE ##################################
##############################################################################################
# Step 1: Crop and mask DEM with province shapefile
print(f"Crop DEM: {PROVINCE_NAME}")
with rasterio.open(DEM_PATH) as src:
    # Project shapefile to match DEM CRS if needed
    if PROVINCE_SHP.crs != src.crs:
        PROVINCE_SHP = PROVINCE_SHP.to_crs(src.crs)

    # Crop DEM to shapefile extent
    out_image, out_transform = mask(src, PROVINCE_SHP.geometry, crop=True)

    # Copy metadata
    out_meta = src.meta.copy()
    out_meta.update({
        "driver": "GTiff",
        "height": out_image.shape[1],
        "width": out_image.shape[2],
        "transform": out_transform
    })

# Step 2: Prepare target raster with desired resolution and projection
print(f"Setting up target raster grid: {PROVINCE_NAME}")
# Project province shapefile to local projection for determining bounds
PROVINCE_SHP_PROJ = PROVINCE_SHP.to_crs(LOCAL_PROJ)
bounds = PROVINCE_SHP_PROJ.total_bounds  # [xmin, ymin, xmax, ymax]

# Calculate dimensions of target raster
width = int((bounds[2] - bounds[0]) / RES)
height = int((bounds[3] - bounds[1]) / RES)

# Create transform for target raster
target_transform = rasterio.transform.from_bounds(
    bounds[0], bounds[1], bounds[2], bounds[3], width, height
)

# Step 3: Reproject using bilinear method
print(f"Bilinear projectRaster DEM: {PROVINCE_NAME}")
bilinear_dem = np.zeros((height, width), dtype=np.float32)

with rasterio.open(DEM_PATH) as src:
    reproject(
        source=src.read(1),
        destination=bilinear_dem,
        src_transform=src.transform,
        src_crs=src.crs,
        dst_transform=target_transform,
        dst_crs=LOCAL_PROJ,
        resampling=Resampling.bilinear
    )

# Step 4: Reproject using nearest neighbor method
print(f"NGb projectRaster DEM: {PROVINCE_NAME}")
ngb_dem = np.zeros((height, width), dtype=np.float32)

with rasterio.open(DEM_PATH) as src:
    reproject(
        source=src.read(1),
        destination=ngb_dem,
        src_transform=src.transform,
        src_crs=src.crs,
        dst_transform=target_transform,
        dst_crs=LOCAL_PROJ,
        resampling=Resampling.nearest
    )

# Step 5: Calculate difference between bilinear and nearest neighbor
diff_dem = bilinear_dem - ngb_dem

# Step 6: Write results to files
print(f"Write rasters DEM: {PROVINCE_NAME}")
# Metadata for output files
out_meta = {
    "driver": "GTiff",
    "height": height,
    "width": width,
    "count": 1,
    "dtype": bilinear_dem.dtype,
    "crs": LOCAL_PROJ,
    "transform": target_transform
}

# Write bilinear resampled DEM
bilinear_path = os.path.join(NEW_DIR, f"DEM_{PROVINCE_NAME}_{RES}m.tif")
with rasterio.open(bilinear_path, "w", **out_meta) as dst:
    dst.write(bilinear_dem, 1)

# Write difference raster
diff_path = os.path.join(NEW_DIR, f"DEM_{PROVINCE_NAME}_{RES}m_diff.tif")
with rasterio.open(diff_path, "w", **out_meta) as dst:
    dst.write(diff_dem, 1)

# Step 7: Calculate and save slope
print(f"Calculate slope DEM: {PROVINCE_NAME}")

# Calculate slope using 3x3 windows (equivalent to terrain with neighbors=8)


def calculate_slope(dem, cell_size=RES):
    # Calculate gradients
    dx = sobel(dem, axis=1) / (8 * cell_size)
    dy = sobel(dem, axis=0) / (8 * cell_size)

    # Calculate slope in radians and convert to percent
    slope_radians = np.arctan(np.sqrt(dx**2 + dy**2))
    slope_percent = np.tan(slope_radians) * 100

    return slope_percent


# Calculate slope
slope = calculate_slope(bilinear_dem)

# Save slope raster
slope_path = os.path.join(INDIR_ELEV, f"Slope_{PROVINCE_NAME}.tif")
with rasterio.open(slope_path, "w", **out_meta) as dst:
    dst.write(slope.astype(rasterio.float32), 1)

# Step 8: Calculate areas with slope less than threshold
# Get slope threshold from params (you'll need to define this)
slope_threshold = 15  # Example value, replace with your actual threshold from params
slope_lower_limit = (slope < slope_threshold).astype(np.uint8)

# Save threshold raster
lower_slope_path = os.path.join(
    INDIR_ELEV, f"Slope_lower_{slope_threshold}perc_{PROVINCE_NAME}.tif")
out_meta.update({"dtype": "uint8"})
with rasterio.open(lower_slope_path, "w", **out_meta) as dst:
    dst.write(slope_lower_limit, 1)

print(f"Processing complete for {PROVINCE_NAME}")

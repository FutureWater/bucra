import os
import numpy as np
import pandas as pd
import rasterio
from rasterio.warp import reproject, Resampling, calculate_default_transform
from rasterio.mask import mask
import geopandas as gpd
from scipy.ndimage import sobel
import matplotlib.pyplot as plt


##############################################################################################
################################### START OF DATA INPUT ######################################
##############################################################################################
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
PARAMETERS = pd.read_csv(os.path.join(current_wd, "Parameters.csv"))
PROVINCES = os.path.join(GIS_DIR, "Shapefiles", "AGO_adm1.shp")

# Set resolution and projection
RES = 250                                       # Set resolution in meters
LOCAL_PROJ = "EPSG:32733"

# Function to ensure directory exists


def ensure_dir(directory):
    os.makedirs(directory, exist_ok=True)


# Create output directory for DEM
DEM_DATA_DIR = os.path.join(RESULTS_DIR, "DEM")
ensure_dir(DEM_DATA_DIR)

# Create output directory for slope
INDIR_ELEV = os.path.join(RESULTS_DIR, "_LS_Results", "Elevation")
ensure_dir(INDIR_ELEV)

# Load original DEM (assuming it's already available)
DEM_PATH = os.path.join(DATA_DIR, "Elevation",
                        "SRTM_30M_Angola.tif")  # Check which dem

# Load province shapefile
PROVINCES_FILEPATH = os.path.join(GIS_DIR, "Shapefiles", "AGO_adm1.shp")
provinces_shp = gpd.read_file(PROVINCES_FILEPATH)
PROVINCE_SHP_SEL = provinces_shp[provinces_shp["NAME_1"] == PROVINCE_NAME]

##############################################################################################
############################# CALCULATE ELEVATION AND SLOPE ##################################
##############################################################################################
# Step 1: Crop and mask DEM with province shapefile.
print(f"Crop DEM: {PROVINCE_NAME}")
with rasterio.open(DEM_PATH) as src:
    # Project shapefile to match DEM CRS if needed
    if PROVINCE_SHP_SEL.crs != src.crs:
        PROVINCE_SHP_SEL = PROVINCE_SHP_SEL.to_crs(src.crs)

    # Crop DEM to shapefile extent.
    out_image, out_transform = mask(src, PROVINCE_SHP_SEL.geometry, crop=True)

    # Copy metadata. Transform matrix is taken from DEM.
    out_meta = src.meta.copy()
    out_meta.update({
        "driver": "GTiff",
        "height": out_image.shape[1],
        "width": out_image.shape[2],
        # (pixel size x, row rotation, x-coordinate of upper-left corner, column rotation, pixel size y, y-coordinate of upper-left corner)
        "transform": out_transform
    })

    # Save cropped DEM to temporary file
    cropped_dem_path = os.path.join(
        TEMP_DIR, f"cropped_dem_{PROVINCE_NAME}_temp.tif")
    with rasterio.open(cropped_dem_path, "w", **out_meta) as dest:
        dest.write(out_image)

# Step 2: Prepare target raster with desired resolution and projection for resampling.
print(f"Setting up target raster grid: {PROVINCE_NAME}")
# Project province shapefile to local projection for determining bounds
PROVINCE_SHP_NEWPROJ = PROVINCE_SHP_SEL.to_crs(LOCAL_PROJ)
bounds = PROVINCE_SHP_NEWPROJ.total_bounds  # [xmin, ymin, xmax, ymax]

# Calculate dimensions of target raster
width = int((bounds[2] - bounds[0]) / RES)
height = int((bounds[3] - bounds[1]) / RES)

# Create transform (local projection) for target raster
target_transform = rasterio.transform.from_bounds(
    bounds[0], bounds[1], bounds[2], bounds[3], width, height
)  # (left, bottom, right, top, width, height)

# Step 3: Reproject using bilinear method
print(f"Bilinear projectRaster DEM: {PROVINCE_NAME}")
bilinear_dem = np.zeros((height, width), dtype=np.float32)

# Fill bilinear_dem array with resampled values.
with rasterio.open(cropped_dem_path) as src:
    reproject(
        source=src.read(1),
        destination=bilinear_dem,
        src_transform=src.transform,
        src_crs=src.crs,
        # Transform matrix from local projection Transform = (pixel size x, row rotation, x-coordinate of upper-left corner, column rotation, pixel size y, y-coordinate of upper-left corner)
        dst_transform=target_transform,
        dst_crs=LOCAL_PROJ,
        resampling=Resampling.bilinear
    )

# # Step 4: Reproject using nearest neighbor method
# print(f"Nearest Neighbor projectRaster DEM: {PROVINCE_NAME}")
# ngb_dem = np.zeros((height, width), dtype=np.float32)

# with rasterio.open(DEM_PATH) as src:
#     reproject(
#         source=src.read(1),
#         destination=ngb_dem,
#         src_transform=src.transform,
#         src_crs=src.crs,
#         dst_transform=target_transform,
#         dst_crs=LOCAL_PROJ,
#         resampling=Resampling.nearest
#     )

# Step 5: Calculate difference between bilinear and nearest neighbor
# Unclear why we need this for now.
# diff_dem = bilinear_dem - ngb_dem

# # Step 6: Write results to files
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
bilinear_path = os.path.join(DEM_DATA_DIR, f"DEM_{PROVINCE_NAME}_{RES}m.tif")
with rasterio.open(bilinear_path, "w", **out_meta) as dst:
    dst.write(bilinear_dem, 1)

# # Write difference raster
# diff_path = os.path.join(DEM_DATA_DIR, f"DEM_{PROVINCE_NAME}_{RES}m_diff.tif")
# with rasterio.open(diff_path, "w", **out_meta) as dst:
#     dst.write(diff_dem, 1)

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
slope_threshold = float(
    PARAMETERS.loc[PARAMETERS['Parameter'] == 'Slope', 'Limit'].values[0])
slope_lower_limit = (slope < slope_threshold).astype(np.uint8)

# Save threshold raster
print(f"Write rasters lower slope limit {PROVINCE_NAME}")
lower_slope_path = os.path.join(
    INDIR_ELEV, f"Slope_lower_{slope_threshold}perc_{PROVINCE_NAME}.tif")
out_meta.update({"dtype": "uint8"})
with rasterio.open(lower_slope_path, "w", **out_meta) as dst:
    dst.write(slope_lower_limit, 1)

# Remove temporary files
os.remove(cropped_dem_path)

print(f"Processing complete for {PROVINCE_NAME}")

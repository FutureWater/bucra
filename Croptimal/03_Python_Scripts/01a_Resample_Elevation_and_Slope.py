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
# Define current and parent working directories
current_wd = os.getcwd()
parent_wd = os.path.dirname(current_wd)
base_wd = os.path.dirname(os.path.dirname(parent_wd))
angola_wd = "/Users/thomasfuturewater/FutureWater Dropbox/Team/Projects/Completed/2019/2019019_G4AW_MavoDiami_Angola/Data/2019019_MavoDiami_LV/2019019_MavoDiami"

# Gets province from subprocess in 000_Run_All.py
PROVINCE_NAME = os.environ.get("PROVINCE")
if not os.environ.get("PROVINCE"):
    PROVINCE_NAME = "Sharkia"                            # dummy variable for testing.

# Define other folders
DATA_DIR = os.path.join(parent_wd, "01_Data")
GIS_DIR = os.path.join(base_wd, "GIS")     # Directory with shapefiles
RESULTS_DIR = os.path.join(parent_wd, "04_Results", PROVINCE_NAME)
os.makedirs(RESULTS_DIR, exist_ok=True)
TEMP_DIR = os.path.join(parent_wd, "05_Temp")

# Create output directory for DEM
DEM_PATH = os.path.join(DATA_DIR, "DEM", "DEM_NileDelta_250m.tif")
NEW_DEM_DIR = os.path.join(RESULTS_DIR, "DEM")
SLOPE_DIR = os.path.join(RESULTS_DIR, "_LS_RESULTS", "Slope")
os.makedirs(SLOPE_DIR, exist_ok=True)
os.makedirs(NEW_DEM_DIR, exist_ok=True)

# Temperature variables and parameters
RES = 250  # Resolution in meters
NO_DATA_VALUE = -9999.0  # No data value
LOCAL_PROJ = "EPSG:32636"  # Local projection
SRC_CRS = "EPSG:4326"  # Projection from NetCDF files: WGS84 projection
PARAMETERS = pd.read_csv(os.path.join(current_wd, "Parameters.csv"))

##############################################################################################
######################################### Mask DEM ###########################################
##############################################################################################
# Step 1: Crop and mask DEM with province shapefile.
print(f"Crop DEM: {PROVINCE_NAME}")
with rasterio.open(DEM_PATH) as src:
    # Get src data and profile
    src_nodata = src.nodata
    profile = src.profile.copy()
    data = src.read(1)
    dem_crs = src.crs

    # Load province shapefile and set to crs of DEM
    provinces_filepath = os.path.join(GIS_DIR, "Nile_delta_bnd_adm1.shp")
    provinces_shp = gpd.read_file(provinces_filepath)
    province_shp_sel = provinces_shp[provinces_shp["ADM1_EN"] == PROVINCE_NAME]
    province_shp_reproj = province_shp_sel.to_crs(LOCAL_PROJ)

    # Crop DEM to shapefile extent.
    out_image, out_transform = mask(
        src, province_shp_reproj.geometry, crop=True)

    # Copy metadata. Transform matrix is taken from DEM.
    out_profile = profile.copy()
    out_profile.update({
        "driver": "GTiff",
        "height": out_image.shape[1],
        "width": out_image.shape[2],
        # Transform: (pixel size x, row rotation, x-coordinate of upper-left corner, column rotation, pixel size y, y-coordinate of upper-left corner)
        "transform": out_transform,
        'nodata': src_nodata,
        'dtype': 'float32',
    })

# Save cropped DEM to temporary file
cropped_dem_path = os.path.join(
    TEMP_DIR, f"cropped_dem_{PROVINCE_NAME}_temp.tif")
with rasterio.open(cropped_dem_path, "w", **out_profile) as dest:
    dest.write(out_image)


##############################################################################################
############################# RESAMPLE ELEVATION AND SLOPE ###################################
##############################################################################################
# Step 2: Prepare target raster with desired resolution and projection for resampling.
print(f"Setting up target raster grid: {PROVINCE_NAME}")
# Determining bounds
bounds = province_shp_reproj.total_bounds  # [xmin, ymin, xmax, ymax]

# Calculate dimensions of target raster
width = int((bounds[2] - bounds[0]) / RES)
height = int((bounds[3] - bounds[1]) / RES)

# Create transform (local projection) for target raster
target_transform = rasterio.transform.from_bounds(
    bounds[0], bounds[1], bounds[2], bounds[3], width, height
)  # (left, bottom, right, top, width, height)

# Step 3: Reproject using bilinear method
print(f"Bilinear Resampling Raster DEM: {PROVINCE_NAME}")
bilinear_dem = np.zeros((height, width), dtype=np.float32)

with rasterio.open(cropped_dem_path) as src:
    reproject(
        source=src.read(1),
        destination=bilinear_dem,
        src_transform=src.transform,
        src_crs=src.crs,
        src_nodata=src.nodata,
        # Transform matrix from local projection Transform = (pixel size x, row rotation, x-coordinate of upper-left corner, column rotation, pixel size y, y-coordinate of upper-left corner)
        dst_transform=target_transform,
        dst_crs=LOCAL_PROJ,
        dst_nodata=NO_DATA_VALUE,
        resampling=Resampling.bilinear
    )

# Metadata for output files
out_profile = {
    "driver": "GTiff",
    "height": height,
    "width": width,
    "count": 1,
    "dtype": bilinear_dem.dtype,
    "crs": LOCAL_PROJ,
    "transform": target_transform,
    "nodata": NO_DATA_VALUE
}

# Create a memory file with the reprojected data for masking
with rasterio.MemoryFile() as memfile:
    with memfile.open(**out_profile) as temp_dst:
        temp_dst.write(bilinear_dem, 1)

        # Perform the masking operation
        masked_data, masked_transform = mask(
            temp_dst,
            province_shp_reproj.geometry,
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

# # Step 6: Write results to files
print(f"Write resampled rasters DEM: {PROVINCE_NAME}")

# Write bilinear resampled DEM
bilinear_path = os.path.join(NEW_DEM_DIR, f"DEM_{PROVINCE_NAME}_{RES}m.tif")
with rasterio.open(bilinear_path, "w", **masked_profile) as dst:
    dst.write(masked_data)

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
slope = calculate_slope(masked_data[0])

# Save slope raster
slope_path = os.path.join(SLOPE_DIR, f"Slope_{PROVINCE_NAME}.tif")
with rasterio.open(slope_path, "w", **out_profile) as dst:
    dst.write(slope.astype(rasterio.float32), 1)

# Step 8: Calculate areas with slope less than threshold
# Get slope threshold from params (you'll need to define this)
slope_threshold = float(
    PARAMETERS.loc[PARAMETERS['Parameter'] == 'Slope', 'Limit'].values[0])
slope_lower_limit = (slope < slope_threshold) & (masked_data[0] != NO_DATA_VALUE)
slope_lower_limit = slope_lower_limit.astype(np.uint8)  # Convert to uint8


# Save threshold raster
print(f"Write rasters slope and slope limit: {PROVINCE_NAME}")
lower_slope_path = os.path.join(
    SLOPE_DIR, f"Slope_lower_{slope_threshold}perc_{PROVINCE_NAME}.tif")
out_profile.update({
    'dtype': 'uint8',
    'nodata': 0,})

with rasterio.open(lower_slope_path, "w", **out_profile) as dst:
    dst.write(slope_lower_limit, 1)

# Remove temporary files
os.remove(cropped_dem_path)

print(f"Processing complete for {PROVINCE_NAME}\n")

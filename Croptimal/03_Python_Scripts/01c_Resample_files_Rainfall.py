import os
import glob
import re
import calendar
import numpy as np
import rasterio
from rasterio.warp import reproject, Resampling
from rasterio.mask import mask
import geopandas as gpd
import matplotlib.pyplot as plt


"""
Python conversion of the R script '01c_Resample_files_Rainfall.R'
This script processes monthly rainfall data by:
1. Reading monthly rainfall TIFs
2. Calculating monthly means across years
3. Cropping to a buffered province boundary
4. Resampling to match a reference DEM
5. Saving the processed data
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

# Set input data directories
DATA_DIR = os.path.join(angola_wd, "01_Data")
GIS_DIR = os.path.join(angola_wd, "02_GIS")     # Directory with shapefiles
RESULTS_DIR = os.path.join(parent_wd, "04_Results", PROVINCE_NAME)
TEMP_DIR = os.path.join(parent_wd, "05_Temp")
PROVINCES_FILEPATH = os.path.join(GIS_DIR, "Shapefiles", "AGO_adm1.shp")
RES = 250  # Set resolution

# Create output directory for processed rainfall data
RESULTS_RAINFALL_DIR = os.path.join(RESULTS_DIR, "Rainfall", "Mean_Monthly")
os.makedirs(RESULTS_RAINFALL_DIR, exist_ok=True)


####################################################################################################
######################### Group monthly rainfall files per month in dictionary #####################
####################################################################################################
# Function to extract month number from filename using regex
def extract_month(filename):
    """Extract month number (6-7th characters) from CHIRPS filename"""
    basename = os.path.basename(filename)
    # In R: substr(var_names_2, 6, 7)
    match = re.search(r'[^X\.]*\.?P_tot\.(\d{2})\.', basename)
    if match:
        return match.group(1)
    else:
        # Fallback: try to extract 6-7th characters directly
        try:
            # Remove 'X' and replace dots with underscores similar to R code
            clean_name = basename.replace('X', '').replace(
                '.', '_').replace('_P_tot', '')
            return clean_name[5:7]  # 6-7th characters (0-indexed)
        except IndexError:
            print(
                f"Warning: Could not extract month from filename: {basename}")
            return None


# Load all rainfall files
INPUT_FILES = glob.glob(os.path.join(
    DATA_DIR, "Rainfall", "Monthly_Sum", "*.tif"))
if not INPUT_FILES:
    raise FileNotFoundError(
        f"No rainfall files found in {os.path.join(DATA_DIR, 'Rainfall', 'Monthly_Sum')}")

print(f"Found {len(INPUT_FILES)} rainfall files")

# Group files by month (equivalent to stackApply in R)
MONTHLY_GROUPS = {}
for file in INPUT_FILES:
    month = extract_month(file)
    if month:
        if month not in MONTHLY_GROUPS:
            MONTHLY_GROUPS[month] = []
        MONTHLY_GROUPS[month].append(file)

####################################################################################################
################################# Calculate monthly mean rainfall #################################
####################################################################################################
# Calculate mean rainfall for each month
print("Calculating monthly mean rainfall...")
MEAN_RAINFALL_DATA = []
MEAN_RAINFALL_PROFILES = []

for month_idx, month in sorted([(int(k), k) for k in MONTHLY_GROUPS.keys()]):
    files = MONTHLY_GROUPS[month]
    print(f"  Processing month {month} ({len(files)} files)")

    # Read all files for this month
    month_data = []
    for file in files:
        with rasterio.open(file) as src:
            # Save profile from first file for this month
            if not month_data:
                profile = src.profile.copy()
            month_data.append(src.read(1))

    # Calculate mean for this month (equivalent to stackApply with mean function)
    if month_data:
        month_mean = np.nanmean(np.stack(month_data), axis=0)
        MEAN_RAINFALL_DATA.append(month_mean)
        MEAN_RAINFALL_PROFILES.append(profile)

    else:
        print(f"Warning: No data found for month {month}")

# Load reference DEM
print("Loading reference DEM...")
DEM_PATH = os.path.join(RESULTS_DIR, "DEM",
                        f"DEM_{PROVINCE_NAME}_{RES}m.tif")  # used to be: "DEM_{PROVINCE_NAME}_{RES}m_diff.tif". But currently not using the difference tif. Maybe needed in the future. Then also change it in 01b file.
with rasterio.open(DEM_PATH) as dem_src:
    DEM_META = dem_src.meta.copy()
    DEM_TRANSFORM = dem_src.transform
    DEM_CRS = dem_src.crs
    DEM_SHAPE = (dem_src.height, dem_src.width)

# Load province shapefile: (EPSG:4326)
print("Loading province shapefile...")
provinces_shp = gpd.read_file(PROVINCES_FILEPATH)
province_shp_sel = provinces_shp[provinces_shp["NAME_1"] == PROVINCE_NAME]

# Reproject to match DEM and have CRS with units in meters. CRS: (EPSG:32733)
province_shp_sel = province_shp_sel.to_crs(DEM_CRS)

# Create a buffer around province shapefile
print("Creating buffer around province...")
province_buffer = province_shp_sel.buffer(500)  # Buffer of 500

# List of month abbreviations (equivalent to month.abb in R)
MONTH_ABBRS = [calendar.month_abbr[i] for i in range(1, 13)]

# Process each month
for i, (month_data, month_profile) in enumerate(zip(MEAN_RAINFALL_DATA, MEAN_RAINFALL_PROFILES)):
    # Get month name and index
    month_idx = i + 1  # 1-based month index
    month_name = MONTH_ABBRS[i]
    print(f"Processing {month_name}...")

    # We need to create a temporary raster to perform masking
    # since we have the data as numpy arrays, not as raster files
    TEMP_RASTER_PATH = os.path.join(
        TEMP_DIR, f"temp_precipitation_{month_name}.tif")
    with rasterio.open(TEMP_RASTER_PATH, 'w', **month_profile) as temp:
        temp.write(month_data, 1)

    # Crop raster to the buffered province boundary
    print(f"  Cropping {month_name} to province boundary.")
    with rasterio.open(TEMP_RASTER_PATH) as src:
        # Ensure CRS province buffer and rainfall raster match: (EPSG:4326)
        if src.crs != province_buffer.crs:
            PROVINCE_BUFFER_PROJECTED = province_buffer.to_crs(src.crs)
        else:
            PROVINCE_BUFFER_PROJECTED = province_buffer

        # Crop raster with correct metadata
        out_image, out_transform = mask(
            src, PROVINCE_BUFFER_PROJECTED.geometry, crop=True)
        out_meta = src.meta.copy()
        out_meta.update({
            "driver": "GTiff",
            "height": out_image.shape[1],
            "width": out_image.shape[2],
            "transform": out_transform,
        })

    # Resample to match DEM resolution and extent
    print(f"  Resampling {month_name}...")
    resampled_data = np.zeros(DEM_SHAPE, dtype=month_data.dtype)

    # Create another temporary file for the cropped data
    TEMP_CROPPED_PATH = os.path.join(
        TEMP_DIR, f"cropped_temp_{month_name}.tif")
    with rasterio.open(TEMP_CROPPED_PATH, 'w', **out_meta) as temp:
        temp.write(out_image[0], 1)

    # Now resample from the cropped temp file to match DEM crs and resolution
    with rasterio.open(TEMP_CROPPED_PATH) as src:
        reproject(
            source=rasterio.band(src, 1),
            destination=resampled_data,
            src_transform=src.transform,
            src_crs=src.crs,
            dst_transform=DEM_TRANSFORM,
            dst_crs=DEM_CRS,
            resampling=Resampling.bilinear
        )

    # Save the final resampled file
    OUTPUT_PATH = os.path.join(RESULTS_RAINFALL_DIR, f"{month_name}.tif")
    out_meta = DEM_META.copy()
    with rasterio.open(OUTPUT_PATH, 'w', **out_meta) as dst:
        dst.write(resampled_data, 1)

    # Clean up temporary files
    try:
        os.remove(TEMP_RASTER_PATH)
        os.remove(TEMP_CROPPED_PATH)
    except Exception as e:
        print(f"  Warning: Could not remove temporary files: {e}")

print("Rainfall resampling complete!")

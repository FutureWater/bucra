import os
import glob
import re
import calendar
import numpy as np
import rasterio
from rasterio.warp import reproject, Resampling
from rasterio.mask import mask
import geopandas as gpd


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
DATA_DIR = "your_data_directory"  # Replace with actual path
RESULTS_DIR = "your_results_directory"  # Replace with actual path
PROVINCE_NAME = "your_province_name"  # Replace with actual province name
RES = 250  # Resolution in meters

# Create output directory for processed rainfall data
NEW_DIR = os.path.join(RESULTS_DIR, PROVINCE_NAME, "Rainfall", "Mean_Monthly")
os.makedirs(NEW_DIR, exist_ok=True)


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
            month_data.append(src.read(1))
            # Save profile from first file for this month
            if not month_data:
                profile = src.profile.copy()

    # Calculate mean for this month (equivalent to stackApply with mean function)
    if month_data:
        month_mean = np.nanmean(np.stack(month_data), axis=0)
        MEAN_RAINFALL_DATA.append(month_mean)
        MEAN_RAINFALL_PROFILES.append(profile)

# Load province shapefile
print("Loading province shapefile...")
PROVINCE_SHP_PATH = os.path.join(
    DATA_DIR, f"{PROVINCE_NAME}.shp")  # Adjust path as needed
PROVINCE_SHP = gpd.read_file(PROVINCE_SHP_PATH)

# Create a buffer around province shapefile (equivalent to gBuffer in R)
print("Creating buffer around province...")
PROVINCE_BUFFER = PROVINCE_SHP.buffer(0.25)  # Buffer of 0.25 degrees

# Load reference DEM
print("Loading reference DEM...")
DEM_PATH = os.path.join(RESULTS_DIR, PROVINCE_NAME, "DEM",
                        f"DEM_{PROVINCE_NAME}_{RES}m_diff.tif")
with rasterio.open(DEM_PATH) as dem_src:
    DEM_META = dem_src.meta.copy()
    DEM_TRANSFORM = dem_src.transform
    DEM_CRS = dem_src.crs
    DEM_SHAPE = (dem_src.height, dem_src.width)

# List of month abbreviations (equivalent to month.abb in R)
MONTH_ABBRS = [calendar.month_abbr[i] for i in range(1, 13)]

# Process each month
for i, (month_data, month_profile) in enumerate(zip(MEAN_RAINFALL_DATA, MEAN_RAINFALL_PROFILES)):
    month_idx = i + 1  # 1-based month index
    month_name = MONTH_ABBRS[i]
    print(f"Processing {month_name}...")

    # We need to create a temporary raster to perform masking
    # since we have the data as numpy arrays, not as raster files
    TEMP_RASTER_PATH = os.path.join(NEW_DIR, f"temp_{month_name}.tif")
    with rasterio.open(TEMP_RASTER_PATH, 'w', **month_profile) as temp:
        temp.write(month_data, 1)

    # Crop raster to the buffered province boundary
    print(f"  Cropping {month_name} to province boundary...")
    with rasterio.open(TEMP_RASTER_PATH) as src:
        # Ensure CRS compatibility
        if src.crs != PROVINCE_SHP.crs:
            PROVINCE_BUFFER_PROJECTED = PROVINCE_BUFFER.to_crs(src.crs)
        else:
            PROVINCE_BUFFER_PROJECTED = PROVINCE_BUFFER

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

    # Need another temporary file for the cropped data
    CROPPED_TEMP_PATH = os.path.join(NEW_DIR, f"cropped_temp_{month_name}.tif")
    with rasterio.open(CROPPED_TEMP_PATH, 'w', **out_meta) as temp:
        temp.write(out_image[0], 1)

    # Now reproject from the cropped temp file to match DEM
    with rasterio.open(CROPPED_TEMP_PATH) as src:
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
    OUTPUT_PATH = os.path.join(NEW_DIR, f"{month_name}.tif")
    out_meta = DEM_META.copy()
    with rasterio.open(OUTPUT_PATH, 'w', **out_meta) as dst:
        dst.write(resampled_data, 1)

    # Clean up temporary files
    try:
        os.remove(TEMP_RASTER_PATH)
        os.remove(CROPPED_TEMP_PATH)
    except Exception as e:
        print(f"  Warning: Could not remove temporary files: {e}")

print("Rainfall resampling complete!")

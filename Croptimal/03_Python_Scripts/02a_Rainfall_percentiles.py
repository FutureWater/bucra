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
# Define current and parent working directories
current_wd = os.getcwd()
parent_wd = os.path.dirname(current_wd)
angola_wd = "/Users/thomasfuturewater/FutureWater Dropbox/Team/Projects/Completed/2019/2019019_G4AW_MavoDiami_Angola/Data/2019019_MavoDiami_LV/2019019_MavoDiami"

# Gets province from subprocess in 000_Run_All.py
PROVINCE_NAME = os.environ.get("PROVINCE")
PROVINCE_NAME = "Zaire"                         # dummy variable for testing.

# Define other folders
DATA_DIR = os.path.join(angola_wd, "01_Data")
GIS_DIR = os.path.join(angola_wd, "02_GIS")     # Directory with shapefiles
RESULTS_DIR = os.path.join(parent_wd, "04_Results", PROVINCE_NAME)
TEMP_DIR = os.path.join(parent_wd, "05_Temp")

# Define constants
RES = 250  # Resolution in meters
P_PERC = [0.05, 0.25]  # Percentiles to calculate (5% and 25%)
NO_DATA_VALUE = -9999.0  # No data value

# Get input rainfall files
input_files = glob.glob(os.path.join(
    DATA_DIR, "Rainfall", "Monthly_Sum", "*.tif"))
if not input_files:
    raise FileNotFoundError(
        f"No rainfall files found in {os.path.join(DATA_DIR, 'Rainfall', 'Monthly_Sum')}")

print(f"Found {len(input_files)} rainfall files")

# Extract month numbers from filenames
months_stack = [int(os.path.basename(file)[5:7]) for file in input_files]


# Load reference DEM
print("Loading reference DEM...")
dem_path = os.path.join(
    RESULTS_DIR, "DEM", f"DEM_{PROVINCE_NAME}_{RES}m.tif")  # used to be: "DEM_{PROVINCE_NAME}_{RES}m_diff.tif"
with rasterio.open(dem_path) as dem_src:
    DEM_META = dem_src.meta.copy()
    DEM_TRANSFORM = dem_src.transform
    DEM_CRS = dem_src.crs
    DEM_SHAPE = (dem_src.height, dem_src.width)
    DEM_PROFILE = dem_src.profile

# Load province shapefile
print("Loading province shapefile...")
provinces_filepath = os.path.join(GIS_DIR, "Shapefiles", "AGO_adm1.shp")
provinces_shp = gpd.read_file(provinces_filepath)
province_shp_sel = provinces_shp[provinces_shp["NAME_1"] == PROVINCE_NAME]
province_shp_sel_reproj = province_shp_sel.to_crs(DEM_CRS)
province_buffer = province_shp_sel.buffer(500)  # Buffer of 0.25 degrees

####################################################################################################
########################## Process each percentile ################################################
####################################################################################################
# Process each percentile
for p_p in P_PERC:
    # Generate directory name by removing decimal point from percentile value
    p_str = str(p_p).replace(".", "")
    new_results_subdir = os.path.join(RESULTS_DIR,
                                      "Rainfall", "Percentiles", f"{p_str}perc")
    os.makedirs(new_results_subdir, exist_ok=True)

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

    # Month abbreviations for naming output files
    month_abbrs = [calendar.month_abbr[i] for i in range(1, 13)]

    # Process each month
    for month, data in out_layers.items():
        # Create a temporary raster file for cropping
        temp_path = os.path.join(TEMP_DIR, f"temp_rainfall_sum_{month}.tif")

        # Write the data to the temporary file
        with rasterio.open(temp_path, 'w', **temp_profile) as dst:
            dst.write(data.astype(rasterio.float32), 1)

        # Open temporary file for resampling and cropping
        with rasterio.open(temp_path) as src:
            temp_data = src.read(1)
            temp_profile = src.profile

            # Create an output array for the resampled data
            resampled_data = np.full(
                DEM_SHAPE, NO_DATA_VALUE, dtype=np.float32)

            # Reproject and resample the data to match the DEM
            reproject(
                source=temp_data,
                destination=resampled_data,
                src_transform=src.transform,
                src_crs=src.crs,  # Assuming all monthly rasters have the same CRS
                dst_transform=DEM_TRANSFORM,
                dst_crs=DEM_CRS,
                resampling=Resampling.bilinear
            )

            # Create a memory file with the reprojected data for masking
            with rasterio.MemoryFile() as memfile:
                with memfile.open(**DEM_PROFILE) as temp_dst:
                    temp_dst.write(resampled_data, 1)

                    # Ensure CRS compatibility
                    province_shp_proj = province_shp_sel_reproj.to_crs(DEM_CRS)

                    # Mask the data to the province boundary
                    masked_data, masked_transform = mask(
                        temp_dst,
                        province_shp_proj.geometry,
                        crop=True,
                        nodata=NO_DATA_VALUE)

                    # Save mask profile
                    masked_profile = temp_dst.profile.copy()
                    masked_profile.update({
                        "height": masked_data.shape[1],
                        "width": masked_data.shape[2],
                        "transform": masked_transform,
                    })

            # Save the masked data to a new file
            output_path = os.path.join(
                new_results_subdir, f"rainfall_{p_str}perc_{month_abbrs[month-1]}.tif")
            with rasterio.open(output_path, 'w', **masked_profile) as dst:
                dst.write(masked_data.astype(rasterio.float32))

            print(f"  Saved {month_abbrs[month-1]} for {p_p*100}% percentile")

            # Clean up temporary file
            os.remove(temp_path)

print("Rainfall percentiles processing complete!")

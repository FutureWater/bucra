import os
import glob
import calendar
import pandas as pd
import rasterio
import numpy as np
import matplotlib.pyplot as plt

"""
This script calculates suitability limits for various parameters:
1) Soil hydraulic properties (Ksat, WCavail)
2) Soil nutrient content
3) NDVI during growing seasons
"""


####################################################################################################
###################### Define directories, constants and file paths ################################
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
HHS_DATA_DIR = os.path.join(DATA_DIR, "Soil_Hydraulic_Properties")
LS_RESULTS = os.path.join(RESULTS_DIR, '_LS_Results')

# Define constants
RES = 250  # Resolution in meters
NO_DATA_VALUE = -9999.0  # No data value
LOCAL_PROJECTION = "EPSG:32733"


# Load cropping calendar and parameters
CROPPING_CAL = pd.read_csv(os.path.join(current_wd, "Cropping_calendar_A.csv"))
PARAMS = pd.read_csv(os.path.join(current_wd, "Parameters.csv"))


####################################################################################################
###################### Process parameters and create limit maps ####################################
####################################################################################################
############## Process each parameter folder ##############
for folder in PARAMS['LS_Results_Folder'].unique():
    "Processing other suitability parameters..."
    # Get parameters for this folder
    folder_mask = PARAMS['LS_Results_Folder'] == folder
    parameters = PARAMS.loc[folder_mask, 'Parameter'].tolist()
    limits = PARAMS.loc[folder_mask, 'Limit'].tolist()
    units = PARAMS.loc[folder_mask, 'Units'].tolist()
    abrevs = PARAMS.loc[folder_mask, 'Abrev'].tolist()

    # Create output directory
    new_results_subdir = os.path.join(RESULTS_DIR, "_LS_Results", folder)
    os.makedirs(new_results_subdir, exist_ok=True)

    # Skip folders handled in other scripts
    if folder in ["Elevation", "Water"]:
        print(f"    Skipping {folder} - calculated in other scripts")
        continue

    ############## Process NDVI ##############
    if folder == "NDVI":
        print(f"    Processing {folder}...")
        # Get NDVI files
        ndvi_files = glob.glob(os.path.join(
            RESULTS_DIR, "NDVI", 'Mean_Monthly', "*.tif"))

        # Process each crop's growing season
        for _, crop_row in CROPPING_CAL.iterrows():
            start_month = crop_row['Start_growing_season']
            end_month = crop_row['End_growing_season']
            crop = crop_row['Crop']

            # Determine months in growing season
            if start_month > end_month:
                months = list(range(start_month, 13)) + \
                    list(range(1, end_month + 1))
            else:
                months = list(range(start_month, end_month + 1))

            month_abbrs = [calendar.month_abbr[m] for m in months]

            # Find NDVI files for these months
            season_ndvi_files = []
            for month_abbr in month_abbrs:
                season_ndvi_files.extend(
                    [f for f in ndvi_files if month_abbr in os.path.basename(f)])

            if not season_ndvi_files:
                print(f"  No NDVI files found for months {month_abbrs}")
                continue

            # Sum all NDVI for the whole season
            ndvi_sum = None
            for ndvi_file in season_ndvi_files:
                with rasterio.open(ndvi_file) as src:
                    ndvi_data = src.read(1)
                    if ndvi_sum is None:
                        # Initialize arrays
                        ndvi_sum = np.zeros(ndvi_data.shape, dtype=np.float32)
                        ndvi_valid = np.zeros(ndvi_data.shape, dtype=bool)
                        ndvi_meta = src.meta.copy()

                    # Only add valid cells. Add zero for non valid cells. Update valid mask.
                    mask = ndvi_data != NO_DATA_VALUE
                    ndvi_data[~mask] = 0
                    ndvi_sum += ndvi_data
                    ndvi_valid |= mask

            # Calculate mean NDVI for the whole season
            ndvi_mean = np.full(
                ndvi_sum.shape, NO_DATA_VALUE, dtype=np.float32)
            ndvi_mean[ndvi_valid] = ndvi_sum[ndvi_valid] / \
                len(season_ndvi_files)

            # Apply limit
            # Assuming NDVI is the last parameter
            limit_value = float(limits[0])
            ndvi_limit = (ndvi_mean > limit_value) & ndvi_valid
            ndvi_limit = ndvi_limit.astype(np.uint8)

            # Save output
            output_path = os.path.join(
                new_results_subdir,
                f"NDVI_limit_higher_{limits[0].replace('.', '')}_{crop}_"
                f"{calendar.month_abbr[start_month]}-{calendar.month_abbr[end_month]}.tif"
            )

            out_meta = ndvi_meta.copy()
            out_meta.update({'dtype': 'uint8',
                             'nodata': 0})

            with rasterio.open(output_path, 'w', **out_meta) as dst:
                dst.write(ndvi_limit, 1)

    ############## Process soil nutrient content ##############
    elif folder == "Soil_Nutrient_Content":
        print(f"Processing {folder}...")

        # Get files
        snc_files = glob.glob(os.path.join(
            RESULTS_DIR, "Soil_Nutrient_Content", "*.tif"))

        # Process each parameter
        for i, param in enumerate(parameters[:1]):
            abrev = abrevs[i]
            limit = limits[i]
            unit = units[i]

            # Find matching file
            matching_files = [
                f for f in snc_files if f"_{abrev}_" in os.path.basename(f).lower()]
            if not matching_files:
                print(f"  No data found for {param}")
                continue

            # Open file and apply limit
            with rasterio.open(matching_files[0]) as src:
                snc_data = src.read(1)
                snc_meta = src.meta.copy()

                # Apply limit
                limit_value = float(limit)
                limit_raster = (snc_data > limit_value).astype(np.uint8)

                # Save output
                output_path = os.path.join(
                    new_results_subdir,
                    f"Extractable_{abrev.upper()}_Soil_higher_than_{limit}{unit}_{PROVINCE_NAME}.tif"
                )

                out_meta = snc_meta.copy()
                out_meta.update({'dtype': 'uint8',
                                'nodata': 0})

                with rasterio.open(output_path, 'w', **out_meta) as dst:
                    dst.write(limit_raster, 1)

    # ############## Process soil hydraulic properties ##############
    # if folder == "Soil_Hydraulic_Properties":
    #   print(f"Processing {folder}...")
    #     # Load reference DEM
    #     dem_path = os.path.join(
    #         RESULTS_DIR, "DEM", f"DEM_{PROVINCE_NAME}_{RES}m.tif")
    #     with rasterio.open(dem_path) as dem_src:
    #         dem_profile = dem_src.profile.copy()
    #         dem_crs = dem_src.crs
    #         dem_transform = dem_src.transform
    #         dem_shape = dem_src.shape

    #     # Create boundary in WGS84 for cropping
    #     boundary_crs = 'EPSG:4326'
    #     boundary_transform, boundary_width, boundary_height = calculate_default_transform(
    #         dem_crs, boundary_crs, dem_shape[1], dem_shape[0],
    #         left = dem_src.bounds.left, bottom = dem_src.bounds.bottom,
    #         right = dem_src.bounds.right, top = dem_src.bounds.top,
    #         resolution=RES
    #     )

    #     # Load HHS data files that match parameters
    #     hhs_files = []
    #     for param in parameters:
    #         hhs_files.extend(glob.glob(os.path.join(
    #             HHS_DATA_DIR, "*.tif")))

    #     # Load and process each HHS parameter
    #     for param in parameters:
    #         print(f"  Processing {param}...")

    #         # Set multiplier based on parameter
    #         MULTIPLIER = 10 if param == "Ksat" else 1000  # mm/d for Ksat, mm/m for WCavail

    #         # Get limit and unit for this parameter
    #         idx = parameters.index(param)
    #         limit = limits[idx]
    #         unit = units[idx]

    #         # Find topsoil and subsoil files
    #         topsoil_file = None
    #         subsoil_file = None
    #         for file in hhs_files:
    #             if param in file.lower():
    #                 if "topsoil" in file.lower():
    #                     topsoil_file = file
    #                 elif "subsoil" in file.lower():
    #                     subsoil_file = file

    #         if not topsoil_file or not subsoil_file:
    #             print(f"  Missing soil data for {param}")
    #             continue

    #         # Process topsoil
    #         with rasterio.open(topsoil_file) as src:
    #             # Division by 10000 as in R script
    #             topsoil_data = src.read(1) / 10000
    #             topsoil_meta = src.meta.copy()

    #             # Reproject to match DEM
    #             topsoil_reproj = np.zeros(dem_shape, dtype=np.float32)
    #             reproject(
    #                 source=rasterio.band(src, 1),
    #                 destination=topsoil_reproj,
    #                 src_transform=src.transform,
    #                 src_crs=src.crs,
    #                 dst_transform=dem_transform,
    #                 dst_crs=dem_crs,
    #                 resampling=Resampling.bilinear,
    #                 src_nodata=NO_DATA_VALUE,
    #                 dst_nodata=NO_DATA_VALUE
    #             )

    #         # Process subsoil
    #         with rasterio.open(subsoil_file) as src:
    #             # Division by 10000 as in R script
    #             subsoil_data = src.read(1) / 10000

    #             # Reproject to match DEM
    #             subsoil_reproj = np.zeros(dem_shape, dtype=np.float32)
    #             reproject(
    #                 source=rasterio.band(src, 1),
    #                 destination=subsoil_reproj,
    #                 src_transform=src.transform,
    #                 src_crs=src.crs,
    #                 dst_transform=dem_transform,
    #                 dst_crs=dem_crs,
    #                 resampling=Resampling.bilinear,
    #                 src_nodata=NO_DATA_VALUE,
    #                 dst_nodata=NO_DATA_VALUE
    #             )

    #         # Calculate weighted average (0.3 * topsoil + 1.7 * subsoil)/2
    #         weighted_topsoil = topsoil_reproj * 0.3
    #         weighted_subsoil = subsoil_reproj * 1.7
    #         weighted_avg = (weighted_topsoil +
    #                         weighted_subsoil) / 2 * MULTIPLIER

    #         # Apply limit
    #         limit_value = float(limit)
    #         limit_raster = (weighted_avg > limit_value).astype(np.uint8)

    #         # Save output
    #         output_path = os.path.join(
    #             new_dir, f"{param}_Soil_higher_{limit}{unit}_{PROVINCE_NAME}.tif")
    #         out_meta = dem_profile.copy()
    #         out_meta.update({
    #             'dtype': 'uint8',
    #             'count': 1
    #         })

    #         with rasterio.open(output_path, 'w', **out_meta) as dst:
    #             dst.write(limit_raster, 1)

print("Parameter limits calculation complete!")

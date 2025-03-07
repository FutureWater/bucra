import os
import glob
import calendar
import pandas as pd
import rasterio
import numpy as np
from rasterio.warp import calculate_default_transform, reproject, Resampling
from rasterio.crs import CRS
from rasterio.features import geometry_mask

"""
This script calculates suitability limits for various parameters:
1) Soil hydraulic properties (Ksat, WCavail)
2) Soil nutrient content
3) NDVI during growing seasons
"""

# Define constants
DATA_DIR = "your_data_directory"  # Replace with actual path
RESULTS_DIR = "your_results_directory"  # Replace with actual path
PROVINCE_NAME = "your_province_name"  # Replace with actual province name
RES = 250  # Resolution in meters
HHS_DATA_DIR = "path_to_hhs_data"  # Replace with actual HHS data directory

# Load parameters table (replace with actual data loading)
# Structure needs: Parameter, LS_Results_Folder, Weight, Limit, Units, Abrev
params = pd.DataFrame({
    'Parameter': ['Ksat', 'WCavail', 'P', 'K', 'N', 'NDVI'],
    'LS_Results_Folder': ['Soil_Hydraulic_Properties', 'Soil_Hydraulic_Properties',
                          'Soil_Nutrient_Content', 'Soil_Nutrient_Content',
                          'Soil_Nutrient_Content', 'NDVI'],
    'Weight': [1, 1, 1, 1, 1, 1],
    'Limit': ['10', '100', '5', '100', '10', '0.6'],
    'Units': ['mm/d', 'mm/m', 'mg/kg', 'mg/kg', 'mg/kg', ''],
    'Abrev': ['ksat', 'wcavail', 'p', 'k', 'n', 'ndvi']
})

# Load crop calendar (replace with actual data loading)
cropping_cal = pd.DataFrame({
    'Crop': ['Wheat'],
    'Start_growing_season': [11],
    'End_growing_season': [4]
})

# Process each parameter folder
for folder in params['LS_Results_Folder'].unique():
    # Get parameters for this folder
    folder_mask = params['LS_Results_Folder'] == folder
    parameters = params.loc[folder_mask, 'Parameter'].tolist()
    limits = params.loc[folder_mask, 'Limit'].tolist()
    units = params.loc[folder_mask, 'Units'].tolist()
    abrevs = params.loc[folder_mask, 'Abrev'].tolist()

    # Create output directory
    new_dir = os.path.join(RESULTS_DIR, PROVINCE_NAME, "_LS_Results", folder)
    os.makedirs(new_dir, exist_ok=True)

    # Skip folders handled in other scripts
    if folder in ["Elevation", "Water"]:
        print(f"Skipping {folder} - calculated in other scripts")
        continue

    # Process soil hydraulic properties
    if folder == "Soil_Hydraulic_Properties":
        print(f"Processing {folder}...")

        # Load reference DEM
        dem_path = os.path.join(
            RESULTS_DIR, PROVINCE_NAME, "DEM", f"DEM_{PROVINCE_NAME}_{RES}m_diff.tif")
        with rasterio.open(dem_path) as dem_src:
            dem_meta = dem_src.meta.copy()
            dem_crs = dem_src.crs
            dem_transform = dem_src.transform
            dem_shape = dem_src.shape

        # Create boundary in WGS84 for cropping
        boundary_crs = CRS.from_epsg(4326)
        boundary_transform, boundary_width, boundary_height = calculate_default_transform(
            dem_crs, boundary_crs, dem_shape[1], dem_shape[0],
            transform=dem_transform
        )

        # Load HHS data files that match parameters
        hhs_files = []
        for param in parameters:
            hhs_files.extend(glob.glob(os.path.join(
                HHS_DATA_DIR, f"**/*{param}*.tif"), recursive=True))

        # Load and process each HHS parameter
        for param in parameters:
            print(f"  Processing {param}...")

            # Set multiplier based on parameter
            multiplier = 10 if param == "Ksat" else 1000  # mm/d for Ksat, mm/m for WCavail

            # Get limit and unit for this parameter
            idx = parameters.index(param)
            limit = limits[idx]
            unit = units[idx]

            # Find topsoil and subsoil files
            topsoil_file = None
            subsoil_file = None
            for file in hhs_files:
                if param in file.lower():
                    if "topsoil" in file.lower():
                        topsoil_file = file
                    elif "subsoil" in file.lower():
                        subsoil_file = file

            if not topsoil_file or not subsoil_file:
                print(f"  Missing soil data for {param}")
                continue

            # Process topsoil
            with rasterio.open(topsoil_file) as src:
                # Division by 10000 as in R script
                topsoil_data = src.read(1) / 10000
                topsoil_meta = src.meta.copy()

                # Reproject to match DEM
                topsoil_reproj = np.zeros(dem_shape, dtype=np.float32)
                reproject(
                    source=rasterio.band(src, 1),
                    destination=topsoil_reproj,
                    src_transform=src.transform,
                    src_crs=src.crs,
                    dst_transform=dem_transform,
                    dst_crs=dem_crs,
                    resampling=Resampling.bilinear
                )

            # Process subsoil
            with rasterio.open(subsoil_file) as src:
                # Division by 10000 as in R script
                subsoil_data = src.read(1) / 10000

                # Reproject to match DEM
                subsoil_reproj = np.zeros(dem_shape, dtype=np.float32)
                reproject(
                    source=rasterio.band(src, 1),
                    destination=subsoil_reproj,
                    src_transform=src.transform,
                    src_crs=src.crs,
                    dst_transform=dem_transform,
                    dst_crs=dem_crs,
                    resampling=Resampling.bilinear
                )

            # Calculate weighted average (0.3 * topsoil + 1.7 * subsoil)/2
            weighted_topsoil = topsoil_reproj * 0.3
            weighted_subsoil = subsoil_reproj * 1.7
            weighted_avg = (weighted_topsoil +
                            weighted_subsoil) / 2 * multiplier

            # Apply limit
            limit_value = float(limit)
            limit_raster = (weighted_avg > limit_value).astype(np.uint8)

            # Save output
            output_path = os.path.join(
                new_dir, f"{param}_Soil_higher_{limit}{unit}_{PROVINCE_NAME}.tif")
            out_meta = dem_meta.copy()
            out_meta.update({
                'dtype': 'uint8',
                'count': 1
            })

            with rasterio.open(output_path, 'w', **out_meta) as dst:
                dst.write(limit_raster, 1)

    # Process soil nutrient content
    elif folder == "Soil_Nutrient_Content":
        print(f"Processing {folder}...")

        # Get files
        snc_files = glob.glob(os.path.join(
            RESULTS_DIR, PROVINCE_NAME, "Soil_Nutrient_Content", "*.tif"))

        # Process each parameter
        for i, param in enumerate(parameters):
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
                    new_dir,
                    f"Extractable_{abrev.upper()}_Soil_higher_{limit}{unit}_{PROVINCE_NAME}.tif"
                )

                out_meta = snc_meta.copy()
                out_meta.update({'dtype': 'uint8'})

                with rasterio.open(output_path, 'w', **out_meta) as dst:
                    dst.write(limit_raster, 1)

    # Process NDVI
    elif folder == "NDVI":
        print(f"Processing {folder}...")

        # Get NDVI files
        ndvi_files = glob.glob(os.path.join(
            RESULTS_DIR, PROVINCE_NAME, "NDVI", "*.tif"))

        # Process each crop's growing season
        for _, crop_row in cropping_cal.iterrows():
            start_month = crop_row['Start_growing_season']
            end_month = crop_row['End_growing_season']

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

            # Calculate mean NDVI
            ndvi_sum = None
            for ndvi_file in season_ndvi_files:
                with rasterio.open(ndvi_file) as src:
                    ndvi_data = src.read(1)
                    if ndvi_sum is None:
                        ndvi_sum = ndvi_data
                        ndvi_meta = src.meta.copy()
                    else:
                        ndvi_sum += ndvi_data

            # Calculate mean
            ndvi_mean = ndvi_sum / len(season_ndvi_files)

            # Apply limit
            # Assuming NDVI is the last parameter
            limit_value = float(limits[-1])
            ndvi_limit = (ndvi_mean > limit_value).astype(np.uint8)

            # Save output
            output_path = os.path.join(
                new_dir,
                f"NDVI_limit_higher_{limits[-1].replace('.', '_')}_"
                f"{calendar.month_abbr[start_month]}-{calendar.month_abbr[end_month]}_{PROVINCE_NAME}.tif"
            )

            out_meta = ndvi_meta.copy()
            out_meta.update({'dtype': 'uint8'})

            with rasterio.open(output_path, 'w', **out_meta) as dst:
                dst.write(ndvi_limit, 1)

print("Parameter limits calculation complete!")

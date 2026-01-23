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
###################### Define fuzzy membership functions ############################################
####################################################################################################
def fuzzy_membership(raster, function_type, parameters):
    result = raster.copy().astype(float)

    if function_type == "Increasing":
        # Two parameters: [a, b]
        # below a = 0, between a-b = 0 to 1, above b = 1
        a = parameters[0]
        b = parameters[1]
        
        result[raster <= a] = 0
        result[raster >= b] = 1
        mask = (a < raster) & (raster < b)
        result[mask] = (raster[mask] - a) / (b-a)
        
    elif function_type == "Decreasing":
        # Two parameters: [a, b]
        # below a = 1, between a-b = 1 to 0, above b = 0
        a = parameters[0]
        b = parameters[1]
        
        result[raster <= a] = 1
        result[raster >= b] = 0
        mask = (a < raster) & (raster < b)
        result[mask] = 1 - (raster[mask] - a) / (b - a)
    
    elif function_type == "Triangular":
        # Three parameters: [a, b, c]
        # below a = 0, a-b = 0 to 1, b-c = 1 to 0, above c = 0
        a = parameters[0]
        b = parameters[1]  # peak point
        c = parameters[2]

        result[raster <= a] = 0
        result[raster >= c] = 0

        mask1 = (a < raster) & (raster < b)
        mask2 = (b < raster) & (raster < c)
        result[mask1] = (raster[mask1] - a) / (b - a)
        result[mask2] = 1 - (raster[mask2] - b) / (c - b)
    
    elif function_type == "Trapezoidal":
        # Four parameters: [a, b, c, d]
        # below a = 0, a-b = 0 to 1, b-c = 1, c-d = 1 to 0, above d = 0
        a = parameters[0]
        b = parameters[1]  # start of plateau
        c = parameters[2]  # end of plateau
        d = parameters[3]

        result[raster <= a] = 0
        result[raster >= d] = 0

        mask1 = (a < raster) & (raster < b)
        mask2 = (b < raster) & (raster < c)
        mask3 = (c < raster) & (raster <d)
        result[mask1] = (raster[mask1] - a) / (b - a)
        result[mask2] = 1
        result[mask3] = 1 - (raster[mask3] -c) / (d - c)
    
    if np.all(result == raster):
        raise ValueError("Result raster is the same as input raster. Fuzzy logic not applied")
    
    return result

def plot_raster(raster, ):
    plt.imshow(raster)
    plt.colorbar()
####################################################################################################
###################### Define directories, constants and file paths ################################
####################################################################################################
# Define current and parent working directories
current_wd = os.getcwd()
parent_wd = os.path.dirname(current_wd)
base_wd = os.path.dirname(os.path.dirname(parent_wd))
angola_wd = "/Users/thomasfuturewater/FutureWater Dropbox/Team/Projects/Completed/2019/2019019_G4AW_MavoDiami_Angola/Data/2019019_MavoDiami_LV/2019019_MavoDiami"

# Gets province from subprocess in 000_Run_All.py
PROVINCE_NAME = os.environ.get("PROVINCE")
if not os.environ.get("PROVINCE"):
    PROVINCE_NAME = "Sharkia"                      # dummy variable for testing.

# Define other folders
DATA_DIR = os.path.join(parent_wd, "01_Data")
GIS_DIR = os.path.join(base_wd, "GIS")     # Directory with shapefiles
RESULTS_DIR = os.path.join(parent_wd, "04_Results", PROVINCE_NAME)
os.makedirs(RESULTS_DIR, exist_ok=True)
TEMP_DIR = os.path.join(parent_wd, "05_Temp")
HHS_DATA_DIR = os.path.join(RESULTS_DIR, "Soil_Hydraulic_Properties")
LS_RESULTS = os.path.join(RESULTS_DIR, '_LS_Results')

# Define constants
RES = 250  # Resolution in meters
NO_DATA_VALUE = -9999.0  # No data value
LOCAL_PROJECTION = "EPSG:32636"


# Load cropping calendar and parameters
CROPPING_CAL = pd.read_csv(os.path.join(current_wd, "Cropping_calendar.csv"), sep = ',')
PARAMS = pd.read_csv(os.path.join(current_wd, "Parameters_fuzzy.csv"))


####################################################################################################
###################### Process parameters and create limit maps ####################################
####################################################################################################
############## Process each parameter folder ##############
print("Processing remaining suitability parameters...")
for folder in PARAMS['LS_Results_Folder'].unique():
    # Get parameters for this folder
    folder_mask = PARAMS['LS_Results_Folder'] == folder
    parameters = PARAMS.loc[folder_mask, 'Parameter'].tolist()
    limits = PARAMS.loc[folder_mask, ['Limit1', 'Limit2','Limit3','Limit4']].values.flatten().tolist()
    limits = [float(limit) for limit in limits]
    units = PARAMS.loc[folder_mask, 'Units'].tolist()
    abrevs = PARAMS.loc[folder_mask, 'Abrev'].tolist()
    fuzzy_functions = PARAMS.loc[folder_mask, 'Function_type'].tolist()

    # Create output directory
    new_results_subdir = os.path.join(RESULTS_DIR, "_LS_Results", folder)
    os.makedirs(new_results_subdir, exist_ok=True)

    # Skip folders handled in other scripts
    if folder in ["Elevation", "Water", "Temperature"]:
        continue

    ############## Process NDVI ##############
    if folder == "NDVI":
        print(f"    Processing {folder}...")
        # Get NDVI files
        ndvi_files = glob.glob(os.path.join(
            RESULTS_DIR, "NDVI", 'Mean_Monthly', "*.tif"))

        # Find NDVI files for the whole year
        ndvi_files = glob.glob(os.path.join(RESULTS_DIR, "NDVI", 'Mean_Monthly', "*.tif"))

        if not ndvi_files:
            print(f"  No NDVI files found")
            continue

        # Sum all NDVI
        all_ndvi = None
        for i, ndvi_file in enumerate(ndvi_files):
            with rasterio.open(ndvi_file) as src:
                ndvi_data = src.read(1)
                
                # Initialize arrays
                if all_ndvi is None:   
                    all_ndvi = np.zeros((len(ndvi_files), ndvi_data.shape[0], ndvi_data.shape[1]), dtype=np.float32)
                    ndvi_meta = src.meta.copy()

                # Add NDVI data to 3d array
                all_ndvi[i] = ndvi_data

        # Calculate max NDVI per cell for the whole year.
        max_ndvi = np.nanmax(all_ndvi, axis=0)

        # Apply fuzzy logic
        ndvi_limit = fuzzy_membership(max_ndvi, fuzzy_functions[0], limits)
        ndvi_limit = ndvi_limit.astype(np.float32)

        # Save output
        output_path = os.path.join(
            new_results_subdir,
            f"NDVI_suitability_score.tif"
        )

        out_profile = ndvi_meta.copy()
        out_profile.update({'dtype': 'float32',
                            'nodata': 0})

        with rasterio.open(output_path, 'w', **out_profile) as dst:
            dst.write(ndvi_limit, 1)

    ############## Process soil nutrient content ##############
    elif folder == "Soil_Nutrient_Content":
        print(f"    Processing {folder}...")

        # Get files
        snc_files = glob.glob(os.path.join(
            RESULTS_DIR, "Soil_Nutrient_Content", "*.tif"))

        limits = [limits[:4], limits[4:]]
        # Process each parameter
        for i, param in enumerate(parameters):
            abrev = abrevs[i]
            limit = limits[i]
            unit = units[i]
            parameter = parameters[i]
            fuzzy_function = fuzzy_functions[i]

            # Find matching file
            matching_files = [
                f for f in snc_files if parameter.lower() in os.path.basename(f).lower()]
            if not matching_files:
                print(f"  No data found for {param}")
                continue

            # Open file and apply limit
            with rasterio.open(matching_files[0]) as src:
                snc_data = src.read(1)
                snc_meta = src.meta.copy()

                # Apply limit
                snc_limit = fuzzy_membership(snc_data, fuzzy_function, limit)
                snc_limit = snc_limit.astype(np.float32)
                
                # Save output
                output_path = os.path.join(
                    new_results_subdir,
                    f"Extractable_{abrev.upper()}_Soil_suitability_score.tif"
                )

                out_profile = snc_meta.copy()
                out_profile.update({'dtype': 'uint8',
                                'nodata': 0})

                with rasterio.open(output_path, 'w', **out_profile) as dst:
                    dst.write(snc_limit, 1)

    # ############## Process soil hydraulic properties ##############
    # if folder == "Soil_Hydraulic_Properties":
    #     print(f"Processing {folder}...")
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

    #         # Find topsoil and subsoil files and read data
    #         topsoil_file = None
    #         subsoil_file = None
    #         for file in hhs_files:
    #             if param in file:
    #                 # Open topsoil and subsoil files
    #                 if "TOPSOIL" in file:
    #                     topsoil_file = file
    #                     with rasterio.open(topsoil_file) as src:
    #                         topsoil_data = src.read(1)
    #                         topsoil_profile = src.profile.copy()
    #                 elif "SUBSOIL" in file:
    #                     subsoil_file = file
    #                     with rasterio.open(subsoil_file) as src:
    #                         subsoil_data = src.read(1)

    #         if not topsoil_file or not subsoil_file:
    #             print(f"  Missing soil data for {param}")
    #             continue

    #         # Calculate weighted average (0.3 * topsoil + 1.7 * subsoil)/2
    #         weighted_topsoil = topsoil_data * 0.3
    #         weighted_subsoil = subsoil_data * 1.7
    #         weighted_avg = (weighted_topsoil +
    #                         weighted_subsoil) / 2 * MULTIPLIER

    #         # Apply limit
    #         limit_value = float(limit)
    #         limit_raster = (weighted_avg > limit_value).astype(np.uint8)

    #         # Save output
    #         output_path = os.path.join(
    #             new_results_subdir, f"{param}_Soil_higher_{limit}{unit}_{PROVINCE_NAME}.tif")
    #         out_profile = topsoil_profile.copy()
    #         out_profile.update({
    #             'dtype': 'uint8',
    #             'count': 1
    #         })

    #         with rasterio.open(output_path, 'w', **out_profile) as dst:
    #             dst.write(limit_raster, 1)

print("Parameter limits calculation complete!")

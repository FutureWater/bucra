# Import necessary libraries
import os
import glob
import datetime
import calendar
import subprocess
from pathlib import Path

import pandas as pd
import numpy as np
import geopandas as gpd
import rasterio
from rasterio.env import Env
import matplotlib.pyplot as plt


##############################################################################################
################################### START OF DATA INPUT ######################################
##############################################################################################
# Get current working directory.
current_wd = os.getcwd()
parent_wd = os.path.dirname(current_wd)
angola_wd = "/Users/thomasfuturewater/FutureWater Dropbox/Team/Projects/Completed/2019/2019019_G4AW_MavoDiami_Angola/Data/2019019_MavoDiami_LV/2019019_MavoDiami"

# General Folder directories
DATA_DIR = os.path.join(angola_wd, "01_Data")
GIS_DIR = os.path.join(angola_wd, "02_GIS")
SCRIPTS_DIR = os.path.join(parent_wd, "03_Python_Scripts")
RESULTS_DIR = os.path.join(parent_wd, "04_Results")
TEMP_DIR = os.path.join(parent_wd, "05_Temp")

# Input data directories
PROVINCES = os.path.join(GIS_DIR, "Shapefiles", "AGO_adm1.shp")
COMMUNES = os.path.join(GIS_DIR, "Shapefiles", "AGO_adm3.shp")
DEM = os.path.join(DATA_DIR, "Elevation", "SRTM_30m_Angola_mask.tif")
HHS_DATA_DIR = os.path.join(DATA_DIR, "__TO_DROPBOX__", "Top_Subsoil")
CROPPING_CALENDER = os.path.join(SCRIPTS_DIR, "Cropping_calendar.csv")
CC_A = os.path.join(SCRIPTS_DIR, "Cropping_calendar_A.csv")
CC_B = os.path.join(SCRIPTS_DIR, "Cropping_calendar_B.csv")
PARAMETERS = os.path.join(SCRIPTS_DIR, "Parameters.csv")
SEASONAL_FORECAST_DIR = os.path.join(parent_wd, "Seasonalforecast")

# Output directories
NEW_DIR = os.path.join(parent_wd, "04_Results", "_LS_Results", "_Communes")
NEW_DROP = os.path.join(parent_wd, "Dropbox (FutureWater)",
                        "FW_VH_RK", "04_Results", "_LS_Results", "_Communes")

##############################################################################################
################################### Parameter setup ##########################################
##############################################################################################
# # If new crops are added, set this switch to 1 (all base data is already present like T, P, etc.)
# Make sure to choose the right cropping_calender file or add lines to the existing one.
SWITCH = 2  # 0: all, 1: new crops only, 2: seasonal forecast only

# Resolution and projection of final rasters
RES = 250  # meter
LOCAL_PRO = "EPSG:32733"

# Temperature input
T_LAPSE_RATE = -0.0065
T_PERC = [0.75, 0.95]
T_vars = ["tavg", "tmin", "tmax"]
T_perc_names = ["Warmer", "Much_Warmer"]

# Rainfall input
P_perc = [0.25, 0.05]
P_perc_names = ["Drier", "Much_Drier"]

##############################################################################################
#################### CREATE SCENARIOS AND CROPPING CALENDARS FOR EACH SCENARIO ###############
##############################################################################################
# Load DEM raster
# with rasterio.open(DEM) as src:
#     DEM_r = src.read(1)
#     dem_profile = src.profile

# Read provinces shapefile
provinces_shp = gpd.read_file(PROVINCES)
provinces_names = provinces_shp['NAME_1'].tolist()

# Read cropping calendars
cropping_cal = pd.read_csv(CROPPING_CALENDER, sep=",")
cropping_cal_a = pd.read_csv(CC_A, sep=",")
cropping_cal_b = pd.read_csv(CC_B, sep=",")

# Read parameters
params = pd.read_csv(PARAMETERS)

# Generate unique combinations of climate scenarios
unique_combis = []

# Add combinations of Average and Temperature scenarios
for t_name in T_perc_names:
    unique_combis.append(f"Average_{t_name}")

# Add combinations of Precipitation scenarios and Average
for p_name in P_perc_names:
    unique_combis.append(f"{p_name}_Average")

# Add combinations of Precipitation and Temperature scenarios
for p_name in P_perc_names:
    for t_name in T_perc_names:
        unique_combis.append(f"{p_name}_{t_name}")

# Add average scenario
unique_combis.append("Average_Average")

# Add prefix to all combinations
unique_combis = [f"Land_Suitability_{combo}" for combo in unique_combis]

##############################################################################################
##### CREATE DATAFRAME FOR FINAL FORECAST WITH CROP SUITABILITY COLUMNS FOR EACH COMMUNE #####
##############################################################################################
communes_shp = gpd.read_file(COMMUNES)
communes_shp_proj = communes_shp.to_crs(LOCAL_PRO)
communes_df = pd.DataFrame({
    "Province": communes_shp_proj['NAME_1'].str.replace(" ", "_"),
    "Municipality": communes_shp_proj['NAME_2'],
    "Commune": communes_shp_proj['NAME_3']
    # "Zone": communes_shp['Zone']  # Uncomment if Zone is in the shapefile
})

# Create month columns and prepare crop names
month_abbr = [calendar.month_abbr[i] for i in range(1, 13)]

# Create final forecast dataframe
final_forecast_df = pd.DataFrame({
    "Province": [None] * (len(communes_df) * len(cropping_cal)),
    "Commune": np.repeat(communes_df['Commune'].values, len(cropping_cal)),
    "Crop": [f"{crop}_{month_abbr[start-1]}_{month_abbr[end-1]}"
             for crop, start, end in zip(
                 cropping_cal['Crop'],
                 cropping_cal['Start_growing_season'],
                 cropping_cal['End_growing_season'])
             for _ in range(len(communes_df))]
})

# Add month columns
for month in month_abbr:
    final_forecast_df[month] = np.nan

# Add additional columns
additional_columns = [
    "Planting_Suitability_M1", "Planting_Suitability_M2", "Planting_Suitability_M3",
    "Crop_Start_Month", "Average_Average", "Threshold"
]
for col in additional_columns:
    final_forecast_df[col] = np.nan

# Fill Province column based on Commune
for i in range(len(final_forecast_df)):
    commune = final_forecast_df.loc[i, 'Commune']
    matching_province = communes_df.loc[communes_df['Commune']
                                        == commune, 'Province'].values
    if len(matching_province) > 0:
        final_forecast_df.loc[i, 'Province'] = matching_province[0]


# Add crop specific columns to communes_df
for _, row in cropping_cal.iterrows():
    crop = row['Crop']
    start_month = row['Start_growing_season']
    end_month = row['End_growing_season']
    name_col = f"{crop}_{month_abbr[start_month-1]}_{month_abbr[end_month-1]}"
    communes_df[name_col] = np.nan

# Create CSV files for each climate scenario
for comb in unique_combis:
    communes_df.to_csv(f"{comb}_Communes.csv")


##############################################################################################
###################### SETUP RASTER OPTIONS AND SWITCH VALUES ################################
##############################################################################################
# Set up raster options to prevent excessive memory or hard disk use.
with Env(
    GDAL_CONFIG_OPTIONS={
        'GDAL_MAX_MEM_ALLOC': '1e+10',  # Set maximum memory allocation
        'GDAL_CACHEMAX': '0.1e+10'      # Set chunk size
    },
    TEMP_DIR=TEMP_DIR
):
    start_time = datetime.datetime.now()

    # Hardcoded switch value (as in the original script)
    switch = 0
    n = 15 if switch == 2 else (10 if switch == 1 else 1)

    # Example: Running seasonal forecast (Script 8) if switch is 2
    # Import helper function for switch == 15
    if n == 15:
        # Implement your Seasonal Forecast function here
        print(
            "Need to convert and import: 01_Download_Seasonal_forecast_WI_API.R to Python")
        # In Python this would be something like:
        # from helper_functions.download_seasonal_forecast import download_seasonal_forecast
        pass

    ##############################################################################################
    ###################### RUN THROUGH ALL SCRIPTS FOR ALL PROVINCES #############################
    ##############################################################################################
    # Create list of all scripts with relative file paths
    list_scripts = [script for script in glob.glob(os.path.join(
        SCRIPTS_DIR, "*.py")) if "000_Run_All.py" not in script]
    scripts_dir = Path(SCRIPTS_DIR)
    list_scripts = [str(script.relative_to(scripts_dir)) for script in scripts_dir.glob(
        "*.py") if script.name != "000_Run_All.py"]

    # Process all scripts
    list_scripts = list_scripts[:1]
    for script in list_scripts:
        print("Script:", script)
        # Run the script for each province.
        province_names = provinces_names[:1]  # When running only one province
        for name in provinces_names:
            # Clean temporary files to have a clean folder to work with.
            temp_files = glob.glob(os.path.join(TEMP_DIR, "*"))
            for f in temp_files:
                try:
                    os.remove(f)
                except:
                    pass

            # Get province name as string
            province_shp = provinces_shp[provinces_shp['NAME_1'] == name]
            province_name = name.replace(" ", "_")

            # Run the script with subprocess
            print(
                f"Run Python : {script} for province:{province_name}.")
            myenv = os.environ.copy()
            myenv["PROVINCE"] = province_name
            # Try to run code with subprocess. When script fails, print error message.
            try:
                result = subprocess.run(["python", script],
                                        env=myenv,
                                        check=True,
                                        capture_output=True,
                                        text=True)
            except subprocess.CalledProcessError as e:
                print(
                    f"Error running script {script} for province {province_name}")
                print(e.stdout)
                print(e.stderr)

            # Clean temporary files
            temp_files = glob.glob(os.path.join(TEMP_DIR, "*.tif"))
            for f in temp_files:
                try:
                    os.remove(f)
                except:
                    pass

# End of for loop. All scripts for all provinces have been run.
end_time = datetime.datetime.now()
print(f"Time elapsed: {end_time - start_time}")

##############################################################################################
#################################### END OF SCRIPT ###########################################
##############################################################################################

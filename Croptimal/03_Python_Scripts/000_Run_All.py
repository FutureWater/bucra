# Import necessary libraries
import os
import glob
import datetime
import sys
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
# Define current and parent working directories
CURRENT_WD = os.getcwd()
PARENT_WD = os.path.dirname(CURRENT_WD)
BASE_WD = os.path.dirname(os.path.dirname(PARENT_WD))
angola_wd = "/Users/thomasfuturewater/FutureWater Dropbox/Team/Projects/Completed/2019/2019019_G4AW_MavoDiami_Angola/Data/2019019_MavoDiami_LV/2019019_MavoDiami"

# Define main folders
DATA_DIR = os.path.join(PARENT_WD, "01_Data")
GIS_DIR = os.path.join(BASE_WD, "GIS")     # Directory with shapefiles
RESULTS_DIR = os.path.join(PARENT_WD, "04_Results")
os.makedirs(RESULTS_DIR, exist_ok=True)
TEMP_DIR = os.path.join(PARENT_WD, "05_Temp")
SCRIPTS_DIR = os.path.join(PARENT_WD, "03_Python_Scripts")


# Input data paths
PROVINCES = os.path.join(GIS_DIR, "Nile_delta_bnd_adm1.shp")
# COMMUNES = os.path.join(GIS_DIR, "Shapefiles", "AGO_adm3.shp")
DEM = os.path.join(DATA_DIR, "DEM\DEM_NileDelta_250m.tif")
HHS_DATA_DIR = os.path.join(DATA_DIR, "__TO_DROPBOX__", "Top_Subsoil")
CROPPING_CALENDER = os.path.join(SCRIPTS_DIR, "Cropping_calendar.csv")
CC_A = os.path.join(SCRIPTS_DIR, "Cropping_calendar_A.csv")
CC_B = os.path.join(SCRIPTS_DIR, "Cropping_calendar_B.csv")
PARAMETERS = os.path.join(SCRIPTS_DIR, "Parameters.csv")
SEASONAL_FORECAST_DIR = os.path.join(PARENT_WD, "Seasonalforecast")

# Output directories
NEW_RESULTS_SUBDIR = os.path.join(RESULTS_DIR, "_LS_Results", "_Communes")
NEW_DROP = os.path.join(PARENT_WD, "Dropbox (FutureWater)",
                        "FW_VH_RK", "04_Results", "_LS_Results", "_Communes")

##############################################################################################
################################### Parameter setup ##########################################
##############################################################################################
# # If new crops are added, set this switch to 1 (all base data is already present like T, P, etc.)
# Make sure to choose the right cropping_calender file or add lines to the existing one.
SWITCH = 2  # 0: all, 1: new crops only, 2: seasonal forecast only

# Resolution and projection of final rasters
RES = 250  # meter
LOCAL_PROJ = "EPSG:32636"

# Temperature input
T_LAPSE_RATE = -0.0065
T_PERC = [0.75, 0.95]
T_vars = ["Tmean", "Tmin", "Tmax"]
T_perc_names = ["Warmer", "Much_Warmer"]

# Rainfall input
P_perc = [0.25, 0.05]
P_perc_names = ["Drier", "Much_Drier"]

# Read provinces shapefile
province_shp = gpd.read_file(PROVINCES)
province_names = province_shp["ADM1_EN"].tolist()



##############################################################################################
###################### SETUP RASTER OPTIONS AND SWITCH VALUES ################################
##############################################################################################
# Set up raster options to prevent excessive memory or hard disk use.
def main():
    with Env(
        GDAL_CONFIG_OPTIONS={
            'GDAL_MAX_MEM_ALLOC': '1e+10',  # Set maximum memory allocation
            'GDAL_CACHEMAX': '0.1e+10'      # Set chunk size
        },
        TEMP_DIR=TEMP_DIR
    ):
        start_time = datetime.datetime.now()

        # Hardcoded switch value
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
            SCRIPTS_DIR, "*.py")) if not ("000_Run_All.py" or '04_Seasonal')in script]
        # Sort the list of scripts from in correct order
        list_scripts.sort()
        # print(list_scripts)

        # Process all scripts
        for script in list_scripts:
            script_name = os.path.basename(script)[:-3]

            # Run the script for each province.
            # When running only one province
            for name in province_names[-1:]:
                # Clean temporary files to have a clean folder to work with.
                temp_files = glob.glob(os.path.join(TEMP_DIR, "*"))

                for f in temp_files:
                    try:
                        os.remove(f)
                    except:
                        pass

                # Get province name as string
                #name = name.replace(" ", "_")

                # Run the script with subprocess. Give all constants as environment variables.
                print(
                    f"Run {script_name} for {name}.")
                myenv = os.environ.copy()
                myenv["PROVINCE"] = name
                #myenv["RESOLUTION"] = RES
                myenv["CWD"] = CURRENT_WD
                myenv["LOCAL_PROJ"] = LOCAL_PROJ
                myenv["CROPPING_CALENDER"] = CROPPING_CALENDER

                # Try to run code with subprocess. When script fails, print error message.
                try:
                    # result = subprocess.run(["python", script],
                    #                         env=myenv,
                    #                         check=True,
                    #                         text=True)
                    
                    process = subprocess.Popen(["python", script],
                             env=myenv,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             universal_newlines=True)
    
                    # Print output in real-time
                    for stdout_line in iter(process.stdout.readline, ""):
                        print(stdout_line, end="")
                        sys.stdout.flush()
                    
                    process.stdout.close()
                    return_code = process.wait()
                    if return_code:
                        raise subprocess.CalledProcessError(return_code, process.args)
                    
                except subprocess.CalledProcessError as e:
                    print(
                        f"Error running script {script} for province {name}")
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

    print("\n\nAll scripts have been run and all temporary files have been deleted.")


if __name__ == "__main__":
    main()
##############################################################################################
#################################### END OF SCRIPT ###########################################
##############################################################################################

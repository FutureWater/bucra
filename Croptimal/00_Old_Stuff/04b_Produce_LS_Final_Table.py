import os
import glob
import calendar
import pandas as pd
import geopandas as gpd
import rasterio
import numpy as np
import rioxarray

"""
This script produces land suitability tables for communes by:
1) Extracting mean suitability values from each Land Suitability raster per commune
2) Generating human-readable recommendations based on suitability thresholds
3) Creating both descriptive recommendation tables and numerical suitability value tables
4) Saving these tables as CSV files for each climate scenario combination
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
LS_RESULTS = os.path.join(RESULTS_DIR, '_LS_Results')
CROPPING_CAL = pd.read_csv(os.path.join(current_wd, "Cropping_calendar_A.csv"))

# Define constants
RES = 250  # Resolution in meters
NO_DATA_VALUE = 0  # No data value
LOCAL_PROJ = "EPSG:32733"

# Define climate scenario combinations
UNIQUE_COMBIS = ["Dry_Warm", "Dry_Much_Warmer", "Dry_Average",
                 "Normal_Warm", "Normal_Much_Warmer", "Normal_Average",
                 "Much_Drier_Warmer", "Much_Drier_Much_Warmer", "Much_Drier_Average",]

# Create output directory
new_results_subdir = os.path.join(RESULTS_DIR, "_LS_Results", "_Communes")
os.makedirs(new_results_subdir, exist_ok=True)

# Load communes shapefile
COMMUNES_FILE = os.path.join(GIS_DIR, "Shapefiles", "AGO_adm3.shp")
commune_shp = gpd.read_file(COMMUNES_FILE)
commune_shp_sel = commune_shp[commune_shp["NAME_1"] == PROVINCE_NAME]
commune_shp_proj = commune_shp_sel.to_crs(LOCAL_PROJ)

# Find all land suitability rasters
all_results = glob.glob(os.path.join(
    LS_RESULTS, "Weighted", f"Land_Suitability*.tif"), recursive=True)

# Initialize dictionaries to store dataframes for each combination
communes_csv_dict = {}
communes_ls_values_dict = {}

# Calculate mean suitability values for each commune across all rasters
print(f"Calculating zonal statistics for {len(all_results)} raster files...")

# Create a combined dataframe to hold all results
mean_all_results_commune = pd.DataFrame(commune_shp_proj["NAME_3"],)

###############################################################################################
##### Process each raster file to get commune-level mean suitability result for each crop #####
###############################################################################################
for raster_file in all_results:
    # Extract the basename without extension for the column name
    file_basename = os.path.basename(raster_file).replace('.tif', '')[17:]

    # Calculate zonal stats
    with rioxarray.open_rasterio(raster_file) as raster:
        column_name = file_basename
        for idx, row in commune_shp_proj.iterrows():
            commune_name = row['NAME_3']

            # Clip raster to commune geometry and calculate mean
            clipped = raster.rio.clip([row.geometry], drop=False)
            mean_value = float(clipped.mean())
            results[commune_name][column_name] = mean_value

    # Extract statistics into a dictionary
    commune_means = {}
    for idx, stat in enumerate(stats):
        commune_id = province_communes_shp.iloc[idx].name
        mean_value = stat['properties']['mean']
        commune_means[commune_id] = mean_value

    # Add to the combined dataframe
    mean_all_results_commune[file_basename] = pd.Series(commune_means)

# Process each climate combination
for comb in UNIQUE_COMBIS:
    print(f"Processing {comb} combination...")

    # Initialize dataframes if first province
    if PROVINCE_NAME == PROVINCES_NAMES[0]:
        # Initialize with commune metadata columns
        communes_csv_dict[comb] = province_communes.copy()
        communes_ls_values_dict[comb] = province_communes.copy()
    else:
        # Load existing dataframes if not first province
        try:
            communes_csv_path = os.path.join(
                new_results_subdir, f"{comb}_Communes.csv")
            communes_ls_values_path = os.path.join(
                new_results_subdir, f"{comb}_Communes_LS_values.csv")

            if os.path.exists(communes_csv_path) and os.path.exists(communes_ls_values_path):
                communes_csv_dict[comb] = pd.read_csv(communes_csv_path)
                communes_ls_values_dict[comb] = pd.read_csv(
                    communes_ls_values_path)
            else:
                # Initialize with commune metadata columns if files don't exist
                communes_csv_dict[comb] = province_communes.copy()
                communes_ls_values_dict[comb] = province_communes.copy()
        except Exception as e:
            print(f"Error loading existing dataframes for {comb}: {e}")
            # Initialize with commune metadata columns if error occurs
            communes_csv_dict[comb] = province_communes.copy()
            communes_ls_values_dict[comb] = province_communes.copy()

    # Filter results for this combination
    comb_results = mean_all_results_commune.filter(regex=comb)

    # Extract crop, start_month and end_month from column names
    # Assuming format is: Land_Suitability_Dry_Cold_Wheat_Nov-Apr.tif
    crop_info = {}
    for col in comb_results.columns:
        parts = col.split('_')
        if len(parts) >= 5:  # Ensure we have enough parts
            crop_name = parts[3]  # Crop name is in position 3 (0-indexed)
            season = parts[4].split('.')[0]  # Remove .tif extension
            start_month, end_month = season.split('-')

            crop_info[col] = {
                'crop': crop_name,
                'start_month': start_month,
                'end_month': end_month
            }

    # Generate recommendations and update dataframes
    for col in comb_results.columns:
        if col in crop_info:
            crop = crop_info[col]['crop']
            start_month = crop_info[col]['start_month']
            end_month = crop_info[col]['end_month']

            # Get crop display name from cropping calendar
            try:
                crop_message = cropping_cal.loc[cropping_cal['Crop']
                                                == crop, 'Name_message'].values[0]
            except:
                # Use the crop name if message is not available
                crop_message = crop

            # Convert month abbreviations to full names
            month_abbrs = list(calendar.month_abbr)[
                1:]  # Skip empty first element
            start_month_num = month_abbrs.index(start_month) + 1
            end_month_num = month_abbrs.index(end_month) + 1
            start_month_full = calendar.month_name[start_month_num]
            end_month_full = calendar.month_name[end_month_num]

            # Generate messages based on suitability values
            for idx, value in comb_results[col].items():
                # Determine suitability level
                if value >= 0.8:
                    suitability = "very suitable"
                elif value >= 0.7:
                    suitability = "suitable"
                else:
                    suitability = ""

                # Generate the recommendation message
                if value < 0.7:
                    message = (f"Growing {crop_message} from {start_month_full} to {end_month_full} "
                               f"is NOT ideal in your location. If you still want to plant this crop: "
                               f"You may have to plant a variety resistant to extreme temperatures. "
                               f"You may have to irrigate. You may have to fertilize a lot. "
                               f"You may have to dig drains.")
                else:
                    message = (f"Growing {crop_message} from {start_month_full} to {end_month_full} "
                               f"is {suitability} in your location.")

                # Update recommendation dataframe
                commune_row = communes_csv_dict[comb].index[communes_csv_dict[comb]
                                                            ['Province'] == PROVINCE_NAME]
                column_name = f"{crop}_{start_month}_{end_month}"

                # Create column if it doesn't exist
                if column_name not in communes_csv_dict[comb].columns:
                    communes_csv_dict[comb][column_name] = ""
                    communes_ls_values_dict[comb][column_name] = np.nan

                # Update values
                communes_csv_dict[comb].loc[commune_row, column_name] = message
                communes_ls_values_dict[comb].loc[commune_row, column_name] = round(
                    value, 3)

# Save final results if this is the last province
if PROVINCE_NAME == PROVINCES_NAMES[-1]:
    print("Saving final results...")

    for comb in UNIQUE_COMBIS:
        # Save recommendation tables
        communes_csv_dict[comb].to_csv(
            os.path.join(new_results_subdir, f"{comb}_Communes.csv"),
            index=False, quoting=3  # QUOTE_NONE in pandas
        )

        # Save numerical value tables
        communes_ls_values_dict[comb].to_csv(
            os.path.join(new_results_subdir, f"{comb}_Communes_LS_values.csv"),
            index=False, quoting=3  # QUOTE_NONE in pandas
        )

    print("All results saved successfully")

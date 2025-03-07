import os
import glob
import calendar
import pandas as pd
import geopandas as gpd
import rasterio
import numpy as np
from rasterstats import zonal_stats

"""
This script produces land suitability tables for communes by:
1) Extracting mean suitability values from each Land Suitability raster per commune
2) Generating human-readable recommendations based on suitability thresholds
3) Creating both descriptive recommendation tables and numerical suitability value tables
4) Saving these tables as CSV files for each climate scenario combination
"""

# Define constants
RESULTS_DIR = "your_results_directory"  # Replace with actual path
PROVINCE_NAME = "your_province_name"  # Replace with actual province name
# Replace with actual province list
PROVINCES_NAMES = ["province1", "province2", "province3"]

# Define climate scenario combinations
UNIQUE_COMBIS = ["Dry_Cold", "Dry_Cool", "Dry_Average",
                 "Normal_Cold", "Normal_Cool", "Normal_Average",
                 "Average_Cold", "Average_Cool", "Average_Average"]

# Create output directory
new_dir = os.path.join(RESULTS_DIR, "_LS_Results", "_Communes")
os.makedirs(new_dir, exist_ok=True)

# Load commune data for the current province
communes_df = pd.read_csv(os.path.join(
    RESULTS_DIR, "communes.csv"))  # Adjust path as needed
province_communes = communes_df[communes_df['Province'] == PROVINCE_NAME]

# Load commune shapefile
communes_shp = gpd.read_file(os.path.join(
    RESULTS_DIR, "communes.shp"))  # Adjust path as needed
province_communes_shp = communes_shp[communes_shp['NAME_1'].str.replace(
    ' ', '_') == PROVINCE_NAME]

# Load crop calendar
cropping_cal = pd.read_csv(os.path.join(
    RESULTS_DIR, "cropping_calendar.csv"))  # Adjust path as needed

# Find all land suitability rasters
all_results = glob.glob(os.path.join(
    RESULTS_DIR, "**", f"Land_Suitability*{PROVINCE_NAME}*.tif"), recursive=True)

# Initialize dictionaries to store dataframes for each combination
communes_csv_dict = {}
communes_ls_values_dict = {}

# Calculate mean suitability values for each commune across all rasters
print(f"Calculating zonal statistics for {len(all_results)} raster files...")

# Create a combined dataframe to hold all results
mean_all_results_commune = pd.DataFrame()

# Process each raster file to get commune-level mean suitability result for each crop
for raster_file in all_results:
    # Extract the basename without extension for the column name
    file_basename = os.path.basename(raster_file).replace('.tif', '')

    # Calculate zonal statistics for this raster
    stats = zonal_stats(
        province_communes_shp,
        raster_file,
        stats="mean",
        geojson_out=True
    )

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
            communes_csv_path = os.path.join(new_dir, f"{comb}_Communes.csv")
            communes_ls_values_path = os.path.join(
                new_dir, f"{comb}_Communes_LS_values.csv")

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
            os.path.join(new_dir, f"{comb}_Communes.csv"),
            index=False, quoting=3  # QUOTE_NONE in pandas
        )

        # Save numerical value tables
        communes_ls_values_dict[comb].to_csv(
            os.path.join(new_dir, f"{comb}_Communes_LS_values.csv"),
            index=False, quoting=3  # QUOTE_NONE in pandas
        )

    print("All results saved successfully")

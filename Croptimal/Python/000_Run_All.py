import os
import pandas as pd
import geopandas as gpd
import rasterio
import numpy as np
from rasterio.enums import Resampling
from datetime import datetime
import glob

# Clean working directory (remove variables if needed)
# Not required in Python as we don't have a global environment like R's workspace

# Define directories
data_drop = "C:/Users/Lisa/Dropbox (FutureWater)/FW_VH_RK/01_Data/"
results_drop = "C:/Users/Lisa/Dropbox (FutureWater)/FW_VH_RK/04_Results/"
provinces = "C:/Users/Lisa/Dropbox (FutureWater)/FW_VH_RK/02_GIS/Shapefiles/AGO_adm1.shp"
DEM = "C:/Users/Lisa/Dropbox (FutureWater)/FW_VH_RK/01_Data/Elevation/SRTM_30m_Angola_mask.tif"
HHS_data_dir = "H:/01_Data/__TO_DROPBOX__/Top_Subsoil/"
cropping_calender = "C:/Users/Lisa/Documents/Projects/2019019_MavoDiami/03_R_Scripts/Cropping_calendar.csv"
cc_a = "C:/Users/Lisa/Documents/Projects/2019019_MavoDiami/03_R_Scripts/Cropping_calendar_A.csv"
cc_b = "C:/Users/Lisa/Documents/Projects/2019019_MavoDiami/03_R_Scripts/Cropping_calendar_B.csv"
parameters = "C:/Users/Lisa/Documents/Projects/2019019_MavoDiami/03_R_Scripts/Parameters.csv"
communes = "C:/Users/Lisa/Documents/Projects/2019019_MavoDiami/Shapefiles/AGO_adm3.shp"
Scripts_dir = "C:/Users/Lisa/Documents/Projects/2019019_MavoDiami/03_R_Scripts"

# Set resolution and projection of final rasters
res = 250  # meter
Local_proj = "EPSG:32733"

# Temperature and Rainfall input parameters
T_vars = ["tavg", "tmin", "tmax"]
T_lapse_rate = -0.0065
T_perc = [0.75, 0.95]
T_perc_names = ["Warmer", "Much_Warmer"]
P_perc = [0.25, 0.05]
P_perc_names = ["Drier", "Much_Drier"]

# Read input files
DEM_r = rasterio.open(DEM)
provinces_shp = gpd.read_file(provinces)
provinces_names = provinces_shp['NAME_1'].unique()
cropping_cal = pd.read_csv(cropping_calender)
cropping_cal_a = pd.read_csv(cc_a)
cropping_cal_b = pd.read_csv(cc_b)
params = pd.read_csv(parameters)

# Create unique combinations for land suitability
unique_combis = [f"Land_Suitability_{x}_{y}" for x in ["Average"] + T_perc_names for y in ["Average"] + P_perc_names]
unique_combis.append("Land_Suitability_Average_Average")

# Read communes shapefile
communes_shp = gpd.read_file(communes)
communes_shp_proj = communes_shp.to_crs(Local_proj)

communes_df = communes_shp[['NAME_1', 'NAME_2', 'NAME_3']].rename(columns={'NAME_1': 'Province', 'NAME_2': 'Municipality', 'NAME_3': 'Commune'})
communes_df['Province'] = communes_df['Province'].str.replace(" ", "_")

# Prepare forecast dataframe
final_forecast_df = pd.DataFrame({
    "Province": [None] * (len(communes_df) * len(cropping_cal)),
    "Commune": np.repeat(communes_df['Commune'], len(cropping_cal)),
    "Crop": [f"{row['Crop']}_{month[:3]}_{month[4:]}" for _, row in cropping_cal.iterrows() for month in zip(cropping_cal['Start_growing_season'], cropping_cal['End_growing_season'])],
    "Planting_Suitability_M1": [None] * (len(communes_df) * len(cropping_cal)),
    "Planting_Suitability_M2": [None] * (len(communes_df) * len(cropping_cal)),
    "Planting_Suitability_M3": [None] * (len(communes_df) * len(cropping_cal)),
    "Crop_Start_Month": [None] * (len(communes_df) * len(cropping_cal)),
    "Average_Average": [None] * (len(communes_df) * len(cropping_cal)),
    "Threshold": [None] * (len(communes_df) * len(cropping_cal)),
})

for i, row in final_forecast_df.iterrows():
    final_forecast_df.loc[i, 'Province'] = communes_df['Province'][communes_df['Commune'] == row['Commune']].values[0]

# Add columns to communes_df
for _, row in cropping_cal.iterrows():
    crop = row['Crop']
    start_month = row['Start_growing_season']
    end_month = row['End_growing_season']
    name_col = f"{crop}_{month.abb[start_month-1]}_{month.abb[end_month-1]}"
    communes_df[name_col] = np.nan

# Loop through and save files for each combination
for comb in unique_combis:
    communes_df.to_csv(f"{comb}_Communes.csv")

# Set options for raster processing
rasterio.env.default_options.update({
    'overwrite': True,
    'max_memory': 1e+10,
    'chunksize': 0.1e+10,
    'tmpdir': temp_dir
})

start_time = datetime.now()

switch = 2  # HARDCODED
n = 15 if switch == 2 else 10

# Example: Running seasonal forecast (Script 8) if switch is 2
if n == 15:
    # Implement your Seasonal Forecast function here
    pass

# Run through scripts for each province
list_scripts = [f for f in glob.glob(os.path.join(Scripts_dir, "*.R")) if "000_Run_All.R" not in f]

for script in list_scripts[n:]:
    for name in provinces_names:
        # Processing for each province
        province_shp = provinces_shp[provinces_shp['NAME_1'] == name]
        province_name = name.replace(" ", "_")
        indir = os.path.join(results_dir, province_name)
        
        print(f"Province: {province_name} & script: {os.path.basename(script)}")
        # Execute the script (you'll need to implement a way to translate R script into Python functions)
        
        # Clean up temporary files after processing
        temp_files = glob.glob(os.path.join(temp_dir, "*.grd")) + glob.glob(os.path.join(temp_dir, "*.gri"))
        for file in temp_files:
            os.remove(file)

end_time = datetime.now()
print(f"Total time taken: {end_time - start_time}")

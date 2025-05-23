import os
import glob
import json
import calendar
import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from rasterio.features import rasterize
from rasterio.warp import reproject, Resampling, calculate_default_transform
from rasterstats import zonal_stats
import ntpath

"""
This script processes seasonal forecast data and creates planting suitability recommendations:
1) Processes temperature and precipitation forecast data from JSON files
2) Converts forecast data to raster grids
3) Calculates average values for each commune
4) Compares forecast values with historical data
5) Creates warnings when conditions are significantly different from normal
"""

# Define constants
RESULTS_DIR = "your_results_directory"  # Replace with actual path
PROVINCE_NAME = "your_province_name"  # Replace with actual province name
SEASONAL_FORECAST_DIR = "your_forecast_directory"  # Replace with actual path
RES = 250  # Resolution in meters
T_LAPSE_RATE = -0.0065  # Temperature lapse rate (°C/m)
# Replace with actual province list
PROVINCES_NAMES = ["province1", "province2"]

# Helper function to extract file parts


def split_path(path):
    """Split a file path into its components"""
    head, tail = ntpath.split(path)
    if not head:
        return [tail]
    return [tail] + split_path(head)


# Find seasonal forecast files
seasonal_temp_files = glob.glob(os.path.join(
    SEASONAL_FORECAST_DIR, "**", "*temperature*.json"), recursive=True)
seasonal_prec_files = glob.glob(os.path.join(
    SEASONAL_FORECAST_DIR, "**", "*rainfall*.json"), recursive=True)

# Load communes data
communes_df = pd.read_csv(os.path.join(
    RESULTS_DIR, "communes.csv"))  # Adjust path as needed
communes_shp = gpd.read_file(os.path.join(
    RESULTS_DIR, "communes.shp"))  # Adjust path as needed
communes_selected = communes_shp[communes_df['Province'] == PROVINCE_NAME]

# Create result dataframe with commune names
result_df = pd.DataFrame(communes_selected['NAME_3'])

# Load crop calendar
cropping_cal = pd.read_csv(os.path.join(
    RESULTS_DIR, "cropping_calendar.csv"))  # Adjust path as needed

# Initialize final forecast dataframe if first province
if PROVINCE_NAME == PROVINCES_NAMES[0]:
    # Create columns for all communes and crops
    all_communes = communes_df['Commune'].tolist()
    all_crops = []
    for _, crop_row in cropping_cal.iterrows():
        crop = crop_row['Crop']
        start_month = crop_row['Start_growing_season']
        end_month = crop_row['End_growing_season']
        crop_season = f"{crop}_{calendar.month_abbr[start_month]}_{calendar.month_abbr[end_month]}"
        all_crops.append(crop_season)

    # Create empty dataframe with structure matching R script's final_forecast_df
    final_forecast_df = pd.DataFrame(columns=[
        'Province', 'Commune', 'Crop',
        'Forecast_M1', 'Planting_Suitability_M1',
        'Forecast_M2', 'Planting_Suitability_M2',
        'Forecast_M3', 'Planting_Suitability_M3',
        'Average_Average', 'Threshold', 'Crop_Start_Month',
        'M1', 'M2', 'M3'
    ])

    # Initialize rows
    for province in PROVINCES_NAMES:
        province_communes = communes_df[communes_df['Province']
                                        == province]['Commune'].tolist()
        for commune in province_communes:
            for crop in all_crops:
                final_forecast_df = final_forecast_df.append({
                    'Province': province,
                    'Commune': commune,
                    'Crop': crop,
                    'Average_Average': np.nan,
                    'Threshold': np.nan,
                    'Crop_Start_Month': np.nan,
                    'M1': 'No',
                    'M2': 'No',
                    'M3': 'No'
                }, ignore_index=True)

###################################
# Process temperature forecast data
###################################
for t_idx, temp_file in enumerate(seasonal_temp_files):
    print(f"Processing temperature file {t_idx+1}/{len(seasonal_temp_files)}")

    # Extract temperature statistic from file path
    t_stat = split_path(temp_file)[1]  # tmin, tmax, or tavg

    # Load JSON data
    with open(temp_file, 'r') as f:
        t_stat_data = json.load(f)

    # Create dataframe with coordinates and temperature data
    temp_df = pd.DataFrame({
        'X': t_stat_data['GridDefinition']['Longitude'],
        'Y': t_stat_data['GridDefinition']['Latitude'],
        'M0': [float(x) if x != '-9999' else np.nan for x in t_stat_data['Data']['0']['Data']],
        'M1': [float(x) if x != '-9999' else np.nan for x in t_stat_data['Data']['1']['Data']],
        'M2': [float(x) if x != '-9999' else np.nan for x in t_stat_data['Data']['2']['Data']],
        # Uncomment if more months are needed
        # 'M3': [float(x) if x != '-9999' else np.nan for x in t_stat_data['Data']['3']['Data']],
        # 'M4': [float(x) if x != '-9999' else np.nan for x in t_stat_data['Data']['4']['Data']],
    })

    # Get dates from data
    dates = [
        t_stat_data['Data']['0']['Date'],
        t_stat_data['Data']['1']['Date'],
        t_stat_data['Data']['2']['Date'],
        # Uncomment if more months are needed
        # t_stat_data['Data']['3']['Date'],
        # t_stat_data['Data']['4']['Date'],
    ]
    dates_t = dates.copy()

    # Rename columns to use dates
    column_mapping = {'X': 'X', 'Y': 'Y'}
    for i, date in enumerate(dates):
        column_mapping[f'M{i}'] = date
    temp_df.rename(columns=column_mapping, inplace=True)

    # Create rasters from points
    x_coords = temp_df['X'].values
    y_coords = temp_df['Y'].values

    # Create base raster
    t_ras = None
    temp_rasters = []

    for date in dates:
        # Get data for this date
        values = temp_df[date].values

        # Determine raster dimensions
        x_min, x_max = min(x_coords), max(x_coords)
        y_min, y_max = min(y_coords), max(y_coords)

        # Assuming regular grid, calculate cell size
        cell_size_x = (x_max - x_min) / (len(np.unique(x_coords)) - 1)
        cell_size_y = (y_max - y_min) / (len(np.unique(y_coords)) - 1)

        # Create transform
        transform = rasterio.transform.from_origin(x_min - cell_size_x/2,
                                                   y_max + cell_size_y/2,
                                                   cell_size_x, cell_size_y)

        # Create raster shape (adjust based on your data structure)
        width = len(np.unique(x_coords))
        height = len(np.unique(y_coords))

        # Create empty raster
        raster_data = np.full((height, width), np.nan)

        # Fill raster with values (basic approach - may need optimization)
        points = np.column_stack((x_coords, y_coords))
        for i, point in enumerate(points):
            # Convert coordinates to raster indices
            col = int((point[0] - (x_min - cell_size_x/2)) / cell_size_x)
            row = int(((y_max + cell_size_y/2) - point[1]) / cell_size_y)
            if 0 <= row < height and 0 <= col < width:
                raster_data[row, col] = values[i]

        # Save to temp file
        output_meta = {
            'driver': 'GTiff',
            'height': height,
            'width': width,
            'count': 1,
            'dtype': rasterio.float32,
            'crs': 'EPSG:4326',
            'transform': transform,
            'nodata': np.nan
        }

        temp_raster_path = os.path.join(
            RESULTS_DIR, "temp", f"temp_{t_stat}_{date}.tif")
        os.makedirs(os.path.dirname(temp_raster_path), exist_ok=True)

        with rasterio.open(temp_raster_path, 'w', **output_meta) as dst:
            dst.write(raster_data.astype(rasterio.float32), 1)

        temp_rasters.append(temp_raster_path)

        # Store first raster as t_ras
        if t_ras is None:
            t_ras = temp_raster_path

    # Load DEM and resample temperature rasters
    dem_path = os.path.join(RESULTS_DIR, PROVINCE_NAME,
                            "DEM", f"DEM_{PROVINCE_NAME}_{RES}m_diff.tif")
    with rasterio.open(dem_path) as dem_src:
        dem_data = dem_src.read(1)
        dem_meta = dem_src.meta
        dem_transform = dem_src.transform
        dem_crs = dem_src.crs
        dem_shape = dem_src.shape

    # Process each month
    months = [calendar.month_abbr[int(date[4:6])] for date in dates]
    months_num = [int(date[4:6]) for date in dates]

    # Get crops that start in relevant months
    crops = cropping_cal[cropping_cal['Start_growing_season'] >= months_num[0]]

    for month_idx, month in enumerate(months):
        # Resample temp raster to DEM resolution
        with rasterio.open(temp_rasters[month_idx]) as src:
            temp_data = src.read(1)

            # Reproject to DEM
            resampled_data = np.zeros(dem_shape, dtype=np.float32)
            reproject(
                source=rasterio.band(src, 1),
                destination=resampled_data,
                src_transform=src.transform,
                src_crs=src.crs,
                dst_transform=dem_transform,
                dst_crs=dem_crs,
                resampling=Resampling.bilinear
            )

            # Apply lapse rate correction
            final_data = resampled_data + dem_data * T_LAPSE_RATE

            # Calculate zonal stats for communes
            communes_mean = zonal_stats(
                communes_selected,
                final_data,
                affine=dem_transform,
                stats="mean",
                nodata=np.nan
            )

            # Extract means
            means = [stat['mean'] for stat in communes_mean]

            # Add to result dataframe
            column_name = f"mean_{t_stat}_{month}"
            result_df[column_name] = means

###################################
# Process precipitation forecast data
###################################
for p_idx, prec_file in enumerate(seasonal_prec_files):
    print(
        f"Processing precipitation file {p_idx+1}/{len(seasonal_prec_files)}")

    # Extract precipitation statistic from file path
    p_stat = split_path(prec_file)[1]  # pmin, pmax, or pavg

    # Load JSON data
    with open(prec_file, 'r') as f:
        p_stat_data = json.load(f)

    # Create dataframe with coordinates and precipitation data
    temp_df = pd.DataFrame({
        'X': p_stat_data['GridDefinition']['Longitude'],
        'Y': p_stat_data['GridDefinition']['Latitude'],
        'M0': [float(x) if x != '-9999' else np.nan for x in p_stat_data['Data']['0']['Data']],
        'M1': [float(x) if x != '-9999' else np.nan for x in p_stat_data['Data']['1']['Data']],
        'M2': [float(x) if x != '-9999' else np.nan for x in p_stat_data['Data']['2']['Data']],
        # Uncomment if more months are needed
        # 'M3': [float(x) if x != '-9999' else np.nan for x in p_stat_data['Data']['3']['Data']],
        # 'M4': [float(x) if x != '-9999' else np.nan for x in p_stat_data['Data']['4']['Data']],
    })

    # Use the same dates as temperature for consistency
    temp_df.rename(columns={
                   'X': 'X', 'Y': 'Y', 'M0': dates_t[0], 'M1': dates_t[1], 'M2': dates_t[2]}, inplace=True)

    # Same process as temperature, but using the temperature raster as reference
    with rasterio.open(t_ras) as src:
        t_meta = src.meta.copy()
        t_transform = src.transform
        t_crs = src.crs

    # Create rasters for each month
    for month_idx, date in enumerate(dates_t):
        # Get data for this month
        temp_matrix = temp_df.loc[:, ['X', 'Y', date]]
        temp_matrix = temp_matrix.values

        # Rasterize points onto temperature raster grid
        with rasterio.open(t_ras) as src:
            # Create points to raster
            temp_path = os.path.join(
                RESULTS_DIR, "temp", f"temp_{p_stat}_{date}.tif")

            # Calculate zonal statistics on precipation raster
            with rasterio.open(dem_path) as dem_src:
                # Resample to match DEM
                resampled_data = np.zeros(dem_shape, dtype=np.float32)

                # Assuming we've already created temp rasters similar to temperature above
                with rasterio.open(temp_path) as p_src:
                    reproject(
                        source=rasterio.band(p_src, 1),
                        destination=resampled_data,
                        src_transform=p_src.transform,
                        src_crs=p_src.crs,
                        dst_transform=dem_transform,
                        dst_crs=dem_crs,
                        resampling=Resampling.bilinear
                    )

                # Calculate zonal stats for communes
                communes_mean = zonal_stats(
                    communes_selected,
                    resampled_data,
                    affine=dem_transform,
                    stats="mean",
                    nodata=np.nan
                )

                # Extract means
                means = [stat['mean'] for stat in communes_mean]

                # Add to result dataframe
                month = months[month_idx]
                column_name = f"mean_{p_stat}_{month}"
                result_df[column_name] = means

# Load historical statistics
p_stats = pd.read_csv(os.path.join(RESULTS_DIR, "_LS_Results",
                      "Rainfall", f"Average_Rainfall_per_Commune_{PROVINCE_NAME}.csv"))
t_stats = pd.read_csv(os.path.join(RESULTS_DIR, "_LS_Results", "Temperature",
                      f"Average_Temperature_per_Commune_{PROVINCE_NAME}.csv"))

# Create results dataframe
df_results_t_p = result_df.copy()

# Compare forecast to historical data
for comm_idx in range(len(result_df)):
    for month_idx, month in enumerate(months):
        # Get columns for this month in stats dataframes
        cols_stats = [col for col in p_stats.columns if month in col]

        # Get values for this commune
        p_values_stats = p_stats.iloc[comm_idx][cols_stats].values
        t_values_stats = t_stats.iloc[comm_idx][cols_stats].values

        # Get forecast values for this commune
        cols_sf = [col for col in df_results_t_p.columns if month in col]
        values_sf = df_results_t_p.iloc[comm_idx][cols_sf].values

        # Compare precipitation (values_sf[4] is pavg)
        if values_sf[3] <= p_values_stats[0]:
            p_res = "Much_Drier"
        elif values_sf[3] <= p_values_stats[1]:
            p_res = "Drier"
        else:
            p_res = "Average"

        # Compare temperature (values_sf[1] is tavg)
        if values_sf[0] >= t_values_stats[1]:
            t_res = "Much_Warmer"
        elif values_sf[0] >= t_values_stats[0]:
            t_res = "Warmer"
        else:
            t_res = "Average"

        # Add forecast result
        forecast_col = f"Forecast_{month}"
        df_results_t_p.loc[comm_idx, forecast_col] = f"{p_res}_{t_res}"

# Process each month
new_dir = os.path.join(RESULTS_DIR, "_LS_Results", "_Communes")
os.makedirs(new_dir, exist_ok=True)

for month_idx, month in enumerate(months):
    print(f"Processing crops for {month}")

    # Create forecast dataframe for this month
    month_cols = [df_results_t_p.columns[0]] + \
        [col for col in df_results_t_p.columns if month in col]
    forecast_df = df_results_t_p[month_cols].copy()

    # Get month number
    month_num = months_num[month_idx]

    # Get crops that are relevant for this month
    all_crops = cropping_cal[(cropping_cal['Start_growing_season'] == month_num) |
                             (cropping_cal['End_growing_season'] >= month_num) |
                             (cropping_cal['End_growing_season'] <= 12)]

    # Create crop column names
    crops = []
    for _, crop_row in all_crops.iterrows():
        crop_name = crop_row['Crop']
        start_month = crop_row['Start_growing_season']
        end_month = crop_row['End_growing_season']
        crop_season = f"{crop_name}_{calendar.month_abbr[start_month]}_{calendar.month_abbr[end_month]}"
        crops.append(crop_season)

    # Add empty columns for each crop
    for crop in crops:
        forecast_df[crop] = np.nan

    # Process each crop
    for crop_idx, crop in enumerate(crops):
        for comm_idx in range(len(forecast_df)):
            # Skip if forecast is NA
            if pd.isna(forecast_df.iloc[comm_idx][7]):
                continue

            # Get commune name
            commune = forecast_df.iloc[comm_idx, 0]

            # Read appropriate Land Suitability file
            forecast_condition = forecast_df.iloc[comm_idx, 7]
            temp_csv_path = os.path.join(
                new_dir, f"Land_Suitability_{forecast_condition}_Communes_LS_values.csv")

            try:
                temp_csv = pd.read_csv(temp_csv_path)

                # Get value for this crop and commune
                value = temp_csv.loc[(temp_csv['Province'] == PROVINCE_NAME) &
                                     (temp_csv['Commune'] == commune),
                                     crop].values[0]

                # Get average value
                avg_avg = pd.read_csv(os.path.join(
                    new_dir, "Land_Suitability_Average_Average_Communes_LS_values.csv"))
                avg_value = avg_avg.loc[(avg_avg['Province'] == PROVINCE_NAME) &
                                        (avg_avg['Commune'] == commune),
                                        crop].values[0]

                # Calculate threshold
                threshold = avg_value - (avg_value * 0.25)

                # Store value in forecast dataframe
                forecast_df.loc[comm_idx, crop] = value

                # Update final forecast dataframe
                mask = ((final_forecast_df['Crop'] == crop) &
                        (final_forecast_df['Commune'] == commune) &
                        (final_forecast_df['Province'] == PROVINCE_NAME))

                # Update values
                final_forecast_df.loc[mask, 'Average_Average'] = avg_value
                final_forecast_df.loc[mask, 'Threshold'] = threshold

                # Set warning if value is below threshold
                month_col = f"M{month_idx+1}"
                if value <= threshold:
                    final_forecast_df.loc[mask, month_col] = "Warning!"
                else:
                    final_forecast_df.loc[mask, month_col] = "No"

                # Update forecast and suitability values
                forecast_col = f"Forecast_M{month_idx+1}"
                suitability_col = f"Planting_Suitability_M{month_idx+1}"
                final_forecast_df.loc[mask, forecast_col] = forecast_condition
                final_forecast_df.loc[mask, suitability_col] = value

                # Set crop start month if applicable
                start_month = all_crops.iloc[crop_idx]['Start_growing_season']
                if month_num == start_month:
                    final_forecast_df.loc[mask, 'Crop_Start_Month'] = month_num

            except Exception as e:
                print(f"Error processing {crop} for {commune}: {e}")

# Process final output if this is the last province
if PROVINCE_NAME == PROVINCES_NAMES[-1]:
    print("Finalizing forecast...")

    # Replace NA with "No Information"
    final_forecast_df['Crop_Start_Month'] = final_forecast_df['Crop_Start_Month'].fillna(
        "No Information")

    # Remove columns that are all NA
    final_forecast_df = final_forecast_df.dropna(axis=1, how='all')

    # Create message column
    final_forecast_df['Message'] = "No Reduced Suitability"
    final_forecast_df.loc[final_forecast_df['Crop_Start_Month']
                          == "No Information", 'Message'] = "No Information"
    final_forecast_df.loc[(final_forecast_df['M1'] == "Warning!") |
                          (final_forecast_df['M2'] == "Warning!") |
                          (final_forecast_df['M3'] == "Warning!"), 'Message'] = "Warning!"

    # Drop Crop_Start_Month column
    final_forecast_df = final_forecast_df.drop('Crop_Start_Month', axis=1)

    # Add reason for warning
    final_forecast_df['dryness'] = ""
    final_forecast_df['temperature'] = ""

    # Extract dryness component from forecast
    for idx, row in final_forecast_df.iterrows():
        if pd.isna(row['Forecast_M1']):
            continue

        forecast = row['Forecast_M1']
        parts = forecast.split('_')

        # Get dryness
        if parts[0] == "Much" and parts[1] == "Drier":
            final_forecast_df.loc[idx, 'dryness'] = "Much Drier"
        elif parts[0] == "Drier":
            final_forecast_df.loc[idx, 'dryness'] = "Drier"
        else:
            final_forecast_df.loc[idx, 'dryness'] = "Average"

        # Get temperature
        if "Much_Warmer" in forecast:
            final_forecast_df.loc[idx, 'temperature'] = "Much Warmer"
        elif "Warmer" in forecast:
            final_forecast_df.loc[idx, 'temperature'] = "Warmer"
        else:
            final_forecast_df.loc[idx, 'temperature'] = "Average"

    # Create reason text
    final_forecast_df['Reason'] = "No Message"

    avg_mask = (final_forecast_df['dryness'] == "Average") & (
        final_forecast_df['temperature'] == "Average")
    final_forecast_df.loc[avg_mask, 'Reason'] = "No Message"

    dry_avg_mask = (final_forecast_df['dryness'] != "Average") & (
        final_forecast_df['temperature'] == "Average")
    final_forecast_df.loc[dry_avg_mask,
                          'Reason'] = final_forecast_df.loc[dry_avg_mask, 'dryness'] + " than normal"

    avg_temp_mask = (final_forecast_df['dryness'] == "Average") & (
        final_forecast_df['temperature'] != "Average")
    final_forecast_df.loc[avg_temp_mask,
                          'Reason'] = final_forecast_df.loc[avg_temp_mask, 'temperature'] + " than normal"

    both_mask = (final_forecast_df['dryness'] != "Average") & (
        final_forecast_df['temperature'] != "Average")
    final_forecast_df.loc[both_mask, 'Reason'] = (final_forecast_df.loc[both_mask, 'dryness'] +
                                                  " and " + final_forecast_df.loc[both_mask, 'temperature'] +
                                                  " than normal")

    # Clean up crop names (remove season info)
    final_forecast_df['Crop'] = final_forecast_df['Crop'].str.split('_').str[0]

    # Create final forecast dataframe
    final_forecast = final_forecast_df[[
        'Province', 'Commune', 'Crop', 'Message', 'Reason']]

    # Save to CSV
    output_path = os.path.join(
        new_dir, f"Forecast_all_communes_{months[0]}_Latest.csv")
    final_forecast.to_csv(output_path, index=False)

    print(f"Saved forecast to {output_path}")

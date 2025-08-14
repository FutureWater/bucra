#!/usr/bin/env python3
"""
Temperature Limits Processing Script
Calculates temperature limits for different crops by:
1) Processing minimum and maximum temperature files.
2) Applying crop-specific base and upper temperature thresholds.
3) Creating temperature suitability maps based on weighted averages
"""

import os
import glob
import calendar
import rasterio
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Import the configuration class
from config import CroptimalConfig
from croptimal_utils import load_reference_dem, fuzzy_membership


def load_cropping_data(config):
    """Load cropping calendar and parameters data."""
    cropping_cal_path = config.scripts_dir / "Cropping_calendar.csv"
    params_path = config.scripts_dir / "Parameters.csv"
    
    if not cropping_cal_path.exists():
        raise FileNotFoundError(f"Cropping calendar not found: {cropping_cal_path}")
    if not params_path.exists():
        raise FileNotFoundError(f"Parameters file not found: {params_path}")
    
    cropping_cal = pd.read_csv(cropping_cal_path, sep=",")
    params = pd.read_csv(params_path)
    
    return cropping_cal, params


def get_temperature_files(config):
    """Get input temperature file paths."""
    temp_dir = config.get_output_path('temperature', '', dir=True)
    
    input_files_tmin = glob.glob(str(temp_dir / "Tmin*" / "**" / "*.tif"), recursive=True)
    input_files_tmax = glob.glob(str(temp_dir / "Tmax*" / "**" / "*.tif"), recursive=True)
    
    if not input_files_tmin:
        raise FileNotFoundError(f"No Tmin files found in {temp_dir}")
    if not input_files_tmax:
        raise FileNotFoundError(f"No Tmax files found in {temp_dir}")
    
    return input_files_tmin, input_files_tmax


def main():
    """Main processing function."""
    print("Processing temperature limits...")

    # Initialize configuration
    config = CroptimalConfig()
    config.validate_inputs()

    print(f"Province: {config.province_name}")

    # Define percentiles to process
    percentiles = ['P75', 'P95', 'Mean']

    # Load cropping calendar and parameters
    cropping_cal, params = load_cropping_data(config)

    # Get input temperature files
    input_files_tmin, input_files_tmax = get_temperature_files(config)

    # Process each percentile
    for percentile in percentiles:
        print(f"Processing temperature limits for {percentile}...")
        
        # Process each crop
        for idx in cropping_cal.index:
            crop_data = cropping_cal.loc[idx]
            crop_name = crop_data["Crop"]
            start_month = int(crop_data["Start_growing_season"])
            end_month = int(crop_data["End_growing_season"])

            # Determine months in growing season of crop
            if start_month > end_month:
                months = list(range(start_month, 13)) + list(range(1, end_month + 1))
            else:
                months = list(range(start_month, end_month + 1))

            # Get weights for each month
            weights = [float(crop_data[f"w_{i+1}"]) for i in range(len(months))]

            # Get temperature thresholds
            t_base = crop_data["T_Base"]
            t_optimal_start = crop_data["T_Optimal_Start"]
            t_optimal_end = crop_data["T_Optimal_End"]
            t_upper = crop_data["T_Upper"]

            # Get month abbreviations for filtering files
            month_abbrs = [calendar.month_abbr[m] for m in months]

            # Filter tmin files for these months with Monthly_Mean
            tmin_files = []
            for month_abbr in month_abbrs:
                tmin_files.extend([f for f in input_files_tmin
                                if month_abbr in os.path.basename(f)])

            # Filter tmax files for these months with specified percentiles
            tmax_files = []
            for month_abbr in month_abbrs:
                tmax_files.extend([f for f in input_files_tmax
                                if month_abbr in os.path.basename(f) and percentile in f])

            # Order files to match months that crop grows in
            tmin_files_ordered = []
            tmax_files_ordered = []
            for month_abbr in month_abbrs:
                for file in tmin_files:
                    if month_abbr in os.path.basename(file):
                        tmin_files_ordered.append(file)
                        break
                for file in tmax_files:
                    if month_abbr in os.path.basename(file):
                        tmax_files_ordered.append(file)
                        break

            if len(tmin_files_ordered) != len(months) or len(tmax_files_ordered) != len(months):
                print(f"  Missing temperature data for {crop_name}, {percentile}")
                # return None

            # Initialize accumulation arrays
            with rasterio.open(tmin_files_ordered[0]) as src:
                shape = src.shape
                tmin_w_data = np.zeros(shape, dtype=np.float32)
                tmax_w_data = np.zeros(shape, dtype=np.float32)
                tmin_valid = np.zeros(shape, dtype=bool)
                tmax_valid = np.zeros(shape, dtype=bool)
                meta = src.meta.copy()

            # Get weighted average of Tmin for whole growing period
            for i, tmin_file in enumerate(tmin_files_ordered):
                with rasterio.open(tmin_file) as src:
                    tmin_data = src.read(1)
                    mask = tmin_data != config.no_data_value
                    tmin_data[mask] *= weights[i]
                    tmin_data[~mask] = 0
                    tmin_w_data += tmin_data
                    tmin_valid |= mask
                        
            tmin_w_data[~tmin_valid] = config.no_data_value

            # Get weighted average of Tmax for whole growing period
            for i, tmax_file in enumerate(tmax_files_ordered):
                with rasterio.open(tmax_file) as src:
                    tmax_data = src.read(1)
                    mask = tmax_data != config.no_data_value
                    tmax_data[mask] *= weights[i]
                    tmax_data[~mask] = 0
                    tmax_w_data += tmax_data
                    tmax_valid |= mask

            tmax_w_data[~tmax_valid] = config.no_data_value

            # Create temperature suitability maps
            if tmin_w_data is not None and tmax_w_data is not None:
                # Apply fuzzy logic
                temperature_lower_limit = fuzzy_membership(tmin_w_data, 'Increasing', [t_base, t_optimal_start])
                temperature_upper_limit = fuzzy_membership(tmax_w_data, 'Decreasing', [t_optimal_end, t_upper])
                temperature_limit = (temperature_lower_limit + temperature_upper_limit) * 0.5
                temperature_limit = temperature_limit.astype(np.float32)
                temperature_limit[~tmax_valid] = config.limit_no_data_value

                # Save output
                season_label = f"{month_abbrs[0]}-{month_abbrs[-1]}"
                folder_name = f"Tmax_{percentile}"
                output_dir = config.get_output_path('suitability', f'Temperature/{folder_name}', dir = True)
                output_path = output_dir / f"Temperature_suitability_{crop_name}_{season_label}.tif"

                # Update metadata for output raster
                meta.update({
                    'dtype': 'float32',
                    'count': 1,
                    'nodata': config.limit_no_data_value,
                    'compress': 'lzw'
                })

                # Write output raster
                with rasterio.open(output_path, 'w', **meta) as dst:
                    dst.write(temperature_limit, 1)

                print(f"  Created temperature suitability map for {crop_name}, {season_label}, {percentile}")

    print("Temperature limits calculation complete!")


if __name__ == "__main__":
    try:
        main()
    except FileNotFoundError as e:
        print(f"Error: Missing required file - {e}")
        exit(1)
    except Exception as e:
        print(f"Error processing temperature limits: {e}")
        exit(1)
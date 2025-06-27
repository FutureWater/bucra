#!/usr/bin/env python3
"""
Land Suitability Processing Script
Produces land suitability maps by:
1) Loading suitability layers from different parameters (NDVI, Temperature, Water, etc.)
2) Applying weightings to each parameter based on its importance
3) Combining all weighted layers to create a final land suitability map
4) Creating variations for different temperature/precipitation conditions
"""

import glob
import calendar
import rasterio
import numpy as np
import pandas as pd

# Import the configuration class
from config import CroptimalConfig


def find_file(file_list, pattern):
    """Find first file matching pattern - only repetitive utility extracted."""
    matching_files = [f for f in file_list if pattern in f]
    return matching_files[0] if matching_files else None


def load_cropping_data(config):
    """Load cropping calendar and parameters data."""
    cropping_cal = pd.read_csv(config.cropping_calendar)
    params_path = config.scripts_dir / "Parameters_fuzzy.csv"
    
    if not params_path.exists():
        raise FileNotFoundError(f"Parameters file not found: {params_path}")
    
    params = pd.read_csv(params_path)
    return cropping_cal, params


def get_suitability_files(config):
    """Get all suitability input files."""
    ls_results_dir = config.get_output_path('suitability', '', dir=True)
    
    files = {
        'ndvi': glob.glob(str(ls_results_dir / "NDVI" / "*.tif")),
        'temperature': glob.glob(str(ls_results_dir / "Temperature" / "**" / "*.tif"), recursive=True),
        'hydraulic': glob.glob(str(ls_results_dir / "Soil_Hydraulic_Properties" / "*.tif")),
        'nutrients': glob.glob(str(ls_results_dir / "Soil_Nutrient_Content" / "*.tif")),
        'slope': glob.glob(str(ls_results_dir / "Elevation" / "*.tif"))
    }
    
    # Validate that we have files
    missing_categories = [k for k, v in files.items() if not v]
    if missing_categories:
        raise FileNotFoundError(f"Missing suitability files for: {missing_categories}")
    
    return files


def main():
    """Main processing function."""
    print("Processing land suitability maps...")
    
    # Initialize configuration
    config = CroptimalConfig()
    config.validate_inputs()
    
    print(f"Province: {config.province_name}")
    
    # Load cropping calendar and parameters
    cropping_cal, params = load_cropping_data(config)
    
    # Get all suitability input files
    files = get_suitability_files(config)
    
    # Create output directory
    output_dir = config.get_output_path('suitability', 'Weighted', dir=True)
    
    # Temperature scenarios
    temp_percentiles = ['75', '95']
    temp_scenario_names = ["Warmer", "Much_Warmer"]
    
    # Process each crop
    for _, crop_row in cropping_cal.iterrows():
        crop = crop_row['Crop']
        start_month = crop_row['Start_growing_season']
        end_month = crop_row['End_growing_season']
        
        # Convert month numbers to abbreviations
        start_month_abbr = calendar.month_abbr[start_month]
        end_month_abbr = calendar.month_abbr[end_month]
        season_label = f"{start_month_abbr}-{end_month_abbr}"
        
        # Process each temperature scenario (including average)
        scenarios = temp_percentiles + [None]  # [75, 95, None]
        
        for t_idx, t_perc in enumerate(scenarios):
            if t_perc is None:
                scenario_name = "Average"
                temp_prefix = "Tmax_Mean"
            else:
                scenario_name = temp_scenario_names[t_idx]
                temp_prefix = f"Tmax_P{t_perc}"
            
            print(f"  Processing {crop} for {scenario_name} conditions...")
            
            # Find temperature file
            temp_pattern = f"{crop}_{season_label}"
            temp_files = [f for f in files['temperature'] if temp_prefix in f and temp_pattern in f]
            if not temp_files:
                print(f"    Temperature file not found for {crop}, {temp_prefix}")
                continue
            
            # Find required parameter files
            ksat_file = find_file(files['hydraulic'], "Ksat")
            wcavail_file = find_file(files['hydraulic'], "WCavail")
            potassium_file = find_file(files['nutrients'], "K_")
            phosphorus_file = find_file(files['nutrients'], "P_")
            slope_file = find_file(files['slope'], "Slope")
            
            # Check if all required files are available
            required_files = [ksat_file, wcavail_file, potassium_file, phosphorus_file, slope_file]
            if not all(required_files):
                print(f"    One or more required files not found for {crop}")
                continue
            
            # Get parameter weights
            ndvi_weight = params.loc[params['Parameter'] == 'NDVI', 'Weight'].values[0]
            temp_weight = params.loc[params['Parameter'] == 'Temperature', 'Weight'].values[0]
            ksat_weight = params.loc[params['Parameter'] == 'Ksat', 'Weight'].values[0]
            wcavail_weight = params.loc[params['Parameter'] == 'WCavail', 'Weight'].values[0]
            potassium_weight = params.loc[params['Parameter'] == 'Potassium', 'Weight'].values[0]
            phosphorus_weight = params.loc[params['Parameter'] == 'Phosphorus', 'Weight'].values[0]
            slope_weight = params.loc[params['Parameter'] == 'Slope', 'Weight'].values[0]
            
            # Read and apply weights to each layer
            with rasterio.open(files['ndvi'][0]) as src:
                no_data_mask = src.read(1) == config.no_data_value
                ndvi_data = src.read(1) * ndvi_weight
                output_meta = src.meta.copy()
            
            with rasterio.open(temp_files[0]) as src:
                temp_data = src.read(1) * temp_weight
            
            with rasterio.open(ksat_file) as src:
                ksat_data = src.read(1) * ksat_weight
            
            with rasterio.open(wcavail_file) as src:
                wcavail_data = src.read(1) * wcavail_weight
            
            with rasterio.open(potassium_file) as src:
                potassium_data = src.read(1) * potassium_weight
            
            with rasterio.open(phosphorus_file) as src:
                phosphorus_data = src.read(1) * phosphorus_weight
            
            with rasterio.open(slope_file) as src:
                slope_data = src.read(1) * slope_weight
            
            # Sum all weighted layers
            suitability_data = (ndvi_data + temp_data + potassium_data + 
                              phosphorus_data + slope_data + ksat_data + wcavail_data)
            
            # Apply no-data mask
            suitability_data[no_data_mask] = config.no_data_value
            
            # Save result
            filename = f"Land_Suitability_{scenario_name}_{crop}_{season_label}.tif"
            output_path = output_dir / filename
            
            output_meta.update({
                'dtype': 'float32',
                'count': 1,
                'nodata': config.no_data_value,
                'compress': 'lzw'
            })
            
            with rasterio.open(output_path, 'w', **output_meta) as dst:
                dst.write(suitability_data.astype(np.float32), 1)
            
            print(f"    Created: {filename}")
    
    print("Land suitability map production complete!")


if __name__ == "__main__":
    try:
        main()
    except FileNotFoundError as e:
        print(f"Error: Missing required file - {e}")
        exit(1)
    except Exception as e:
        print(f"Error processing land suitability: {e}")
        exit(1)
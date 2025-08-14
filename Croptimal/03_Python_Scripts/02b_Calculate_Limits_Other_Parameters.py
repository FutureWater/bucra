#!/usr/bin/env python3
"""
Parameter Limits Processing Script
Calculates suitability limits for various parameters:
1) Soil hydraulic properties (Ksat, WCavail)
2) Soil nutrient content
3) NDVI during growing seasons
4) Slope suitability
"""

import os
import glob
import rasterio
import numpy as np
import pandas as pd

# Import the configuration class
from config import CroptimalConfig
from croptimal_utils import load_reference_dem, fuzzy_membership


def load_parameters_data(config):
    """Load parameters configuration data."""
    params_path = config.scripts_dir / "Parameters_fuzzy.csv"
    
    if not params_path.exists():
        raise FileNotFoundError(f"Parameters file not found: {params_path}")
    
    params = pd.read_csv(params_path)
    return params


def process_slope_suitability(config, folder_params):
    """Process slope suitability limits."""    
    # Get slope file
    slope_files = glob.glob(str(config.slope_output_dir / "*.tif"))
    if not slope_files:
        print("    No slope files found")
        return
    
    slope_file = slope_files[0]
    
    # Get parameters
    limits = folder_params[['Limit1', 'Limit2', 'Limit3', 'Limit4']].iloc[0].values.astype(float)
    limits = limits[~pd.isna(limits)]  # Remove NaN values
    fuzzy_function = folder_params['Function_type'].iloc[0]
    
    # Process slope data
    with rasterio.open(slope_file) as src:
        slope_data = src.read(1)
        slope_meta = src.meta.copy()
        
        # Apply fuzzy membership
        slope_limit = fuzzy_membership(slope_data, fuzzy_function, limits)
        slope_limit[slope_data == config.no_data_value] = config.limit_no_data_value
        
        # Save output
        output_dir = config.get_output_path('suitability', 'Elevation', dir=True)
        output_path = output_dir / "Slope_suitability_score.tif"
        
        out_profile = slope_meta.copy()
        out_profile.update({
            'dtype': 'float32',
            'nodata': config.limit_no_data_value,
            'compress': 'lzw'
        })
        
        with rasterio.open(output_path, 'w', **out_profile) as dst:
            dst.write(slope_limit.astype(np.float32), 1)
        
        print(f"    Created slope suitability map.")


def process_ndvi_suitability(config, folder_params):
    """Process NDVI suitability limits."""    
    # Get NDVI files
    ndvi_dir = config.get_output_path('ndvi', 'Mean_Monthly',dir=True)
    ndvi_files = glob.glob(str(ndvi_dir / "*.tif"))
    
    if not ndvi_files:
        print("    No NDVI files found")
        return
    
    # Get parameters
    limits = folder_params[['Limit1', 'Limit2', 'Limit3', 'Limit4']].iloc[0].values.astype(float)
    limits = limits[~pd.isna(limits)]  # Remove NaN values
    fuzzy_function = folder_params['Function_type'].iloc[0]
    
    # Process all NDVI files to find maximum
    all_ndvi = None
    ndvi_meta = None
    
    for i, ndvi_file in enumerate(ndvi_files):
        with rasterio.open(ndvi_file) as src:
            ndvi_data = src.read(1)
            
            # Initialize arrays on first iteration
            if all_ndvi is None:
                all_ndvi = np.zeros((len(ndvi_files), ndvi_data.shape[0], ndvi_data.shape[1]), dtype=np.float32)
                ndvi_meta = src.meta.copy()
            
            # Add NDVI data to 3D array
            all_ndvi[i] = ndvi_data
    
    # Calculate maximum NDVI per cell for the whole year
    max_ndvi = np.nanmax(all_ndvi, axis=0)
    
    # Apply fuzzy membership
    ndvi_limit = fuzzy_membership(max_ndvi, fuzzy_function, limits)
    ndvi_limit[max_ndvi == config.no_data_value] = config.limit_no_data_value
    
    # Save output
    output_dir = config.get_output_path('suitability', 'NDVI', dir = True)
    output_path = output_dir / "NDVI_suitability_score.tif"
    
    out_profile = ndvi_meta.copy()
    out_profile.update({
        'dtype': 'float32',
        'nodata': config.limit_no_data_value,
        'compress': 'lzw'
    })
    
    with rasterio.open(output_path, 'w', **out_profile) as dst:
        dst.write(ndvi_limit.astype(np.float32), 1)
    
    print(f"    Created NDVI suitability map.")


def process_soil_nutrient_content(config, folder_params):
    """Process soil nutrient content suitability limits."""    
    # Get soil nutrient content files
    snc_dir = config.get_output_path('soil_nutrient_content', '',dir=True)
    snc_files = glob.glob(str(snc_dir / "*.tif"))
    
    if not snc_files:
        print("    No soil nutrient content files found")
        return
    
    # Process each parameter in this folder
    for idx, row in folder_params.iterrows():
        param = row['Parameter']
        abrev = row['Abrev']
        fuzzy_function = row['Function_type']
        
        # Get limits for this parameter
        limits = row[['Limit1', 'Limit2', 'Limit3', 'Limit4']].values.astype(float)
        limits = limits[~pd.isna(limits)]  # Remove NaN values
        
        # Find matching file
        matching_files = [f for f in snc_files if param.lower() in os.path.basename(f).lower()]
        if not matching_files:
            print(f"    No data found for {param}")
            continue
        
        # Process the parameter
        with rasterio.open(matching_files[0]) as src:
            snc_data = src.read(1)
            snc_meta = src.meta.copy()
            
            # Apply fuzzy membership
            snc_limit = fuzzy_membership(snc_data, fuzzy_function, limits)
            snc_limit[snc_data == config.no_data_value] = config.limit_no_data_value
            
            # Save output
            output_dir = config.get_output_path('suitability', 'Soil_Nutrient_Content', dir = True)
            output_path = output_dir / f"{abrev.upper()}_suitability_score.tif"
            
            out_profile = snc_meta.copy()
            out_profile.update({
                'dtype': 'float32',
                'nodata': config.limit_no_data_value,
                'compress': 'lzw'
            })
            
            with rasterio.open(output_path, 'w', **out_profile) as dst:
                dst.write(snc_limit.astype(np.float32), 1)
            
            print(f"    Created {param} suitability map.")


def process_soil_hydraulic_properties(config, folder_params):
    """Process soil hydraulic properties suitability limits."""    
    # Get soil hydraulic properties files
    hhs_dir = config.get_output_path('soil_hydraulic_properties', '', dir = True)
    hhs_files = glob.glob(str(hhs_dir / "*.tif"))
    
    if not hhs_files:
        print("    No soil hydraulic properties files found")
        return
    
    # Process each parameter in this folder
    for idx, row in folder_params.iterrows():
        param = row['Parameter']
        fuzzy_function = row['Function_type']
        
        # Get limits for this parameter
        limits = row[['Limit1', 'Limit2', 'Limit3', 'Limit4']].values.astype(float)
        limits = limits[~pd.isna(limits)]  # Remove NaN values
                
        # Find topsoil and subsoil files
        topsoil_file = None
        subsoil_file = None
        topsoil_data = None
        subsoil_data = None
        topsoil_profile = None
        
        for file in hhs_files:
            if param in file:
                if "TOPSOIL" in file:
                    topsoil_file = file
                    with rasterio.open(topsoil_file) as src:
                        topsoil_data = src.read(1)
                        topsoil_profile = src.profile.copy()
                elif "SUBSOIL" in file:
                    subsoil_file = file
                    with rasterio.open(subsoil_file) as src:
                        subsoil_data = src.read(1)
        
        if topsoil_file is None or subsoil_file is None:
            print(f"    Missing soil data for {param}")
            continue
        
        # Calculate weighted average (0.3 * topsoil + 1.7 * subsoil) / 2
        weighted_topsoil = topsoil_data * 0.3
        weighted_subsoil = subsoil_data * 1.7
        weighted_avg = (weighted_topsoil + weighted_subsoil) / 2
        
        # Apply fuzzy membership
        shp_limit = fuzzy_membership(weighted_avg, fuzzy_function, limits)
        shp_limit[weighted_avg == config.no_data_value] = config.limit_no_data_value
        
        # Save suitability output
        output_dir = config.get_output_path('suitability', 'Soil_Hydraulic_Properties', dir = True)
        output_path = output_dir / f"{param}_suitability_score_{config.province_name}.tif"
        
        out_profile = topsoil_profile.copy()
        out_profile.update({
            'dtype': 'float32',
            'count': 1,
            'nodata': config.limit_no_data_value,
            'compress': 'lzw'
        })
        
        with rasterio.open(output_path, 'w', **out_profile) as dst:
            dst.write(shp_limit.astype(np.float32), 1)
        
        # Save weighted average for reference
        hhs_output_dir = config.get_output_path('soil_hydraulic_properties', '')
        weighted_avg_path = hhs_output_dir / f"weighted_avg_{param}.tif"
        
        with rasterio.open(weighted_avg_path, 'w', **out_profile) as dst:
            dst.write(weighted_avg.astype(np.float32), 1)
        
        print(f"    Created {param} suitability map.")


def main():
    """Main processing function."""
    print("Processing parameter limits...")

    # Initialize configuration
    config = CroptimalConfig()
    config.validate_inputs()

    print(f"Province: {config.province_name}")

    # Load parameters data
    params = load_parameters_data(config)

    # Process each parameter folder
    folders_to_process = params['LS_Results_Folder'].unique()

    for folder in folders_to_process:
        # Skip temperature folder (handled by separate script)
        if folder == "Temperature":
            continue
        
        print(f"Processing {folder}...")
        
        # Get parameters for this folder
        folder_params = params[params['LS_Results_Folder'] == folder]
        
        # Process based on folder type
        if folder == "Elevation":
            process_slope_suitability(config, folder_params)
        elif folder == "NDVI":
            process_ndvi_suitability(config, folder_params)
        elif folder == "Soil_Nutrient_Content":
            process_soil_nutrient_content(config, folder_params)
        elif folder == "Soil_Hydraulic_Properties":
            process_soil_hydraulic_properties(config, folder_params)
        else:
            print(f"    Unknown folder type: {folder}")

    print("Parameter limits calculation complete!")


if __name__ == "__main__":
    try:
        main()
    except FileNotFoundError as e:
        print(f"Error: Missing required file - {e}")
        exit(1)
    except Exception as e:
        print(f"Error processing parameter limits: {e}")
        exit(1)
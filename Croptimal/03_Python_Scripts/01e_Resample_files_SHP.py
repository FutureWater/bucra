#!/usr/bin/env python3
"""
Soil Hydraulic Properties Processing Script
Processes soil hydraulic properties data by resampling to DEM resolution and applying unit conversions.
"""

import glob
import rasterio
import numpy as np
from rasterio.warp import reproject, Resampling
import geopandas as gpd

# Import the configuration class
from config import CroptimalConfig
from croptimal_utils import load_reference_dem


def get_soil_hydraulic_files(config):
    """Find and validate soil hydraulic properties files."""
    # Define variable names and their directories
    var_names = ['WCavail', 'Ksat']
    
    # Define the base path for soil hydraulic data
    base_path = r"C:\Users\thomv\FutureWater Dropbox\Team\Data\Global\Soil\HiHydroSoil_250m\Top_Subsoil"
    
    # Find files for each variable
    input_files = []
    for var in var_names:
        var_files = glob.glob(f"{base_path}/{var}/*.tif")
        input_files.extend(var_files)
    
    if not input_files:
        raise FileNotFoundError(f"No soil hydraulic files found in {base_path}")
    
    return input_files, var_names


def apply_unit_conversion(data, var_name, no_data_value):
    """Apply variable-specific unit conversions to soil hydraulic data."""
    converted_data = data.copy()
    
    if var_name == 'Ksat':
        # Convert from original units to mm/day
        converted_data[converted_data != no_data_value] *= 0.001
    else:
        # Convert WCavail to proper units
        converted_data[converted_data != no_data_value] *= 0.0001
    
    return converted_data

def main():
    """Main processing function."""
    print(f"Processing Soil Hydraulic Properties: {' ' * 5}")
    
    # Initialize configuration
    config = CroptimalConfig()
    config.validate_inputs()
        
    # Load reference DEM
    reference_dem = load_reference_dem(config)
    
    # Load province boundary with buffer
    provinces_shp = gpd.read_file(config.provinces_shapefile)
    province_boundary = provinces_shp[provinces_shp["ADM1_EN"] == config.province_name]
    
    if province_boundary.empty:
        raise ValueError(f"Province '{config.province_name}' not found in shapefile")
    
    province_boundary_reproj = province_boundary.to_crs(reference_dem['crs'])
    province_boundary_buffer = province_boundary_reproj.buffer(20000)
    
    # Get soil hydraulic files
    input_files, var_names = get_soil_hydraulic_files(config)
    
    # Create output directory
    output_dir = config.get_output_path('soil_hydraulic_properties', '', dir=True)
    
    # Process each variable
    output_paths = []
    for var_name in var_names:
        print(f"    Processing {var_name} data...")
        
        # Get files for this variable
        var_files = [file for file in input_files if var_name in file]
        
        
        # Process each file for this variable
        var_results_dict = {}
        for file_path in var_files:
            with rasterio.open(file_path) as src:
                # Calculate window that covers the buffered province extent
                province_window_crs = province_boundary_buffer.to_crs(src.crs)
                window = src.window(*province_window_crs.total_bounds)
                
                # Read only the subset covering the province
                subset_data = src.read(1, window=window)
                subset_transform = src.window_transform(window)
                
                # Create destination array
                destination_array = np.full(
                    (reference_dem['height'], reference_dem['width']), 
                    config.no_data_value, 
                    dtype=np.float32
                )
                
                # Resample to match DEM resolution
                reproject(
                    source=subset_data,
                    destination=destination_array,
                    src_transform=subset_transform,
                    src_crs=src.crs,
                    src_nodata=src.nodata,
                    dst_nodata=config.no_data_value,
                    dst_transform=reference_dem['transform'],
                    dst_crs=reference_dem['crs'],
                    resampling=Resampling.bilinear
                )
                
                # Mask areas where DEM has no data
                destination_array[reference_dem['data'] == reference_dem['nodata']] = config.no_data_value
                
                # Apply unit conversion
                converted_data = apply_unit_conversion(destination_array, var_name, config.no_data_value)
                
                # Get layer name: topsoil or subsoil. And save to results dict.
                layer_name = file_path.split("_")[-1].split(".")[0]  # Extract "Top" or "Sub"
                var_results_dict[layer_name] = converted_data

        # Average top and subsoil layers
        weighted_topsoil = var_results_dict['TOPSOIL'] * 0.3
        weighted_subsoil = var_results_dict['SUBSOIL'] * 1.7
        weighted_avg_data = (weighted_topsoil + weighted_subsoil) / 2.0

        output_name = f"{var_name}_{config.province_name}.tif"
        output_path = output_dir / output_name
        
        # Save processed data of variable
        with rasterio.open(output_path, 'w', **reference_dem['profile']) as dst:
            dst.write(weighted_avg_data.astype(rasterio.float32), 1)


    print("Soil hydraulic properties processing complete!")


if __name__ == "__main__":
    try:
        main()
    except FileNotFoundError as e:
        print(f"Error: Missing required file - {e}")
        exit(1)
    except Exception as e:
        print(f"Error processing soil hydraulic data: {e}")
        exit(1)
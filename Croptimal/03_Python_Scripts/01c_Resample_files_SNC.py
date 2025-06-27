#!/usr/bin/env python3
"""
Soil Nutrient Content (SNC) Processing Script
Processes soil nutrient content data by cropping, resampling to match DEM resolution,
and masking to province boundaries.
"""

import os
import rasterio
import numpy as np
from rasterio.warp import reproject, Resampling
from rasterio.mask import mask
import geopandas as gpd
from pathlib import Path

# Import the configuration class
from config import CroptimalConfig
from croptimal_utils import load_reference_dem


def load_province_boundary(config, target_crs):
    """Load and reproject province boundary shapefile."""
    provinces_shp = gpd.read_file(config.provinces_shapefile)
    province_boundary = provinces_shp[provinces_shp["ADM1_EN"] == config.province_name]
    
    if province_boundary.empty:
        raise ValueError(f"Province '{config.province_name}' not found in shapefile")
    
    return province_boundary.to_crs(target_crs)

def main():
    """Main processing function."""
    # Initialize configuration
    config = CroptimalConfig()
    
    # Validate inputs
    config.validate_inputs()
    
    print(f"Processing SNC files: {config.province_name}")
    
    # Set up paths
    snc_input_dir = config.data_dir / "Soil_Nutrient_Content"
    snc_output_dir = config.get_output_path('soil_nutrient_content', '', dir=True)
    
    # Find input files
    input_files = list(snc_input_dir.glob("*.tif"))
    if not input_files:
        raise FileNotFoundError(f"No SNC files found in {snc_input_dir}")
    
    # Load reference DEM
    reference_dem = load_reference_dem(config)
    
    # Load province boundary
    province_boundary = load_province_boundary(config, reference_dem['crs'])
    
    # Process each SNC file
    for input_file in input_files:
        # Generate output filename
        output_name = input_file.name.replace("0to20cm", config.province_name)
        output_name = f"Extractable_{output_name}"
        output_path = snc_output_dir / output_name
        
        # Process the file
        with rasterio.open(input_file) as src:
            # Crop to province extent first to reduce processing
            province_window_crs = province_boundary.to_crs(src.crs)
            window = src.window(*province_window_crs.total_bounds)
            
            # Read windowed data
            windowed_data = src.read(1, window=window)
            windowed_transform = src.window_transform(window)
            
            # Prepare destination array
            destination_array = np.full(
                (reference_dem['height'], reference_dem['width']), 
                config.no_data_value, 
                dtype=np.float32
            )
            
            # Resample to match DEM resolution
            reproject(
                source=windowed_data,
                destination=destination_array,
                src_transform=windowed_transform,
                src_crs=src.crs,
                dst_transform=reference_dem['transform'],
                dst_crs=reference_dem['crs'],
                resampling=Resampling.bilinear,
                src_nodata=src.nodata,
                dst_nodata=config.no_data_value
            )
            
            # Mask areas where DEM has no data
            dem_mask = reference_dem['data'] == reference_dem['profile']['nodata']
            destination_array[dem_mask] = config.no_data_value
        
            # Save result
            output_profile = reference_dem['profile'].copy()
            output_profile.update({
                'dtype': rasterio.float32,
                'nodata': config.no_data_value,
                'compress': 'lzw'
            })
            
            with rasterio.open(output_path, 'w', **output_profile) as dst:
                dst.write(destination_array.astype(rasterio.float32), 1)
    
    print(f"SNC processing complete for {config.province_name}!\n")


if __name__ == "__main__":
    try:
        main()
    except FileNotFoundError as e:
        print(f"Error: Missing required file - {e}")
        exit(1)
    except Exception as e:
        print(f"Error processing SNC data: {e}")
        exit(1)
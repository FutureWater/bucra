#!/usr/bin/env python3
"""
NDVI Processing Script
Resamples and crops monthly NDVI data to match DEM resolution and extent.
"""

import calendar
import rasterio
import numpy as np
from rasterio.warp import reproject, Resampling

# Import the configuration class
from config import CroptimalConfig
from croptimal_utils import load_reference_dem


def main():
    """Main processing function."""
    print(f"Processing NDVI data: {' ' * 25}")
    
    # Initialize configuration
    config = CroptimalConfig()
    config.validate_inputs()
    
    print(f"Province: {config.province_name}")
    
    # Set up paths
    ndvi_input_path = config.data_dir / "NDVI" / "NDVI_mean_monthly_stack.tif"
    ndvi_output_dir = config.get_output_path('ndvi', 'Mean_Monthly', dir=True)
    
    if not ndvi_input_path.exists():
        raise FileNotFoundError(f"NDVI input file not found: {ndvi_input_path}")
    
    # Load reference DEM using the dedicated function
    reference_dem = load_reference_dem(config)
    
    # Generate month abbreviations
    month_abbrs = [calendar.month_abbr[i] for i in range(1, 13)]
    
    # Process NDVI stack
    with rasterio.open(ndvi_input_path) as src:
        num_bands = src.count
        print(f"Processing {num_bands} monthly NDVI bands...")
        
        for band_idx in range(1, num_bands + 1):
            month = month_abbrs[band_idx - 1]
            print(f"  Processing NDVI for {month}...")
            
            # Read band data
            ndvi_data = src.read(band_idx)
            
            # Create destination array
            resampled_data = np.full(
                (reference_dem['height'], reference_dem['width']), 
                config.no_data_value, 
                dtype=np.float32
            )
            
            # Resample NDVI to match DEM
            reproject(
                source=ndvi_data,
                destination=resampled_data,
                src_transform=src.transform,
                src_crs=src.crs,
                dst_transform=reference_dem['transform'],
                dst_crs=reference_dem['crs'],
                resampling=Resampling.bilinear,
                src_nodata=src.nodata,
                dst_nodata=config.no_data_value
            )
            
            # Mask areas where DEM has no data
            resampled_data[reference_dem['data'] == reference_dem['nodata']] = config.no_data_value
            
            # Save resampled data
            output_path = ndvi_output_dir / f"NDVI_{month}.tif"
            
            output_profile = reference_dem['profile'].copy()
            output_profile.update({
                'dtype': rasterio.float32,
                'nodata': config.no_data_value,
                'compress': 'lzw'
            })
            
            with rasterio.open(output_path, 'w', **output_profile) as dst:
                dst.write(resampled_data.astype(rasterio.float32), 1)
    
    print(f"NDVI processing complete for {config.province_name}!\n")


if __name__ == "__main__":
    try:
        main()
    except FileNotFoundError as e:
        print(f"Error: Missing required file - {e}")
        exit(1)
    except Exception as e:
        print(f"Error processing NDVI data: {e}")
        exit(1)
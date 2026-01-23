#!/usr/bin/env python3
"""
Temperature Processing Script
Processes temperature data by resampling to DEM resolution and applying lapse rate correction.
"""

import glob
import calendar
import rasterio
import numpy as np
from pathlib import Path
from rasterio.warp import reproject, Resampling

# Import the configuration class
from config import CroptimalConfig
from croptimal_utils import load_reference_dem
from IPython import embed

#def main():
"""Main processing function."""
print(f"Processing Temperature data: {' ' * 17}")

# Initialize configuration
config = CroptimalConfig()
config.validate_inputs()
    
# Temperature lapse rate (-0.0065 °C/m)
lapse_rate = config.lapse_rate

# Find temperature input files
input_files = [config.data_dir / "Temperature" / "Tmax.tif", config.data_dir / "Temperature" / "Tmin.tif"]
if not input_files:
    raise FileNotFoundError(f"No temperature files found in {config.data_dir / 'Temperature'}")

# Extract variable names from filenames
temp_vars = [Path(file).stem for file in input_files] 

# Load reference DEM using the dedicated function
reference_dem = load_reference_dem(config)

def main():
    # Process each temperature variable
    for file_name, var in zip(input_files, temp_vars):
        #process_temperature_variable(file_name, var, reference_dem, config, lapse_rate)
        """Process a single temperature variable file: resample and apply lapse rate correction."""

        print(f"Processing {var}")
        
        # Create output directory for this variable
        var_output_dir = config.get_output_path('temperature', var, dir = True)

        # Generate month abbreviations
        month_abbrs = [calendar.month_abbr[i] for i in range(1, 13)]
        
        # Process each month in the temperature file
        with rasterio.open(file_name) as src:
            num_bands = src.count
            
            for band_idx in range(1, num_bands + 1):
                month = month_abbrs[band_idx - 1]
                
                # Read temperature data
                temperature_data = src.read(band_idx)
                
                # Create destination array
                destination_array = np.full(
                    (reference_dem['height'], reference_dem['width']), 
                    config.no_data_value, 
                    dtype=np.float32
                )
                
                # Resample temperature to match DEM
                reproject(
                    source=temperature_data,
                    destination=destination_array,
                    src_transform=src.transform,
                    src_crs=src.crs,
                    dst_transform=reference_dem['transform'],
                    dst_crs=reference_dem['crs'],
                    resampling=Resampling.bilinear,
                    src_nodata=src.nodata,
                    dst_nodata=config.no_data_value
                )
                
                # Mask areas where DEM has no data
                destination_array[reference_dem['data'] == reference_dem['nodata']] = config.no_data_value
                
                # Apply lapse rate correction based on elevation
                valid_mask = destination_array != config.no_data_value
                corrected_data = np.where(
                    valid_mask, 
                    destination_array + reference_dem['data'] * lapse_rate, 
                    config.no_data_value
                )
                
                # Save processed temperature data
                output_path = var_output_dir / f"{var}_{month}_{config.province_name}.tif"
                
                # Update profile with correct data type and nodata value
                output_profile = reference_dem['profile'].copy()
                output_profile.update({
                    'dtype': rasterio.float32,
                    'nodata': config.no_data_value,
                    'compress': 'lzw'
                })
                
                with rasterio.open(output_path, 'w', **output_profile) as dst:
                    dst.write(corrected_data.astype(rasterio.float32), 1)

    print("Temperature data processing complete!")


if __name__ == "__main__":
    try:
        main()
    except FileNotFoundError as e:
        print(f"Error: Missing required file - {e}")
        exit(1)
    except Exception as e:
        print(f"Error processing temperature data: {e}")
        exit(1)
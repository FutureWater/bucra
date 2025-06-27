#!/usr/bin/env python3
"""
DEM and Slope Processing Script
Processes digital elevation models by cropping, resampling, and calculating slope.
Optimized for simplicity and efficiency while maintaining readability.
"""

import numpy as np
import rasterio
from rasterio.warp import reproject, Resampling
from rasterio.mask import mask
import geopandas as gpd
from scipy.ndimage import sobel, distance_transform_edt

# Import the shared configuration class
from config import CroptimalConfig


def load_province_boundary(config):
    """Load and reproject province boundary shapefile."""
    provinces = gpd.read_file(config.provinces_shapefile)
    target_province = provinces[provinces["ADM1_EN"] == config.province_name]
    
    if target_province.empty:
        available = provinces["ADM1_EN"].tolist()
        raise ValueError(f"Province '{config.province_name}' not found. Available: {available}")
    
    return target_province.to_crs(config.local_projection)


def calculate_slope_sobel(dem_array, cell_size_meters, no_data_value):
    """Calculate slope using Sobel edge detection filters."""
    # Prepare data for processing
    valid_mask = dem_array != no_data_value
    dem_processed = dem_array.astype(np.float64)
    
    # Fill no-data areas with nearest valid values to prevent edge artifacts
    invalid_mask = ~valid_mask
    if np.any(invalid_mask):
        indices = distance_transform_edt(invalid_mask, return_distances=False, return_indices=True)
        dem_processed[invalid_mask] = dem_processed[tuple(indices[:, invalid_mask])]
    
    # Calculate gradients using Sobel filters
    dx = sobel(dem_processed, axis=1) / (8 * cell_size_meters)
    dy = sobel(dem_processed, axis=0) / (8 * cell_size_meters)
    
    # Convert gradients to slope percentage
    slope_radians = np.arctan(np.sqrt(dx**2 + dy**2))
    slope_percent = np.tan(slope_radians) * 100
    
    # Restore original no-data mask
    slope_percent[~valid_mask] = no_data_value
    
    return slope_percent


def main():
    """Main processing function."""
    print(f"Processing DEM and Slope: {' ' * 20}")
    
    # Initialize configuration
    config = CroptimalConfig()
    config.validate_inputs()
    
    print(f"Province: {config.province_name}")
    
    # Load province boundary
    province_boundary = load_province_boundary(config)
    
    # Process DEM data
    print("Processing DEM data...")
    
    with rasterio.open(config.dem_file) as src:
        # Step 1: Crop DEM to province extent
        cropped_data, cropped_transform = mask(src, province_boundary.geometry, crop=True)
        
        # Step 2: Set up target grid for resampling
        bounds = province_boundary.total_bounds
        target_width = int((bounds[2] - bounds[0]) / config.resolution_meters)
        target_height = int((bounds[3] - bounds[1]) / config.resolution_meters)
        
        target_transform = rasterio.transform.from_bounds(
            bounds[0], bounds[1], bounds[2], bounds[3], 
            target_width, target_height
        )
        
        # Step 3: Resample to target resolution
        resampled_dem = np.zeros((target_height, target_width), dtype=np.float32)
        
        reproject(
            source=cropped_data[0],
            destination=resampled_dem,
            src_transform=cropped_transform,
            src_crs=src.crs,
            src_nodata=src.nodata or config.no_data_value,
            dst_transform=target_transform,
            dst_crs=config.local_projection,
            dst_nodata=config.no_data_value,
            resampling=Resampling.bilinear
        )
        
        # Step 4: Apply final precise mask
        profile = {
            "driver": "GTiff",
            "height": target_height,
            "width": target_width,
            "count": 1,
            "dtype": resampled_dem.dtype,
            "crs": config.local_projection,
            "transform": target_transform,
            "nodata": config.no_data_value
        }
        
        # Use memory file for final masking
        with rasterio.MemoryFile() as memfile:
            with memfile.open(**profile) as temp_dst:
                temp_dst.write(resampled_dem, 1)
                
                final_data, final_transform = mask(
                    temp_dst,
                    province_boundary.geometry,
                    crop=True,
                    nodata=config.no_data_value
                )
        
        # Update profile for final output
        final_profile = profile.copy()
        final_profile.update({
            "height": final_data.shape[1],
            "width": final_data.shape[2],
            "transform": final_transform
        })
    
    # Calculate slope
    print("Calculating slope...")
    slope_data = calculate_slope_sobel(
        final_data[0], 
        config.resolution_meters, 
        config.no_data_value
    )
    
    # Save results
    print("Saving results...")
    
    # Generate output paths
    dem_filename = f"DEM_{config.province_name}_{config.resolution_meters}m.tif"
    slope_filename = f"Slope_{config.province_name}.tif"
    
    dem_path = config.get_output_path('dem', dem_filename)
    slope_path = config.get_output_path('slope', slope_filename)
    
    # Save DEM
    with rasterio.open(dem_path, "w", **final_profile) as dst:
        dst.write(final_data[0].astype(rasterio.float32), 1)
    
    # Save Slope
    with rasterio.open(slope_path, "w", **final_profile) as dst:
        dst.write(slope_data.astype(rasterio.float32), 1)
    
    print(f"Processing complete for {config.province_name}\n")


if __name__ == "__main__":
    try:
        main()
    except FileNotFoundError as e:
        print(f"Error: Missing required file - {e}")
        exit(1)
    except Exception as e:
        print(f"Error processing DEM data: {e}")
        exit(1)
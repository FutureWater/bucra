#!/usr/bin/env python3
"""
Shared utility functions for Croptimal processing.

This module contains common functions used across multiple processing scripts
to avoid code duplication and ensure consistency.
"""

import numpy as np
import rasterio
import matplotlib.pyplot as plt

def load_reference_dem(config):
    """
    Load DEM data to use as reference for resampling.
    
    Parameters:
    -----------
    config : CroptimalConfig
        Configuration object containing paths and settings
        
    Returns:
    --------
    dict
        Dictionary containing DEM data and metadata:
        - 'data': DEM elevation values
        - 'profile': Rasterio profile for output files
        - 'transform': Geospatial transform
        - 'crs': Coordinate reference system
        - 'height': Raster height in pixels
        - 'width': Raster width in pixels
        - 'nodata': No-data value
    """
    dem_path = config.get_output_path('dem', f"DEM_{config.province_name}_{config.resolution_meters}m.tif")
    
    if not dem_path.exists():
        raise FileNotFoundError(f"Reference DEM not found: {dem_path}")
    
    with rasterio.open(dem_path) as src:
        return {
            'data': src.read(1),
            'profile': src.profile.copy(),
            'transform': src.transform,
            'crs': src.crs,
            'height': src.height,
            'width': src.width,
            'nodata': src.nodata
        }
    

def fuzzy_membership(raster, function_type, parameters):
    """
    Apply fuzzy membership functions to raster data.

    """
    result = raster.copy().astype(float)

    if function_type == "Increasing":
        # Two parameters: [a, b]
        # below a = 0, between a-b = 0 to 1, above b = 1
        a = parameters[0]
        b = parameters[1]
        
        result[raster <= a] = 0
        result[raster >= b] = 1
        mask = (a < raster) & (raster < b)
        result[mask] = (raster[mask] - a) / (b - a)
        
    elif function_type == "Decreasing":
        # Two parameters: [a, b]
        # below a = 1, between a-b = 1 to 0, above b = 0
        a = parameters[0]
        b = parameters[1]
        
        result[raster <= a] = 1
        result[raster >= b] = 0
        mask = (a < raster) & (raster < b)
        result[mask] = 1 - (raster[mask] - a) / (b - a)
    
    elif function_type == "Triangular":
        # Three parameters: [a, b, c]
        # below a = 0, a-b = 0 to 1, b-c = 1 to 0, above c = 0
        a = parameters[0]
        b = parameters[1]  # peak point
        c = parameters[2]

        result[raster <= a] = 0
        result[raster >= c] = 0

        mask1 = (a < raster) & (raster < b)
        mask2 = (b < raster) & (raster < c)
        result[mask1] = (raster[mask1] - a) / (b - a)
        result[mask2] = 1 - (raster[mask2] - b) / (c - b)
    
    elif function_type == "Trapezoidal":
        # Four parameters: [a, b, c, d]
        # below a = 0, a-b = 0 to 1, b-c = 1, c-d = 1 to 0, above d = 0
        a = parameters[0]
        b = parameters[1]  # start of plateau
        c = parameters[2]  # end of plateau
        d = parameters[3]

        result[raster <= a] = 0
        result[raster >= d] = 0

        mask1 = (a < raster) & (raster < b)
        mask2 = (b < raster) & (raster < c)
        mask3 = (c < raster) & (raster < d)
        result[mask1] = (raster[mask1] - a) / (b - a)
        result[mask2] = 1
        result[mask3] = 1 - (raster[mask3] - c) / (d - c)
    
    else:
        raise ValueError(f"Unknown function type: {function_type}. "
                        "Must be one of: 'Increasing', 'Decreasing', 'Triangular', 'Trapezoidal'")
    
    if np.all(result == raster):
        raise ValueError("Result raster is the same as input raster. Fuzzy logic not applied")
    
    return result

def save_suitability_raster(data, template_meta, output_path, no_data_value):
    """Save suitability raster with standardized profile."""
    out_profile = template_meta.copy()
    out_profile.update({
        'dtype': 'float32',
        'nodata': no_data_value,
        'compress': 'lzw'
    })
    
    with rasterio.open(output_path, 'w', **out_profile) as dst:
        dst.write(data.astype(np.float32), 1)

        
def plot_raster(raster, no_data_value, title=None, cmap='viridis', save_path=None):
    """
    Plot raster data with masked no-data values.
    
    Parameters:
    -----------
    raster : numpy.ndarray
        Raster data to plot
    no_data_value : float
        No-data value to mask
    title : str, optional
        Title for the plot
    cmap : str, optional
        Colormap to use (default: 'viridis')
    save_path : str or Path, optional
        Path to save the plot
    """
    plt.figure(figsize=(10, 8))
    masked_raster = np.ma.masked_where(raster == no_data_value, raster)
    im = plt.imshow(masked_raster, cmap=cmap)
    plt.colorbar(im)
    
    if title:
        plt.title(title)
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close()
    else:
        plt.show()


def calculate_weighted_average(raster_files, weights, no_data_value):
    """
    Calculate weighted average of multiple raster files.
    
    Parameters:
    -----------
    raster_files : list
        List of raster file paths
    weights : list
        List of weights corresponding to each raster
    no_data_value : float
        No-data value to handle
        
    Returns:
    --------
    tuple
        (weighted_average_array, validity_mask, metadata)
    """
    import rasterio
    
    if len(raster_files) != len(weights):
        raise ValueError("Number of raster files must match number of weights")
    
    # Initialize with first raster
    with rasterio.open(raster_files[0]) as src:
        shape = src.shape
        weighted_data = np.zeros(shape, dtype=np.float32)
        valid_mask = np.zeros(shape, dtype=bool)
        meta = src.meta.copy()
    
    # Process each raster
    for i, raster_file in enumerate(raster_files):
        with rasterio.open(raster_file) as src:
            data = src.read(1)
            mask = data != no_data_value
            data[mask] *= weights[i]
            data[~mask] = 0
            weighted_data += data
            valid_mask |= mask
    
    # Set invalid cells to no-data value
    weighted_data[~valid_mask] = no_data_value
    
    return weighted_data, valid_mask, meta


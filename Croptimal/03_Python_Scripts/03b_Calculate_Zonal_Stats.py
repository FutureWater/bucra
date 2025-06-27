import glob
import os
from pathlib import Path
import rasterio
import rasterio.mask
import fiona
import geopandas as gpd
import pandas as pd
import numpy as np
from shapely.geometry import shape

from IPython import embed

# Import the centralized configuration
from config import CroptimalConfig

# Setup functions
def load_communes(shapefile_path, province_name):
    """Load and filter communes by province name."""
    with fiona.open(shapefile_path) as shapefile:
        communes = [f for f in shapefile if f['properties']['ADM1_EN'] == province_name]
    
    if not communes:
        raise ValueError(f"No communes found for province: {province_name}")
    
    return communes

def get_variable_name(raster_path):
    """Extract variable name from raster filename."""
    filename = Path(raster_path).stem
    if 'Temperature' in raster_path:
        filename = filename.replace('_', ' ')
    elif 'Land_Suitability' not in raster_path:
        filename = filename.split('_')[0]
    return filename

def calculate_zonal_stats(shapefile_path, raster_folder, nodata_value, province_name, max_communes=None, max_rasters = None):
    """Calculate zonal statistics for raster files within commune boundaries."""
    # Load communes
    communes = load_communes(shapefile_path, province_name)
    
    # Find all raster files
    raster_files = list(Path(raster_folder).rglob('*.tif'))
    if not raster_files:
        raise ValueError(f"No raster files found in {raster_folder}")
    
    # Limit for testing if specified
    if max_communes:
        communes = communes[:max_communes]
    if max_rasters:
        raster_files = raster_files[:max_rasters]

    commune_names = [f['properties']['ADM3_PCODE'] for f in communes]
    commune_geoms = [shape(f['geometry']) for f in communes]

    
    # Initialize results
    results = {'Commune': commune_names}
    
    i = 1
    # Process each raster file
    for raster_path in raster_files:
        variable_name = get_variable_name(str(raster_path))

        # Progress update
        print(f'Processing {variable_name }. Raster {i}/{len(raster_files)}')
        i += 1

        # Open raster files
        with rasterio.open(raster_path) as src:
            commune_means = []
            for geom in commune_geoms:
                try:
                    # Mask and calculate mean
                    out_image, _ = rasterio.mask.mask(src, [geom], crop=True, nodata=nodata_value)
                    valid_pixels = out_image[out_image != nodata_value]
                    
                    mean_val = np.mean(valid_pixels).round(3) if len(valid_pixels) > 0 else np.nan
                    commune_means.append(mean_val)
                    
                except Exception:
                    commune_means.append(np.nan)
            
            results[variable_name] = commune_means
    
    return pd.DataFrame(results)

def reproject_shapefile(input_path, output_path, target_crs):
    """Reproject shapefile to target CRS."""
    gdf = gpd.read_file(input_path)
    gdf_reprojected = gdf.to_crs(target_crs)
    gdf_reprojected.to_file(output_path)
    return output_path

def main():
    # Initialize configuration
    config = CroptimalConfig()
    config.validate_inputs()
    
    print(f"Province: {config.province_name}")
    
    # File paths using config
    suitability_files = config.province_results_dir / 'Suitability'
    communes_file = config.gis_dir / 'egy_admbnda_adm3.shp'
    
    # Reproject communes shapefile
    communes_output_path = config.province_results_dir / 'egy_admbnda_adm3.shp'
    reproject_shapefile(communes_file, communes_output_path, config.local_projection)
    
    # Calculate zonal statistics (limit to 3 communes for testing)
    zonal_stats_df = calculate_zonal_stats(communes_file, 
                                           suitability_files, 
                                           config.no_data_value,
                                           config.province_name,
                                           max_communes=None,
                                           max_rasters= None)
    
    # Save final results to csv
    output_csv = config.province_results_dir / 'Suitability' / 'suitability_per_commune.csv'
    zonal_stats_df.to_csv(output_csv)
    print(f"Processed {len(zonal_stats_df)} communes")

# Main execution
if __name__ == "__main__":
    try:
        main()
    except FileNotFoundError as e:
        print(f"Error: Missing required file - {e}")
        exit(1)
    except Exception as e:
        print(f"Error processing parameter limits: {e}")
        exit(1)
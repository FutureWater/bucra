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


# Import the centralized configuration
from config import CroptimalConfig

# Setup functions
def load_communes(config, shapefile_path, province_name):
    """Load and filter communes by province name."""
    gdf = gpd.read_file(shapefile_path)
    communes_gdf = gdf[gdf['ADM1_EN'] == province_name]
    
    if communes_gdf.empty:
        raise ValueError(f"No communes found for province: {province_name}")
    
    # Reproject
    communes_gdf = communes_gdf.to_crs(config.local_projection)
    return communes_gdf

def get_variable_name(raster_path):
    """Extract variable name from raster filename."""
    # Get the filename without extension
    filename = Path(raster_path).stem

    # Remove the province name
    filename = filename.rsplit('_',1)[0]

    return filename

def calculate_zonal_stats(config, communes_gdf, max_communes=None, max_rasters=None):
    """Calculate zonal statistics for raster files within commune boundaries."""
    # Load data from config
    nodata_value = config.no_data_value
    province_results_dir = config.province_results_dir

    # Find all raster files
    raster_files = list(Path(province_results_dir).rglob('*.tif'))
    if not raster_files:
        raise ValueError(f"No raster files found in {province_results_dir}")
    
    # Limit for testing if specified
    if max_communes:
        communes_gdf = communes_gdf.iloc[:max_communes]
    if max_rasters:
        raster_files = raster_files[:max_rasters]

    commune_names = communes_gdf['ADM3_PCODE'].tolist()
    commune_geoms = communes_gdf.geometry.tolist()
    
    # Initialize results
    results = {'Commune': commune_names}
    
    # Process each raster file
    i = 1
    for raster_path in raster_files:
        variable_name = get_variable_name(str(raster_path))
        # Progress update
        print(f'Processing {variable_name}. Raster {i}/{len(raster_files)}')
        i += 1

        # Open raster files
        with rasterio.open(raster_path) as src:
            commune_means = []
            for idx, geom in enumerate(commune_geoms):
                try:
                    # Mask and calculate mean
                    out_image, _ = rasterio.mask.mask(src, [geom], crop=True, nodata=nodata_value)
                    valid_pixels = out_image[out_image != nodata_value]
                    
                    mean_val = round(np.mean(valid_pixels), 2) if len(valid_pixels) > 0 else np.nan
                    commune_means.append(mean_val)
                    
                except Exception as e:
                    print(f"  - Warning: Could not process commune {commune_names[idx]} for raster {variable_name}.")
                    print(f"  - Error details: {e}")
                    commune_means.append(np.nan)
            
            results[variable_name] = commune_means
    
    # Convert results to DataFrame and join with gdf
    results_df = pd.DataFrame(results).round(2)
    communes_gdf_with_zonal_stats = communes_gdf.merge(results_df, left_on="ADM3_PCODE", right_on="Commune")

    # Save gdf with new variables (for inspection if needed)
    output_dir = config.get_output_path('zonal_stats_communes', '', dir=True)
    output_geopackage = output_dir / (f'communes_with_zonal_stats_{config.province_name}.gpkg')
    communes_gdf_with_zonal_stats.to_file(output_geopackage, driver='GPKG')

    # Print number of cells with no data
    print("\nNumber of missing values per variable:")
    print(results_df.isna().sum())
    print("Total number of missing values:", results_df.isna().sum().sum())
    print("\nPercentage of missing values:")
    print((results_df.isna().sum() / len(results_df) * 100).round(1))
    print("Total percentage of missing values:", (results_df.isna().sum().sum() / (len(results_df) * (len(results_df.columns)-1)) * 100).round(1))

    return results_df

def main():
    # Initialize configuration
    config = CroptimalConfig()
    config.validate_inputs()
    
    print(f"Province: {config.province_name}")
    
    # Reproject communes shapefile
    communes_filepath = config.communes_shapefile
    communes_gdf = load_communes(config, communes_filepath, config.province_name)
    
    # Calculate zonal statistics (limit to 3 communes for testing)
    zonal_stats_df = calculate_zonal_stats(config, 
                                            communes_gdf, 
                                            max_communes=None,
                                            max_rasters= None)
    
    # Save final results to csv
    output_dir = config.get_output_path('zonal_stats_communes', '', dir=True)
    output_csv = output_dir / 'zonal_stats_per_commune.csv'
    zonal_stats_df.to_csv(output_csv)
    print(f"Processed {len(zonal_stats_df)} communes")

# Main execution
if __name__ == "__main__":
    # try:
    main()
    # except FileNotFoundError as e:
    #     print(f"Error: Missing required file - {e}")
    #     exit(1)
    # except Exception as e:
    #     print(f"Error processing parameter limits: {e}")
    #     exit(1)
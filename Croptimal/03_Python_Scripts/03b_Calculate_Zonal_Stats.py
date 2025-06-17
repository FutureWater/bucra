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

# Configuration
PROVINCE_NAME = os.environ.get("PROVINCE", "Sharkia")
NO_DATA_VALUE = 0
LOCAL_PROJECTION = "EPSG:32636"

# Setup paths using pathlib
current_dir = Path.cwd()
parent_dir = current_dir.parent
base_dir = parent_dir.parent.parent

RESULTS_DIR = parent_dir / "04_Results" / PROVINCE_NAME
GIS_DIR = base_dir / "GIS"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)

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
    return filename if 'Land_Suitability' in raster_path else filename.split('_')[0]

def calculate_zonal_stats(shapefile_path, raster_folder, nodata_value, max_communes=None):
    """Calculate zonal statistics for raster files within commune boundaries."""
    # Load communes
    communes = load_communes(shapefile_path, PROVINCE_NAME)
    
    # Limit for testing if specified
    if max_communes:
        communes = communes[:max_communes]
    
    commune_names = [f['properties']['ADM3_PCODE'] for f in communes]
    commune_geoms = [shape(f['geometry']) for f in communes]
    
    # Find all raster files
    raster_files = list(Path(raster_folder).rglob('*.tif'))
    if not raster_files:
        raise ValueError(f"No raster files found in {raster_folder}")
    
    # Initialize results
    results = {'Commune': commune_names}
    
    # Process each raster file
    for raster_path in raster_files:
        variable_name = get_variable_name(str(raster_path))
        
        with rasterio.open(raster_path) as src:
            commune_means = []
            
            for geom in commune_geoms:
                try:
                    # Mask and calculate mean
                    out_image, _ = rasterio.mask.mask(src, [geom], crop=True, nodata=nodata_value)
                    valid_pixels = out_image[out_image != nodata_value]
                    
                    mean_val = np.mean(valid_pixels) if len(valid_pixels) > 0 else np.nan
                    commune_means.append(mean_val)
                    
                except Exception:
                    commune_means.append(np.nan)
            
            results[variable_name] = commune_means
    
    return pd.DataFrame(results)

def reproject_shapefile(input_path, output_path, target_crs=LOCAL_PROJECTION):
    """Reproject shapefile to target CRS."""
    gdf = gpd.read_file(input_path)
    gdf_reprojected = gdf.to_crs(target_crs)
    gdf_reprojected.to_file(output_path)
    return output_path

# Main execution
if __name__ == "__main__":
    # File paths
    suitability_files = RESULTS_DIR / '_LS_Results'
    communes_file = GIS_DIR / 'egy_admbnda_adm2.shp'
    
    # Reproject communes shapefile
    communes_output_path = RESULTS_DIR / 'egy_admbnda_adm3.shp'
    reproject_shapefile(communes_file, communes_output_path)
    
    # Calculate zonal statistics (limit to 3 communes for testing)
    zonal_stats_df = calculate_zonal_stats(communes_file, suitability_files, NO_DATA_VALUE, max_communes=10)
    zonal_stats_df.to_csv(RESULTS_DIR / '_LS_Results' / 'suitability_per_commune.csv')
    print(f"Processed {len(zonal_stats_df)} communes")

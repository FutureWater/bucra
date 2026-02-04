import os
from pathlib import Path
from typing import Optional, Dict, List

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
import rasterio.mask

from config import CroptimalConfig


def load_communes(config: CroptimalConfig, shapefile_path: str, province_name: str) -> gpd.GeoDataFrame:
    """Load and filter communes by province name."""
    gdf = gpd.read_file(shapefile_path)
    communes_gdf = gdf[gdf['ADM1_EN'] == province_name]
    
    if communes_gdf.empty:
        raise ValueError(f"No communes found for province: {province_name}")
    
    return communes_gdf.to_crs(config.local_projection)


def get_variable_name(raster_path: str) -> str:
    """Extract variable name from raster filename."""
    return Path(raster_path).stem.rsplit('_', 1)[0]


def _calculate_commune_mean(src, geom, nodata_value: float) -> float:
    """Calculate mean value for a geometry, excluding nodata."""
    try:
        out_image, _ = rasterio.mask.mask(src, [geom], crop=True, nodata=nodata_value)
        valid_pixels = out_image[out_image != nodata_value]
        return round(np.mean(valid_pixels), 2) if len(valid_pixels) > 0 else np.nan
    except Exception as e:
        return np.nan


def calculate_zonal_stats(
    config: CroptimalConfig,
    communes_gdf: gpd.GeoDataFrame,
    max_communes: Optional[int] = None,
    max_rasters: Optional[int] = None
) -> pd.DataFrame:
    """Calculate zonal statistics for raster files within commune boundaries."""
    province_results_dir = Path(config.province_results_dir)
    raster_files = sorted(province_results_dir.rglob('*.tif'))
    
    if not raster_files:
        raise ValueError(f"No raster files found in {province_results_dir}")
    
    communes_gdf = communes_gdf.iloc[:max_communes] if max_communes else communes_gdf
    raster_files = raster_files[:max_rasters] if max_rasters else raster_files
    
    commune_ids = communes_gdf['ADM3_PCODE'].tolist()
    geoms = communes_gdf.geometry.tolist()
    
    results = {'Commune_ID': commune_ids}
    
    for idx, raster_path in enumerate(raster_files, 1):
        var_name = get_variable_name(str(raster_path))
        print(f'Processing {var_name}. Raster {idx}/{len(raster_files)}')
        
        with rasterio.open(raster_path) as src:
            results[var_name] = [
                _calculate_commune_mean(src, geom, config.no_data_value)
                for geom in geoms
            ]
    
    results_df = pd.DataFrame(results).round(2)
    communes_with_stats = communes_gdf.merge(results_df, left_on="ADM3_PCODE", right_on="Commune_ID")
    
    # Save geopackage
    output_dir = config.get_output_path('zonal_stats_communes', '', dir=True)
    output_gpkg = output_dir / f'communes_with_zonal_stats_{config.province_name}.gpkg'
    communes_with_stats.to_file(output_gpkg, driver='GPKG')
    
    # Report missing values
    missing_values_mask = results_df.isna().sum()
    missing_values = missing_values_mask[missing_values_mask > 0]
    
    if not missing_values.empty:
        total_missing = missing_values.sum()
        total_cells = len(results_df) * (len(results_df.columns) - 1)
        
        print("\nMissing values per variable:")
        print(missing_values)
        print(f"Total missing values: {total_missing}")
        print("\nPercentage of missing values:")
        print((missing_values / len(results_df) * 100).round(1))
        print(f"Total percentage: {(total_missing / total_cells * 100).round(1)}%")

        # Save missing values report
        missing_values_df = missing_values.reset_index(headers=['Variable', 'Missing_Values'])
        missing_values_df["Percentage_Missing"] = (missing_values_df["Missing_Values"] / len(results_df) * 100).round(1)
        missing_values_df.to_csv(output_dir / f'missing_values_report_{config.province_name}.csv', index=False)
    
    return results_df


def main():
    config = CroptimalConfig()
    config.validate_inputs()
    
    print(f"Province: {config.province_name}")
    
    # Load communes
    communes_gdf = load_communes(config, config.communes_shapefile, config.province_name)
    
    # Calculate zonal statistics
    zonal_stats_df = calculate_zonal_stats(config, communes_gdf)
    
    # Save to CSV
    output_dir = config.get_output_path('zonal_stats_communes', '', dir=True)
    output_csv = output_dir / 'zonal_stats_per_commune.csv'
    zonal_stats_df.to_csv(output_csv, index=False)
    print(f"Processed {len(zonal_stats_df)} communes")


if __name__ == "__main__":
    main()
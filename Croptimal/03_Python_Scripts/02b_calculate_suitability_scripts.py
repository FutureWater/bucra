#!/usr/bin/env python3
"""
Crop Suitability Calculation Script
Calculates crop suitability scores for each commune using fuzzy membership functions.
Processes temperature parameters against crop-specific optimal ranges from the cropping calendar.
"""

import calendar
from pathlib import Path

import geopandas as gpd
import pandas as pd

# Import the centralized configuration and utilities
from config import CroptimalConfig
from croptimal_utils import fuzzy_membership


def get_zonal_stats_gdf(config):
    """Load zonal statistics GeoDataFrame for the province."""
    zonal_stats_file = (
        config.province_results_dir 
        / "Zonal_Stats_Communes" 
        / f"communes_with_zonal_stats_{config.province_name}.gpkg"
    )
    
    if not zonal_stats_file.exists():
        raise FileNotFoundError(f"Zonal stats file not found: {zonal_stats_file}")
    
    return gpd.read_file(zonal_stats_file)


def get_cropping_calendar_df(config):
    """Load cropping calendar containing crop parameters and limits."""
    if not config.cropping_calendar.exists():
        raise FileNotFoundError(f"Cropping calendar not found: {config.cropping_calendar}")
    
    return pd.read_csv(config.cropping_calendar)


def get_parameters_df(config):
    """Load parameters file containing fuzzy membership function definitions."""
    if not config.parameter_limits_file.exists():
        raise FileNotFoundError(f"Parameters file not found: {config.parameter_limits_file}")
    
    return pd.read_csv(config.parameter_limits_file)


def get_parameter_list():
    """Generate list of temperature parameters to process."""
    month_abbrs = [calendar.month_abbr[i] for i in range(1, 13)]
    
    tmax_params = [f'Tmax_{month}' for month in month_abbrs]
    tmin_params = [f'Tmin_{month}' for month in month_abbrs]
    
    return tmax_params + tmin_params


def get_crop_suitability_df(config, crop_name):
    """
    Load existing crop suitability DataFrame or create empty one.
    
    Returns:
        tuple: (DataFrame, Path to output CSV)
    """
    output_dir = config.get_output_path('suitability', '', dir=True)
    output_csv = output_dir / f'suitability_{crop_name}.csv'
    
    if output_csv.exists():
        return pd.read_csv(output_csv), output_csv
    
    return pd.DataFrame(), output_csv


def get_parameter_limits(crop_row, parameters_df, parameter):
    """
    Extract parameter limits and membership function type for a crop-parameter combination.
    
    Args:
        crop_row: Row from cropping calendar containing crop-specific limits
        parameters_df: DataFrame with parameter metadata (function types)
        parameter: Parameter name (e.g., 'Tmax_Jan')
    
    Returns:
        dict: Contains 'l1'-'l4' limits and 'function_type'
    """
    # Extract base parameter name (Tmax/Tmin) for lookup
    base_parameter = parameter.split('_', 1)[0] if parameter.startswith(('Tmax', 'Tmin')) else parameter
    
    # Collect limits from cropping calendar
    param_limits = {}
    for limit_type in ['l1', 'l2', 'l3', 'l4']:
        col_name = f"{base_parameter}_{limit_type}"
        if col_name in crop_row.index:
            param_limits[limit_type] = crop_row[col_name]
    
    # Get membership function type from parameters file
    param_row = parameters_df[parameters_df['Parameter'] == base_parameter]
    if param_row.empty:
        raise ValueError(f"Parameter '{base_parameter}' not found in parameters file")
    
    param_limits['function_type'] = param_row.iloc[0]['Function_type']
    
    return param_limits


def calculate_crop_suitability(zonal_stats_gdf, crop_row, parameters_df, parameter_list):
    """
    Calculate suitability scores for all parameters for a single crop.
    
    Args:
        zonal_stats_gdf: GeoDataFrame with zonal statistics per commune
        crop_row: Row from cropping calendar with crop-specific limits
        parameters_df: DataFrame with parameter metadata
        parameter_list: List of parameters to process
    
    Returns:
        DataFrame: Suitability scores for each commune and parameter
    """
    # Extract administrative columns for output
    adm_columns = [col for col in zonal_stats_gdf.columns if col.startswith('ADM')] + ['Commune']
    results = zonal_stats_gdf[adm_columns].copy()
    
    # Calculate suitability for each parameter
    for parameter in parameter_list:
        if parameter not in zonal_stats_gdf.columns:
            raise ValueError(f"Parameter '{parameter}' not found in zonal stats data")
        
        # Get observed values and parameter limits
        observed_values = zonal_stats_gdf[parameter]
        param_limits = get_parameter_limits(crop_row, parameters_df, parameter)
        
        # Apply fuzzy membership function
        results[parameter] = fuzzy_membership(observed_values, param_limits)
    
    return results


def update_crop_suitability_df(existing_df, new_results):
    """
    Merge new suitability results with existing DataFrame, avoiding duplicates.
    
    Args:
        existing_df: Existing crop suitability DataFrame (may be empty)
        new_results: New suitability results to add
    
    Returns:
        tuple: (Updated DataFrame, count of new and updated communes)
    """
    if existing_df.empty:
        return new_results, (len(new_results), 0)
    
    # Update existing DataFrame with communes already present
    existing_communes = existing_df['Commune'].values
    existing_mask = new_results['Commune'].isin(existing_communes)
    existing_df.update(new_results[existing_mask])

    # Find communes not already in existing data
    existing_communes = existing_df['Commune'].values
    new_communes_mask = ~new_results['Commune'].isin(existing_communes)
    new_rows = new_results[new_communes_mask]
    
    if new_rows.empty:
        return existing_df, (0, len(existing_mask))
    
    return pd.concat([existing_df, new_rows], ignore_index=True), (len(new_rows), len(existing_mask))


def main():
    """Main processing function."""
    print(f"Calculating Crop Suitability: {' ' * 15}")
    
    # Initialize configuration
    config = CroptimalConfig()
    config.validate_inputs()
    
    print(f"Province: {config.province_name}")
    
    # Load input data
    zonal_stats_gdf = get_zonal_stats_gdf(config)
    cropping_calendar_df = get_cropping_calendar_df(config)
    parameters_df = get_parameters_df(config)
    parameter_list = get_parameter_list()
    
    print(f"Processing {len(cropping_calendar_df)} crops...")
    
    # Process each crop
    for _, crop_row in cropping_calendar_df.iterrows():
        if crop_row.empty:
            continue
        
        crop_name = crop_row['Crop']
        
        # Load or create suitability DataFrame for this crop
        crop_suitability_df, output_path = get_crop_suitability_df(config, crop_name)
        
        # Calculate suitability scores
        new_results = calculate_crop_suitability(
            zonal_stats_gdf, crop_row, parameters_df, parameter_list
        )
        
        # Update and save results
        crop_suitability_df, new_counts = update_crop_suitability_df(crop_suitability_df, new_results)
        crop_suitability_df.to_csv(output_path, encoding='utf-8-sig', index=False)
        
        if new_counts[0] > 0:
            print(f"  {crop_name}: added {new_counts[0]} communes, updated {new_counts[1]} communes")
        else:
            print(f"  {crop_name}: no new communes to add")
    
    print(f"Crop suitability calculation complete for {config.province_name}!\n")


if __name__ == "__main__":
    try:
        main()
    except FileNotFoundError as e:
        print(f"Error: Missing required file - {e}")
        exit(1)
    except Exception as e:
        print(f"Error calculating crop suitability: {e}")
        exit(1)
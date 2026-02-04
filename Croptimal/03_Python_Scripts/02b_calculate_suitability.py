#!/usr/bin/env python3
"""
Crop Suitability Calculation Script

Calculates crop suitability scores for each commune using fuzzy membership functions.
Processes temperature parameters against crop-specific optimal ranges from the cropping calendar.
"""

import calendar
from pathlib import Path
from typing import Dict, Tuple

import geopandas as gpd
import numpy as np
import pandas as pd

# Import the centralized configuration and utilities
from config import CroptimalConfig
from croptimal_utils import fuzzy_membership


def get_zonal_stats_gdf(config: CroptimalConfig) -> gpd.GeoDataFrame:
    """Load zonal statistics GeoDataFrame for the province."""
    zonal_stats_file = (
        config.province_results_dir 
        / "Zonal_Stats_Communes" 
        / f"communes_with_zonal_stats_{config.province_name}.gpkg"
    )
    
    if not zonal_stats_file.exists():
        raise FileNotFoundError(f"Zonal stats file not found: {zonal_stats_file}")
    
    return gpd.read_file(zonal_stats_file)


def get_cropping_calendar_df(config: CroptimalConfig) -> pd.DataFrame:
    """Load cropping calendar containing crop parameters and limits."""
    if not config.cropping_calendar.exists():
        raise FileNotFoundError(f"Cropping calendar not found: {config.cropping_calendar}")
    
    return pd.read_csv(config.cropping_calendar)


def get_parameters_df(config: CroptimalConfig) -> pd.DataFrame:
    """Load parameters file containing fuzzy membership function definitions."""
    if not config.parameter_limits_file.exists():
        raise FileNotFoundError(f"Parameters file not found: {config.parameter_limits_file}")
    
    return pd.read_csv(config.parameter_limits_file)


def get_parameter_list() -> list[str]:
    """Generate list of temperature parameters to process."""
    month_abbrs = [calendar.month_abbr[i] for i in range(1, 13)]
    
    return (
        [f'Tmax_{month}' for month in month_abbrs]
        + [f'Tmin_{month}' for month in month_abbrs]
    )


def get_crop_suitability_df(config: CroptimalConfig, crop_name: str) -> Tuple[pd.DataFrame, Path]:
    """Load existing crop suitability DataFrame or create empty one."""
    output_dir = config.get_output_path('suitability', '', dir=True)
    output_csv = output_dir / f'suitability_{crop_name}.csv'
    
    if output_csv.exists():
        return pd.read_csv(output_csv), output_csv
    
    return pd.DataFrame(), output_csv


def get_parameter_limits(
    crop_row: pd.Series, 
    parameters_df: pd.DataFrame, 
    parameter: str
) -> Dict[str, any]:
    """Extract parameter limits and membership function type for a crop-parameter combination.
    
    Args:
        crop_row: Row from cropping calendar containing crop-specific limits
        parameters_df: DataFrame with parameter metadata (function types)
        parameter: Parameter name (e.g., 'Tmax_Jan')
    
    Returns:
        Dictionary containing 'l1'-'l4' limits and 'function_type'
    """
    # Extract base parameter name (Tmax/Tmin) for lookup
    base_parameter = parameter.split('_', 1)[0] if parameter.startswith(('Tmax', 'Tmin')) else parameter
    
    # Collect limits from cropping calendar
    param_limits = {
        limit_type: crop_row[f"{base_parameter}_{limit_type}"]
        for limit_type in ['l1', 'l2', 'l3', 'l4']
        if f"{base_parameter}_{limit_type}" in crop_row.index
    }
    
    # Get membership function type from parameters file
    param_row = parameters_df[parameters_df['Parameter'] == base_parameter]
    if param_row.empty:
        raise ValueError(f"Parameter '{base_parameter}' not found in parameters file")
    
    param_limits['function_type'] = param_row.iloc[0]['Function_type']
    
    return param_limits


def get_adm_columns(
    zonal_stats_gdf: gpd.GeoDataFrame, 
    config: CroptimalConfig
) -> pd.DataFrame:
    """Extract administrative columns from zonal stats GeoDataFrame.
    
    Updates ADM csv with new communes.
    """
    output_dir = config.get_output_path('suitability', '', dir=True)
    adm_filepath = output_dir / f'ADM_{config.country_name}.csv'
    
    # Load existing ADM file or create from zonal stats
    columns = (
            [col for col in zonal_stats_gdf.columns if col.startswith('ADM') and col.endswith('EN')]
            + ['Commune_ID', 'ADM3_AR']
        )
    if adm_filepath.exists():
        adm_df = pd.read_csv(adm_filepath)
    else:
        adm_df = zonal_stats_gdf[columns].copy()
    
    # Add new communes if they don't exist yet
    existing_communes = set(adm_df['Commune_ID'].values)
    new_communes_mask = ~zonal_stats_gdf['Commune_ID'].isin(existing_communes)
    
    if new_communes_mask.any():
        new_rows = zonal_stats_gdf.loc[new_communes_mask, columns]
        adm_df = pd.concat([adm_df, new_rows], ignore_index=True)

    # Set country as ADM0_EN if not present
    if 'ADM0_EN' not in adm_df.columns:
        adm_df['ADM0_EN'] = config.country_name

    # Save to csv
    adm_df.to_csv(adm_filepath, encoding='utf-8-sig', index=False)

    return adm_df


def calculate_crop_suitability(
    zonal_stats_gdf: gpd.GeoDataFrame,
    crop_row: pd.Series,
    parameters_df: pd.DataFrame,
    parameter_list: list[str]
) -> pd.DataFrame:
    """Calculate suitability scores for all parameters for a single crop.
    
    Args:
        zonal_stats_gdf: GeoDataFrame with zonal statistics per commune
        crop_row: Row from cropping calendar with crop-specific limits
        parameters_df: DataFrame with parameter values
        parameter_list: List of parameters to process
    
    Returns:
        DataFrame with suitability scores for each commune and parameter
    """
    # Extract commune ids for output
    results = zonal_stats_gdf[['Commune_ID']].copy()
    
    # Calculate suitability for each parameter
    for parameter in parameter_list:
        if parameter not in zonal_stats_gdf.columns:
            raise ValueError(f"Parameter '{parameter}' not found in zonal stats data")
        
        # Get observed values and parameter limits
        observed_values = zonal_stats_gdf[parameter]
        param_limits = get_parameter_limits(crop_row, parameters_df, parameter)
        
        # Calculate suitability value: Apply fuzzy membership function
        results[parameter] = fuzzy_membership(observed_values, param_limits)
    
    return results


def update_crop_suitability_df(
    existing_df: pd.DataFrame, 
    new_results: pd.DataFrame
) -> Tuple[pd.DataFrame, Tuple[int, int]]:
    """Merge new suitability results with existing DataFrame, avoiding duplicates.
    
    Args:
        existing_df: Existing crop suitability DataFrame (may be empty)
        new_results: New suitability results to add
    
    Returns:
        Tuple of (Updated DataFrame, (count of new communes, count of updated communes))
    """
    if existing_df.empty:
        return new_results, (len(new_results), 0)
    
    # Identify which communes already exist
    existing_communes = set(existing_df['Commune_ID'].values)
    existing_mask = new_results['Commune_ID'].isin(existing_communes)
    updated_count = existing_mask.sum()
    
    # Update existing rows by Commune_ID
    if updated_count > 0:
        for idx, row in new_results[existing_mask].iterrows():
            commune_id = row['Commune_ID']
            existing_df.loc[existing_df['Commune_ID'] == commune_id] = row.values

    # Find and add new communes
    new_communes_mask = ~new_results['Commune_ID'].isin(existing_communes)
    new_rows = new_results[new_communes_mask]
    
    if new_rows.empty:
        return existing_df, (0, updated_count)
    
    return pd.concat([existing_df, new_rows], ignore_index=True), (len(new_rows), updated_count)


def update_final_suitability_df(final_results: Dict[str, pd.DataFrame], config: CroptimalConfig) -> None:
    """Create/update a total suitability Excel file with tabs for each crop suitability.
    
    Includes cropping calendar, parameter list, and zonal stats.
    
    Args:
        final_results: Dictionary with all final results DataFrames
        config: CroptimalConfig object
    """
    output_dir = config.get_output_path('suitability', '', dir=True)
    output_xlsx = output_dir / 'total_suitability_results.xlsx'
    
    # Write all sheets to the Excel file
    with pd.ExcelWriter(output_xlsx, engine='openpyxl') as writer:
        for sheet_name, dataframe in final_results.items():
            dataframe.to_excel(writer, sheet_name=sheet_name, index=False)

    print("  Total suitability results excel created/updated.")


def main() -> None:
    """Main processing function."""
    # Initialize configuration
    config = CroptimalConfig()
    config.validate_inputs()
    
    print(f"Calculating Crop Suitability for {config.province_name}")
    print(f"Province: {config.province_name}")
    
    # Load input data
    zonal_stats_gdf = get_zonal_stats_gdf(config)
    cropping_calendar_df = get_cropping_calendar_df(config)
    parameters_df = get_parameters_df(config)
    parameter_list = get_parameter_list()
    
    print(f"Processing {len(cropping_calendar_df)} crops...")

    # Get administrative columns and dataframe with all commune info
    adm_df = get_adm_columns(zonal_stats_gdf, config)

    # Create dictionary to hold final results
    final_results = {
        'zonal_stats': zonal_stats_gdf,
        'cropping_calendar': cropping_calendar_df,
        'parameters': parameters_df,
        'ADM_info': adm_df,
    }

    # Process each crop
    for _, crop_row in cropping_calendar_df.iterrows():
        if crop_row.empty:
            continue
        
        crop_name = crop_row['Crop']
        
        # Load/create suitability DataFrame for this crop
        crop_suitability_df, output_path = get_crop_suitability_df(config, crop_name)
        
        # Calculate suitability scores
        new_results = calculate_crop_suitability(
            zonal_stats_gdf, crop_row, parameters_df, parameter_list
        )
        
        # Update and save results
        crop_suitability_df, (new_count, updated_count) = update_crop_suitability_df(
            crop_suitability_df, new_results
        )
        crop_suitability_df.to_csv(output_path, encoding='utf-8-sig', index=False)

        # Store crop suitability results in final results dictionary
        final_results[crop_name] = crop_suitability_df
        
        # Print summary of communes added/updated
        if new_count > 0:
            print(f"  {crop_name}: added {new_count} communes, updated {updated_count} communes")
        else:
            print(f"  {crop_name}: no new communes. Updated {updated_count} communes.")    

    # Save final results dictionary as Excel file with multiple sheets
    update_final_suitability_df(final_results, config)
    print(f"Crop suitability calculation complete for {config.province_name}!\n")


if __name__ == "__main__":
    main()
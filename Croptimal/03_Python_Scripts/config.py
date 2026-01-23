"""
Centralized configuration for Croptimal processing.

This module should be imported by all processing scripts to ensure
consistent configuration across the entire system.
"""

import os
from pathlib import Path


class CroptimalConfig:
    """
    Centralized configuration class for all Croptimal processing scripts.
    
    This class is defined ONCE in this file and imported by all processing scripts.
    When you need to change configuration logic, you only change it here.
    """
    
    def __init__(self):
        # Get province and country from environment (set by main script) or use default for testing
        self.province_name = os.environ.get("PROVINCE", "Sharkia")
        self.country_name = os.environ.get("COUNTRY", "Egypt")
        
        # Calculate base directory paths
        current_dir = Path.cwd()
        parent_dir = current_dir.parent
        base_dir = parent_dir.parent.parent
        
        # Set up main directories - use environment variables when available
        self.data_dir = Path(os.environ.get("DATA_DIR", parent_dir / "01_Data" / "Nile_Delta"))
        self.results_dir = Path(os.environ.get("RESULTS_DIR", parent_dir / "04_Results"))
        self.temp_dir = Path(os.environ.get("TEMP_DIR", parent_dir / "05_Temp"))
        self.gis_dir = base_dir / "GIS"
        self.scripts_dir = current_dir
        
        # Processing parameters with environment variable overrides
        self.resolution_meters = int(os.environ.get("RESOLUTION", "250"))
        self.local_projection = os.environ.get("LOCAL_PROJ", "EPSG:32636")
        self.no_data_value = -9999.0
        self.limit_no_data_value = 0.0
        self.lapse_rate = -0.0065
        
        # Input file paths
        self.dem_file = self.data_dir / "DEM" / "DEM_NileDelta_250m.tif"
        self.provinces_shapefile = self.gis_dir / "Nile_delta_bnd_adm1.shp"
        self.communes_shapefile = self.gis_dir / "egy_admbnda_adm3.shp"
        self.cropping_calendar = Path(os.environ.get("CROPPING_CALENDER", current_dir / "Cropping_calendar.csv"))
        self.parameter_limits_file = self.scripts_dir / "Parameters_fuzzy.csv"
        
        # Province-specific output directories
        self.country_results_dir = self.results_dir / self.country_name
        self.province_results_dir = self.country_results_dir / self.province_name
        self.dem_output_dir = self.province_results_dir / "DEM"
        self.slope_output_dir = self.province_results_dir / "Slope"
        self.temperature_output_dir = self.province_results_dir / "Temperature"
        self.shp_output_dir = self.province_results_dir / "Soil_Hydraulic_Properties"
        self.snc_output_dir = self.province_results_dir / "Soil_Nutrient_Content"
        self.ndvi_output_dir = self.province_results_dir / "NDVI"
        self.new_communes_shapefile = self.province_results_dir / "Zonal_Stats_Communes"

        # Country-specific output directories
        self.crop_suitability_dir = self.country_results_dir / "Crop_Suitability"
        
        # Ensure all output directories exist
        for output_dir in [self.dem_output_dir, 
                           self.slope_output_dir, 
                          self.temperature_output_dir,
                          self.crop_suitability_dir,
                          self.shp_output_dir,
                          self.snc_output_dir,
                          self.ndvi_output_dir,
                          self.new_communes_shapefile]:
            output_dir.mkdir(parents=True, exist_ok=True)
    
    def validate_inputs(self):
        """Verify that all required input files exist."""
        required_files = {
            'DEM file': self.dem_file,
            'Provinces shapefile': self.provinces_shapefile,
            'Cropping calendar': self.cropping_calendar
        }
        
        missing_files = []
        for description, filepath in required_files.items():
            if not filepath.exists():
                missing_files.append(f"  - {description}: {filepath}")
        
        if missing_files:
            error_message = "Missing required input files:\n" + "\n".join(missing_files)
            raise FileNotFoundError(error_message)
        
        print(f"✓ All required input files found for {self.province_name}")
    
    def get_output_path(self, data_type, filename, dir = False):
        """Generate output file paths with automatic directory creation."""
        output_dirs = {
            'dem': self.dem_output_dir,
            'slope': self.slope_output_dir,
            'temperature': self.temperature_output_dir,
            'suitability': self.crop_suitability_dir,
            'shp': self.shp_output_dir,
            'snc': self.snc_output_dir,
            'ndvi': self.ndvi_output_dir,
            'zonal_stats_communes': self.new_communes_shapefile
        }
        
        # If it is a folder in dict: retrieve it or create it.
        if (data_type.lower() in output_dirs) and dir:
            output_dir = output_dirs[data_type.lower()] / filename
            output_dir.mkdir(parents=True, exist_ok=True)
            return output_dir
        # If it is a file: retrieve file path
        elif data_type.lower() in output_dirs:
            return output_dirs[data_type.lower()] / filename
        # If it is a new folder, outside dict: create it.
        else:
            # Handle new data types by creating appropriate subdirectories
            new_output_dir = self.province_results_dir / data_type.title()
            new_output_dir.mkdir(parents=True, exist_ok=True)
            return new_output_dir / filename
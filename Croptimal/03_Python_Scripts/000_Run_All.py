#!/usr/bin/env python3
"""
Main Orchestration Script
Runs all Croptimal processing scripts for specified provinces in the correct order.
Manages environment variables and handles script execution with proper error handling.
"""

import os
import glob
import datetime
import sys
import subprocess
from pathlib import Path
import geopandas as gpd
from rasterio.env import Env

# Import the configuration class
from config import CroptimalConfig


def get_processing_scripts(scripts_dir):
    """Get list of processing scripts in correct execution order."""
    # Find all Python scripts except this orchestration script
    script_files = glob.glob(str(scripts_dir / "*.py"))
    
    # Exclude orchestration and seasonal forecast scripts
    excluded_patterns = ["000_Run_All.py", "04_Seasonal"]
    processing_scripts = [
        script for script in script_files 
        if not any(pattern in script for pattern in excluded_patterns)
    ]
    
    # Sort to ensure correct execution order
    processing_scripts.sort()
    
    if not processing_scripts:
        raise FileNotFoundError("No processing scripts found")
    
    return processing_scripts


def clean_temp_directory(temp_dir):
    """Clean temporary files from the temp directory."""
    temp_files = glob.glob(str(temp_dir / "*"))
    
    for temp_file in temp_files:
        try:
            os.remove(temp_file)
        except OSError:
            pass  # File might be in use or already deleted


def run_script_for_province(script_path, province_name, config):
    """Run a single processing script for a specific province."""
    script_name = Path(script_path).stem
    print(f"Run {script_name} for {province_name}.")
    
    # Set up environment variables for the subprocess
    script_env = os.environ.copy()
    script_env["PROVINCE"] = province_name
    script_env["RESOLUTION"] = str(config.resolution_meters)
    script_env["LOCAL_PROJ"] = config.local_projection
    script_env["DATA_DIR"] = str(config.data_dir)
    script_env["RESULTS_DIR"] = str(config.results_dir)
    script_env["TEMP_DIR"] = str(config.temp_dir)
    script_env["CROPPING_CALENDER"] = str(config.cropping_calendar)
    
    try:
        # Run script with real-time output
        process = subprocess.Popen(
            ["python", script_path],
            env=script_env,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True
        )
        
        # Print output in real-time
        for stdout_line in iter(process.stdout.readline, ""):
            print(stdout_line, end="")
            sys.stdout.flush()
        
        process.stdout.close()
        return_code = process.wait()
        
        if return_code:
            raise subprocess.CalledProcessError(return_code, process.args)
            
    except subprocess.CalledProcessError as e:
        print(f"Error running script {script_name} for province {province_name}")
        if hasattr(e, 'stderr') and e.stderr:
            print(f"Error details: {e.stderr}")
        raise


def load_province_names(config):
    """Load list of province names from shapefile."""
    provinces_gdf = gpd.read_file(config.provinces_shapefile)
    province_names = provinces_gdf["ADM1_EN"].tolist()
    
    if not province_names:
        raise ValueError("No provinces found in shapefile")
    
    return province_names


def main():
    """Main orchestration function."""
    print("Starting Croptimal processing pipeline...")
    start_time = datetime.datetime.now()
    
    # Initialize configuration
    config = CroptimalConfig()
    
    # Set up raster processing environment
    with Env(
        GDAL_CONFIG_OPTIONS={
            'GDAL_MAX_MEM_ALLOC': '1e+10',
            'GDAL_CACHEMAX': '0.1e+10'
        },
        TEMP_DIR=str(config.temp_dir)
    ):
        # Load province names
        province_names = load_province_names(config)
        print(f"Found {len(province_names)} provinces: {', '.join(province_names)}")
        
        # Get processing scripts
        processing_scripts = get_processing_scripts(config.scripts_dir)
        print(f"Found {len(processing_scripts)} processing scripts")
        
        # Process each province (currently set to last province only for testing)
        for province_name in province_names[:]:
            print(f"\n{'='*60}")
            print(f"Processing province: {province_name}")
            print(f"{'='*60}")
            
            # Clean temp directory before starting
            clean_temp_directory(config.temp_dir)
            
            # Run each processing script for this province
            for script_path in processing_scripts:
                try:
                    run_script_for_province(script_path, province_name, config)
                    print("\n")
                except subprocess.CalledProcessError:
                    print(f"Failed to process {province_name} with {Path(script_path).name}")
                    # Continue with next script rather than stopping entirely
                    continue
                
                # Clean temp files after each script
                clean_temp_directory(config.temp_dir)
        
        # Final cleanup
        clean_temp_directory(config.temp_dir)
    
    # Calculate and display total processing time
    end_time = datetime.datetime.now()
    elapsed_time = end_time - start_time
    print(f"\n{'='*60}")
    print(f"Processing complete!")
    print(f"Total time elapsed: {elapsed_time}")
    print(f"All scripts have been run and temporary files cleaned.")


if __name__ == "__main__":
    try:
        main()
    except FileNotFoundError as e:
        print(f"Error: Missing required file - {e}")
        exit(1)
    except Exception as e:
        print(f"Error in main orchestration: {e}")
        exit(1)
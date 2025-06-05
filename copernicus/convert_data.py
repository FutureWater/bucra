from utils import load_copernicus_data, reshape_projections

##############################
### CONVERT DATA TO NetCDF ###
##############################
# --- SORTED BY DATASET ---

# Reanalysis Hindcast
reanalysis_hindcast_path = "GRIB_data/reanalysis_hindcast_copernicus.grib"
reanalysis_hindcast_netcdf_path = "NetCDF_data/reanalysis_hindcast_copernicus.nc"
reanalysis_hindcast_ds = load_copernicus_data(reanalysis_hindcast_path)
reanalysis_hindcast_ds = reanalysis_hindcast_ds.drop_isel(time=range(12))
reanalysis_hindcast_ds.to_netcdf(reanalysis_hindcast_netcdf_path)

# Reanalysis Validation
reanalysis_validation_path = "GRIB_data/reanalysis_validation_copernicus.grib"
reanalysis_validation_netcdf_path = "NetCDF_data/reanalysis_validation_copernicus.nc"
reanalysis_validation_ds = load_copernicus_data(reanalysis_validation_path)
reanalysis_validation_ds = reanalysis_validation_ds.drop_isel(time=range(12))
reanalysis_validation_ds.to_netcdf(reanalysis_validation_netcdf_path)

# ECMWF Projections Validation
ECMWF_projections_validation_path = "GRIB_data/ECMWF_validation_projections.grib"
ECMWF_projections_validation_netcdf_path = "NetCDF_data/ECMWF_validation_projections.nc"
ECMWF_projections_validation_ds = load_copernicus_data(ECMWF_projections_validation_path)
ECMWF_projections_validation_ds = reshape_projections(ECMWF_projections_validation_ds)
ECMWF_projections_validation_ds = ECMWF_projections_validation_ds.drop_isel(time=range(6))
ECMWF_projections_validation_ds.to_netcdf(ECMWF_projections_validation_netcdf_path)

# ECMWF Hindcast Projections
ECMWF_hindcast_projections_path = "GRIB_data/ECMWF_hindcast_projections.grib"
ECMWF_hindcast_projections_netcdf_path = "NetCDF_data/ECMWF_hindcast_projections.nc"
ECMWF_hindcast_projections_ds = load_copernicus_data(ECMWF_hindcast_projections_path)
ECMWF_hindcast_projections_ds = reshape_projections(ECMWF_hindcast_projections_ds)
ECMWF_hindcast_projections_ds = ECMWF_hindcast_projections_ds.drop_isel(time=range(6))
ECMWF_hindcast_projections_ds.to_netcdf(ECMWF_hindcast_projections_netcdf_path)

# ECCC Projections Validation
projections_validation_path = "GRIB_data/ECCC_validation_projections.grib"
projections_validation_netcdf_path = "NetCDF_data/ECCC_validation_projections.nc"
projections_validation_ds = load_copernicus_data(projections_validation_path)
projections_validation_ds = reshape_projections(projections_validation_ds)
projections_validation_ds = projections_validation_ds.drop_isel(time=range(6))
projections_validation_ds.to_netcdf(projections_validation_netcdf_path)

# ECCC Hindcast Projections
projections_hindcast_path = "GRIB_data/ECCC_hindcast_projections.grib"
projections_hindcast_netcdf_path = "NetCDF_data/ECCC_hindcast_projections.nc"
projections_hindcast_ds = load_copernicus_data(projections_hindcast_path)
projections_hindcast_ds = reshape_projections(projections_hindcast_ds)
projections_hindcast_ds = projections_hindcast_ds.drop_isel(time=range(6))
projections_hindcast_ds.to_netcdf(projections_hindcast_netcdf_path)

from utils import load_copernicus_data, reshape_projections

##############################
### CONVERT DATA TO NetCDF ###
##############################

reanalysis_hindcast_path = "GRIB_data/reanalysis_hindcast_copernicus.grib"
reanalysis_validation_path = "GRIB_data/reanalysis_validation_copernicus.grib"
projections_hindcast_path = "GRIB_data/ECCC_hindcast_projections.grib"
projections_validation_path = "GRIB_data/ECCC_validation_projections.grib"

reanalysis_hindcast_netcdf_path = "NetCDF_data/reanalysis_hindcast_copernicus.nc"
reanalysis_validation_netcdf_path = "NetCDF_data/reanalysis_validation_copernicus.nc"
projections_hindcast_netcdf_path = "NetCDF_data/ECCC_hindcast_projections.nc"
projections_validation_netcdf_path = "NetCDF_data/ECCC_validation_projections.nc"

# Load Data
reanalysis_hindcast_ds = load_copernicus_data(reanalysis_hindcast_path)
reanalysis_validation_ds = load_copernicus_data(reanalysis_validation_path)
projections_hindcast_ds = load_copernicus_data(projections_hindcast_path)
projections_validation_ds = load_copernicus_data(projections_validation_path)

# Reshape the projections to align the lead-time dimension, this drops the first six entries
projections_hindcast_ds = reshape_projections(projections_hindcast_ds)
projections_validation_ds = reshape_projections(projections_validation_ds)

# Drop entries from the reanalysis and projection data to make sure they are the same length
# and every month is represented the same number of times (full years in data).
reanalysis_hindcast_ds = reanalysis_hindcast_ds.drop_isel(time=range(12))
reanalysis_validation_ds = reanalysis_validation_ds.drop_isel(time=range(12))
projections_hindcast_ds = projections_hindcast_ds.drop_isel(time=range(6))
projections_validation_ds = projections_validation_ds.drop_isel(time=range(6))

# Save as NetCDF
reanalysis_hindcast_ds.to_netcdf(reanalysis_hindcast_netcdf_path)
reanalysis_validation_ds.to_netcdf(reanalysis_validation_netcdf_path)
projections_hindcast_ds.to_netcdf(projections_hindcast_netcdf_path)
projections_validation_ds.to_netcdf(projections_validation_netcdf_path)

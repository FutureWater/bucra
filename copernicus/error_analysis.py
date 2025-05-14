from utils import load_copernicus_data

## LOAD DATA
reanalysis_path = r"GRIB_data/reanalysis_copernicus.grib"
projections_path = r"GRIB_data/ECCC_hindcast_projections.grib"

reanalysis_ds = load_copernicus_data(reanalysis_path)
projections_ds = load_copernicus_data(projections_path)

## CALCULATE ERROR PER MONTH

# Collapse ensemble member dimension in projection data
projections_ds = projections_ds.mean(dim="number")

# Interpolate reanalysis_ds to the (courser) dimension of projections_ds
reanalysis_ds = reanalysis_ds.interp_like(projections_ds)

# Calculate the absolute error
abs_error = abs(projections_ds - reanalysis_ds)

# Reduce error dimension to only forecastMonth
abs_error_leadtime = abs_error.mean(dim=("latitude", "longitude", "time"))

abs_error_leadtime.t2m.plot()

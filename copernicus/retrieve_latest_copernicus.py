# %%
import cdsapi
import requests
from datetime import datetime
import xarray as xr
from utils import load_copernicus_data, reshape_projections
import numpy as np
from dateutil.relativedelta import relativedelta
from ibicus.debias import ISIMIP

# Get last month and year
metadata = "https://cds.climate.copernicus.eu/api/catalogue/v1/collections/seasonal-monthly-single-levels"
r = requests.get(metadata)
time_interval = r.json()["extent"]["temporal"]["interval"][0]
end_time = time_interval[1]
dt = datetime.fromisoformat(end_time.replace("Z", "+00:00"))
year = dt.year
month = dt.month
month_str = f"{month:02d}"

# %%
## Download the data
dataset = "seasonal-monthly-single-levels"
request = {
    "originating_centre": "ecmwf",
    "system": "51",
    "variable": [
        "2m_temperature",
        "minimum_2m_temperature_in_the_last_24_hours",
        "maximum_2m_temperature_in_the_last_24_hours",
    ],
    "year": [year],
    "month": [month_str],
    "leadtime_month": ["1", "2", "3", "4", "5", "6"],
    "data_format": "grib",
    "product_type": ["monthly_mean"],
    "area": [31.8, 28.4, 28, 31.4],
}
target = f"GRIB_data/ECMWF_latest_projections_{year}_{month_str}.grib"

# client = cdsapi.Client()
# client.retrieve(dataset, request, target)
# %%
## Convert data to preferred format
ds = load_copernicus_data(target)
ds = ds.mean(dim="number")

# Add an array of months to an array of timestamps, which apparently takes many lines:
array = np.full(
    (ds.forecastMonth.size), datetime.fromisoformat(np.datetime_as_string(ds.time.values))
)
month_range = np.arange(1, ds.forecastMonth.size + 1)

adjusted_array = np.empty_like(array)
for i, months in enumerate(month_range):
    adjusted_array[i] = array[i] + relativedelta(months=months)

ds = ds.assign_coords(predictedTime=("forecastMonth", adjusted_array))


## Bias correction
tas_debiaser = ISIMIP.from_variable("tas")
tas_debiaser.running_window_mode = True
tas_debiaser.detrending = False

reanalysis_hindcast_netcdf_path = "NetCDF_data/reanalysis_hindcast_copernicus.nc"
projections_hindcast_netcdf_path = "NetCDF_data/ECMWF_hindcast_projections.nc"

reanalysis_hindcast_ds = xr.open_dataset(reanalysis_hindcast_netcdf_path)
projections_hindcast_ds = xr.open_dataset(projections_hindcast_netcdf_path)

projections_hindcast_ds = projections_hindcast_ds.mean(dim="number")
reanalysis_hindcast_ds = reanalysis_hindcast_ds.interp_like(projections_hindcast_ds)

reanalysis_t2m = reanalysis_hindcast_ds["t2m"].to_numpy()  # (time, lat, lon)
projections_t2m = projections_hindcast_ds["t2m"].to_numpy()  # (leadtime, time, lat, lon)


def debias_variable(
    var_name, ds, reanalysis_hindcast_ds, projections_hindcast_ds, tas_debiaser
):
    debiased = [None] * ds.predictedTime.size
    for i, time in enumerate(ds.predictedTime.values):
        future = ds[var_name].isel(forecastMonth=i).to_numpy()
        future = np.reshape(future, (1, future.shape[0], future.shape[1]))
        debiased[i] = tas_debiaser.apply(
            reanalysis_hindcast_ds[var_name].to_numpy(),
            projections_hindcast_ds[var_name].isel(forecastMonth=i).to_numpy(),
            future,
            time_obs=reanalysis_hindcast_ds.time.values,
            time_cm_hist=projections_hindcast_ds.time.values,
            time_cm_future=np.array([time]),
        )
    debiased_array = np.array(debiased)
    debiased_xr = xr.DataArray(
        debiased_array.squeeze(),
        dims=("forecastMonth", "latitude", "longitude"),
        coords={
            "forecastMonth": ds.forecastMonth,
            "latitude": ds.coords["latitude"],
            "longitude": ds.coords["longitude"],
            "predictedTime": ds.predictedTime,
        },
        name=f"{var_name}_debiased",
    )
    return debiased_xr


# Debias all variables
t2m_debiased = debias_variable(
    "t2m", ds, reanalysis_hindcast_ds, projections_hindcast_ds, tas_debiaser
)
mn2t24_debiased = debias_variable(
    "mn2t24", ds, reanalysis_hindcast_ds, projections_hindcast_ds, tas_debiaser
)
mx2t24_debiased = debias_variable(
    "mx2t24", ds, reanalysis_hindcast_ds, projections_hindcast_ds, tas_debiaser
)

# # Bias should be corrected per leadtime
# debiased = [None] * ds.predictedTime.size
# for i, time in enumerate(ds.predictedTime.values):
#     future = ds["t2m"].isel(forecastMonth=i).to_numpy()
#     future = np.reshape(future, (1, future.shape[0], future.shape[1]))
#     debiased[i] = tas_debiaser.apply(
#         reanalysis_hindcast_ds["t2m"].to_numpy(),
#         projections_hindcast_ds["t2m"].isel(forecastMonth=i).to_numpy(),
#         future,
#         time_obs=reanalysis_hindcast_ds.time.values,
#         time_cm_hist=projections_hindcast_ds.time.values,
#         time_cm_future=np.array([time]),
#     )

# debiased_array = np.array(debiased)
# debiased_xr = xr.DataArray(
#     debiased_array.squeeze(),
#     dims=("forecastMonth", "latitude", "longitude"),
#     coords={
#         "forecastMonth": ds.forecastMonth,
#         "latitude": ds.coords["latitude"],
#         "longitude": ds.coords["longitude"],
#         "predictedTime": ds.predictedTime,
#     },
#     name="t2m_debiased"
# )


# TODO:
# - Add projections?
# %%

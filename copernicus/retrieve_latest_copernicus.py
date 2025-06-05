import requests
import cdsapi
import xarray as xr
import numpy as np
from datetime import datetime
from utils import load_copernicus_data
from dateutil.relativedelta import relativedelta
from ibicus.debias import QuantileDeltaMapping
from os import path


## Get latest available month from Copernicus metadata
metadata = "https://cds.climate.copernicus.eu/api/catalogue/v1/collections/seasonal-monthly-single-levels"
r = requests.get(metadata)
time_interval = r.json()["extent"]["temporal"]["interval"][0]
end_time = time_interval[1]
dt = datetime.fromisoformat(end_time.replace("Z", "+00:00"))
year = dt.year
month = dt.month
month_str = f"{month:02d}"

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
    "area": [31.7, 24.5, 22, 37],
}
target = f"GRIB_data/ECMWF_latest_projections_{year}_{month_str}.grib"

## (UN)COMMENT THESE TO ENABLE/DISABLE DOWNLOADING
client = cdsapi.Client()
client.retrieve(dataset, request, target)

## Convert data to preferred format
ds = load_copernicus_data(target)

## Average out ensemble members
ds = ds.mean(dim="number")

## Adding lead-time corrected time coordinates to the data array
array = np.full(
    (ds.forecastMonth.size), datetime.fromisoformat(np.datetime_as_string(ds.time.values))
)
month_range = np.arange(ds.forecastMonth.size)
adjusted_array = np.empty_like(array)
for i, months in enumerate(month_range):
    adjusted_array[i] = array[i] + relativedelta(months=months)

ds = ds.assign_coords(predictedTime=("forecastMonth", adjusted_array))


## Bias correction
tas_debiaser = QuantileDeltaMapping.from_variable("tas")

reanalysis_hindcast_netcdf_path = "NetCDF_data/reanalysis_hindcast_copernicus.nc"
projections_hindcast_netcdf_path = "NetCDF_data/ECMWF_hindcast_projections.nc"

reanalysis_hindcast_ds = xr.open_dataset(reanalysis_hindcast_netcdf_path)
projections_hindcast_ds = xr.open_dataset(projections_hindcast_netcdf_path)

projections_hindcast_ds = projections_hindcast_ds.mean(dim="number")
reanalysis_hindcast_ds = reanalysis_hindcast_ds.interp_like(projections_hindcast_ds)


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


## Debiasing
t2m_debiased = debias_variable(
    "t2m", ds, reanalysis_hindcast_ds, projections_hindcast_ds, tas_debiaser
)

## Uncomment to debias the min/max temperatures in the following code once observation data is available ###

# mn2t24_debiased = debias_variable(
#     "mn2t24", ds, reanalysis_hindcast_ds, projections_hindcast_ds, tas_debiaser
# )
# mx2t24_debiased = debias_variable(
#     "mx2t24", ds, reanalysis_hindcast_ds, projections_hindcast_ds, tas_debiaser
# )

# Get mn2t24 and mx2t24 from original ds (instead of debiased data)
mn2t24 = ds["mn2t24"]
mx2t24 = ds["mx2t24"]


## K to C
t2m_debiased -= 273.15
# mn2t24_debiased -= 273.15
# mx2t24_debiased -= 273.15
mn2t24 -= 273.15
mx2t24 -= 273.15

## Writing NetCDF
netcdf_path = "NetCDF_data/"


def write_nc(ds, target, metadata=True):
    if metadata:
        if "time" in ds.coords:
            ds = ds.drop_vars("time")

        # rename the forecastMonth dim to time
        ds = ds.rename({"forecastMonth": "time"})

        # replace the time coordinate values with your actual dates
        ds = ds.assign_coords(time=ds["predictedTime"])
        ds = ds.drop_vars("predictedTime")  # no longer needed as a separate var

        # label the time axis according to CF
        ds.time.attrs.update(
            {
                "standard_name": "time",
                "long_name": "forecast valid time",
                "axis": "T",
            }
        )

        # label the desired projection
        ds = ds.rio.set_spatial_dims(x_dim="longitude", y_dim="latitude")
        ds = ds.rio.write_crs("EPSG:4326")

    # 5. write out
    ds.to_netcdf(target)
    return ds


## To disable writing metadata for GIS-software, add metadata=False argument
write_nc(t2m_debiased, path.join(netcdf_path, "latest_t2m_debiased.nc"))
# write_nc(mn2t24_debiased, path.join(netcdf_path, "latest_mn2t24_debiased.nc"))
# write_nc(mx2t24_debiased, path.join(netcdf_path, "latest_mx2t24_debiased.nc"))
write_nc(mn2t24, path.join(netcdf_path, "latest_mn2t24.nc"))
write_nc(mx2t24, path.join(netcdf_path, "latest_mx2t24.nc"))

# %%

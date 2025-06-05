import xarray as xr

map_path = "NetCDF_data/map.nc"
concencus_path = "NetCDF_data/consensus.nc"
ds = xr.open_dataset(map_path)
concensus = xr.open_dataset(concencus_path)


import cdsapi

dataset = "projections-climate-atlas"
request = {
    "origin": "cordex",
    "experiment": "rcp_8_5",
    "domain": "africa",
    "period": "2006-2100",
    "variable": "monthly_mean_of_daily_mean_temperature",
}

target = "GRIB_data/climate-atlas.nc"
# client = cdsapi.Client()
# client.retrieve(dataset, request, target)

# climate = xr.open_dataset(
#     target,
#     engine="cfgrib",
# )

climate = xr.open_dataset(target)

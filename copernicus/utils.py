import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature


def load_copernicus_data(data_path, st_dim_name="time"):
    """
    Loads the Copernicus dataset using xarray with specified time dimensions.

    Parameters:
    -----------
    data_path : str
        Path to the GRIB file containing the dataset.
    st_dim_name : str, optional
        Name of the time dimension to use. Default is "time".

    Returns:
    --------
    xarray.Dataset
        The loaded dataset.
    """
    return xr.open_dataset(
        data_path,
        engine="cfgrib",
        backend_kwargs=dict(time_dims=("forecastMonth", st_dim_name)),
    )


# Reshape Algorithm
# Start from 7th index, as the 6 months before do not have 6 months of leadtime data
# End date can the month after the last index, but this might be hard to process with datetimes,
# so end at last index.
# Loop through leadtimes (enumerate ds.forecastMonth)
#   Loop through months (from 6th index to last)
#       Copy numlonlat data from (current month - leadtime) in ds and overwrite in copy of ds.


def reshape_projections(data):
    """
    # Reshape Dataset for Forecast Alignment
    ### **! Important !** Removes a number of entries from the start of the data equal to the number of leadtimes! \n
    Aligns forecast months with the corresponding time indices in the dataset, so selecting a month
    and a lead-time (x) will always yield the selected month predicted from x months into the past.
    ## Parameters:
    - **data** (*xarray.Dataset*):
      The input dataset containing forecast data.

    ## Returns:
    - **xarray.Dataset**:
      The reshaped dataset with adjusted time indices.
    """
    # This code cuts off part of the data (and a bit too much).
    # This does not matter as the full first year of data is discarded later regardless.
    ds_reshape = data.copy(deep=True)
    fcsize = data.forecastMonth.size
    for var in data.data_vars:
        for i, fcmonth in enumerate(data.forecastMonth):
            for j in range(len(data.time[fcsize:])):
                ds_reshape[var][:, i, j + fcsize, :, :] = data[var][
                    :, i, j + fcsize - i, :, :
                ]

    ds_reshape = ds_reshape.drop_isel(time=range(0, fcsize))
    return ds_reshape


def plot_temperature_with_borders(ds_reshape, time_idx, leadtime, ensemble_num=None):
    """
    # Plot Temperature with Country Borders

    Plots temperature data with country borders and coastlines using a given dataset.

    ## Parameters:
    - **ds_reshape** (*xarray.Dataset*):
      The reshaped dataset containing temperature data (`t2m`) and associated dimensions.
    - **time_idx** (*int*):
      The index of the time dimension to select for plotting.
    - **leadtime** (*int*):
      The forecast month to select for plotting.
    - **ensemble_num** (*int, optional*):
      The ensemble member number to select. If None, the number dimension is averaged out.

    ## Notes:
    - The function uses the `cartopy` library for map projections and features.
    - The temperature data is plotted in Kelvin with a color bar labeled accordingly.
    - Country borders and coastlines are added to the map for geographical context.

    ## Returns:
    - **None**:
      The function displays the plot but does not return any value.
    """

    if ensemble_num:
        try:
            data = ds_reshape.t2m.sel(number=ensemble_num)
        except Exception:
            print(
                "No 'number' dimension found! Is your input data correct or "
                "is the 'number' dimension already averaged out? "
                "In the last case, leave the 'number' argument 'None'"
            )
    else:
        try:
            data = ds_reshape.t2m.mean(dim="number")
        except Exception:
            print(
                "No 'number' dimension found! Is your input data correct or"
                " is the 'number' dimension already averaged out?"
            )

    plt.figure(figsize=(10, 6))
    ax = plt.axes(projection=ccrs.PlateCarree())
    data.sel(time=ds_reshape.time[time_idx], forecastMonth=leadtime).plot(
        ax=ax, transform=ccrs.PlateCarree(), cbar_kwargs={"label": "Temperature (K)"}
    )
    ax.add_feature(cfeature.BORDERS, linestyle="--", edgecolor="black")
    ax.add_feature(cfeature.COASTLINE, edgecolor="black")
    plt.title("Temperature with Country Borders")
    plt.show()


if __name__ == "__main__":
    data_path = r".\eb1625b16b16f60dd24693f27390bc0e.grib"
    ds = load_copernicus_data(data_path)
    null_count = ds.t2m.isnull().any(dim=("latitude", "longitude")).sum()
    print(null_count.values)
    ds_reshape = reshape_projections(ds)
    plot_temperature_with_borders(ds_reshape, 1, 1)

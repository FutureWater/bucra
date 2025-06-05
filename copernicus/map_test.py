# %%
from utils import load_copernicus_data
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

#################
### LOAD DATA ###
#################
reanalysis_path = "GRIB_data/reanalysis_validation_copernicus.grib"
projections_path = "GRIB_data/ECMWF_hindcast_projections.grib"

reanalysis_ds = load_copernicus_data(reanalysis_path)
projections_ds = load_copernicus_data(projections_path)

#######################
### MAE CALCULATION ###
#######################

# Collapse ensemble member dimension in projection data
projections_ds = projections_ds.mean(dim="number")

# Interpolate reanalysis_ds to the (courser) dimension of projections_ds
reanalysis_ds = reanalysis_ds.interp_like(projections_ds)

# Calculate the absolute error
abs_error = abs(projections_ds - reanalysis_ds)

################
### PLOTTING ###
################
# %%
## MAE per leadtime
# Reduce error dimension to only forecastMonth
abs_error_leadtime = abs_error.mean(dim=("latitude", "longitude", "time"))
abs_error_leadtime.t2m.plot()

# %%
## Spatial MAE for all leadtimes
abs_error_map = abs_error.mean(dim=("time"))


# TODO: Turn into function, feed non-bias-adjusted and bias-adjusted
def plot_spatial_mae(data, suptitle):
    fig, axes = plt.subplots(
        2, 3, figsize=(25, 10), subplot_kw={"projection": ccrs.PlateCarree()}
    )
    forecast_months = range(1, 7)
    fig.suptitle(suptitle, fontsize=30, fontweight="bold")
    for i, forecast_month in enumerate(forecast_months):
        ax = axes[i // 3, i % 3]
        cbar = (
            data.t2m.sel(forecastMonth=forecast_month)
            .plot(
                ax=ax,
                transform=ccrs.PlateCarree(),
                vmin=0,
                vmax=12,
            )
            .colorbar
        )
        cbar.ax.tick_params(labelsize=16)
        cbar.ax.set_ylabel("")
        ax.set_title(f"Lead-time (months): {forecast_month}", fontsize=24)
        ax.add_feature(cfeature.BORDERS, linestyle="--", edgecolor="black")
        ax.add_feature(cfeature.COASTLINE, edgecolor="black")
    plt.show()


# %%

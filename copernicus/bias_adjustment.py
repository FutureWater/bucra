import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from ibicus.debias import QuantileMapping
import xarray as xr
import numpy as np
import pandas as pd
import seaborn as sns


#################
### LOAD DATA ###
#################
reanalysis_hindcast_netcdf_path = "NetCDF_data/reanalysis_hindcast_copernicus.nc"
reanalysis_validation_netcdf_path = "NetCDF_data/reanalysis_validation_copernicus.nc"
projections_hindcast_netcdf_path = "NetCDF_data/ECCC_hindcast_projections.nc"
projections_validation_netcdf_path = "NetCDF_data/ECCC_validation_projections.nc"

reanalysis_hindcast_ds = xr.open_dataset(reanalysis_hindcast_netcdf_path)
reanalysis_validation_ds = xr.open_dataset(reanalysis_validation_netcdf_path)
projections_hindcast_ds = xr.open_dataset(projections_hindcast_netcdf_path)
projections_validation_ds = xr.open_dataset(projections_validation_netcdf_path)

# Collapse ensemble member dimension in projection data
projections_hindcast_ds = projections_hindcast_ds.mean(dim="number")
projections_validation_ds = projections_validation_ds.mean(dim="number")

# Interpolate reanalysis_ds to the (courser) dimension of projections_ds
reanalysis_hindcast_ds = reanalysis_hindcast_ds.interp_like(projections_hindcast_ds)
reanalysis_validation_ds = reanalysis_validation_ds.interp_like(projections_validation_ds)

# Convert to numpy arrays to work with ibicus
reanalysis_hindcast_np = reanalysis_hindcast_ds["t2m"].to_numpy()  # (time, lat, lon)
reanalysis_validation_np = reanalysis_validation_ds["t2m"].to_numpy()  # (time, lat, lon)
projections_hindcast_np = projections_hindcast_ds[
    "t2m"
].to_numpy()  # (leadtime, time, lat, lon)
projections_validation_np = projections_validation_ds[
    "t2m"
].to_numpy()  # (lead-time, time, lat, lon)


#######################
### BIAS ADJUSTMENT ###
#######################
QM_debiaser = QuantileMapping.from_variable("tas")

QM_val = [None] * len(projections_hindcast_ds.forecastMonth)
QM_val = np.empty_like(projections_validation_np)
for i, _ in enumerate(projections_hindcast_ds.forecastMonth):
    QM_val[i] = QM_debiaser.apply(
        reanalysis_hindcast_np, projections_hindcast_np[i], projections_validation_np[i]
    )

# %%
# Bias evaluation
# TODO:
# - Threshold analysis for crop thesholds spatially (per leadtime)
# - Bias evaluation for full dataset spatially (per leadtime)
# - Previous two bias evaluations but done per calendar month so see temporal distribution of bias
#       If a one-dimensional metric can be devised, you can plot the calendar month v leadtime bias in a heatmap

# %%
# TODO:
# - Compare more methods using ibicus evaluation and find the best
# - By hand threshold analysis, apply over corrected and base ds
#       - Per location true false for threshold,
#       - Select cell with farm: sum per calendar month, show result in heatmap leadtime v month)
# - Plot like above but for MAE throughout the year for different leadtimes

# %%
######################
### ERROR ANALYSIS ###
######################


def calculate_threshold_exceedances_per_leadtime_month(df, threshold):
    threshold_map = (df > threshold[0]) & (df < threshold[1])
    data = np.zeros((6, 12))
    for i, _ in enumerate(threshold_map[:, 0, 0, 0]):
        for j, _ in enumerate(threshold_map[0, :, 0, 0]):
            data[i, j % 12] = data[i, j % 12] + threshold_map[i, j].sum()
    return data


def calculate_threshold_exceedances_per_month(obs_df, threshold):
    threshold_map = (obs_df > threshold[0]) & (obs_df < threshold[1])
    # obs_df shape: (time, lat, lon)
    data = np.zeros(12)
    for i in range(obs_df.shape[0]):
        data[i % 12] += threshold_map[i].sum()
    return data


# Create a DataFrame for the data
def plot_heatmap(data, title, colorbar_label, decimals=1):
    df = pd.DataFrame(
        data,
        columns=[
            "Jan",
            "Feb",
            "Mar",
            "Apr",
            "May",
            "Jun",
            "Jul",
            "Aug",
            "Sep",
            "Oct",
            "Nov",
            "Dec",
        ],
    )
    df.index = [f"Leadtime {i + 1}" for i in range(data.shape[0])]

    plt.figure(figsize=(10, 6))
    sns.heatmap(
        df,
        annot=True,
        fmt=f".{decimals}f",
        cmap="coolwarm",
        cbar_kws={"label": colorbar_label},
    )
    plt.title(title)
    plt.xlabel("Month")
    plt.ylabel("Leadtime")
    plt.show()


sample_threshold = [293, 308]
QM_threshold = calculate_threshold_exceedances_per_leadtime_month(
    QM_val, sample_threshold
)
raw_threshold = calculate_threshold_exceedances_per_leadtime_month(
    projections_validation_np, sample_threshold
)
obs_threshold = calculate_threshold_exceedances_per_month(
    reanalysis_validation_np, sample_threshold
)
diff_threshold = QM_threshold - obs_threshold

plot_heatmap(
    diff_threshold,
    "Threshold Exceedance differences (Projected - Observed) per Leadtime and Month",
    "Threshold Exceedances",
    0,
)
# %%


## MAE per leadtime per month
def calculate_MAE_per_leadtime_month(df, obs):
    AE_map = abs(df - obs)
    data = np.zeros((6, 12))
    for i, _ in enumerate(AE_map[:, 0, 0, 0]):
        for j, _ in enumerate(AE_map[0, :, 0, 0]):
            data[i, j % 12] = data[i, j % 12] + AE_map[i, j, 1, 3]
    # Since means can be summed over multiple years, divide by number of years in the data
    data = data / (obs.shape[0] / 12)
    return data


MAE_leadtime_month = calculate_MAE_per_leadtime_month(
    projections_validation_np, reanalysis_validation_np
)
plot_heatmap(MAE_leadtime_month, "Raw AE per Leadtime and Month for Qabunah", "AE (K)")

MAE_leadtime_month = calculate_MAE_per_leadtime_month(QM_val, reanalysis_validation_np)
plot_heatmap(MAE_leadtime_month, "QM AE per Leadtime and Month for Qabunah", "AE (K)")


# %%
# MAE Maps
## Spatial MAE for all leadtimes


def plot_spatial_mae(data, suptitle):
    fig, axes = plt.subplots(
        2, 3, figsize=(25, 10), subplot_kw={"projection": ccrs.PlateCarree()}
    )
    forecast_months = range(1, 7)
    fig.suptitle(suptitle, fontsize=30, fontweight="bold")
    for i, forecast_month in enumerate(forecast_months):
        ax = axes[i // 3, i % 3]
        cbar = (
            data.sel(forecastMonth=forecast_month)
            .plot(
                ax=ax,
                transform=ccrs.PlateCarree(),
                vmin=0,
                vmax=5,
            )
            .colorbar
        )
        cbar.ax.tick_params(labelsize=16)
        cbar.ax.set_ylabel("")
        ax.set_title(f"Lead-time (months): {forecast_month}", fontsize=24)
        ax.add_feature(cfeature.BORDERS, linestyle="--", edgecolor="black")
        ax.add_feature(cfeature.COASTLINE, edgecolor="black")
    plt.show()


abs_error = abs(QM_val - reanalysis_validation_np)
abs_error_map = abs_error.mean(axis=(1))
abs_error_map = xr.DataArray(
    abs_error_map,
    dims=("forecastMonth", "latitude", "longitude"),
    coords={
        "forecastMonth": projections_validation_ds.forecastMonth,
        "latitude": projections_validation_ds.latitude,
        "longitude": projections_validation_ds.longitude,
    },
    name="t2m",
)

abs_error_raw = abs(projections_validation_np - reanalysis_validation_np)
abs_error_raw_map = abs_error_raw.mean(axis=(1))
abs_error_raw_map = xr.DataArray(
    abs_error_raw_map,
    dims=("forecastMonth", "latitude", "longitude"),
    coords={
        "forecastMonth": projections_validation_ds.forecastMonth,
        "latitude": projections_validation_ds.latitude,
        "longitude": projections_validation_ds.longitude,
    },
    name="t2m",
)

plot_spatial_mae(abs_error_raw_map, "Raw MAE")
plot_spatial_mae(abs_error_map, "Bias corrected MAE")

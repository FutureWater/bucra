# %%
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from ibicus.debias import QuantileMapping, ISIMIP, QuantileDeltaMapping
import xarray as xr
import numpy as np
import pandas as pd
import seaborn as sns


#################
### LOAD DATA ###
#################
reanalysis_hindcast_netcdf_path = "NetCDF_data/reanalysis_hindcast_copernicus.nc"
reanalysis_validation_netcdf_path = "NetCDF_data/reanalysis_validation_copernicus.nc"
projections_hindcast_netcdf_path = "NetCDF_data/ECMWF_hindcast_projections.nc"
projections_validation_netcdf_path = "NetCDF_data/ECMWF_validation_projections.nc"

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
def apply_debiaser(debiaser):
    output = np.empty_like(projections_validation_np)
    for i, _ in enumerate(projections_hindcast_ds.forecastMonth):
        output[i] = debiaser.apply(
            reanalysis_hindcast_np,
            projections_hindcast_np[i],
            projections_validation_np[i],
        )
    return output


# Quartile Mapping adjustment (good base)
QM_debiaser = QuantileMapping.from_variable("tas")
QM_val = apply_debiaser(QM_debiaser)


# ISIMIP (advanced method)
ISIMIP_debiaser = ISIMIP.from_variable("tas")
ISIMIP_val = apply_debiaser(ISIMIP_debiaser)

# QuantileDeltaMapping (advanced version of ECDFM, uses a running window)
QDM_debiaser = QuantileDeltaMapping.from_variable("tas")
QDM_val = apply_debiaser(QDM_debiaser)


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
            data[i, j % 12] = data[i, j % 12] + threshold_map[i, j, 1, 3]
    return data


def calculate_threshold_exceedances_per_month(obs_df, threshold):
    threshold_map = (obs_df > threshold[0]) & (obs_df < threshold[1])
    # obs_df shape: (time, lat, lon)
    data = np.zeros(12)
    for i in range(obs_df.shape[0]):
        data[i % 12] += threshold_map[i, 1, 3]
    return data


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


sample_threshold = [294, 313]
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
def calculate_MAE_per_leadtime_month(df):
    AE_map = abs(df - reanalysis_validation_np)
    data = np.zeros((6, 12))
    for i, _ in enumerate(AE_map[:, 0, 0, 0]):
        for j, _ in enumerate(AE_map[0, :, 0, 0]):
            data[i, j % 12] = data[i, j % 12] + AE_map[i, j, 1, 3]
    # Since means can be summed over multiple years, divide by number of years in the data
    data = data / (reanalysis_validation_np.shape[0] / 12)
    return data


MAE_leadtime_month = calculate_MAE_per_leadtime_month(projections_validation_np)
plot_heatmap(MAE_leadtime_month, "Raw AE per Leadtime and Month for Qabunah", "AE (K)")

MAE_leadtime_month = calculate_MAE_per_leadtime_month(QM_val)
plot_heatmap(MAE_leadtime_month, "QM AE per Leadtime and Month for Qabunah", "AE (K)")

MAE_leadtime_month = calculate_MAE_per_leadtime_month(ISIMIP_val)
plot_heatmap(MAE_leadtime_month, "ISIMIP AE per Leadtime and Month for Qabunah", "AE (K)")

MAE_leadtime_month = calculate_MAE_per_leadtime_month(QDM_val)
plot_heatmap(MAE_leadtime_month, "QDM AE per Leadtime and Month for Qabunah", "AE (K)")
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


def abs_error_map(debiased_val):
    abs_error = abs(debiased_val - reanalysis_validation_np)
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
    return abs_error_map


raw_abs_error_map = abs_error_map(projections_validation_np)
QM_abs_error_map = abs_error_map(QM_val)
ISIMIP_abs_error_map = abs_error_map(ISIMIP_val)
QDM_abs_error_map = abs_error_map(QDM_val)

plot_spatial_mae(raw_abs_error_map, "Raw MAE")
plot_spatial_mae(QM_abs_error_map, "QM Bias corrected MAE")
plot_spatial_mae(ISIMIP_abs_error_map, "ISIMIP Bias corrected MAE")
plot_spatial_mae(QDM_abs_error_map, "QDM Bias corrected MAE")

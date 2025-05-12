# %%
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# %%
# LOAD DATA
data_path = ".\eb1625b16b16f60dd24693f27390bc0e.grib"
ds = xr.open_dataset(data_path, engine="cfgrib")
# %%
# PRELIMINARY PLOT

# INVESTIGATE MISSING DATA
missing_data = ds.sst.isnull()
missing_count = missing_data.sum().item()
print(f"Number of missing data points: {missing_count}")

# Plot missing data
plt.figure(figsize=(10, 6))
missing_data.mean(axis=0).plot()
plt.title("Missing Data Distribution")
plt.show()
# %%
# TODO: Figure out if leadtime is time before month, or time after...
ds.t2m[0, -1, 1].plot()
# %%
# Plot with country borders
plt.figure(figsize=(10, 6))
ax = plt.axes(projection=ccrs.PlateCarree())
ds.t2m[0, -1, 2].plot(
    ax=ax, transform=ccrs.PlateCarree(), cbar_kwargs={"label": "Temperature (K)"}
)
ax.add_feature(cfeature.BORDERS, linestyle="--", edgecolor="black")
ax.add_feature(cfeature.COASTLINE, edgecolor="black")
plt.title("Temperature with Country Borders")
plt.show()

# %%
# FIND COMBINATIONS WITH NO MISSING VALUES
valid_combinations = ds.sst.where(~ds.sst.isnull(), drop=True)
unique_combinations = valid_combinations.groupby(["time", "number", "step"]).count()

# FILTER COMBINATIONS WITH NO MISSING VALUES
no_missing_combinations = unique_combinations.where(
    unique_combinations
    == valid_combinations.sizes["latitude"] * valid_combinations.sizes["longitude"],
    drop=True,
)

print("Number/Timeframe/Leadtime combinations with no missing values:")
print(no_missing_combinations)

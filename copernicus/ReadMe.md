# Copernicus Bias Adjustment & Evaluation Toolkit

This repository provides a workflow for downloading, converting, bias-adjusting, and evaluating Copernicus seasonal forecast and reanalysis data, with a focus on temperature (2m temperature, `t2m`). The workflow is designed for agricultural applications, such as crop threshold analysis and spatial/temporal error evaluation.

---

## Workflow Overview

1. **Download Data**  
   Download reanalysis and forecast (projections) data for both hindcast and validation periods from Copernicus using the provided scripts.

2. **Convert Data to NetCDF**  
   Convert the downloaded GRIB files to NetCDF format and align/reshape the datasets for analysis.

3. **Bias Adjustment**  
   Apply bias adjustment methods (Quantile Mapping, ISIMIP, QuantileDeltaMapping) to the forecast data.

4. **Error & Bias Evaluation**  
   Evaluate the performance of bias adjustment using threshold exceedance analysis, MAE (Mean Absolute Error), and spatial/temporal plots.

---

## File Descriptions & Usage Order

### 1. Download Data

- **`download_hindcast_reanalysis_copernicus.py`**  
  Downloads ERA5 reanalysis data for the hindcast period (historical years).

- **`download_validation_reanalysis_copernicus.py`**  
  Downloads ERA5 reanalysis data for the validation period (recent years).

- **`download_hindcast_projections_copernicus.py`**  
  Downloads ECMWF seasonal forecast (projections) for the hindcast period.

- **`download_validation_projections_copernicus.py`**  
  Downloads ECMWF seasonal forecast (projections) for the validation period.

- **`retrieve_latest_copernicus.py`**  
  Downloads the latest available ECMWF seasonal forecast and processes the time coordinates.

**Usage:**  
Run these scripts to populate the `GRIB_data/` directory with the required `.grib` files.

---

### 2. Convert Data to NetCDF

- **`convert_data.py`**  
  Loads the downloaded GRIB files, reshapes the forecast data for leadtime alignment, and saves all datasets as NetCDF files in `NetCDF_data/`.

- **`utils.py`**  
  Contains utility functions for loading Copernicus data (`load_copernicus_data`), reshaping projections (`reshape_projections`), and plotting.

**Usage:**  
Run `convert_data.py` after downloading all GRIB files.  
You may use functions from `utils.py` in your own scripts for data loading and visualization.

---

### 3. Bias Adjustment

- **`bias_adjustment.py`**  
  Loads the NetCDF datasets, applies bias adjustment methods (using the `ibicus` library), and prepares bias-adjusted forecast arrays.  
  Also includes functions for threshold analysis, MAE calculation, and plotting (heatmaps, spatial maps).

**Usage:**  
Run or adapt `bias_adjustment.py` to perform bias correction and analyze results.  
Edit threshold values or analysis parameters as needed for your application.

---

### 4. Error & Bias Evaluation

- **`bias_adjustment.py`**  
  (Continued) Contains code for:
  - Threshold exceedance analysis (per leadtime and month)
  - MAE calculation (per leadtime, month, and spatially)
  - Visualization of results

- **`error_analysis.py`**  
  Example script for calculating and plotting absolute errors between projections and reanalysis data.

- **`legacy_bias_eval.py`**  
  Example of using custom metrics and marginal bias analysis (legacy code, may require adaptation).

---

## Typical Usage Order

1. **Download all required data**  
   Run the four `download_*.py` scripts and `retrieve_latest_copernicus.py` as needed.

2. **Convert and align data**  
   Run `convert_data.py` to produce NetCDF files.

3. **Bias adjustment and evaluation**  
   Run `bias_adjustment.py` to apply bias correction and generate evaluation plots.

4. **Further analysis**  
   Use `error_analysis.py` and `legacy_bias_eval.py` for additional or custom analyses.

---

## Dependencies

- Python 3.x
- [xarray](https://xarray.pydata.org/)
- [cdsapi](https://cds.climate.copernicus.eu/api-how-to)
- [cartopy](https://scitools.org.uk/cartopy/docs/latest/)
- [matplotlib](https://matplotlib.org/)
- [seaborn](https://seaborn.pydata.org/)
- [ibicus](https://github.com/ibicus-org/ibicus) (for bias adjustment)
- [pandas](https://pandas.pydata.org/)
- [numpy](https://numpy.org/)
- [requests](https://docs.python-requests.org/)
- [dateutil](https://dateutil.readthedocs.io/)

Install dependencies with:
```sh
pip install xarray cdsapi cartopy matplotlib seaborn ibicus pandas numpy requests python-dateutil
```

---

## Notes

- **Copernicus API Key:**  
  You must set up your [CDS API key](https://ads.atmosphere.copernicus.eu/how-to-api) for data downloads.
- **Data Paths:**  
  Adjust file paths in scripts as needed for your directory structure.
- **Thresholds & Locations:**  
  Update threshold values and grid cell indices in analysis scripts for your specific crop or region.

---

## References

- [Copernicus Climate Data Store](https://cds.climate.copernicus.eu/)
- [ibicus bias adjustment library](https://github.com/ibicus-org/ibicus)

---

**Contact:**  
For questions or contributions, please open an issue or pull request on this repository.
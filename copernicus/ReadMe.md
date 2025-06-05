# Copernicus Seasonal Forecasting Scripts – BUCRA

This folder contains scripts and data for downloading, converting, bias-adjusting, and analyzing Copernicus seasonal forecast data. The workflow is divided into several steps, but for convenience, pre-converted data is provided so you can start with analysis immediately if desired.

## 1. Data Download Scripts

- `download_hindcast_projections_copernicus.py`
- `download_hindcast_reanalysis_copernicus.py`
- `download_validation_projections_copernicus.py`
- `download_validation_reanalysis_copernicus.py`

These scripts download the necessary Copernicus data for error analysis and bias correction. **Before running these scripts, set up access to the Copernicus API** (see [Copernicus API instructions](https://cds.climate.copernicus.eu/how-to-api)).

## 2. Data Conversion Script

- `convert_data.py`

This script converts downloaded GRIB data to NetCDF format. The original data has both a time and a lead-time dimension (plus longitude and latitude). The script reshapes the data so that the time dimension indicates the month for which the data is valid, and the lead-time dimension indicates how many months in advance the forecast was made. This makes it easier to compare forecasts for the same month from different lead-times. Note: The first six months of data are lost in reshaping, so the first year is removed from all datasets to keep them aligned and starting in January.

## 3. Bias Adjustment and Error Analysis

- `bias_adjustment_analysis.ipynb`

This Jupyter Notebook contains code for bias adjustment and analysis of the converted hindcast and validation data. You can run individual cells to generate specific figures or analyses. All error-analysis figures in the report and annex can be generated using this notebook as a template. The variable analyzed (e.g., mean, min, or max temperature) can be changed at the start of the script. **You do not need to run this script before retrieving the latest forecasts.**

## 4. Retrieving Latest Forecasts

- `retrieve_latest_copernicus.py`

This script downloads, debiases (where possible), and converts the latest six-month forecast available from the Copernicus Data Store via API. If min/max temperature observations are available, you can uncomment relevant code to enable debiasing for all variables. **Note:** If you change the selected area for Copernicus downloads, you must also redownload and convert the hindcast projections and reanalysis data for the same area for debiasing to work correctly. Follow the [Copernicus API instructions](https://cds.climate.copernicus.eu/how-to-api) to set up access.

## 5. Utilities

- `utils.py`

Contains utility functions for loading and reshaping GRIB data, used by several scripts.

## 6. Data

- `GRIB_data/` (folder)
- `NetCDF_data/` (folder)

These folders contain sample downloaded and converted Copernicus data. They are also available on the GitHub branch. This allows you to use `bias_adjustment_analysis.ipynb` and `retrieve_latest_copernicus.py` without first downloading and converting data.

---

**Typical Workflow:**

1. Run the data download scripts to obtain raw Copernicus data (unless you use the provided data).
2. Use `convert_data.py` to convert and reshape the data to NetCDF.
3. Analyze and bias-adjust the data using `bias_adjustment_analysis.ipynb`.
4. Retrieve and process the latest forecasts with `retrieve_latest_copernicus.py`.

## Dependencies

- Python 3.x
- requirements.txt

Install dependencies with:
```sh
pip install -r requirements.txt
```
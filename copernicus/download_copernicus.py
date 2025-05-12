# %%
# Get Copernicus data
import cdsapi

dataset = "seasonal-monthly-single-levels"
request = {
    "originating_centre": "ecmwf",
    "system": "51",
    "variable": ["2m_temperature"],
    "product_type": ["monthly_mean"],
    "year": ["2024", "2025"],
    "month": ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"],
    "leadtime_month": ["1", "2", "3", "4", "5", "6"],
    "data_format": "grib",
    "area": [35, 25, 28, 35],
}

client = cdsapi.Client()
client.retrieve(dataset, request).download()

# %%

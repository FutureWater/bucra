# %%

import cdsapi

dataset = "seasonal-monthly-single-levels"
request = {
    "originating_centre": "eccc",
    "system": "5",
    "variable": ["2m_temperature"],
    "year": [
        "2017",
        "2018",
    ],
    "month": ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"],
    "leadtime_month": ["1", "2", "3", "4", "5", "6"],
    "data_format": "grib",
    "product_type": ["monthly_mean"],
    "area": [31.8, 28.4, 28, 31.4],
}
target = "GRIB_data/ECCC_validation_projections.grib"

client = cdsapi.Client()
client.retrieve(dataset, request, target)

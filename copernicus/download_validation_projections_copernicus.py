# %%

import cdsapi

dataset = "seasonal-monthly-single-levels"
request = {
    "originating_centre": "ecmwf",
    "system": "5",
    "variable": [
        "2m_temperature",
        "minimum_2m_temperature_in_the_last_24_hours",
        "maximum_2m_temperature_in_the_last_24_hours",
    ],
    "year": [
        "2017",
        "2018",
        "2019",
        "2020",
    ],
    "month": ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"],
    "leadtime_month": ["1", "2", "3", "4", "5", "6"],
    "data_format": "grib",
    "product_type": ["monthly_mean"],
    "area": [31.8, 28.4, 28, 31.4],
}
target = "GRIB_data/ECMWF_validation_projections.grib"

client = cdsapi.Client()
client.retrieve(dataset, request, target)

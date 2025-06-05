# %%

import cdsapi

dataset = "seasonal-monthly-single-levels"
request = {
    "originating_centre": "eccc",
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
        "2021",
        "2022",
        "2023",
        "2024",
    ],
    "month": ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"],
    "leadtime_month": ["1", "2", "3", "4", "5", "6"],
    "data_format": "grib",
    "product_type": ["monthly_mean"],
    "area": [31.7, 24.5, 22, 37],
}
target = "GRIB_data/ECCC_validation_projections.grib"

client = cdsapi.Client()
client.retrieve(dataset, request, target)

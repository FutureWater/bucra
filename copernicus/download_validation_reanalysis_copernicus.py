import cdsapi

dataset = "reanalysis-era5-single-levels-monthly-means"
request = {
    "product_type": ["monthly_averaged_reanalysis"],
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
    "time": ["00:00"],
    "data_format": "grib",
    "download_format": "unarchived",
    "variable": [
        "2m_temperature",
        "minimum_2m_temperature_in_the_last_24_hours",
        "maximum_2m_temperature_in_the_last_24_hours",
    ],
    "area": [31.7, 24.5, 22, 37],
}
target = "GRIB_data/reanalysis_validation_copernicus.grib"

client = cdsapi.Client()
client.retrieve(dataset, request, target)

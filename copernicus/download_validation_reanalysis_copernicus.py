import cdsapi

dataset = "reanalysis-era5-single-levels-monthly-means"
request = {
    "product_type": ["monthly_averaged_reanalysis"],
    "year": [
        "2017",
        "2018",
        "2019",
        "2020",
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
    "area": [31.8, 28.4, 28, 31.4],
}
target = "GRIB_data/reanalysis_validation_copernicus.grib"

client = cdsapi.Client()
client.retrieve(dataset, request, target)

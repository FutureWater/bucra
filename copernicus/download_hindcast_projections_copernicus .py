# %%
# # Get Copernicus data
# import cdsapi

# dataset = "seasonal-monthly-single-levels"
# request = {
#     "originating_centre": "ecmwf",
#     "system": "51",
#     "variable": ["2m_temperature"],
#     "product_type": ["monthly_mean"],
#     "year": ["2024", "2025"],
#     "month": ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"],
#     "leadtime_month": ["1", "2", "3", "4", "5", "6"],
#     "data_format": "grib",
#     "area": [35, 25, 28, 35],
# }

# client = cdsapi.Client()
# client.retrieve(dataset, request).download()

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
        "1993",
        "1994",
        "1995",
        "1996",
        "1997",
        "1998",
        "1999",
        "2000",
        "2001",
        "2002",
        "2003",
        "2004",
        "2005",
        "2006",
        "2007",
        "2008",
        "2009",
        "2010",
        "2011",
        "2012",
        "2013",
        "2014",
        "2015",
        "2016",
    ],
    "month": ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"],
    "leadtime_month": ["1", "2", "3", "4", "5", "6"],
    "data_format": "grib",
    "product_type": ["monthly_mean"],
    "area": [31.7, 24.5, 22, 37],
}
target = "GRIB_data/ECCC_hindcast_projections.grib"

client = cdsapi.Client()
client.retrieve(dataset, request, target)

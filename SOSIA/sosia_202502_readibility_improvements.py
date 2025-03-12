"""
Smart Irrigation Advice System using Google Earth Engine and TAHMO data
-----------------------------------------------------------------------

This script calculates irrigation recommendations for farmers based on:
1. Field data from a farmer API (dimensions, crop type, planting dates)
2. Evapotranspiration data from WAPOR and other sources
3. Weather station data from TAHMO (when available)
4. Climate forecast data from CFSv2 (when TAHMO data is unavailable)

The script then posts irrigation advice back to the API for farmer access.
"""

import datetime as datetime
from datetime import date
import ee
import pandas as pd
import requests

# Initialize Google Earth Engine with authentication
ee.Authenticate()
ee.Initialize()

# Fetch farmer field data from the API
url_api = "https://sosia.tahmo.org/api/fields/"
payload_api = {}
headers_api = {
    "Authorization": "Basic bC52ZXJzY2h1cmVuQGZ1dHVyZXdhdGVyLm5sOnA3XlE1OTdNNmx3Wg=="
}
# Get the response and load all farmer data as a dictionary
response_api = requests.request(
    "GET", url_api, headers=headers_api, data=payload_api)
print(response_api.json())

# Process each farmer in the response data
farmer = response_api.json()[0]  # Initialize with first farmer
for farmer in response_api.json():
    # Extract farmer field information from JSON
    field_id = farmer["id"]
    print(field_id)
    field_name = farmer["name"]
    latitude = farmer["latitude"]
    longitude = farmer["longitude"]
    if latitude == None or longitude == None:
        continue
    lat_float = float(latitude)
    long_float = float(longitude)
    field_point = ee.Geometry.Point(long_float, lat_float)
    point_coords = str(longitude) + "," + str(latitude)

    # Skip farmers with missing location data (e.g. Holland Green Tech)

    # Extract irrigation system parameters
    drip_lines_count = farmer["numberOfDriplines"]
    drip_line_length = farmer["lengthOfDriplines"]
    bed_width = farmer["bedWidth"]
    emitter_spacing = farmer["emitterSpacing"]
    emitter_flow_rate = farmer["emitterFlowRate"]
    init_time = farmer["initialisationTime"]

    # Calculate derived irrigation system parameters
    total_flow_rate = (drip_line_length / emitter_spacing) * emitter_flow_rate
    irrigated_area = drip_line_length * (bed_width / drip_lines_count)
    loss_rate = 1.1  # Loss rate of drip irrigation system (10% losses)
    application_efficiency = 0.9  # Application efficiency (90% efficiency)

    # Extract crop information
    crop_type = farmer["cropSpecific"]["cropType"]
    planting_date = ee.Date(farmer["cropSpecific"]["plantingDate"])
    harvest_date = ee.Date(farmer["cropSpecific"]
                           ["lastExpectedHarvestingDate"])

    # Set crop-specific parameters based on crop type
    # Each crop type has specific kc values and growth stage durations
    if crop_type == "Habanero Peppers":
        # Initial, mid-season, and late season crop coefficients
        kc_initial = 0.6
        kc_mid = 1.05
        kc_late = 1.05

        # Define growth stages (days from planting)
        season_dev_start = planting_date.advance(0, "day")
        season_dev_end = planting_date.advance(65, "day")
        season_mid_start = planting_date.advance(66, "day")
        season_mid_end = season_dev_end.advance(40, "day")
        season_end_start = season_dev_end.advance(41, "day")
        season_end = harvest_date.advance(0, "day")

    elif crop_type == "French Beans":
        kc_initial = 0.5
        kc_mid = 1.05
        kc_late = 1.05

        season_dev_start = planting_date.advance(0, "day")
        season_dev_end = planting_date.advance(65, "day")
        season_mid_start = planting_date.advance(66, "day")
        season_mid_end = season_dev_end.advance(40, "day")
        season_end_start = season_dev_end.advance(41, "day")
        season_end = harvest_date.advance(0, "day")

    elif crop_type == "Lettuce":
        kc_initial = 0.5
        kc_mid = 1.05
        kc_late = 1.05

        season_dev_start = planting_date.advance(0, "day")
        season_dev_end = planting_date.advance(85, "day")
        season_mid_start = planting_date.advance(86, "day")
        season_mid_end = season_dev_end.advance(40, "day")
        season_end_start = season_dev_end.advance(41, "day")
        season_end = harvest_date.advance(0, "day")

    elif crop_type == "Brassica":
        kc_initial = 0.7
        kc_mid = 1.05
        kc_late = 1.05

        season_dev_start = planting_date.advance(0, "day")
        season_dev_end = planting_date.advance(65, "day")
        season_mid_start = planting_date.advance(66, "day")
        season_mid_end = season_dev_end.advance(40, "day")
        season_end_start = season_dev_end.advance(41, "day")
        season_end = harvest_date.advance(0, "day")

    elif crop_type == "Cucumbers":
        kc_initial = 0.60
        kc_mid = 1.00
        kc_late = 1.00

        season_dev_start = planting_date.advance(0, "day")
        season_dev_end = planting_date.advance(60, "day")
        season_mid_start = planting_date.advance(61, "day")
        season_mid_end = season_dev_end.advance(50, "day")
        season_end_start = season_dev_end.advance(51, "day")
        season_end = harvest_date.advance(0, "day")

    elif crop_type == "Okra":
        kc_initial = 0.30
        kc_mid = 1.00
        kc_late = 0.90

        season_dev_start = planting_date.advance(0, "day")
        season_dev_end = planting_date.advance(41, "day")
        season_mid_start = planting_date.advance(42, "day")
        season_mid_end = season_dev_end.advance(25, "day")
        season_end_start = season_dev_end.advance(26, "day")
        season_end = harvest_date.advance(0, "day")

    # Time variables for analysis
    months = ee.List.sequence(1, 12)
    days_of_year = ee.List.sequence(1, 365)
    list_days = ee.List.sequence(1, 9)
    historical_start_date = ee.Date("2010-01-01")
    historical_end_date = ee.Date("2022-12-31")
    current_day = date.today()
    current_day_str = str(current_day)
    today_ee = ee.Date(str(current_day))
    today_formatted = today_ee.format("YYYY-MM-dd")
    past_date = today_ee.advance(-10, "day")
    past_date_formatted = past_date.format("YYYY-MM-dd")
    beginning_date = ee.Date("2021-12-31")
    millis_in_day = 24 * 60 * 60 * 1000  # Milliseconds in a day

    # Load required Earth Engine datasets
    # Reference Evapotranspiration
    wapor_ret = ee.ImageCollection("FAO/WAPOR/2/L1_RET_E")
    # Global Precipitation
    gpm = ee.ImageCollection("NASA/GPM_L3/IMERG_V06")
    # Climate Forecast System
    cfsv2 = ee.ImageCollection("NOAA/CFSV2/FOR6H")
    # Digital Elevation Model
    dem = ee.Image("NASA/NASADEM_HGT/001")

    # Create buffer around field point for spatial analysis (30m radius)
    clip_geometry = field_point.buffer(30)

    # Function to clip images to the field buffer area
    def clip_image(image):
        """Clip an image to the field buffer area."""
        clipped = image.clip(clip_geometry)
        return clipped

    # Filter WAPOR data for historical date range and clip to field area
    wapor_filtered = wapor_ret.filterDate(
        historical_start_date, historical_end_date).map(clip_image)

    # Function to correct WAPOR values (divide by 10) and add as a new band
    def correct_wapor_values(img):
        """Correct WAPOR evapotranspiration values by dividing by 10."""
        corrected = img.select("L1_RET_E").divide(10).rename("corrected")
        return img.addBands(corrected)

    # Apply correction to all WAPOR images
    wapor_corrected = wapor_filtered.map(correct_wapor_values)

    # Function to calculate daily average WAPOR values for a specific day of year
    def calculate_wapor_daily_avg(day_of_year):
        """
        Calculate the mean WAPOR value for a specific day of year across multiple years.
        Returns a single value representing the average for that day of year.
        """
        # Get all images for this day of year across all years
        day_images = wapor_corrected.select("corrected").filter(
            ee.Filter.calendarRange(start=day_of_year, field="day_of_year")
        )
        # Calculate mean and format to 2 decimal places
        mean_wapor = (
            ee.Image(day_images.mean())
            .multiply(100)
            .round()
            .divide(100)
            .set("DOY", day_of_year)
        )
        return ee.Number(mean_wapor)

    # Create a collection of daily average WAPOR values for each day of the year
    wapor_by_day = ee.ImageCollection(
        days_of_year.map(calculate_wapor_daily_avg))

    print(crop_type)  # Print crop type for debugging

    # Historical Crop Schedule - Calculate irrigation needs for entire growing season
    # -----------------------------------------------------------------------------

    # Calculate reference ET for the entire growing season
    season_start_millis = season_dev_start.millis()
    season_end_millis = season_end.millis()

    # Function to convert date in milliseconds to day of year
    def milliseconds_to_day_of_year(date_millis):
        """Convert date in milliseconds to day of year (1-365)."""
        return ee.Number.parse(ee.Date(date_millis).format("DDD"))

    # Get all days in the growing season as days of year
    days_in_season = ee.List.sequence(
        season_start_millis, season_end_millis, millis_in_day).map(milliseconds_to_day_of_year)

    # Calculate historical average reference ET for each day in the growing season
    wapor_ref_seasonal = ee.ImageCollection(
        days_in_season.map(calculate_wapor_daily_avg))

    # Function to process ET reference values (round to 2 decimals)
    def process_et_ref(image):
        """Process ET reference values by rounding to 2 decimal places."""
        return (
            image.multiply(100)
            .round()
            .divide(100)
            .copyProperties(image)
            .set("DOY", image.get("DOY"))
        )

    et_ref_processed = wapor_ref_seasonal.map(process_et_ref)

    # Function to add the Date property to each image (days from planting date)
    def add_date_from_planting(img):
        """Add Date property to image based on days from planting date."""
        day_index = ee.Number.parse(img.get("system:index"))
        date_val = planting_date.advance(day_index, "day")
        return img.set("Date", date_val)

    et_ref_with_date = et_ref_processed.map(add_date_from_planting)

    # Select the corrected band and rename it
    et_ref_renamed = et_ref_with_date.select(
        ["corrected"], ["1 ETref in mm/day"])

    # Calculate ET for each growth stage (development, mid-season, end)
    # ---------------------------------------------------------------------

    # DEVELOPMENT STAGE ET CALCULATION
    dev_start_millis = season_dev_start.millis()
    dev_end_millis = season_dev_end.millis()

    # Get days of year for development stage
    days_in_dev_stage = ee.List.sequence(
        dev_start_millis, dev_end_millis, millis_in_day).map(milliseconds_to_day_of_year)

    # Get WAPOR values for development stage days
    wapor_dev_stage = ee.ImageCollection(
        days_in_dev_stage.map(calculate_wapor_daily_avg))

    # Apply kc_initial to the development stage ET values
    def apply_kc_initial(image):
        """Apply the initial stage crop coefficient (kc_initial) to ET values."""
        return (
            image.multiply(kc_initial)
            .multiply(100)
            .round()
            .divide(100)
            .copyProperties(image)
            .set("DOY", image.get("DOY"))
        )

    et_dev_stage = wapor_dev_stage.map(apply_kc_initial)

    # Add Date property to development stage images
    def add_dev_stage_date(img):
        """Add Date property for development stage images."""
        day_index = ee.Number.parse(img.get("system:index"))
        date_val = planting_date.advance(day_index, "day")
        return img.set("Date", date_val)

    et_dev_with_date = et_dev_stage.map(add_dev_stage_date)

    # MID-SEASON STAGE ET CALCULATION
    mid_start_millis = season_mid_start.millis()
    mid_end_millis = season_mid_end.millis()

    # Get days of year for mid-season stage
    days_in_mid_stage = ee.List.sequence(
        mid_start_millis, mid_end_millis, millis_in_day).map(milliseconds_to_day_of_year)

    # Get WAPOR values for mid-season stage days
    wapor_mid_stage = ee.ImageCollection(
        days_in_mid_stage.map(calculate_wapor_daily_avg))

    # Apply kc_mid to the mid-season stage ET values
    def apply_kc_mid(image):
        """Apply the mid-season crop coefficient (kc_mid) to ET values."""
        return (
            image.multiply(kc_mid)
            .multiply(100)
            .round()
            .divide(100)
            .copyProperties(image)
            .set("DOY", image.get("DOY"))
        )

    et_mid_stage = wapor_mid_stage.map(apply_kc_mid)

    # Add Date property to mid-season stage images
    def add_mid_stage_date(img):
        """Add Date property for mid-season stage images."""
        day_index = ee.Number.parse(img.get("system:index"))
        date_val = season_mid_start.advance(day_index, "day")
        return img.set("Date", date_val)

    et_mid_with_date = et_mid_stage.map(add_mid_stage_date)

    # END STAGE ET CALCULATION
    end_start_millis = season_end_start.millis()
    end_end_millis = season_end.millis()

    # Get days of year for end stage
    days_in_end_stage = ee.List.sequence(
        end_start_millis, end_end_millis, millis_in_day).map(milliseconds_to_day_of_year)

    # Get WAPOR values for end stage days
    wapor_end_stage = ee.ImageCollection(
        days_in_end_stage.map(calculate_wapor_daily_avg))

    # Apply kc_late to the end stage ET values
    def apply_kc_late(image):
        """Apply the late season crop coefficient (kc_late) to ET values."""
        return (
            image.multiply(kc_late)
            .multiply(100)
            .round()
            .divide(100)
            .copyProperties(image)
            .set("DOY", image.get("DOY"))
        )

    et_end_stage = wapor_end_stage.map(apply_kc_late)

    # Add Date property to end stage images
    def add_end_stage_date(img):
        """Add Date property for end stage images."""
        day_index = ee.Number.parse(img.get("system:index"))
        date_val = season_end_start.advance(day_index, "day")
        return img.set("Date", date_val)

    et_end_with_date = et_end_stage.map(add_end_stage_date)

    # Merge ET from all three stages
    et_merged = et_dev_with_date.merge(et_mid_with_date)
    et_all_stages = et_merged.merge(et_end_with_date)  # mm per day
    et_all_renamed = et_all_stages.select(["corrected"], ["2 ETc in mm/day"])

    # Calculate irrigation needs based on ET values
    # ---------------------------------------------

    def calculate_irrigation_needs(image):
        """
        Calculate irrigation volume (m³/day) and time (minutes/day) from ET data.

        Args:
            image: Image containing ET data in mm/day

        Returns:
            Image with added bands for irrigation volume and time
        """
        # Calculate irrigation volume in m³/day
        irr_volume = (
            image.divide(1000)                     # Convert mm to m
            # Multiply by area to get volume
            .multiply(irrigated_area)
            .multiply(loss_rate)                   # Apply loss rate
            # Account for application efficiency
            .divide(application_efficiency)
            .multiply(100).round().divide(100)     # Round to 2 decimals
            .copyProperties(image)
            .set("DOY", image.get("DOY"))
        )

        # Calculate irrigation time in minutes/day
        irr_time = (
            image.divide(1000)                     # Convert mm to m
            # Multiply by area to get volume
            .multiply(irrigated_area)
            .multiply(loss_rate)                   # Apply loss rate
            # Account for application efficiency
            .divide(application_efficiency)
            .multiply(60000)                       # Convert hours to minutes
            .divide(total_flow_rate)               # Divide by flow rate
            .round()
            .add(init_time)                        # Add initialization time
            .copyProperties(image)
            .set("DOY", image.get("DOY"))
        )

        # Create new bands for volume and time
        volume_band = ee.Image([irr_volume])
        time_band = ee.Image([irr_time])

        # Name the bands appropriately
        volume_band_renamed = volume_band.select(
            ["2 ETc in mm/day"], ["3 Irrigation needs in m3/day"]
        )
        time_band_renamed = time_band.select(
            ["2 ETc in mm/day"], ["4 Irrigation time in min/day"]
        )

        # Add the new bands to the original image
        return image.addBands(volume_band_renamed).addBands(time_band_renamed)

    # Apply the irrigation calculation to the ET data
    total_irrigation = et_all_renamed.map(calculate_irrigation_needs)

    # Set up filters and joins for combining datasets
    date_filter = ee.Filter.equals(leftField="Date", rightField="Date")
    doy_filter = ee.Filter.equals(leftField="DOY", rightField="DOY")
    simple_join = ee.Join.inner()

    # Join reference ET with crop ET and irrigation needs
    inner_join = ee.ImageCollection(simple_join.apply(
        et_ref_renamed, total_irrigation, date_filter))

    # Function to merge joined images
    def merge_joined_images(feature):
        """Merge primary and secondary images from a join result."""
        return ee.Image.cat(feature.get("primary"), feature.get("secondary"))

    # Merge joined images (125 days of the season)
    historical_schedule = inner_join.map(merge_joined_images)

    # Convert EE images to tabular data for API
    # -----------------------------------------

    def extract_field_means(img):
        """
        Extract mean values for all required bands at the field location.

        Args:
            img: Image with ET and irrigation bands

        Returns:
            Image with properties set for all required values
        """
        # Get day of year
        doy = img.get("DOY")
        date_string = ee.Number(doy).format()
        model_run_id = ee.String(date_string).cat("_SAT").cat(current_day_str)

        # Calculate mean values for each band at the field location
        et_ref_mean = img.reduceRegion(
            reducer=ee.Reducer.mean(), geometry=clip_geometry, scale=30
        ).get("1 ETref in mm/day")

        et_crop_mean = img.reduceRegion(
            reducer=ee.Reducer.mean(), geometry=clip_geometry, scale=30
        ).get("2 ETc in mm/day")

        irr_vol_mean = img.reduceRegion(
            reducer=ee.Reducer.mean(), geometry=clip_geometry, scale=30
        ).get("3 Irrigation needs in m3/day")

        irr_time_mean = img.reduceRegion(
            reducer=ee.Reducer.mean(), geometry=clip_geometry, scale=30
        ).get("4 Irrigation time in min/day")

        # Set all calculated values as properties
        return (
            img.set("Etref", et_ref_mean)
            .set("Etc", et_crop_mean)
            .set("Irrvol", irr_vol_mean)
            .set("Irrtime", irr_time_mean)
            .set("Date", doy)
            .set("ModelRun", model_run_id)
        )

    # Apply extraction to all images
    field_data_images = historical_schedule.map(extract_field_means)

    # Convert to a list format for API
    field_data_list = (
        field_data_images.reduceColumns(
            ee.Reducer.toList(5), ["Date", "ModelRun",
                                   "Etc", "Irrvol", "Irrtime"]
        )
        .values()
        .get(0)
    )

    # Convert to pandas DataFrame
    historical_df = pd.DataFrame(
        field_data_list.getInfo(),
        columns=[
            "date",
            "modelRun",
            "evaporation",
            "advisedWaterVolume",
            "advisedIrrigationTime",
        ],
    )

    # Format dates and add required columns
    current_year = 2024  # Hardcoded for the example
    historical_df["date"] = pd.to_datetime(
        current_year * 1000 + historical_df["date"], format="%Y%j")
    historical_df["modelRun"] = historical_df["date"].dt.strftime(
        "%Y%m%d") + "_SAT" + current_day_str
    historical_df["dateDisplay"] = historical_df["date"].dt.strftime(
        "%d-%m").map(lambda x: str(x)[-5:])
    historical_df["date"] = pd.to_datetime(
        historical_df["date"], format="%d-%m-%Y").dt.strftime("%Y-%m-%d")
    historical_df["field"] = field_id

    # Select and order columns
    historical_df = historical_df[
        [
            "date",
            "dateDisplay",
            "modelRun",
            "evaporation",
            "advisedWaterVolume",
            "advisedIrrigationTime",
            "field",
        ]
    ]

    # Convert to JSON for API
    historical_json = historical_df.to_json(orient="records")
    print(historical_json)

    # Post historical irrigation schedule to API
    url_post = "https://sosia.tahmo.org/api/seasonal_schedule/"
    headers_post = {
        "Content-type": "application/json",
        "Authorization": "Basic bC52ZXJzY2h1cmVuQGZ1dHVyZXdhdGVyLm5sOnA3XlE1OTdNNmx3Wg==",
    }

    response_post = requests.post(
        url_post, headers=headers_post, data=historical_json)
    print(response_post.json())

    # HINDCAST CALCULATION - Recent and current irrigation needs
    # --------------------------------------------------------
    # Use recent data to provide current irrigation recommendations

    # Get TAHMO weather station data if available
    url_tahmo = "https://smartirrigation.tahmo.org/sosia/nearestwetnessreport/" + point_coords
    payload_tahmo = {}
    headers_tahmo = {
        "Authorization": "Basic ZnV0dXJld2F0ZXI6R2MzYVdMN3kyckRkR2Y3RQ=="}

    response_tahmo = requests.request(
        "GET", url_tahmo, headers=headers_tahmo, data=payload_tahmo
    )

    tahmo_result = response_tahmo.json()

    # Process TAHMO data if available, otherwise use CFSv2 data
    if tahmo_result:
        # TAHMO data is available - process it
        print("Processing TAHMO weather station data")

        # Normalize and prepare TAHMO results
        tahmo_data = pd.json_normalize(tahmo_result, "results", "distance")
        tahmo_cleaned = tahmo_data.drop(
            columns=["_id", "Wetness", "FC", "RAM", "WP", "Eact"])
        tahmo_cleaned["Day"] = pd.to_datetime(
            tahmo_cleaned["Time"]).dt.strftime("%m/%d/%y")
        tahmo_cleaned.drop(columns=["Time"])

        # Pivot data by day and station
        tahmo_days = tahmo_cleaned.pivot_table("Eref", "Day", "Station")
        tahmo_distance = tahmo_cleaned.pivot_table("distance", "Station")

        # Apply Inverse Distance Weighting based on number of stations
        if tahmo_distance.count()[0] == 1:
            print("One station retrieved, using nearest station value")
            tahmo_point_value = tahmo_days

        elif tahmo_distance.count()[0] == 2:
            print("Two stations retrieved, applying IDW")
            d1 = tahmo_distance["distance"][0] + 28
            d2 = tahmo_distance["distance"][1]
            v1 = tahmo_days.iloc[:, 0]
            v2 = tahmo_days.iloc[:, 1]
            tahmo_point_value = ((v1 / d1) + (v2 / d2)) / ((1 / d1) + (1 / d2))
            print(tahmo_point_value)

        else:
            print("Three stations retrieved, applying IDW")
            d1 = tahmo_distance["distance"][0]
            d2 = tahmo_distance["distance"][1]
            d3 = tahmo_distance["distance"][2]
            v1 = tahmo_days.iloc[:, 0]
            v2 = tahmo_days.iloc[:, 1]
            v3 = tahmo_days.iloc[:, 2]
            tahmo_point_value = ((v1 / d1) + (v2 / d2) + (v3 / d3)) / (
                (1 / d1) + (1 / d2) + (1 / d3)
            )
            print(tahmo_point_value)

        # Get the TAHMO values for processing
        tahmo_values = tahmo_days.iloc[:, 0]
        tahmo_values_list = tahmo_values.reset_index(drop=True).tolist()

        # Create EE data structure with TAHMO values
        week_ago = today_ee.advance(-6, "day")
        today_millis = today_ee.millis()
        week_ago_millis = week_ago.millis()

        # Get recent days as day of year values
        hindcast_days = ee.List.sequence(
            week_ago_millis, today_millis, millis_in_day).map(milliseconds_to_day_of_year)

        # Create images with the structure of WAPOR but values from TAHMO
        tahmo_hindcast = ee.ImageCollection(
            hindcast_days.map(calculate_wapor_daily_avg))

        # Zip WAPOR structure with TAHMO values
        max_elements = 1000
        zipped_list = tahmo_hindcast.toList(
            max_elements).zip(tahmo_values_list)

        # Function to create images with TAHMO values
        def create_tahmo_image(list_item):
            """Create an image with TAHMO ET values using WAPOR image structure."""
            list_item = ee.List(list_item)
            wapor_img = ee.Image(list_item.get(0))
            tahmo_value = list_item.getNumber(1)

            # Create a new image with the TAHMO value
            tahmo_image = wapor_img.multiply(0).add(tahmo_value)
            tahmo_image = tahmo_image.double().rename("1 ETref in mm/day")

            # Copy all properties from the original image
            return tahmo_image.copyProperties(wapor_img, wapor_img.propertyNames())

        # Create TAHMO ET reference collection
        tahmo_images = zipped_list.map(create_tahmo_image)
        tahmo_et_collection = ee.ImageCollection.fromImages(tahmo_images)

    else:
        # No TAHMO data available - use CFSv2 forecast data
        print("No TAHMO data available, using CFSv2 forecast data")

        # Define date range for hindcast (from 9 to 2 days ago)
        week_ago = today_ee.advance(-9, "day")
        week_ago_formatted = ee.Date(week_ago.format("YYYY-MM-dd"))
        week_ago_doy = ee.Number.parse(week_ago.format("DDD"))

        yesterday = today_ee.advance(-2, "day")
        yesterday_formatted = ee.Date(yesterday.format("YYYY-MM-dd"))
        yesterday_doy = ee.Number.parse(yesterday.format("DDD"))

        # Get days in the hindcast period
        hindcast_days = ee.List.sequence(week_ago_doy, yesterday_doy, 1)

        # Filter CFSv2 data for the hindcast period
        cfsv2_filtered_initial = cfsv2.filterDate(
            week_ago_formatted, yesterday_formatted)
        cfsv2_filtered = ee.ImageCollection(
            cfsv2_filtered_initial.map(clip_image))

        # Function to calculate net radiation and wind from CFSv2 data
        def calculate_radiation_and_wind(image):
            """Calculate net radiation and wind components from CFSv2 data."""
            # Get radiation components
            downward_shortwave = image.select(
                "Downward_Short-Wave_Radiation_Flux_surface_6_Hour_Average"
            )
            downward_longwave = image.select(
                "Downward_Long-Wave_Radp_Flux_surface_6_Hour_Average")
            upward_longwave = image.select(
                "Upward_Long-Wave_Radp_Flux_surface_6_Hour_Average")

            # Calculate upward shortwave (albedo * downward)
            upward_shortwave = downward_shortwave.multiply(0.23).rename(
                "Upward_Short-Wave_Radiation_Flux_surface_6_Hour_Average"
            )

            # Calculate net radiation
            net_radiation = (
                downward_shortwave.subtract(upward_shortwave)
                .add(downward_longwave.subtract(upward_longwave))
                .rename("Net_Radiation_6_Hour_Average")
            )

            # Calculate wind speed from components
            u_wind = image.select("u-component_of_wind_height_above_ground")
            v_wind = image.select("v-component_of_wind_height_above_ground")
            u_wind_squared = u_wind.pow(2)
            v_wind_squared = v_wind.pow(2)
            wind_speed = (u_wind_squared.add(v_wind_squared)
                          ).sqrt().rename("Wind_component")

            # Add new bands to the image
            return image.addBands(net_radiation).addBands(wind_speed)

        # Calculate radiation and wind for all images
        cfsv2_with_radiation = cfsv2_filtered.map(calculate_radiation_and_wind)
        cfsv2_temperature = cfsv2_filtered.select(
            "Temperature_height_above_ground")

        # Calculate daily means over the hindcast period
        days_count = yesterday_formatted.difference(week_ago_formatted, "days")

        # Function to calculate daily means
        def calculate_daily_mean(day_offset):
            """Calculate daily mean values for a specific day in the hindcast period."""
            start_date = week_ago_formatted.advance(day_offset, "days")
            end_date = start_date.advance(1, "days")
            date_formatted = start_date.format("YYYY-MM-dd")
            doy = ee.Number.parse(start_date.format("DDD"))

            return (
                cfsv2_with_radiation.filterDate(start_date, end_date)
                .mean()
                .set("Date", date_formatted)
                .set("DOY", doy)
            )

        # Calculate daily means for each day in the hindcast period
        daily_means = ee.ImageCollection(
            ee.List.sequence(0, days_count.subtract(1)
                             ).map(calculate_daily_mean)
        )

        # Function to calculate daily maximum temperature
        def calculate_daily_max_temp(day_offset):
            """Calculate daily maximum temperature for a specific day."""
            start_date = week_ago_formatted.advance(day_offset, "days")
            end_date = start_date.advance(1, "days")
            return (
                cfsv2_temperature.filterDate(start_date, end_date)
                .max()
                .set("Date", start_date.format("YYYY-MM-dd"))
            )

        # Calculate daily maximum temperatures
        daily_max_temps = ee.ImageCollection(
            ee.List.sequence(0, days_count.subtract(1)).map(
                calculate_daily_max_temp)
        )

        temp_max_renamed = daily_max_temps.map(
            lambda image: image.select("Temperature_height_above_ground").rename(
                "Max Temperature"
            )
        )

        # Function to calculate daily minimum temperature
        def calculate_daily_min_temp(day_offset):
            """Calculate daily minimum temperature for a specific day."""
            start_date = week_ago_formatted.advance(day_offset, "days")
            end_date = start_date.advance(1, "days")
            return (
                cfsv2_temperature.filterDate(start_date, end_date)
                .min()
                .set("Date", start_date.format("YYYY-MM-dd"))
            )

        # Calculate daily minimum temperatures
        daily_min_temps = ee.ImageCollection(
            ee.List.sequence(0, days_count.subtract(1)).map(
                calculate_daily_min_temp)
        )

        temp_min_renamed = daily_min_temps.map(
            lambda image: image.select("Temperature_height_above_ground").rename(
                "Min Temperature"
            )
        )

        # Join daily min and max temperatures
        min_max_join = ee.ImageCollection(
            simple_join.apply(temp_max_renamed, temp_min_renamed, date_filter)
        )

        temp_combined = min_max_join.map(
            lambda feature: ee.Image.cat(
                feature.get("primary"), feature.get("secondary")
            )
        )

        # Join temperatures with daily means
        all_data_join = ee.ImageCollection(
            simple_join.apply(daily_means, temp_combined, date_filter)
        )

        cfsv2_combined = all_data_join.map(
            lambda feature: ee.Image.cat(
                feature.get("primary"), feature.get("secondary")
            )
        )

        # Function to calculate ET0 using Penman-Monteith equation
        def calculate_et0_penman_monteith(image):
            """
            Calculate reference evapotranspiration (ET0) using the Penman-Monteith equation.

            Args:
                image: Image with temperature, humidity, pressure, wind, and radiation data

            Returns:
                Image with added ET0 band
            """
            # Extract required inputs
            min_temp_kelvin = image.select("Min Temperature")
            max_temp_kelvin = image.select("Max Temperature")
            wind_speed_10m = image.select("Wind_component")
            specific_humidity = image.select(
                "Specific_humidity_height_above_ground")
            pressure_pascals = image.select("Pressure_surface")
            radiation_net_watts = image.select("Net_Radiation_6_Hour_Average")

            # Convert units
            min_temp_celsius = min_temp_kelvin.subtract(273.15).rename("Tmin")
            max_temp_celsius = max_temp_kelvin.subtract(273.15).rename("Tmax")
            pressure_kpa = pressure_pascals.divide(1000).rename("Atm pressure")
            wind_speed_2m = wind_speed_10m.multiply(
                0.75).rename("Wind")  # Convert 10m to 2m height
            radiation_net_mjm2 = radiation_net_watts.multiply(
                0.0864).rename("Rnet")  # W/m² to MJ/m²/day

            # Base for exponential calculations
            exp_base = pressure_kpa.multiply(0).add(2.71828)  # e constant

            # Calculate mean temperature
            temp_mean = min_temp_celsius.add(max_temp_celsius).divide(2)

            # Calculate slope of saturation vapor pressure curve (delta)
            delta_term1 = temp_mean.multiply(
                17.27).divide(temp_mean.add(237.3))
            delta_term2 = exp_base.pow(delta_term1).multiply(0.6108)
            delta_term3 = temp_mean.add(237.3)
            delta_term4 = delta_term3.pow(2)
            delta = delta_term2.multiply(4098).divide(delta_term4)

            # Calculate saturation vapor pressure (es)
            e0_min_term1 = min_temp_celsius.multiply(
                17.27).divide(min_temp_celsius.add(237.3))
            e0_min = exp_base.pow(e0_min_term1).multiply(0.6108)
            e0_max_term1 = max_temp_celsius.multiply(
                17.27).divide(max_temp_celsius.add(237.3))
            e0_max = exp_base.pow(e0_max_term1).multiply(0.6108)
            es_term1 = e0_max.add(e0_min)
            es = es_term1.divide(2)

            # Calculate actual vapor pressure (ea)
            rh_max = (
                specific_humidity.multiply(pressure_kpa)
                .multiply(1.6077717)
                .divide(e0_min)
                .rename("RHmax")
            )
            rh_min = (
                specific_humidity.multiply(pressure_kpa)
                .multiply(1.6077717)
                .divide(e0_max)
                .rename("RHmin")
            )
            ea_term1 = e0_min.multiply(rh_max)
            ea_term2 = e0_max.multiply(rh_min)
            ea_term3 = ea_term1.add(ea_term2)
            ea = ea_term3.divide(2)

            # Calculate ET0 using Penman-Monteith equation
            part_1 = delta.multiply(radiation_net_mjm2).multiply(0.408)
            psychrometric_constant = pressure_kpa.multiply(
                0.001).divide(1.53634)
            part_2b = temp_mean.add(273)
            part_2 = psychrometric_constant.multiply(900).divide(part_2b)
            part_3 = wind_speed_2m.multiply(es.subtract(ea))
            part_4a = wind_speed_2m.multiply(0.34).add(1)
            part_4 = psychrometric_constant.multiply(part_4a).add(delta)
            et0_part1 = part_1.divide(part_4)
            et0_part2 = part_2.multiply(part_3).divide(part_4)
            et0 = et0_part1.add(et0_part2).rename("1 ETref in mm/day")

            # Add ET0 band to the original image
            return image.addBands(et0)

        # Calculate ET0 for all days in the hindcast period
        cfsv2_et = cfsv2_combined.map(calculate_et0_penman_monteith)
        cfsv2_et_reference = cfsv2_et.select("1 ETref in mm/day")

    # Create Kc values for each growth stage to multiply with ET reference
    # These are used for both TAHMO and CFSv2 methods

    # Development stage Kc values
    def create_kc_dev(image):
        """Create an image with the initial stage Kc value."""
        return (
            image.multiply(0)
            .add(kc_initial)
            .copyProperties(image)
            .set("DOY", image.get("DOY"))
        )

    kc_dev = wapor_dev_stage.map(create_kc_dev)

    def add_date_to_kc_dev(image):
        """Add Date property to development stage Kc image."""
        doy = ee.Number.parse(image.get("system:index"))
        date_formatted = (planting_date.advance(
            doy, "day")).format("YYYY-MM-dd")
        return image.set("Date", date_formatted)

    kc_dev_with_date = kc_dev.map(add_date_to_kc_dev)

    # Mid-season stage Kc values
    def create_kc_mid(image):
        """Create an image with the mid-season Kc value."""
        return (
            image.multiply(0)
            .add(kc_mid)
            .copyProperties(image)
            .set("DOY", image.get("DOY"))
        )

    kc_mid_stage = wapor_mid_stage.map(create_kc_mid)

    def add_date_to_kc_mid(image):
        """Add Date property to mid-season stage Kc image."""
        doy = ee.Number.parse(image.get("system:index"))
        date_formatted = (season_mid_start.advance(
            doy, "day")).format("YYYY-MM-dd")
        return image.set("Date", date_formatted)

    kc_mid_with_date = kc_mid_stage.map(add_date_to_kc_mid)

    # End stage Kc values
    def create_kc_end(image):
        """Create an image with the late season Kc value."""
        return (
            image.multiply(0)
            .add(kc_late)
            .copyProperties(image)
            .set("DOY", image.get("DOY"))
        )

    kc_end_stage = wapor_end_stage.map(create_kc_end)

    def add_date_to_kc_end(image):
        """Add Date property to end stage Kc image."""
        doy = ee.Number.parse(image.get("system:index"))
        date_formatted = (season_end_start.advance(
            doy, "day")).format("YYYY-MM-dd")
        return image.set("Date", date_formatted)

    kc_end_with_date = kc_end_stage.map(add_date_to_kc_end)

    # Merge all Kc values for the entire growing season
    kc_merged = kc_dev_with_date.merge(kc_mid_with_date)
    kc_all_stages = kc_merged.merge(kc_end_with_date)
    kc_all_renamed = kc_all_stages.select(["corrected"], ["2 Daily Kc"])

    # Join ET reference with Kc values - different join methods for TAHMO vs CFSv2
    if tahmo_result:
        # Join TAHMO ET with Kc by DOY (day of year)
        kc_et_join = ee.ImageCollection(simpleJoin.apply(
            tahmo_et_collection, kc_all_renamed, doy_filter))
        print("TAHMO data used")
    else:
        # Join CFSv2 ET with Kc by Date
        kc_et_join = ee.ImageCollection(simple_join.apply(
            cfsv2_et_reference, kc_all_renamed, date_filter))
        print("CFSv2 data used")

    # Merge joined images
    def merge_kc_et_images(feature):
        """Merge ET reference and Kc images."""
        return ee.Image.cat(feature.get("primary"), feature.get("secondary"))

    # 5-7 days of hindcast data
    joined_kc_et = kc_et_join.map(merge_kc_et_images)

    # Calculate crop ET (ETc) by multiplying ET reference with Kc
    def calculate_crop_et(image):
        """Calculate crop evapotranspiration by multiplying ET reference with Kc."""
        et_ref = image.select("1 ETref in mm/day")
        kc = image.select("2 Daily Kc")
        et_crop = et_ref.multiply(kc).rename("3 ETc in mm/day")
        return image.addBands(et_crop)

    et_crop_hindcast = joined_kc_et.map(calculate_crop_et)

    # Calculate irrigation needs for the hindcast period
    def calculate_hindcast_irrigation(image):
        """
        Calculate irrigation volume and time for the hindcast period.
        Similar to calculate_irrigation_needs but adapted for hindcast data structure.
        """
        # Select the ETc band
        et_crop = image.select("3 ETc in mm/day")

        # Calculate irrigation volume (m³/day)
        irr_volume = (
            et_crop.divide(1000)                      # Convert mm to m
            .multiply(irrigated_area)                 # Multiply by area
            .multiply(loss_rate)                      # Apply loss rate
            .divide(application_efficiency)           # Account for efficiency
            .multiply(100).round().divide(100)        # Round to 2 decimals
            .copyProperties(image)
        )

        # Calculate irrigation time (minutes/day)
        irr_time = (
            et_crop.divide(1000)                      # Convert mm to m
            .multiply(irrigated_area)                 # Multiply by area
            .multiply(loss_rate)                      # Apply loss rate
            .divide(application_efficiency)           # Account for efficiency
            .multiply(60000)                          # Convert to minutes
            .divide(total_flow_rate)                  # Divide by flow rate
            .round()
            .add(init_time)                           # Add initialization time
            .copyProperties(image)
        )

        # Create new bands for volume and time
        new_bands = ee.Image([irr_volume, irr_time])
        adjusted_bands = new_bands.select(
            ["3 ETc in mm/day", "3 ETc in mm/day_1"],
            ["4 Irrigation needs in m3/day", "5 Irrigation time in min/day"]
        )

        return image.addBands(adjusted_bands)

    # Apply irrigation calculation to hindcast data
    hindcast_irrigation = et_crop_hindcast.map(calculate_hindcast_irrigation)

    # Extract field means for hindcast data
    def extract_hindcast_means(img):
        """Extract mean values at field location for hindcast data."""
        doy = img.get("DOY")
        date_string = ee.Number(doy).format()
        model_run_id = ee.String(date_string).cat("_SAT").cat(current_day_str)

        # Calculate mean values for each band at the field location
        et_ref_mean = img.reduceRegion(
            reducer=ee.Reducer.mean(), geometry=clip_geometry, scale=30
        ).get("1 ETref in mm/day")

        et_crop_mean = img.reduceRegion(
            reducer=ee.Reducer.mean(), geometry=clip_geometry, scale=30
        ).get("3 ETc in mm/day")

        irr_vol_mean = img.reduceRegion(
            reducer=ee.Reducer.mean(), geometry=clip_geometry, scale=30
        ).get("4 Irrigation needs in m3/day")

        irr_time_mean = img.reduceRegion(
            reducer=ee.Reducer.mean(), geometry=clip_geometry, scale=30
        ).get("5 Irrigation time in min/day")

        # Set all calculated values as properties
        return (
            img.set("Etref", et_ref_mean)
            .set("Etc", et_crop_mean)
            .set("Irrvol", irr_vol_mean)
            .set("Irrtime", irr_time_mean)
            .set("Date", doy)
            .set("ModelRun", model_run_id)
        )

    # Apply extraction to hindcast images
    hindcast_data_images = hindcast_irrigation.map(extract_hindcast_means)

    # Convert to a list format for API
    hindcast_data_list = (
        hindcast_data_images.reduceColumns(
            ee.Reducer.toList(5), ["Date", "ModelRun",
                                   "Etc", "Irrvol", "Irrtime"]
        )
        .values()
        .get(0)
    )

    # Convert to pandas DataFrame
    hindcast_df = pd.DataFrame(
        hindcast_data_list.getInfo(),
        columns=[
            "date",
            "modelRun",
            "evaporation",
            "advisedWaterVolume",
            "advisedIrrigationTime",
        ],
    )

    # Round evaporation values
    hindcast_df["evaporation"] = hindcast_df["evaporation"].round(2)

    # Format dates and add required columns
    current_year = 2024  # Hardcoded for the example
    hindcast_df["date"] = pd.to_datetime(
        current_year * 1000 + hindcast_df["date"], format="%Y%j")
    hindcast_df["modelRun"] = hindcast_df["date"].dt.strftime(
        "%Y%m%d") + "_SAT" + current_day_str

    # Add display date with offset
    hindcast_df["datedisp"] = hindcast_df["date"] + pd.DateOffset(days=0)
    hindcast_df["dateDisplay"] = hindcast_df["datedisp"].dt.strftime(
        "%d-%m").map(lambda x: str(x)[-5:])

    # Format date column for API
    hindcast_df["date"] = pd.to_datetime(
        hindcast_df["date"] + pd.DateOffset(days=0), format="%d-%m-%Y").dt.strftime("%Y-%m-%d")

    # Add field ID
    hindcast_df["field"] = field_id

    # Select and order columns
    hindcast_df = hindcast_df[
        [
            "date",
            "dateDisplay",
            "modelRun",
            "evaporation",
            "advisedWaterVolume",
            "advisedIrrigationTime",
            "field",
        ]
    ]

    # Convert to JSON for API
    hindcast_json = hindcast_df.to_json(orient="records")
    print(hindcast_json)

    # Post hindcast data to API
    url_hindcast = "https://sosia.tahmo.org/api/hindcast/"
    headers_hindcast = {
        "Content-type": "application/json",
        "Authorization": "Basic bC52ZXJzY2h1cmVuQGZ1dHVyZXdhdGVyLm5sOnA3XlE1OTdNNmx3Wg==",
    }

    response_hindcast = requests.post(
        url_hindcast, headers=headers_hindcast, data=hindcast_json)
    print(response_hindcast.json())

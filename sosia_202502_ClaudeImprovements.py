import datetime
import ee
import pandas as pd
import requests

# Configuration constants
API_BASE_URL = "https://sosia.tahmo.org/api"
AUTH_HEADER = "Basic bC52ZXJzY2h1cmVuQGZ1dHVyZXdhdGVyLm5sOnA3XlE1OTdNNmx3Wg=="
TAHMO_AUTH_HEADER = "Basic ZnV0dXJld2F0ZXI6R2MzYVdMN3kyckRkR2Y3RQ=="
BUFFER_RADIUS = 30  # meters
# Constants for irrigation calculations
LOSS_RATE = 1.1  # Loss rate of drip irrigation system
APPLICATION_EFFICIENCY = 0.9  # Application efficiency
# Time constants
MILLIS_IN_DAY = 24 * 60 * 60 * 1000
HISTORICAL_START_DATE = "2010-01-01"
HISTORICAL_END_DATE = "2022-12-31"


def initialize_earth_engine():
    """Initialize the Google Earth Engine API"""
    ee.Authenticate()
    ee.Initialize()


def fetch_farmer_data():
    """Fetch farmer field data from the API"""
    url = f"{API_BASE_URL}/fields/"
    headers = {"Authorization": AUTH_HEADER}
    response = requests.get(url, headers=headers)
    return response.json()


def get_earth_engine_datasets():
    """Load required Earth Engine datasets"""
    return {
        "wapor_ret": ee.ImageCollection("FAO/WAPOR/2/L1_RET_E"),
        "gpm": ee.ImageCollection("NASA/GPM_L3/IMERG_V06"),
        "cfsv2": ee.ImageCollection("NOAA/CFSV2/FOR6H"),
        "dem": ee.Image("NASA/NASADEM_HGT/001")
    }


def create_point_geometry(longitude, latitude):
    """Create an Earth Engine point geometry from coordinates"""
    return ee.Geometry.Point(longitude, latitude)


def create_buffer(geometry, radius=BUFFER_RADIUS):
    """Create a buffer around a geometry"""
    return geometry.buffer(radius)


def clip_to_geometry(image, geometry):
    """Clip an image to a geometry"""
    return image.clip(geometry)


def get_crop_parameters(crop_type):
    """Get crop-specific parameters based on crop type"""
    crop_params = {
        "Habanero Peppers": {"kc": [0.6, 1.05, 1.05], "dev_days": 65, "mid_days": 40},
        "French Beans": {"kc": [0.5, 1.05, 1.05], "dev_days": 65, "mid_days": 40},
        "Lettuce": {"kc": [0.5, 1.05, 1.05], "dev_days": 85, "mid_days": 40},
        "Brassica": {"kc": [0.7, 1.05, 1.05], "dev_days": 65, "mid_days": 40},
        "Cucumbers": {"kc": [0.60, 1.00, 1.00], "dev_days": 60, "mid_days": 50},
        "Okra": {"kc": [0.30, 1.00, 0.90], "dev_days": 41, "mid_days": 25}
    }

    if crop_type not in crop_params:
        raise ValueError(f"Unknown crop type: {crop_type}")

    return crop_params[crop_type]


def calculate_growing_seasons(planting_date, harvest_date, crop_params):
    """Calculate the start and end dates for each growing season stage"""
    dev_days = crop_params["dev_days"]
    mid_days = crop_params["mid_days"]

    seasons = {
        "dev_start": planting_date.advance(0, "day"),
        "dev_end": planting_date.advance(dev_days, "day"),
        "mid_start": planting_date.advance(dev_days + 1, "day"),
        "mid_end": planting_date.advance(dev_days + mid_days, "day"),
        "end_start": planting_date.advance(dev_days + mid_days + 1, "day"),
        "end_end": harvest_date
    }

    return seasons


def process_wapor_data(wapor_collection, buffer_geometry, start_date, end_date):
    """Process WAPOR reference evapotranspiration data"""
    # Filter and clip WAPOR data
    wapor_filtered = wapor_collection.filterDate(start_date, end_date).map(
        lambda img: clip_to_geometry(img, buffer_geometry)
    )

    # Correct WAPOR data (divide by 10)
    wapor_corrected = wapor_filtered.map(
        lambda img: img.addBands(
            img.select("L1_RET_E").divide(10).rename("corrected")
        )
    )

    return wapor_corrected


def calculate_daily_average(collection, doys):
    """Calculate daily averages for a collection over a list of days of year"""
    def wapor_day_avg(doy):
        day_images = collection.select("corrected").filter(
            ee.Filter.calendarRange(start=doy, field="day_of_year")
        )
        mean_value = (
            ee.Image(day_images.mean())
            .multiply(100)
            .round()
            .divide(100)
            .set("DOY", doy)
        )
        return mean_value

    return ee.ImageCollection(doys.map(wapor_day_avg))


def calculate_stage_etc(wapor_days, kc_value, start_date, days_sequence):
    """Calculate crop evapotranspiration for a specific growth stage"""
    # Convert date to milliseconds
    start_millis = start_date.millis()
    end_millis = start_date.advance(days_sequence, "day").millis()

    # Get days in sequence
    def get_doy(date_millis):
        return ee.Number.parse(ee.Date(date_millis).format("DDD"))

    days = ee.List.sequence(start_millis, end_millis,
                            MILLIS_IN_DAY).map(get_doy)

    # Get WAPOR data for these days
    wapor_stage = ee.ImageCollection(days.map(
        lambda doy: wapor_days.filter(ee.Filter.eq("DOY", doy)).first()
    ))

    # Apply kc value
    etc_stage = wapor_stage.map(
        lambda img: img.multiply(kc_value)
                       .multiply(100)
                       .round()
                       .divide(100)
                       .copyProperties(img)
                       .set("DOY", img.get("DOY"))
    )

    # Add date property
    def add_date(img):
        doy = ee.Number.parse(img.get("system:index"))
        date = start_date.advance(doy, "day")
        return img.set("Date", date)

    return etc_stage.map(add_date)


def calculate_irrigation_needs(etc_collection, field_area, drip_params):
    """Calculate irrigation needs in volume and time"""
    flow_rate = drip_params["flow_rate"]
    init_time = drip_params["init_time"]

    def add_irrigation_bands(img):
        # Calculate irrigation volume in m³
        irr_m3 = (
            img.divide(1000)
            .multiply(field_area)
            .multiply(LOSS_RATE)
            .divide(APPLICATION_EFFICIENCY)
            .multiply(100)
            .round()
            .divide(100)
            .copyProperties(img)
            .set("DOY", img.get("DOY"))
        )

        # Calculate irrigation time in minutes
        irr_min = (
            img.divide(1000)
            .multiply(field_area)
            .multiply(LOSS_RATE)
            .divide(APPLICATION_EFFICIENCY)
            .multiply(60000)
            .divide(flow_rate)
            .round()
            .add(init_time)
            .copyProperties(img)
            .set("DOY", img.get("DOY"))
        )

        # Add new bands
        new_bands = ee.Image([irr_m3, irr_min])
        renamed_bands = new_bands.select(
            ["2 ETc in mm/day", "2 ETc in mm/day_1"],
            ["3 Irrigation needs in m3/day", "4 Irrigation time in min/day"]
        )

        return img.addBands(renamed_bands)

    return etc_collection.map(add_irrigation_bands)


def extract_point_data(collection, buffer_geometry):
    """Extract point data from an image collection"""
    def extract_values(img):
        doy = img.get("DOY")
        date_string = ee.Number(doy).format()
        today_string = datetime.date.today().strftime("%Y-%m-%d")
        model_run = ee.String(date_string).cat("_SAT").cat(today_string)

        # Extract values at point
        values = img.reduceRegion(
            reducer=ee.Reducer.mean(),
            geometry=buffer_geometry,
            scale=30
        )

        # Set properties with extracted values
        return img.set({
            "Etref": values.get("1 ETref in mm/day"),
            "Etc": values.get("2 ETc in mm/day"),
            "Irrvol": values.get("3 Irrigation needs in m3/day"),
            "Irrtime": values.get("4 Irrigation time in min/day"),
            "Date": doy,
            "ModelRun": model_run
        })

    reduced_imgs = collection.map(extract_values)

    nested_list = (
        reduced_imgs.reduceColumns(
            ee.Reducer.toList(5),
            ["Date", "ModelRun", "Etc", "Irrvol", "Irrtime"]
        )
        .values()
        .get(0)
    )

    return nested_list.getInfo()


def format_dataframe(data, field_id):
    """Format extracted data into a pandas DataFrame"""
    df = pd.DataFrame(
        data,
        columns=[
            "date", "modelRun", "evaporation",
            "advisedWaterVolume", "advisedIrrigationTime"
        ]
    )

    # Format date columns
    current_year = datetime.date.today().year
    df["date"] = pd.to_datetime(
        current_year * 1000 + df["date"], format="%Y%j")
    df["modelRun"] = df["date"].dt.strftime(
        "%Y%m%d") + "_SAT" + datetime.date.today().strftime("%Y-%m-%d")
    df["dateDisplay"] = df["date"].dt.strftime(
        "%d-%m").map(lambda x: str(x)[-5:])
    df["date"] = pd.to_datetime(
        df["date"], format="%d-%m-%Y").dt.strftime("%Y-%m-%d")
    df["field"] = field_id

    # Reorder columns
    return df[[
        "date", "dateDisplay", "modelRun", "evaporation",
        "advisedWaterVolume", "advisedIrrigationTime", "field"
    ]]


def post_data_to_api(endpoint, data):
    """Post data to API endpoint"""
    url = f"{API_BASE_URL}/{endpoint}/"
    headers = {
        "Content-type": "application/json",
        "Authorization": AUTH_HEADER
    }

    json_data = data.to_json(orient="records")
    response = requests.post(url, headers=headers, data=json_data)

    return response.json()


def calculate_field_area(drip_lines, drip_length, bed_width):
    """Calculate field area based on drip line parameters"""
    return drip_length * (bed_width / drip_lines)


def get_flow_rate(drip_length, emitter_spacing, flow_rate):
    """Calculate total flow rate"""
    return (drip_length / emitter_spacing) * flow_rate


def process_tahmo_data(point_string):
    """Fetch and process TAHMO weather station data"""
    url = f"https://smartirrigation.tahmo.org/sosia/nearestwetnessreport/{point_string}"
    headers = {"Authorization": TAHMO_AUTH_HEADER}

    response = requests.get(url, headers=headers)
    result = response.json()

    if not result:
        print("No TAHMO data available")
        return None

    # Process TAHMO data
    results = pd.json_normalize(result, "results", "distance")
    results2 = results.drop(
        columns=["_id", "Wetness", "FC", "RAM", "WP", "Eact"])
    results2["Day"] = pd.to_datetime(results2["Time"]).dt.strftime("%m/%d/%y")

    days = results2.pivot_table("Eref", "Day", "Station")
    distance = results.pivot_table("distance", "Station")

    # Calculate weighted values based on distance (IDW - Inverse Distance Weighting)
    if distance.count()[0] == 1:
        print("One station retrieved, value of nearest station retrieved")
        point_value = days
    elif distance.count()[0] == 2:
        print("Two stations retrieved")
        d1 = distance["distance"][0] + 28
        d2 = distance["distance"][1]
        v1 = days.iloc[:, 0]
        v2 = days.iloc[:, 1]
        point_value = ((v1 / d1) + (v2 / d2)) / ((1 / d1) + (1 / d2))
    else:
        print("Three stations retrieved")
        d1 = distance["distance"][0]
        d2 = distance["distance"][1]
        d3 = distance["distance"][2]
        v1 = days.iloc[:, 0]
        v2 = days.iloc[:, 1]
        v3 = days.iloc[:, 2]
        point_value = ((v1 / d1) + (v2 / d2) + (v3 / d3)) / \
                       ((1 / d1) + (1 / d2) + (1 / d3))

    tahmo_closest = days.iloc[:, 0].reset_index(drop=True)
    return tahmo_closest.tolist()


def calculate_hargreaves_et(cfsv2_collection, buffer_geometry, start_date, end_date):
    """Calculate reference evapotranspiration using Hargreaves method"""
    # Filter CFSv2 data
    cfs_filtered = cfsv2_collection.filterDate(start_date, end_date).map(
        lambda img: clip_to_geometry(img, buffer_geometry)
    )

    # Calculate net radiation
    def calc_rnet(image):
        dsw = image.select(
            "Downward_Short-Wave_Radiation_Flux_surface_6_Hour_Average")
        dlw = image.select(
            "Downward_Long-Wave_Radp_Flux_surface_6_Hour_Average")
        ulw = image.select("Upward_Long-Wave_Radp_Flux_surface_6_Hour_Average")
        usw = dsw.multiply(0.23).rename(
            "Upward_Short-Wave_Radiation_Flux_surface_6_Hour_Average")

        rnet = (
            dsw.subtract(usw)
            .add(dlw.subtract(ulw))
            .rename("Net_Radiation_6_Hour_Average")
        )

        # Calculate wind components
        u = image.select("u-component_of_wind_height_above_ground")
        v = image.select("v-component_of_wind_height_above_ground")
        wind_tot = u.pow(2).add(v.pow(2)).sqrt().rename("Wind_component")

        return image.addBands(rnet).addBands(wind_tot)

    cfs_with_radiation = cfs_filtered.map(calc_rnet)

    # Daily aggregations
    def aggregate_daily(collection, start, num_days, aggregator):
        def process_day(day_offset):
            day_start = start.advance(day_offset, "days")
            day_end = day_start.advance(1, "days")
            date_formatted = day_start.format("YYYY-MM-dd")
            doy = ee.Number.parse(day_start.format("DDD"))

            return (
                collection.filterDate(day_start, day_end)
                [aggregator]()
                .set("Date", date_formatted)
                .set("DOY", doy)
            )

        day_sequence = ee.List.sequence(0, num_days.subtract(1))
        return ee.ImageCollection(day_sequence.map(process_day))

    # Calculate daily means, min and max
    num_days = end_date.difference(start_date, "days")
    daily_mean = aggregate_daily(
        cfs_with_radiation, start_date, num_days, "mean")

    temp_collection = cfs_filtered.select("Temperature_height_above_ground")
    temp_max = aggregate_daily(temp_collection, start_date, num_days, "max").map(
        lambda img: img.select(
            "Temperature_height_above_ground").rename("Max Temperature")
    )
    temp_min = aggregate_daily(temp_collection, start_date, num_days, "min").map(
        lambda img: img.select(
            "Temperature_height_above_ground").rename("Min Temperature")
    )

    # Join temperature and radiation data
    join_filter = ee.Filter.equals(leftField="Date", rightField="Date")
    simple_join = ee.Join.inner()

    inner_join_min = ee.ImageCollection(
        simple_join.apply(temp_max, temp_min, join_filter))
    temp_max_min = inner_join_min.map(
        lambda feature: ee.Image.cat(feature.get(
            "primary"), feature.get("secondary"))
    )

    inner_join_max = ee.ImageCollection(
        simple_join.apply(daily_mean, temp_max_min, join_filter))
    cfs_total = inner_join_max.map(
        lambda feature: ee.Image.cat(feature.get(
            "primary"), feature.get("secondary"))
    )

    # Calculate ET using modified Penman-Monteith
    def calculate_et0(image):
        # Extract required variables
        t_min = image.select("Min Temperature").subtract(273.15).rename("Tmin")
        t_max = image.select("Max Temperature").subtract(273.15).rename("Tmax")
        pressure = image.select("Pressure_surface").divide(
            1000).rename("Atm pressure")
        wind = image.select("Wind_component").multiply(0.75).rename("Wind")
        r_net = image.select("Net_Radiation_6_Hour_Average").multiply(
            0.0864).rename("Rnet")

        # Calculate mean temperature
        t_mean = t_min.add(t_max).divide(2)

        # Constants
        exp_base = pressure.multiply(0).add(2.71828)

        # Calculate delta (slope of saturation vapor pressure curve)
        delta_term1 = t_mean.multiply(17.27).divide(t_mean.add(237.3))
        delta_term2 = exp_base.pow(delta_term1).multiply(0.6108)
        delta_term3 = t_mean.add(237.3).pow(2)
        delta = delta_term2.multiply(4098).divide(delta_term3)

        # Calculate saturation vapor pressure
        e0_min_term = t_min.multiply(17.27).divide(t_min.add(237.3))
        e0_min = exp_base.pow(e0_min_term).multiply(0.6108)
        e0_max_term = t_max.multiply(17.27).divide(t_max.add(237.3))
        e0_max = exp_base.pow(e0_max_term).multiply(0.6108)
        e_s = e0_min.add(e0_max).divide(2)

        # Calculate actual vapor pressure
        rh_max = (
            image.select("Specific_humidity_height_above_ground")
            .multiply(pressure)
            .multiply(1.6077717)
            .divide(e0_min)
            .rename("RHmax")
        )
        rh_min = (
            image.select("Specific_humidity_height_above_ground")
            .multiply(pressure)
            .multiply(1.6077717)
            .divide(e0_max)
            .rename("RHmin")
        )
        e_a = e0_min.multiply(rh_max).add(e0_max.multiply(rh_min)).divide(2)

        # Calculate ET0 using modified Penman-Monteith
        part1 = delta.multiply(r_net).multiply(0.408)
        psy = pressure.multiply(0.001).divide(1.53634)
        part2 = psy.multiply(900).divide(t_mean.add(273))
        part3 = wind.multiply(e_s.subtract(e_a))
        part4 = psy.multiply(wind.multiply(0.34).add(1)).add(delta)

        et0 = part1.add(part2.multiply(part3)).divide(
            part4).rename("1 ETref in mm/day")

        return image.addBands(et0)

    # Calculate ET0
    return cfs_total.map(calculate_et0).select("1 ETref in mm/day")


def main():
    """Main function to process irrigation advice"""
    # Initialize Earth Engine
    initialize_earth_engine()

    # Fetch farmer data
    farmers = fetch_farmer_data()

    # Get Earth Engine datasets
    datasets = get_earth_engine_datasets()

    # Process for each farmer
    for farmer in farmers:
        print(f"Processing field: {farmer['id']}")

        # Extract farmer data
        field_id = farmer["id"]
        field_name = farmer["name"]
        lat = float(farmer["latitude"])
        lon = float(farmer["longitude"])

        # Create point geometry
        point = create_point_geometry(lon, lat)
        point_string = f"{lon},{lat}"
        buffer_geo = create_buffer(point)

        # Drip irrigation parameters
        drip_lines = farmer["numberOfDriplines"]
        drip_length = farmer["lengthOfDriplines"]
        bed_width = farmer["bedWidth"]
        emitter_spacing = farmer["emitterSpacing"]
        emitter_flow_rate = farmer["emitterFlowRate"]
        init_time = farmer["initialisationTime"]

        # Calculate derived parameters
        flow = get_flow_rate(drip_length, emitter_spacing, emitter_flow_rate)
        field_area = calculate_field_area(drip_lines, drip_length, bed_width)

        drip_params = {
            "flow_rate": flow,
            "init_time": init_time
        }

        # Crop information
        crop_type = farmer["cropSpecific"]["cropType"]
        planting_date = ee.Date(farmer["cropSpecific"]["plantingDate"])
        harvest_date = ee.Date(
            farmer["cropSpecific"]["lastExpectedHarvestingDate"])

        # Get crop parameters and calculate seasons
        crop_params = get_crop_parameters(crop_type)
        seasons = calculate_growing_seasons(
            planting_date, harvest_date, crop_params)

        # Process WAPOR data
        wapor_corrected = process_wapor_data(
            datasets["wapor_ret"],
            buffer_geo,
            ee.Date(HISTORICAL_START_DATE),
            ee.Date(HISTORICAL_END_DATE)
        )

        # Calculate daily averages for all days of year
        doys = ee.List.sequence(1, 365)
        wapor_daily_avg = calculate_daily_average(wapor_corrected, doys)

        # Get days in growing season
        season_start_millis = seasons["dev_start"].millis()
        season_end_millis = seasons["end_end"].millis()

        def get_doy(date_millis):
            return ee.Number.parse(ee.Date(date_millis).format("DDD"))

        season_days = ee.List.sequence(
            season_start_millis,
            season_end_millis,
            MILLIS_IN_DAY
        ).map(get_doy)

        # Calculate reference ET for growing season
        wapor_ref = ee.ImageCollection(season_days.map(
            lambda doy: wapor_daily_avg.filter(
                ee.Filter.eq("DOY", doy)).first()
        ))

        # Format reference ET
        etc_ref = wapor_ref.map(
            lambda img: img.multiply(100)
                           .round()
                           .divide(100)
                           .copyProperties(img)
                           .set("DOY", img.get("DOY"))
        )

        # Add date property
        def add_date_from_planting(img):
            doy = ee.Number.parse(img.get("system:index"))
            date = planting_date.advance(doy, "day")
            return img.set("Date", date)

        total_ref_date = etc_ref.map(add_date_from_planting).select(
            ["corrected"], ["1 ETref in mm/day"]
        )

        # Calculate ETc for each growth stage
        kc_values = crop_params["kc"]

        # Development stage
        etc_dev = calculate_stage_etc(
            wapor_daily_avg,
            kc_values[0],
            seasons["dev_start"],
            crop_params["dev_days"]
        )

        # Mid stage
        etc_mid = calculate_stage_etc(
            wapor_daily_avg,
            kc_values[1],
            seasons["mid_start"],
            crop_params["mid_days"]
        )

        # End stage
        end_days = harvest_date.difference(seasons["end_start"], "day")
        etc_end = calculate_stage_etc(
            wapor_daily_avg,
            kc_values[2],
            seasons["end_start"],
            end_days
        )

        # Merge stages
        etc_merged = etc_dev.merge(etc_mid).merge(etc_end)
        etc = etc_merged.select(["corrected"], ["2 ETc in mm/day"])

        # Calculate irrigation needs
        total = calculate_irrigation_needs(etc, field_area, drip_params)

        # Join with reference ET
        join_filter = ee.Filter.equals(leftField="Date", rightField="Date")
        simple_join = ee.Join.inner()

        joined = ee.ImageCollection(simple_join.apply(
            total_ref_date, total, join_filter
        ))

        # Create merged collection
        historical_crop_schedule = joined.map(
            lambda feature: ee.Image.cat(
                feature.get("primary"), feature.get("secondary")
            )
        )

        # Extract point data
        historical_data = extract_point_data(
            historical_crop_schedule, buffer_geo)

        # Format data
        historical_df = format_dataframe(historical_data, field_id)

        # Post historical data
        historical_response = post_data_to_api(
            "seasonal_schedule", historical_df)
        print(f"Historical data posted for {field_id}")

        # Process hindcast data
        today = ee.Date(datetime.date.today().strftime("%Y-%m-%d"))
        week_ago = today.advance(-9, "day")
        yesterday = today.advance(-2, "day")

        # Try to get TAHMO data
        tahmo_data = process_tahmo_data(point_string)

        if tahmo_data:
            # Use TAHMO data for hindcast
            week_ago_recent = today.advance(-6, "day")
            today_millis = today.millis()
            week_ago_millis = week_ago_recent.millis()

            recent_days = ee.List.sequence(
                week_ago_millis, today_millis, MILLIS_IN_DAY
            ).map(get_doy)

            tahmo_hind = ee.ImageCollection(recent_days.map(
                lambda doy: wapor_daily_avg.filter(
                    ee.Filter.eq("DOY", doy)).first()
            ))

            # Create images from TAHMO data
            zipped_list = tahmo_hind.toList(1000).zip(tahmo_data)

            def create_scaled_image(list_item):
                list_item = ee.List(list_item)
                img = ee.Image(list_item.get(0))
                scale = list_item.getNumber(1)
                scaled_image = img.multiply(0).add(
                    scale).double().rename("1 ETref in mm/day")
                return scaled_image.copyProperties(img, img.propertyNames())

            tahmo_list = zipped_list.map(create_scaled_image)
            et0_collection = ee.ImageCollection.fromImages(tahmo_list)

            print("Using TAHMO data for hindcast")
        else:
            # Calculate ET0 using Hargreaves when no TAHMO data available
            et0_collection = calculate_hargreaves_et(
                datasets["cfsv2"],
                buffer_geo,
                week_ago,
                yesterday
            )

            print("Using CFSv2 data for hindcast")

        # Create Kc collections for each growth stage
        kc_dev = etc_dev.map(
            lambda img: img.multiply(0)
                           .add(kc_values[0])
                           .copyProperties(img)
                           .set("DOY", img.get("DOY"))
        )

        kc_mid = etc_mid.map(
            lambda img: img.multiply(0)
                           .add(kc_values[1])
                           .copyProperties(img)
                           .set("DOY", img.get("DOY"))
        )

        kc_end = etc_end.map(
            lambda img: img.multiply(0

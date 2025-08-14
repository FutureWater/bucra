"""
Smart Irrigation Advice System using Google Earth Engine and TAHMO data
--------------------------------------------------------------------

This script calculates irrigation recommendations for farmers based on field characteristics,
crop type, climate data, and evapotranspiration models. It uses Google Earth Engine (GEE) for
remote sensing data and TAHMO weather station data when available.

Main steps of the program:
1. Initialize GEE and authenticate
2. Fetch farmer field data from API (field dimensions, crop type, planting dates)
3. For each farmer field:
   a. Set crop-specific parameters (kc values and growth stage durations)
   b. Calculate field-specific parameters (irrigation system details)
   c. Process historical evapotranspiration data from WAPOR dataset
   d. Calculate irrigation needs for entire growing season ("Historical Crop Schedule")
   e. Calculate recent and current irrigation needs ("Hindcast"):
      i. Try to use TAHMO weather station data if available
      ii. Fall back to CFSv2 forecast data if no TAHMO stations nearby
   f. Post results back to API
"""

import datetime as dt
from datetime import date
import ee
import pandas as pd
import requests


def initialize_gee():
    """Initialize Google Earth Engine with authentication."""
    cloud_project = 'ee-ameliafdezrodriguez2' # replace with your project id

    try:
        ee.Initialize(
            project=cloud_project,
            opt_url='https://earthengine-highvolume.googleapis.com')
        
    except:
        ee.Authenticate()
        ee.Initialize(
            project=cloud_project,
            opt_url='https://earthengine-highvolume.googleapis.com')



def fetch_farmer_data(api_url, headers):
    """Fetch farmer field data from the API."""
    response = requests.request("GET", api_url, headers=headers, data={})
    return response.json()


def create_point_geometry(longitude, latitude):
    """Create an EE point geometry from coordinates."""
    return ee.Geometry.Point(longitude, latitude)


def get_crop_parameters(crop_type):
    """Return crop-specific coefficients and growth stages."""
    # Default values (will be overwritten)
    params = {
        "kc1": 0.5,  # Initial stage crop coefficient
        "kc2": 1.0,  # Mid-season stage crop coefficient
        "kc3": 0.9,  # Late season crop coefficient
        "dev_days": 60,   # Days in development stage
        "mid_days": 40,   # Days in mid-season stage
    }

    # Crop-specific parameters
    if crop_type == "Habanero Peppers":
        params.update({"kc1": 0.6, "kc2": 1.05, "kc3": 1.05,
                      "dev_days": 65, "mid_days": 40})
    elif crop_type == "French Beans":
        params.update({"kc1": 0.5, "kc2": 1.05, "kc3": 1.05,
                      "dev_days": 65, "mid_days": 40})
    elif crop_type == "Lettuce":
        params.update({"kc1": 0.5, "kc2": 1.05, "kc3": 1.05,
                      "dev_days": 85, "mid_days": 40})
    elif crop_type == "Brassica":
        params.update({"kc1": 0.7, "kc2": 1.05, "kc3": 1.05,
                      "dev_days": 65, "mid_days": 40})
    elif crop_type == "Cucumbers":
        params.update({"kc1": 0.60, "kc2": 1.00, "kc3": 1.00,
                      "dev_days": 60, "mid_days": 50})
    elif crop_type == "Okra":
        params.update({"kc1": 0.30, "kc2": 1.00, "kc3": 0.90,
                      "dev_days": 41, "mid_days": 25})

    return params


def calculate_growth_stages(planting_date, harvest_date, dev_days, mid_days):
    """Calculate the start and end dates for each growth stage."""
    stages = {
        "dev_start": planting_date.advance(0, "day"),
        "dev_end": planting_date.advance(dev_days, "day"),
        "mid_start": planting_date.advance(dev_days + 1, "day"),
        "mid_end": planting_date.advance(dev_days + mid_days, "day"),
        "end_start": planting_date.advance(dev_days + mid_days + 1, "day"),
        "end_end": harvest_date
    }
    return stages


def calculate_irrigation_params(drip_lines, drip_length, bed_width, emitter_spacing, flow_rate, init_time):
    """Calculate irrigation system parameters."""
    flow = (drip_length / emitter_spacing) * \
        flow_rate  # Total flow rate of the system
    area = drip_length * (bed_width / drip_lines)       # Irrigated area
    loss_rate = 1.1    # Loss rate of drip irrigation system (constant)
    efficiency = 0.9   # Application efficiency (constant)

    return {
        "flow": flow,
        "area": area,
        "loss_rate": loss_rate,
        "efficiency": efficiency,
        "init_time": init_time
    }


def load_datasets():
    """Load required Earth Engine datasets."""
    return {
        "wapor_ret": ee.ImageCollection("FAO/WAPOR/2/L1_RET_E"),
        "gpm": ee.ImageCollection("NASA/GPM_L3/IMERG_V06"),
        "cfsv2": ee.ImageCollection("NOAA/CFSV2/FOR6H"),
        "dem": ee.Image("NASA/NASADEM_HGT/001")
    }


def clip_image_to_buffer(buffer_geometry, image):
    """Function that clips images to a buffer geometry."""
    return image.clip(buffer_geometry)


def wapor_correction(img):
    """Function to correct WAPOR RET values (divide by 10)."""
    corrected = img.select("L1_RET_E").divide(10).rename("corrected")
    return img.addBands(corrected)


def wapor_daily_average(wapor_corr,day_of_year):
    """Function to calculate daily averages of WAPOR data for specific days of year."""
    day_images = wapor_corr.select("corrected").filter(
            ee.Filter.calendarRange(start=day_of_year, field="day_of_year")
        )
    mean_wapor = (
            ee.Image(day_images.mean())
            .multiply(100)
            .round()
            .divide(100)
            .set("DOY", day_of_year)
        )
    return mean_wapor


def calculate_historical_et(stages, wapor_day_avg_func, kc_values, planting_date, millis_per_day):
    """Calculate historical evapotranspiration for the entire growing season."""
    # Calculate days in the growing season (from planting to harvest)
    start_millis = stages["dev_start"].millis()
    end_millis = stages["end_end"].millis()

    # Convert date millis to days of year
    def millis_to_doy(date_millis):
        return ee.Number.parse(ee.Date(date_millis).format("DDD"))

    # Get all days of year in the growing season
    growing_season_days = ee.List.sequence(
        start_millis, end_millis, millis_per_day).map(millis_to_doy)

    # Get WAPOR reference ET for all days in the growing season
    wapor_ref = ee.ImageCollection(growing_season_days.map(wapor_day_avg_func))

    # Process reference ET values
    def process_et_ref(image):
        return (
            image.multiply(100)
            .round()
            .divide(100)
            .copyProperties(image)
            .set("DOY", image.get("DOY"))
        )

    et_ref = wapor_ref.map(process_et_ref)

    # Add date property to each image
    def add_date_property(img):
        day_offset = ee.Number.parse(img.get("system:index"))
        date_val = planting_date.advance(day_offset, "day")
        return img.set("Date", date_val)

    et_ref_with_date = et_ref.map(add_date_property)

    # Select and rename the band
    et_ref_renamed = et_ref_with_date.select(
        ["corrected"], ["1 ETref in mm/day"])

    # Calculate ET for each growth stage using the appropriate kc value
    stages_et = calculate_stage_et(
        stages, wapor_day_avg_func, kc_values, planting_date, millis_per_day)

    return {
        "et_ref": et_ref_renamed,
        "et_stages": stages_et
    }


def calculate_stage_et(stages, wapor_day_avg_func, kc_values, planting_date, millis_per_day):
    """Calculate ET for each growth stage."""
    # Development stage
    dev_et = calculate_stage_specific_et(
        stages["dev_start"],
        stages["dev_end"],
        wapor_day_avg_func,
        kc_values["kc1"],
        planting_date,
        millis_per_day
    )

    # Mid-season stage
    mid_et = calculate_stage_specific_et(
        stages["mid_start"],
        stages["mid_end"],
        wapor_day_avg_func,
        kc_values["kc2"],
        stages["mid_start"],
        millis_per_day
    )

    # End stage
    end_et = calculate_stage_specific_et(
        stages["end_start"],
        stages["end_end"],
        wapor_day_avg_func,
        kc_values["kc3"],
        stages["end_start"],
        millis_per_day
    )

    # Merge the stages
    et_merged = dev_et.merge(mid_et).merge(end_et)
    et_final = et_merged.select(["corrected"], ["2 ETc in mm/day"])

    return et_final


def calculate_stage_specific_et(start_date, end_date, wapor_day_avg_func, kc_value, reference_date, millis_per_day):
    """Calculate ET for a specific growth stage."""
    # Get time range in milliseconds
    start_millis = start_date.millis()
    end_millis = end_date.millis()

    # Convert to days of year
    def millis_to_doy(date_millis):
        return ee.Number.parse(ee.Date(date_millis).format("DDD"))

    days = ee.List.sequence(start_millis, end_millis,
                            millis_per_day).map(millis_to_doy)

    # Get WAPOR data for these days
    wapor_data = ee.ImageCollection(days.map(wapor_day_avg_func))

    # Apply kc value to WAPOR data
    def apply_kc(image):
        return (
            image.multiply(kc_value)
            .multiply(100)
            .round()
            .divide(100)
            .copyProperties(image)
            .set("DOY", image.get("DOY"))
        )

    et_with_kc = wapor_data.map(apply_kc)

    # Add date property
    def add_date_to_image(img):
        day_offset = ee.Number.parse(img.get("system:index"))
        date_val = reference_date.advance(day_offset, "day")
        return img.set("Date", date_val)

    return et_with_kc.map(add_date_to_image)


def calculate_irrigation_needs(et_data, irr_params):
    """Calculate irrigation needs from ET data."""
    def calculate_water_and_time(image):
        # Calculate irrigation volume (m³/day)
        irr_m3 = (
            image.divide(1000)                # Convert mm to m
            .multiply(irr_params["area"])     # Multiply by area to get volume
            .multiply(irr_params["loss_rate"])  # Account for losses
            .divide(irr_params["efficiency"])  # Account for efficiency
            .multiply(100).round().divide(100)  # Round to 2 decimals
            .copyProperties(image)
            .set("DOY", image.get("DOY"))
        )

        # Calculate irrigation time (minutes/day)
        irr_min = (
            image.divide(1000)                # Convert mm to m
            .multiply(irr_params["area"])     # Multiply by area to get volume
            .multiply(irr_params["loss_rate"])  # Account for losses
            .divide(irr_params["efficiency"])  # Account for efficiency
            .multiply(60000)                  # Convert hours to minutes
            .divide(irr_params["flow"])       # Divide by flow rate
            .round()
            .add(irr_params["init_time"])     # Add initialization time
            .copyProperties(image)
            .set("DOY", image.get("DOY"))
        )

        # Create new bands with the calculated values
        vol_band = ee.Image([irr_m3])
        time_band = ee.Image([irr_min])

        vol_band_renamed = vol_band.select(
            ["2 ETc in mm/day"], ["3 Irrigation needs in m3/day"]
        )
        time_band_renamed = time_band.select(
            ["2 ETc in mm/day"], ["4 Irrigation time in min/day"]
        )

        return image.addBands(vol_band_renamed).addBands(time_band_renamed)

    # Apply the calculations to the ET data
    return et_data.map(calculate_water_and_time)


def process_tahmo_data(point_coords, headers):
    """Process TAHMO weather station data if available."""
    url = f"https://smartirrigation.tahmo.org/sosia/nearestwetnessreport/{point_coords}"
    response = requests.get(url, headers=headers, data={})
    result = response.json()

    if not result:
        return None

    # Process TAHMO data if available
    results = pd.json_normalize(result, "results", "distance")
    results = results.drop(
        columns=["_id", "Wetness", "FC", "RAM", "WP", "Eact"])
    results["Day"] = pd.to_datetime(results["Time"]).dt.strftime("%m/%d/%y")
    results = results.drop(columns=["Time"])

    days = results.pivot_table("Eref", "Day", "Station")
    distance = results.pivot_table("distance", "Station")

    # Inverse Distance Weighting based on number of stations
    if distance.count()[0] == 1:
        print("One station retrieved, value of nearest station used")
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

    tahmo_data = days.iloc[:, 0]
    return tahmo_data.reset_index(drop=True).tolist()


def calculate_hindcast_et(today, datasets, clip_func, wapor_day_avg_func, tahmo_data=None):
    """Calculate hindcast ET using either TAHMO data or CFSv2 data."""
    millis_per_day = 24 * 60 * 60 * 1000

    if tahmo_data:
        # Use TAHMO data
        week_ago = today.advance(-6, "day")
        today_millis = today.millis()
        week_ago_millis = week_ago.millis()

        def millis_to_doy(date_millis):
            return ee.Number.parse(ee.Date(date_millis).format("DDD"))

        hindcast_days = ee.List.sequence(
            week_ago_millis, today_millis, millis_per_day).map(millis_to_doy)
        tahmo_hindcast = ee.ImageCollection(
            hindcast_days.map(wapor_day_avg_func))

        # Combine WAPOR structure with TAHMO values
        max_elements = 1000
        zipped_list = tahmo_hindcast.toList(max_elements).zip(tahmo_data)

        def create_tahmo_image(list_item):
            list_item = ee.List(list_item)
            img = ee.Image(list_item.get(0))
            scale = list_item.getNumber(1)
            scaled_image = img.multiply(0).add(scale)
            scaled_image = scaled_image.double().rename("1 ETref in mm/day")
            return scaled_image.copyProperties(img, img.propertyNames())

        tahmo_list = zipped_list.map(create_tahmo_image)
        return ee.ImageCollection.fromImages(tahmo_list)
    else:
        # Use CFSv2 data
        print("No TAHMO data available, using CFSv2")

        # Dates for hindcast (from 9 days ago to 2 days ago due to data availability)
        week_ago = today.advance(-9, "day")
        week_ago_fmt = ee.Date(week_ago.format("YYYY-MM-dd"))
        yesterday = today.advance(-2, "day")
        yesterday_fmt = ee.Date(yesterday.format("YYYY-MM-dd"))

        # Filter CFSv2 data for the hindcast period
        cfsv2_data = datasets["cfsv2"].filterDate(week_ago_fmt, yesterday_fmt)
        cfsv2_clipped = ee.ImageCollection(cfsv2_data.map(clip_func))

        # Calculate net radiation
        def calculate_radiation(image):
            dsw = image.select(
                "Downward_Short-Wave_Radiation_Flux_surface_6_Hour_Average")
            dlw = image.select(
                "Downward_Long-Wave_Radp_Flux_surface_6_Hour_Average")
            ulw = image.select(
                "Upward_Long-Wave_Radp_Flux_surface_6_Hour_Average")
            usw = dsw.multiply(0.23).rename(
                "Upward_Short-Wave_Radiation_Flux_surface_6_Hour_Average")

            rnet = dsw.subtract(usw).add(dlw.subtract(ulw)).rename(
                "Net_Radiation_6_Hour_Average")

            # Calculate wind component
            u = image.select("u-component_of_wind_height_above_ground")
            v = image.select("v-component_of_wind_height_above_ground")
            wind_tot = u.pow(2).add(v.pow(2)).sqrt().rename("Wind_component")

            return image.addBands(rnet).addBands(wind_tot)

        # Calculate all components for ET calculation
        cfsv2_with_rad = cfsv2_clipped.map(calculate_radiation)
        cfsv2_temp = cfsv2_clipped.select("Temperature_height_above_ground")

        # Calculate daily values
        num_days = yesterday_fmt.difference(week_ago_fmt, "days")

        def calculate_daily_mean(day_offset):
            start = week_ago_fmt.advance(day_offset, "days")
            end = start.advance(1, "days")
            return (
                cfsv2_with_rad.filterDate(start, end)
                .mean()
                .set("Date", start.format("YYYY-MM-dd"))
                .set("DOY", ee.Number.parse(start.format("DDD")))
            )

        daily_mean = ee.ImageCollection(
            ee.List.sequence(0, num_days.subtract(1)).map(calculate_daily_mean)
        )

        # Calculate daily max temperature
        def calculate_daily_max_temp(day_offset):
            start = week_ago_fmt.advance(day_offset, "days")
            end = start.advance(1, "days")
            return (
                cfsv2_temp.filterDate(start, end)
                .max()
                .set("Date", start.format("YYYY-MM-dd"))
            )

        daily_max_temp = ee.ImageCollection(
            ee.List.sequence(0, num_days.subtract(1)).map(
                calculate_daily_max_temp)
        )

        max_temp = daily_max_temp.map(
            lambda image: image.select(
                "Temperature_height_above_ground").rename("Max Temperature")
        )

        # Calculate daily min temperature
        def calculate_daily_min_temp(day_offset):
            start = week_ago_fmt.advance(day_offset, "days")
            end = start.advance(1, "days")
            return (
                cfsv2_temp.filterDate(start, end)
                .min()
                .set("Date", start.format("YYYY-MM-dd"))
            )

        daily_min_temp = ee.ImageCollection(
            ee.List.sequence(0, num_days.subtract(1)).map(
                calculate_daily_min_temp)
        )

        min_temp = daily_min_temp.map(
            lambda image: image.select(
                "Temperature_height_above_ground").rename("Min Temperature")
        )

        # Join all the data
        join_filter = ee.Filter.equals(leftField="Date", rightField="Date")
        simple_join = ee.Join.inner()

        # Join min and max temperatures
        temp_join = ee.ImageCollection(
            simple_join.apply(max_temp, min_temp, join_filter))
        temp_combined = temp_join.map(
            lambda feature: ee.Image.cat(feature.get(
                "primary"), feature.get("secondary"))
        )

        # Join with daily mean data
        final_join = ee.ImageCollection(simple_join.apply(
            daily_mean, temp_combined, join_filter))
        combined_data = final_join.map(
            lambda feature: ee.Image.cat(feature.get(
                "primary"), feature.get("secondary"))
        )

        # Calculate ET0 using the Penman-Monteith equation
        def calculate_et0(image):
            # Get required bands
            t_min_k = image.select("Min Temperature")
            t_max_k = image.select("Max Temperature")
            wind = image.select("Wind_component")
            spec_humidity = image.select(
                "Specific_humidity_height_above_ground")
            pressure_pa = image.select("Pressure_surface")
            rad_net = image.select("Net_Radiation_6_Hour_Average")

            # Convert units
            t_min = t_min_k.subtract(273.15).rename("Tmin")
            t_max = t_max_k.subtract(273.15).rename("Tmax")
            pressure = pressure_pa.divide(1000).rename("Atm pressure")
            wind_adj = wind.multiply(0.75).rename("Wind")
            rad_net_adj = rad_net.multiply(0.0864).rename("Rnet")

            # For exponential calculations
            exp_base = pressure.multiply(0).add(2.71828)

            # Calculate mean temperature
            t_mean = t_min.add(t_max).divide(2)

            # Calculate delta (slope of vapor pressure curve)
            delta_term1 = t_mean.multiply(17.27).divide(t_mean.add(237.3))
            delta_term2 = exp_base.pow(delta_term1).multiply(0.6108)
            delta_term3 = t_mean.add(237.3).pow(2)
            delta = delta_term2.multiply(4098).divide(delta_term3)

            # Calculate saturation vapor pressure
            e0_min_term = t_min.multiply(17.27).divide(t_min.add(237.3))
            e0_min = exp_base.pow(e0_min_term).multiply(0.6108)
            e0_max_term = t_max.multiply(17.27).divide(t_max.add(237.3))
            e0_max = exp_base.pow(e0_max_term).multiply(0.6108)
            es = e0_min.add(e0_max).divide(2)

            # Calculate actual vapor pressure
            rh_max = (
                spec_humidity.multiply(pressure)
                .multiply(1.6077717)
                .divide(e0_min)
                .rename("RHmax")
            )
            rh_min = (
                spec_humidity.multiply(pressure)
                .multiply(1.6077717)
                .divide(e0_max)
                .rename("RHmin")
            )
            ea = e0_min.multiply(rh_max).add(e0_max.multiply(rh_min)).divide(2)

            # Calculate ET0 using Penman-Monteith
            part1 = delta.multiply(rad_net_adj).multiply(0.408)
            psy = pressure.multiply(0.001).divide(1.53634)
            part2 = psy.multiply(900).divide(t_mean.add(273))
            part3 = wind_adj.multiply(es.subtract(ea))
            part4 = psy.multiply(wind_adj.multiply(0.34).add(1)).add(delta)

            et0 = part1.divide(part4).add(part2.multiply(
                part3).divide(part4)).rename("1 ETref in mm/day")

            return image.addBands(et0)

        # Calculate ET0 for all days
        et_data = combined_data.map(calculate_et0)
        return et_data.select("1 ETref in mm/day")


def create_kc_images(stages, wapor_collections, kc_values):
    """Create images with Kc values for each growth stage."""
    # Development stage
    def create_kc_dev(image):
        return (
            image.multiply(0)
            .add(kc_values["kc1"])
            .copyProperties(image)
            .set("DOY", image.get("DOY"))
        )

    kc_dev = wapor_collections["dev"].map(create_kc_dev)

    def add_date_to_dev(image):
        doy = ee.Number.parse(image.get("system:index"))
        date = stages["dev_start"].advance(doy, "day").format("YYYY-MM-dd")
        return image.set("Date", date)

    kc_dev_with_date = kc_dev.map(add_date_to_dev)

    # Mid stage
    def create_kc_mid(image):
        return (
            image.multiply(0)
            .add(kc_values["kc2"])
            .copyProperties(image)
            .set("DOY", image.get("DOY"))
        )

    kc_mid = wapor_collections["mid"].map(create_kc_mid)

    def add_date_to_mid(image):
        doy = ee.Number.parse(image.get("system:index"))
        date = stages["mid_start"].advance(doy, "day").format("YYYY-MM-dd")
        return image.set("Date", date)

    kc_mid_with_date = kc_mid.map(add_date_to_mid)

    # End stage
    def create_kc_end(image):
        return (
            image.multiply(0)
            .add(kc_values["kc3"])
            .copyProperties(image)
            .set("DOY", image.get("DOY"))
        )

    kc_end = wapor_collections["end"].map(create_kc_end)

    def add_date_to_end(image):
        doy = ee.Number.parse(image.get("system:index"))
        date = stages["end_start"].advance(doy, "day").format("YYYY-MM-dd")
        return image.set("Date", date)

    kc_end_with_date = kc_end.map(add_date_to_end)

    # Merge all stages
    kc_merged = kc_dev_with_date.merge(
        kc_mid_with_date).merge(kc_end_with_date)
    return kc_merged.select(["corrected"], ["2 Daily Kc"])


def calculate_hindcast_irrigation(et_ref, kc_images, irr_params, join_by_date=True):
    """Calculate irrigation needs for the hindcast period."""
    # Join ET reference with Kc values
    join_filter = ee.Filter.equals(
        leftField="Date" if join_by_date else "DOY",
        rightField="Date" if join_by_date else "DOY"
    )
    simple_join = ee.Join.inner()

    kc_join = ee.ImageCollection(
        simple_join.apply(et_ref, kc_images, join_filter))

    def combine_images(feature):
        return ee.Image.cat(feature.get("primary"), feature.get("secondary"))

    joined_images = kc_join.map(combine_images)

    # Calculate crop ET
    def calculate_etc(image):
        et_ref = image.select("1 ETref in mm/day")
        kc = image.select("2 Daily Kc")
        etc = et_ref.multiply(kc).rename("3 ETc in mm/day")
        return image.addBands(etc)

    etc_images = joined_images.map(calculate_etc)

    # Calculate irrigation needs
    def calculate_irrigation(image):
        etc = image.select("3 ETc in mm/day")

        # Calculate irrigation volume (m³/day)
        irr_m3 = (
            etc.divide(1000)                 # Convert mm to m
            .multiply(irr_params["area"])    # Multiply by area to get volume
            .multiply(irr_params["loss_rate"])  # Account for losses
            .divide(irr_params["efficiency"])  # Account for efficiency
            .multiply(100).round().divide(100)  # Round to 2 decimals
            .copyProperties(image)
        )

        # Calculate irrigation time (minutes/day)
        irr_min = (
            etc.divide(1000)                 # Convert mm to m
            .multiply(irr_params["area"])    # Multiply by area to get volume
            .multiply(irr_params["loss_rate"])  # Account for losses
            .divide(irr_params["efficiency"])  # Account for efficiency
            .multiply(60000)                 # Convert hours to minutes
            .divide(irr_params["flow"])      # Divide by flow rate
            .round()
            .add(irr_params["init_time"])    # Add initialization time
            .copyProperties(image)
        )

        # Create new bands with the calculated values
        bands = ee.Image([irr_m3, irr_min])
        renamed_bands = bands.select(
            ["3 ETc in mm/day", "3 ETc in mm/day_1"],
            ["4 Irrigation needs in m3/day", "5 Irrigation time in min/day"]
        )

        return image.addBands(renamed_bands)

    return etc_images.map(calculate_irrigation)


def format_data_for_api(image_collection, clip_geo, field_id, current_date):
    """Convert Earth Engine image collection to pandas DataFrame format for API."""
    def extract_data(img):
        doy = img.get("DOY")
        date_string = ee.Number(doy).format()
        model_run = ee.String(date_string).cat("_SAT").cat(str(current_date))

        etref = img.reduceRegion(reducer=ee.Reducer.mean(
        ), geometry=clip_geo, scale=30).get("1 ETref in mm/day")
        etc = img.reduceRegion(reducer=ee.Reducer.mean(
        ), geometry=clip_geo, scale=30).get("3 ETc in mm/day")
        irr_vol = img.reduceRegion(reducer=ee.Reducer.mean(
        ), geometry=clip_geo, scale=30).get("4 Irrigation needs in m3/day")
        irr_time = img.reduceRegion(reducer=ee.Reducer.mean(
        ), geometry=clip_geo, scale=30).get("5 Irrigation time in min/day")

        return (
            img.set("Etref", etref)
            .set("Etc", etc)
            .set("Irrvol", irr_vol)
            .set("Irrtime", irr_time)
            .set("Date", doy)
            .set("ModelRun", model_run)
        )

    processed_imgs = image_collection.map(extract_data)

    nested_list = (
        processed_imgs.reduceColumns(
            ee.Reducer.toList(5), ["Date", "ModelRun",
                                   "Etc", "Irrvol", "Irrtime"]
        )
        .values()
        .get(0)
    )

    # Convert to pandas DataFrame
    df = pd.DataFrame(
        nested_list.getInfo(),
        columns=["date", "modelRun", "evaporation",
                 "advisedWaterVolume", "advisedIrrigationTime"]
    )

    # Format dates
    current_year = dt.datetime.now().year
    df["date"] = pd.to_datetime(
        current_year * 1000 + df["date"], format="%Y%j")
    df["modelRun"] = df["date"].dt.strftime(
        "%Y%m%d") + "_SAT" + str(current_date)
    df["dateDisplay"] = df["date"].dt.strftime(
        "%d-%m").map(lambda x: str(x)[-5:])
    df["date"] = pd.to_datetime(
        df["date"], format="%d-%m-%Y").dt.strftime("%Y-%m-%d")
    df["field"] = field_id

    # Select and order columns
    return df[[
        "date", "dateDisplay", "modelRun", "evaporation",
        "advisedWaterVolume", "advisedIrrigationTime", "field"
    ]]


def post_data_to_api(data, url, headers):
    """Post data to API endpoint."""
    data_json = data.to_json(orient="records")
    # response = requests.post(url, headers=headers, data=data_json)
    return #response.json()


def main():

    """
    Main execution function.
    """
    # Initialize Earth Engine
    initialize_gee()

    # API endpoints and authentication
    api_config = {
        "field_url": "https://sosia.tahmo.org/api/fields/",
        "schedule_url": "https://sosia.tahmo.org/api/seasonal_schedule/",
        "hindcast_url": "https://sosia.tahmo.org/api/hindcast/",
        "tahmo_api_auth": "Basic ZnV0dXJld2F0ZXI6R2MzYVdMN3kyckRkR2Y3RQ==",
        "sosia_api_auth": "Basic bC52ZXJzY2h1cmVuQGZ1dHVyZXdhdGVyLm5sOnA3XlE1OTdNNmx3Wg=="
    }

    # Headers for SOSIA API
    sosia_headers = {
        "Authorization": api_config["sosia_api_auth"]
    }

    # Headers for TAHMO API
    tahmo_headers = {
        "Authorization": api_config["tahmo_api_auth"]
    }

    # Headers for POST requests
    post_headers = {
        "Content-type": "application/json",
        "Authorization": api_config["sosia_api_auth"]
    }

    # Load Earth Engine datasets
    datasets = load_datasets()

    # Get farmer field data
    farmer_data = fetch_farmer_data(api_config["field_url"], sosia_headers)

    # Constants
    millis_per_day = 24 * 60 * 60 * 1000
    current_date = date.today()
    today_ee = ee.Date(str(current_date))

    # Process each farmer field
    for farmer in farmer_data:
        # Extract basic field information
        field_id = farmer["id"]
        field_name = farmer["name"]
        print(f"Processing field: {field_id} - {field_name}")

        # Extract coordinates
        latitude = farmer["latitude"]
        longitude = farmer["longitude"]
        coords_float = (float(longitude), float(latitude))
        point = create_point_geometry(*coords_float)
        point_string = f"{longitude},{latitude}"

        # Create buffer for spatial analysis
        buffer_geom = point.buffer(30)
        clipped_image = clip_image_to_buffer(buffer_geom)

        # Extract field parameters
        drip_lines = farmer["numberOfDriplines"]
        drip_length = farmer["lengthOfDriplines"]
        bed_width = farmer["bedWidth"]
        emitter_spacing = farmer["emitterSpacing"]
        flow_rate = farmer["emitterFlowRate"]
        init_time = farmer["initialisationTime"]

        # Calculate irrigation system parameters
        irr_params = calculate_irrigation_params(
            drip_lines, drip_length, bed_width, emitter_spacing, flow_rate, init_time
        )

        # Extract crop information
        crop_type = farmer["cropSpecific"]["cropType"]
        planting_date = ee.Date(farmer["cropSpecific"]["plantingDate"])
        harvest_date = ee.Date(
            farmer["cropSpecific"]["lastExpectedHarvestingDate"])

        # Get crop-specific parameters
        crop_params = get_crop_parameters(crop_type)

        # Calculate growth stages
        stages = calculate_growth_stages(
            planting_date, harvest_date,
            crop_params["dev_days"], crop_params["mid_days"]
        )

        # Process WAPOR data
        wapor_filtered = datasets["wapor_ret"].filterDate(
            ee.Date("2010-01-01"), ee.Date("2022-12-31")
        ).map(clipped_image)

        wapor_correct_func = wapor_correction()
        wapor_corrected = wapor_filtered.map(wapor_correct_func)

        wapor_day_avg = wapor_daily_average(wapor_corrected)

        # Define day sequences for each stage
        days_dev = ee.List.sequence(
            stages["dev_start"].millis(),
            stages["dev_end"].millis(),
            millis_per_day
        ).map(lambda ms: ee.Number.parse(ee.Date(ms).format("DDD")))

        days_mid = ee.List.sequence(
            stages["mid_start"].millis(),
            stages["mid_end"].millis(),
            millis_per_day
        ).map(lambda ms: ee.Number.parse(ee.Date(ms).format("DDD")))

        days_end = ee.List.sequence(
            stages["end_start"].millis(),
            stages["end_end"].millis(),
            millis_per_day
        ).map(lambda ms: ee.Number.parse(ee.Date(ms).format("DDD")))

        # Create collections for each stage
        wapor_collections = {
            "dev": ee.ImageCollection(days_dev.map(wapor_day_avg)),
            "mid": ee.ImageCollection(days_mid.map(wapor_day_avg)),
            "end": ee.ImageCollection(days_end.map(wapor_day_avg))
        }

        # Calculate historical ET and irrigation needs
        historical_et = calculate_historical_et(
            stages, wapor_day_avg, crop_params, planting_date, millis_per_day
        )

        # Calculate irrigation needs
        total_irrigation = calculate_irrigation_needs(
            historical_et["et_stages"], irr_params)

        # Join reference ET with crop ET and irrigation needs
        join_filter = ee.Filter.equals(leftField="Date", rightField="Date")
        simple_join = ee.Join.inner()

        inner_join = ee.ImageCollection(
            simple_join.apply(
                historical_et["et_ref"], total_irrigation, join_filter)
        )

        historical_schedule = inner_join.map(
            lambda feature: ee.Image.cat(feature.get(
                "primary"), feature.get("secondary"))
        )

        # Format historical data for API
        ################## Bug here ##################
        historical_df = format_data_for_api(
            historical_schedule, buffer_geom, field_id, current_date)

        # Post historical data to API
        historical_response = post_data_to_api(
            historical_df, api_config["schedule_url"], post_headers
        )
        print(f"Historical data posted for field {field_id}")

        # --- Hindcast (recent data) processing ---

        # Try to get TAHMO data
        tahmo_data = process_tahmo_data(point_string, tahmo_headers)

        # Calculate hindcast ET
        hindcast_et_ref = calculate_hindcast_et(
            today_ee, datasets, clipped_image, wapor_day_avg, tahmo_data
        )

        # Create Kc images for hindcast period
        kc_images = create_kc_images(stages, wapor_collections, crop_params)

        # Calculate irrigation for hindcast period
        join_by_doy = tahmo_data is not None
        hindcast_irrigation = calculate_hindcast_irrigation(
            hindcast_et_ref, kc_images, irr_params, not join_by_doy
        )

        # Format hindcast data for API
        hindcast_df = format_data_for_api(
            hindcast_irrigation, buffer_geom, field_id, current_date)
        hindcast_df["evaporation"] = hindcast_df["evaporation"].round(2)

        # Post hindcast data to API
        hindcast_response = post_data_to_api(
            hindcast_df, api_config["hindcast_url"], post_headers
        )
        print(f"Hindcast data posted for field {field_id}")


# if __name__ == "__main__":
#     main()

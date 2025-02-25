import ee
import geemap
import geemap.chart as chart
from ipyleaflet import *
from ipywidgets import *
import ipywidgets as widgets
from ipywidgets import IntSlider, Text, HTML
import json

import os
from branca.element import *
from IPython.display import Image, display

import numpy as np
import pandas as pd
import matplotlib as plt

import datetime as datetime 
from datetime import date

import requests

ee.Initialize()

#This part loads in the farmer data from the API from the dashboard. 
#We need add data of the farmer to the specific points here. 

url_api = "https://sosia.tahmo.org/api/fields/"

payload_api={}
headers_api={
    'Authorization': 'Basic bC52ZXJzY2h1cmVuQGZ1dHVyZXdhdGVyLm5sOnA3XlE1OTdNNmx3Wg=='
}

response_api = requests.request("GET", url_api, headers=headers_api, data=payload_api)
print(response_api.json())

farmer = response_api.json()[0]
for farmer in response_api.json():                                     
    #info from json file                                               
    name = farmer['id']
    print(name)
    fieldname = farmer['name']
    Lat = farmer['latitude']
    Long = farmer['longitude']
    latcoords = float(Lat)
    longcoords = float(Long)
    point = ee.Geometry.Point(longcoords, latcoords)
    pointstring = str(Long) + ',' + str(Lat)
    drip_nr = farmer['numberOfDriplines']
    dripl = farmer['lengthOfDriplines']
    wd = farmer['bedWidth']
    emitterspacing = farmer['emitterSpacing']
    flowrate = farmer['emitterFlowRate']
    intime = farmer['initialisationTime']
    flow = (dripl / emitterspacing) * flowrate
    s = dripl * (wd/drip_nr)
    LR = 1.1  #loss rate of drip irrigation system
    AE = 0.9  #application efficiency

    crops = farmer['cropSpecific']['cropType']
    pldt = ee.Date(farmer['cropSpecific']['plantingDate'])
    hvdt = ee.Date(farmer['cropSpecific']['lastExpectedHarvestingDate'])
                                                                        
    #This was the farmer specific part, I think we can just have a dataframe for the specific points. 

    if crops == 'Habanero Peppers':                     #this part is crop specific, this is still needed in the magda file
        kc1 = 0.6                                              # but  here also the ndvi algorithm needs to be added
        kc2 = 1.05
        kc3 = 1.05

        #season days
        season_devstart = pldt.advance(0, "day")
        season_dev = pldt.advance(65, "day")
        season_midstart = pldt.advance(66, "day") 
        season_mid = season_dev.advance(40, "day")
        season_endstart = season_dev.advance(41, 'day') 
        season_end = hvdt.advance(0, "day") 

    if crops == 'French Beans':
        kc1 = 0.5
        kc2 = 1.05
        kc3 = 1.05

        #season days
        season_devstart = pldt.advance(0, "day")
        season_dev = pldt.advance(65, "day")
        season_midstart = pldt.advance(66, "day") 
        season_mid = season_dev.advance(40, "day")
        season_endstart = season_dev.advance(41, 'day') 
        season_end = hvdt.advance(0, "day") 

    if crops == 'Lettuce':
        kc1 = 0.5
        kc2 = 1.05
        kc3 = 1.05

        #season days
        season_devstart = pldt.advance(0, "day")
        season_dev = pldt.advance(85, "day")
        season_midstart = pldt.advance(86, "day") 
        season_mid = season_dev.advance(40, "day")
        season_endstart = season_dev.advance(41, 'day') 
        season_end = hvdt.advance(0, "day") 

    if crops == "Brassica":
        kc1 = 0.7
        kc2 = 1.05
        kc3 = 1.05

        #season days
        season_devstart = pldt.advance(0, "day")
        season_dev = pldt.advance(65, "day")
        season_midstart = pldt.advance(66, "day") 
        season_mid = season_dev.advance(40, "day")
        season_endstart = season_dev.advance(41, 'day')
        season_end = hvdt.advance(0, "day") 

    if crops == "Cucumbers":
        kc1 = 0.60
        kc2 = 1.00
        kc3 = 1.00

        #season days
        season_devstart = pldt.advance(0, "day")
        season_dev = pldt.advance(60, "day")
        season_midstart = pldt.advance(61, "day") 
        season_mid = season_dev.advance(50, "day")
        season_endstart = season_dev.advance(51, 'day') 
        season_end = hvdt.advance(0, "day")

    if crops == "Okra":
        kc1 = 0.30
        kc2 = 1.00
        kc3 = 0.90

        #season days
        season_devstart = pldt.advance(0, "day")
        season_dev = pldt.advance(41, "day")
        season_midstart = pldt.advance(42, "day") 
        season_mid = season_dev.advance(25, "day")
        season_endstart = season_dev.advance(26, 'day') 
        season_end = hvdt.advance(0, "day") 


    # All Time variables and functions (a few variables can probably be removed but haven't cleaned it up here)
        
    months = ee.List.sequence(1,12)
    doys = ee.List.sequence(1,365)
    list_days = ee.List.sequence(1,9)
    startdate = ee.Date('2010-01-01') # for historical analysis
    enddate = ee.Date('2022-12-31') # for historical analysis
    day1 = date.today()
    stringday = str(day1)
    today = ee.Date(str(day1))
    todayformat = today.format('YYYY-MM-dd')
    pastdate = today.advance(-10,'day')
    pastdateformat = pastdate.format('YYYY-MM-dd')
    beginning = ee.Date('2021-12-31')
    millisInDay = 24*60*60*1000

    #load datasets from google earth engine
    #here we need to load the datasets from Martina and the SPHY statistics and the NDVI data
    wapor_RET = ee.ImageCollection("FAO/WAPOR/2/L1_RET_E")
    GPM = ee.ImageCollection("NASA/GPM_L3/IMERG_V06")
    CFSv2 = ee.ImageCollection("NOAA/CFSV2/FOR6H")
    DEM = ee.Image("NASA/NASADEM_HGT/001")

    #point is defined as long/lat coordinates 
    clip_geo = point.buffer(30)   # here 30 means average of the point + surrounding pixels 
    def Clipping(image):
        clipped = image.clip(clip_geo)
        return clipped
    
    #clip the point you want with the wapor data for the corresponding dates, the end result is the average daily ET0 based on the last 10 years
    # for the point that is selected

    wapor_filt = wapor_RET.filterDate(startdate,enddate).map(Clipping)

    def func_cri(img):
        corrected = img.select('L1_RET_E').divide(10).rename('corrected')
        return img.addBands(corrected)

    wapor_corr = wapor_filt.map(func_cri)

    def waporDayAve(doyx):
        dayImages = wapor_corr.select('corrected').filter(
            ee.Filter.calendarRange(start=doyx, field='day_of_year'))
        #dayID = ee.Number(doyx)
        meanWapor = ee.Image(dayImages.mean()).multiply(100).round().divide(100).set('DOY',doyx)
        return ee.Number(meanWapor)

    wapor_byDay = ee.ImageCollection(doys.map(waporDayAve))

    print(crops) #output test 

#Historical Crop Schedule --> We don't need an historical crop schedule, only the hindcast so you can ignore this part

# --------------------------------------------------------------------------------------------------------------------

    #ET ref for the season
    sMillis = season_devstart.millis()
    eMillis = season_end.millis()

    def func_iyz(dateMillis):
        return ee.Number.parse(ee.Date(dateMillis).format('DDD'))

    days_plant = ee.List.sequence(sMillis, eMillis, millisInDay).map(func_iyz)

    wapor_ref = ee.ImageCollection(days_plant.map(waporDayAve))
    
    def func_fks(image):
        return image.multiply(100).round().divide(100).copyProperties(image).set('DOY', image.get('DOY'))

    Etc_ref = wapor_ref.map(func_fks)

    def func_ldq(img):
        doyIm = ee.Number.parse(img.get('system:index'))
        dateIm = pldt.advance(doyIm, 'day')
        return img.set('Date', dateIm)

    total_refDatex = Etc_ref.map(func_ldq)

    total_refDate = total_refDatex.select(['corrected'], ['1 ETref in mm/day'])

    # ---------------------------------------------   development stage --------------------------------------------------------
    sMillis1 = season_devstart.millis()
    eMillis1 = season_dev.millis()

    def func_bsn(dateMillis):
        return ee.Number.parse(ee.Date(dateMillis).format('DDD'))

    days_dev = ee.List.sequence(sMillis1, eMillis1, millisInDay).map(func_bsn)

    wapor_dev = ee.ImageCollection(days_dev.map(waporDayAve))

    ### --------------- Convert WaPoR data according to the stage of the season where we are. Development ----> kc1 * ET0 (Wapor)
    def func_uds(image):
        return image.multiply(kc1).multiply(100).round().divide(100).copyProperties(image).set('DOY', image.get('DOY'))
                                                                                            
    Etc_dev = wapor_dev.map(func_uds)

    def func_jeu(img):
        doyIm = ee.Number.parse(img.get('system:index'))
        dateIm = pldt.advance(doyIm, 'day')
        return img.set('Date', dateIm)

    total_devDate = Etc_dev.map(func_jeu)



    # -------------------------------------------  mid stage ------------------------------------------------------------------
    sMillis2 = season_midstart.millis()
    eMillis2 = season_mid.millis()

    def func_jlx(dateMillis):
        return ee.Number.parse(ee.Date(dateMillis).format('DDD'))

    days_mid = ee.List.sequence(sMillis2, eMillis2, millisInDay).map(func_jlx)

    wapor_mid = ee.ImageCollection(days_mid.map(waporDayAve))

    ### --------------- Convert WaPoR data according to the stage of the season where we are. Mid ----> kc2* ET0 (Wapor)
    def func_tma(image):
        return image.multiply(kc2).multiply(100).round().divide(100).copyProperties(image).set('DOY', image.get('DOY'))

    Etc_mid = wapor_mid.map(func_tma)

    def func_rgc(img):
        doyIm = ee.Number.parse(img.get('system:index'))
        dateIm = season_midstart.advance(doyIm, 'day')
        return img.set('Date', dateIm)

    total_midDate = Etc_mid.map(func_rgc)

    # -------------------------------------------- end stage ----------------------------------------------------------------
    sMillis3 = season_endstart.millis()
    eMillis3 = season_end.millis()

    def func_sud(dateMillis):
        return ee.Number.parse(ee.Date(dateMillis).format('DDD'))

    days_end = ee.List.sequence(sMillis3, eMillis3 , millisInDay).map(func_sud)

    wapor_end = ee.ImageCollection(days_end.map(waporDayAve))

    ### --------------- Convert WaPoR data according to the stage of the season where we are. End ----> kc3* ET0 (Wapor)
    def func_kll(image):
        return image.multiply(kc3).multiply(100).round().divide(100).copyProperties(image).set('DOY', image.get('DOY'))

    Etc_end = wapor_end.map(func_kll)

    def func_rdn(img):
        doyIm = ee.Number.parse(img.get('system:index'))
        dateIm = season_endstart.advance(doyIm, 'day')
        return img.set('Date', dateIm)

    total_endDate = Etc_end.map(func_rdn)

    # ------------------------------------ Merge Etc computed from the three stages: dev, mid and end ------------------------------------------------------
    #create results
    Etc_merge = total_devDate.merge(total_midDate)
    Etc = Etc_merge.merge(total_endDate);                # mm per dag
    ETc = Etc.select(['corrected'],['2 ETc in mm/day'])


    # ----------------------------------- Computing the irrigation needs ------------------------------------------------------------------------------------
    def func_wpd(image):
        irr_m3 = image.divide(1000).multiply(s).multiply(LR).divide(AE).multiply(100).round().divide(100).copyProperties(image).set('DOY', image.get('DOY'))
        irr_min = image.divide(1000).multiply(s).multiply(LR).divide(AE).multiply(60000).divide(flow).round().add(intime).copyProperties(image).set('DOY', image.get('DOY'))
        newbands = ee.Image([irr_m3])
        adjustbands = newbands.select(['2 ETc in mm/day'], ['3 Irrigation needs in m3/day'])
        newbands1 = ee.Image([irr_min])
        adjustbands1 = newbands1.select(['2 ETc in mm/day'], ['4 Irrigation time in min/day'])
        return image.addBands(adjustbands).addBands(adjustbands1)

    # ------------------ one image per day in the season including 3 bands: Etc, irr (m3/day) and irr (min/day) ----------------------------------------------
    total = ETc.map(func_wpd)



    ## ------------------------------------------- formatting for the final historical crop table -------------------------------------------------------------- 
    # FILTERS 
    #Combine two image collections by date
    joinFilter = ee.Filter.equals(leftField='Date', rightField='Date')

    #Combine two image collections by DOY
    doyFilter = ee.Filter.equals(leftField='DOY', rightField='DOY')

    # Create the join.
    simpleJoin = ee.Join.inner()

    # Inner join
    innerJoin = ee.ImageCollection(simpleJoin.apply(total_refDate, total, joinFilter))

    def func_nlb(feature):
        return ee.Image.cat(feature.get('primary'), feature.get('secondary'))



    ## One image per day of the season, including Etref (WaPOR), Etc (WaPOR*kcx), irr (m3) and irr (minutes) from the formulas  
    
    Historical_crop_schedule = innerJoin.map(func_nlb)  #joined is csv of 125 days of the season 

    ## ---------------- Translate from EE images to csv tables (by averaging around the buffer created at the begining) --------------------------
    def poi_mean(img):
            #Kc = img.reduceRegion(reducer=ee.Reducer.mean(), geometry=clip_geo, scale=30).get('0 Daily Kc')
            DOY = img.get('DOY')
            datestring = ee.Number(DOY).format()
            ModelRun1 = ee.String(datestring).cat("_SAT")
            ModelRun = ee.String(ModelRun1).cat(stringday)#d2)
            Etref = img.reduceRegion(reducer=ee.Reducer.mean(), geometry=clip_geo, scale=30).get('1 ETref in mm/day')
            Etc = img.reduceRegion(reducer=ee.Reducer.mean(), geometry=clip_geo, scale=30).get('2 ETc in mm/day')
            irrvol = img.reduceRegion(reducer=ee.Reducer.mean(), geometry=clip_geo, scale=30).get('3 Irrigation needs in m3/day')
            irrtime = img.reduceRegion(reducer=ee.Reducer.mean(), geometry=clip_geo, scale=30).get('4 Irrigation time in min/day')
            return img.set('Etref', Etref).set('Etc', Etc).set('Irrvol', irrvol).set('Irrtime', irrtime).set('Date', DOY).set('ModelRun', ModelRun)

    poi_reduced_imgs = Historical_crop_schedule.map(poi_mean)

    nested_list = poi_reduced_imgs.reduceColumns(ee.Reducer.toList(5), ['Date', 'ModelRun', 'Etc', 'Irrvol', 'Irrtime']).values().get(0)

    df = pd.DataFrame(nested_list.getInfo(), columns=["date", "modelRun", "evaporation", "advisedWaterVolume", "advisedIrrigationTime"])
    df["date"] = pd.to_datetime(2024 * 1000 + df["date"], format='%Y%j')
    df["modelRun"] = df["date"].dt.strftime("%Y%m%d") + '_SAT' + stringday
    df["dateDisplay"] = df["date"].dt.strftime('%d-%m').map(lambda x: str(x)[-5:])
    df['date'] = pd.to_datetime(df['date'], format='%d-%m-%Y').dt.strftime('%Y-%m-%d')
    df["field"] =  farmer['id'] #farmer['name']
    df = df[["date", "dateDisplay", "modelRun", "evaporation", "advisedWaterVolume", "advisedIrrigationTime", "field"]]
    data_json = df.to_json(orient='records')
    print(data_json)
    
    url_post = "https://sosia.tahmo.org/api/seasonal_schedule/"
    headers_post = {'Content-type': 'application/json',
            'Authorization': 'Basic bC52ZXJzY2h1cmVuQGZ1dHVyZXdhdGVyLm5sOnA3XlE1OTdNNmx3Wg=='   
            }

    response_post = requests.post(url_post, headers=headers_post, data=data_json)
    print(response_post.json())
    
    ### ------------------------------------------------------------------------------------------------------------------------

# -------------------- END OF THE HISTORICAL CROP SCHEDULE ---------------------------------------------------------------------


    #From here we are going to make the hindcasted ET0, based on both hindcasted data from gee and tahmo 

    ### TAHMO ### # Here the tahmo data is called

    #load in the tahmo wetness reports (multiple variable, we only use ETact)
    url_tahmo = "https://smartirrigation.tahmo.org/sosia/nearestwetnessreport/" + pointstring 

    payload_tahmo={}
    headers_tahmo={
        'Authorization': 'Basic ZnV0dXJld2F0ZXI6R2MzYVdMN3kyckRkR2Y3RQ=='
    }

    response_tahmo = requests.request("GET", url_tahmo, headers=headers_tahmo, data=payload_tahmo)

    result = response_tahmo.json() #.get('result')

    #Either 0, 1, 2, or 3 stations will be close to a specific point. 
    #if we have a tahmo output (so 1,2 or 3 close tahmo stations), the if result loop will generate results. otherwise the else part is chosen.
   
    if result:
        results = pd.json_normalize(result, "results", "distance")  #SOLUTION TO DISTANCE PROBLEM :) 
        results2 = results.drop(columns=["_id", "Wetness", "FC", "RAM",  "WP", "Eact"])
        results2["Day"] = pd.to_datetime(results2["Time"]).dt.strftime('%m/%d/%y')
        results2.drop(columns=["Time"])

        days = results2.pivot_table("Eref", "Day", "Station")
        distance = results.pivot_table("distance", "Station")
  
        # IDW
        if distance.count()[0] == 1:
            print('one station retrieved, value of nearest station retrieved')
            point_value = days
            
        elif distance.count()[0] == 2: 
            print('two stations retrieved')
            point_value = days
            d1 = distance['distance'][0] + 28
            d2 = distance['distance'][1]
            v1 = days.iloc[:,0]
            v2 = days.iloc[:,1]
            point_value = ((v1/d1) + (v2/d2)) / ((1/d1)+(1/d2))
            print(point_value)

        else:
            print('three stations retrieved')
            d1 = distance['distance'][0]
            d2 = distance['distance'][1]
            d3 = distance['distance'][2]
            v1 = days.iloc[:,0]
            v2 = days.iloc[:,1]
            v3 = days.iloc[:,2]
            statdist1 = v1/d1
            statdist2 = v2/d2
            statdist3 = v3/d3
            point_value = ((v1/d1) + (v2/d2) + (v3/d3)) / ((1/d1)+(1/d2)+(1/d3))
            print(point_value)

        ### IF POINT VALUE = smaller than 1: get CFSv2 waarde

        TAHMO = days.iloc[:, 0]
        TAHMO_Closest = TAHMO.reset_index(drop=True)
        tahmolist = TAHMO_Closest.tolist()


        # going one week backwards to compute ETref
        # tahmo data converted to imagecollection
        weekago = today.advance(-6, 'day')
        todaymillis = today.millis()
        weekagomillis = weekago.millis()

        def func_iyz(dateMillis):
            return ee.Number.parse(ee.Date(dateMillis).format('DDD'))

        days_hind = ee.List.sequence(weekagomillis, todaymillis, millisInDay).map(func_iyz)

        tahmo_hind = ee.ImageCollection(days_hind.map(waporDayAve))

        max_elements = 1000
        zippedList = tahmo_hind.toList(max_elements).zip(tahmolist)

        def listmaker(list):
            list = ee.List(list)
            img = ee.Image(list.get(0))
            scale = list.getNumber(1)
            scaledImage = img.multiply(0).add(scale)
            scaledImage = scaledImage.double().rename('1 ETref in mm/day')
            return scaledImage.copyProperties(img, img.propertyNames())

        tahmo_list = zippedList.map(listmaker)

        TAHMO_Eref = ee.ImageCollection.fromImages(tahmo_list)

    ### ---------------------------------- This is the important part, when no stations are available -------------------------------------------------------------------------
    
    else:    # If there is no tahmo data available for a point, this loop will follow. this is the loop we want, except that here we use penman 
             # Instead of Hargreaves. 
        print('No  tahmo available')
        # Hindcast CFSv2
        weekago = today.advance(-9, 'day')   #-9 and -2 are chosen, because the hindcast only gives information of 2 days ago, so we are using a moving window loop here
        weekago2 = ee.Date(weekago.format('YYYY-MM-dd'))
        weekagoformat = ee.Number.parse(weekago.format('DDD'))
        yesterdate = today.advance(-2, 'day')
        yesterdate2 = ee.Date(yesterdate.format('YYYY-MM-dd'))
        yesterdateformat = ee.Number.parse(yesterdate.format('DDD'))
        doys7 = ee.List.sequence(weekagoformat, yesterdateformat, 1)

        #Filter data
        CFSR_filt_first  = CFSv2.filterDate(weekago2,yesterdate2)
        CFSR_filt = ee.ImageCollection(CFSR_filt_first.map(Clipping))

        # Compute radiation from CSFv2 
        def calc_rnet(image):
            dsw = image.select('Downward_Short-Wave_Radiation_Flux_surface_6_Hour_Average')
            dlw = image.select('Downward_Long-Wave_Radp_Flux_surface_6_Hour_Average')
            ulw = image.select('Upward_Long-Wave_Radp_Flux_surface_6_Hour_Average')
            usw = dsw.multiply(0.23).rename('Upward_Short-Wave_Radiation_Flux_surface_6_Hour_Average')
            rnet = dsw.subtract(usw).add(dlw.subtract(ulw)).rename('Net_Radiation_6_Hour_Average')
            u = image.select('u-component_of_wind_height_above_ground')
            v = image.select('v-component_of_wind_height_above_ground')
            Uwind = u.pow(2)
            Vwind = v.pow(2)
            Wind_tot = (Uwind.add(Vwind)).sqrt().rename('Wind_component')
            rnetband = image.addBands(rnet).addBands(Wind_tot)
            return rnetband

        #input for daily mean, min and max
        CFSR_tot = CFSR_filt.map(calc_rnet);
        CFSR_Temp = CFSR_filt.select('Temperature_height_above_ground')

        #Daily mean of images 
        numberOfDays = yesterdate2.difference(weekago2, 'days')

        def daily_mean(dayOffset): 
            start = weekago2.advance(dayOffset, 'days')
            end = start.advance(1, 'days')
            dateIm = start.format('YYYY-MM-dd')
            doyIm = ee.Number.parse(start.format('DDD'))
            return CFSR_tot.filterDate(start, end).mean().set('Date', dateIm).set('DOY', doyIm)   #hier ergens bij deze dingen gaat het mis; string is 085, je wil 85

        daily_mean = ee.ImageCollection(ee.List.sequence(0, numberOfDays.subtract(1)).map(daily_mean))

        #daily max of temp images

        def CFSR_Temp_Max(dayOffset): 
            start = weekago2.advance(dayOffset, 'days')
            end = start.advance(1, 'days')
            return CFSR_Temp.filterDate(start, end).max().set('Date', start.format('YYYY-MM-dd'))#.set('DOY', start.format('DDD'))

        CFSR_Temp_Max = ee.ImageCollection(ee.List.sequence(0, numberOfDays.subtract(1)).map(CFSR_Temp_Max))

        Temp_max = CFSR_Temp_Max.map(lambda image: image.select('Temperature_height_above_ground').rename('Max Temperature'))

        #daily min of temp images

        def CFSR_Temp_Min(dayOffset): 
            start = weekago2.advance(dayOffset, 'days')
            end = start.advance(1, 'days')
            return CFSR_Temp.filterDate(start, end).min().set('Date', start.format('YYYY-MM-dd'))#.set('DOY', start.format('DDD'))

        CFSR_Temp_Min = ee.ImageCollection(ee.List.sequence(0, numberOfDays.subtract(1)).map(CFSR_Temp_Min))

        Temp_min = CFSR_Temp_Min.map(lambda image: image.select('Temperature_height_above_ground').rename('Min Temperature'))

        # join daily mean, max en min

        innerJoin_min = ee.ImageCollection(simpleJoin.apply(Temp_max, Temp_min, joinFilter))

        tempmaxmin = innerJoin_min.map(lambda feature: ee.Image.cat(feature.get('primary'), feature.get('secondary')))

        innerJoin_max = ee.ImageCollection(simpleJoin.apply(daily_mean, tempmaxmin, joinFilter))

        CFSR_total = innerJoin_max.map(lambda feature: ee.Image.cat(feature.get('primary'), feature.get('secondary')))


        # Compute ET from Hargreaves 
        def allbands(image):
            Tmin_d_k = image.select('Min Temperature')
            Tmax_d_k = image.select('Max Temperature')
            Wind_d_10m = image.select('Wind_component')
            SHum_d = image.select('Specific_humidity_height_above_ground')
            Press_d_pa = image.select('Pressure_surface')
            Rnet_w = image.select('Net_Radiation_6_Hour_Average')
            #correct Units
            Tmin_d = Tmin_d_k.subtract(273.15).rename('Tmin')
            Tmax_d = Tmax_d_k.subtract(273.15).rename('Tmax')
            Press_d = Press_d_pa.divide(1000).rename('Atm pressure')
            Wind_d = Wind_d_10m.multiply(0.75).rename('Wind')
            Rnet_d = Rnet_w.multiply(0.0864).rename('Rnet')
            exps = Press_d.multiply(0).add(2.71828)
            #calculate delta
            tmean_1 = Tmin_d.add(Tmax_d)
            tmean = tmean_1.divide(2)
            delta_1 = tmean.multiply(17.27).divide(tmean.add(237.3))
            delta_2 = exps.pow(delta_1).multiply(0.6108)
            delta_3 = tmean.add(237.3)
            delta_4 = delta_3.pow(2)
            delta_d = delta_2.multiply(4098).divide(delta_4)
            #calculate es
            e0_min_1 = Tmin_d.multiply(17.27).divide(Tmin_d.add(237.3))
            e0_min = exps.pow(e0_min_1).multiply(0.6108)
            e0_max_1 = Tmax_d.multiply(17.27).divide(Tmax_d.add(237.3))
            e0_max = exps.pow(e0_max_1).multiply(0.6108)
            es_d_1 = e0_max.add(e0_min)
            es_d = es_d_1.divide(2)
            #calculate ea
            rhmax = SHum_d.multiply(Press_d).multiply(1.6077717).divide(e0_min).rename('RHmax')
            rhmin = SHum_d.multiply(Press_d).multiply(1.6077717).divide(e0_max).rename('RHmin')
            ea_1 = e0_min.multiply(rhmax)
            ea_2 = e0_max.multiply(rhmin)
            ea_3 = ea_1.add(ea_2)
            ea_d = ea_3.divide(2)
            #-------------    calculate ET0 --------------------------------------------------
            part_1 = delta_d.multiply(Rnet_d).multiply(0.408)
            psy_c = Press_d.multiply(0.001).divide(1.53634)
            part_2b = tmean.add(273)
            part_2 = psy_c.multiply(900).divide(part_2b)
            part_3 = Wind_d.multiply(es_d.subtract(ea_d))
            part_4a = Wind_d.multiply(0.34).add(1)
            part_4 = psy_c.multiply(part_4a).add(delta_d)
            ET0_1 = part_1.divide(part_4)
            ET0_2 = part_2.multiply(part_3).divide(part_4)
            ET0 = ET0_1.add(ET0_2).rename('1 ETref in mm/day') 
            #newbands = ee.Image([ET0]);#, irr_min]);
            return image.addBands(ET0)

        ## 7 images with all bands above + 1 band EtRef mm per day (Hargreaves) ---------------- 
        ET = CFSR_total.map(allbands)
        ET0 = ET.select('1 ETref in mm/day')
        
    # --------- Go from ET0/ETref to ETc -----> multiplying by kc  -----------------------------
    
    #create kc values to multiply with the ETref. The Kc values are dependent by day. For our script, we need to add and 
    # IF ELSE statement with these FAO56 values
    
    ##### WAPOR here again? ------------------------
    def kcdev(image): 

        return image.multiply(0).add(kc1).copyProperties(image).set('DOY', image.get('DOY'))
        
    Kc_dev = wapor_dev.map(kcdev)

    def kcdevdate(image):
        doyIm = ee.Number.parse(image.get('system:index'))
        dateIm = (pldt.advance(doyIm, 'day')).format('YYYY-MM-dd')
        return image.set('Date', dateIm)
    
    Kc_devDate = Kc_dev.map(kcdevdate)

    def kcmid(image):
        return image.multiply(0).add(kc2).copyProperties(image).set('DOY', image.get('DOY'))
    Kc_mid = wapor_mid.map(kcmid)

    def kcmiddate(image):
        doyIm = ee.Number.parse(image.get('system:index'))
        dateIm = (season_midstart.advance(doyIm,'day')).format('YYYY-MM-dd')
        return image.set('Date', dateIm)
    Kc_midDate = Kc_mid.map(kcmiddate)

    def kcend(image):
        return image.multiply(0).add(kc3).copyProperties(image).set('DOY', image.get('DOY'))
    Kc_end = wapor_end.map(kcend)

    def kcenddate(image):
        doyIm = ee.Number.parse(image.get('system:index'))
        dateIm = (season_endstart.advance(doyIm, 'day')).format('YYYY-MM-dd')
        return image.set('Date', dateIm)
    Kc_endDate = Kc_end.map(kcenddate)

    Kc_merge = Kc_devDate.merge(Kc_midDate)
    Kc = Kc_merge.merge(Kc_endDate)
    Kc_season = Kc.select(['corrected'], ['2 Daily Kc'])


    #---- HERE all KC values are put in a dataframe and now we can multiply the dataframe with KC values with the ETref values 
    
    if result:
        KcJoin = ee.ImageCollection(simpleJoin.apply(TAHMO_Eref, Kc_season, doyFilter))
        print('tahmodata used')
    else:
        KcJoin = ee.ImageCollection(simpleJoin.apply(ET0, Kc_season, joinFilter))
        print('CFSv2 data used')

    def func_kc(feature):
        return ee.Image.cat(feature.get('primary'), feature.get('secondary'))

    joined_kc_cfs = KcJoin.map(func_kc)  #joined is csv of 5 days hindcast

    def etc7(image):
        ETref = image.select('1 ETref in mm/day')
        Kc = image.select('2 Daily Kc')
        ETc = ETref.multiply(Kc).rename('3 ETc in mm/day')
        return image.addBands(ETc)

    ETc_7 = joined_kc_cfs.map(etc7)

    # ETC7 = ETc_7.select(['1 ETref in mm/day'], ['2 ETc in mm/day'])

    def tot7(image):
        irr_m3 = image.divide(1000).multiply(s).multiply(LR).divide(AE).multiply(100).round().divide(100).copyProperties(image)
        irr_min = image.divide(1000).multiply(s).multiply(LR).divide(AE).multiply(60000).divide(flow).round().add(intime).copyProperties(image)
        newbands = ee.Image([irr_m3, irr_min])
        adjustbands = newbands.select(['3 ETc in mm/day', '3 ETc in mm/day_1'], ['4 Irrigation needs in m3/day', '5 Irrigation time in min/day'])
        return image.addBands(adjustbands)
    
    ### Irrigation for the hindcast period ------ (same structure as at the begining)
    Hindcasted_crop_schedule = ETc_7.map(tot7) # deze moet DOY bevatten

    # Spatial average 
    def poi2_mean(img):
        DOY = img.get('DOY')
        datestring = ee.Number(DOY).format()
        ModelRun1 = ee.String(datestring).cat("_SAT")
        ModelRun = ee.String(ModelRun1).cat(stringday)#d2)
        Etref = img.reduceRegion(reducer=ee.Reducer.mean(), geometry=clip_geo, scale=30).get('1 ETref in mm/day')
        #Kc = img.reduceRegion(reducer=ee.Reducer.mean(), geometry=clip_geo, scale=30).get('1 ETref in mm/day_1')
        Etc = img.reduceRegion(reducer=ee.Reducer.mean(), geometry=clip_geo, scale=30).get('3 ETc in mm/day')
        irrvol = img.reduceRegion(reducer=ee.Reducer.mean(), geometry=clip_geo, scale=30).get('4 Irrigation needs in m3/day')
        irrtime = img.reduceRegion(reducer=ee.Reducer.mean(), geometry=clip_geo, scale=30).get('5 Irrigation time in min/day')
        return img.set('Etref', Etref).set('Etc', Etc).set('Irrvol', irrvol).set('Irrtime', irrtime).set('Date', DOY).set('ModelRun', ModelRun)
    
    poi_reduced_imgs2 = Hindcasted_crop_schedule.map(poi2_mean)

    nested_list2 = poi_reduced_imgs2.reduceColumns(ee.Reducer.toList(5), ['Date', 'ModelRun', 'Etc', 'Irrvol', 'Irrtime']).values().get(0)

    df2 = pd.DataFrame(nested_list2.getInfo(), columns=["date", "modelRun", "evaporation", "advisedWaterVolume", "advisedIrrigationTime"])#
    df2['evaporation']  = df2['evaporation'].round(2)
    df2["date"] = pd.to_datetime(2024 * 1000 + df2["date"], format='%Y%j')
    df2["modelRun"] = df2["date"].dt.strftime("%Y%m%d") + '_SAT' + stringday
    df2['datedisp'] = df2["date"] + pd.DateOffset(days=0)#.map(lambda x: str(x)[-5:])
    df2["dateDisplay"] = df2["datedisp"].dt.strftime('%d-%m').map(lambda x: str(x)[-5:])
    df2["date"] = pd.to_datetime(df2['date'] + pd.DateOffset(days=0), format='%d-%m-%Y').dt.strftime('%Y-%m-%d')
    df2["field"] =  farmer['id'] #farmer['name']i\

    df2 = df2[["date", "dateDisplay", "modelRun", "evaporation", "advisedWaterVolume", "advisedIrrigationTime", "field"]]
    #print(df2)
    data_json_hind = df2.to_json(orient='records')
    print(data_json_hind)

    ### ook command om te kijken of er data aanwezig is voor het juiste timeframe! 

    url_hind = "https://sosia.tahmo.org/api/hindcast/"
    headers_hind = {'Content-type': 'application/json',
            'Authorization': 'Basic bC52ZXJzY2h1cmVuQGZ1dHVyZXdhdGVyLm5sOnA3XlE1OTdNNmx3Wg=='   
            }

    response_hind = requests.post(url_hind, headers=headers_hind, data=data_json_hind)

    response_hind.json()


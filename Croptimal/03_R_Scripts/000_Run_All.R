#### Clean working directory
rm(list = ls())

## Load necessary libraries
library(raster)
library(ncdf4)
library(rjson)
library(geojsonio)
library(jsonlite)
library(Rcurl)
library(httr)
library(RJSONIO)
library(dplyr)
library(rgdal)
library(rgeos)
#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

##############################################################################################
##################################### START OF INPUT #########################################
##############################################################################################

### Directories dropbox

data_drop = "C:\\Users\\Lisa\\Dropbox (FutureWater)\\FW_VH_RK\\01_Data\\"
results_drop = "C:\\Users\\Lisa\\Dropbox (FutureWater)\\FW_VH_RK\\04_Results\\"
#temp_dir = "C:\\Users\\Lisa\\Dropbox (FutureWater)\\FW_VH-RK\\05_Temp\\"
provinces = "C:\\Users\\Lisa\\Dropbox (FutureWater)\\FW_VH_RK\\02_GIS\\Shapefiles\\AGO_adm1.shp"
#communes = "C:\\Users\\Lisa\\Dropbox (FutureWater)\\FW_VH_RK\\02_GIS\\Shapefiles\\AGO_adm3.shp"
DEM = "C:\\Users\\Lisa\\Dropbox (FutureWater)\\FW_VH_RK\\01_Data\\Elevation\\SRTM_30m_Angola_mask.tif"
# cropping_calender = paste0(getwd(), "\\Cropping_calender.csv")
# parameters = paste0(getwd(), "\\Parameters.csv")
HHS_data_dir = "H:\\01_Data\\__TO_DROPBOX__\\Top_Subsoil\\"
#Seasonal_Forecast_dir = "D:\\Dropbox (FutureWater)\\FW_VH_RK\\01_Data\\Seasonal_Forecast\\"

### Directories local
#data_dir = "./01_Data/"
results_dir = "C:\\Users\\Lisa\\Documents\\Projects\\2019019_MavoDiami\\04_Results\\"
temp_dir = "C:\\Users\\Lisa\\Documents\\Projects\\2019019_MavoDiami\\05_Temp"
#provinces = "./02_GIS/Shapefiles/AGO_adm1.shp"
communes = "C:\\Users\\Lisa\\Documents\\Projects\\2019019_MavoDiami\\Shapefiles\\AGO_adm3.shp" #Communes_AB
#communes = "./02_GIS/Shapefiles/AGO_adm3.shp"
#DEM = "./01_Data/Elevation/SRTM_30m_Angola_mask.tif"
cropping_calender = "C:\\Users\\Lisa\\Documents\\Projects\\2019019_MavoDiami\\03_R_Scripts\\Cropping_calendar.csv"
cc_a = "C:\\Users\\Lisa\\Documents\\Projects\\2019019_MavoDiami\\03_R_Scripts\\Cropping_calendar_A.csv"
cc_b ="C:\\Users\\Lisa\\Documents\\Projects\\2019019_MavoDiami\\03_R_Scripts\\Cropping_calendar_B.csv"
parameters = "C:\\Users\\Lisa\\Documents\\Projects\\2019019_MavoDiami\\03_R_Scripts\\Parameters.csv"
#HHS_data_dir = "./01_Data/__TO_DROPBOX__/Top_Subsoil/"
Seasonal_Forecast_dir = "C:\\Users\\Lisa\\Documents\\Projects\\2019019_MavoDiami\\Seasonalforecast\\"
Scripts_dir = "C:\\Users\\Lisa\\Documents\\Projects\\2019019_MavoDiami\\03_R_Scripts" 


new_dir = "C:\\Users\\Lisa\\Documents\\Projects\\2019019_MavoDiami\\04_Results\\_LS_Results\\_Communes\\"
new_drop = "C:\\Users\\Lisa\\Dropbox (FutureWater)\\FW_VH_RK\\04_Results\\_LS_Results\\_Communes\\"

### If new crops are added, set this switch to 1 (all base data is already present like T, P, etc.)
### Make sure to choose the right cropping_calender file or add lines to the existing one.
switch = 2 #0: all, 1: new crops only, 2: seasonal forecast only 

### Resolution and projection of final rasters
res = 250 #meter
Local_proj  = "EPSG:32733"

### Temperature input
T_vars = c("tavg", "tmin", "tmax")
T_lapse_rate <- -0.0065
T_perc = c(0.75, 0.95)
T_perc_names = c("Warmer", "Much_Warmer")

### Rainfall input
P_perc = c(0.25, 0.05)
P_perc_names = c("Drier", "Much_Drier")

##############################################################################################
##################################### END OF INPUT ###########################################
##############################################################################################

list.scripts = list.files(Scripts_dir, pattern = "\\.R$", full.names = T)
list.scripts = list.scripts[-grep("000_Run_All.R", list.scripts)]

DEM_r = raster(DEM)
provinces_shp = readOGR(provinces)
provinces_names = provinces_shp$NAME_1
cropping_cal = read.csv(cropping_calender, check.names=FALSE, sep = ",")
cropping_cal_a = read.csv(cc_a, check.names=FALSE, sep = ",")
cropping_cal_b = read.csv(cc_b, check.names=FALSE, sep = ",")

params = read.csv(parameters, check.names=FALSE)

unique_combis <- c(outer("Average",T_perc_names, paste, sep = "_"), 
                   outer(P_perc_names, "Average" , paste, sep = "_"),
                   outer(P_perc_names, T_perc_names, paste, sep = "_"),
                   paste0("Average_Average"))

unique_combis = paste0("Land_Suitability_", unique_combis)

communes_shp = readOGR(communes)
communes_shp_proj = spTransform(communes_shp, Local_proj)
communes_df = data.frame("Province" = gsub(" " , "_", communes_shp$NAME_1), "Municipality" = communes_shp$NAME_2, 
                         "Commune" = communes_shp$NAME_3)#, "Zone" = communes_shp$Zone)

final_forecast_df = data.frame("Province" = NA,
                                 "Commune" = rep(communes_df$Commune, each=nrow(cropping_cal)), 
                                 "Crop" = paste0(cropping_cal$Crop, "_", month.abb[cropping_cal$Start_growing_season], "_",
                                                                                  month.abb[cropping_cal$End_growing_season]),
                                 matrix(NA, ncol = 12, nrow = nrow(cropping_cal)*nrow(communes_df), dimnames = list(NULL, month.abb[seq(1,12,1)])),
                               "Planting_Suitability_M1" = NA,
                               "Planting_Suitability_M2" = NA,
                               "Planting_Suitability_M3" = NA,
                               "Crop_Start_Month" = NA,
                               #"Planting_Suitability_M5" = NA,
                               "Average_Average" = NA, 
                               "Threshold" = NA) 
                               

for(i in 1:nrow(final_forecast_df)){
  final_forecast_df$Province[i] = communes_df$Province[final_forecast_df$Commune[i] == communes_df$Commune]
}


for(c in 1:nrow(cropping_cal)){
  crop = cropping_cal$Crop[c]
  start_month = cropping_cal$Start_growing_season[c]
  end_month = cropping_cal$End_growing_season[c]
  name_col = paste0(crop, "_", month.abb[start_month], "_", month.abb[end_month])
  communes_df[, ncol(communes_df) + 1] <- NA
  names(communes_df)[ncol(communes_df)] <- name_col
}

for(comb in unique_combis){
  assign(paste0(comb, "_Communes.csv"), communes_df)
}

#rasterOptions(overwrite=TRUE,maxmemory = 1e+10, chunksize=0.1e+10, progress="text", tmpdir=temp_dir) 
rasterOptions(overwrite=TRUE,maxmemory = 1e+10, chunksize=0.1e+10, progress="", tmpdir=temp_dir) 
start_time <- Sys.time()

## no. of script: switch = 1: 05_Calculate_Limits_Parameters.R = 11, switch = 2: 08_Seasonal_Forecast.R

switch = 2 # HARDCODED 
n = ifelse(switch == 1, 10, ifelse(switch == 2, 15, 1))

if(n == 15) {source("C:\\Users\\Lisa\\Documents\\Projects\\2019019_MavoDiami\\03_R_Scripts\\Helper_functions\\01_Download_Seasonal_forecast_WI_API.R")}

###alle provincies
for(script in n:length(list.scripts)){
  print(script)
  for(name in provinces_names){
    rasterTmpFile("clean_this_after_")
    
    province_shp = provinces_shp[provinces_shp$NAME_1 == name,]
    province_name = gsub(" ", "_", name)
    
    indir = paste0(results_dir, province_name, "\\")
  
    print(paste0("Province: ",province_name, " & script: ", basename(list.scripts[script])))
    source(list.scripts[script])
    
    temp_files = list.files(temp_dir, pattern = "\\.grd$|\\.gri$", full.names = T)
    unlink(temp_files, recursive = T)
  }
}
end_time <- Sys.time()
end_time - start_time


#### Run script voor alleen 1 provincie;

for(script in n:length(list.scripts)){
  print(script)
  name = provinces_names[2]
  rasterTmpFile("clean_this_after_")
  
  province_shp = provinces_shp[provinces_shp$NAME_1 == name,]
  province_name = gsub(" ", "_", name)
  
  indir = paste0(results_dir, province_name, "\\")
  
  print(paste0("Province: ",province_name, " & script: ", basename(list.scripts[script])))
  source(list.scripts[script])
  
  temp_files = list.files(temp_dir, pattern = "\\.grd$|\\.gri$", full.names = T)
  unlink(temp_files, recursive = T)
}
end_time <- Sys.time()
end_time - start_time
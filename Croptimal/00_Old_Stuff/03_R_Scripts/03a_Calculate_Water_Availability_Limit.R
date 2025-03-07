#This is a script to obtain the following maps:
#1) Total ETc for a specific growth period and crop type 
#2) Total Precipitation for a specific growth period and crop type
#3) water availability index for a specific growth period and crop type
#Kc values are obtained from FAO: "Crop evapotranspiration - Guidelines for computing crop water requirements - FAO Irrigation and drainage paper 56"

### Select directories abd files
input_files_ETc = list.files(paste0(results_dir, province_name, "\\Ref_ET\\Mean_Monthly"), pattern = "\\.tif$", full.names = T)
input_files_P = list.files(paste0(results_dir, province_name, "\\Rainfall"), pattern = "\\.tif$", full.names = T, recursive =T)
new_dir = paste0(results_dir,province_name, "\\_LS_Results\\")

for (f in c("Water", "ETc", "Rainfall")){
  dir.create(paste0(new_dir,f))
}


for(p in 1:(length(P_perc)+1)){
  if(p > length(P_perc)){
    per = paste0("Mean_Monthly")
    dir.create(paste0(new_dir,"Rainfall\\", per))
    dir.create(paste0(new_dir,"Water\\", per))
  }else{
    per = paste0(gsub("\\.", "" ,P_perc[p]), "perc")
    dir.create(paste0(new_dir,"Rainfall\\", per))
    dir.create(paste0(new_dir,"Water\\", per))
  }
  
  for(c in 1:nrow(cropping_cal)){
    crop_name = cropping_cal$Crop[c]
    start_month = cropping_cal$Start_growing_season[c]
    end_month = cropping_cal$End_growing_season[c]
    if(start_month > end_month){
      months = c(seq(start_month,12,1), seq(1, end_month, 1))
      kc = as.numeric(unlist(c(cropping_cal[c,4 + seq(start_month,12,1)], cropping_cal[c, 4 + seq(1, end_month, 1)])))
    }else{
      months = seq(start_month, end_month, 1) 
      kc = as.numeric(cropping_cal[c, 4 + seq(start_month, end_month, 1)])
    }
    
    ET_files = stack(input_files_ETc[grep(paste(month.abb[months],collapse="|"), input_files_ETc)])
    P_files_all = input_files_P[grep(paste(month.abb[months], collapse="|"), input_files_P)]
    P_files_per = stack(P_files_all[grep(per, P_files_all)])
    ET_sum = sum(ET_files)
    P_sum = sum(P_files_per)
    Water = P_sum/ET_sum
    
    ### Calculate the limits for water availability
    limit = as.numeric(params$Limit[params$Parameter == "Water"])
    Water_limit <- Water>as.numeric(limit)
    
    ### Write rasters
    writeRaster(ET_sum, paste0(new_dir, "ETc\\ETc_",crop_name, "_", month.abb[start_month], "-", month.abb[end_month], ".tif"), overwrite = T)
    writeRaster(P_sum, paste0(new_dir, "Rainfall\\", per, "\\P_",crop_name, "_", month.abb[start_month], "-", month.abb[end_month], ".tif"), overwrite = T)
    writeRaster(Water_limit, paste0(new_dir, "Water\\", per, "\\Water_rain_",per, "_",crop_name, "_", month.abb[start_month], "-", month.abb[end_month], "_higher_", limit, "_", province_name, ".tif"), overwrite = T)
  }
} 
  

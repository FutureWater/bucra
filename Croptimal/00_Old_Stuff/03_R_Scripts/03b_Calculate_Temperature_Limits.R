### Select directories abd files
input_files_tmin = list.files(paste0(results_dir, province_name, "\\Temperature\\tmin"), pattern = "\\.tif$", full.names = T, recursive =T)
input_files_tmax = list.files(paste0(results_dir, province_name, "\\Temperature\\tmax"), pattern = "\\.tif$", full.names = T, recursive =T)
new_dir = paste0(results_dir,province_name, "\\_LS_Results\\Temperature\\")


for(p in 1:(length(T_perc)+1)){
  if(p > length(T_perc)){
    per <- paste0("Monthly_Mean")
    folder_name <- "Tbase_Tupper\\"
    dir.create(paste0(new_dir,folder_name), recursive = T)
    start_name <- "T_two_limits_W_"
  }else{
    per <- paste0(gsub("\\.", "" ,T_perc[p]), "perc")
    folder_name <- paste0("Tmax_",per,"\\")
    dir.create(paste0(new_dir,folder_name), recursive = T)
    start_name <- paste0("T_two_limits_Tmax_",per,"_W_")
  }
  
  for(c in 1:nrow(cropping_cal)){
    crop_name <- cropping_cal$Crop[c]
    start_month <- cropping_cal$Start_growing_season[c]
    end_month <- cropping_cal$End_growing_season[c]
    if(start_month > end_month){
      months <- c(seq(start_month,12,1), seq(1, end_month, 1))
    }else{
      months <- seq(start_month, end_month, 1) 
    }
    weights <- as.numeric(cropping_cal[c,match("w_1", names(cropping_cal)):match(paste0("w_",length(months)), names(cropping_cal))])
    T_Base <- cropping_cal$T_Base[c]
    T_Upper <- cropping_cal$T_Upper[c]
    tmin_files <- input_files_tmin[grep(paste(month.abb[months],collapse="|"), input_files_tmin)]
    tmin_files <- stack(tmin_files[grep("Monthly_Mean", tmin_files)])
    tmax_files <- input_files_tmax[grep(paste(month.abb[months],collapse="|"), input_files_tmax)]
    tmax_files <- stack(tmax_files[grep(per, tmax_files)])
    
    Tmin_W <- sum(tmin_files * weights)
    Tmax_W <- sum(tmax_files * weights)
    
    Tmin_W_higher_Tbase <- Tmin_W>T_Base
    Tmax_W_lower_Tupper <- Tmax_W<T_Upper
    Two_limits_W <- Tmax_W_lower_Tupper*Tmin_W_higher_Tbase
    
    writeRaster(Two_limits_W, paste0(new_dir, folder_name, start_name, crop_name, "_", month.abb[start_month], "-", month.abb[end_month], ".tif"), overwrite=TRUE)
  
  }
}

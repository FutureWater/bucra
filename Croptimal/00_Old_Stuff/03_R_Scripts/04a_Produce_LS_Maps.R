new_dir = paste0(results_dir,province_name, "\\_LS_Results\\_Weighted\\")
dir.create(new_dir)

input_files_NDVI = list.files(paste0(results_dir,province_name, "\\_LS_Results\\NDVI"), pattern = "\\.tif$", full.names = T, recursive = T)
input_files_Temperature = list.files(paste0(results_dir,province_name, "\\_LS_Results\\Temperature"), pattern = "\\.tif$", full.names = T, recursive = T)
input_files_Water = list.files(paste0(results_dir,province_name, "\\_LS_Results\\Water"), pattern = "\\.tif$", full.names = T, recursive = T)
input_files_HHS = list.files(paste0(results_dir,province_name, "\\_LS_Results\\Soil_Hydraulic_Properties"), pattern = "\\.tif$", full.names = T, recursive = T)
input_files_SNC = list.files(paste0(results_dir,province_name, "\\_LS_Results\\Soil_Nutrient_Content"), pattern = "\\.tif$", full.names = T, recursive = T)
input_files_Slope = list.files(paste0(results_dir,province_name, "\\_LS_Results\\Elevation"), pattern = "\\.tif$", full.names = T, recursive = T)

for(c in 1:nrow(cropping_cal)){
  crop = cropping_cal$Crop[c]
  start_month = month.abb[cropping_cal$Start_growing_season[c]]
  end_month = month.abb[cropping_cal$End_growing_season[c]]
  for(p_T in 1:(length(T_perc)+1)){
    if(p_T > length(T_perc)){
      warm = "Average"
      per_T = "Tbase_Tupper"
    }else{
      warm = T_perc_names[p_T]
      per_T = paste0("Tmax_", gsub("\\.", "", T_perc[p_T]), "perc")
    }
    for(p_P in 1:(length(P_perc)+1)){
      if(p_P > length(P_perc)){
        dry = "Average"
        per_P = "Mean_Monthly"
      }else{
        dry = P_perc_names[p_P]
        per_P = paste0(gsub("\\.", "", P_perc[p_P]), "perc")
      }
      file_name = paste0("Land_Suitability_",dry, "_", warm,"_", crop, "_", start_month, "-", end_month, ".tif")
      
      NDVI_weighted = raster(input_files_NDVI[grep(paste0(start_month, "-", end_month), input_files_NDVI)]) * params$Weight[grep("NDVI", params$Parameter)]
      Temperature_fin = input_files_Temperature[grep(per_T, input_files_Temperature)]
      Temperature_weighted = raster(Temperature_fin[grep(paste0(crop, "_", start_month, "-", end_month), Temperature_fin)]) * params$Weight[grep("Temperature", params$Parameter)]
      Water_fin = input_files_Water[grep(per_P, input_files_Water)]
      Water_weighted = raster(Water_fin[grep(paste0(crop, "_", start_month, "-", end_month, "_higher_"), Water_fin)]) * params$Weight[grep("Water", params$Parameter)]
      Ksat_weighted = raster(input_files_HHS[grep("Ksat", input_files_HHS)]) * params$Weight[grep("Ksat", params$Parameter)]
      WCavail_weighted = raster(input_files_HHS[grep("WCavail", input_files_HHS)]) * params$Weight[grep("WCavail", params$Parameter)]
      Potasium_weighted = raster(input_files_SNC[grep("_K_", input_files_SNC)]) * params$Weight[grep("Potasium", params$Parameter)]
      Phosphorus_weighted = raster(input_files_SNC[grep("_P_", input_files_SNC)]) * params$Weight[grep("Phosphorus", params$Parameter)]
      Slope_weighted = raster(input_files_Slope[grep("lower", input_files_Slope)]) * params$Weight[grep("Slope", params$Parameter)]
      
      LS_r = sum(NDVI_weighted, Temperature_weighted, Water_weighted, Ksat_weighted, WCavail_weighted, Potasium_weighted, Phosphorus_weighted, Slope_weighted)
      
      writeRaster(LS_r, paste0(new_dir, file_name), overwrite =T)
    }
  }
  
}
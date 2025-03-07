for(par in unique(params$LS_Results_Folder)){
  parameter = params$Parameter[grep(par, params$LS_Results_Folder)]
  weight = params$Weight[grep(par, params$LS_Results_Folder)]
  limit = params$Limit[grep(par, params$LS_Results_Folder)]
  units = params$Units[grep(par, params$LS_Results_Folder)]
  folder = par
  abrev = params$Abrev[grep(par, params$LS_Results_Folder)]
  
  new_dir = paste0(results_dir,province_name, "\\_LS_Results\\", folder, "\\")
  dir.create(new_dir)
  if(folder == "Elevation" | folder == "Water"){next} ## Slope and Water limits are calculated in script 00_XXX.R and 03_XXX.R
  if(folder == "Soil_Hydraulic_Properties"){
    
    resampled_dem_diff_ras = raster(paste0(results_dir,province_name,"\\DEM\\DEM_",province_name,"_",res,"m_diff.tif"))
    boundary_proj = projectRaster(resampled_dem_diff_ras, crs = "EPSG:4326")
    boundary_ext = extent(boundary_proj)
    ext_ras = raster(boundary_ext, res = res(boundary_proj), crs = crs(boundary_proj))
    HHS_data = list.files(HHS_data_dir, pattern = "\\.tif$", full.names = T, recursive = T)
    HHS_data = stack(HHS_data[grep(paste(parameter, collapse = "|"), HHS_data)])
    
    HHS_data_cropped_10000 = crop(HHS_data,ext_ras)/10000
    HHS_data_cropped_10000_proj = projectRaster(HHS_data_cropped_10000, resampled_dem_diff_ras, res = res)
    HHS_data_weighted_Top = HHS_data_cropped_10000_proj[[grep("TOPSOIL", names(HHS_data_cropped_10000_proj))]]*0.3
    HHS_data_weighted_Sub = HHS_data_cropped_10000_proj[[grep("SUBSOIL", names(HHS_data_cropped_10000_proj))]]*1.7
    
    for(HHS in parameter){
      if(HHS == "Ksat"){
        multiplier = 10 #mm/d
      }
      if(HHS == "WCavail"){
        multiplier = 1000 # mm/m
      }
      lim = limit[grep(HHS, parameter)]
      unit = units[grep(HHS, parameter)]
      
      HHS_weighted = (HHS_data_weighted_Top[[grep(HHS, names(HHS_data_weighted_Top))]] + 
                        HHS_data_weighted_Sub[[grep(HHS, names(HHS_data_weighted_Sub))]])/2*multiplier
      HHS_limit = HHS_weighted > as.numeric(lim)
      
      writeRaster(HHS_limit, paste0(new_dir, HHS, "_Soil_higher_", lim, unit,"_", province_name,".tif"), overwrite =T)
      
    }
  }
  if(folder == "Soil_Nutrient_Content"){
    SNC_data = list.files(paste0(results_dir, province_name, "\\Soil_Nutrient_Content"), pattern = "\\.tif$", full.names = T, recursive = T)
    SNC_data = stack(SNC_data[grep(paste0("_", abrev, "_", collapse = "|"), SNC_data)])
    
    SNC_data_limit = SNC_data > as.numeric(limit)
    
    for(SNC in parameter){
      abr = abrev[grep(SNC, parameter)] 
      lim = limit[grep(SNC, parameter)]
      unit = units[grep(SNC, parameter)]
      writeRaster(SNC_data_limit[[grep(SNC, parameter)]], paste0(new_dir, "Extractable_", toupper(abr),"_Soil_higher_", lim, unit,"_", province_name,".tif"), overwrite =T)
    }
  }
  if(folder == "NDVI"){
    NDVI_data = list.files(paste0(results_dir, province_name, "\\NDVI"), pattern = "\\.tif$", full.names = T, recursive = T)
    
    for(c in 1:nrow(cropping_cal)){
      start_month = cropping_cal$Start_growing_season[c]
      end_month = cropping_cal$End_growing_season[c]
      if(start_month > end_month){
        months = month.abb[c(seq(start_month,12,1), seq(1, end_month, 1))]
      }else{
        months = month.abb[seq(start_month, end_month, 1)] 
      }
      NDVI = mean(stack(NDVI_data[grep(paste0(months, collapse = "|"), NDVI_data)]))
      NDVI_limit = NDVI > as.numeric(limit)
      
      writeRaster(NDVI_limit, paste0(new_dir, "NDVI_limit_higher_", gsub("\\.", "_", limit),"_", 
                                     month.abb[start_month], "-", month.abb[end_month], "_", province_name,".tif"), overwrite =T)
      
    }
  }
}

 
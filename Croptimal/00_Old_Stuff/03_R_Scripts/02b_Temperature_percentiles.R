for (var in T_vars){
  
  for (p_T in T_perc){
  
    new_dir = paste0(results_dir, province_name, "\\Temperature\\", var, "\\", gsub("\\.","", as.character(p_T)),"perc\\")
    dir.create(new_dir, recursive = T)
    input_files = list.files(paste0(data_dir, "Temperature\\",var) , pattern = "\\.nc$", full.names = T)
    
    if(exists(paste0("sstStack_T_", var, "_", gsub("\\.","", as.character(p_T)))) == F){
      
      sstStack_temp = stack(input_files)
      assign(paste0("sstStack_T_", var, "_", gsub("\\.","", as.character(p_T))), sstStack_temp)
      months_sstStack = rep(seq(1,12,1),length(input_files))
    }
    if(exists(paste0("outSt_T_", var, "_", gsub("\\.","", as.character(p_T)))) == F)   {
      ### Because StackApply doesnt work with percentiles, this is some sort of replacement 
      outSt <- stack()
      for (mn in 1:12){
        st <- subset(get(paste0("sstStack_T_", var, "_", gsub("\\.","", as.character(p_T)))), which(months_sstStack == mn))
        mn_perc <- calc(st, fun=function(x) raster::quantile(x, probs=p_T, na.rm=T))
        outSt <- addLayer(outSt, mn_perc)
      }
      assign(paste0("outSt_T_", var, "_", gsub("\\.","", as.character(p_T))), outSt)
    }
    
    ### crop raster to extent of the province
    crop_var_raster = crop(get(paste0("outSt_T_", var, "_", gsub("\\.","", as.character(p_T)))), gBuffer(province_shp,0.25, byid = T))
    
    ### import Diff_DEM of the province
    resampled_dem_diff_ras = raster(paste0(results_dir,province_name,"\\DEM\\DEM_",province_name,"_",res,"m_diff.tif"))
    
    ### resample var_raster to match dem_raster
    resampled_var_ras = projectRaster(crop_var_raster, resampled_dem_diff_ras, res = res, method = "bilinear")
    
    ### use lapse rate to correct temperature
    final_ras <- resampled_var_ras 
    
    ### write result to folder
    for (x in 1:nlayers(final_ras)) {
      writeRaster(final_ras[[x]], paste0(new_dir, month.abb[x],".tif"), overwrite = T)
    }
  }
}

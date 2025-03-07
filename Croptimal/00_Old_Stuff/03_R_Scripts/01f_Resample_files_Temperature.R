for (var in T_vars){
  input_files = list.files(paste0(data_dir, "Temperature\\",var) , pattern = "\\.nc$", full.names = T)
  new_dir = paste0(results_dir, province_name, "\\Temperature\\", var, "\\Monthly_Mean\\")
  dir.create(new_dir, recursive = T)
  
  stack_data = stack(input_files)
  index_names <- ifelse(nchar(rep(seq(1,12,1),length(input_files)))==1, paste0("0",rep(seq(1,12,1),length(input_files))),rep(seq(1,12,1),length(input_files)))
  data_mean <- stackApply(stack_data, index_names, fun = mean)
  
  crop_var_raster = crop(data_mean, gBuffer(province_shp,0.25, byid = T))
  ### import Diff_DEM of the province
  resampled_dem_diff_ras = raster(paste0(results_dir,province_name,"\\DEM\\DEM_",province_name,"_",res,"m_diff.tif"))
  
  ### resample var_raster to match dem_raster
  resampled_var_ras = projectRaster(crop_var_raster, resampled_dem_diff_ras, res = res, method = "bilinear")
  
  ### use lapse rate to correct temperature
  final_ras <- resampled_var_ras + resampled_dem_diff_ras*T_lapse_rate 
  
  lapply(1:12,function(x){
    writeRaster(final_ras[[x]], paste0(new_dir, month.abb[x],".tiff"),overwrite=T)
  })  
}

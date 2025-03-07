
input_files = list.files(paste0(data_dir, "Rainfall\\Monthly_Sum\\") , pattern = "\\.tif$", full.names = T)
new_dir = paste0(results_dir, province_name, "\\Rainfall\\Mean_Monthly\\")
dir.create(new_dir, recursive = T)
#################################################################################################

### Create output folder if that does not exist yet based on output_dir, var and province name
ifelse(dir.exists(new_dir)==F, dir.create(new_dir),paste0(""))

if(exists("var_stack_P")==F){
  var_stack_P = stack(input_files)
  var_names = names(var_stack_P)
  var_names_2 = gsub("X", "", var_names)
  var_names_2 = gsub("\\.", "_", var_names_2)
  var_names_2 = gsub("_P_tot", "", var_names_2)
  
  index_names = substr(var_names_2,6,7)
  data_mean_rainfall <- stackApply(var_stack_P, index_names, fun = mean)
}

### crop raster to extent of the province
crop_var_raster = crop(data_mean_rainfall, gBuffer(province_shp,0.25, byid = T))

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






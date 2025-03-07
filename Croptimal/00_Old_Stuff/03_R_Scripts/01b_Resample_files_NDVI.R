### Select directories abd files
input_files = list.files(paste0(data_dir, "NDVI\\",province_name,"\\Mean_Monthly\\") , pattern = "\\.tif$", full.names = T)
new_dir = paste0(results_dir, "\\", province_name, "\\NDVI\\Mean_Monthly\\")
dir.create(new_dir, recursive = T)

### Read raster file and select band (month)
var_raster = stack(input_files)
names_raster = names(var_raster)
### import Diff_DEM of the province
resampled_dem_diff_ras = raster(paste0(results_dir,province_name,"\\DEM\\DEM_",province_name,"_",res,"m_diff.tif"))

### crop raster to extent of the province
crop_var_raster = crop(var_raster, resampled_dem_diff_ras)

### resample var_raster to match dem_raster
resampled_var_ras = projectRaster(crop_var_raster, resampled_dem_diff_ras, res = res, method = "bilinear")

### use lapse rate to correct temperature
final_ras <- resampled_var_ras 

### plot result
lapply(1:12,function(x){
  writeRaster(final_ras[[x]], paste0(new_dir, names_raster[x],".tiff"),overwrite=T)
})

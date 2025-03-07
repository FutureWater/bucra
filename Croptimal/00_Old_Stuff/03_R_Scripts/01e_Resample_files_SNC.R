### Select directories abd files
input_files = list.files(paste0(data_dir, "Soil_Nutrient_Content\\") , pattern = "\\.tif$", full.names = T)
new_dir = paste0(results_dir, province_name, "\\Soil_Nutrient_Content\\")
dir.create(new_dir, recursive = T)

var_raster = stack(input_files)
var_names = names(var_raster)

var_names_2 = gsub("af", paste0(province_name,"_"), var_names)

### Create output folder if that does not exist yet based on output_dir, var and province name
ifelse(dir.exists(new_dir)==F, dir.create(new_dir),paste0(""))

### crop raster to extent of the province
crop_var_raster = crop(var_raster, gBuffer(province_shp,0.25, byid = T))
### import Diff_DEM of the province
resampled_dem_diff_ras = raster(paste0(results_dir,province_name,"\\DEM\\DEM_",province_name,"_",res,"m_diff.tif"))

### resample var_raster to match dem_raster
resampled_var_ras = projectRaster(crop_var_raster, resampled_dem_diff_ras, res = res, method = "bilinear")
    
### use lapse rate to correct temperature
final_ras <- resampled_var_ras 
    
### write result to folder
print(paste0("write SNC ", province_name))
for (x in 1:length(var_names)) {
  writeRaster(final_ras[[x]], paste0(new_dir,"\\",var_names_2[x],".tif"), overwrite = T)
} 




### Select directories abd files
input_files = list.files(paste0(data_dir, "Ref_ET\\") , pattern = "ymonmean.nc$", full.names = T)
new_dir = paste0(results_dir, "\\", province_name, "\\Ref_ET\\Mean_Monthly\\")
dir.create(new_dir, recursive = T)

multiplier = c(31,28,31,30,31,30,31,31,30,31,30,31) # This is used to multiply the daily ET_ref data with the
                                                    # number of days of the month
nc = nc_open(input_files)
time_nc = ncvar_get(nc,"time")
st_date = as.Date(time_nc,origin="1900-01-01")

index_names = substr(st_date, 5,8)

brick_data = brick(input_files)

data_mean<- stackApply(brick_data, index_names, fun = mean)

nc_close(nc)

### Clip var_raster to each province shape seperately
### Create output folder if that does not exist yet based on output_dir, var and province name
ifelse(dir.exists(new_dir)==F, dir.create(new_dir),paste0(""))

### crop raster to extent of the province
crop_var_raster = crop(data_mean, gBuffer(province_shp,0.25, byid = T))
### import Diff_DEM of the province
resampled_dem_diff_ras = raster(paste0(results_dir,province_name,"\\DEM\\DEM_",province_name,"_",res,"m_diff.tif"))

### resample var_raster to match dem_raster
resampled_var_ras = projectRaster(crop_var_raster, resampled_dem_diff_ras, res = res, method = "bilinear")

### use lapse rate to correct temperature
final_ras <- resampled_var_ras 

lapply(1:12,function(x){
  writeRaster(final_ras[[x]]*multiplier[x], paste0(new_dir, month.abb[x],".tiff"),overwrite=T)
})





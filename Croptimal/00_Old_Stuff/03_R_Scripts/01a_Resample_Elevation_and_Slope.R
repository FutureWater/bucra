print(paste0("crop DEM: ", province_name))
DEM_crop = crop(DEM_r, province_shp)
print(paste0("mask DEM: ", province_name))
DEM_mask = mask(DEM_crop, province_shp)

new_dir = paste0(results_dir, "\\", province_name, "\\DEM")
dir.create(new_dir, recursive = T)

province_shp_proj = spTransform(province_shp, Local_proj)
province_shp_proj_ext = extent(province_shp_proj)

r <- raster(xmn= province_shp_proj_ext[1], ymn= province_shp_proj_ext[3], xmx = province_shp_proj_ext[2],ymx = province_shp_proj_ext[4], 
            resolution = res,
            crs = Local_proj)

print(paste0("bilinear projectRaster DEM: ", province_name))
prov_proj_res = projectRaster(DEM_mask, r, res= res, method ='bilinear')

print(paste0("ngb projectRaster DEM: ", province_name))
prov_proj_res_ngb = projectRaster(DEM_mask, r, res= res, method ='ngb')

prov_proj_diff = prov_proj_res - prov_proj_res_ngb

print(paste0("write rasters DEM: ", province_name))
writeRaster(prov_proj_res,paste0(new_dir, "\\DEM_", province_name, "_", res, "m.tif"), overwrite = T)
writeRaster(prov_proj_diff,paste0(new_dir, "\\DEM_", province_name, "_", res, "m_diff.tif"), overwrite = T)


print(paste0("calculate slope DEM:", province_name))

indir_elev = paste0(indir, "_LS_Results\\Elevation\\")
dir.create(indir_elev, recursive = T)

x <- terrain(prov_proj_res, opt="slope", unit="tangent", neighbors=8)
x <- x*100
writeRaster(x, paste0(indir_elev, 'Slope_',province_name,".tif"), overwrite =T)
Slope = as.numeric(params$Limit[params$Parameter ==  "Slope"])
Slope_lower_limit <- x<Slope
writeRaster(Slope_lower_limit, paste0(indir_elev, 'Slope_lower_',Slope,'perc_',province_name,".tif"), overwrite =T)


new_dir_T = paste0(results_dir, "\\_LS_Results\\Temperature\\")
new_dir_P = paste0(results_dir, "\\_LS_Results\\Rainfall\\")
dir.create(new_dir_T, recursive = T)
dir.create(new_dir_P, recursive = T)

### Compare to temperature
T_avg_files = list.files(paste0(results_dir,province_name,"\\Temperature\\tavg\\"), pattern = "\\.tif$", full.names = T, recursive = T)
T_max_files = list.files(paste0(results_dir,province_name,"\\Temperature\\tmax\\"), pattern = "\\.tif$", full.names = T, recursive = T)
T_min_files = list.files(paste0(results_dir,province_name,"\\Temperature\\tmin\\"), pattern = "\\.tif$", full.names = T, recursive = T)

T_files_all = c(T_avg_files, T_max_files, T_min_files)

P_avg_files = list.files(paste0(results_dir,province_name,"\\Rainfall\\005perc\\"), pattern = "\\.tif$", full.names = T, recursive = T)
P_max_files = list.files(paste0(results_dir,province_name,"\\Rainfall\\025perc\\"), pattern = "\\.tif$", full.names = T, recursive = T)
P_min_files = list.files(paste0(results_dir,province_name,"\\Rainfall\\Mean_Monthly\\"), pattern = "\\.tif$", full.names = T, recursive = T)

P_files_all = c(P_avg_files, P_max_files, P_min_files)

communes_selected = communes_shp[communes_df$Province == province_name,]

for (T_file in 1:length(T_files_all)){
  print(paste0("file ", T_file, " of ", length(T_files_all)))
  temp <- raster(T_files_all[T_file])
  Tempi = basename(dirname(dirname(T_files_all[T_file])))
  perci = basename(dirname(T_files_all[T_file]))
  month = substr(basename(T_files_all[T_file]),0,3)
  
  name = paste0(Tempi, "_", perci, "_", month)
  if(T_file ==1){
    df_res = extract(temp, communes_selected, fun=mean, na.rm=TRUE)
    names(df_res) <- name
  }else{
    df2_res = extract(temp, communes_selected, fun=mean, na.rm=TRUE)
    names(df2_res) <- name
    df_res = cbind(df_res, df2_res)
  }
  names(df_res)[T_file] <-name
  
}

df_res2 = data.frame(df_res)
for (T_file in 1:length(T_files_all)){
  print(paste0("file ", T_file, " of ", length(T_files_all)))
  temp <- raster(T_files_all[T_file])
  Tempi = basename(dirname(dirname(T_files_all[T_file])))
  perci = basename(dirname(T_files_all[T_file]))
  month = substr(basename(T_files_all[T_file]),0,3)
  
  name = paste0(Tempi, "_", perci, "_", month)
  names(df_res2)[T_file] <-name
  
}

write.csv(df_res2, paste0(new_dir_T, "Average_Temperature_per_Commune_", province_name, ".csv"),
          col.names = T, row.names = communes_selected$NAME_3)


for (P_file in 1:length(P_files_all)){
  print(paste0("file ", P_file, " of ", length(P_files_all)))
  temp <- raster(P_files_all[P_file])
  Tempi = "P"
  perci = basename(dirname(P_files_all[P_file]))
  month = substr(basename(P_files_all[P_file]),0,3)
  
  name = paste0(Tempi, "_", perci, "_", month)
  if(P_file ==1){
    df_res = extract(temp, communes_selected, fun=mean, na.rm=TRUE)
    names(df_res) <- name
  }else{
    df2_res = extract(temp, communes_selected, fun=mean, na.rm=TRUE)
    names(df2_res) <- name
    df_res = cbind(df_res, df2_res)
  }
  names(df_res)[P_file] <-name
  
}

df_res2 = data.frame(df_res)
for (P_file in 1:length(P_files_all)){
  print(paste0("file ", P_file, " of ", length(P_files_all)))
  temp <- raster(P_files_all[P_file])
  Tempi = "P"
  perci = basename(dirname(P_files_all[P_file]))
  month = substr(basename(P_files_all[P_file]),0,3)
  
  name = paste0(Tempi, "_", perci, "_", month)
  names(df_res2)[P_file] <-name
  
}

write.csv(df_res2, paste0(new_dir_P, "Average_Rainfall_per_Commune_", province_name, ".csv"),
          col.names = T, row.names = communes_selected$NAME_3)

new_dir = paste0(results_dir, "\\_LS_Results\\_Communes\\")
dir.create(new_dir, recursive = T)

df = communes_df[communes_df$Province == province_name,]

all_results = list.files(results_dir, pattern = "Land_Suitability", recursive = T, full.names = T)
all_results = all_results[grep( "\\.tif$", all_results)]
all_results = all_results[grep(province_name, all_results)]
all_results_stack = stack(all_results)
mean_all_results_commune =  extract(all_results_stack, communes_shp_proj[gsub(" " , "_", communes_shp$NAME_1) == province_name, ],
                                     fun=mean, na.rm=TRUE, df=TRUE)

for(comb in unique_combis){
  df_temp = get(paste0(comb, "_Communes.csv"))
  mean_all_results_commune_comb = mean_all_results_commune[, grep(comb, names(mean_all_results_commune))]
  crops = unlist(strsplit(names(df_temp[4:ncol(df_temp)]), "_"))[seq(1,3*(ncol(df_temp)-3),3)]
  crop_message = cropping_cal$Name_message
  start_months = unlist(strsplit(names(df_temp[4:ncol(df_temp)]), "_"))[seq(2,3*(ncol(df_temp)-3),3)]
  end_months = unlist(strsplit(names(df_temp[4:ncol(df_temp)]), "_"))[seq(3,3*(ncol(df_temp)-3),3)]
  
  for(x in 1:length(crops)){
    col = grep(paste0(crops[x], "_", start_months[x], ".", end_months[x]),
                names(mean_all_results_commune_comb))
    row = grep(province_name, df_temp$Province)
    mean_all_results_commune_comb_message = ifelse(mean_all_results_commune_comb >= 0.8, "very suitable",
                                                   ifelse(mean_all_results_commune_comb >= 0.7, "suitable",
                                                          ""))
    final_result = ifelse(mean_all_results_commune_comb[col] < 0.7,
      sprintf("Growing %s from %s to %s is NOT ideal in your location. If you still want to plant this crop: You may have to plant a variety resistant to extreme temperatures. You may have to irrigate. You may have to fertilize a lot. You may have to dig drains.",
                          crop_message[col],
                          month.name[match(start_months[col], month.abb)],
                          month.name[match(end_months[col], month.abb)]),
      sprintf("Growing %s from %s to %s is %s in your location.",
                          crop_message[col],
                          month.name[start_month],
                          month.name[end_month],
                          mean_all_results_commune_comb_message[col])
    )

    df_temp = get(paste0(comb, "_Communes.csv"))
    df_temp[row, col+3] <- final_result
    assign(paste0(comb, "_Communes.csv"), df_temp)

  }
  if(province_name == provinces_names[1]){
    df_temp2 = get(paste0(comb, "_Communes.csv"))
    df_temp2[grep(province_name, df_temp2$Province),4:ncol(df_temp2)] <- round(mean_all_results_commune_comb,3)
    assign(paste0(comb, "_Communes_LS_values.csv"), df_temp2)
  }else{
    df_temp2 = get(paste0(comb, "_Communes_LS_values.csv"))
    df_temp2[grep(province_name, df_temp2$Province),4:ncol(df_temp2)] <- round(mean_all_results_commune_comb,3)
    assign(paste0(comb, "_Communes_LS_values.csv"), df_temp2)
  }
  
}



if(province_name == provinces_shp$NAME_1[length(provinces_names)]){
  for(comb in unique_combis){
    write.csv(get(paste0(comb, "_Communes.csv")), paste0(new_dir, comb, "_Communes.csv"),
              row.names = F, col.names = T, quote = F)
    write.csv(get(paste0(comb, "_Communes_LS_values.csv")), paste0(new_dir, comb, "_Communes_LS_values.csv"),
              row.names = F, col.names = T, quote = F)
  }
}else{next}

      



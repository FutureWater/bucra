#Script 000_Run_All.R heeft een loop waarin dit scriptje wordt ingeladen de loop is met for name in province name --> hij loopt over dit scriptje per provincie! 
# In dit geval loopt ie eerst over Bengo. het script is dus een output per provincie tot ie bij de laatste provincie is en dan schrijft ie een csv file! 

split_path <- function(x) if (dirname(x)==x) x else c(basename(x),split_path(dirname(x)))

Seasonal_Temp_files = list.files(Seasonal_Forecast_dir, pattern = "temperature", full.names = T, recursive = T)   #leest JSON temperatuur files uit
Seasonal_Prec_files = list.files(Seasonal_Forecast_dir, pattern = "rainfall", full.names = T, recursive = T)      #leest JSON Precipitation files uit
communes_selected = communes_shp[communes_df$Province == province_name,]                                          #attribute table van communes shapefile (dataframe)
result_df <- data.frame(communes_selected$NAME_3)                                                                 #Maakt dataframe met als enige kolom de communes

#Loop 1 --> eindresultaat is result_df waarbij voor iedere commune een gemiddelde per maand is van de Tavg, Tmin, Tmax

for(t in 1:length(Seasonal_Temp_files)){                                                                          # loopt over drie files (Tavg, Tmax en Tmin)
  Temp_file = Seasonal_Temp_files[t]                                                                              # in iedere loop wordt een andere temp file ingeladen als je temp_file doet laat ie de laatste zien (=Tmin)
  T_stat = split_path(Temp_file)[2]                                                                               # T_stat = "tmin"
  T_stat_data = RJSONIO::fromJSON(Temp_file)                                                                      # leest tmin json uit en geeft een tabel terug met alle waarden
  # https://stackoverflow.com/questions/2991514/prevent-unlist-to-drop-null-values
  temp = data.frame("X" = T_stat_data$GridDefinition$Longitude,                                                   #dataframe met X en Y coordinaten
                    "Y" = T_stat_data$GridDefinition$Latitude,
                    "M0"= as.numeric(as.character(T_stat_data$Data$`0`$Data)),
                    "M1"= as.numeric(as.character(T_stat_data$Data$`1`$Data)),
                    "M2"= as.numeric(as.character(T_stat_data$Data$`2`$Data)))
                    #"M3"= as.numeric(as.character(T_stat_data$Data$`3`$Data))
                    #"M4"= as.numeric(as.character(T_stat_data$Data$`4`$Data)))
  
  dates = c(T_stat_data$Data$`0`$Date,                                                                            # 5 datums 
            T_stat_data$Data$`1`$Date,
            T_stat_data$Data$`2`$Date)#,      
            #T_stat_data$Data$`3`$Date,
            #T_stat_data$Data$`4`$Date)
  dates_T = dates                                                                                               # diezelfde dates 
  names(temp) <- c("X", "Y", dates)                                                                             # weer diezelfde dates "20220501" "20220601" "20220701" "20220801" "20220901"
  coordinates(temp) <- ~ X + Y                                                                                  #  X en Y coordinaten dataframe
  gridded(temp) <- TRUE                                                                                         # true
  temp2 = stack(temp)                                                                                           # temp2 is rasterbrick met temperatuurwaarden met kolommen voor iedere maand voor heel angola 
  crs(temp2) <- "EPSG:4326"                                                                                     # in goede projectie zetten
  temp2[temp2==-9999] <- NA                                                                                     # -9999 waarden verwijderen
  T_ras = temp2[[1]]                                                                                            # Raster met waarden tussen de 15.9 en 26.8 voor de eerste maand

  
  ### import Diff_DEM of the province
  resampled_dem_diff_ras = raster(paste0(results_drop,province_name,"\\DEM\\DEM_",province_name,"_",res,"m_diff.tif"))  #DEM per provincie importeren
  
  ### resample var_raster to match dem_raster

  resampled_var_ras = projectRaster(temp2, resampled_dem_diff_ras, res = res, method = "bilinear")                      #raster resample --> project values of a raster object to a new raster object with another projection

  ### use lapse rate to correct temperature
  final_ras <- resampled_var_ras + resampled_dem_diff_ras*T_lapse_rate                                                  #raster van de hele provincie met temperatuurwaarden per pixel en dan zijn er 5 kolommen die ieder een maand weergeven
  #plot(final_ras[[1]])
  
  communes_selected = communes_shp[communes_df$Province == province_name,]                                              #shapefile van alle communes in de geselecteerde provincie
 
  months = month.abb[as.numeric(substr(dates, 5,6))]                                                                    # May, Jun, Jul, Aug, Sep
  months_num = as.numeric(substr(dates, 5,6))                                                                           # 5,6,7,8,9
  
  crops = cropping_cal[cropping_cal$Start_growing_season >= months_num[1],]                                             # inladen cropping calendar
  

  for(mon in 1:length(months)){                                                                                         # for loop voor de x maanden die meegenomen worden
    month = months[mon]                                                                                                 # de maand waar de loop zich op dat moment in bevind; bijv "Sep"
    #Monthly_Mean

    communes_mean = raster::extract(temp2[[mon]], communes_selected, fun=mean, na.rm=TRUE)                              # dataframe van de gemiddelde/minimum/maximum temperatuur per commune in de desbetreffende maand

    result_df[ncol(result_df)+1] <- communes_mean                                                                       # voegt telkens de nieuw communes_column (tmax/min/avg en m 1,2,3,4 of 5) toe aan de dataframe iedere loop
    colnames(result_df)[ncol(result_df)] <- paste0("mean_", T_stat, "_", month)                                         # geeft naam aan de columnheader
    }
}

#loop 2

for(p in 1:length(Seasonal_Prec_files)){                                                                                # loopt over drie files (Pavg, Pmax en Pin)
  Prec_file = Seasonal_Prec_files[p]                                                                                    # locatie j.son rainfall.min file
  P_stat = split_path(Prec_file)[2]                                                                                     # "pmin"
  P_stat_data = RJSONIO::fromJSON(Prec_file)                                                                            # leest j.son file en maakt er een list van met de (in dit geval) Pmin waarden per pixel
  temp = data.frame("X" = P_stat_data$GridDefinition$Longitude,                                                         # dataframe met x en y coordinaten en 5 lege colums voor de datums (datums zijn aparte columns)
                    "Y" = P_stat_data$GridDefinition$Latitude,
                    "M0"= as.numeric(as.character(P_stat_data$Data$`0`$Data)),
                    "M1"= as.numeric(as.character(P_stat_data$Data$`1`$Data)),
                    "M2"= as.numeric(as.character(P_stat_data$Data$`2`$Data)))
                    #"M3"= as.numeric(as.character(P_stat_data$Data$`3`$Data)),
                    #"M4"= as.numeric(as.character(P_stat_data$Data$`4`$Data)))
  dates = c(P_stat_data$Data$`0`$Date,                                                                                  # 5 strings die de maanden weergeven 
            P_stat_data$Data$`1`$Date,
            P_stat_data$Data$`2`$Date)#,
            #P_stat_data$Data$`3`$Date,
            #P_stat_data$Data$`4`$Date)
  # names(temp) <- c("X", "Y", dates)
  names(temp) <- c("X", "Y", dates_T)                                                                                   #namen van de headers
  temp4 <- data.matrix(temp)                                                                                            # xy ingevuld data nog niet
  temp3 <- rasterize(temp4[, 1:2], T_ras, temp4[,3:ncol(temp4)], fun=mean)                                              # assign XY points to grid cells van temperatuur raster
  temp2 = stack(temp3)                                                                                                  # rasterbrick met regendata voor de 5 datums(maanden) per pixel 
  crs(temp2) <- "EPSG:4326"                                                                                             # zet projectie goed
  temp2[temp2==-9999] <- NA                                                                                             # verwijder alle -9999 waarden

  ### import Diff_DEM of the province
  resampled_dem_diff_ras = raster(paste0(results_drop,province_name,"\\DEM\\DEM_",province_name,"_",res,"m_diff.tif"))  # DEM per provincie importeren
  
  ### resample var_raster to match dem_raster
  final_ras = projectRaster(temp2, resampled_dem_diff_ras, res = res, method = "bilinear")                              #final_ras met waarden van de regenval voor de provincie waar je doorheen loopt op dat moment per pixel voor 5 mnd
  
  #plot(final_ras[[1]])
   
  communes_selected = communes_shp[communes_df$Province == province_name,]                                              # shapes van de communes in desbetreffende provincie 
  
  # months = month.abb[as.numeric(substr(dates, 5,6))]
  # months_num = as.numeric(substr(dates, 5,6))
  months = month.abb[as.numeric(substr(dates_T, 5,6))]                                                                  #may jun jul aug sep 
  months_num = as.numeric(substr(dates_T, 5,6))                                                                         # 5,6,7,8,9
  
  for(mon in 1:length(months)){                                                                                         # loop through 5 maanden
    month = months[mon]                                                                                                 # geeft de maand terug waar de loop zich nu in bevind
    #Monthly_Mean
    communes_mean = raster::extract(temp2[[mon]], communes_selected, fun=mean, na.rm=TRUE)                              # dataframe van de gemiddelde/minimum/maximum (waar de loop zich bevind) regenval per commune in de desbetreffende maand
    result_df[ncol(result_df)+1] <- communes_mean                                                                       # voegt aan de dataframe (ook met temp waarden) ook de nieuwe precipitation kolomen toe, maar ligt eraan waar je je bevind in de loop (bijv Tavg_sep)
    colnames(result_df)[ncol(result_df)] <- paste0("mean_", P_stat, "_", month)                                         # voegt header toe aan net toegevoegde column
  }
}


P_stats = read.csv(paste0(results_drop, "_LS_Results\\Rainfall\\Average_Rainfall_per_Commune_", province_name, ".csv"))                # inlezen data average rainfall per commune (005per en 025 perc en Pmean_monthly)
T_stats = read.csv(paste0(results_drop, "_LS_Results\\Temperature\\Average_Temperature_per_Commune_", province_name, ".csv"))          # inlezen data average temp per communce per maand (tmin: MM, 95perc, 75perc (hetzelfde voor tmax en tavg))

df_results_T_P = result_df                                                                                              #eindresultaat: dataframe: per commune per maand; tavg, tmax, tmin, pavg, pmax, pmin (volgens de laatste WI voorspellingen)

#loop 3

for(comm in 1:nrow(result_df)){                                                                                         # Loopen door de communes van de provincie waar je nu in bent 
  row = comm                                                                                                            # rij selectie
  ncolums = ncol(result_df)                                                                                             # aantal kolommen (=31) van df_results
  for(month in 1:length(months)){                                                                                       # door de maanden heen loopen
    cols_stats = grep(months[month], names(P_stats))                                                                    # zoekt naar de kolommen in P_stats die relevant zijn voor de maand waar je in bent (bijv voor sept; 13, 25 en 37) kolommen zijn t zelfde in T_stats
    P_values_stats = P_stats[row, cols_stats]                                                                           # pakt voor de commune waar je in bent de waarden voor de drie kolommen die net zijn geselecteerd
    T_values_stats = T_stats[row, cols_stats]                                                                           # pakt voor de commune waar je in bent de waarden voor de drie kolommen die net zijn geselecteerd
    
    cols_SF = grep(months[month], names(df_results_T_P))                                                                # in de grote df_results_T_P dataframe worden hier de kolommen van de goede maand geselecteerd (zijn er 6)
    values_SF = df_results_T_P[row, cols_SF]                                                                            # pakt per commune voor de goede kolommen de juiste waarden
    
    P_res = ifelse(values_SF[4] <= P_values_stats[1], "Much_Drier",                                                     # vergelijkt de Values_SF waarden van de kolom Pavg met P_values_stat kolom P_005perc (als true; much drier)
                   ifelse(values_SF[4] <= P_values_stats[2], "Drier", "Average"))                                       # anders drier en anders average
    T_res = ifelse(values_SF[1] >= T_values_stats[2], "Much_Warmer",                                                    # t zelfde maar dan voor de temperatuur (je krijg dat een Much_drier_Much_warmer (of andere variatie) voor die maand voor die commune qua weerssituatie)
                   ifelse(values_SF[1] >= T_values_stats[2], "Warmer", "Average"))
    
    df_results_T_P[row, ncolums + month] <- paste0(P_res, "_", T_res)                                                   # extra kolom toegevoegd aan df_results_T_P net de voorspellingen voor alle communes van de weersverwachting (drier_average of variatie) voor de maand waar de loop is
    colnames(df_results_T_P)[ncol(df_results_T_P)] <- paste0("Forecast_", months[month])                                # header een naam geven

  }
}

#output loop 3 is een hele grote dataframe met alle communes en de verwachting van het weer (1 van de 9 categorien) voor iedere maand
#Loop 4 

for(month in 1:length(months)){                                                                                       # loopen door de maanden! 
  print(month)
  forecast_df = data.frame(df_results_T_P[,1],df_results_T_P[,grep(months[month], names(df_results_T_P))])            # dataframe voor de maand waar je in je loop zit met de P+T waarden, de forecast_month en daarachter de planting suitability waarde voor de crops voor die weerscategorie 
  month_num = match(months[month],month.abb)                                                                          # getal van de maand waar we in zitten (sep = 9)
  all_crops = cropping_cal[cropping_cal$Start_growing_season == month_num |                                           # cropping calendar inladen
                           cropping_cal$End_growing_season >= month_num   |
                           cropping_cal$End_growing_season <= 12,]
  crops = paste0(all_crops$Crop, "_", month.abb[all_crops$Start_growing_season], "_",                                 # headers alle crops 
                 month.abb[all_crops$End_growing_season])
  crop_col_name = paste0(all_crops$Crop)
  
  #if(crops == "__"){next}                                                                                             # als false; dan kunnen we naar de volgende loop
 
  #loop 4.2 --> voor maand x gaan we nu loopen door alle crops 
  for(crop in 1:length(crops)){                                                                                  
        forecast_df[,ncol(forecast_df)+1] = NA                                                                        # extra kolom met NA values aanmaken voor crop x
        colnames(forecast_df)[ncol(forecast_df)] <- paste0(crops[crop])                                               # naamgeving kolom values naar crop x in het df
        #print(forecast_df) # we hebben nu een dataframe met iedere loop een rij NA waarden voor de crops

  #loop 4.3 --> voor crop x en maand x gaan we nu door de communes loopen in de provincie waar we in zijn
        for(comm in 1:nrow(forecast_df)){   
          # for commune X
          if(length(grep("NA", forecast_df[comm, 8]))>0){next}                                                           
            temp_csv = read.csv(paste0(new_dir, "Land_Suitability_", forecast_df[comm, 8], "_Communes_LS_values.csv")) # avg_warmer tabel of iets dergelijks lezen
            value = temp_csv[temp_csv$Province == province_name & temp_csv$Commune == forecast_df[comm,1],             # leest crop column uit de tabel voor commune
                            grep(crops[crop], names(temp_csv))][1]
            
            avg_avg = read.csv(paste0(new_dir, "Land_Suitability_Average_Average_Communes_LS_values.csv"))             # average_average tabel uitlezen
            avg_value = temp_csv[temp_csv$Province == province_name & temp_csv$Commune == forecast_df[comm,1],         # leest crop column uit de tabel voor commune
                                 grep(crops[crop], names(temp_csv))][1]
            
            threshold = avg_value - (avg_value*0.25)
            
            
            #print(final_forecast_df)
            #print(value)}}}

            #print(forecast_df[comm, ncol(forecast_df)]) is hier nog NA
            forecast_df[comm, ncol(forecast_df)] <- as.numeric(value)                                                  # value aan laatste kolom van forecast_df geven       
            
          
            #final_forecast_df vullen
            #average kolom aanmaken
            final_forecast_df$Average_Average[final_forecast_df$Crop == crops[crop] & 
                                                              final_forecast_df$Commune == forecast_df[comm,1] &
                                                              final_forecast_df$Province == province_name] <- as.numeric(avg_value)
            #threshold kolom aanmaken
            final_forecast_df$Threshold[final_forecast_df$Crop == crops[crop] & 
                                                final_forecast_df$Commune == forecast_df[comm,1] &
                                                final_forecast_df$Province == province_name] <- as.numeric(threshold)
            #warning kolom aanmaken M1
            final_forecast_df <- final_forecast_df %>% 
              mutate(M1 = if_else(Planting_Suitability_M1 <= Threshold, "Warning!", "No"))
    
            #final_forecast_df$Reduced_Suitability_M1 <- "No"
            #final_forecast_df$Reduced_Suitability_M1[final_forecast_df$Planting_Suitability_M1<=final_forecast_df$Threshold] <- "Warning!"
            #final_forecast_df$Reduced_Suitability_Start <- ifelse(final_forecast_df$Planting_Suitability_M1 <= final_forecast_df$Threshold, "No", "Warning!") #alt
          
            #voorspelling kolom vullen
            final_forecast_df[final_forecast_df$Crop == crops[crop] & 
                                final_forecast_df$Commune == forecast_df[comm,1] &
                                final_forecast_df$Province == province_name, 
                              grep(month.abb[month_num], names(final_forecast_df))] <- forecast_df[comm,8]            
            
            #planting suitability kolom per maand vullen
            final_forecast_df[final_forecast_df$Crop == crops[crop] &
                                final_forecast_df$Commune == forecast_df[comm,1] &
                                final_forecast_df$Province == province_name,
                              grep(month, names(final_forecast_df))] <- as.numeric(value)
            
            #zorgen dat alleen de goede crops geselecteerd worden
            if(month_num == all_crops$Start_growing_season[crop]){ # Note Corjan: this throws a FALSE sometimes --> What to do?  
              final_forecast_df$Crop_Start_Month[final_forecast_df$Crop == crops[crop] & 
                                                                final_forecast_df$Commune == forecast_df[comm,1] &
                                                                final_forecast_df$Province == province_name] <- month_num #as.numeric(value)
            }
        }
  }
}


#wegschrijven naar csv bestand en laatste aanpassingen dataframe
if(province_name == provinces_names[length(provinces_names)]){
  
  # NA is no_information
  final_forecast_df$Crop_Start_Month <- replace(final_forecast_df$Crop_Start_Month, 
                                                is.na(final_forecast_df$Crop_Start_Month), "No Information")
  
  final_forecast_df <- final_forecast_df %>% 
    mutate(M2 = if_else(Planting_Suitability_M2 <= Threshold, "Warning!", "No")) 
  
  final_forecast_df <- final_forecast_df %>% 
    mutate(M3 = if_else(Planting_Suitability_M3 <= Threshold, "Warning!", "No"))
  
  final_forecast_df <- final_forecast_df[ , colSums(is.na(final_forecast_df)) <nrow(final_forecast_df)]   #delete NA rows
  
  final_forecast_df$Message = ifelse(final_forecast_df$Crop_Start_Month == "No Information", "No Information", ifelse(final_forecast_df$M1 == "Warning!" | final_forecast_df$M2  == "Warning!" | final_forecast_df$M3 == "Warning!", "Warning!", "No Reduced Suitability"))

  final_forecast_df <- subset(final_forecast_df, select= -c(Crop_Start_Month))
  
  #after the message, also provide the reason
  
  final_forecast_df$dryness <- gsub("_", " ", ifelse(final_forecast_df[[4]] == "Much_Drier_Much_Warmer" | final_forecast_df[[4]] == "Much_Drier_Warmer" | final_forecast_df[[4]] == "Much_Drier_Average", "Much_Drier", 
                                      ifelse(final_forecast_df[[4]] == "Drier_Much_Warmer" | final_forecast_df[[4]] == "Drier_Warmer" | final_forecast_df[[4]] == "Drier_Average", "Drier", "Average")))
  
  final_forecast_df$temperature <- gsub("_", " ", ifelse(final_forecast_df[[4]] == "Much_Drier_Much_Warmer" | final_forecast_df[[4]] == "Drier_Much_Warmer" | final_forecast_df[[4]] == "Average_Much_Warmer", "Much_Warmer", 
                                   ifelse(final_forecast_df[[4]] == "Much_Drier_Warmer" | final_forecast_df[[4]] == "Drier_Warmer" | final_forecast_df[[4]] == "Average_Warmer", "Warmer", "Average")))
  
  final_forecast_df$Reason <- ifelse(final_forecast_df$dryness == "Average" & final_forecast_df$temperature == "Average", "No Message", 
                                     ifelse(final_forecast_df$dryness == "Average", paste(final_forecast_df$temperature, "than normal"), 
                                            ifelse(final_forecast_df$temperature == "Average", paste(final_forecast_df$dryness, "than normal"), paste(final_forecast_df$dryness, "and", final_forecast_df$temperature, "than normal"))))
  
  final_forecast_df$Crop <- sub("_.*", "", final_forecast_df$Crop)
  
  final_forecast <- select(final_forecast_df, Province, Commune, Crop, Message, Reason)    #create final dataframe

  # zorgen dat reshape 2 als package wordt geinstalleerd en library
  #M1 <- dcast(final_forecast, fun.aggregate = list, Province + Commune ~ Crop)
  
  write.csv(final_forecast, paste0(new_dir, "Forecast_all_communes_", months[1], "_Latest.csv"))#,row.names = F, col.names(M1) = T, quote = F)
  
  # setup FTP parameters
  #file_in <- ftp_outname
  #file_out <- str_sub(file_in, 14)
  #path_out <- paste0("ftp://213.127.133.58/MavoDiami/Dynamic_Planting_Suitability/",file_out)
  #userpwd <- "MavoDiami:FW@MavoDiami2020!"
  
  # Push to FTP
  #ftpUpload(file_in, path_out, userpwd=userpwd)
  
}




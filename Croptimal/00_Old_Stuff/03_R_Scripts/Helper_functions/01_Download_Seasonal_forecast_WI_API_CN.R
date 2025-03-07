library(httr)
library(jsonlite)
library(R.utils)

import_info = read.csv(paste0(Seasonal_Forecast_dir, "Seasonal_structure.csv"))
url = 'https://service.weatherimpact.com/api/data/mavo_diami_seasonal/'
authkey = '0a970d7bcb110aaa2c1f5afcea1129ef'


## Import seasonal files
getMet <- function(url, endpoint, authkey, import_info, attempts){
  
  url_st = paste0(url, endpoint)
  
  request <- NULL
  attempt <- 0
  while(is.null(request) && attempt <= attempts) {
    attempt <- attempt + 1
    print(paste0("Attempt ", attempt, " to download data for: ", import_info$Par[import_info$Endpoint == endpoint]))
    try(
      request <- withTimeout(GET(url = url_st, add_headers(.headers = c('authkey'= authkey))),
                               timeout = 10, onTimeout = "silent")
    )
  } 
  response <- httr::content(request)
}

## Run the function for all parameters

for(info in 1:nrow(import_info)){
  print(import_info$Par[info])
  get_res <- getMet(url, import_info$Endpoint[info], authkey, import_info, 50)
  print(paste0("successfully imported ", import_info$Par[info]))
  outfolder <- ifelse(length(grep("rainfall", import_info$Output[info])==1), "Rainfall/", "Temperature/") 
  write_json(get_res, paste0(Seasonal_Forecast_dir, outfolder,  import_info$Par[info], "/", import_info$Output[info]))
}


extract_pixelValues <- function(targetLong = 23.61, targetLat = -13.45){
  
  ##  
  ## libs
  ##
  library(httr)
  library(RJSONIO)
  
  targetLong = 23.61
  targetLat = -13.45
  
  ##
  ## init
  ## 
  apiCallVector <- c("forecast_class", "forecast_text", 
                     "rainfall_min", "rainfall_median", "rainfall_max")
                     # "temperature_min", "temperature_median", "temperature_max")
  dateValue <- NULL
  resultContainer <- data.frame()
  
  
  
  ## 
  ## collect datalayers
  ##
  for(nn in seq(from=1, to=length(apiCallVector))){
    
    print(sprintf("... now processing %s", apiCallVector[nn]))
    
    
    query <- sprintf("https://service.weatherimpact.com/api/data/mavo_diami_seasonal/%s?datetime=latest", apiCallVector[nn])
    fc.param <- GET(url=query, add_headers(Authorization="bearer 0a970d7bcb110aaa2c1f5afcea1129ef"))
    fc.param.j <- fromJSON(content(fc.param,type="text"))
    
    
    if (nn==1){
    ## 
    ## calculate which element of the grid vector are closest
    ##
    dLong <-abs(fc.param.j$GridDefinition$Longitude-targetLong)
    indexLong <- which(dLong == min(dLong))
    dLat <- abs(fc.param.j$GridDefinition$Latitude-targetLat)
    indexLat <- which(dLat == min(dLat))
    
    
    ## 
    ## identify the element which is closest in both the lat and long vectors
    ##
    indexCombined <- intersect(indexLat, indexLong)
    }
    
    
    dataValue <- NULL
    for (mm in seq(from=1, to=5)){
      ## 
      ## extract data value @ of closest pixel
      ##
      if (nn==1){
        dateValue <- c(dateValue, fc.param.j$Data[[mm]]$Date)
      }
      dataValue <- c(dataValue, fc.param.j$Data[[mm]]$Data[indexCombined])
    }
    
    if (nn==1){
      resultContainer <- data.frame(date=dateValue)
      resultContainer <- cbind(resultContainer, dataValue)
    }
    else {
      resultContainer <- cbind(resultContainer, dataValue)
    }
    
    dd <- 1
  }
  
  colnames(resultContainer) <- c("date", apiCallVector)
  rownames(resultContainer) <- as.character(seq(from=1, to=nrow(resultContainer)))
  
  return(resultContainer)
  
}








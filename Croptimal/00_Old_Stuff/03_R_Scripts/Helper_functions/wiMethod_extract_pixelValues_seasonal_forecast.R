extract_pixelValues_seasonal_forecast <- function(targetLong = 23.61, targetLat = -13.45, forecastType="rain"){
  
  
  ##
  ## forecastType can be either "rain" or "temperature"
  ##  
  
  
  
  ##
  ## libs =========================================================================================
  ##
  library(httr)
  library(RJSONIO)
  
  
  
  ##
  ## init =========================================================================================
  ## 
  if (forecastType=="rain"){
    apiCallVector <- c("forecast_class", "forecast_text", 
                       "rainfall_min", "rainfall_median", "rainfall_max")
  }
  else if (forecastType=="temperature"){
    apiCallVector <- c("temperature_min", "temperature_median", "temperature_max")
  }
  else {
    stop("Dummy! You entered an incorrect forecast type!")
  }
  
  dateValue <- NULL
  resultContainer <- data.frame()
  
  
  
  ## 
  ## collect datalayers ===========================================================================
  ##
  for(nn in seq(from=1, to=length(apiCallVector))){
    
    print(sprintf("... now processing %s", apiCallVector[nn]))
    
    query <- sprintf("https://service.weatherimpact.com/api/data/mavo_diami_seasonal/%s?datetime=latest", apiCallVector[nn])
    fc.param <- GET(url=query, add_headers(Authorization="bearer 0a970d7bcb110aaa2c1f5afcea1129ef"))
    fc.param.j <- fromJSON(content(fc.param,type="text"))
    
    
    if (nn==1){
      ##
      ## for the first datalayer ==================================================================
      ## calculate which elements of the grid vectors are closest
      ##
      dLong <-abs(fc.param.j$GridDefinition$Longitude-targetLong)
      indexLong <- which(dLong == min(dLong))
      dLat <- abs(fc.param.j$GridDefinition$Latitude-targetLat)
      indexLat <- which(dLat == min(dLat))
      
      
      ## 
      ## identify the element which is closest in both the lat and long vectors ===================
      ##
      indexCombined <- intersect(indexLat, indexLong)
      
      print(sprintf("Target location: (%f, %f)", targetLong, targetLat))
      print(sprintf("Closest gridcell: (%f, %f)", fc.param.j$GridDefinition$Longitude[indexCombined], 
                    fc.param.j$GridDefinition$Latitude[indexCombined]))
    }
    
    
    dataValue <- NULL
    for (mm in seq(from=1, to=5)){
      ## 
      ## extract data value @ of closest pixel ====================================================
      ##
      if (nn==1){
        dateValue <- c(dateValue, fc.param.j$Data[[mm]]$Date)
      }
      if (forecastType=="rain"){
        dataValue <- c(dataValue, fc.param.j$Data[[mm]]$Data[indexCombined])
      }
      else {
        dataValue <- c(dataValue, fc.param.j$Data[[mm]]$Data[indexCombined][[1]])
      }
    }
    
    
    ##
    ## create or append collected data  ===========================================================
    ##
    
    if (nn==1){
      resultContainer <- data.frame(date=dateValue)
      resultContainer <- cbind(resultContainer, dataValue)
    }
    else {
      resultContainer <- cbind(resultContainer, dataValue)
    }
    
  }
  
  ##
  ## tidy and return ==============================================================================
  colnames(resultContainer) <- c("date", apiCallVector)
  rownames(resultContainer) <- as.character(seq(from=1, to=nrow(resultContainer)))
  
  return(resultContainer)
  
}









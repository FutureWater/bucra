extract_pixelValues_seasonal_forecast <- function(targetLongVector = c(23.61, 21.66, 20.22), targetLatVector = c(-13.45, -13.95, -13.15),  forecastType="rain"){
  
  
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
  # resultContainer.location <- data.frame()
  resultContainer.final <- list()
  
  
  
  ## 
  ## collect datalayers ===========================================================================
  ##
  
  ## nn counts the data layers
  for(nn in seq(from=1, to=length(apiCallVector))){
    
    print(sprintf("... now processing %s", apiCallVector[nn]))

    query <- sprintf("https://service.weatherimpact.com/api/data/mavo_diami_seasonal/%s?datetime=latest", apiCallVector[nn])
    fc.param <- GET(url=query, add_headers(Authorization="bearer 0a970d7bcb110aaa2c1f5afcea1129ef"))
    fc.param.j <- fromJSON(content(fc.param,type="text"))

    resultContainer.location <- data.frame()
    

    ## pp counts the locations
    for (pp in seq(from=1, to=length(targetLongVector))){

            print(sprintf("... now processing location %0.0f", pp))
      
      
      targetLong <- targetLongVector[pp]
      targetLat <- targetLatVector[pp]
      
      
    
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
    
    
    dateValue <- NULL
    dataValue <- NULL
    
    ## mm counts the time steps
    for (mm in seq(from=1, to=5)){
      ## 
      ## extract data value @ of closest pixel ====================================================
      ##
      if (pp==1){
        dateValue <- c(dateValue, fc.param.j$Data[[mm]]$Date)
      }

      if (forecastType=="rain"){
        dataValue <- c(dataValue, fc.param.j$Data[[mm]]$Data[indexCombined])
      }
      else {
        dataValue <- c(dataValue, fc.param.j$Data[[mm]]$Data[indexCombined][[1]])
      }

    } ## end mm
    
    
    ##
    ## create or append collected data  ===========================================================
    ##
    
    if (pp==1){
      resultContainer.location <- data.frame(date=dateValue)
      resultContainer.location <- cbind(resultContainer.location, dataValue)
    }
    else {
      resultContainer.location <- cbind(resultContainer.location, dataValue)
    }
    
    

  } ## end pp
  
  ##
  ## tidy and return ==============================================================================
  colnames(resultContainer.location) <- c("date", sprintf("location%02.0f", seq(from=1,to=length(targetLongVector))))
    rownames(resultContainer.location) <- as.character(seq(from=1, to=nrow(resultContainer.location)))
    

    resultContainer.final[[nn]] <- list(parameter=apiCallVector[nn],
                                        targetLongVector=targetLongVector,
                                        targetLatVector=targetLatVector,
                                        resultContainer.location=resultContainer.location)
    
    remove(resultContainer.location, dateValue, dataValue)
  } ## end nn
    
  
  return(resultContainer.final)
  
}









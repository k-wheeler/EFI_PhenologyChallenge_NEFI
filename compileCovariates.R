###Script for covariates
##ERA5 data
#' Compile Covariates
#'
#' @param forecastStartDate 
#' @param siteID 
#' @input PhenoForecast
#'
#' @return
#' @export
#'
#' @examples
compileCovariates <- function(forecastStartDate,siteID){
  dates <- seq(as.Date("2021-01-01"),forecastStartDate,"day")
  
  datTairEns <- load_ERA5_Tair_New2(ERA5dataFolder=ERA5dataFolder,endDate=forecastStartDate,stacked=TRUE)
  TairMu <- colMeans(datTairEns)
  TairPrec <- 1/(apply(X=datTairEns,FUN=sd,MARGIN = 2)**2)
  #TairMu <- numeric()
  #TairPrec <- numeric()
  
  ##GEFS for gap filling
  for(d in (length(TairMu)+1):length(dates)){
    #for(d in (24+1):length(dates)){ ##Need to change 
    dte <- dates[d]
    input <- GEFS_Data(dte=dte,siteID=siteID)
    TairMu <- c(TairMu,colMeans(input,na.rm = TRUE)[1])
    TairPrec <- c(TairPrec,1/(apply(X=input,FUN=sd,MARGIN = 2,na.rm=TRUE)**2)[1])
  }
  
  ##GEFS forecast 
  input <- GEFS_Data(dte=forecastStartDate,siteID=siteID)
  TairMu <- c(TairMu,colMeans(input,na.rm = TRUE))
  TairPrec <- c(TairPrec,1/(apply(X=input,FUN=sd,MARGIN = 2,na.rm=TRUE)**2))
  
  output <- list(TairMu=TairMu,TairPrec=TairPrec)
}


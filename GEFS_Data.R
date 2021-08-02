###GEFS Download
# library('devtools')
# library(ncdf4)

#' GEFS Data
#'
#' @param dte 
#' @param siteID 
#' @import devtools
#' @import ncdf4
#'
#' @return
#' @export
#'
#' @examples
GEFS_Data <- function(dte,siteID){
  #source_url('https://github.com/eco4cast/neon4cast/blob/main/R/noaa.R?raw=TRUE')
  time_interval <- "1hr"
  cycle <- "00"
  download_noaa(siteID=siteID,interval=time_interval,date=dte,cycle=cycle,dir='data')
  allMus <- matrix(nrow=30,ncol=35)
  for(i in 1:30){
    if(i<10){
      i <- paste("0",as.character(i),sep="")
    }

    fileName <- dir(path=paste("data/noaa/noaa/NOAAGEFS_",time_interval,"/",siteID,"/",dte,"/",cycle,sep=""),
                              pattern=paste(i,".nc",sep=""))
    fileName <- paste("data/noaa/noaa/NOAAGEFS_",time_interval,"/",siteID,"/",dte,"/",cycle,"/",fileName,sep="")
                          
    nc <- nc_open(fileName)
    time <- as.integer(ncdf4::ncvar_get(nc, "time"))
    tustr <- lubridate::as_datetime(strsplit(ncdf4::ncatt_get(nc, varid = "time", "units")$value
                                             , " ")[[1]][3])
    time <- as.POSIXct.numeric((time*60*60), origin = tustr,tz = "UTC")
    temperature <- ncvar_get(nc,"air_temperature")-273 #Convert from Kelvin to Celcius
    avgTemp <- numeric()
    for(d in seq(forecastStartDate+1,lubridate::date(time[length(time)]),"day")){
      avgTemp <- c(avgTemp,mean(temperature[lubridate::date(time)==d],na.rm = TRUE)) #Averaged over UTC day (not local day)
    }
    allMus[as.numeric(i),1:length(avgTemp)] <- avgTemp
  }
  
  return(allMus)
}
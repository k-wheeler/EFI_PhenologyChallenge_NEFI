#' Download NOAA Weather forecasts for NEON sites from the EFI server
#'  
#' @param siteID vector of 4-character NEON siteIDs
#' @param interval Time interval for the forecast
#' @param date start date for the forecast
#' @param cycle NOAA hour-cycle (first hour of the forecast),
#'  options are "00", "06", "12", "18"; only the "00" forecasts
#'  run 35 days into future.
#' @param dir storage location.  Use tempdir unless you want to keep this 
#' data around on your computer, in which case, `neonstore::neon_dir()` might
#' be a convenient choice.
#' @export
#' @examples 
#' download_noaa("ABBY")
download_noaa <- function(siteID, 
                          interval = "6hr",
                          date = Sys.Date()-2, 
                          cycle = "00", 
                          dir = tempdir()){
  lapply(siteID, download_noaa_, interval, date, cycle, dir)
  invisible(dir)
}
download_noaa_ <- function(siteID, 
                           interval = "6hr",
                           date = Sys.Date()-2, 
                           cycle = "00", 
                           dir = tempdir()){
  
  noaadir <- file.path(dir, "noaa")
  dir.create(noaadir, FALSE, TRUE)
  prefix <- paste("noaa", paste0("NOAAGEFS_", interval), 
                  siteID, date, cycle, sep="/")
  object <- aws.s3::get_bucket("drivers",
                               prefix = prefix,
                               region = "data",
                               base_url = "ecoforecast.org")
  
  #data <- purrr::map_chr(object, ~ .x$Key)
  
  for(i in seq_along(object)){
    aws.s3::save_object(object[[i]], 
                        bucket = "drivers", 
                        file = file.path(noaadir, object[[i]]$Key),
                        region = "data",
                        base_url = "ecoforecast.org")
  }
}

#' Stack downloaded NOAA files
#' 
#' @inheritParams download_noaa
#' @param forecast_date Include only forecasts issued on this date
#' @examples 
#' stack_noaa()
#' @export
stack_noaa <- function(dir = tempdir(), forecast_date = NULL) {
  files <- list.files(file.path(dir, "noaa"), pattern = "[.]nc",
                      recursive = TRUE, full.names = TRUE)
  names(files) <- basename(files)
  if(!is.null(forecast_date)){
    files <- files[stringr::str_detect(files, forecast_date)]
  }
  
  
  out <- purrr::map_dfr(files, function(x){
    tidync::hyper_tibble(tidync::tidync(x))
  }, .id = "file")
  
  ## Add metadata from filename as column...
  out <- tidyr::separate(out, file, "_",
                         into=c("model","interval","siteID",
                                "runStartDate", "runEndDate", "ensemble"))
  
  start_time <- stringr::str_split_fixed(out$runStartDate, pattern = "T", n = 2)
  
  out$time <- lubridate::as_datetime(start_time[, 1]) + lubridate::hours(start_time[, 2]) + lubridate::hours(out$time)
  
  out$ensemble <- stringr::str_split_fixed(out$ensemble, ".nc", 2)[, 1]
  
  return(out)
}
############ and we're ready to go:


# download_noaa("ABBY")
# df <- stack_noaa()




###GEFS Download
# library('devtools')
# library(ncdf4)

#' GEFS Data
#'
#' @param dte 
#' @param siteID 
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
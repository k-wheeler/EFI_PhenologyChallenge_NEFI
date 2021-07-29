library("reticulate")

#' Download ERA5 Temp for Forecast
#'
#' @param end_date 
#' @param siteName 
#' @param lat 
#' @param long 
#' @import reticulate
#'
#' @return
#' @export
#'
#' @examples
downloadERA5Temp <- function(end_date,siteName,lat,long){
  start_date <- as.Date("2021-07-01") ##
  
  setwd("/projectnb/dietzelab/kiwheel/ERA5")
  outfolder <- paste("Data/",siteName,sep="")
  fname <- file.path(outfolder, paste(siteName,"_",start_date,"_",end_date,"_era5AirTemperatureMembers.nc", sep =""))
  
  cdsapi <- reticulate::import("cdsapi")
  cclient <- cdsapi$Client(timeout=600,quiet=FALSE)
  
  variables <- tibble::tribble(
    ~cf_name, ~units, ~api_name, ~ncdf_name,
    "air_temperature", "Kelvin", "2m_temperature", "t2m",
    "air_pressure", "Pa", "surface_pressure", NA_character_
  )
  var <- variables[["api_name"]][[1]]
  area <- rep(round(c(lat, long) * 4) / 4, 2)
  
  do_next <- tryCatch({
    cclient$retrieve(
      "reanalysis-era5-single-levels",
      list(
        variable = var,
        product_type = 'ensemble_members',
        date = paste(start_date, end_date, sep = "/"),
        time = "00/to/23/by/1",
        area = area,
        grid = c(0.25, 0.25),
        format = "netcdf"
      ),
      fname
    )
    FALSE
  }, error = function(e) {
    print("Failed to download variable Mean")
    TRUE
  })
}

library(doParallel)

n.cores <- 8
registerDoParallel(cores=n.cores)

siteData <- read.csv("data/phenologyForecastSites2.csv",header=TRUE)
siteData <- siteData[seq(13,20),]
forecastStartDate <- Sys.Date()

#for(s in 1:nrow(siteData)){
foreach(s=1:nrow(siteData)) %dopar% {
  ##Load Calibration Data: 
  siteName <- as.character(siteData$siteName[s])
  print(siteName)
  lat <- as.numeric(siteData[s,2])
  long <- as.numeric(siteData[s,3])
  
  downloadERA5Temp(end_date=forecastStartDate,siteName=siteName,lat=lat,long=long)

}
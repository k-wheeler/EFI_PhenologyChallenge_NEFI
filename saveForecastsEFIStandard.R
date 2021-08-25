setwd("/projectnb/dietzelab/kiwheel/EFI_PhenologyChallenge_NEFI") 
myPaths <- .libPaths()   # get the paths
myPaths <- c("/projectnb/dietzelab/kiwheel/Rlibrary/4.0.5",myPaths[2]) # switch them
.libPaths(myPaths)  

library('rjags')
library('runjags')

siteData <- read.csv("data/phenologyForecastSites2.csv",header=TRUE)

siteData <- siteData[seq(13,20),]

forecastStartDate <- Sys.Date()
output <- matrix(nrow=0,ncol=6)
for(s in 1:nrow(siteData)){
  ##Load Calibration Data: 
  siteName <- as.character(siteData$siteName[s])
  siteID <- strsplit(siteName,"[.]")[[1]][3]
  print(siteID)
  load(paste0("forecastOutputs/",siteName,"_",forecastStartDate,"_EFI_ForecastChallenge_varBurn.RData"))
  pred.mat <- as.matrix(out.burn$predict)
  
  colnames(output) <- c("time","siteID","forecast","data_assimilation","statistic","gcc_90")

  colNums <- seq(ncol(pred.mat)-34,ncol(pred.mat))
  
  for(d in 1:35){
    output <- rbind(output,c(c(as.character(forecastStartDate + d -1),siteID,1,0,"mean",mean(pred.mat[,colNums[d]]))))
    output <- rbind(output,c(as.character(forecastStartDate + d -1),siteID,1,0,"sd",sd(pred.mat[,colNums[d]])))
  }
}

csvFileName <- paste0("phenology-",forecastStartDate,"-NEFIpheno.csv")
write.csv(file=csvFileName,output,row.names=FALSE)

neon4cast::submit(forecast_file = csvFileName)



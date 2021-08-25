library('rjags')
library('runjags')

siteData <- read.csv("data/phenologyForecastSites2.csv",header=TRUE)

siteData <- siteData[seq(13,20),]

forecastStartDate <- Sys.Date()

for(s in 1:nrow(siteData)){
  ##Load Calibration Data: 
  siteName <- as.character(siteData$siteName[s])
  siteID <- strsplit(siteName,"[.]")[[1]][3]
  load(paste0("forecastOutputs/",siteName,"_",forecastStartDate,"_EFI_ForecastChallenge_varBurn.RData"))
  pred.mat <- as.matrix(out.burn$predict)
  
  output <- matrix(nrow=0,ncol=6)
  colnames(output) <- c("time","sideID","forecast","data_assimilation","statistic","gcc_90")

  colNums <- seq(ncol(pred.mat)-34,ncol(pred.mat))
  
  for(d in 1:35){
    output <- rbind(output,c(c(as.character(forecastStartDate + d -1),siteID,1,0,"mean",mean(pred.mat[,colNums[d]]))))
    output <- rbind(output,c(as.character(forecastStartDate + d -1),siteID,1,0,"sd",sd(pred.mat[,colNums[d]])))
  }
}

csvFileName <- paste0("phenology-",forecastStartDate,"-NEFIpheno.csv")
write.csv(file=csvFileName,output,row.names="FALSE")



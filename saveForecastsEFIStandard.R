library('rjags')
library('runjags')

siteData <- read.csv("data/phenologyForecastSites2.csv",header=TRUE)

siteData <- siteData[seq(13,20),]

forecastStartDate <- Sys.Date()

for(s in 1:nrow(siteData)){
  ##Load Calibration Data: 
  siteName <- as.character(siteData$siteName[s])
  siteID <- strsplit(siteName,"[.]")[[1]][3]
  
  pred.mat <- as.matrix(out.burn$predict)
  
  output <- matrix(nrow=0,ncol=6)
  colnames(output) <- c("time","sideID","forecast","data_assimilation","statistic","gcc_90")

  colNums <- seq(ncol(pred.mat)-34,ncol(pred.mat))
  
  for(d in 1:35){
    output <- rbind(output,c(c(as.character(forecastStartDate + d -1),siteID,1,0,"mean",mean(pred.mat[,colNums[d]]))))
    output <- rbind(output,c(as.character(forecastStartDate + d -1),siteID,1,0,"sd",sd(pred.mat[,colNums[d]])))
  }
}

csvFileName <- paste("phenology-",forecastStartDate,"-NEFIpheno.csv",sep="")
write.csv(file=csvFileName,output,row.names="FALSE")



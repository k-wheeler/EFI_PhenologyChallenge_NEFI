library(phenopix)
library(zoo)


estimateDatesAndScale <- function(p,siteName,year){
  newP <- zoo(na.approx(p))
  #out <- KlostermanFit(newP)
  out <- ElmoreFit(newP)
  output <- c(year,out$fit$sf)
  x=out$fit$predicted
  dts <- PhenoDeriv(x=x,fit=out$fit,uncert=TRUE)
  output <- c(output,dts[4])
  xd <- c(NA,diff(x))
  xd2 <- c(NA,diff(xd))
  #xd3 <- c(NA,diff(xd2))
  sof <- 200+which.min(xd2[200:length(xd2)])-3
  #LastMax <- which.max(xd2[200:length(xd2)])
  #sof <- 200+which.min(xd2[200:(200+LastMax)])-3
  
  output <- c(output,sof)
  PhenoPlot(out,metrics="GCC",main=paste(siteName,year))
  points(newP,col="green",pch=20)
  abline(v=output[c(4,5)],col="red")
  return(output)
}

calculateFits <- function(siteName,URLs,startDate,endDate,year){
  ###Download PhenoCam data and format
  phenoData <- matrix(nrow=0,ncol=32)
  print(URLs[1])
  for(u in 1:length(URLs)){
    phenoDataSub <- download.phenocam(URLs[u])
    phenoData <- rbind(phenoData,phenoDataSub)
  }
  
  ##Order and remove duplicate PC data
  phenoData2 <- phenoData[order(phenoData$date),]
  phenoData3 <- phenoData2[!duplicated(phenoData2$date),]
  phenoData <- phenoData3
  
  phenoData <- phenoData[phenoData$date<=endDate,]
  p.old <- phenoData$gcc_90
  time.old <-  as.Date(phenoData$date)
  
  days <- seq(as.Date(startDate),(as.Date(endDate)),"day")
  p <- rep(NA,length(days))
  
  for(i in 1:length(p.old)){
    p[which(days==time.old[i])] <- p.old[i]
  }
  
  months <- lubridate::month(days)
  years <- lubridate::year(days)
  
  dat2 <- data.frame(dates=days,years=years,months=months,p=p)
  
  dat2 <- dat2[as.numeric(format(dat2$dates,"%j"))%in% seq(1,365),] #Remove final day in leap years
  allDat <- matrix(nrow=0,ncol=5)
  
  for(i in (lubridate::year(as.Date(dat2$dates[1]))):lubridate::year(as.Date(dat2$dates[length(dat2$dates)]))){##I know this includes the forecasted stuff, but it shouldn't really matter because of the JAGS model setup
    print(i)
    subDat <- dat2[lubridate::year(as.Date(dat2$dates))==i,]
    allDat <- rbind(allDat,estimateDatesAndScale(p=subDat$p,siteName,year=i))
  }
  colnames(allDat) <- c("Year","Low","High","PeakDay","FallStartDay")
  return(allDat)
}

siteData <- read.csv("data/phenologyForecastSites2.csv",header=TRUE)
siteData <- siteData[seq(13,20),] #Thinned for NEON sites
dataDirectory <- "data/"
endDate <- as.Date("2020-12-31")

pdf(file="EFIForecastChallenge_AutumnForecastPreviousCurveElmoreFits.pdf",height=5,width=8)
for(i in 1:nrow(siteData)){
  #for(i in 1:5){
  siteName <- as.character(siteData[i,1])
  print(siteName)
  
  startDate <- (as.Date(siteData[i,7]))
  URL <- as.character(siteData$URL[i])
  URL2 <- as.character(siteData$URL2[i])
  URL3 <- as.character(siteData$URL3[i])
  if(nchar(URL2)>0){
    URL <- c(URL,URL2)
    if(nchar(URL3)>0){
      URL <- c(URL,URL3)
    }
  }
  TZ <- as.numeric(siteData[i,6])
  URLs <- URL
  allDat <- calculateFits(siteName,URLs,startDate,endDate)
  save(allDat,file=paste(dataDirectory,siteName,"_phenopixOutputs.RData",sep=""))
}
dev.off()

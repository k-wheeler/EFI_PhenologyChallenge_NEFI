setwd("/projectnb/dietzelab/kiwheel/EFI_PhenologyChallenge_NEFI") 
myPaths <- .libPaths()   # get the paths
myPaths <- c("/projectnb/dietzelab/kiwheel/Rlibrary/4.0.5",myPaths[2]) # switch them
.libPaths(myPaths)  

library(PhenoForecast)
library(PhenologyBayesModeling)
library(rjags)
library(runjags)
library(doParallel)
library(ncdf4)
library(devtools)

##' Calculates cummulative Tair (within one year for one ensemble)
##'
##' @param lat The site latitude
##' @param long The site longitude
##' @param endDate
##' @param calDatesT
##' @param ERA5dataFolder
##' @param TZ_name
##' @param stacked
##' @import xts
##' @import ncdf4
##' @export
load_ERA5_Tair_New2 <- function(lat="",long="",endDate="",calDatesT=TRUE,ERA5dataFolder,TZ_name="America/New_York",stacked=TRUE) {
  if(calDatesT){ #Need to include calibration data
    calFileName <- dir(path=ERA5dataFolder,pattern=paste(endDate,"_era5AirTemperatureMembers",sep=""))
    
    if(length(calFileName)==0){
      downloadERA5Calibration(var="ensemble_members") ##***Need to add this function***
      calFileName <- dir(path=ERA5dataFolder,pattern=paste(endDate,"_era5AirTemperatureMembers",sep=""))
    }
    #print(calFileName)
    #Load data
    ensembleFile <- nc_open(paste(ERA5dataFolder,calFileName,sep=""))
    
    Tairs <- ncvar_get(ensembleFile)-273 #Convert from Kelvin to C
    
    if(stacked){
      Tairs <- Tairs[,1,]
    }
    
    timeHours <- ensembleFile$dim$time$vals #Hours since 1900-01-01 00:00:00.0
    
    ##Convert times to actual times
    times <- as.POSIXct(timeHours*3600, origin = "1900-01-01",tz = "GMT")
    
    attributes(times)$tzone <- TZ_name
    #Daily average
    allDates <- lubridate::date(times)
    dates <- seq(lubridate::date(times[5]),lubridate::date(times[length(times)]),"day")
    TairsDaily <- matrix(nrow=10,ncol=length(dates))
    print(dim(Tairs))
    for(d in 1:length(dates)){
      subTairs <- Tairs[,allDates==dates[d]]
      TairsDaily[,d] <- apply(subTairs,MARGIN=1,mean)
    }
    Tairs <- ncvar_get(ensembleFile)-273 #Convert from Kelvin to C
    
    if(stacked){
      Tairs <- Tairs[,2,]
    }
    
    attributes(times)$tzone <- TZ_name
    #Daily average
    TairsDaily2 <- matrix(nrow=10,ncol=length(dates))
    for(d in 1:length(dates)){
      subTairs <- Tairs[,allDates==dates[d]]
      TairsDaily2[,d] <- apply(subTairs,MARGIN=1,mean)
    }
    TairsDaily[is.na(TairsDaily)] <- TairsDaily2[is.na(TairsDaily)]
    TairsDaily[,is.na(TairsDaily[1,])] <- TairsDaily[,seq(1,ncol(TairsDaily))[is.na(TairsDaily[1,])]-1] 

  }
  ##Current Year (already downloaded)
  ##Will fill in once I get the calibration done
  return(TairsDaily)
  
}

source('compileCovariates.R')
source('GEFS_Data.R')
#source('downloadERA5Temp.R')
#source('saveForecastsEFIStandard.R')

n.cores <- 8

##Put output in EFI standard

#naming: phenology-year-month-day-team_name_ID.csv #First day of the submitted forecast
##statistic column:
  ##mean, sd (predictive distribution), optional: Conf_interv_02.5, Conf_interv_97.5, Pred_interval_02.5,Pred_interval_97.5
siteData <- read.csv("data/phenologyForecastSites2.csv",header=TRUE)

siteData <- siteData[seq(13,20),]
forecastStartDate <- Sys.Date()
startDate <- as.Date("2021-01-01")

baseTemp <- 20
nchain=5
forecastLength <- 35

variableNames <- c("p.PC","p.proc","x","b0","b1","b2","a","CDDtrigger")#,"Dtrigger")

generalModelYear = "
model {
### Data Models for complete years
for(i in 1:n){
p[i] ~ dnorm(x[i],p.PC)
}

#### Process Model
for(i in 2:n){
Tair[i] ~ dnorm(TairMu[i],TairPrec[i])
CDDs[i] <- ifelse(Tair[i]<baseTemp,CDDs[(i-1)]+baseTemp - Tair[i],CDDs[(i-1)])
xmu[i] <- max(min(x[(i-1)] + ifelse(CDDs[i]>CDDtrigger,(b0 + (b1 * x[(i-1)]) + (b2 * x[(i-1)] ** 2)),a),x[1]),0)
x[i] ~ dnorm(xmu[i],p.proc)
}

#### Priors
x[1] ~ dbeta(x1.a,x1.b)
CDDs[1] <- 0
p.PC ~ dgamma(s1.PC,s2.PC)
p.proc ~ dgamma(s1.proc,s2.proc)
CDDtrigger ~ dnorm(CDDtrigger_mu,CDDtrigger_prec)

a ~ dnorm(a_mu,a_prec)
b0 ~ dnorm(b0_mu,b0_prec)
b1 ~ dnorm(b1_mu,b1_prec) 
b2 ~ dnorm(b2_mu,b2_prec)
}
"

for(s in 1:nrow(siteData)){
  ##Load Calibration Data: 
  siteName <- as.character(siteData$siteName[s])
  siteID <- strsplit(siteName,"[.]")[[1]][3]
  print(siteID)
  
  calFileName <- paste0("partial2_",siteName,"_EFI_ForecastChallenge_calibration_varBurn.RData")
  outputFileName <- paste0("forecastOutputs/",siteName,"_",forecastStartDate,"_EFI_ForecastChallenge_varBurn.RData")
  ERA5dataFolder <- paste0("/projectnb/dietzelab/kiwheel/ERA5/Data/",siteName,"/")
  
  lat <- as.numeric(siteData[s,2])
  long <- as.numeric(siteData[s,3])
  #startDate <- (as.Date(siteData[s,7]))
  startDate <- as.Date("2021-01-01")
  URL <- as.character(siteData$URL[s])
  URL2 <- as.character(siteData$URL2[s])
  URL3 <- as.character(siteData$URL3[s])
  if(nchar(URL2)>0){
    URL <- c(URL,URL2)
    if(nchar(URL3)>0){
      URL <- c(URL,URL3)
    }
  }
  TZ <- as.numeric(siteData[s,6])
  URLs <- URL
  
  phenoData <- matrix(nrow=0,ncol=32)
  print(URLs)
  for(u in 1:length(URLs)){
    phenoDataSub <- download.phenocam(URLs[u])
    phenoData <- rbind(phenoData,phenoDataSub)
  }

  ##Order and remove duplicate PC data
  phenoData2 <- phenoData[order(phenoData$date),]
  phenoData3 <- phenoData2[!duplicated(phenoData2$date),]
  phenoData <- phenoData3
  
  phenoData <- phenoData[phenoData$date<forecastStartDate,]
  p.old <- phenoData$gcc_90
  time.old <-  as.Date(phenoData$date)
  
  days <- seq(as.Date(startDate),(as.Date(forecastStartDate))+forecastLength,"day")
  p <- rep(NA,length(days))
  
  for(i in 1:length(p.old)){
    p[which(days==time.old[i])] <- p.old[i]
  }
  
  months <- lubridate::month(days)
  years <- lubridate::year(days)
  
  dat2 <- data.frame(dates=days,years=years,months=months,p=p)
  datTairs <- compileCovariates(forecastStartDate=forecastStartDate,siteID=siteID)
  
  print("Finished loading met")
  
  ICsdat <- dat2[as.numeric(format(dat2$dates,"%j"))%in% seq(203,212),]
  dat2$TairMu <- datTairs$TairMu
  dat2$TairPrec <- datTairs$TairPrec
  dat2 <- dat2[as.numeric(format(dat2$dates,"%j"))%in% seq(213,365),]
  
  #nrowNum <- as.numeric(format(forecastStartDate,'%j')) -212 + forecastLength
  nrowNum <- length(dat2$p)
  
  dataFinal <- list(p=dat2$p,TairMu=dat2$TairMu,TairPrec=dat2$TairPrec,baseTemp=baseTemp)
  dataFinal$n <- nrowNum
  ICs <- ICsdat$p
  mu <- mean(ICs,na.rm=TRUE)
  vr <- var(ICs,na.rm = TRUE)
  x1a <- (mu**2-mu**3-mu*vr)/(vr)
  x1b <- (mu-2*mu**2+mu**3-vr+mu*vr)/(vr)
  dataFinal$x1.a <- x1a
  dataFinal$x1.b <- x1b
  
  load(calFileName)
  out.burn <- partialOutput
  out.mat <- data.frame(as.matrix(out.burn$params))
  dataFinal$CDDtrigger_mu <- mean(out.mat$CDDtrigger)
  dataFinal$CDDtrigger_prec <- 1/var(out.mat$CDDtrigger)*100
  dataFinal$a_mu <- mean(out.mat$a)
  dataFinal$a_prec <- 1/var(out.mat$a)
  dataFinal$b0_mu <- mean(out.mat$b0)
  dataFinal$b0_prec <- 1/var(out.mat$b0)
  dataFinal$b1_mu <- mean(out.mat$b1)
  dataFinal$b1_prec <- 1/var(out.mat$b1)
  dataFinal$b2_mu <- mean(out.mat$b2)
  dataFinal$b2_prec <- 1/var(out.mat$b2)
  
  dataFinal$s1.PC <- mean(out.mat$p.PC)**2/var(out.mat$p.PC)
  dataFinal$s2.PC <- mean(out.mat$p.PC)/var(out.mat$p.PC)
  dataFinal$s1.proc <- mean(out.mat$p.proc)**2/var(out.mat$p.proc)
  dataFinal$s2.proc <- mean(out.mat$p.proc)/var(out.mat$p.proc)

  inits <- list()
  
  for(i in 1:nchain){
    inits[[i]] <- list(CDDtrigger = rnorm(1,mean(out.mat$CDDtrigger),sd(out.mat$CDDtrigger)),
                       a=rnorm(1,mean(out.mat$a),sd(out.mat$a)),
                       b0=rnorm(1,mean(out.mat$b0),sd(out.mat$b0)),
                       b1=rnorm(1,mean(out.mat$b1),sd(out.mat$b1)),
                       b2=rnorm(1,mean(out.mat$b2),sd(out.mat$b2)))
  }
  
  j.model   <- jags.model(file = textConnection(generalModelYear),
                          data = dataFinal,
                          inits = inits,
                          n.chains = nchain,
                          n.adapt = 1000)
  
  out.burn <- runForecastIter(j.model=j.model,variableNames=variableNames,
                              baseNum = 10000,iterSize = 10000,effSize = 5000,partialFile = paste("partial_",outputFileName,sep=""),
                              maxIter = 200000)
  out.mat <- as.matrix(out.burn$params)
  thinAmount <- round(nrow(out.mat)/1000,digits=0)
  out.burn2 <- list()
  out.burn2$params <- window(out.burn$params,thin=thinAmount)
  out.burn2$predict <- window(out.burn$predict,thin=thinAmount)
  out.burn <- out.burn2
  save(out.burn,file = outputFileName)
}


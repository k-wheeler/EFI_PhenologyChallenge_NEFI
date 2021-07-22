library(PhenoForecast)
library(PhenologyBayesModeling)
library(rjags)
library(runjags)
#siteData <- read.csv("PhenologyForecastData/phenologyForecastSites.csv",header=TRUE)
siteData <- read.csv("data/phenologyForecastSites2.csv",header=TRUE)
dataDirectory <- "data/"
endDate <- as.Date("2020-12-31")

baseTemp <- 20
nchain=5

siteData <- siteData[seq(13,20),] #Thinned for NEON sites

for(i in 1:nrow(siteData)){
  
  siteName <- as.character(siteData[i,1])
  print(siteName)
  ERA5dataFolder <- paste("/projectnb/dietzelab/kiwheel/ERA5/Data/",siteName,"/",sep="")
  outputFileName <- paste(siteName,"_EFI_ForecastChallenge_calibration_varBurn.RData",sep="")
  load(file=paste(dataDirectory,siteName,"_phenopixOutputs.RData",sep=""))
  fittedDat=allDat
  
  lat <- as.numeric(siteData[i,2])
  long <- as.numeric(siteData[i,3])
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
  
  ###Download PhenoCam data and format
  phenoData <- matrix(nrow=0,ncol=32)
  print(URLs[1])
  for(u in 1:length(URLs)){
    phenoDataSub <- download.phenocam(URLs[u])
    phenoData <- rbind(phenoData,phenoDataSub)
  }
  #print(phenoData)
  ##Order and remove duplicate PC data
  phenoData2 <- phenoData[order(phenoData$date),]
  phenoData3 <- phenoData2[!duplicated(phenoData2$date),]
  phenoData <- phenoData3
  
  phenoData <- phenoData[phenoData$date<endDate,]
  p.old <- phenoData$gcc_90
  time.old <-  as.Date(phenoData$date)
  #print("passed time.old")
  
  days <- seq(as.Date(startDate),(as.Date(endDate)),"day")
  #print("past days")
  p <- rep(NA,length(days))
  
  for(i in 1:length(p.old)){
    p[which(days==time.old[i])] <- p.old[i]
  }
  
  months <- lubridate::month(days)
  years <- lubridate::year(days)
  
  dat2 <- data.frame(dates=days,years=years,months=months,p=p)
  datTairEns <- load_ERA5_Tair_New(ERA5dataFolder=ERA5dataFolder,endDate=endDate)
  print(head(datTairEns))
  print("Finished loading ERA5")
  TairMu <- apply(X=datTairEns,MARGIN=2,FUN=mean)
  TairPrec <- 1/apply(X=datTairEns,MARGIN=2,FUN=var)
  dat2$TairMu <- TairMu + (0-min(TairMu)) ##Done to make sure all temperatures are >= 0 
  baseTempOrig <- baseTemp
  baseTemp <- baseTemp + (0-min(TairMu)) 
  dat2$TairPrec<- TairPrec
  
  ICsdat <- dat2[as.numeric(format(dat2$dates,"%j"))%in% seq(203,212),]
  dat2 <- dat2[as.numeric(format(dat2$dates,"%j"))%in% seq(213,365),]
  
  nrowNum <- 365-212
  p <- matrix(nrow=nrowNum,ncol=0)
  TairMu <- matrix(nrow=nrowNum,ncol=0)
  D <- matrix(nrow=nrowNum,ncol=0)
  ICs <- matrix(nrow=10,ncol=0)
  TairPrec <- matrix(nrow=nrowNum,ncol=0)
  valNum <- 0
  days2 <- matrix(nrow=nrowNum,ncol=0)
  
  finalYrs <- numeric()
  sofs <- numeric()
  
  for(i in (lubridate::year(as.Date(dat2$dates[1]))):lubridate::year(as.Date(dat2$dates[length(dat2$dates)]))){##I know this includes the forecasted stuff, but it shouldn't really matter because of the JAGS model setup
    #print(i)
    subDat <- dat2[lubridate::year(as.Date(dat2$dates))==i,]
    valNum <- valNum + 1
    newCol <- subDat$p
    p <- cbind(p,newCol)
    ICs <- cbind(ICs,ICsdat[lubridate::year(as.Date(ICsdat$dates))==i,]$p)
    days2 <- cbind(days2,as.Date(subDat$dates))
    finalYrs <- c(finalYrs,i)
    TairMu <- cbind(TairMu,subDat$TairMu)
    TairPrec <- cbind(TairPrec,subDat$TairPrec)
    sofs <- c(sofs,(fittedDat[valNum,'FallStartDay']-212))
  }
  
  dataFinal <- list(p=p,years=finalYrs,sofMean=mean(sofs))
  dataFinal$n <- nrowNum
  dataFinal$N <- ncol(dataFinal$p)
  dataFinal$CDDtrigger.lower <- 0
  dataFinal$CDDtrigger.upper <- 500
  dataFinal$s1.PC <- 1.56#1262.626 ## Very roughly based off of what I think are reasonable and uninformed priors
  dataFinal$s2.PC <- 0.016
  dataFinal$s1.proc <- 1.56
  dataFinal$s2.proc <- 0.016
  
  x1a <- numeric()
  x1b <- numeric()
  for(yr in 1:dataFinal$N){
    mu <- mean(ICs[,yr],na.rm=TRUE)
    vr <- var(ICs[,yr],na.rm = TRUE)
    x1a <- c(x1a,(mu**2-mu**3-mu*vr)/(vr))
    x1b <- c(x1b,(mu-2*mu**2+mu**3-vr+mu*vr)/(vr))
  }
  
  dataFinal$x1.a <- x1a
  dataFinal$x1.b <- x1b
  
  dataFinal$TairMu <- TairMu
  dataFinal$TairPrec <- TairPrec
  dataFinal$baseTemp <- baseTemp
  
  inits <- list()
  
  dataFinal$b0_lower <- -1
  dataFinal$b0_upper <- 1
  dataFinal$b1_lower <- -1
  dataFinal$b1_upper <- 1
  dataFinal$b2_lower <- -1
  dataFinal$b2_upper <- 1
  dataFinal$a_upper <- 0
  dataFinal$a_lower <- -0.01
  dataFinal$CDDtrigger.lower <- 0
  dataFinal$CDDtrigger.upper <- 500
  transitionCDD <- numeric()
  agings <- numeric()
  ##Determine Inits:
  
  for(yr in 1:dataFinal$N){
    CDD <- 0
    for(i in 1:sofs[yr]){
      print(dim(dataFinal$baseTemp))
      print(dataFinal$baseTemp)
      if(dataFinal$TairMu[i,yr]<dataFinal$baseTemp){
        CDD <- CDD + (dataFinal$baseTemp - dataFinal$TairMu[i,yr])
      }
    }
    transitionCDD <- c(transitionCDD,CDD)
    agings <- c(agings,
                as.numeric(lm(dataFinal$p[,yr]~seq(1,length(dataFinal$p[,yr])))$coefficients[2]))
  }
  
  for(i in 1:nchain){
    inits[[i]] <- list(CDDtrigger = rnorm(1,mean(transitionCDD),10),
                       a=rnorm(1,mean(agings),0.0005),
                       b0=rnorm(1,-0.25,0.05),
                       b1=rnorm(1,0.6,0.08),
                       b2=rnorm(1,-0.5,0.1))
  }
  
  variableNames <- c("p.PC","p.proc","x","b0","b1","b2","a","CDDtrigger")#,"Dtrigger")
  
  generalModel = "
    model {
  ### Data Models for complete years
  for(yr in 1:(N)){
  for(i in 1:n){
  p[i,yr] ~ dnorm(x[i,yr],p.PC)
  }
  }
  
  #### Process Model
  for(yr in 1:(N)){
  for(i in 2:n){
  Tair[i,yr] ~ dnorm(TairMu[i,yr],TairPrec[i,yr])
  CDDs[i,yr] <- ifelse(TairMu[i,yr]<baseTemp,CDDs[(i-1),yr]+baseTemp - Tair[i,yr],CDDs[(i-1),yr])
  xmu[i,yr] <- x[(i-1),yr] + ifelse(CDDs[i,yr]>CDDtrigger,(b0 + (b1 * x[(i-1),yr]) + (b2 * x[(i-1),yr] ** 2)),a)
  x[i,yr] ~ dnorm(xmu[i,yr],p.proc) T(0,0.999)
  }
  }
  
  #### Priors
  for(yr in 1:N){ ##Initial Conditions
  x[1,yr] ~ dbeta(x1.a[yr],x1.b[yr])
  CDDs[1,yr] <- 0
  }
  p.PC ~ dgamma(s1.PC,s2.PC)
  p.proc ~ dgamma(s1.proc,s2.proc)
  CDDtrigger ~ dunif(CDDtrigger.lower,CDDtrigger.upper)
  
  a ~ dunif(a_lower,a_upper)
  b0 ~ dunif(b0_lower,b0_upper)
  b2 ~ dunif(b2_lower,b2_upper)
  b1 ~ dunif(b1_lower,b1_upper) 
  }
  "
  save(dataFinal,file=paste(siteName,"_EFIChallengeCalibration_dataFinal.RData",sep=""))
  
  j.model   <- jags.model(file = textConnection(generalModel),
                          data = dataFinal,
                          n.chains = nchain,
                          n.adapt = 1500)
  
  out.burn <- runForecastIter(j.model=j.model,variableNames=variableNames,
                              baseNum = 15000,iterSize = 5000,effSize = 5000,partialFile = paste("partial_",outputFileName,sep=""))
  ##Thin the data:
  out.mat <- as.matrix(out.burn$params)
  thinAmount <- round(nrow(out.mat)/5000,digits=0)
  out.burn2 <- list()
  out.burn2$params <- window(out.burn$params,thin=thinAmount)
  out.burn2$predict <- window(out.burn$predict,thin=thinAmount)
  out.burn <- out.burn2
  save(out.burn,file = outputFileName)
  
}


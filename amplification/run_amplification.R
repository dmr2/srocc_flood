#!/opt/local/bin/R
# run_amplification.R
# 
# last updated: Tue Apr  3 14:11:05 EDT 2018
# 
#  Calculate flood amplification factors for
#  given return periods and selected years for
#  multiple global mean surface temperature
#  scenarios.
#
#  Inputs:
#   - a file with GPD parameters for each site
#   - 5/50/95th percentile sea-level rise projections 
#
#  Outputs:
#   - writes amplification factors to a table

rm(list=ls(all=TRUE))

mainDir = "/Users/dmr/Dropbox/SROCC/srocc_flood/amplification" # current directory
setwd(mainDir) 

source("routines/GPDsample.R")
source("routines/GPDNExceed.R")
source("routines/calcAF.R")

# Location of local sea level rise projections
dir = "../"

type = "srocc" # label for files
targ_years = c(2050,2100) # Years to generate AFs for...

scenarios = c("RCP2.6","RCP4.5","RCP8.5") # e.g., RCPs or GMST targets
slab = c("rcp26","rcp45","rcp85")

# generate AFs for these return periods, e.g., 10-, 20-, 50-, 100-, 500-yr
rp_list = c(10,20,50,100,500)

###
#targ_years = c(2100) # Years to generate AFs for...

#scenarios = c("RCP8.5") # e.g., RCPs or GMST targets
#slab = c("rcp85")

# generate AFs for these return periods, e.g., 10-, 20-, 50-, 100-, 500-yr
#rp_list = c(100) 
###

# Open GPD parameters and historical flood data
dat <- read.csv("GPDfits_uhawaii_projectLSLfmt17.tsv",header=T,sep="\t")
sites <- dat$Site # tide gauge sites

AF <- array(NaN,dim=c(length(sites),length(targ_years),
                      length(scenarios),length(rp_list),3))

histFloodHeight <- array(NaN,dim=c(length(sites),length(targ_years),
                                  length(scenarios),length(rp_list)))
FutureFreq <- array(NaN,dim=c(length(sites),length(targ_years),
                              length(scenarios),length(rp_list),3))
slr <- array(NaN,dim=c(length(sites),length(scenarios),length(targ_years),3))

for(j in 1:length(sites)){  

  site <- dat$Site[j] # tide gauge site name
  scale <- dat$Scale50[j] # median scale parameter
  shape <- dat$Shape50[j] # median shape parameter
  UHid <- dat$UHAWAII_ID[j] # U Hawaii site identifier
  threshold <- dat$Q99[j] # GPD threshold
  lambda <- dat$Lambda[j] # mean Poisson arrival rate of threshold
  shapescaleV <- dat$Vscaleshape[j] # covariance of GPD scale and shape parameters
  shapeV <- dat$Vshape[j] # variance of shape parameter
  scaleV <- dat$Vscale[j] # variance of scale parameter
  psmslid <- dat$PSMSL_ID[j] # Tide gauge ID
  basin <- dat$Basin[j] # Ocean
  
 # Account for GPD parameter uncertainty by making draws from a
 # bivariate normal distribution using Latin hypercube sampling
  GPD <- GPDsample(1000, scale, shape, shapeV, scaleV, shapescaleV)

  z <- seq(0,10,.01) # some flood heights (meters above tide gauge MHHW)
  
  # Expected historical flood height curve (No SLR) (GPD uncertainty)
  histCurveSamps <- matrix(NaN, nrow=length(GPD[,2]), ncol=length(z))
  for(iii in 1:length(GPD[,2]) ){
    histCurveSamps[iii,] <- GPDNExceed(z-threshold,lambda,-threshold,GPD[iii,2],GPD[iii,1])
  }
  
  Ne_hist <- apply(histCurveSamps,2,mean,na.rm=T) 
  
  for( s in 1:length(scenarios) ){
    
    for( t in 1:length(targ_years) ){
          
      # Get sea level rise projections
      fil <- Sys.glob(paste( dir,"/SROCC_",slab[s],"_",targ_years[t],".txt",sep=""))
      if (length(fil)==0){
         print(paste("Skipping... No SLR projections for site: ",site,sep=""))
         next
      }
      
      x = read.table( paste(fil, sep=""), skip=0, sep=" ", header=TRUE)
    
      print(paste("[",dat$Site[j],"] Calculating AFs for year: ",targ_years[t]," Scenario: ", scenarios[s]))
      
      slr[j,s,t,] <- as.numeric(x[x$Station==site,5:7])/100 # centimeters to meters
      
      for (r in 1:length(rp_list)){
        
       print(paste("  ",rp_list[r]))
       q <- 1/rp_list[r]
        
        for (jj in 1:3){ # each SLR percenilte (5th, 50th, 95th)
          af <- get_AF(rp_list[r],z,Ne_hist,histCurveSamps,slr[j,s,t,jj])
          AF[j,t,s,r,jj] <- round(mean(af$AF,2),2)
          FutureFreq[j,t,s,r,jj] <- round(AF[j,t,s,r,jj] * 1/rp_list[r],2)
          
          # Mask AF if dominated by tidal events (i.e., every other day)
          limit <- (365.25/2)/(1/rp_list[r])
          if( AF[j,t,s,r,jj] >= limit ){
            AF[j,t,s,r,jj] = NA
            FutureFreq[j,t,s,r,jj] <- NA
          }
          
        }
      
       histFloodHeight[j,t,s,r] <- af$HistHeight
      }
      
    } # each target year
  } # each scenario
  }# each tide gauge

# Write all AFs for all sites to a table

dir.create(file.path(mainDir, "output"))

for( r in 1:length(rp_list) ){
for( s in 1:length(scenarios) ){
  for( t in 1:length(targ_years) ){
    
    # put the data we want into data frames for writing to disk
    
    df <- data.frame(Site=dat$Site,Country=dat$Region,Region=dat$Basin,
                     ID=dat$PSMSL_ID,UHAWAII_ID=dat$UHAWAII_ID,
                     Lat=dat$Lat,Lon=dat$Lon,
                     Scenario=rep(slab[s],length(dat$Site)),
                     Year=rep(targ_years[t],length(dat$Site)),
                     SLR_5_m=slr[,s,t,1],SLR_50_m=slr[,s,t,2],SLR_95=slr[,s,t,3],
                     GPD_threshold_m=dat$Q99,
                     Hist_Freq=rep(round(1/rp_list[r],3),length(dat$Site)),
                     Hist_Height_m=histFloodHeight[,t,s,r],Fut_Freq_5=FutureFreq[,t,s,r,1],
                     Fut_Freq_50=FutureFreq[,t,s,r,2],Fut_Freq_95=FutureFreq[,t,s,r,3],
                     AF_5=AF[,t,s,r,1],AF_50=AF[,t,s,r,2],AF_95=AF[,t,s,r,3])
                  
    outf <- paste("output/af_",type,"_",slab[s],"_rp",rp_list[r],"_",targ_years[t],".csv",sep="")
    write.table(df,outf,sep=",",row.names = FALSE)
  }
}
}
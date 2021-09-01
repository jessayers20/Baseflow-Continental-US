###################################################################
#########           MK trends for the entire USA           ########
#########                                                  ########
###################################################################


# This script outputs the data from the model's predictive probability distribution
# the goal of these tables is to use to recreate the MK trends
# however, could also plot this data (use for centiles plots)


rm(list=ls(all=TRUE)) ; cat("\014")

library(ggplot2); library(gamlss); library(data.table); library(raster); library(modifiedmk)
library(data.table); library(sp); library(rgeos); library(rgdal); library(maptools)
library(maps); library(plyr); library(gridExtra); library(tidyr); library(corrplot)
library(gtools); library(tools); library(doParallel);  library(foreach); library(Kendall)

#mon.con = read.csv("/Users/jesayers/USGS_Streamflow/Data/Stations/USA_Monthly_Contribution.csv")
mon.con = read.csv("T:/USGS_Streamflow/Data/Stations/USA_Monthly_Contribution.csv")

ref.gages <- read.csv("T:/Bf-Climate (US)/Data/Reference Sites Data/USA_BF_refgages_2021.csv",header=T)
sites <- ref.gages$site
sites <- sites[-1029]

#stations = bf.sites = list.files("T:/Historic Baseflow Data/Climate")

trend.results <- list()

for(s in 1:length(sites)){
  #s=1000
  station = as.numeric(sites[s])
  #station = 1010000
 
  #temp = read.csv(paste0("T:/USGS_Streamflow/Data/Discharge/No_Filter_Monthly_Baseflow/Monthly_BF_",station,".csv"),sep=",", header = TRUE)
  temp = read.csv(paste0("T:/Data/Baseflow_Monthly_Onepar/",station,"_monthly_BF.csv"),header=T)
  temp = temp[,c("Year","Month","Baseflow")]
  
  temp$Month = factor(temp$Month,levels =c(1,2,3,4,5,6,7,8,9,10,11,12),
                         labels=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))
  
  #Subset for which months contribute more
  
  pathS = paste0("T:/Bf-Climate (US)/Data/Model_Selection/Eck_Prob_Dist_1940/",station,"")
  files = list.files(paste0("T:/Bf-Climate (US)/Data/Model_Selection/Eck_Prob_Dist_1940/",station,"/"))
  
  if(length(files)==0) {
    next
  }
  
  GET <- function(number){ # skips the first 5 lines of each csv file and tab deleniates them
    mods=read.csv(paste0(pathS,"/",files[number]),header=T, sep=",")
    m = substr(unlist(strsplit(files[number],"/"))[[1]],1,30) 
    mt = lapply(m, function(x) substring(x,nchar(m)-6))
    mon = unlist(mt)
    month = substr(mon,1,nchar(mon)-4)
    mods$month = paste0(month)
    return(mods)
  }
  
  # Data frames separated by year (in list), so need to merge them together
  each = lapply(1:length(files), function(i){ GET(i)  }) #applies function to length of list
  all.betvals = rbindlist(each)
  mnths = paste(unique(all.betvals$month))
  
  con.t <- mon.con[which(mon.con$site == station),]
  
  temp = temp[which(temp$Month %in% mnths),]
  tmon = spread(temp,key="Month",value="Baseflow")
  all = na.omit(tmon)
  
  library(quantreg);library(modifiedmk); library(Kendall)
  
  all = all[which(all$Year >= 1989 & all$Year <=2019),]
  #Subset for predictions or observations

  library(quantreg);library(modifiedmk); library(Kendall)
  allT = all[,-1]
  #allT <- allT %>% dplyr::select_if( ~ length(unique(.)) > 1)
  runmo = colnames(allT)
  if(length(runmo>0) & nrow(all)>25){
    pw<-apply(allT[,1:length(allT)],2,function(x){ mkttest(as.vector(x))}) #because zero values in columns
    pw = as.data.frame(pw)#grabs tau from output results
    zval = sapply(pw,"[[",1)
    tau = sapply(pw,"[[",6)
    sl= sapply(pw,"[[",2)
    pval = sapply(pw,"[[",5)
    pos <- ifelse(sl >= 0, "Positive", "Negative")
    five <- ifelse(abs(pval) <= 0.05, "Significant", "Nonsignificant")
    res1 = data.frame(zval,pval,tau,sl,pos,five)
    res1$month = row.names(res1)
    res1$site = paste(station)
 write.csv(res1,file=paste0("T:/Bf-Climate (US)/Data/Trends/Onepar_Observed/BF_Trends_1989_",station,".csv"), row.names=F)
  }
    }
  
  
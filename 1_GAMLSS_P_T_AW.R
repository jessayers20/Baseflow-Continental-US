######################################################
####### GAMLSS Selection and Write CSV Results #######
#######      P, T, and AW, for entire US       #######
######################################################

# This Script Runs the selection for the GAMLSS model with antecedent wetness (sum of prev 3 months)
# output is an excel file with all the functions, correlations, beta values, and p-vals

# Load libraries, data ----------------------------------------------------

rm(list=ls(all=TRUE)) ; cat("\014")

library(sp); library(rgdal); library(gamlss); library(sp); library(tools); library(data.table); library(ggplot2); 
library(tidyr); library(foreach); library(doParallel); library(gamlss)

months = c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")

#sites = read.csv("T:/USGS_Streamflow/Contiguous US Baseflow/Sites_Information.csv")
#stations = sites$site_no
sites = list.files("/Users/jesayers/Bf-Climate (US)/Data/Readied_CSVs_for_Analysis/Eckhardt_2/")
#sites = list.files("T:/Bf-Climate (US)/Data/Readied_CSVs_for_Analysis/Eckhardt_2/")


mon.con = read.csv("/Users/jesayers/USGS_Streamflow/Data/Stations/USA_Monthly_Contribution.csv")
#mon.con = read.csv("T:/USGS_Streamflow/Data/Stations/USA_Monthly_Contribution.csv")
mon.con = mon.con[,-14]
mon.con=gather(mon.con, key="month",value="cont",-site)

mon.con = mon.con[which(mon.con$site %in% sites),]
length(unique(mon.con$site))

#cl = makeCluster(12)
#registerDoParallel(cl)  

# Begin loop through months and stations ----------------------------------------
month_list = list()
for(m in 1:length(months))
{
  #m=4
  month = months[m]
  
  #Call stations list for each month with >5% monthly contribution (unique for each)
  cont = mon.con[which(mon.con$month == month),]
  stations = cont[which(cont$cont > 5),]$site
  #stations = stations[1:20]

  station_summary = data.frame(matrix(NA,length(stations),12))
  colnames(station_summary)=c("Station", "Mu", "Function", "R", "RS", "Intercetp","beta_1","beta_2","beta_3","pValue_1","pValue_2","pValue_3")

  #stat_list = list()  

for(s in 1:length(stations)){
   #stat_list = foreach (station = stations) %dopar%
 tryCatch({     

    #station = 3237500
    #s=100
    station = stations[s]
    library(magrittr); library(tidyr)
    #s=match(station,stations)
    print(station)

    #clim = read.csv(paste0("/Users/jesayers/Bf-Climate (US)/Data/Readied_CSVs_for_Analysis/Eckhardt/",station,"/Readied_",station,"_",month,".csv"), sep=",",header = T)
    clim = read.csv(paste0("T:/Bf-Climate (US)/Data/Readied_CSVs_for_Analysis/Eckhardt/",station,"/Readied_",station,"_",month,".csv"), sep=",",header = T)
    
    clim = clim[,2:6]
    colnames(clim) = c("Year","Bf","Pr","Te","Aw")
    
    temp = read.csv(paste0("/Users/jesayers/Data/Baseflow_Monthly_Onepar/",station,"_monthly_BF.csv"),sep=",", header = TRUE)
    temp$Month = factor(temp$Month, levels = c(1,2,3,4,5,6,7,8,9,10,11,12),
                        labels = c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))
    temp = temp[which(temp$Month == month),]
    library(dplyr)
    c2 = inner_join(temp,clim, by = c("Year"))
    clim=c2[,c(1,3,7:9)]
    colnames(clim) = c("Year","Bf","Pr","Te","Aw")
    
    table = clim[which(clim$Year >=1940),]
    
    #Years for computation; Y those which do not have observations; removal must be after extract years vector
    #ns = min(which(!is.na(table$Bf)))
    #last = nrow(table)
    #table = table[ns:last,]
    
    table$Year = as.numeric(table$Year)
    minyear<-min(table$Year);maxyear<-max(table$Year)
    years = minyear:maxyear
    
    table[is.na(table[,3]),2] = NA
    table[is.na(table[,4]),2] = NA
    table[is.na(table[,5]),2] = NA
    table[is.na(table[,2]),3] = NA
    table[is.na(table[,2]),4] = NA
    table[is.na(table[,2]),5] = NA
      
    if(any(is.na(table$Bf))){years <- years[-which(is.na(table$Bf))]}
    if(any(is.na(table$Bf))){table <- table[-which(is.na(table$Bf)),]}
    
    years <- table$Year
      
    Bf  = table$Bf
    Pr = table$Pr
    Aw = table$Aw
    Te = table$Te
      
    # Create Dataset to input into the model after NAs have been removed from original 'table'
    mod.sel <- data.frame(Bf,Pr,Te,Aw)
    
    drivers = colnames(mod.sel[,-1]) #name of drivers used in model selection
    
    test.year = years[which(years>1985)]
    
    if(any(Bf == 0)){next}
    
    if(length(test.year)>=30 & sum(table$Bf >0)){
      
    library(gamlss)
    mod1 <- gamlss(mod.sel$Bf ~ ., data = mod.sel, family = GA(mu.link="log"))
    
    #Model to selection for Mu parameter
    mod2<-stepGAICAll.B(mod1,what="mu", k=log(length(Bf)))
    mod.mu <- summary(mod2)
    
    mu.pars = substring(mod2$call[2],14)

    sregr = data.frame(summary(mod2))
    
    SBC = mod2$sbc
    
    betas = sregr$Estimate[2:length(sregr$Estimate)] # Take all beta coefficients (excluding intercept)
    pValues = sregr[,4][2:length(sregr[,4])] # Take all pValues 
    
    source("/Users/jesayers/BF_Factors/centiles1.R") # I modified the function to create a dataframe containing the values of centiles
    #source("T:/BF_Factors/centiles1.R") 
    
    a <- centiles1(mod2, xvar = years, cent = c(5,25,50,75,95),save=T, plot=F)
    a$formula =  paste0(substring(mod2$call[2],14))

    head(a);
    head(table)
    #dir.create(paste0("/Users/jesayers/Bf-Climate (US)/Data/Model_Selection/Log_Prob_Dist_1940/",station,""))
    write.csv(a,paste0("/Users/jesayers/Bf-Climate (US)/Data/Model_Selection/Log_Prob_Dist_1940/",station,"/Centiles_",station,"_",month,".csv"),row.names = F)
    
    intercept = sregr$Estimate[1]
    R = cor.test(a$Obs,a$n_50)$est 
    RS = cor.test(a$Obs,a$n_50, method = "spearman")$est
   
    ##### write table #####
    station_summary[s,1] = station #write results to a table
    station_summary[s,2] = paste0(mu.pars)
    station_summary[s,3] = deparse(mod2$call$formula)
    station_summary[s,4] = R
    station_summary[s,5] = RS
    station_summary[s,6] = intercept
    station_summary[s,7] = betas[1]
    station_summary[s,8] = betas[2]
    station_summary[s,9] = betas[3]
    station_summary[s,10] = pValues[1]
    station_summary[s,11] = pValues[2]
    station_summary[s,12] = pValues[3]
      
    }
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
  #a = table_empty
  write.csv(station_summary,paste0("/Users/jesayers/Bf-Climate (US)/Data/Model_Selection/Predictors_Formula_1940/Model_results_",month,".csv"),row.names = F)

}



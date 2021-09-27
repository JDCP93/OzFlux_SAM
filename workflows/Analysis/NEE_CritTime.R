rm(list=ls())

Var = "Tair"
#Transect = c("NATT","SAWS")
#Transect = c("NATT")
Transect = c("SAWS")

Sites = c("AU-ASM","AU-Cpr","AU-Cum","AU-DaS","AU-Dry","AU-Gin","AU-GWW",
                  "AU-How","AU-Stp","AU-TTE","AU-Tum","AU-Whr","AU-Wom")

library(lubridate)
library(magrittr)
library(dplyr)
library(coda)

# Collect the analysis outputs and name them with each site
for (Site in Sites){
  File = list.files("analysis/RTPV/",pattern = paste0("NEE_analysis_RTPV_",Site))
  message("Loading analysis output for ",Site," where file is ",File)
  load(paste0("analysis/RTPV/",File))
  assign(Site,output)
  rm(output)
}

message("Plotting weights for sites...")
# Do some shit with it
CumWeightParams = c(sprintf("cum_weightA[%d,%d]",rep(1:4,14),rep(1:14,each=4)),
                    sprintf("cum_weightAP[%d]",seq(1:8)))

# Extract the cumulative weights
CumWeights = data.frame("Sites"=rep(Sites,each = 64),
                        "Variable" = rep(rownames(eval(as.name(Sites[1]))$CumWeights),length(Sites)),
                        "Lag" = rep(0,each = 64*length(Sites)),
                        "Low" = unlist(lapply(Sites,function(x) eval(as.name(x))$CumWeights$WeightsLow)),
                        "Med" = unlist(lapply(Sites,function(x) eval(as.name(x))$CumWeights$WeightsMedian)),
                        "High" = unlist(lapply(Sites,function(x) eval(as.name(x))$CumWeights$WeightsHigh)))


# Assign lags
CumWeights$Lag[substr(CumWeights$Variable,15,18)=="1]"] = 0
CumWeights$Lag[substr(CumWeights$Variable,15,18)=="2]"] = 1
CumWeights$Lag[substr(CumWeights$Variable,15,18)=="3]"] = 2
CumWeights$Lag[substr(CumWeights$Variable,15,18)=="4]"] = 3
CumWeights$Lag[substr(CumWeights$Variable,15,18)=="5]"] = 4
CumWeights$Lag[substr(CumWeights$Variable,15,18)=="6]"] = 5
CumWeights$Lag[substr(CumWeights$Variable,15,18)=="7]"] = 6
CumWeights$Lag[substr(CumWeights$Variable,15,18)=="8]"] = 7
CumWeights$Lag[substr(CumWeights$Variable,15,18)=="9]"] = 8
CumWeights$Lag[substr(CumWeights$Variable,15,18)=="10]"] = 9
CumWeights$Lag[substr(CumWeights$Variable,15,18)=="11]"] = 10
CumWeights$Lag[substr(CumWeights$Variable,15,18)=="12]"] = 11
CumWeights$Lag[substr(CumWeights$Variable,15,18)=="13]"] = 12
CumWeights$Lag[substr(CumWeights$Variable,15,18)=="14]"] = 13
CumWeights$Lag[substr(CumWeights$Variable,14,15)=="1]"] = 0
CumWeights$Lag[substr(CumWeights$Variable,14,15)=="2]"] = 20
CumWeights$Lag[substr(CumWeights$Variable,14,15)=="3]"] = 29
CumWeights$Lag[substr(CumWeights$Variable,14,15)=="4]"] = 59
CumWeights$Lag[substr(CumWeights$Variable,14,15)=="5]"] = 119
CumWeights$Lag[substr(CumWeights$Variable,14,15)=="6]"] = 179
CumWeights$Lag[substr(CumWeights$Variable,14,15)=="7]"] = 269
CumWeights$Lag[substr(CumWeights$Variable,14,15)=="8]"] = 365

# Rename variables
CumWeights$Variable[substr(CumWeights$Variable,13,13)==1] = "Tair"
CumWeights$Variable[substr(CumWeights$Variable,13,13)==2] = "Fsd"
CumWeights$Variable[substr(CumWeights$Variable,13,13)==3] = "VPD"
CumWeights$Variable[substr(CumWeights$Variable,13,13)==4] = "PPTshort"
CumWeights$Variable[substr(CumWeights$Variable,11,12)=="AP"] = "PPTlong"

# Limit the dataframe to the variables requested
CumWeights = CumWeights[CumWeights$Variable == Var,]

# Rename variables to nice names
CumWeights$Variable[CumWeights$Variable == "Fsd"] = "Shortwave Radiation"
CumWeights$Variable[CumWeights$Variable == "Tair"] = "Air Temperature"
CumWeights$Variable[CumWeights$Variable == "PPTshort"] = "Short-term Precipitation"
CumWeights$Variable[CumWeights$Variable == "PPTlong"] = "Long-term Precipitation"

# Assign levels to Variable
CumWeights$Variable = factor(CumWeights$Variable,levels = sort(unique(CumWeights$Variable)))

# Assign vector length
if (Var == "PPTlong"){
  veclen = 8
} else {
  veclen = 14
}

# Summarise the weights
CritTime <- CumWeights %>%
                  group_by(Sites,Variable) %>%
                  mutate(PostLag = max(min(Lag[Med>0.5]),0),
                         PreLag = max(max(Lag[Med<0.5]),0),
                         Grad = max((Med[Lag==PostLag]-Med[Lag==PreLag])/(PostLag-PreLag),0),
                         Intercept = min(Med),
                         Lag50 = max(c((0.5-Med[Lag==PreLag])/Grad+PreLag,1-Intercept),na.rm=TRUE)) %>%
                select(Sites,Lag50) %>%
                unique() %>%
                as.data.frame()


load("site_data/SiteMetrics_worldclim_0.5res.Rdata")


DF = merge(CritTime,WorldClimMetrics,by="Sites")


DF = DF[DF$Transect %in% Transect,]

df = data.frame("Metric" = colnames(DF)[7:ncol(DF)],
                "rho" = NA,
                "p" = NA)

for (i in 7:ncol(DF)){
  test = cor.test(DF$Lag50,DF[,i],method="spearman")
  df$rho[i-6] = test$estimate
  df$p[i-6] = test$p.value
  if (df$p[i-6] < 0.05){
  print(colnames(DF)[i])
  print(df$rho[i-6])
  print(df$p[i-6])
  }
}

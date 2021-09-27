
rm(list=ls())

library(cowplot)

Sites = c("AU-ASM","AU-Cpr","AU-Cum","AU-DaS","AU-Dry","AU-Gin","AU-GWW",
          "AU-How","AU-Stp","AU-TTE","AU-Tum","AU-Whr","AU-Wom")

Transects = c("NATT","SAWS","SAWS","NATT","NATT","SAWS","SAWS",
              "NATT","NATT","NATT","SAWS","SAWS","SAWS")

#*******************************************************************************
# NEE Model Performance
#*******************************************************************************

source("functions/NEE_R2BarPlot_function_RTPV.R")
NEE_R2 = NEE_R2BarPlot_RTPV(Sites,Transects,"AnnualPPT", Clusters = 0)
#NEE_R2

#*******************************************************************************
# LE Model Performance
#*******************************************************************************

source("functions/LE_R2BarPlot_function_RTPV.R")
LE_R2 = LE_R2BarPlot_RTPV(Sites,Transects,"AnnualPPT", Clusters = 0)
#LE_R2
#
load("R2BarPlotLegend.Rdata")
plot = plot_grid(NEE_R2,LE_R2,labels = c("(a)","(b)"),label_size = 20,label_x = -0.01)
png("All_R2_legend.png",width = 1200, height = 400)
plot_grid(plot,legend,nrow=2,rel_heights = c(8,1))
dev.off()
#*******************************************************************************
# Flux Model Performance
#*******************************************************************************

# source("functions/Flux_R2BarPlot_function_RTPV.R")
# Plot = Flux_R2BarPlot_RTPV(Sites,Transects,"AnnualPPT", Clusters = 0)
# Plot

#*******************************************************************************
# Grouping Model Performance
#*******************************************************************************
NATT = c("AU-ASM","AU-DaS","AU-Dry","AU-How","AU-Stp","AU-TTE")
SAWS = c("AU-Cpr","AU-Cum","AU-Gin","AU-GWW","AU-Tum","AU-Whr","AU-Wom")

# Plot = NEE_R2BarPlot_RTPV(NATT,"NATT","PPTSeasonality", Clusters = 0)
# Plot
# 
# Plot = NEE_R2BarPlot_RTPV(SAWS,"SAWS","AnnualMeanTemp", Clusters = 0)
# Plot

# Plot = LE_R2BarPlot_RTPV(NATT,"NATT","AnnualPPT", Clusters = 0)
# Plot
# 
# Plot = LE_R2BarPlot_RTPV(SAWS,"SAWS","PPTDryQtr", Clusters = 0)
# Plot

#*******************************************************************************
# Climate Sensitivity
#*******************************************************************************

source("functions/NEE_SDScaledSensitivity_RTPV.R")
NEE_Sen = NEE_SDScaledSensitivity_RTPV(Sites,Vars=c("Tair","Fsd","VPD","PPTshort","PPTlong"))
#NEE_Sen

source("functions/LE_SDScaledSensitivity_RTPV.R")
LE_Sen = LE_SDScaledSensitivity_RTPV(Sites,Vars=c("Tair","Fsd","VPD","PPTshort","PPTlong"))
#LE_Sen

png("All_Sen.png",width = 1200, height = 800)
plot_grid(NEE_Sen,LE_Sen,labels = c("(a)","(b)"),label_size = 20, label_x = -0.01)
dev.off()

#*******************************************************************************
# Climate Sensitivity
#*******************************************************************************

source("functions/SensitivityBoxPlot_RTPV.R")
out = SensitivityBoxPlot_RTPV (Sites,Transects,Vars=c("Tair","Fsd","VPD","PPTshort","PPTlong"))

png("Sen_Point.png",width = 800, height = 800)
plot(out$PointPlot)
dev.off()

png("Sen_Box.png",width = 1200, height = 800)
plot(out$BoxPlot)
dev.off()

#*******************************************************************************
# NEE stacked weights with different colours per transect
#*******************************************************************************

source("functions/NEE_StackedWeightPlot_Transects_RTPV.R")
NEETairSites = c("AU-Cpr","AU-Cum","AU-DaS","AU-Gin","AU-GWW",
          "AU-How","AU-Stp","AU-Tum","AU-Whr","AU-Wom")

NEETairTransects = c("SAWS","SAWS","NATT","SAWS","SAWS",
              "NATT","NATT","SAWS","SAWS","SAWS")

NEE_Tair = NEE_StackedWeightPlot_RTPV(NEETairSites,NEETairTransects,"Tair","AnnualPPT")
#NEE_Tair

NEE_PPT = NEE_StackedWeightPlot_RTPV(Sites,Transects,"PPTlong","AnnualPPT")
#NEE_PPT
png("NEE_lag_TairPPTlong.png",width = 1200, height = 600)
plot_grid(NEE_Tair,NEE_PPT,labels = c("(a)","(b)"),label_size = 20, hjust = 0)
dev.off()
#*******************************************************************************
# LE stacked weights with different colours per transect
#*******************************************************************************

source("functions/LE_StackedWeightPlot_Transects_RTPV.R")
LETairSites = c("AU-ASM","AU-Cpr","AU-Cum","AU-DaS","AU-Gin","AU-GWW",
              "AU-How","AU-Stp","AU-TTE","AU-Tum","AU-Whr","AU-Wom")

LETairTransects = c("NATT","SAWS","SAWS","NATT","SAWS","SAWS",
                  "NATT","NATT","NATT","SAWS","SAWS","SAWS")

LE_Tair = LE_StackedWeightPlot_RTPV(LETairSites,LETairTransects,"Tair","AnnualPPT")
#LE_Tair

LEPPTSites = c("AU-ASM","AU-Cpr","AU-Cum","AU-DaS","AU-Dry","AU-Gin","AU-GWW",
            "AU-How","AU-Stp","AU-TTE","AU-Whr","AU-Wom")

LEPPTTransects = c("NATT","SAWS","SAWS","NATT","NATT","SAWS","SAWS",
                "NATT","NATT","NATT","SAWS","SAWS")

LE_PPT = LE_StackedWeightPlot_RTPV(LEPPTSites,LEPPTTransects,"PPTlong","AnnualPPT")
#LE_PPT
png("LE_lag_TairPPTlong.png",width = 1200, height = 600)
plot_grid(LE_Tair,LE_PPT,labels = c("(a)","(b)"),label_size = 20, label_x = -0.02)
dev.off()
#*******************************************************************************
# Metric Improvements - SUPPLEMENTARY
#*******************************************************************************

source("functions/NEE_MetricsPlot_function_RTPV.R")
png("NEE_Metrics.png",width = 1200, height = 1200)
NEE_Metrics = NEE_MetricsPlot_function_RTPV(Sites)
NEE_Metrics
dev.off()

source("functions/LE_MetricsPlot_function_RTPV.R")
png("LE_Metrics.png",width = 1200, height = 1200)
LE_Metrics = LE_MetricsPlot_function_RTPV(Sites)
LE_Metrics
dev.off()
#*******************************************************************************
# Model Performance Time Series - SUPPLEMENTARY
#*******************************************************************************

source("functions/NEE_DailyObsVsPred_MA_RTPV.R")
for (Site in Sites){
  Plots = NEE_DailyObsVsPred_MA_RTPV(Site,k=30)
  plot(Plots$StackedPlot)
}

#*******************************************************************************
# Individual Weight Plots - SUPPLEMENTARY
#*******************************************************************************

source("functions/NEE_WeightPlot_RTPV.R")
png("NEE_NATT_Weights.png",width = 1200, height = 700)
NEE_NATT_Weights = NEE_WeightPlot_RTPV(NATT)
NEE_NATT_Weights
dev.off()

png("NEE_SAWS_Weights.png",width = 1200, height = 800)
NEE_SAWS_Weights = NEE_WeightPlot_RTPV(SAWS)
NEE_SAWS_Weights
dev.off()

source("functions/LE_WeightPlot_RTPV.R")
png("LE_NATT_Weights.png",width = 1200, height = 700)
LE_NATT_Weights = LE_WeightPlot_RTPV(NATT)
LE_NATT_Weights
dev.off()

png("LE_SAWS_Weights.png",width = 1200, height = 800)
LE_SAWS_Weights = LE_WeightPlot_RTPV(SAWS)
LE_SAWS_Weights
dev.off()


#*******************************************************************************
# TEST - Stacked Weights
#*******************************************************************************

source("functions/NEE_StackedWeightPlot_TransectsSplitPlot_RTPV.R")
NEETairSites = c("AU-Cpr","AU-Cum","AU-DaS","AU-Gin","AU-GWW",
                 "AU-How","AU-Stp","AU-Tum","AU-Whr","AU-Wom")

NEETairTransects = c("SAWS","SAWS","NATT","SAWS","SAWS",
                     "NATT","NATT","SAWS","SAWS","SAWS")

NEE_Tair = NEE_StackedWeightPlot_split_RTPV(NEETairSites,NEETairTransects,"Tair","AnnualPPT")
#NEE_Tair

NEE_PPT = NEE_StackedWeightPlot_split_RTPV(Sites,Transects,"PPTlong","AnnualPPT")
#NEE_PPT
png("NEE_lag_TairPPTlong_split.png",width = 1200, height = 800)
plot_grid(NEE_Tair,NEE_PPT,
          nrow = 2,
          labels = c("(a) Air Temperature",
                     "(b) Long-term Precipitation"),
          label_size = 20,
          hjust = 0)
dev.off()


source("functions/LE_StackedWeightPlot_TransectsSplitPlot_RTPV.R")
LETairSites = c("AU-ASM","AU-Cpr","AU-Cum","AU-DaS","AU-Gin","AU-GWW",
                "AU-How","AU-Stp","AU-TTE","AU-Tum","AU-Whr","AU-Wom")

LETairTransects = c("NATT","SAWS","SAWS","NATT","SAWS","SAWS",
                    "NATT","NATT","NATT","SAWS","SAWS","SAWS")

LE_Tair = LE_StackedWeightPlot_split_RTPV(LETairSites,LETairTransects,"Tair","AnnualPPT")
#LE_Tair

LEPPTSites = c("AU-ASM","AU-Cpr","AU-Cum","AU-DaS","AU-Dry","AU-Gin","AU-GWW",
               "AU-How","AU-Stp","AU-TTE","AU-Whr","AU-Wom")

LEPPTTransects = c("NATT","SAWS","SAWS","NATT","NATT","SAWS","SAWS",
                   "NATT","NATT","NATT","SAWS","SAWS")

LE_PPT = LE_StackedWeightPlot_split_RTPV(LEPPTSites,LEPPTTransects,"PPTlong","AnnualPPT")
#LE_PPT
png("LE_lag_TairPPTlong_split.png",width = 1200, height = 800)
plot_grid(LE_Tair,LE_PPT,
          nrow = 2,
          labels = c("(a) Air Temperature",
                     "(b) Long-term Precipitation"),
          label_size = 20,
          hjust = 0)
dev.off()



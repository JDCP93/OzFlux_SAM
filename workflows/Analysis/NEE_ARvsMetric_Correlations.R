
# Make sure everything is clean
rm(list=ls())

# Load the required packages
library(tidyverse)
library(lubridate)
library(viridis)
library(ggpubr)
library(gridExtra)

# List all sites
Sites = c("AU-ASM"
          ,"AU-Cpr"
          ,"AU-Cum"
          ,"AU-DaS"
          ,"AU-Dry"
          ,"AU-Gin"
          ,"AU-GWW"
          ,"AU-How"
          ,"AU-Stp"
          ,"AU-TTE"
          ,"AU-Tum"
          ,"AU-Whr"
          ,"AU-Wom"
)

# Assign transects
Transects = c("NATT",
              "SAWS",
              "SAWS",
              "NATT",
              "NATT",
              "SAWS",
              "SAWS",
              "NATT",
              "NATT",
              "NATT",
              "SAWS",
              "SAWS",
              "SAWS"
)

# Source worldclim correlations and climate metrics
load("site_data/SiteMetrics_worldclim_0.5res_AI.Rdata")

# Initialise the dataframe of model R2 values
R2 = data.frame("Site" = Sites,
                "Transect" = Transects,
                "R2.CUR" = 0,
                "R2.SAM" = 0,
                "R2.AR1" = 0)

# For each site
for (Site in Sites){
  
  # Collect the R2 values from the analysis scripts
  message("Collating R2 values for ",Site)
  
  # Load the analysis results
  File = list.files("analysis/RTPV/",pattern = paste0("NEE_analysis_RTPV_",Site))
  load(paste0("analysis/RTPV/",File))
  R2$R2.SAM[R2$Site==Site] = output$SAM.R2
  
  File = list.files("analysis/RTPV/",pattern = paste0("NEE_current_analysis_RTPV_",Site))
  load(paste0("analysis/RTPV/",File))
  R2$R2.CUR[R2$Site==Site] = output$CUR.R2
  
  File = list.files("analysis/RTPV/",pattern = paste0("NEE_AR1_analysis_RTPV_",Site))
  load(paste0("analysis/RTPV/",File))
  R2$R2.AR1[R2$Site==Site] = output$AR1.R2
}

# Initialise the dataframe for the correlation values and p values
Correlations = data.frame("Metric" = colnames(WorldClimMetrics[-(1:4)]),
                          "AbsVal" = 0,
                          "AbsP" = 0,
                          "AbsImpVal" = 0,
                          "AbsImpP" = 0,
                          "RelImpVal" = 0,
                          "RelImpP" = 0,
                          "NATTAbsVal" = 0,
                          "NATTAbsP" = 0,
                          "NATTAbsImpVal" = 0,
                          "NATTAbsImpP" = 0,
                          "NATTRelImpVal" = 0,
                          "NATTRelImpP" = 0,
                          "SAWSAbsVal" = 0,
                          "SAWSAbsP" = 0,
                          "SAWSAbsImpVal" = 0,
                          "SAWSAbsImpP" = 0,
                          "SAWSRelImpVal" = 0,
                          "SAWSRelImpP" = 0)

# For each metric
for (i in 5:ncol(WorldClimMetrics)){
  
  # Calculate the correlations across every site
  # For absolute AR1 R2 value
  metricR2 = cor.test(x = WorldClimMetrics[,i], 
                      y = R2$R2.AR1, 
                      method = "spearman")$estimate
  metricPvalue = cor.test(x = WorldClimMetrics[,i], 
                          y = R2$R2.AR1, 
                          method = "spearman")$p.value
  message("Absolute memory strength and ",
          colnames(WorldClimMetrics)[i],
          " are correlated with R2 value ",
          round(metricR2,3),
          " and p value ",
          round(metricPvalue,3))
  Correlations$AbsVal[i-4] = metricR2
  Correlations$AbsP[i-4] = metricPvalue
  # For absolute R2 improvement between AR1 and EM
  metricR2 = cor.test(x = WorldClimMetrics[,i], 
                      y = (R2$R2.AR1-R2$R2.SAM), 
                      method = "spearman")$estimate
  metricPvalue = cor.test(x = WorldClimMetrics[,i], 
                          y = (R2$R2.AR1-R2$R2.SAM), 
                          method = "spearman")$p.value
  message("Absolute memory improvement and ",
          colnames(WorldClimMetrics)[i],
          " are correlated with R2 value ",
          round(metricR2,3),
          " and p value ",
          round(metricPvalue,3))
  Correlations$AbsImpVal[i-4] = metricR2
  Correlations$AbsImpP[i-4] = metricPvalue
  # For relative R2 improvement between AR1 and EM
  metricR2 = cor.test(x = WorldClimMetrics[,i],
                      y = (R2$R2.AR1-R2$R2.SAM)/R2$R2.SAM,
                      method = "spearman")$estimate
  metricPvalue = cor.test(x = WorldClimMetrics[,i],
                          y = (R2$R2.AR1-R2$R2.SAM)/R2$R2.SAM,
                          method = "spearman")$p.value
  message("Relative memory improvement and ",
          colnames(WorldClimMetrics)[i],
          " are correlated with R2 value ",
          round(metricR2,3),
          " and p value ",
          round(metricPvalue,3),"\n")
  Correlations$RelImpVal[i-4] = metricR2
  Correlations$RelImpP[i-4] = metricPvalue
  
  # Calculate correlations for just NATT sites
  # For absolute AR1 R2 value
  metricR2 = cor.test(x = WorldClimMetrics[WorldClimMetrics$Transect=="NATT",i],
                      y = (R2$R2.AR1[R2$Transect=="NATT"]),
                      method = "spearman")$estimate
  metricPvalue = cor.test(x = WorldClimMetrics[WorldClimMetrics$Transect=="NATT",i],
                          y = (R2$R2.AR1[R2$Transect=="NATT"]),
                          method = "spearman")$p.value
  Correlations$NATTAbsVal[i-4] = metricR2
  Correlations$NATTAbsP[i-4] = metricPvalue
  # For absolute R2 improvement between AR1 and EM  
  metricR2 = cor.test(x = WorldClimMetrics[WorldClimMetrics$Transect=="NATT",i],
                      y = (R2$R2.AR1[R2$Transect=="NATT"]-R2$R2.SAM[R2$Transect=="NATT"]),
                      method = "spearman")$estimate
  metricPvalue = cor.test(x = WorldClimMetrics[WorldClimMetrics$Transect=="NATT",i],
                          y = (R2$R2.AR1[R2$Transect=="NATT"]-R2$R2.SAM[R2$Transect=="NATT"]),
                          method = "spearman")$p.value
  Correlations$NATTAbsImpVal[i-4] = metricR2
  Correlations$NATTAbsImpP[i-4] = metricPvalue
  # For relative R2 improvement between AR1 and EM
  metricR2 = cor.test(x = WorldClimMetrics[WorldClimMetrics$Transect=="NATT",i],
                      y = (R2$R2.AR1[R2$Transect=="NATT"]-R2$R2.SAM[R2$Transect=="NATT"])/R2$R2.SAM[R2$Transect=="NATT"],
                      method = "spearman")$estimate
  metricPvalue = cor.test(x = WorldClimMetrics[WorldClimMetrics$Transect=="NATT",i],
                          y = (R2$R2.AR1[R2$Transect=="NATT"]-R2$R2.SAM[R2$Transect=="NATT"])/R2$R2.SAM[R2$Transect=="NATT"],
                          method = "spearman")$p.value
  Correlations$NATTRelImpVal[i-4] = metricR2
  Correlations$NATTRelImpP[i-4] = metricPvalue
  
  
  # Calculate correlations for just SAWS sites
  # For absolute AR1 R2 value
  metricR2 = cor.test(x = WorldClimMetrics[WorldClimMetrics$Transect=="SAWS",i],
                      y = (R2$R2.AR1[R2$Transect=="SAWS"]),
                      method = "spearman")$estimate
  metricPvalue = cor.test(x = WorldClimMetrics[WorldClimMetrics$Transect=="SAWS",i],
                          y = (R2$R2.AR1[R2$Transect=="SAWS"]),
                          method = "spearman")$p.value
  Correlations$SAWSAbsVal[i-4] = metricR2
  Correlations$SAWSAbsP[i-4] = metricPvalue
  # For absolute R2 improvement between AR1 and EM  
  metricR2 = cor.test(x = WorldClimMetrics[WorldClimMetrics$Transect=="SAWS",i],
                      y = (R2$R2.AR1[R2$Transect=="SAWS"]-R2$R2.SAM[R2$Transect=="SAWS"]),
                      method = "spearman")$estimate
  metricPvalue = cor.test(x = WorldClimMetrics[WorldClimMetrics$Transect=="SAWS",i],
                          y = (R2$R2.AR1[R2$Transect=="SAWS"]-R2$R2.SAM[R2$Transect=="SAWS"]),
                          method = "spearman")$p.value
  Correlations$SAWSAbsImpVal[i-4] = metricR2
  Correlations$SAWSAbsImpP[i-4] = metricPvalue
  # For relative R2 improvement between AR1 and EM
  metricR2 = cor.test(x = WorldClimMetrics[WorldClimMetrics$Transect=="SAWS",i],
                      y = (R2$R2.AR1[R2$Transect=="SAWS"]-R2$R2.SAM[R2$Transect=="SAWS"])/R2$R2.SAM[R2$Transect=="SAWS"],
                      method = "spearman")$estimate
  metricPvalue = cor.test(x = WorldClimMetrics[WorldClimMetrics$Transect=="SAWS",i],
                          y = (R2$R2.AR1[R2$Transect=="SAWS"]-R2$R2.SAM[R2$Transect=="SAWS"])/R2$R2.SAM[R2$Transect=="SAWS"],
                          method = "spearman")$p.value
  Correlations$SAWSRelImpVal[i-4] = metricR2
  Correlations$SAWSRelImpP[i-4] = metricPvalue
}

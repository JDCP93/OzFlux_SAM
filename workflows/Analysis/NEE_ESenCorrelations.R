
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
Vars = c("Tair","Fsd","VPD","PPTshort","PPTlong")
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
WorldClimMetrics = WorldClimMetrics[WorldClimMetrics$Sites %in% Sites,]

# Initialise dataframe for weights to apply to sensitivities
Weights = data.frame("Site" = Sites,
                     "Tair" = NA,
                     "Fsd" = NA,
                     "VPD" = NA,
                     "PPTshort" = NA,
                     "PPTlong" = NA,
                     "PPT" = NA)

# Collect the analysis outputs and name them with each site
for (Site in Sites){
  # Load analysis for ESen 
  File = list.files("analysis/RTPV/",pattern = paste0("NEE_analysis_RTPV_",Site))
  load(paste0("analysis/RTPV/",File))
  assign(Site,output)
  rm(output)
  
  # Load the coefficients
  File = list.files("analysis/RTPV/",pattern = paste0("NEE_summary_RTPV_",Site))
  load(paste0("analysis/RTPV/",File))
  
  # Load the input file and extract required data
  load(paste0("inputs/RTPV/",Site,"_Input_RTPV.Rdata"))
  input = eval(as.name(paste0(Site,"_Input")))
  
  # We now calculate the weighted climate timesteps on a site-by-site basis
  # Combine climate data and precip data
  climend = nrow(input$clim)
  climate = cbind(input$clim,
                  rbind(matrix(NA,1,4),input$clim[-climend,]),
                  rbind(matrix(NA,2,4),input$clim[-((climend-1):climend),]),
                  rbind(matrix(NA,3,4),input$clim[-((climend-2):climend),]),
                  rbind(matrix(NA,4,4),input$clim[-((climend-3):climend),]),
                  rbind(matrix(NA,5,4),input$clim[-((climend-4):climend),]),
                  rbind(matrix(NA,6,4),input$clim[-((climend-5):climend),]),
                  rbind(matrix(NA,7,4),input$clim[-((climend-6):climend),]),
                  rbind(matrix(NA,8,4),input$clim[-((climend-7):climend),]),
                  rbind(matrix(NA,9,4),input$clim[-((climend-8):climend),]),
                  rbind(matrix(NA,10,4),input$clim[-((climend-9):climend),]),
                  rbind(matrix(NA,11,4),input$clim[-((climend-10):climend),]),
                  rbind(matrix(NA,12,4),input$clim[-((climend-11):climend),]),
                  rbind(matrix(NA,13,4),input$clim[-((climend-12):climend),]),
                  input$ppt_multiscale)
  # Rename columns for ease
  colnames(climate) = c("Tair_1",
                        "Fsd_1",
                        "VPD_1",
                        "PPTshort_1",
                        "Tair_2",
                        "Fsd_2",
                        "VPD_2",
                        "PPTshort_2",
                        "Tair_3",
                        "Fsd_3",
                        "VPD_3",
                        "PPTshort_3",
                        "Tair_4",
                        "Fsd_4",
                        "VPD_4",
                        "PPTshort_4",
                        "Tair_5",
                        "Fsd_5",
                        "VPD_5",
                        "PPTshort_5",
                        "Tair_6",
                        "Fsd_6",
                        "VPD_6",
                        "PPTshort_6",
                        "Tair_7",
                        "Fsd_7",
                        "VPD_7",
                        "PPTshort_7",
                        "Tair_8",
                        "Fsd_8",
                        "VPD_8",
                        "PPTshort_8",
                        "Tair_9",
                        "Fsd_9",
                        "VPD_9",
                        "PPTshort_9",
                        "Tair_10",
                        "Fsd_10",
                        "VPD_10",
                        "PPTshort_10",
                        "Tair_11",
                        "Fsd_11",
                        "VPD_11",
                        "PPTshort_11",
                        "Tair_12",
                        "Fsd_12",
                        "VPD_12",
                        "PPTshort_12",
                        "Tair_13",
                        "Fsd_13",
                        "VPD_13",
                        "PPTshort_13",
                        "Tair_14",
                        "Fsd_14",
                        "VPD_14",
                        "PPTshort_14",
                        "PPTlong_1",
                        "PPTlong_2",
                        "PPTlong_3",
                        "PPTlong_4",
                        "PPTlong_5",
                        "PPTlong_6",
                        "PPTlong_7",
                        "PPTlong_8")
  
  # Remove first year, which has no PPT data
  climate = climate[-(1:365),]
  NEE = input$NEE[-(1:365)]
  
  # Extract climate weights
  weightnames = c(sprintf("weightA[%d,%d]",rep(1:4,14),rep(1:14,each=4)),
                  sprintf("weightAP[%d]",seq(1:8)))
  weights = output.summary$statistics[rownames(output.summary$statistics) %in% weightnames,1]
  # Rename weights
  names(weights)=str_replace_all(names(weights),pattern="weightA\\[1,","Tair_")
  names(weights)=str_replace_all(names(weights),pattern="weightA\\[2,","Fsd_")
  names(weights)=str_replace_all(names(weights),pattern="weightA\\[3,","VPD_")
  names(weights)=str_replace_all(names(weights),pattern="weightA\\[4,","PPTshort_")
  names(weights)=str_replace_all(names(weights),pattern="weightAP\\[","PPTlong_")
  names(weights)=str_replace_all(names(weights),pattern="\\]","")
  
  # Initialise matrix
  wts = data.frame("Timestep" = 1:nrow(climate),
                   "Tair" = NA,
                   "Fsd" = NA,
                   "VPD" = NA,
                   "PPTshort" = NA,
                   "PPTlong" = NA,
                   "PPT" = NA)
  
  # Calculated the weighted climate term at each timestep
  for (t in wts$Timestep){
    wts$Tair[t] = sum(weights[substr(names(weights),1,4)=="Tair"]*climate[t,substr(colnames(climate),1,4)=="Tair"])
    wts$Fsd[t] = sum(weights[substr(names(weights),1,3)=="Fsd"]*climate[t,substr(colnames(climate),1,3)=="Fsd"])
    wts$VPD[t] = sum(weights[substr(names(weights),1,3)=="VPD"]*climate[t,substr(colnames(climate),1,3)=="VPD"])
    wts$PPTshort[t] = sum(weights[substr(names(weights),1,4)=="PPTs"]*climate[t,substr(colnames(climate),1,4)=="PPTs"])
    wts$PPTlong[t] = sum(weights[substr(names(weights),1,4)=="PPTl"]*climate[t,substr(colnames(climate),1,4)=="PPTl"])
    wts$PPT = wts$PPTshort+wts$PPTlong
  }
  
  # Calculate the standard deviation of the weighted climate
  Weights[Weights$Site==Site,2:7] = sqrt(diag(var(wts)))[2:7]
  rm(output.summary)
}

# Extract the sensitivity covariates
NEE_ESen = data.frame("Site"=rep(Sites,each = 6),
                      "Transect" = rep(Transects,each=6),
                      "Variable" = rep(rownames(eval(as.name(Sites[1]))$ESen),length(Sites)),
                      "Low" = unlist(lapply(Sites,function(x) eval(as.name(x))$ESen$ESenLow)),
                      "Med" = unlist(lapply(Sites,function(x) eval(as.name(x))$ESen$ESenMedian)),
                      "High" = unlist(lapply(Sites,function(x) eval(as.name(x))$ESen$ESenHigh)))


# Weighted ESen
NEE_ESen[NEE_ESen$Variable=="Tair",4:6] = NEE_ESen[NEE_ESen$Variable=="Tair",4:6]*Weights$Tair
NEE_ESen[NEE_ESen$Variable=="Fsd",4:6] = NEE_ESen[NEE_ESen$Variable=="Fsd",4:6]*Weights$Fsd
NEE_ESen[NEE_ESen$Variable=="VPD",4:6] = NEE_ESen[NEE_ESen$Variable=="VPD",4:6]*Weights$VPD
NEE_ESen[NEE_ESen$Variable=="PPTshort",4:6] = NEE_ESen[NEE_ESen$Variable=="PPTshort",4:6]*Weights$PPTshort
NEE_ESen[NEE_ESen$Variable=="PPTlong",4:6] = NEE_ESen[NEE_ESen$Variable=="PPTlong",4:6]*Weights$PPTlong
NEE_ESen[NEE_ESen$Variable=="PPT",4:6] = NEE_ESen[NEE_ESen$Variable=="PPT",4:6]*Weights$PPT

# Limit the dataframe to the variables requested
NEE_ESen = NEE_ESen[NEE_ESen$Variable %in% Vars,]

#NEE_ESen$Med = abs(NEE_ESen$Med)

ESen = NEE_ESen %>% pivot_wider(id_cols = c("Site","Transect"),names_from = "Variable",values_from = "Med")

# Initialise the dataframe for the correlation values and p values
Correlations = data.frame(matrix(nrow=5*3*ncol(WorldClimMetrics[-(1:4)]),ncol=6))
colnames(Correlations) = c("Flux","Metric","Variable","Transect","CCO","Pval")
Correlations["Metric"] = rep(colnames(WorldClimMetrics[-(1:4)]),each = 15)
Correlations["Variable"] = rep(Vars,each=3)
Correlations["Transect"] = c("All","NATT","SAWS")
Correlations["Flux"] = "NEE"

# For each metric
for (i in 5:ncol(WorldClimMetrics)){
  for (Var in Vars){
    Loopmetric = colnames(WorldClimMetrics)[i]
    # Calculate the correlations across every site
    test = cor.test(x = WorldClimMetrics[,i], 
                    y = pull(ESen[,Var]), 
                    method = "spearman")
    R2 = round(test$estimate,2)
    P = round(test$p.value,3)
    message(Var,
            " and ",
            Loopmetric,
            " are correlated with R2 value ",
            round(R2,3),
            " and p value ",
            round(P,3))
    Correlations[Correlations$Metric==Loopmetric & Correlations$Variable==Var & Correlations$Transect == "All","CCO"] = R2
    Correlations[Correlations$Metric==Loopmetric & Correlations$Variable==Var & Correlations$Transect == "All","Pval"] = P
    
    # Calculate correlations for just NATT sites
    test = cor.test(x = WorldClimMetrics[WorldClimMetrics$Transect=="NATT",i],
                    y = pull(ESen[ESen$Transect=="NATT",Var]),
                    method = "spearman")
    R2 = round(test$estimate,2)
    P = round(test$p.value,3)
    message(Var,
            " and ",
            Loopmetric,
            " are correlated in the NATT with R2 value ",
            round(R2,3),
            " and p value ",
            round(P,3))
    Correlations[Correlations$Metric==Loopmetric & Correlations$Variable==Var & Correlations$Transect == "NATT","CCO"] = R2
    Correlations[Correlations$Metric==Loopmetric & Correlations$Variable==Var & Correlations$Transect == "NATT","Pval"] = P
    
    # Calculate correlations for just SAWS sites
    test = cor.test(x = WorldClimMetrics[WorldClimMetrics$Transect=="SAWS",i],
                    y = pull(ESen[ESen$Transect=="SAWS",Var]),
                    method = "spearman")
    R2 = round(test$estimate,2)
    P = round(test$p.value,3)
    message(Var,
            " and ",
            Loopmetric,
            " are correlated in the SAWS with R2 value ",
            round(R2,3),
            " and p value ",
            round(P,3))
    Correlations[Correlations$Metric==Loopmetric & Correlations$Variable==Var & Correlations$Transect == "SAWS","CCO"] = R2
    Correlations[Correlations$Metric==Loopmetric & Correlations$Variable==Var & Correlations$Transect == "SAWS","Pval"] = P
  }
}

NEECorrelations = Correlations

save(NEECorrelations, file = "NEE_ESenMetricCorrelations.Rdata")
write_excel_csv(NEECorrelations, file = "NEE_ESenMetricCorrelations.csv")

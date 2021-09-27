## Workflow to find the correlation between MAP and model performance for NEE and LE
rm(list=ls())

Sites = c("AU-ASM","AU-Cpr","AU-Cum","AU-DaS","AU-Dry","AU-Gin","AU-GWW",
          "AU-How","AU-Stp","AU-TTE","AU-Tum","AU-Whr","AU-Wom")

# Source worldclim correlations and climate metrics
load("site_data/SiteMetrics_worldclim_0.5res.Rdata")

# Initialise dataframe
NEEdf = data.frame("Site" = rep(Sites,each = 15),
                "Model" = rep(c("Current Environment (SAM)",
                                "Environmental Memory (SAM)",
                                "Biological Memory (AR1)"),
                              each = 5,
                              times = length(Sites)),
                "Metric" = rep(c("R^2",
                                 "Mean Bias Error",
                                 "Normalised Mean Error",
                                 "Std. Dev. Difference",
                                 "Correlation Coeff."),
                               times = 3*length(Sites)),
                "Value" = NA)

# For all sites
for (Site in Sites){
  # Load the metrics from each of the 3 model runs
  File = list.files("analysis/RTPV/metrics/",
                    pattern = paste0("NEE_current_metrics_RTPV_",Site))
  load(paste0("analysis/RTPV/metrics/",File))
  NEEdf$Value[NEEdf$Site == Site & NEEdf$Model == "Current Environment (SAM)"] = unlist(output[c(1:5)])
  
  File = list.files("analysis/RTPV/metrics/",
                    pattern = paste0("NEE_metrics_RTPV_",Site))
  load(paste0("analysis/RTPV/metrics/",File))
  NEEdf$Value[NEEdf$Site == Site & NEEdf$Model == "Environmental Memory (SAM)"] = unlist(output[c(1:5)])
  
  File = list.files("analysis/RTPV/metrics/",
                    pattern = paste0("NEE_AR1_metrics_RTPV_",Site))
  load(paste0("analysis/RTPV/metrics/",File))
  NEEdf$Value[NEEdf$Site == Site & NEEdf$Model == "Biological Memory (AR1)"] = unlist(output[c(1:5)])
}

# Turn Model into Factor
NEEdf$Model = factor(NEEdf$Model,
                  levels = c("Current Environment (SAM)",
                             "Environmental Memory (SAM)",
                             "Biological Memory (AR1)"))
# Turn Metric into Factor
NEEdf$Metric = factor(NEEdf$Metric, 
                   levels = c("Normalised Mean Error",
                              "Std. Dev. Difference",
                              "Mean Bias Error",
                              "R^2",
                              "Correlation Coeff."),
                   labels = c(expression("Normalised~Mean~Error"),
                              expression("Std.~Dev.~Difference"),
                              expression("Mean~Bias~Error"),
                              expression("R^2"),
                              expression("Correlation~Coeff.")))

# Correlation test
cor.test(NEEdf$Value[NEEdf$Model=="Current Environment (SAM)" & NEEdf$Metric=="R^2"],WorldClimMetrics$AnnualPPT, method = "spearman")
cor.test(NEEdf$Value[NEEdf$Model=="Environmental Memory (SAM)" & NEEdf$Metric=="R^2"],WorldClimMetrics$AnnualPPT, method = "spearman")
cor.test(NEEdf$Value[NEEdf$Model=="Biological Memory (AR1)" & NEEdf$Metric=="R^2"],WorldClimMetrics$AnnualPPT, method = "spearman")

# Test just for NATT
cor.test(NEEdf$Value[NEEdf$Model=="Current Environment (SAM)" & NEEdf$Metric=="R^2" & NEEdf$Site %in% c("AU-ASM","AU-DaS","AU-Dry","AU-How","AU-Stp","AU-TTE")],WorldClimMetrics$AnnualPPT[WorldClimMetrics$Transect=="NATT"], method = "spearman")
cor.test(NEEdf$Value[NEEdf$Model=="Environmental Memory (SAM)" & NEEdf$Metric=="R^2" & NEEdf$Site %in% c("AU-ASM","AU-DaS","AU-Dry","AU-How","AU-Stp","AU-TTE")],WorldClimMetrics$AnnualPPT[WorldClimMetrics$Transect=="NATT"], method = "spearman")
cor.test(NEEdf$Value[NEEdf$Model=="Biological Memory (AR1)" & NEEdf$Metric=="R^2" & NEEdf$Site %in% c("AU-ASM","AU-DaS","AU-Dry","AU-How","AU-Stp","AU-TTE")],WorldClimMetrics$AnnualPPT[WorldClimMetrics$Transect=="NATT"], method = "spearman")

# Test just for SAWS
cor.test(NEEdf$Value[NEEdf$Model=="Current Environment (SAM)" & NEEdf$Metric=="R^2" & NEEdf$Site %in% c("AU-Cpr","AU-Cum","AU-Gin","AU-GWW","AU-Tum","AU-Whr","AU-Wom")],WorldClimMetrics$AnnualPPT[WorldClimMetrics$Transect=="SAWS"], method = "spearman")
cor.test(NEEdf$Value[NEEdf$Model=="Environmental Memory (SAM)" & NEEdf$Metric=="R^2" & NEEdf$Site %in% c("AU-Cpr","AU-Cum","AU-Gin","AU-GWW","AU-Tum","AU-Whr","AU-Wom")],WorldClimMetrics$AnnualPPT[WorldClimMetrics$Transect=="SAWS"], method = "spearman")
cor.test(NEEdf$Value[NEEdf$Model=="Biological Memory (AR1)" & NEEdf$Metric=="R^2" & NEEdf$Site %in% c("AU-Cpr","AU-Cum","AU-Gin","AU-GWW","AU-Tum","AU-Whr","AU-Wom")],WorldClimMetrics$AnnualPPT[WorldClimMetrics$Transect=="SAWS"], method = "spearman")

##
## Repeat the above for latent heat
##
 
LEdf = data.frame("Site" = rep(Sites,each = 15),
                "Model" = rep(c("Current Environment (SAM)",
                                "Environmental Memory (SAM)",
                                "Biological Memory (AR1)"),
                              each = 5,
                              times = length(Sites)),
                "Metric" = rep(c("R^2",
                                 "Mean Bias Error",
                                 "Normalised Mean Error",
                                 "Std. Dev. Difference",
                                 "Correlation Coeff."),
                               times = 3*length(Sites)),
                "Value" = NA)


for (Site in Sites){
  
  File = list.files("analysis/RTPV/",
                    pattern = paste0("LE_current_analysis_RTPV_",Site))
  load(paste0("analysis/RTPV/",File))
  LEdf$Value[LEdf$Site == Site & LEdf$Model == "Current Environment (SAM)"] = unlist(output[c(4:8)])
  
  File = list.files("analysis/RTPV/",
                    pattern = paste0("LE_analysis_RTPV_",Site))
  load(paste0("analysis/RTPV/",File))
  LEdf$Value[LEdf$Site == Site & LEdf$Model == "Environmental Memory (SAM)"] = unlist(output[c(4:8)])
  
  File = list.files("analysis/RTPV/",
                    pattern = paste0("LE_AR1_analysis_RTPV_",Site))
  load(paste0("analysis/RTPV/",File))
  LEdf$Value[LEdf$Site == Site & LEdf$Model == "Biological Memory (AR1)"] = unlist(output[c(4:8)])
}

LEdf$Model = factor(LEdf$Model,
                  levels = c("Current Environment (SAM)",
                             "Environmental Memory (SAM)",
                             "Biological Memory (AR1)"))
LEdf$Metric = factor(LEdf$Metric, 
                   levels = c("Normalised Mean Error",
                              "Std. Dev. Difference",
                              "Mean Bias Error",
                              "R^2",
                              "Correlation Coeff."),
                   labels = c(expression("Normalised~Mean~Error"),
                              expression("Std.~Dev.~Difference"),
                              expression("Mean~Bias~Error"),
                              expression("R^2"),
                              expression("Correlation~Coeff.")))

# Correlation test
cor.test(LEdf$Value[LEdf$Model=="Current Environment (SAM)" & LEdf$Metric=="R^2"],WorldClimMetrics$AnnualPPT, method = "spearman")
cor.test(LEdf$Value[LEdf$Model=="Environmental Memory (SAM)" & LEdf$Metric=="R^2"],WorldClimMetrics$AnnualPPT, method = "spearman")
cor.test(LEdf$Value[LEdf$Model=="Biological Memory (AR1)" & LEdf$Metric=="R^2"],WorldClimMetrics$AnnualPPT, method = "spearman")


cor.test(LEdf$Value[LEdf$Model=="Current Environment (SAM)" & LEdf$Metric=="R^2" & LEdf$Site %in% c("AU-ASM","AU-DaS","AU-Dry","AU-How","AU-Stp","AU-TTE")],WorldClimMetrics$AnnualPPT[WorldClimMetrics$Transect=="NATT"], method = "spearman")
cor.test(LEdf$Value[LEdf$Model=="Environmental Memory (SAM)" & LEdf$Metric=="R^2" & LEdf$Site %in% c("AU-ASM","AU-DaS","AU-Dry","AU-How","AU-Stp","AU-TTE")],WorldClimMetrics$AnnualPPT[WorldClimMetrics$Transect=="NATT"], method = "spearman")
cor.test(LEdf$Value[LEdf$Model=="Biological Memory (AR1)" & LEdf$Metric=="R^2" & LEdf$Site %in% c("AU-ASM","AU-DaS","AU-Dry","AU-How","AU-Stp","AU-TTE")],WorldClimMetrics$AnnualPPT[WorldClimMetrics$Transect=="NATT"], method = "spearman")


cor.test(LEdf$Value[LEdf$Model=="Current Environment (SAM)" & LEdf$Metric=="R^2" & LEdf$Site %in% c("AU-Cpr","AU-Cum","AU-Gin","AU-GWW","AU-Tum","AU-Whr","AU-Wom")],WorldClimMetrics$AnnualPPT[WorldClimMetrics$Transect=="SAWS"], method = "spearman")
cor.test(LEdf$Value[LEdf$Model=="Environmental Memory (SAM)" & LEdf$Metric=="R^2" & LEdf$Site %in% c("AU-Cpr","AU-Cum","AU-Gin","AU-GWW","AU-Tum","AU-Whr","AU-Wom")],WorldClimMetrics$AnnualPPT[WorldClimMetrics$Transect=="SAWS"], method = "spearman")
cor.test(LEdf$Value[LEdf$Model=="Biological Memory (AR1)" & LEdf$Metric=="R^2" & LEdf$Site %in% c("AU-Cpr","AU-Cum","AU-Gin","AU-GWW","AU-Tum","AU-Whr","AU-Wom")],WorldClimMetrics$AnnualPPT[WorldClimMetrics$Transect=="SAWS"], method = "spearman")



NEE_SDScaledSensitivity_RTPV = function(Sites,Vars = c("Tair","Fsd","VPD","PPTshort","PPTlong","PPT"),Metric = "AnnualPPT"){
  
  # Source packages needed
  library(lubridate)
  library(magrittr)
  library(dplyr)
  library(coda)
  
  # Make the title nice
  source("functions/FindMetric.R")
  TitleUnits = FindMetric(Metric)
  
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
  ESen = data.frame("Site"=rep(Sites,each = 6),
                    "Variable" = rep(rownames(eval(as.name(Sites[1]))$ESen),length(Sites)),
                    "Low" = unlist(lapply(Sites,function(x) eval(as.name(x))$ESen$ESenLow)),
                    "Med" = unlist(lapply(Sites,function(x) eval(as.name(x))$ESen$ESenMedian)),
                    "High" = unlist(lapply(Sites,function(x) eval(as.name(x))$ESen$ESenHigh)))
  # Check whether the climate variable is significant - does the CI contain 0?
  ESen$Significant = sign(ESen$Low*ESen$High)==1
  
  
  # Weighted ESen
  ESen[ESen$Variable=="Tair",3:5] = ESen[ESen$Variable=="Tair",3:5]*Weights$Tair
  ESen[ESen$Variable=="Fsd",3:5] = ESen[ESen$Variable=="Fsd",3:5]*Weights$Fsd
  ESen[ESen$Variable=="VPD",3:5] = ESen[ESen$Variable=="VPD",3:5]*Weights$VPD
  ESen[ESen$Variable=="PPTshort",3:5] = ESen[ESen$Variable=="PPTshort",3:5]*Weights$PPTshort
  ESen[ESen$Variable=="PPTlong",3:5] = ESen[ESen$Variable=="PPTlong",3:5]*Weights$PPTlong
  ESen[ESen$Variable=="PPT",3:5] = ESen[ESen$Variable=="PPT",3:5]*Weights$PPT
  
  # Limit the dataframe to the variables requested
  ESen = ESen[ESen$Variable %in% Vars,]
  
  # Rename variables to nice names
  ESen$Variable[ESen$Variable == "Fsd"] = "Shortwave Radiation"
  ESen$Variable[ESen$Variable == "Tair"] = "Air Temperature"
  ESen$Variable[ESen$Variable == "PPTshort"] = "Short-term Precipitation"
  ESen$Variable[ESen$Variable == "PPTlong"] = "Long-term Precipitation"
  ESen$Variable[ESen$Variable == "PPT"] = "All Precipitation"
  
  # Assign levels to Variable
  levels = c("Shortwave Radiation",
             "Air Temperature",
             "VPD",
             "Short-term Precipitation",
             "Long-term Precipitation",
             "All Precipitation")
  ESen$Variable = factor(ESen$Variable,levels)
  
  # Order sites by chosen metric
  load("site_data/SiteMetrics_worldclim_0.5res.Rdata")
  metric = WorldClimMetrics[WorldClimMetrics$Sites %in% Sites,c("Sites",Metric)]
  colnames(metric) = c("Sites","Metric")
  # SiteOrder = paste0(metric[order(metric[,2]),1]," - ",metric[order(metric[,2]),2],TitleUnits$Units)
  # ESen$Site = paste0(ESen$Site,
  #                    " - ",
  #                    rep(metric[unique(match(ESen$Site,metric$Sites)),2],
  #                        each = nrow(ESen)/length(Sites)),TitleUnits$Units)
  # ESen$Site = factor(ESen$Site,levels=SiteOrder)
  ESen$Site = factor(ESen$Site, levels = metric[order(metric[,2]),1])
  
  # Plot the sensitivity covariates
  library(ggplot2)
  library(viridisLite)
  ESenPlot = ggplot(ESen) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_pointrange(aes(x = Site, 
                        ymin = Low, 
                        y = Med, 
                        ymax = High, 
                        color = Significant),
                    size = 0.25) +
    geom_errorbar(aes(x = Site, 
                      ymin = Low, 
                      ymax = High, 
                      color = Significant), 
                  width = 0.5,
                  size = 1) +
    facet_wrap(Variable~.,
             #   scales = "free_y",
             # ncol = (length(Vars)>=4)*2+(length(Vars)<4*1)) +
             ncol = 1) +
    #facet_grid(Variable~.) +
    scale_color_viridis_d(name="",
                          labels=c("FALSE"="Non-Significant",
                                   "TRUE"="Significant"),
                          guide="legend",
                          begin=0,
                          end=0.75,
                          direction = -1) +
    ylab("Sensitivity") +
    xlab("") +
    theme_bw() +
    theme(legend.position = "none",
          text = element_text(size=20),
          axis.text.x = element_text(angle=45, hjust=1),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.major.x = element_blank(),
          plot.title = element_text(hjust=0.5)) +
    ggtitle("NEE")
}

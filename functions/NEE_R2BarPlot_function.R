NEE_R2BarPlot_RTPV = function(Sites,Transects,Metric,Clusters = 0){

  
  # Load the required packages
  library(tidyverse)
  library(lubridate)
  library(viridis)
  library(ggpubr)
  library(gridExtra)
  
  # Make a nice metric title
  Titles =c("Annual Mean Temp",
            "Mean Diurnal Range",
            "Isothermality",
            "Temp Seasonality",
            "Max Temp of Hottest Month",
            "Min Temp of Coldest Month",
            "Temp Annual Range",
            "Mean Temp in Wettest Qtr",
            "Mean Temp in Driest Qtr",
            "Mean Temp in Summer",
            "Mean Temp in Winter",
            "Mean Annual PPT",
            "PPT in Wettest Month",
            "PPT in Driest Month",
            "PPT Seasonality",
            "PPT in Wettest Qtr",
            "PPT in Driest Qtr",
            "PPT in Summer",
            "PPT in Winter")
  
  PotMetric = c("AnnualMeanTemp",
                "MeanDiurnalRange",
                "Isothermality",
                "TempSeasonality",
                "MaxTempHotMon",
                "MinTempColdMon",
                "TempAnnualRange",
                "MeanTempWetQtr",
                "MeanTempDryQtr",
                "MeanTempHotQtr",
                "MeanTempColdQtr",
                "AnnualPPT",
                "PPTWetMon",
                "PPTDryMon",
                "PPTSeasonality",
                "PPTWetQtr",
                "PPTDryQtr",
                "PPTHotQtr",
                "PPTColdQtr")
  
  Title = Titles[Metric == PotMetric]

  # Initiliase R2 data frame
  R2 = data.frame("Site" = Sites,
                  "Transect" = Transects,
                  "R2.CUR" = 0,
                  "R2.SAM" = 0,
                  "R2.AR1" = 0,
                  "R2.KMN" = 0,
                  "R2.KMC" = 0)
  
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
    if (Clusters > 0){
    load(paste0("alternate/RTPV/results/NEE_output_",Clusters,"cluster_kmean_current_NDVI_RTPV_",Site,".Rdata"))
    R2$R2.KMN[R2$Site==Site] = output$r.squared
    
    load(paste0("alternate/RTPV/results/NEE_output_",Clusters,"cluster_kmean_current_RTPV_",Site,".Rdata"))
    R2$R2.KMC[R2$Site==Site] = output$r.squared
    } else {
      R2$R2.KMN = 0
      R2$R2.KMC = 0
    }
  }
  
  # Source worldclim correlations and climate metrics
  load("site_data/SiteMetrics_worldclim_0.5res.Rdata")
  
  # Plot the data!
  # Create the plot dataframe
  Site = rep(Sites,5)
  Site = factor(Site, levels = WorldClimMetrics[order(WorldClimMetrics[colnames(WorldClimMetrics)==Metric]),1])
  
  Transect = rep(Transects,5)
  
  # Order the models as we want
  Model = rep(c("Current Environment (k-means with no NDVI)",
                "Current Environment (k-means with NDVI)",
                "Current Environment (SAM)",
                "Environmental Memory (SAM)",
                "Biological Memory (AR1)"),
              each=length(Sites))
  Model = factor(Model,
                 levels=c("Current Environment (k-means with no NDVI)",
                          "Current Environment (k-means with NDVI)",
                          "Current Environment (SAM)",
                          "Environmental Memory (SAM)",
                          "Biological Memory (AR1)"))
  
  Value = c(R2$R2.KMC,
            R2$R2.KMN,
            R2$R2.CUR,
            R2$R2.SAM,
            R2$R2.AR1)
  
  Fig = data.frame(Site,
                   Transect,
                   Model,
                   Value)
  
  # Here we make sure that the bars are ordered by their size
  Fig = Fig[order(Fig$Value, decreasing = TRUE),]
  Fig$ValueFactor<-factor(Fig$Value, levels = unique(Fig$Value))
  
  # We then replace value with the difference so that the bars are "cumulative"
  for(site in Sites){
    Fig$Value[Fig$Site==site][1:4] = rev(diff(rev(Fig$Value[Fig$Site==site])))
  }
  
  # If we aren't plotting k-means, remove this data
  if (Clusters == 0){
    Fig = Fig[!(Fig$Model %in% c("Current Environment (k-means with no NDVI)",
                                "Current Environment (k-means with NDVI)")),]
  }
  
  # Plot for every site based on supplied metric
  Plot = ggplot(Fig,aes(fill=Model,y=Value,x=Site,group = ValueFactor)) +
    geom_bar(position="stack",stat="identity") +
    geom_bar(stat = "identity",color = "black",size = 1) +
    geom_point(data = Fig[Fig$Transect=="NATT",],aes(x=Site,y = rep(1,length(Site))),shape=8,size = 2,show.legend=FALSE) +
    geom_point(data = Fig[Fig$Transect=="NATT",],aes(x=Site,y = rep(-0.025,length(Site))),shape=8,size = 2,show.legend=FALSE) +
    scale_fill_viridis_d(direction=-1,begin=0,end=1) +
    guides(fill=guide_legend(nrow=ceiling(length(unique(Fig$Model))/3),byrow=TRUE)) +
    coord_flip(ylim=c(0,1)) +
    ylab(parse(text="R^2")) +
    xlab("") +
    theme_bw() +
    theme(legend.position = "none", 
          text = element_text(size = 20),
          legend.title = element_blank(),
          panel.grid.major.y = element_blank(),
          plot.title = element_text(hjust=0.5)) +
    ggtitle("NEE")
  Plot

}

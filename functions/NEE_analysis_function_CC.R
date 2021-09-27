NEE_analysis_current_RTPV <- function(Site){
  
  # A function to take the output from a R2jags model run for an OzFlux site and
  # turn it into something useful and interesting and possibly, hopefully, 
  # insightful. Fingers crossed!
  # 
  # INPUTS:
  # - Site: A character vector of length 6 containing the offical OzFlux code
  #             for the site, including the 'AU-' part
  # 
  # OUTPUTS:
  # - We'll have to wait and see what I come up with!
  # 
  # HERE WE GO!
  # 
  
  # ##################
  # Let's ACTIVATE!
  # ##################
  
  # Let the user know which site the function is looking at
  message("*** Analysing current R2jags output for ",Site," ***")
  # 
  # Load in the output data we are analysing
  # Look in folder "results" for the data
  File = list.files("output/RTPV/",pattern = paste0("NEE_current_output_RTPV_",Site))
  # Read the data into R - note that if multiple results are available for a 
  # site, we take the most recent
  message("File is ",File[length(File)])
  load(paste0("output/RTPV/",File[length(File)]))
  
  # Source the necessary packages
  library(coda)
  library(ggplot2)
  library(dplyr)
  library(rjags)
  library(R2jags)
  library(mcmcplots)
  library(lubridate)
  library(magrittr)
  library(zoo)
  source('functions/DBDA2E-utilities.R')
  
  # ##################
  # Convergence checks
  # ##################
  
  # List the "fundamental" parameters - e.g. those that are assigned priors and
  # are not a function of other parameters. Stochastic parameters? Maybe.
  stochastic.params = c("phi0",
                        "sig_y",
                        sprintf("an[%d]",seq(1:16)),
                        sprintf("ag[%d]",seq(1:16)),
                        sprintf("deltaXAP[%d]",seq(1:8)),
                        sprintf("deltaXA[%d,%d]",rep(1:4,10),rep(1:10,each=4))
                        )

  # Convert output to an mcmc object
  # Either take the object already saved as an mcmc object for the current 
  # workflows or, to maintain compatibility with older workflows, calculate it
  # from the rjags object
  if (class(output) == "list"){
    output.mcmc = output$output.mcmc
  }else{
    output.mcmc = as.mcmc.rjags(output)
  }
  # Make sure we can analyse this as an mcmc list
  output.mcmc = as.mcmc.list(output.mcmc)
  rm(output)
  
  # Produce plots of each parameter to assess convergence.
  for (param in stochastic.params){
    # Output the variable
    print(param)
    diagMCMC(output.mcmc,param,saveName = paste0(Site,"_current_",Sys.Date()))
    Sys.sleep(1)
    graphics.off()
  }
  
  message("Running Gelman Diagnostics for ",Site)
  # We find the Gelman diagnostic (it has a proper name but I'm a hack)
  # I think it's the shrink factor or something lol
  Gelman = gelman.diag(output.mcmc,multivariate=FALSE)
  # Find how many, and which, parameters fall outside of the acceptable limits
  # which here is set to 1.1
  Gelman.Fail = Gelman$psrf[Gelman$psrf[,2]>1.1,]
  # This picks up NA values for those parameters that are fixed to a certain 
  # value - we should exclude these
  Gelman.Fail = Gelman.Fail[complete.cases(Gelman.Fail),]
  
  message("Running Effective Sample Size for ",Site)
  # We find the effective sample size for each parameter
  ESS.raw = effectiveSize(output.mcmc)
  # Where parameters are forced to 0, then the ESS is also 0. Therefore we exclude
  # these when considering the fit. In general, higher ESS is better, with 10,000+
  # being ideal
  ESS = ESS.raw[ESS.raw > 0]
  # Plot a histogram to visualise how the ESS distribution breaks down
  ESSPlot = ggplot(data.frame(ESS)) +
            geom_histogram(aes(ESS),binwidth=250)
  # See which parameters are way below 10,000 ESS
  ESS.Fail = ESS[ESS<10000] # & names(ESS) %in% stochastic.params]
  
  message("Running Geweke Diagnostics for ",Site)
  # We calculate the Geweke diagnostic - this should fall within the confidence 
  # bounds of -2 and 2. 
  Geweke = geweke.diag(output.mcmc)
  # I think this is less important - or at least, it depends on the length of the
  # burn-in period
  # Count how many elements are outside the bounds
  GewekeCount = unlist(lapply(Geweke, function(i) sum((i$z>2 | i$z<(-2)), #& names(i$z) %in% stochastic.params,
                                                      na.rm=TRUE)))
  GewekeNames = (lapply(Geweke, function(i) names(i$z)[(i$z>2 | i$z<(-2))])) # & names(i$z) %in% stochastic.params]))
  Geweke.Fail = mean(GewekeCount)
  
  # ##################
  # Model Performance
  # ##################
  message("Running Model Performance for ",Site)
  
  # Load the observations
  name = paste0(Site,"_Input")
  load(paste0("inputs/RTPV/",name,"_RTPV.Rdata"))
  assign("obs",eval(as.name(name)))
  
  # Create dataframe of observed vs modelled with confidence intervals
#  NEE_pred = output$BUGSoutput$median$NEE_pred
#  NEE_pred_min = output$BUGSoutput$summary[substr(rownames(output$BUGSoutput$summary),1,3)=="NEE",3]
#  NEE_pred_max = output$BUGSoutput$summary[substr(rownames(output$BUGSoutput$summary),1,3)=="NEE",7]
  summary = summary(output.mcmc)
  NEE_pred = summary$statistics[substr(rownames(summary$statistics),1,5)=="NEE_p",1]
  NEE_pred_min = summary$quantiles[substr(rownames(summary$quantiles),1,5)=="NEE_p",1]
  NEE_pred_max = summary$quantiles[substr(rownames(summary$quantiles),1,5)=="NEE_p",5]
  NEE_obs = obs$NEE[-(1:365)]
  
  df = data.frame("Date"=obs$DailyData$TIMESTAMP[-(1:365)],
                  "Pred"=NEE_pred,
                  "Min"=NEE_pred_min,
                  "Max"=NEE_pred_max,
                  "Obs"=NEE_obs)
  
  # Plot the daily data
  ObsVsNEE_daily = ggplot(df) +
    geom_ribbon(aes(x=Date,ymin=Min,ymax=Max, fill="Pred"),alpha=0.5) +
    geom_line(aes(x=Date,y=Obs,color="Obs", fill="Obs")) +
    geom_line(aes(x=Date,y=Pred,color="Pred", fill="Obs")) +
    scale_color_viridis_d(name="Data",
                          labels=c("Obs"="Observations","Pred"="Predicted"),
                          guide="legend",
                          option="magma",
                          direction=-1,
                          begin=0.2,
                          end=0.8) +
    scale_fill_viridis_d(name="Data",
                         labels=c("Obs"="Observations","Pred"="Predicted"),
                         guide="legend",
                         option="magma",
                         direction=-1,
                         begin=0.2,
                         end=0.8) +
    theme_bw()
  
  # We also calculate a dataframe of moving averages to smooth out the plot
  # Set the window for moving average
  k = 15
  # Create the dataframe
  df_ma = data.frame("Date"=rollmedian(df$Date,k),
                    "Pred"=rollmedian(df$Pred,k),
                    "Min"=rollmedian(df$Min,k),
                    "Max"=rollmedian(df$Max,k),
                    "Obs"=rollmedian(df$Obs,k))
  
  # Plot the moving average data
  ObsVsNEE_ma = ggplot(df_ma) +
    geom_ribbon(aes(x=Date,ymin=Min,ymax=Max, fill="Pred"),alpha=0.5) +
    geom_line(aes(x=Date,y=Obs,color="Obs", fill="Obs")) +
    geom_line(aes(x=Date,y=Pred,color="Pred", fill="Obs")) +
    scale_color_viridis_d(name="Data",
                          labels=c("Obs"="Observations","Pred"="Predicted"),
                          guide="legend",
                          option="magma",
                          direction=-1,
                          begin=0.2,
                          end=0.8) +
    scale_fill_viridis_d(name="Data",
                         labels=c("Obs"="Observations","Pred"="Predicted"),
                         guide="legend",
                         option="magma",
                         direction=-1,
                         begin=0.2,
                         end=0.8) +
    theme_bw()

  # Since the daily data is likely to be very noisy, aggregate into monthly
  # data with sums 
  # I HAVE NO IDEA IF SUMMING CONFIDENCE INTERVALS IS LEGIT
  df$year = year(df$Date)
  df$month = month(df$Date)
  
  df_monthly = df %>%
    group_by(year,month) %>%               
    summarise(Pred=sum(Pred),
              Min=sum(Min),
              Max=sum(Max),
              Obs=sum(Obs))
  
  df_monthly$Date = as.yearmon(paste(df_monthly$year, df_monthly$month), "%Y %m")
  # Plot monthly data
  ObsVsNEE_monthly = ggplot(df_monthly) +
    geom_ribbon(aes(x=Date,ymin=Min,ymax=Max, fill="Pred"),alpha=0.5) +
    geom_line(aes(x=Date,y=Obs,color="Obs", fill="Obs")) +
    geom_line(aes(x=Date,y=Pred,color="Pred", fill="Obs")) +
    scale_color_viridis_d(name="Data",
                          labels=c("Obs"="Observations","Pred"="Predicted"),
                          guide="legend",
                          option="magma",
                          direction=-1,
                          begin=0.2,
                          end=0.8) +
    scale_fill_viridis_d(name="Data",
                         labels=c("Obs"="Observations","Pred"="Predicted"),
                         guide="legend",
                         option="magma",
                         direction=-1,
                         begin=0.2,
                         end=0.8) +
    theme_bw()
  
  # Calculate the r squared value for the SAM model
  CUR.R2 = summary(lm(NEE_pred ~ NEE_obs))$r.squared
  Phi = summary$statistics["phi0",]
  # Mean Bias Error
  CUR.MBE = sum(NEE_pred-NEE_obs,na.rm=TRUE)/length(NEE_pred)
  # Normalised Mean Error
  CUR.NME = sum(abs(NEE_pred-NEE_obs),na.rm=TRUE)/sum(abs(mean(NEE_obs,na.rm=TRUE)-NEE_obs),na.rm=TRUE)
  # Standard Deviation Difference
  CUR.SDD = abs(1-sd(NEE_pred,na.rm=TRUE)/sd(NEE_obs,na.rm=TRUE))
  # Correlation Coefficient
  CUR.CCO = cor(NEE_pred,NEE_obs,use = "complete.obs", method = "pearson")
  
  # Extract climate sensitivities
  
  # Define params
  SensitivityCovariates = c(sprintf("ESen[%d]",seq(1:7)))
  # Extract 2.5%, median and 97.5%
  ESenMedian = summary$statistics[rownames(summary$statistics) %in% SensitivityCovariates,1]
  ESenLow = summary$quantiles[rownames(summary$quantiles) %in% SensitivityCovariates,1]
  ESenHigh = summary$quantiles[rownames(summary$quantiles) %in% SensitivityCovariates,5]
  # Place in dataframe
  ESen = data.frame(ESenLow,ESenMedian,ESenHigh)
  rownames(ESen) = c("Tair",
                     "Fsd",
                     "VPD",
                     "PPTshort",
                     "PPTlong",
                     "Precip")
  
  # Extract cumulative weights
  # Define params
  CumWeightParams = c(sprintf("cum_weightA[%d,%d]",rep(1:5,14),rep(1:14,each=5)),
                      sprintf("cum_weightAP[%d]",seq(1:8)))
  # Extract 2.5%, median and 97.5%
  WeightsMedian = summary$statistics[rownames(summary$statistics) %in% CumWeightParams,1]
  WeightsLow = summary$quantiles[rownames(summary$quantiles) %in% CumWeightParams,1]
  WeightsHigh = summary$quantiles[rownames(summary$quantiles) %in% CumWeightParams,5]
  # Place in dataframe
  CumWeights = data.frame(WeightsLow,WeightsMedian,WeightsHigh)
  # rownames(CumWeights) = c(rep("Tair",10),
  #                    rep("Fsd",10),
  #                    rep("VPD",10),
  #                    rep("curSWC",10),
  #                    rep("antSWC",10),
  #                    rep("Precip",8))
  
  # Create a nice output and save it
  output = list("Gelman.Fail" = Gelman.Fail,
                "ESS.Fail" = ESS.Fail,
                "Geweke.Fail" = Geweke.Fail,
                "CUR.R2" = CUR.R2,
                "CUR.MBE" = CUR.MBE,
                "CUR.NME" = CUR.NME,
                "CUR.SDD" = CUR.SDD,
                "CUR.CCO" = CUR.CCO,
                "df" = df,
                "ESen" = ESen,
                "CumWeights" = CumWeights,
                "Phi0" = Phi)
  
  save(output,file = paste0("NEE_current_analysis_RTPV_",Site,"_",Sys.Date(),".Rdata"))
}


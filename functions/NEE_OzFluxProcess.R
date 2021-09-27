OzFluxProcess_noSWC = function(Site){

  # A function to extract the necessary data from the daily OzFlux netcdf for a 
  # given site. This is designed to work with the model from Liu et al, 2019. 
  # To work correctly, the OzFlux netcdf for the site must be saved in a
  # subfolder named as Site_raw_data. 
  # 
  # The below variables are extracted, quality-checked and return in a nice, 
  # usable fashion:
  # 
  #   NEE_LL          Net Ecosystem Exchange
  #   SW_IN_F         Incoming Shortwave Radiation
  #   TA_F            Air Temperature
  #   VPD_F           Vapour Pressure Deficit
  #   SWC_F_MDS_1     Top-layer Soil Moisture Content
  #   P_F             Precipitation
  # 
  # 
  # ############################################################################
  # Function inputs and outputs
  # ############################################################################
  #
  # INPUTS:
  #  - Site: A character vector with the FluxNet siteID
  #  
  #  OUTPUTS:
  #  - "AU-Site_Input_noSWC.Rdata": a .Rdata file containing a list of the same name. 
  #                           The list includes all necessary inputs for the 
  #                           NEEModel.R workflows as well as a dataframe of
  #                           daily data and a list of QC metrics

  # Let the user know which site the function is looking at
  message("*** Extracting data for ",Site," ***")
  
  # ##################
  # Extract raw data
  # ##################
  # Source the ncdf package
  library(ncdf4)
  # Load in the OzFlux data we need:
  # Look in folder "Site_raw_data" for the data
  File = list.files("Site_data",pattern = Site)
  # Read the data into R 
  NCDF = nc_open(paste0("Site_data/",File))
  # Change timestamps into dates
  TIMESTAMP = ncvar_get(NCDF,"time")
  # Convert time into datetimes
  # Find the origin of the times
  time_from = substr(ncatt_get(NCDF, "time")$units, 12, 30)
  # times are in days since origin
  Days_to_Secs = 24*3600
  TIMESTAMP = as.POSIXct(TIMESTAMP*Days_to_Secs, origin=time_from, tz = ncatt_get(NCDF, 0)$time_zone)

  # List the variables we want to extract as well as their quality control
  # Note some sites are using L5 PI data as recommended by site PIs - as such
  # we must check and take different variables for this data
  if (substr(File,8,9)=="L5"){
  Variables = c("NEE",
                "Fsd",
                "Ta",
                "VPD",
                "Precip")
  } else {
    Variables = c("Fc",
                  "Fsd",
                  "Ta",
                  "VPD",
                  "Precip")
  }
  # Extract the variables we require
  Data = data.frame(TIMESTAMP)
  for (Var in Variables){
    Data[Var] = ncvar_get(NCDF,Var)
    QC = paste0(Var,"_QCFlag")
    Data[QC] = ncvar_get(NCDF,QC)
  }
  
  # Give common names to NEE/Fc columns
  colnames(Data)[2] = "NEE"
  colnames(Data)[3] = "NEE_QCFlag"
  
  # ##################
  # Quality check data
  # ##################
  
  # A QC flag that is not mod10 is bad data - i.e. should be masked
  # A QC flag of greater than 0 but mod10 is gap-filled in a robust manner - but
  # still gap-filled! 
  
  # Identify QC columns
  QCcols = grep("QC",colnames(Data))
  
  # First, we deal with BAD data - i.e data that is not robust/reliable
  # Remove first row if any data is poor - repeat as necessary
  count = 0
  while(any(Data[1,QCcols]%%10!=0)){
    Data = Data[-1,]
    count = count + 1
  }
  if (count > 0){
  message("Info! ",count," rows removed from start of record due to poor data (QC flag not mod 10).")
  }
  # Remove last row if any data is poor - repeat as necessary
  count = 0
  while(any(Data[nrow(Data),QCcols]%%10!=0)){
    Data = Data[-nrow(Data),]
    count = count + 1
  }
  if (count > 0){
  message("Info! ",count," rows removed from end of record due to poor data (QC flag not mod 10).")
  }
  
  # Find any remaining bad data and mask as NA:
  QC = Data[,QCcols]%%10 != 0
  # Find the corresponding rows
  BadDataRows = which(apply(QC,MARGIN=1,function(x) any(x == TRUE)) %in% TRUE)
  # Create outputs for QC reference
  BadObs = sum(QC)
  Bad30Min = length(BadDataRows)
  
  # Calculate percentage of data that is bad.
  PercentBad30Min = round(Bad30Min*100/nrow(Data),digits = 0)
  PercentBadObs = round(100*BadObs/(nrow(Data)*5),digits = 0)
  
  if (BadObs>0){
    message("Info! ",
            BadObs,
            " individual observations (",
            PercentBadObs,
            "% of all observations) are bad data and will be masked as NA!")
    message("These observations occur in ",
            Bad30Min,
            " different 30-min periods, equal to ",
            PercentBad30Min,
            "% of all 30 min periods")
  } else {
    message("Info! All data are measured or good quality gap-filled!")
  }
  ## NOTE: The below loops are quite slow... attempt to vectorise if possible
  #for all rows with bad data
  for(i in BadDataRows){
    # for each QC flag
    for (j in QCcols){
      # if the flag indicates bad data
      if (Data[i,j] %%10 !=0){
        # set the data to NA
        Data[i,(j-1)] = NA
      }
    }
  }
  
  # Now check for any data that is GOOD but GAP-FILLED i.e. > 0 && mod10 = 0
  
  # Arbitrarily decide that less than 95% measured/good data is 
  # worrying
  # Count rows that have one or more gap-fill QC flag
  QC = sum(apply(Data[,QCcols],MARGIN=1,function(x) any(x > 0 && x%%10 == 0)))
  GapFilledData = Data[,QCcols] > 0 & Data[,QCcols]%%10 == 0  
  # Calculate percentage of data that is gap-filled.
  PercentFilled30Min = round(QC*100/nrow(Data),digits = 0)
  PercentFilledObs = round(100*sum(GapFilledData)/(nrow(Data)*5),digits = 0)
  # Create outputs for QC reference
  FilledObs = sum(GapFilledData)
  Filled30Min = QC
  # If more than 5% of the data is poor, print a warning
  if (PercentFilledObs > 0){
    message("Info! ",
            PercentFilled30Min,
            "% of the 30-minute data intervals contain at least one gap-filled observation!")
    message("This is ",
            FilledObs,
            " individual observations, equivalent to ",
            PercentFilledObs,
            "% of all observations.")
  } else {
    "All data are good measurements! Rejoice!"
  }
  # Check for excessive consecutive streaks of poor data
  # Find the sequences of poor/good data
  Seq = rle(apply(Data[,QCcols],MARGIN=1,function(x) any(x != 0)))
  # Find the lengths of these sequences for the poor data
  Lengths = Seq$lengths[Seq$values==TRUE]
  # If a run of 5 or more days of poor data exists, print a warning
  if (length(Lengths)>0){
    if (max(Lengths)>=5){
      message("Info! There is a run of ",
                 max(Lengths),
                 " consecutive half-hours with at least one gap-filled observation!")
    }
  }
  
  # #####################
  # Retime data to daily
  # #####################
  # Source required packages
  library(lubridate)
  library(magrittr)
  library(tidyverse)
  library(zoo)
  
  # Create dataframe of daily values and QC counts
  Data_day <- Data %>%
    mutate(TIMESTAMP=as.Date(TIMESTAMP, 
                             format="%Y-%m-%d %H:%M:%S", 
                             tz = ncatt_get(NCDF, 0)$time_zone)) %>%
    group_by(TIMESTAMP) %>%               # group by the day column
    # Count the amount of bad and gap-filled observations each day
    summarise(NEE=mean(NEE),
              NEE_BD=sum(NEE_QCFlag%%10 != 0),
              NEE_GF = sum(NEE_QCFlag%%10 == 0 & NEE_QCFlag > 0),
              Fsd=mean(Fsd),
              Fsd_BD=sum(Fsd_QCFlag%%10 != 0),
              Fsd_GF = sum(Fsd_QCFlag%%10 == 0 & Fsd_QCFlag > 0),
              Ta=mean(Ta),
              Ta_BD=sum(Ta_QCFlag%%10 != 0),
              Ta_GF = sum(Ta_QCFlag%%10 == 0 & Ta_QCFlag > 0),
              VPD=mean(VPD),
              VPD_BD=sum(VPD_QCFlag%%10 != 0),
              VPD_GF = sum(VPD_QCFlag%%10 == 0 & VPD_QCFlag > 0),
              Precip=sum(Precip),
              Precip_BD=sum(Precip_QCFlag%%10 != 0),
              Precip_GF = sum(Precip_QCFlag%%10 == 0 & Precip_QCFlag > 0))
  
  # ####################
  # Trim to full years
  # ####################
  
  # First remove any days at the beginning or end of the record with more than
  # 6 instances (12.5%) of bad data
  # Remove first row if too much data is poor - repeat as necessary
  # Define QC cols
  QCcols_day = seq(3,16,by=3)
  count_daystart = 0
  while(any(Data_day[1,QCcols_day]>6)){
    Data_day = Data_day[-1,]
    count_daystart = count_daystart + 1
  }
  if (count_daystart > 0){
    message("Info! ",
            count_daystart,
            " days removed from start of daily record due to more than 12.5% bad data existing")
  }
  # Remove last row if too much data is poor - repeat as necessary
  count_dayend = 0
  while(any(Data_day[nrow(Data_day),QCcols_day]>6)){
    Data_day = Data_day[-nrow(Data_day),]
    count_dayend = count_dayend + 1
  }
  if (count_dayend > 0){
    message("Info! ",
            count_dayend,
            " days removed from end of daily record due to more than 12.5% bad data existing")
  }
  
  #
  DaysRemoved = count_daystart+count_dayend
  # Find the years that are complete in the daily data
  AllYears = unique(year(Data_day$TIMESTAMP))
  FullYears = rle(year(Data_day$TIMESTAMP))$values[rle(year(Data_day$TIMESTAMP))$lengths > 364]
  CutYears = AllYears[!(AllYears %in% FullYears)]
  
  # Trim the daily data to these years
  Data_day = Data_day[year(Data_day$TIMESTAMP) %in% FullYears,]
  
  # Report the cut years
  message("Info! The years ",paste(shQuote(CutYears), collapse=", "),
          " have been cut for being incomplete!")
  message(FullYears[length(FullYears)]+1-FullYears[1],
          " full years remain between ",
          FullYears[1],
          " to ",
          FullYears[length(FullYears)])
  
  
  
  # ####################
  # Create inputs required for modelling
  # ####################
  
  # First, create the fixed parameters:
  
  # Number of short-term climatic predictors
  # Tave, SW, VPD, PPTshort
  Nv = 4
  # Days into past considered for short-term predictors
  Nlag = 14
  # Time lag for PPTlong(number of different time periods)
  NlagP = 8
  # Total number of climate covariates (see paper for info)
  Ns = 16
  # Number of blocks for PPTlong
  NblocksP = NlagP
  # Time blocks i.e. how the Nlag days are grouped together
  block = c(1:7, rep(8:9, each = 2), rep(10, 3))
  # The size of the block that each day is included in
  BlockSize = c(rep(1, 7), rep(2, 4),rep(3, 3))
  # Number of blocks
  Nblocks = max(block)
  
  # Calculate other parameters:
  
  # Number of days that antecedent conditions can be calculated for
  Nmem = nrow(Data_day)-365
  # Indices of days that can be modelled
  Mem_records = 366:nrow(Data_day)
  
  # Create the climate predictor matrix
  # See Model_Liu inputs for correct order
  # SWC is repeated to account for current and antecedent.
  clim = matrix(c(Data_day$Ta,
                  Data_day$Fsd,
                  Data_day$VPD,
                  Data_day$Precip),
                ncol = Nv)
  # Fill any climate NAs by interpolating over the gap
  clim = na.approx(clim)
  # Mean centre the climatic variables
  clim = scale(clim,scale=FALSE)
  
  # Create the NEE vector
  NEE = Data_day$NEE
  
  ## NDVI
  # NDVI was extracted from Google Earth Engine on 17/09/20 at 12:00 
  # Product was MCD43A3.006 MODIS Albedo Daily 500m
  load("inputs/VegIndex.Rdata")
  NDVI_df = VegIndex[VegIndex$site==Site,]
  NDVI_df = NDVI_df[as.Date(NDVI_df$date) %in% as.Date(Data_day$TIMESTAMP),c("date","ndvi_sg")]
  NDVI_df$date = as.Date(NDVI_df$date)
  # Check for missing data
  for (i in Data_day$TIMESTAMP[!(Data_day$TIMESTAMP %in% NDVI_df$date)]){
    NewDate = as.Date(i,origin="1970-01-01")
    NewNDVI = mean(c(NDVI_df$ndvi_sg[NDVI_df$date==(i-1)],NDVI_df$ndvi_sg[NDVI_df$date==(i+1)]))
    df = data.frame("date" = NewDate, "ndvi_sg" = NewNDVI)
    NDVI_df = rbind(NDVI_df,df) 
  }
  NDVI_df = NDVI_df[order(NDVI_df$date),]
  
  NDVI = NDVI_df$ndvi_sg
  
  # Scale NDVI to vary from 0 to 1 - lowest NDVI is "not growing season" and
  # highest is "growing season"
  NDVI = (NDVI-min(NDVI))/(max(NDVI)-min(NDVI))
  
  ## PPTlong
  
  # Source the ppt processing function
  source("PPTProcess.R")
  # Extract the ppt matrix
  ppt_multiscale = PPTProcess(Data_day)
  # Mean centre ppt
  ppt_multiscale = scale(ppt_multiscale,scale=FALSE)
  # Note this gives precip values of less than 0 which is intriguing.
  
  ## Create QC list for reference
  QCList = list("BadObs" = BadObs,
                "Bad30Min" = Bad30Min,
                "PercentBadObs" = PercentBadObs,
                "PercentBad30Min" = PercentBad30Min,
                "FilledObs" = FilledObs,
                "Filled30Min" = Filled30Min,
                "PercentFilledObs" = PercentFilledObs,
                "PercentFilled30Min" = PercentFilled30Min,
                "DaysRemoved" = DaysRemoved,
                "CutYears" = CutYears)
  
  
  # ###############
  # Create output list
  # ###############
  
  output = list("Nv"=Nv,
                "Ns"=Ns,
                "Nlag"=Nlag,
                "Nmem"=Nmem,
                "NlagP"=NlagP,
                "Mem_records"=Mem_records,
                "clim"=clim,
                "ppt_multiscale"=ppt_multiscale,
                "NEE"=NEE,
                "NDVI"=NDVI,
                "NblocksP" = NblocksP,
                "block" = block,
                "BlockSize" = BlockSize,
                "Nblocks" = Nblocks,
                "DailyData" = Data_day,
                "QCList" = QCList)
  name = paste0(Site,"_input")
  assign(name,output)
  save(list=c(name),file=paste0(name,"_RTPV.Rdata"))
  
  
##### STILL TO BE DONE - Check QC for L6 data
}
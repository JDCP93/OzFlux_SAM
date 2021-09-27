CABLEoutput = function(Site){
  # Let the user know which site the function is looking at
  message("*** Extracting CABLE output data for ",Site," ***")
  
  # ##################
  # Extract raw data
  # ##################
  # Source the ncdf package
  library(ncdf4)
  # Load in the OzFlux data we need:
  # Look in folder "Site_raw_data" for the data
  File = list.files("CABLE",pattern = Site)
  # Read the data into R 
  NCDF = nc_open(paste0("CABLE/",File))
  # Change timestamps into dates
  TIMESTAMP = ncvar_get(NCDF,"time")
  # Convert time into datetimes
  # Find the origin of the times
  time_from = substr(ncatt_get(NCDF, "time")$units, 15, 33)
  # times are in days since origin
  Days_to_Secs = 24*3600
  TIMESTAMP = as.POSIXct(TIMESTAMP, origin=time_from, tz = "UTC")

  # We want the NEE and LE

  Variables = c("NEE",
                "Qle")
  Data = data.frame(TIMESTAMP)
  for (Var in Variables){
    Data[Var] = ncvar_get(NCDF,Var)
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
                             tz = "UTC")) %>%
    group_by(TIMESTAMP) %>%               # group by the day column
    # Count the amount of bad and gap-filled observations each day
    summarise(NEE=mean(NEE),
              LE=mean(Qle))
  # ###############
  # Create output list
  # ###############
  
  output = Data_day
  save(output,file=paste0("CABLE/processed/",Site,"_extracted.Rdata"))
  
  name = paste0(Site,"_Input")
  load(paste0("inputs/RTPV/",name,"_RTPV.Rdata"))
  assign("NEEobs",eval(as.name(name)))
  
  name = paste0(Site,"_input")
  load(paste0("inputs/RTPV/LE/",name,"_LE_RTPV.Rdata"))
  assign("LEobs",eval(as.name(name)))
  
  Data = merge(NEEobs$DailyData,Data_day,by.x="TIMESTAMP",by.y="TIMESTAMP") %>%
          merge(LEobs$DailyData,by.x="TIMESTAMP",by.y="TIMESTAMP") %>%
          select(TIMESTAMP,"NEE_obs"=NEE.x,"NEE_pred"=NEE.y,"LE_obs"=Fe,"LE_pred"=LE)
  
  ggplot(Data) + geom_line(aes(x=TIMESTAMP,y=NEE_obs),color="blue") + geom_line(aes(x=TIMESTAMP,y=NEE_pred),color="red")
  
  output = list()
  output$CABLE.NEE.R2 = summary(lm(Data$NEE_pred ~ Data$NEE_obs))$r.squared
  output$CABLE.NEE.MBE = sum(Data$NEE_pred-Data$NEE_obs,na.rm=TRUE)/length(Data$NEE_pred)
  output$CABLE.NEE.NME = sum(abs(Data$NEE_pred-Data$NEE_obs),na.rm=TRUE)/sum(abs(mean(Data$NEE_obs,na.rm=TRUE)-Data$NEE_obs),na.rm=TRUE)
  output$CABLE.NEE.SDD = abs(1-sd(Data$NEE_pred,na.rm=TRUE)/sd(Data$NEE_obs,na.rm=TRUE))
  output$CABLE.NEE.CCO = cor(Data$NEE_pred,Data$NEE_obs,use = "complete.obs", method = "pearson")
  
  output$CABLE.LE.R2 = summary(lm(Data$LE_pred ~ Data$LE_obs))$r.squared
  output$CABLE.LE.MBE = sum(Data$LE_pred-Data$LE_obs,na.rm=TRUE)/length(Data$LE_pred)
  output$CABLE.LE.NME = sum(abs(Data$LE_pred-Data$LE_obs),na.rm=TRUE)/sum(abs(mean(Data$LE_obs,na.rm=TRUE)-Data$LE_obs),na.rm=TRUE)
  output$CABLE.LE.SDD = abs(1-sd(Data$LE_pred,na.rm=TRUE)/sd(Data$LE_obs,na.rm=TRUE))
  output$CABLE.LE.CCO = cor(Data$LE_pred,Data$LE_obs,use = "complete.obs", method = "pearson")
  
  save(output,file=paste0("CABLE/processed/",Site,"_metrics.Rdata"))
}
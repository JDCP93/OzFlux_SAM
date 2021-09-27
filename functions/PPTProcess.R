PPTProcess = function(Data){
  
  # A function to sit within FluxNetProcess. The function takes the current data
  # from the FluxNetProcess function and extracts mean antecedent precipitation
  # for each day from day 366 onwards.
  # 
  # The ppt_multiscale output calculates the mean daily rainfall over past time
  # windows for each day k, corresponding to row k. Column j for row k contains
  # the mean daily rainfall over the below days:
  #   j=1 - Set to 0 within model code so not required (set to 0 here)
  #   j=2 - k-14 to k-20
  #   j=3 - k-21 to k-29
  #   j=4 - k-30 to k-59
  #   j=5 - k-60 to k-119
  #   j=6 - k-120 to k-179
  #   j=7 - k-180 to k-269
  #   j=8 - k-270 to k-365
  # 
  # ############################################################################
  # Function inputs and outputs
  # ############################################################################
  #
  # INPUTS:
  #  - Data: A working dataframe from within FluxNetProcess
  #  
  #  OUTPUTS:
  #  - ppt_multiscale: A nrow(Data)x8 matrix of antecedent rainfall calculated 
  #                    as daily means over varying time windows 
  

  # Extract the precipitation data from the input
  PPT = unname(unlist(Data[,7]))
  # Initiliase the precipitation dataframe
  ppt_multiscale = data.frame()

  # for each day in the main data series, calculate the antecedent rainfall
  for (i in 366:nrow(Data)){
  
    ppt_multiscale[i,1] = 0
    ppt_multiscale[i,2] = mean(PPT[(i-14):(i-20)],na.rm=TRUE)
    ppt_multiscale[i,3] = mean(PPT[(i-21):(i-29)],na.rm=TRUE)
    ppt_multiscale[i,4] = mean(PPT[(i-30):(i-59)],na.rm=TRUE)
    ppt_multiscale[i,5] = mean(PPT[(i-60):(i-119)],na.rm=TRUE)
    ppt_multiscale[i,6] = mean(PPT[(i-120):(i-179)],na.rm=TRUE)
    ppt_multiscale[i,7] = mean(PPT[(i-180):(i-269)],na.rm=TRUE)
    ppt_multiscale[i,8] = mean(PPT[(i-270):(i-365)],na.rm=TRUE)
  } 
  
  # Turn into matrix
  ppt_multiscale = as.matrix(ppt_multiscale)
}
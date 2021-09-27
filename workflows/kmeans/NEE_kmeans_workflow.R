
# Tidy up
rm(list=ls())
# List sites
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

############################################################################
# All sites, variable k amount with optimisation
############################################################################


# # Source the function
# source("alternate/RTPV/NEE_variable_kmean_function_current_RTPV.R")
# source("alternate/RTPV/NEE_variable_kmean_function_current_NDVI_RTPV.R")
# source("alternate/RTPV/NEE_variable_kmean_function_singlePPT_RTPV.R")
# source("alternate/RTPV/NEE_variable_kmean_function_singlePPT_NDVI_RTPV.R")
# source("alternate/RTPV/NEE_variable_kmean_function_allPPT_RTPV.R")
# source("alternate/RTPV/NEE_variable_kmean_function_allPPT_NDVI_RTPV.R")
# source("alternate/RTPV/NEE_variable_kmean_function_alllags_RTPV.R")
# source("alternate/RTPV/NEE_variable_kmean_function_alllags_NDVI_RTPV.R")
# 
# # Calculate the values for the variable functions
# for (Site in Sites){
#   NEE_current_variable_kmean_RTPV(Site)
#   NEE_current_NDVI_variable_kmean_RTPV(Site)
#   for (k in 2:8){
#     NEE_singlePPT_variable_kmean_RTPV(Site,k)
#     NEE_singlePPT_NDVI_variable_kmean_RTPV(Site,k)
#   }
#   NEE_allPPT_variable_kmean_RTPV(Site)
#   NEE_allPPT_NDVI_variable_kmean_RTPV(Site)
#   NEE_alllags_variable_kmean_RTPV(Site)
#   NEE_alllags_NDVI_variable_kmean_RTPV(Site)
# }

############################################################################
# All sites, one cluster amount
############################################################################
# Source the function
# source("alternate/RTPV/NEE_kmean_function_current_RTPV.R")
# source("alternate/RTPV/NEE_kmean_function_current_NDVI_RTPV.R")
# source("alternate/RTPV/NEE_kmean_function_singlePPT_RTPV.R")
# source("alternate/RTPV/NEE_kmean_function_singlePPT_NDVI_RTPV.R")
# source("alternate/RTPV/NEE_kmean_function_allPPT_RTPV.R")
# source("alternate/RTPV/NEE_kmean_function_allPPT_NDVI_RTPV.R")
# source("alternate/RTPV/NEE_kmean_function_alllags_RTPV.R")
# source("alternate/RTPV/NEE_kmean_function_alllags_NDVI_RTPV.R")
# 
# k = 5
# # Calculate the values for the variable functions
# for (Site in Sites){
#   NEE_current_kmean_RTPV(Site,k)
#   NEE_current_NDVI_kmean_RTPV(Site,k)
#   for (Lag in 2:8){
#     NEE_singlePPT_kmean_RTPV(Site,k,Lag)
#     NEE_singlePPT_NDVI_kmean_RTPV(Site,k,Lag)
#   }
#   NEE_allPPT_kmean_RTPV(Site,k)
#   NEE_allPPT_NDVI_kmean_RTPV(Site,k)
#   NEE_alllags_kmean_RTPV(Site,k)
#   NEE_alllags_NDVI_kmean_RTPV(Site,k)
# }

############################################################################
# For a single Site but many different clusters!!!
############################################################################

# # Source the function
# source("alternate/RTPV/NEE_kmean_function_current_RTPV.R")
# source("alternate/RTPV/NEE_kmean_function_current_NDVI_RTPV.R")
# source("alternate/RTPV/NEE_kmean_function_singlePPT_RTPV.R")
# source("alternate/RTPV/NEE_kmean_function_singlePPT_NDVI_RTPV.R")
# source("alternate/RTPV/NEE_kmean_function_allPPT_RTPV.R")
# source("alternate/RTPV/NEE_kmean_function_allPPT_NDVI_RTPV.R")
# source("alternate/RTPV/NEE_kmean_function_alllags_RTPV.R")
# source("alternate/RTPV/NEE_kmean_function_alllags_NDVI_RTPV.R")
# 
# # Set the site
# Sites = "AU-How"
# 
# ks = c(60,70,80,90,100)
# # Calculate the values for the variable functions
# for (k in ks){
#   for (Site in Sites){
# #    NEE_current_kmean_RTPV(Site,k)
#     NEE_current_NDVI_kmean_RTPV(Site,k)
#     for (Lag in 2:8){
# #      NEE_singlePPT_kmean_RTPV(Site,k,Lag)
#       NEE_singlePPT_NDVI_kmean_RTPV(Site,k,Lag)
#     }
# #    NEE_allPPT_kmean_RTPV(Site,k)
#     NEE_allPPT_NDVI_kmean_RTPV(Site,k)
# #    NEE_alllags_kmean_RTPV(Site,k)
#     NEE_alllags_NDVI_kmean_RTPV(Site,k)
#   }
# }


# k = 4
# Type = "allPPT"
# source("alternate/RTPV/NEE_kmean_function_allPPT_NDVI_RTPV.R")
# source("functions/NEE_kmeans_CoeffPlot_function.R")
# for (Site in Sites){
#   NEE_allPPT_NDVI_kmean_RTPV(Site,k)
#   Plot = ClusterCoefficients(Site,k,Type)
#   plot(Plot)
#   Sys.sleep(5)
# }
# 

source("alternate/RTPV/NEE_kmean_function_current_NDVI_RTPV.R")
source("alternate/RTPV/NEE_kmean_function_Fsdlags_NDVI_RTPV.R")
source("alternate/RTPV/NEE_kmean_function_Tairlags_NDVI_RTPV.R")
source("alternate/RTPV/NEE_kmean_function_VPDlags_NDVI_RTPV.R")
source("alternate/RTPV/NEE_kmean_function_PPTlags_NDVI_RTPV.R")
source("alternate/RTPV/NEE_kmean_function_allPPT_NDVI_RTPV.R")
source("alternate/RTPV/NEE_kmean_function_alllags_NDVI_RTPV.R")
source("alternate/RTPV/NEE_kmean_function_alllags_scaledNEE_NDVI_RTPV.R")
for (k in 2:8){
  for (Site in Sites){
    #NEE_alllags_scaledNEE_NDVI_kmean_RTPV(Site,k)
    #NEE_alllags_NDVI_kmean_RTPV(Site,k)
    NEE_allPPT_NDVI_kmean_RTPV(Site,k)
    NEE_current_NDVI_kmean_RTPV(Site,k)
    #NEE_Fsdlags_NDVI_kmean_RTPV(Site,k)
    #NEE_Tairlags_NDVI_kmean_RTPV(Site,k)
    #NEE_VPDlags_NDVI_kmean_RTPV(Site,k)
    #NEE_PPTlags_NDVI_kmean_RTPV(Site,k)
  }
}


source("alternate/RTPV/NEE_kmean_function_Fsdlags_NDVI_RTPV.R")
source("alternate/RTPV/NEE_kmean_function_Tairlags_NDVI_RTPV.R")
source("alternate/RTPV/NEE_kmean_function_VPDlags_NDVI_RTPV.R")
source("alternate/RTPV/NEE_kmean_function_PPTlags_NDVI_RTPV.R")
source("alternate/RTPV/NEE_kmean_function_allPPT_NDVI_RTPV.R")
for (k in 167:200){
    NEE_PPTlags_NDVI_kmean_RTPV("AU-How",k)
}


NEE_alllags_NDVI_kmean_RTPV = function(Site,k){
  
# Load the required packages
library(cluster)
library(tidyverse)
library(factoextra)
library(NbClust)

starttime = Sys.time()  

message("Performing ",k,"-cluster k-means clustering with NDVI and all lags for ",Site," at ",Sys.time())
# Load the input file and extract required data
load(paste0("inputs/RTPV/",Site,"_Input_RTPV.Rdata"))
input = eval(as.name(paste0(Site,"_Input")))
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
                input$NDVI,
                input$ppt_multiscale)
colnames(climate) = c("Ta",
                      "Fsd",
                      "VPD",
                      "PPT",
                      "Ta_1",
                      "Fsd_1",
                      "VPD_1",
                      "PPT_1",
                      "Ta_2",
                      "Fsd_2",
                      "VPD_2",
                      "PPT_2",
                      "Ta_3",
                      "Fsd_3",
                      "VPD_3",
                      "PPT_3",
                      "Ta_4",
                      "Fsd_4",
                      "VPD_4",
                      "PPT_4",
                      "Ta_5",
                      "Fsd_5",
                      "VPD_5",
                      "PPT_5",
                      "Ta_6",
                      "Fsd_6",
                      "VPD_6",
                      "PPT_6",
                      "Ta_7",
                      "Fsd_7",
                      "VPD_7",
                      "PPT_7",
                      "Ta_8",
                      "Fsd_8",
                      "VPD_8",
                      "PPT_8",
                      "Ta_9",
                      "Fsd_9",
                      "VPD_9",
                      "PPT_9",
                      "Ta_10",
                      "Fsd_10",
                      "VPD_10",
                      "PPT_10",
                      "Ta_11",
                      "Fsd_11",
                      "VPD_11",
                      "PPT_11",
                      "Ta_12",
                      "Fsd_12",
                      "VPD_12",
                      "PPT_12",
                      "Ta_13",
                      "Fsd_13",
                      "VPD_13",
                      "PPT_13",
                      "NDVI",
                      "PPT0",
                      "PPT_14_20",
                      "PPT_21_29",
                      "PPT_30_59",
                      "PPT_60_119",
                      "PPT_120_179",
                      "PPT_180_269",
                      "PPT_270_365")
# Remove PPT intercept
climate = climate[,-c(58)]
# Remove first year, which has no PPT data and scale
climate = scale(climate[-(1:365),])
NEE = input$NEE[-(1:365)]

# Find the cluster allocations for recommended number of clusters
kmean.output = kmeans(climate,k,iter.max = 100, nstart = 50)

# Initialise the comparison dataframe
compare = data.frame("NEE_obs" = NEE,"NEE_pred" = 0, "cluster" = kmean.output$cluster)

output = list()

for (i in 1:k){
  climate_cluster = climate[kmean.output$cluster==i,]
  NEE_cluster = NEE[kmean.output$cluster==i]
  lin.mod = lm(NEE_cluster ~ climate_cluster,na.action = na.exclude)
  r.squared = summary(lin.mod)$r.squared
  # Place the k-means fitted NEE into the data frame
  compare$NEE_pred[kmean.output$cluster==i] = fitted(lin.mod)
  # Assign and output the cluster info
  name = paste0("cluster_",i)
  assign(name,list("climate" = climate_cluster,
                   "NEE" = NEE_cluster,
                   "model" = lin.mod,
                   "r.squared" = r.squared,
                   "mod.r" = r.squared*nrow(climate_cluster)))
  
  output[[paste0("cluster_",i)]] = eval(as.name(name))
  # Tidy up
  rm(list = c("climate_cluster",
              "NEE_cluster",
              "name",
              "lin.mod",
              "r.squared"))
}

NEE_obs = compare$NEE_obs
NEE_pred = compare$NEE_pred

if (any(kmean.output$size<50)){
  message("                     ##**## WARNING! ##**##\n",
          sum(kmean.output$size<50)," clusters have too few observations for a reliable regression!\n",
          "                     ##**## WARNING! ##**##")
  Sys.sleep(1)
}

output[["r.squared"]] = summary(lm(compare$NEE_obs ~ compare$NEE_pred))$r.squared
output[["MBE"]] = sum(NEE_pred-NEE_obs,na.rm=TRUE)/length(NEE_pred)
output[["NME"]] = sum(abs(NEE_pred-NEE_obs),na.rm=TRUE)/sum(abs(mean(NEE_obs,na.rm=TRUE)-NEE_obs),na.rm=TRUE)
output[["SDD"]] = abs(1-sd(NEE_pred,na.rm=TRUE)/sd(NEE_obs,na.rm=TRUE))
output[["CCO"]] = cor(NEE_pred,NEE_obs,use = "complete.obs", method = "pearson")
output[["series"]] = compare
output[["totwithinss"]] = kmean.output$tot.withinss
output[["bet.tot.ratio"]] = kmean.output$betweenss/kmean.output$totss

ss <- silhouette(kmean.output$cluster, dist(climate))
ss = mean(ss[, 3])
output[["avg.sil"]] = ss

runtime = Sys.time()-starttime
output[["runtime"]] = runtime
save(output,file = paste0("alternate/RTPV/results/NEE_output_",k,"cluster_kmean_alllags_NDVI_RTPV_",Site,".Rdata"))

}

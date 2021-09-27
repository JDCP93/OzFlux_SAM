# LE SAM memory exploration
# Based on Liu et al, 2019
# Uses the same model

# Start the workflow
rm(list = ls())
message("Start the workflow at ",Sys.time())

# load needed packages
nee_packs <- c('rjags', 'coda', 'stats', 'R2jags', 'parallel','runjags','mcmcplots')
lapply(nee_packs, require, character.only = T)

# variables to monitor
monitor_vars <- c("b0","b1","sig.res","LE.res","LE.res_rep","mu.res","deviance")

# Load the site inputs
load('inputs/RTPV/LE/AU-ASM_input_LE_RTPV.Rdata') 
data = `AU-ASM_input`

# Load in the output data 
# Look in folder "results" for the data
File = list.files("output/RTPV/",pattern = "LE_output_RTPV_AU-ASM")
# Read the data into R - note that if multiple results are available for a 
# site, we take the most recent
message("File is ",File[length(File)])
load(paste0("output/RTPV/",File[length(File)]))

# Either take the object already saved as an mcmc object for the current 
# workflows or, to maintain compatibility with older workflows, calculate it
# from the rjags object
if (class(output) == "list"){
  output.mcmc = output$output.mcmc
}else{
  output.mcmc = as.mcmc.rjags(output)
}
rm(output)

# Make sure the file is in mcmc format
output.mcmc = as.mcmc.list(output.mcmc)

message("Extracting observed LE at ",Sys.time())
# Extract predicted and observed LE
summary = summary(output.mcmc)
LE_pred = summary$statistics[substr(rownames(summary$statistics),1,7)=="LE_pred",1]
LE_obs = data$LE[-(1:365)]
LE_res = LE_pred - LE_obs

inputdata = list("LE.res" = LE_res,
                 "Nmem"=data$Nmem)


# parallelize using runjags
message("Begin model run at ",Sys.time())
# run model in parallel with 6 chains and cores
output <- jags.parallel(model.file = 'LEModel_AR1_RTPV_r2jags.R',
                            parameters.to.save = monitor_vars,
                            data = inputdata,
                            n.chains = 6, 
                            n.burnin = 500000, 
                            n.iter = 1000000,
                            jags.module = c('glm','dic'),
                            n.thin = 50)

message("Save model output at ",Sys.time())
# Transform output into mcmc object to save space
output.mcmc = as.mcmc.rjags(output)
DIC = output$BUGSoutput$DIC
pD = output$BUGSoutput$pD
output = list("output.mcmc"=output.mcmc,
              "DIC" = DIC,
              "pD" = pD)
# Save the results
save(output, file=paste('LE_AR1_output_RTPV_AU-ASM_', Sys.Date(),'.Rdata', sep = ''))

# Tidy up
rm(list=ls())

# Finish
message("Finished the model run at ",Sys.time())
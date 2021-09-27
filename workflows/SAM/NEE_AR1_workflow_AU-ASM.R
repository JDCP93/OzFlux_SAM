# NEE SAM memory exploration
# Based on Liu et al, 2019
# Uses the same model

# Start the workflow
rm(list = ls())
message("Start the workflow at ",Sys.time())

# load needed packages
nee_packs <- c('rjags', 'coda', 'stats', 'R2jags', 'parallel','runjags','mcmcplots')
lapply(nee_packs, require, character.only = T)

# variables to monitor
monitor_vars <- c("b0","b1","sig.res","NEE.res","NEE.res_rep","mu.res","deviance")

# Load the site inputs
load('inputs/RTPV/AU-ASM_Input_RTPV.Rdata') 
data = `AU-ASM_Input`

# Load in the output data 
# Look in folder "results" for the data
File = list.files("output/RTPV/",pattern = "NEE_output_RTPV_AU-ASM")
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

output.mcmc = as.mcmc.list(output.mcmc)

message("Extracting observed NEE at ",Sys.time())
# Extract predicted and observed NEE
summary = summary(output.mcmc)
NEE_pred = summary$statistics[substr(rownames(summary$statistics),1,8)=="NEE_pred",1]
NEE_obs = data$NEE[-(1:365)]
NEE_res = NEE_pred - NEE_obs

inputdata = list("NEE.res" = NEE_res,
                 "Nmem"=data$Nmem)


# parallelize using runjags
message("Begin model run at ",Sys.time())
# run model in parallel with 6 chains and cores
output <- jags.parallel(model.file = 'NEEModel_AR1_RTPV_r2jags.R',
                            parameters.to.save = monitor_vars,
                            data = inputdata,
                            n.chains = 6, 
                            n.burnin = 100000, 
                            n.iter = 500000,
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
save(output, file=paste('NEE_AR1_output_RTPV_AU-ASM_', Sys.Date(),'.Rdata', sep = ''))

# Tidy up
rm(list=ls())

# Finish
message("Finished the model run at ",Sys.time())
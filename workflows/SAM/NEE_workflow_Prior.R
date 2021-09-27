# NEE SAM memory exploration
# Based on Liu et al, 2019
# Uses the same model with SWC removed and PPT instead. NDVI is also 0-1 normalised

# Start the workflow
rm(list = ls())
message("Start the workflow at ",Sys.time())

# load needed packages
nee_packs <- c('rjags', 'coda', 'stats', 'R2jags', 'parallel','runjags','mcmcplots')
lapply(nee_packs, require, character.only = T)

# variables to monitor
monitor_vars <- c("an", "ag", "phi0", "deltaXA", "weightA", "weightAP", "deltaXAP", 
                  "cum_weightA", "cum_weightAP", "sig_y", "NEE_pred","muNEE",
                  "ESen")

# Load the site inputs
load('./inputs/RTPV/AU-ASM_Input_RTPV.Rdata') 
data = `AU-ASM_Input`
inputdata = list("Nv"=data$Nv,
                 "Ns"=data$Ns,
                 "Nlag"=data$Nlag,
                 "Nmem"=data$Nmem,
                 "NlagP"=data$NlagP,
                 "Mem_records"=data$Mem_records,
                 "clim"=data$clim,
                 "ppt_multiscale"=data$ppt_multiscale,
                 "NEE"=NULL,
                 "NDVI"=data$NDVI,
                 "NblocksP"=data$NblocksP,
                 "block"=data$block,
                 "BlockSize"=data$BlockSize,
                 "Nblocks"=data$Nblocks)


# parallelize using runjags
message("Begin model run at ",Sys.time())

# For some reason, jags.parallel errors "subscript out of bounds" when trying
# to calculate a prior. As such, we use jags.model plus coda.samples.
jags = jags.model("NEEModel_RTPV_r2jags.R", data=inputdata, n.chains=6, n.adapt=100000) 
output = coda.samples(jags, n.iter=1000000, n.burnin=100000, thin=1000,
                      variable.names=monitor_vars)


message("Save model output at ",Sys.time())
# Transform output into mcmc object to save space
output.mcmc = output
output = list("output.mcmc"=output.mcmc)
# Save the results
save(output, file=paste('NEE_output_RTPV_Prior_', Sys.Date(),'.Rdata', sep = ''))

# Tidy up
rm(list=ls())

# Finish
message("Finished the model run at ",Sys.time())

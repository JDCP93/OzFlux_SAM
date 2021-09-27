# LE SAM memory exploration
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
                  "cum_weightA", "cum_weightAP", "sig_y", "LE_pred","muLE",
                  "ESen", 'deviance')

# Load the site inputs
load('./inputs/RTPV/LE/AU-ASM_input_LE_RTPV.Rdata') 
data = `AU-ASM_input`
inputdata = list("Nv"=data$Nv,
                 "Ns"=data$Ns,
                 "Nlag"=data$Nlag,
                 "Nmem"=data$Nmem,
                 "NlagP"=data$NlagP,
                 "Mem_records"=data$Mem_records,
                 "clim"=data$clim,
                 "ppt_multiscale"=data$ppt_multiscale,
                 "LE"=data$LE,
                 "NDVI"=data$NDVI,
                 "NblocksP"=data$NblocksP,
                 "block"=data$block,
                 "BlockSize"=data$BlockSize,
                 "Nblocks"=data$Nblocks)


# parallelize using runjags
message("Begin model run at ",Sys.time())
# run model in parallel with 6 chains and cores
output <- jags.parallel(model.file = 'LEModel_current_RTPV_r2jags.R',
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
save(output, file=paste('LE_current_output_RTPV_AU-ASM_', Sys.Date(),'.Rdata', sep = ''))

# Tidy up
rm(list=ls())

# Finish
message("Finished the model run at ",Sys.time())

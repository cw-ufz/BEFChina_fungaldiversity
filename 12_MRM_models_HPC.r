######### Parallel R script to calculate all mrm model combinations on EVE ####

# R script to be run on cluster computer needed for 
# R script (12/13) for analyses in Weissbecker et al. 2018
# version: April 2018

#######################

library(parallel)
library("ecodist", lib.loc = "/gpfs0/home/lachmann/R/x86_64-pc-linux-gnu-library/3.4/")
##########

cl <- makeForkCluster(outfile = "")

# Initialize data
input<-"/gpfs1/work/lachmann/Rprojects/China_beta_NP/Data"
output<-"/gpfs1/work/lachmann/Rprojects/China_beta_NP/output"

RDS.object<-readRDS(file= paste0(input,"/10_AM_beta_mrm_var_avreps.RDS"))
all_model_combinations<-readRDS(file= paste0(input,"/11_AM_model_combinations_avreps.RDS"))

# Make variables available to cluster 
clusterExport(cl=cl, varlist=c("RDS.object", "all_model_combinations", "input", "output"))

results<-list()
for(i in 1:(length(all_model_combinations)-1)){
  results[[i]] <- parApply(cl, X = all_model_combinations[[i]], MARGIN = 2, function(x) MRM(as.formula(paste("RDS.object$Fungal_beta ~ ", paste(paste0("RDS.object$",x),collapse="+"))), nperm=10000,mrank=FALSE))
}

saveRDS(results, file = paste0(output,"/12_AM_EVE_mrm_bestsubsets_avreps.RDS"))

stopCluster(cl)

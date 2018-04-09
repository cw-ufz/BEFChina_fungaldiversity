#### Best Subset Model selection: Multiple regression of distance matrices ####

# R script (11_12/13) for analyses in Weissbecker et al. 2018
# version: April 2018

# actual MRM calculation done on a more powerful computer (cluster computer) 
################################

library(dplyr)      #dplyr_0.7.4
library(ecodist)    #ecodist_2.0.1 
library(vegan)      #vegan_2.4-4  
#R version 3.4.2 (2017-09-28)
###############################

source(file="./R_functions/11_MRM_subset_combinations_NP.R")

# Preparation of complete subset model lists 
Sapro_beta_mrm_var<-readRDS(file="./output/10_Sapro_beta_mrm_var_avreps.RDS")
Patho_beta_mrm_var<-readRDS(file="./output/10_Patho_beta_mrm_var_avreps.RDS")
EcM_beta_mrm_var<-readRDS(file="./output/10_EcM_beta_mrm_var_avreps.RDS")
AM_beta_mrm_var<-readRDS(file="./output/10_AM_beta_mrm_var_avreps.RDS")

Sapro_mrm_names<-as.list(names(Sapro_beta_mrm_var)[-1])
Sapro_model_combinations<-mrm.subset.list(Sapro_mrm_names)

EcM_mrm_names<-as.list(names(EcM_beta_mrm_var)[-1])
EcM_model_combinations<-mrm.subset.list(EcM_mrm_names)

Patho_mrm_names<-as.list(names(Patho_beta_mrm_var)[-1])
Patho_model_combinations<-mrm.subset.list(Patho_mrm_names)

AM_mrm_names<-as.list(names(AM_beta_mrm_var)[-1])
AM_model_combinations<-mrm.subset.list(AM_mrm_names)

############## script output 

# MRM model combinations list
#saveRDS(Sapro_model_combinations, file="./output/11_Sapro_model_combinations_avreps.RDS")
#saveRDS(EcM_model_combinations, file="./output/11_EcM_model_combinations_avreps.RDS")
#saveRDS(Patho_model_combinations, file="./output/11_Patho_model_combinations_avreps.RDS")
#saveRDS(AM_model_combinations, file="./output/11_AM_model_combinations_avreps.RDS")

############# ........ Run analysis on EVE with extra script ......... #################


############# ........ Load results obtained by cluster computer       #################  

results_dopar<-readRDS(file="./output/12_AM_EVE_mrm_bestsubsets_avreps.RDS")
RDS.object<-readRDS("./output/10_AM_beta_mrm_var_avreps.RDS")
model_combinations<-readRDS("./output/11_AM_model_combinations_avreps.RDS")

dim_EVE<-length(results_dopar[])

dopar_mrm_pre<-list()
dopar_mrm_pre2<-list()
list_results_converted<-list()
for(i in 1:dim_EVE){
  dopar_mrm_pre[[i]]<-lapply(results_dopar[[i]], function(x) c(length(rownames(x$coef))-1, 
                                                               x$F.test, x$r, x$coef[,2]))
  dopar_mrm_pre2[[i]]<-lapply(dopar_mrm_pre[[i]], function(x) t(cbind(x)))
  list_results_converted[[i]]<-ldply(dopar_mrm_pre2[[i]], data.frame)
  }
unlist_results_dopar <- do.call(rbind, list_results_converted)


# Calculate full model locally
full_model_comb<-length(model_combinations[])
full_model<-MRM(as.formula(paste("RDS.object$Fungal_beta ~ ", 
              paste(paste0("RDS.object$",model_combinations[[full_model_comb]]),
                    collapse="+"))), nperm=10000,mrank=FALSE)

# Merge results
full_mrm_pre<-c(length(rownames(full_model$coef))-1, full_model$F.test,
                full_model$r, full_model$coef[,2])
all_unlist<-rbind(unlist_results_dopar,full_mrm_pre)
colnames(all_unlist)[] <- lapply(colnames(all_unlist), gsub, 
                                 pattern = "RDS.object.", replacement = "",
                                 fixed = TRUE)

# Extract the best subsets from the result file
results_ordered<-all_unlist[order(all_unlist$V1, -as.numeric(all_unlist$F)),]
no_sets<-3
best_sets<-results_ordered %>% group_by(V1) %>% top_n(no_sets, F)

### script output

#---type the right file manually:
#write.table(results_ordered, file="./output/12_AM_EVE_mrm_readable_avreps.csv")
#write.table(best_sets, file="./output/12_AM_EVE_mrm_readable_bestsets_avreps.csv")





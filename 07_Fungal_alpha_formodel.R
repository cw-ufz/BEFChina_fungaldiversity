############################### Fungal alpha diversity analysis ###############

# R script (07/13) for analyses in Weissbecker et al. 2018 New Phytologist
# version: April 2018

###############################
library(vegan)      #vegan_2.4-4
library(sp)         #sp_1.2-5            #spDists function
library(car)        #car_2.1-5           #Anova
library(lme4)       #lme4_1.1-14         #lmer
library(nlme)       #nlme_3.1-131        #lme
library(MuMIn)      #MuMIn_1.40.4        #corr.AIC, provides R^2 for lmer models
library(data.table) #data.table_1.10.4-1 #merge list to dataframe
library(dplyr)      #dplyr_0.7.4         
#R version 3.4.2
################

### Warning message in forward slection procedure: 
#"In selected_cols != which(names(model_results_ordered) %in% retained_vars) :
#longer object length is not a multiple of shorter object length"
#-->This is due to the R recycling habit and can be ignored here!

###############

source(file="./R_functions/07_Vif_function.r")

Sapro_otu_sum_avreps<-readRDS(file="./output/03_Sapro_otu_sum_avreps.RDS")
EcM_otu_sum_avreps<-readRDS(file="./output/03_EcM_otu_sum_avreps.RDS")
Patho_otu_sum_avreps<-readRDS(file="./output/03_Patho_otu_sum_avreps.RDS")
AM_otu_sum_avreps<-readRDS(file="./output/03_AM_otu_sum_avreps.RDS")

Sapro_meta_av<-readRDS(file="./output/03_Sapro_metadata_trans_avreps.RDS")
EcM_meta_av<-readRDS(file="./output/03_EcM_metadata_trans_avreps.RDS")
Patho_meta_av<-readRDS(file="./output/03_Patho_metadata_trans_avreps.RDS")
AM_meta_av<-readRDS(file="./output/03_AM_metadata_trans_avreps.RDS")

# Fungal richness calculation

Fungal_otudata<-AM_otu_sum_avreps
Fungal_metadata<-AM_meta_av

fungal_richness<-specnumber(Fungal_otudata)
fungal_rich_residuals<-residuals(lm(fungal_richness~
                                      sqrt(rowSums(Fungal_otudata))))

# Generate Prinicpal Components Neighborhod Matrix (PCNM) eigenvectors

geo_distance_m <- spDists(as.matrix(Fungal_metadata$location), 
                          as.matrix(Fungal_metadata$location),  
                          longlat=TRUE)*1000
pcnm_alpha<-pcnm(geo_distance_m) # auto threshold, sapro=98.06807, weights=1
pcnm_alpha<-pcnm_alpha$vectors
pcnm_alpha<-data.frame(pcnm_alpha)

# Transform to factor type: Sample_Tree, Plot_ID

Sample_tree_factor<-data.frame()
for (i in seq(1,nrow(Fungal_metadata$Sample_tree),by=1)){
  Sample_tree_factor[i,1]<-names(Fungal_metadata$Sample_tree)[which(Fungal_metadata$Sample_tree[i,]==1)]
}
names(Sample_tree_factor)<-"Sample_tree_factor"
Sample_tree_factor<-as.factor(Sample_tree_factor$Sample_tree_factor)

Plot_factor<-data.frame()
for (i in seq(1,nrow(Fungal_metadata$plot),by=1)){
 Plot_factor[i,1]<-names(Fungal_metadata$plot)[which(Fungal_metadata$plot[i,]==1)]
}
names(Plot_factor)<-"Plot_ID"
Plot_factor<-as.factor(Plot_factor$Plot_ID)


Fungal_data_alpha<-list(Fungal_metadata$Env_var, Sample_tree_factor,Plot_factor,
                        Fungal_metadata$richness_vars, pcnm_alpha)

names(Fungal_data_alpha)<-c("Env_var", "Sample_tree_factor", "Plot_factor",
                            "richness_vars", "pcnm_alpha")

### vif (variance inflation) analysis: stepwise vif function in 07_vif_function

fungal_rich_residuals<-readRDS(file="./output/07_AM_richness_avreps_residuals.RDS")
Fungal_data_alpha<-readRDS(file="./output/07_AM_alpha_avreps_metadata.RDS")

vif_test<-vif_func(list(Fungal_data_alpha$pcnm_alpha, Fungal_data_alpha$Env_var,
                        Fungal_data_alpha$richness_vars),thresh=10,trace=T)
vif_test[length(vif_test)+1]<-"Plot_factor"
vif_test[length(vif_test)+1]<-"Sample_tree_factor"

# Full model formulation 
fungal_rich_residuals<-readRDS(file="./output/07_EcM_richness_avreps_residuals.RDS")
vif_test<-readRDS(file="./output/07_EcM_alpha_vifallmodelvars_avreps.RDS")
Fungal_data_alpha<-readRDS(file="./output/07_EcM_alpha_avreps_metadata.RDS")

Sample_tree_factor<-data.frame(Fungal_data_alpha$Sample_tree_factor)
names(Sample_tree_factor)<-"Sample_tree_factor"

Plot_factor<-data.frame(Fungal_data_alpha$Plot_factor)
names(Plot_factor)<-"Plot_factor"

### mixed effect model selection                               

# Set variables for first run
var_set<-vif_test[-c(which(vif_test %in% c("Plot_factor")))]
all_vars<-var_set
retained_vars<-vector()

model_select<-list()
model_results_comb<-list()
for (i in seq(from=1, to=length(all_vars), by=1)){
  model_select[[i]]<-list()
}

for(i in seq(1,length(all_vars),by=1)){
  print("Anzahl Variablen:");print(no_vars<-length(retained_vars)+1)  
  print("potenzielle Variablen:"); print(var_set)
  print("Selected vars:"); print(retained_vars)
  for (iteration in seq(from=1, to=length(var_set), by=1)){
    no_vars<-length(retained_vars)+1
    
    if(length(retained_vars)==0){
      lmer_model<-lmer(as.formula(paste("fungal_rich_residuals ~ ", 
            paste(paste0(c(var_set[iteration],"(1|Plot_factor)"),collapse="+")))), 
            data=c(Fungal_data_alpha$pcnm_alpha,Fungal_data_alpha$Env_var,
                   Fungal_data_alpha$richness_vars,Plot_factor, Sample_tree_factor),
            contrasts=list(Plot_factor="contr.sum", Sample_tree_factor="contr.sum"), 
            REML=FALSE)
      if(var_set[iteration]=="Sample_tree_factor"){
        lme_model<-lme(as.formula(paste("fungal_rich_residuals ~ ",
            paste(paste0(var_set[iteration],collapse="+")))),random=~1|Plot_factor,
            data=c(Fungal_data_alpha$pcnm_alpha,Fungal_data_alpha$Env_var,
                  Fungal_data_alpha$richness_vars,Plot_factor, Sample_tree_factor),
            contrasts=list(Plot_factor="contr.sum", Sample_tree_factor="contr.sum"),
            method="REML")
      }else{
        lme_model<-lme(as.formula(paste("fungal_rich_residuals ~ ",
            paste(paste0(var_set[iteration],collapse="+")))),random=~1|Plot_factor, 
            data=c(Fungal_data_alpha$pcnm_alpha,Fungal_data_alpha$Env_var,
                   Fungal_data_alpha$richness_vars,Plot_factor, Sample_tree_factor),
            contrasts=list(Plot_factor="contr.sum"), method="REML")
      }
    }else if(length(retained_vars)==length(all_vars)){
      lmer_model<-lmer(as.formula(paste("fungal_rich_residuals ~ ",
            paste(paste0(c(retained_vars,"(1|Plot_factor)"),collapse="+")))), 
            data=c(Fungal_data_alpha$pcnm_alpha,Fungal_data_alpha$Env_var,
                   Fungal_data_alpha$richness_vars,Plot_factor, Sample_tree_factor), 
            contrasts=list(Plot_factor="contr.sum", Sample_tree_factor="contr.sum"),
            REML=FALSE)
      if("Sample_tree_factor" %in% retained_vars){
        lme_model<-lme(as.formula(paste("fungal_rich_residuals ~ ",
            paste(paste0(retained_vars,collapse="+")))),random=~1|Plot_factor, 
            data=c(Fungal_data_alpha$pcnm_alpha,Fungal_data_alpha$Env_var,
                   Fungal_data_alpha$richness_vars,Plot_factor, Sample_tree_factor),
            contrasts=list(Plot_factor="contr.sum", Sample_tree_factor="contr.sum"), 
            method="REML")
      }else{
        lme_model<-lme(as.formula(paste("fungal_rich_residuals ~ ", 
            paste(paste0(retained_vars,collapse="+")))),random=~1|Plot_factor, 
            data=c(Fungal_data_alpha$pcnm_alpha,Fungal_data_alpha$Env_var,
                   Fungal_data_alpha$richness_vars,Plot_factor, Sample_tree_factor),
            contrasts=list(Plot_factor="contr.sum"), method="REML")
      }
    }else{
      lmer_model<-lmer(as.formula(paste("fungal_rich_residuals ~ ", 
            paste(paste0(c(retained_vars,var_set[iteration],"(1|Plot_factor)"),
                         collapse="+")))), 
            data=c(Fungal_data_alpha$pcnm_alpha,Fungal_data_alpha$Env_var,
                   Fungal_data_alpha$richness_vars,Plot_factor, Sample_tree_factor),
            contrasts=list(Plot_factor="contr.sum", Sample_tree_factor="contr.sum"),
            REML=FALSE)
      if("Sample_tree_factor" %in% c(retained_vars, var_set[iteration])){
        lme_model<-lme(as.formula(paste("fungal_rich_residuals ~ ", 
            paste(paste0(c(retained_vars,var_set[iteration]),collapse="+")))),
            random=~1|Plot_factor, data=c(Fungal_data_alpha$pcnm_alpha,
            Fungal_data_alpha$Env_var,Fungal_data_alpha$richness_vars,
            Plot_factor, Sample_tree_factor),
            contrasts=list(Plot_factor="contr.sum", Sample_tree_factor="contr.sum"), 
            method="REML")
      }else{
        lme_model<-lme(as.formula(paste("fungal_rich_residuals ~ ",
            paste(paste0(c(retained_vars,var_set[iteration]),collapse="+")))),
            random=~1|Plot_factor, 
            data=c(Fungal_data_alpha$pcnm_alpha,Fungal_data_alpha$Env_var,
                   Fungal_data_alpha$richness_vars,Plot_factor, Sample_tree_factor),
            contrasts=list(Plot_factor="contr.sum"), method="REML")
      }
    }
    
    corrAIC<-AICc(lmer_model)
    r2_info<-r.squaredGLMM(lmer_model)
    
    eval_lme<-Anova(lme_model, white.adjust=TRUE, type=2)
    shapiro_norm<-shapiro.test(residuals(lme_model))
    
    model_select[[no_vars]][[iteration]]<-list(lmer_model,lme_model,corrAIC,
                                               r2_info,eval_lme,shapiro_norm)
    names(model_select[[no_vars]][[iteration]])<-c("lmer_model", "lme_model",
                                                   "AIC", "R2", "Anova", "Shapiro")
   
  }
  print("Modelle gebaut")
  
  # Extract model information to table
  model_results <- matrix(data=NA, nrow = length(var_set) , 
                          ncol = length(all_vars)+5)
  colnames(model_results)<-c("Vars", "AICc", "R2m", "R2c", "norm_p",
                             all_vars)
  
  for(no_model in seq(from=1, to=length(var_set), by=1)){
    model_results[no_model,1]<-length(retained_vars)+1
    model_results[no_model,2]<-model_select[[no_vars]][[no_model]]$AIC    
    model_results[no_model,3]<-model_select[[no_vars]][[no_model]]$R2[1]
    model_results[no_model,4]<-model_select[[no_vars]][[no_model]]$R2[2]
    model_results[no_model,5]<-model_select[[no_vars]][[no_model]]$Shapiro[2]$p.value
    
    result_col<-which(colnames(model_results) %in% 
                        rownames(model_select[[no_vars]][[no_model]]$Anova))
    model_results[no_model,result_col]<-model_select[[no_vars]][[no_model]]$Anova[3]$`Pr(>Chisq)`
  }
  print("Modellkennzahlen extrahiert")
  
  # Save ordered model_results in list element and rearrange variable sets
  model_results<-data.frame(model_results)
  model_results_ordered<-model_results[order(model_results$AICc),]
  model_results_comb[[no_vars]]<-model_results_ordered
  
  current_models<-which(model_results_ordered[,"Vars"]==no_vars)
  selected_cols<-which(!is.na(model_results_ordered[min(current_models),]))
  
  if(length(retained_vars)<1){
    selected_vars<-names(model_results_ordered)[selected_cols[selected_cols>5]]
    retained_vars<-selected_vars
    var_set<-var_set[var_set!=retained_vars]
    print("Modellkennzahlen gespeichert und neue var_set und retained_vars bestimmt")
  }else if (length(retained_vars)== length(all_vars)){
    print("Full model reached")
  }else {
    selected_vars<-names(model_results_ordered)[selected_cols[selected_cols>5 & 
                        selected_cols != which(names(model_results_ordered) %in% 
                                                 retained_vars)]]
    retained_vars<-unique(c(retained_vars,selected_vars))
    var_set<-setdiff(var_set, retained_vars)
    print("Modellkennzahlen gespeichert und neue var_set und retained_vars bestimmt")}
  
  if (!length(which(model_results_ordered[1,-(1:5)] > 0.05)) == 0) {
    stop("p-values in best model above 0.05, no further variables are added: 
         stop of forward selection")}
}

# Inspect results
alpha_forwmodels_plain<-rbindlist(lapply(model_results_comb, as.data.frame.list),
                                  fill=TRUE)
no_sets<-3
best_sets<-alpha_forwmodels_plain %>% group_by(Vars) %>% top_n(no_sets, -AICc)

# Inspect final model
final_vars<-retained_vars[1:(length(retained_vars)-1)]

lmer.full<- lmer(as.formula(paste("fungal_rich_residuals ~ ",
                paste(paste0(c(final_vars,"(1|Plot_factor)"),collapse="+")))), 
                data=c(Fungal_data_alpha$pcnm_alpha,Fungal_data_alpha$Env_var,
                       Fungal_data_alpha$richness_vars,Plot_factor, Sample_tree_factor),
                contrasts=list(Plot_factor="contr.sum", Sample_tree_factor="contr.sum"), 
                REML=FALSE)
if("Sample_tree_factor" %in% final_vars){
  lme.full<- lme(as.formula(paste("fungal_rich_residuals ~ ",
                paste(paste0(final_vars,collapse="+")))),random=~1|Plot_factor, 
                data=c(Fungal_data_alpha$pcnm_alpha,Fungal_data_alpha$Env_var,
                       Fungal_data_alpha$richness_vars,Sample_tree_factor, Plot_factor), 
                method="REML")
}else{
  lme.full<- lme(as.formula(paste("fungal_rich_residuals ~ ",
                paste(paste0(final_vars,collapse="+")))),random=~1|Plot_factor, 
                data=c(Fungal_data_alpha$pcnm_alpha,Fungal_data_alpha$Env_var,
                       Fungal_data_alpha$richness_vars,Plot_factor), 
                method="REML")
}

############## 

#---choose the right file manually:

# Fungal residual richness + full metadata 

#saveRDS(fungal_rich_residuals, file="./output/07_Sapro_richness_avreps_residuals.RDS")
#saveRDS(Fungal_data_alpha, file="./output/07_Sapro_alpha_avreps_metadata.RDS")

#saveRDS(fungal_rich_residuals, file="./output/07_Patho_richness_avreps_residuals.RDS")
#saveRDS(Fungal_data_alpha, file="./output/07_Patho_alpha_avreps_metadata.RDS")

#saveRDS(fungal_rich_residuals, file="./output/07_EcM_richness_avreps_residuals.RDS")
#saveRDS(Fungal_data_alpha, file="./output/07_EcM_alpha_avreps_metadata.RDS")

#saveRDS(fungal_rich_residuals, file="./output/07_AM_richness_avreps_residuals.RDS")
#saveRDS(Fungal_data_alpha, file="./output/07_AM_alpha_avreps_metadata.RDS")

# Retained vif model vars metadata

#saveRDS(vif_test, file="./output/07_Sapro_alpha_vifallmodelvars_avreps.RDS")
#saveRDS(vif_test, file="./output/07_EcM_alpha_vifallmodelvars_avreps.RDS")
#saveRDS(vif_test, file="./output/07_AM_alpha_vifallmodelvars_avreps.RDS")
#saveRDS(vif_test, file="./output/07_Patho_alpha_vifallmodelvars_avreps.RDS")

# Model selection results

# write.csv2(best_sets, file="./output/07_AM_alpha_forward_bestsets_avreps.csv")
# write.csv2(alpha_forwmodels_plain, file="./output/07_AM_alpha_forward_models_avreps.csv")
# saveRDS(final_vars, file="./output/07_AM_forward_model_retained_vars_avreps.RDS")
# write.csv2(Anova_bestmodel, file="./output/07_AM_alpha_forward_best_Anova_avreps.csv")



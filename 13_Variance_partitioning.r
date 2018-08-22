############################## Variance partitioning ##########################

# R script (13/13) for analyses in Weissbecker et al. 2018
# version: August 2018

# variation partitioning alpha diversity  
# lmer model --> rand() + VarCorr for random factor contribution, 
# hier.part for fixed factor contribution, Procedure descibed in:
# https://benwhalley.github.io/just-enough-r/icc-and-vpc.html
#################
library(vegan)      #vegan_2.5-2    #variance partitioning distance class data
library(lme4)       #lme4_1.1-17   #lmer, VarCorr
library(hier.part)  #hier.part_1.0-4#variance partitioning linear regressions
library(dplyr)      #dplyr_0.7.6     #can handle %>%, mutate
library(lmerTest)   #lmerTest_3.0-1#
library(MuMIn)      #MuMIn_1.42.1  #provides R^2 for lmer models
library(ggplot2)    #ggplot2_3.0.0
library(extrafont)  #extrafont_0.17
#font_import()
#fonts()
library(scales)     #scales_0.5.0
require(grid)     
#sessionInfo()
#R version 3.5.1 (2018-07-02)
#################

alpha_model_variances <- data.frame(frac=rep(c("location", "Plot","shared","soil",
                                               "unexplained"), each=4),
                                   group=rep(c("Saprotrophs", "Pathogens",
                                               "EcM fungi", "AM fungi"),5),
                                   perc=NA)
levels(alpha_model_variances$frac)<-c("unexplained","shared",
                                      "Plot", "location", "soil")

fungal_rich_residuals<-readRDS(file="./output/07_Sapro_richness_avreps_residuals.RDS")
Fungal_data_alpha<-readRDS(file="./output/07_Sapro_alpha_avreps_metadata.RDS")

fungal_rich_residuals<-readRDS(file="./output/07_Patho_richness_avreps_residuals.RDS")
Fungal_data_alpha<-readRDS(file="./output/07_Patho_alpha_avreps_metadata.RDS")

fungal_rich_residuals<-readRDS(file="./output/07_EcM_richness_avreps_residuals.RDS")
Fungal_data_alpha<-readRDS(file="./output/07_EcM_alpha_avreps_metadata.RDS")

fungal_rich_residuals<-readRDS(file="./output/07_AM_richness_avreps_residuals.RDS")
Fungal_data_alpha<-readRDS(file="./output/07_AM_alpha_avreps_metadata.RDS")

Plot_factor<-data.frame(Fungal_data_alpha$Plot_factor)
names(Plot_factor)<-"Plot_factor"

# Sapro
i=0
fungal_alpha_best<- lmer(fungal_rich_residuals~Fungal_data_alpha$Env_var$Ntot+
                    Fungal_data_alpha$pcnm_alpha$PCNM8+Fungal_data_alpha$pcnm_alpha$PCNM21+
                    Fungal_data_alpha$pcnm_alpha$PCNM3+(1|Plot_factor), 
                    data=c(Fungal_data_alpha$pcnm_alpha,Fungal_data_alpha$Env_var,
                    Fungal_data_alpha),REML=TRUE)

Fungal_best_envvar<-data.frame(Fungal_data_alpha$Env_var$Ntot,
                    Fungal_data_alpha$pcnm_alpha$PCNM8,
                    Fungal_data_alpha$pcnm_alpha$PCNM21,
                    Fungal_data_alpha$pcnm_alpha$PCNM3)


# Patho
i=1

fungal_alpha_best<- lmer(fungal_rich_residuals~Fungal_data_alpha$pcnm_alpha$PCNM36+
                    Fungal_data_alpha$pcnm_alpha$PCNM24+Fungal_data_alpha$pcnm_alpha$PCNM40+
                    (1|Plot_factor), 
                    data=c(Fungal_data_alpha$pcnm_alpha,Fungal_data_alpha), REML=TRUE)

Fungal_best_envvar<-data.frame(Fungal_data_alpha$pcnm_alpha$PCNM36,
                               Fungal_data_alpha$pcnm_alpha$PCNM24,
                               Fungal_data_alpha$pcnm_alpha$PCNM40)

# EcM
i=2

fungal_alpha_best<- lmer(fungal_rich_residuals~Fungal_data_alpha$Env_var$Soil_moisture+
                    Fungal_data_alpha$pcnm_alpha$PCNM5+Fungal_data_alpha$pcnm_alpha$PCNM12+
                    Fungal_data_alpha$pcnm_alpha$PCNM1+Fungal_data_alpha$pcnm_alpha$PCNM29+
                    Fungal_data_alpha$pcnm_alpha$PCNM14+(1|Plot_factor), 
                    data=c(Fungal_data_alpha$pcnm_alpha,Fungal_data_alpha$Env_var,Fungal_data_alpha),
                    REML=TRUE)

Fungal_best_envvar<-data.frame(Fungal_data_alpha$Env_var$Soil_moisture,
                               Fungal_data_alpha$pcnm_alpha$PCNM5,
                               Fungal_data_alpha$pcnm_alpha$PCNM12,
                               Fungal_data_alpha$pcnm_alpha$PCNM1,
                               Fungal_data_alpha$pcnm_alpha$PCNM29,
                               Fungal_data_alpha$pcnm_alpha$PCNM14)

# AM
i=3

fungal_alpha_best<- lmer(fungal_rich_residuals~Fungal_data_alpha$Env_var$Soil_moisture+
                    Fungal_data_alpha$Env_var$SLOPE+Fungal_data_alpha$Env_var$KAK.eff+
                    Fungal_data_alpha$pcnm_alpha$PCNM25+Fungal_data_alpha$pcnm_alpha$PCNM13+
                    Fungal_data_alpha$pcnm_alpha$PCNM34+Fungal_data_alpha$pcnm_alpha$PCNM3+
                    Fungal_data_alpha$pcnm_alpha$PCNM4+Fungal_data_alpha$pcnm_alpha$PCNM5+
                    Fungal_data_alpha$pcnm_alpha$PCNM9+Fungal_data_alpha$pcnm_alpha$PCNM22+
                    Fungal_data_alpha$pcnm_alpha$PCNM35+Fungal_data_alpha$pcnm_alpha$PCNM32+
                    (1|Plot_factor), 
                    data=c(Fungal_data_alpha$pcnm_alpha,Fungal_data_alpha$Env_var,Fungal_data_alpha),
                    REML=TRUE)

Fungal_best_envvar<-data.frame(Fungal_data_alpha$Env_var$Soil_moisture,
                               Fungal_data_alpha$Env_var$SLOPE,
                               Fungal_data_alpha$Env_var$KAK.eff,
                               Fungal_data_alpha$pcnm_alpha$PCNM25,
                               Fungal_data_alpha$pcnm_alpha$PCNM13,
                               Fungal_data_alpha$pcnm_alpha$PCNM34,
                               Fungal_data_alpha$pcnm_alpha$PCNM3,
                               Fungal_data_alpha$pcnm_alpha$PCNM4,
                               Fungal_data_alpha$pcnm_alpha$PCNM5,
                               Fungal_data_alpha$pcnm_alpha$PCNM9,
                               Fungal_data_alpha$pcnm_alpha$PCNM22,
                               Fungal_data_alpha$pcnm_alpha$PCNM35)
# !current implementation does not allow >=13 variables, 
# thus the last location vector was omitted in the %individual contribtion
# thus the result is not perfectly right but still reasonable compared to a
# complete ind.contr. calculation

Random_factor<-VarCorr(fungal_alpha_best) %>% 
  as.data.frame() %>% 
  mutate(icc=vcov/sum(vcov)) %>% 
  select(grp, icc)
Fixed_factor<-hier.part(fungal_rich_residuals, Fungal_best_envvar, gof = "Rsqu",
                        barplot = FALSE)

R2m<-r.squaredGLMM(fungal_alpha_best)[1]
R2c<-r.squaredGLMM(fungal_alpha_best)[2]
R_random<-Random_factor[1,2] #Sapro, Patho
#R_random<-0 #EcM #AM

Soil_perc<-Fixed_factor$I.perc$I[1]/100 #Sapro, #EcM
#Soil_perc<-0 #Pathp
#Soil_perc<-sum(Fixed_factor$I.perc$I[1:3]/100)#AM

Location_perc<-sum(Fixed_factor$I.perc$I[2:4])/100 #sapro
#Location_perc<-sum(Fixed_factor$I.perc$I[1:3])/100 #patho
#Location_perc<-sum(Fixed_factor$I.perc$I[2:6])/100 #EcM
#Location_perc<-sum(Fixed_factor$I.perc$I[4:12])/100 #AM

# Collect results
Result_order<-c(17+i,13+i,9+i, 5+i, 1+i) 
#soil,location,plot, shared(random+fixed), unexplained
#soil
alpha_model_variances[Result_order[1], 3]<-R2m*Soil_perc
#location
alpha_model_variances[Result_order[2], 3]<-R2m*Location_perc
#plot
alpha_model_variances[Result_order[3], 3]<-R_random*R2c
#shared fixed and random
alpha_model_variances[Result_order[4], 3]<-R2c-R2m-R_random*R2c
# unexplained
alpha_model_variances[Result_order[5], 3]<-1-R2c

sum(alpha_model_variances[Result_order[1:5],3])

alpha_model_variances$group<-gsub("Sapro", "Saprotrophs", alpha_model_variances$group)
alpha_model_variances$group<-gsub("Patho", "Pathogens", alpha_model_variances$group)
alpha_model_variances$group<-gsub("EcM", "EcM fungi", alpha_model_variances$group)
alpha_model_variances$group<-gsub("AM", "AM fungi", alpha_model_variances$group)

############## script output
#write.csv2(alpha_model_variances, file="./output/13_varpart_alpha_avreps.csv")

# Variance partitioning for beta diversity analysis
Sapro_beta_mrm_var<-readRDS(file="./output/10_Sapro_beta_mrm_var_avreps.RDS")
Patho_beta_mrm_var<-readRDS(file="./output/10_Patho_beta_mrm_var_avreps.RDS")
EcM_beta_mrm_var<-readRDS(file="./output/10_EcM_beta_mrm_var_avreps.RDS")
AM_beta_mrm_var<-readRDS(file="./output/10_AM_beta_mrm_var_avreps.RDS")

# Sapro models to test in variance partitioning
Sapro_stochastic<-cbind(Sapro_beta_mrm_var$Plot_beta, Sapro_beta_mrm_var$location_beta)
Sapro_soil<-cbind(Sapro_beta_mrm_var$Ctot_beta, Sapro_beta_mrm_var$CN_beta, Sapro_beta_mrm_var$CEC_beta,
                  Sapro_beta_mrm_var$BS_beta, Sapro_beta_mrm_var$moisture_beta)
Sapro_variance<-varpart(Sapro_beta_mrm_var$Fungal_beta,Sapro_stochastic,Sapro_soil)

Sapro_stochastic<-cbind(Sapro_beta_mrm_var$Plot_beta)
Sapro_soil<-cbind(Sapro_beta_mrm_var$CN_beta, Sapro_beta_mrm_var$BS_beta)
Sapro_variance<-varpart(Sapro_beta_mrm_var$Fungal_beta,Sapro_stochastic,Sapro_soil)

Patho_stochastic<-cbind(Patho_beta_mrm_var$location_beta)
Patho_tree<-cbind(Patho_beta_mrm_var$Tree_beta)
Patho_soil<-cbind(Patho_beta_mrm_var$CN_beta)
Patho_variance<-varpart(Patho_beta_mrm_var$Fungal_beta,Patho_tree,Patho_stochastic,Patho_soil)

Patho_tree<-cbind(Patho_beta_mrm_var$Tree_beta)
Patho_soil<-cbind(Patho_beta_mrm_var$CN_beta)
Patho_variance<-varpart(Patho_beta_mrm_var$Fungal_beta,Patho_tree,Patho_soil)

EcM_stochastic<-cbind(EcM_beta_mrm_var$Plot_beta,EcM_beta_mrm_var$location_beta)
EcM_tree1<-cbind(EcM_beta_mrm_var$Sample_Tree_beta, EcM_beta_mrm_var$Myco_beta, EcM_beta_mrm_var$EcM_ab_beta)
EcM_tree2<-cbind(EcM_beta_mrm_var$Sample_Tree_beta, EcM_beta_mrm_var$Myco_beta, EcM_beta_mrm_var$EcM_richness_beta)
EcM_variance1<-varpart(EcM_beta_mrm_var$Fungal_beta,EcM_stochastic,EcM_tree1)
EcM_variance2<-varpart(EcM_beta_mrm_var$Fungal_beta,EcM_stochastic,EcM_tree2)

EcM_stochastic<-cbind(EcM_beta_mrm_var$Plot_beta)
EcM_tree<-cbind(EcM_beta_mrm_var$Myco_beta)
EcM_variance<-varpart(EcM_beta_mrm_var$Fungal_beta,EcM_stochastic,EcM_tree)

AM_stochastic<-cbind(AM_beta_mrm_var$location_beta)
AM_tree<-cbind(AM_beta_mrm_var$Tree_beta)
AM_soil<-cbind(AM_beta_mrm_var$pH_beta, AM_beta_mrm_var$CN_beta, AM_beta_mrm_var$CEC_beta)
AM_variance<-varpart(AM_beta_mrm_var$Fungal_beta,AM_stochastic,AM_tree,AM_soil)

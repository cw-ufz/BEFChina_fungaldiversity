################################# Mantel test for distance matrices ###########

# R script (10/13) for analyses in Weissbecker et al. 2017
# version: April 2018

####################################
library(ecodist)  #ecodist_2.0.1  #Mantel and partial Mantel tests
library(vegan)    #vegan_2.4-4    #calculate distance matrices with vegdist
library(gdata)    #gdata_2.18.0   #extracting lower triangular matrix
library(sp)       #sp_1.2-5       #calculation geographic distances (km) 
                                  #from latitude and logitude point information
#R version 3.4.2 (2017-09-28)
####################################

source(file="./R_functions/09_Plot_beta_fungi_localtree_function.R")
source(file="./R_functions/07_Vif_function.r")

Sapro_otu<-readRDS(file="./output/03_Sapro_otu_hell_avreps.RDS")
EcM_otu<-readRDS(file="./output/03_EcM_otu_hell_avreps.RDS")
Patho_otu<-readRDS(file="./output/03_Patho_otu_hell_avreps.RDS")
AM_otu<-readRDS(file="./output/03_AM_otu_hell_avreps.RDS")

Sapro_metadata<-readRDS(file="./output/03_Sapro_metadata_trans_avreps.RDS")
EcM_metadata<-readRDS(file="./output/03_EcM_metadata_trans_avreps.RDS")
Patho_metadata<-readRDS(file="./output/03_Patho_metadata_trans_avreps.RDS")
AM_metadata<-readRDS(file="./output/03_AM_metadata_trans_avreps.RDS")

# Create dissimilarity matrices metadata

Sapro_beta<-create_beta(Sapro_otu, "none" ,"bray")
EcM_beta<-create_beta(EcM_otu, "none" ,"bray")
AM_beta<-create_beta(AM_otu, "none" ,"bray")
Patho_beta<-create_beta(Patho_otu, "none" ,"bray")


#---choose the right file manually:
Fungal_metadata<-Sapro_metadata
Fungal_beta<-Sapro_beta
#---

# Tree data (T)
Tree_beta<-create_beta(Fungal_metadata$local_tree, "none" , "bray")
Sample_Tree_beta<-create_beta(Fungal_metadata$Sample_tree, "none", "bray")
Myco_beta<-create_beta(Fungal_metadata$Myco_tree, "none", "bray")
Richness_beta<-create_beta(Fungal_metadata$richness_vars$local_richness, "none" , "euclidean")
Shannon_beta<-create_beta(Fungal_metadata$richness_vars$local_shannon, "none" , "euclidean")
Simpson_beta<-create_beta(Fungal_metadata$richness_vars$local_simpson, "none" , "euclidean")
EcM_ab_beta<-create_beta(Fungal_metadata$richness_vars$local_EcM_ab, "none" , "euclidean")
non_EcM_ab_beta<-create_beta(Fungal_metadata$richness_vars$local_non_EcM_ab, "none" , "euclidean")
EcM_richness_beta<-create_beta(Fungal_metadata$richness_vars$local_EcM_rich, "none" , "euclidean")
non_EcM_richness_beta<-create_beta(Fungal_metadata$richness_vars$local_non_EcM_rich, "none" , "euclidean")

# Environmental soil and topographic data (E)
pH_beta<-create_beta(Fungal_metadata$Env_var$pH, "none", "euclidean")
Ntot_beta<-create_beta(Fungal_metadata$Env_var$Ntot, "none", "euclidean")
Ctot_beta<-create_beta(Fungal_metadata$Env_var$Ctot, "none", "euclidean")
CN_beta<-create_beta(Fungal_metadata$Env_var$C.N, "none", "euclidean")
CEC_beta<-create_beta(Fungal_metadata$Env_var$KAK.eff, "none", "euclidean")
BS_beta<-create_beta(Fungal_metadata$Env_var$BS, "none", "euclidean")
moisture_beta<-create_beta(Fungal_metadata$Env_var$Soil_moisture, "none", "euclidean")
altitude_beta<-create_beta(Fungal_metadata$Env_var$ALTITUDE, "none", "euclidean")
slope_beta<-create_beta(Fungal_metadata$Env_var$SLOPE, "none", "euclidean")

# Design data (sample locations and plot assignment) (D)
Plot_beta<-create_beta(Fungal_metadata$plot, "none", "bray")

coordinates.sorted<-Fungal_metadata$location 
Fungal_metadata$location[ , c(1,2)] <- Fungal_metadata$location[ , c(2,1)]
colnames(Fungal_metadata$location)[c(1,2)] <- colnames(Fungal_metadata$location)[c(2,1)]
geo_distance_m <- spDists(as.matrix(Fungal_metadata$location), 
                          as.matrix(Fungal_metadata$location),  
                          longlat=TRUE)*1000
location_beta<-create_beta(geo_distance_m, "log(x+1)", "none")

fungal_beta_meta<-list(Fungal_beta,
                       pH_beta, Ntot_beta,Ctot_beta,CN_beta,CEC_beta,BS_beta, 
                       moisture_beta,altitude_beta,slope_beta,
                       Tree_beta,Sample_Tree_beta,Myco_beta,Richness_beta,Shannon_beta,
                       Simpson_beta,EcM_ab_beta,non_EcM_ab_beta,EcM_richness_beta,
                       non_EcM_richness_beta,
                       Plot_beta,location_beta)

names(fungal_beta_meta)<-list("Fungal_beta",
                       "pH_beta", "Ntot_beta","Ctot_beta","CN_beta","CEC_beta",
                       "BS_beta", "moisture_beta","altitude_beta","slope_beta",
                       "Tree_beta","Sample_Tree_beta","Myco_beta","Richness_beta",
                       "Shannon_beta", "Simpson_beta","EcM_ab_beta","non_EcM_ab_beta",
                       "EcM_richness_beta","non_EcM_richness_beta",
                       "Plot_beta","location_beta")

############################# create Mantel result table

mantel_conv_results <- matrix(data=NA, nrow = 4, ncol = 42)
row.names(mantel_conv_results)<-c("Saprotroph", "EcM-Fungi", 
                                  "Pathogen", "AM-Fungi")
mantel_conv_results<-data.frame(mantel_conv_results)

mantel_names<-c("Plot", "Geographic distance", 
                "Tree community composition", "Sampling tree identity",
                "Sampling tree mycorrhiza type","Tree richness", 
                "Tree Shannon diversity", "Tree Simpson diversity",
                "EcM tree abundance", "AM tree abundance", "EcM tree richness", 
                "AM tree richness",
                "pH", "Ntot", "Ctot", "C:N", "CEC", "BS", "Moisture", "Elevation", 
                "Slope")

pos=0
nam=1

for (cols in 1:21){
 colnames(mantel_conv_results)[(2*pos+1)]<-mantel_names[nam]
 pos=pos+1
 nam=nam+1
 colnames(mantel_conv_results)[(2*pos)]<-"p"
}

#############################

# Mantel and partial Mantel tests

#---choose the right file+index manually:
Fungal_data<-list(Sapro_beta)
fung<-1 
#---

TED_data<-list(Plot_beta,location_beta, 
               Tree_beta, Sample_Tree_beta, Myco_beta,
               Richness_beta, Shannon_beta, Simpson_beta, 
               EcM_ab_beta, non_EcM_ab_beta, EcM_richness_beta,non_EcM_richness_beta,  
               pH_beta, Ntot_beta, Ctot_beta, CN_beta, CEC_beta, BS_beta, moisture_beta,
               altitude_beta, slope_beta)

detach(package:vegan) # vegan also has a Mantel function that needs to be deactivated

# Ensure datasets and col names result matrix are the same order
for(groups in Fungal_data){
 evar<-1
  for(dats in TED_data){
   print("Start Mantel calculation"); print(evar)
   mantel_conventional<-mantel(groups~dats,nperm=10000,mrank=FALSE)
   print("Start partial Mantel calculation")
   mantel_partial<-mantel(groups~dats+location_beta,nperm=10000,mrank=FALSE)
   print("Start MRM calculation")
   mrm.single<-MRM(groups ~ dats, nperm=10000,mrank=FALSE)
   
   mantel_conv_results[fung,evar]<-mantel_conventional[1]
   mantel_partial_results[fung,evar]<-mantel_partial[1]
   MRM_single_results[fung,evar]<-mrm.single$F.test[1]
   
   evar=evar+1
   mantel_conv_results[fung,evar]<-mantel_conventional[2]
   mantel_partial_results[fung,evar]<-mantel_partial[2]
   MRM_single_results[fung,evar]<-mrm.single$F.test[2]
   evar=evar+1
   }
 }

# Check for variance inflation factor= highly correlated parameters
Sapro_beta_meta_data<-readRDS(file="./output/10_Sapro_beta_meta_avreps.RDS")
Sapro_mantel_var<-list( Sapro_beta_meta_data$Fungal_beta,
                        Sapro_beta_meta_data$Plot_beta,
                        Sapro_beta_meta_data$location_beta,
                        Sapro_beta_meta_data$pH_beta,
                        Sapro_beta_meta_data$Ctot_beta, 
                        Sapro_beta_meta_data$CN_beta,
                        Sapro_beta_meta_data$CEC_beta,
                        Sapro_beta_meta_data$BS_beta,
                        Sapro_beta_meta_data$moisture_beta)
names(Sapro_mantel_var)<-c("Fungal_beta","Plot_beta","location_beta","pH_beta",
                           "Ctot_beta", "CN_beta", "CEC_beta","BS_beta",
                           "moisture_beta")

Sapro_vif_result<-vif_func(Sapro_mantel_var[-1],thresh=10,trace=T)

Patho_beta_meta_data<-readRDS(file="./output/10_Patho_beta_meta_avreps.RDS")
Patho_mantel_var<-list(Patho_beta_meta_data$Fungal_beta,
                        Patho_beta_meta_data$Plot_beta,
                        Patho_beta_meta_data$location_beta, 
                        Patho_beta_meta_data$Tree_beta,
                        Patho_beta_meta_data$Ctot_beta,
                        Patho_beta_meta_data$CN_beta)
names(Patho_mantel_var)<-c("Fungal_beta","Plot_beta", "location_beta", "Tree_beta",
                           "Ctot_beta", "CN_beta")
Patho_vif_result<-vif_func(Patho_mantel_var[-1],thresh=10,trace=T)

EcM_beta_meta_data<-readRDS(file="./output/10_EcM_beta_meta_avreps.RDS")
EcM_mantel_var<-list(EcM_beta_meta_data$Fungal_beta,
                      EcM_beta_meta_data$Plot_beta,
                      EcM_beta_meta_data$location_beta,
                      EcM_beta_meta_data$Sample_Tree_beta,
                      EcM_beta_meta_data$Myco_beta,
                      EcM_beta_meta_data$EcM_ab_beta,
                      EcM_beta_meta_data$EcM_richness_beta)
names(EcM_mantel_var)<-c("Fungal_beta","Plot_beta", "location_beta","Sample_Tree_beta",
                         "Myco_beta", "EcM_ab_beta", "EcM_richness_beta")
EcM_vif_result<-vif_func(list(EcM_mantel_var[-1]),thresh=10,trace=T)

AM_beta_meta_data<-readRDS(file="./output/10_AM_beta_meta_avreps.RDS")
AM_mantel_var<-list(AM_beta_meta_data$Fungal_beta,
                    AM_beta_meta_data$Plot_beta,
                    AM_beta_meta_data$location_beta,
                    AM_beta_meta_data$Tree_beta,
                    AM_beta_meta_data$pH_beta,
                    AM_beta_meta_data$Ctot_beta,
                    AM_beta_meta_data$CN_beta,
                    AM_beta_meta_data$CEC_beta,
                    AM_beta_meta_data$BS_beta)
names(AM_mantel_var)<-c("Fungal_beta","Plot_beta", "location_beta", "Tree_beta",
                        "pH_beta", "Ctot_beta",
                        "CN_beta", "CEC_beta", "BS_beta")
AM_vif_result<-vif_func(AM_mantel_var[-1],thresh=10,trace=T) 


############## script output 

#---choose the right file manually:

# Save beta variables
#saveRDS(fungal_beta_meta, file="./output/10_Sapro_beta_meta_avreps.RDS")
#saveRDS(fungal_beta_meta, file="./output/10_Patho_beta_meta_avreps.RDS")
#saveRDS(fungal_beta_meta, file="./output/10_EcM_beta_meta_avreps.RDS")
#saveRDS(fungal_beta_meta, file="./output/10_AM_beta_meta_avreps.RDS")

# Mantel results
#saveRDS(mantel_conv_results, file="./output/10_mantel_conv_avrep.RDS")
#saveRDS(mantel_partial_results, file="./output/10_mantel_part_avrep.RDS")
#saveRDS(MRM_single_results, file="./output/10_mrm_single_avrep.RDS")

#write.csv2(mantel_conv_results, file="./output/10_mantel_conv_avrep.csv")
#write.csv2(mantel_partial_results, file="./output/10_mantel_part_avrep.csv")
#write.csv2(MRM_single_results, file="./output/10_mrm_single_avrep.csv")

# Variables for mrm analysis
#saveRDS(AM_mantel_var, "./output/10_AM_beta_mrm_var_avreps.RDS")
#saveRDS(Sapro_mantel_var, "./output/10_Sapro_beta_mrm_var_avreps.RDS")
#saveRDS(EcM_mantel_var, "./output/10_EcM_beta_mrm_var_avreps.RDS")
#saveRDS(Patho_mantel_var, "./output/10_Patho_beta_mrm_var_avreps.RDS")


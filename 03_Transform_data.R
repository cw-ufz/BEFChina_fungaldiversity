############################### Transform data ################################

# R script (03/13) for analyses in Weissbecker et al. 2018
# version: August 2018

###############################
library(dplyr)      # dplyr_0.7.6   #data transformation
#sessionInfo()
#R version 3.5.1 (2018-07-02)
################

BEF_Sapro<-readRDS(file="./Cleaned_Data/02_BEF_dom_Sapro_min_250.RDS")
BEF_EcM<-readRDS(file="./Cleaned_Data/02_BEF_dom_EcM_min_210.RDS")
BEF_Patho<-readRDS(file="./Cleaned_Data/02_BEF_dom_Patho_min_40.RDS")
BEF_AM<-readRDS(file="./Cleaned_Data/02_BEF_dom_AM_min_20.RDS")

Sapro_var<-readRDS(file="./Cleaned_Data/02_BEF_dom_Sapro_varlist_min_250.RDS")
EcM_var<-readRDS(file="./Cleaned_Data/02_BEF_dom_EcM_varlist_min_210.RDS")
Patho_var<-readRDS(file="./Cleaned_Data/02_BEF_dom_Patho_varlist_min_40.RDS")
AM_var<-readRDS(file="./Cleaned_Data/02_BEF_dom_AM_varlist_min_20.RDS")

### Transformation OTU abundance counts
Sapro_hell_otu<-decostand(data.frame(t(otu_table(BEF_Sapro))), 
                          method="hellinger")
Sapro_otu_pa<-specnumber(data.frame(t(otu_table(BEF_Sapro))))
range(Sapro_otu_pa) #25 136
mean(Sapro_otu_pa)  #78.18663

EcM_hell_otu<-decostand(data.frame(t(otu_table(BEF_EcM))), 
                        method="hellinger")
EcM_otu_pa<-specnumber(data.frame(t(otu_table(BEF_EcM))))
range(EcM_otu_pa) #1 25
mean(EcM_otu_pa)  #9.81383

Patho_hell_otu<-decostand(data.frame(t(otu_table(BEF_Patho))), 
                          method="hellinger")
Patho_otu_pa<-specnumber(data.frame(t(otu_table(BEF_Patho))))
range(Patho_otu_pa) #5 44
mean(Patho_otu_pa)  #20.53107

AM_hell_otu<-decostand(data.frame(t(otu_table(BEF_AM))), 
                       method="hellinger")
AM_otu_pa<-specnumber(data.frame(t(otu_table(BEF_AM))))
range(AM_otu_pa) #2 43
mean(AM_otu_pa)  #14.7243

### Transform metadata

#---need to put applicable dataset manually here:
Data_var<-Sapro_var
#---

var_trans<-Data_var$Env_var %>%
  mutate_at(vars(Ntot:SLOPE), log1p) #log1p() calculates log(1+x)

var_rich_trans<-Data_var$richness_vars %>%
  mutate_at(vars(-local_shannon,-local_simpson), log1p) 

Data_single_trans<-list(var_trans,Data_var$plot,Data_var$Sample_tree,
                        Data_var$Myco_tree, Data_var$local_tree, 
                        var_rich_trans, Data_var$location)

names(Data_single_trans)<-c("Env_var", "plot", "Sample_tree", 
                            "Myco_tree", "local_tree", "richness_vars", 
                            "location")

#### Averaging sample replicates per tree species for each plot

#---need to put applicable dataset manually here:
Var_Dataset<-AM_var
Dataset<-AM_hell_otu
otu_read_datase<-BEF_AM
#---

#Extract plot ID and tree ID information for proper averaging
plot_matrix<-Var_Dataset$plot
tree_matrix<-Var_Dataset$Sample_tree

plot_ID_prep<-cbind(rownames(plot_matrix),plot_matrix)
plot_pre_red<-melt(plot_ID_prep, id.vars ="rownames(plot_matrix)")
plot_ID<-plot_pre_red[plot_pre_red$value>0,-3]
names(plot_ID)<-c("Sample", "Plot_ID")

tree_ID_prep<-cbind(rownames(tree_matrix),tree_matrix)
tree_pre_red<-melt(tree_ID_prep, id.vars ="rownames(tree_matrix)")
tree_ID<-tree_pre_red[tree_pre_red$value>0,-3]
names(tree_ID)<-c("Sample", "Tree_ID")

plot_tree_ID<-merge(plot_ID,tree_ID, by="Sample", all.x=TRUE)
rownames(plot_tree_ID)<-plot_tree_ID$Sample

 otu_hell<-Dataset
 
 plot_tree_ID.sorted<-plot_tree_ID[match(rownames(otu_hell), 
                                         rownames(plot_tree_ID)),-1]
 merged_otu<-merge(plot_tree_ID.sorted,otu_hell,by="row.names",all.x=TRUE)
 
 otu_mean<-merged_otu[,-1] %>%
  group_by(Plot_ID, Tree_ID) %>%
  summarise_all(mean)

#Determine read sums sample replicates
 otu_reads<-data.frame(t(otu_table(otu_read_datase)))
 
 plot_tree_ID.sorted<-plot_tree_ID[match(rownames(otu_reads), 
                                         rownames(plot_tree_ID)),-1]
 merged_otu<-merge(plot_tree_ID.sorted,otu_reads,by="row.names",all.x=TRUE)
 
 otu_sum<-merged_otu[,-1] %>%
   group_by(Plot_ID, Tree_ID) %>%
   summarise_all(sum)
 
# Average local tree community
local_hell<-decostand(Var_Dataset$local_tree, method="hellinger")

plot_tree_ID.sorted<-plot_tree_ID[match(rownames(local_hell), 
                                        rownames(plot_tree_ID)),-1]
merged_env<-merge(plot_tree_ID.sorted,local_hell,by="row.names",all.x=TRUE)

var_mean<-merged_env[,-1] %>%
  group_by(Plot_ID, Tree_ID) %>% 
  summarise_all(mean)
local_tree<-var_mean[,c(-1,-2)]

# Average local soil + typographic variables
plot_tree_ID.sorted<-plot_tree_ID[match(rownames(Var_Dataset$Env_var),
                                        rownames(plot_tree_ID)),-1]
merged_env<-merge(plot_tree_ID.sorted,Var_Dataset$Env_var,by="row.names",all.x=TRUE)

var_trans<-merged_env %>%
  mutate_at(vars(Ntot:SLOPE), log1p)

var_mean<-var_trans[,-1] %>%
  group_by(Plot_ID, Tree_ID) %>%
  summarise_all(mean)

Env_var<-var_mean[,c(-1,-2,-3)]

# Average plot information
plot_tree_ID.sorted<-plot_tree_ID[match(rownames(Var_Dataset$plot), 
                                        rownames(plot_tree_ID)),-1]
merged_plot<-merge(plot_tree_ID.sorted,Var_Dataset$plot,by="row.names",all.x=TRUE)

plot_mean<-merged_plot[,-1] %>%
  group_by(Plot_ID, Tree_ID) %>% 
  summarise_all(mean)
plot<-plot_mean[,c(-1,-2)]

# Average Sample_tree information
plot_tree_ID.sorted<-plot_tree_ID[match(rownames(Var_Dataset$Sample_tree), 
                                        rownames(plot_tree_ID)),-1]
merged_sample_tree<-merge(plot_tree_ID.sorted,Var_Dataset$Sample_tree,by="row.names",all.x=TRUE)

sample_tree_mean<-merged_sample_tree[,-1] %>%
  group_by(Plot_ID, Tree_ID) %>% 
  summarise_all(mean)
Sample_tree<-sample_tree_mean[,c(-1,-2)]

# Average Myco tree information
plot_tree_ID.sorted<-plot_tree_ID[match(rownames(Var_Dataset$Myco_tree),
                                        rownames(plot_tree_ID)),-1]
merged_myco_tree<-merge(plot_tree_ID.sorted,Var_Dataset$Myco_tree,by="row.names",all.x=TRUE)

myco_tree_mean<-merged_myco_tree[,-1] %>%
  group_by(Plot_ID, Tree_ID) %>% 
  summarise_all(mean)
Myco_tree<-myco_tree_mean[,c(-1,-2)]

# Average sample location information
plot_tree_ID.sorted<-plot_tree_ID[match(rownames(Var_Dataset$location), 
                                        rownames(plot_tree_ID)),-1]
merged_location<-merge(plot_tree_ID.sorted,Var_Dataset$location,by="row.names",all.x=TRUE)

location_mean<-merged_location[,-1] %>%
  group_by(Plot_ID, Tree_ID) %>% 
  summarise_all(mean)
location<-location_mean[,c(-1,-2)]

# Average tree richness information
plot_tree_ID.sorted<-plot_tree_ID[match(rownames(Var_Dataset$richness_vars), 
                                        rownames(plot_tree_ID)),-1]
merged_richness_tree<-merge(plot_tree_ID.sorted,Var_Dataset$richness_vars,
                            by="row.names",all.x=TRUE)

richness_trans<-merged_richness_tree %>%
  mutate_at(vars(local_richness,local_non_EcM_ab:local_non_EcM_rich), log1p) 

richness_mean<-richness_trans[,-1] %>%
  group_by(Plot_ID, Tree_ID) %>% 
  summarise_all(mean) 

richness_vars<-richness_mean[,c(-1,-2)]

# Summarize to new variable list for averaged replicates per tree species pre plot
Var_Dataset_avrep<-list(Env_var,plot,Sample_tree,
                        Myco_tree, local_tree, richness_vars,location)

names(Var_Dataset_avrep)<-c("Env_var", "plot", "Sample_tree", 
                            "Myco_tree", "local_tree", "richness_vars", "location")

############## script output

#---choose the right file manually:

# saveRDS(otu_mean[,c(-1,-2)], file="./output/03_Sapro_otu_hell_avreps.RDS")
# saveRDS(otu_mean[,c(-1,-2)], file="./output/03_EcM_otu_hell_avreps.RDS")
# saveRDS(otu_mean[,c(-1,-2)], file="./output/03_Patho_otu_hell_avreps.RDS")
# saveRDS(otu_mean[,c(-1,-2)], file="./output/03_AM_otu_hell_avreps.RDS")
# 
# saveRDS(otu_sum[c(-1,-2)], file="./output/03_Sapro_otu_sum_avreps.RDS")
# saveRDS(otu_sum[c(-1,-2)], file="./output/03_EcM_otu_sum_avreps.RDS")
# saveRDS(otu_sum[c(-1,-2)], file="./output/03_Patho_otu_sum_avreps.RDS")
# saveRDS(otu_sum[c(-1,-2)], file="./output/03_AM_otu_sum_avreps.RDS")
#
#saveRDS(Var_Dataset_avrep, file="./output/03_Sapro_metadata_trans_avreps.RDS")
#saveRDS(Var_Dataset_avrep, file="./output/03_EcM_metadata_trans_avreps.RDS")
#saveRDS(Var_Dataset_avrep, file="./output/03_Patho_metadata_trans_avreps.RDS")
#saveRDS(Var_Dataset_avrep, file="./output/03_AM_metadata_trans_avreps.RDS")


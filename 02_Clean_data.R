############################### Clean data ####################################

# R script (02/13) for analyses in Weissbecker et al. 2018 New Phytologist
# version: April 2018

###############################

library(vegan)  # vegan_2.4-4 #decostand function
#R version 3.4.2 (2017-09-28)
###############################

BEF_dom<-readRDS(file="./Cleaned_Data/01_dominant_phyloseq.RDS")
all_var_sorted<-readRDS(file="./Cleaned_Data/01_dominant_metadata.RDS")

BEF_dom_fungi <- subset_taxa(BEF_dom, Kingdom=="Fungi") 

# 2) remove samples with unknown tree neighborhood (3 samples)
rm_meta_samples<-unique(row.names(which(is.na(all_var_sorted$local_tree),
                                        arr.ind=TRUE)))
BEF_fungi_valid<-subset_samples(BEF_dom_fungi,
                                !rownames(sample_data(BEF_dom_fungi)) %in% 
                                  rm_meta_samples)

# 3) split dataset according to fungal functional groups
tax<-data.frame(tax_table(BEF_dom))
unique(tax)

BEF_Sapro_pre = subset_taxa(BEF_fungi_valid, Fungal_guild=="Saprotroph")#1768
BEF_EcM_pre = subset_taxa(BEF_fungi_valid, Fungal_guild=="Ectomycorrhizal")#410
BEF_AM_pre = subset_taxa(BEF_fungi_valid, Fungal_guild=="Arbuscular Mycorrhizal") #320
BEF_Plant_Patho_pre = subset_taxa(BEF_fungi_valid, Fungal_guild=="Plant Pathogen")#310
BEF_Unknown = subset_taxa(BEF_fungi_valid, Fungal_guild=="Unknown") #2604
BEF_Erico = subset_taxa(BEF_fungi_valid,Fungal_guild=="Ericoid Mycorrhizal")#34
BEF_Endo = subset_taxa(BEF_fungi_valid, Fungal_guild=="Endophyte")#75
BEF_Lich = subset_taxa(BEF_fungi_valid, Fungal_guild=="Lichenized")#36
BEF_Fungal_Patho = subset_taxa(BEF_fungi_valid,Fungal_guild=="Fungal Parasite")# 46
BEF_Lichen_Patho = subset_taxa(BEF_fungi_valid,Fungal_guild=="Lichen Parasite")# 1
BEF_Animal_Patho = subset_taxa(BEF_fungi_valid,Fungal_guild=="Animal Pathogen")#48
BEF_Epiphyt = subset_taxa(BEF_fungi_valid, Fungal_guild=="Epiphyte")# 8
BEF_OrchidM = subset_taxa(BEF_fungi_valid, Fungal_guild=="Orchid Mycorrhizal")# 5

#OTUs annotated to fungal functional group:
#1-2604/5665=0.5403354 classified to any functional group
#1768/5665=0.31 saprotroph
#410/5665=0.07 EcM
#320/5665=0.05 AM
#310/5665=0.05 Patho

#Pathogen composition
#310+46+1+48 = 405
#310/405=0.7654321 #plant pathogen
#46/405=0.1135802 # fungal parasite
#1/405=0.002469136 #lichen parasite
#48/405=0.1185185 #animal pathogen

# 4)Extract samples with a minimum threshold community of target fungal group
# >= appr. 10% of max sequence reads achieved for respective groups
# extract metadata subset for each fungal main group

max(sample_sums(BEF_Sapro_pre))       # 2764 reads
BEF_Sapro = prune_samples(sample_sums(BEF_Sapro_pre)>=250, BEF_Sapro_pre)           #361 samples
max(sample_sums(BEF_EcM_pre))         # 2180 reads
BEF_EcM = prune_samples(sample_sums(BEF_EcM_pre)>=210, BEF_EcM_pre)                 #190 samples
max(sample_sums(BEF_Plant_Patho_pre)) # 430 reads
BEF_Patho = prune_samples(sample_sums(BEF_Plant_Patho_pre)>=40,BEF_Plant_Patho_pre) #178 samples
max(sample_sums(BEF_AM_pre))          # 201 reads
BEF_AM = prune_samples(sample_sums(BEF_AM_pre)>=20, BEF_AM_pre)                     #215 samples

Sapro_var<-lapply(all_var_sorted, function(x) x[rownames(x) %in% rownames(sample_data(BEF_Sapro)),])
EcM_var<-lapply(all_var_sorted, function(x) x[rownames(x) %in% rownames(sample_data(BEF_EcM)),])
Patho_var<-lapply(all_var_sorted, function(x) x[rownames(x) %in% rownames(sample_data(BEF_Patho)),])
AM_var<-lapply(all_var_sorted, function(x) x[rownames(x) %in% rownames(sample_data(BEF_AM)),])


############## script output

saveRDS(BEF_Sapro, file="./Cleaned_Data/02_BEF_dom_Sapro_min_250.RDS")
saveRDS(BEF_EcM, file="./Cleaned_Data/02_BEF_dom_EcM_min_210.RDS")
saveRDS(BEF_Patho, file="./Cleaned_Data/02_BEF_dom_Patho_min_40.RDS")
saveRDS(BEF_AM, file="./Cleaned_Data/02_BEF_dom_AM_min_20.RDS")

saveRDS(Sapro_var, file="./Cleaned_Data/02_BEF_dom_Sapro_varlist_min_250.RDS")
saveRDS(EcM_var, file="./Cleaned_Data/02_BEF_dom_EcM_varlist_min_210.RDS")
saveRDS(Patho_var, file="./Cleaned_Data/02_BEF_dom_Patho_varlist_min_40.RDS")
saveRDS(AM_var, file="./Cleaned_Data/02_BEF_dom_AM_varlist_min_20.RDS")

write.csv2(data.frame(tax_table(BEF_Sapro)), file="./output/02_BEF_dom_Sapro_taxtable.csv")
write.csv2(data.frame(tax_table(BEF_EcM)), file="./output/02_BEF_dom_EcM_taxtable.csv")
write.csv2(data.frame(tax_table(BEF_Patho)), file="./output/02_BEF_dom_Patho_taxtable.csv")
write.csv2(data.frame(tax_table(BEF_AM)), file="./output/02_BEF_dom_AM_taxtable.csv")

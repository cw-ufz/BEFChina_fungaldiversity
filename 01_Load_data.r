############################### Load data #####################################

# R script (01/13) for analyses in Weissbecker et al. 2018 New Phytologist
# version: April 2018

#########################

library(phyloseq)         #phyloseq_1.20.0        #merge otu,tax,metadata info
library(vegan)            #vegan_2.4-4            #diversity calculations
library(data.table)       #data.table_1.10.4-1    #dcast function
library(biomformat)       #biomformat_1.4.0       #read biom file
#R version 3.4.2 (2017-09-28)
#Platform: x86_64-w64-mingw32/x64 (64-bit)
#Running under: Windows 7 x64 (build 7601) Service Pack 1

######################### install biomformat package
# source("https://bioconductor.org/biocLite.R")
# biocLite("rhdf5")
# biocLite("biomformat")
######################### 

#1) Load otu_table, prune to dominant OTU dataset (reads>=4)
x = read_biom("./Data/BEF_vipA_454_ITS_1283866816.json_Nov2017.biom") 

otu_matrix<-as(biom_data(x), "matrix")                               
BEF_all_otu<- phyloseq(otu_table(otu_matrix, taxa_are_rows = TRUE))
BEF_dom_otu = prune_taxa(taxa_sums(BEF_all_otu) > 3, BEF_all_otu)

#2) Load taxonomic information and fungal group assignment, create phyloseq object
fungal_taxonomy<-read.csv2(file="./Data/BEF_vipA_454_ITS_2011_curated_taxonomy_guild.csv",
                           row.names = 1)
taxonomy_fungi_pre<-as.matrix(fungal_taxonomy)
row.names(taxonomy_fungi_pre)<-rownames(fungal_taxonomy)
taxonomy_fungi<-taxonomy_fungi_pre 

BEF_dom<- phyloseq(otu_table(BEF_dom_otu, taxa_are_rows = TRUE),
                   tax_table(taxonomy_fungi))

#3) Load metadata
Env_var<-read.csv("./Data/BEF_vipA_2011_soil_typo_metadata.csv",row.names=1, 
                  dec=",", sep=";") 
Sample_tree<- read.csv("./Data/BEF_vipA_2011_Sampling_tree.csv",row.names=1, 
                       na.strings = "NA", sep=";")
Myco_tree<- read.csv("./Data/BEF_vipA_2011_Tree_mycotype_composition.csv",
                     row.names=1, na.strings = "NA", sep=";")
local_tree<-read.csv("./Data/BEF_vipA_2011_local_9_trees.csv", na.strings = "NA",
                     row.names =1, sep=";")
location <- read.csv("./Data/BEF_vipA_2011_sample_location.csv",row.names=1,
                     na.strings = "NA")

# get wide format plot information to create binary plot information
plot_long<-data.frame(cbind(row.names(Env_var), as.character(Env_var$Plot), 
                            rep(1,length(Env_var$Plot)))) 
names(plot_long)<-c("Sample_ID", "Plot_ID", "value")
plot_long$value<-as.numeric(plot_long$value)
plot_wide <- dcast(plot_long, Sample_ID ~ Plot_ID ,value.var="value", fill=0)
rownames(plot_wide)<-plot_wide$Sample_ID
plot_wide[]<-data.frame(sapply(plot_wide, as.numeric))
plot<-plot_wide[,-1]

# load calculate tree diversity information
local_richness<-data.frame(specnumber(local_tree, MARGIN=1))
local_shannon<-data.frame(diversity(local_tree, index="shannon", MARGIN = 1))
local_simpson<-data.frame(diversity(local_tree, index="simpson", MARGIN = 1))
local_non_EcM_ab<-data.frame(rowSums(local_tree[,c(2,5,3,10,11,7,12,8)]))
local_EcM_ab<-data.frame(rowSums(local_tree[,c(1,4,6,7,9,13,14,15,16,17)]))
local_EcM_rich<-data.frame(specnumber(local_tree[,c(1,4,6,7,9,13,14,15,16,17)],
                                      MARGIN=1))
local_non_EcM_rich<-data.frame(specnumber(local_tree[,c(2,5,3,10,11,7,12,8)],
                                          MARGIN=1))

#Assure correct data order
reference<-t(otu_matrix)
richness_vars<-cbind.data.frame(local_richness,local_shannon,local_simpson,
                                local_non_EcM_ab,local_EcM_ab,local_EcM_rich,
                                local_non_EcM_rich)
colnames(richness_vars)<-c("local_richness","local_shannon","local_simpson",
                           "local_non_EcM_ab","local_EcM_ab","local_EcM_rich",
                           "local_non_EcM_rich")

var_list<-list(Env_var,plot,Sample_tree,Myco_tree,local_tree,richness_vars,
               location)
names(var_list)<-list("Env_var","plot","Sample_tree","Myco_tree","local_tree",
                      "richness_vars","location")
all_var_sorted<-lapply(var_list, function(x) x[match(rownames(reference),
                                                     rownames(x)),]) 

BEF_dom_addmeta<- phyloseq(otu_table(BEF_dom_otu, taxa_are_rows = TRUE), 
                           tax_table(taxonomy_fungi), 
                           sample_data(all_var_sorted$Env_var))

#### Summary statistics

BEF_dom_fungi <- subset_taxa(BEF_dom, Kingdom=="Fungi") 
taxonomy_readable<-data.frame(tax_table(BEF_dom_fungi))
fungi_taxa<-length(row.names(taxonomy_readable))

Kingdom<-which(colnames(taxonomy_readable)== "Kingdom")
Phylum<-which(colnames(taxonomy_readable)== "Phylum")
Class<-which(colnames(taxonomy_readable)== "Class")
Order<-which(colnames(taxonomy_readable)== "Order")
Family<-which(colnames(taxonomy_readable)== "Family")
Genus<-which(colnames(taxonomy_readable)== "Genus")

Fungi_Phylum<-taxonomy_readable[which(taxonomy_readable$Kingdom=="Fungi"),
                                Phylum]
out_Phylum<-(which(is.na(Fungi_Phylum)))
Taxa_Phylum<-unique(Fungi_Phylum[-out_Phylum])
phyla_taxs<-as.list(Taxa_Phylum)
phyla_otu_numbers<-lapply(phyla_taxs, 
                          function(x)length(which(taxonomy_readable[,Phylum]==x)))
names(phyla_otu_numbers)<-c(Taxa_Phylum)

Fungi_Class<-taxonomy_readable[which(taxonomy_readable$Kingdom=="Fungi"),Class]
out_Class<-(which(is.na(Fungi_Class)))
Taxa_Class<-unique(Fungi_Class[-out_Class])
class_taxs<-as.list(Taxa_Class)
class_otu_numbers<-lapply(class_taxs, 
                          function(x)length(which(taxonomy_readable[,Class]==x)))
names(class_otu_numbers)<-c(Taxa_Class)

Fungi_Order<-taxonomy_readable[which(taxonomy_readable$Kingdom=="Fungi"),Order]
out_Order<-(which(is.na(Fungi_Order)))
Taxa_Order<-unique(Fungi_Order[-out_Order])
order_taxs<-as.list(Taxa_Order)
order_otu_numbers<-lapply(order_taxs, 
                          function(x)length(which(taxonomy_readable[,Order]==x)))
names(order_otu_numbers)<-c(Taxa_Order)

Fungi_Family<-taxonomy_readable[,Family]
out_Family<-(which(is.na(Fungi_Family)))
Taxa_Family<-unique(Fungi_Family[-out_Family])
family_taxs<-as.list(Taxa_Family)
family_otu_numbers<-lapply(family_taxs, 
                           function(x)length(which(taxonomy_readable[,Family]==x)))
names(family_otu_numbers)<-c(Taxa_Family)

Fungi_Genus<-taxonomy_readable[,Genus]
out_Genus<-(which(is.na(Fungi_Genus)))
Taxa_Genus<-unique(Fungi_Genus[-out_Genus])
genus_taxs<-as.list(Taxa_Genus)
genus_otu_numbers<-lapply(genus_taxs, 
                          function(x)length(which(taxonomy_readable[,Genus]==x)))
names(genus_otu_numbers)<-c(Taxa_Genus)

# Percentage taxonomic identification at each taxonomic level
n_taxa<-ntaxa(BEF_dom_addmeta)

pKingdom<-length(which(taxonomy_readable$Kingdom=="Fungi"))/n_taxa*100 
pPhylum<-(fungi_taxa-length(out_Phylum))/n_taxa*100 
pClass<-(fungi_taxa-length(out_Class))/n_taxa*100 
pOrder<-(fungi_taxa-length(out_Order))/n_taxa*100
pFamily<-(fungi_taxa-length(out_Family))/n_taxa*100
pGenus<-(fungi_taxa-length(out_Genus))/n_taxa*100



############## script output

saveRDS(BEF_dom_addmeta, 
        file="./Cleaned_Data/01_dominant_phyloseq.RDS")
saveRDS(all_var_sorted, 
        file="./Cleaned_Data/01_dominant_metadata.RDS")

write.csv2(data.frame(unlist(phyla_otu_numbers)), 
           file="./output/01_taxonomy_phylum_sums.csv")
write.csv2(data.frame(unlist(class_otu_numbers)), 
           file="./output/01_taxonomy_class_sums.csv")
write.csv2(data.frame(unlist(order_otu_numbers)), 
           file="./output/01_taxonomy_order_sums.csv")
write.csv2(data.frame(unlist(family_otu_numbers)), 
           file="./output/01_taxonomy_family_sums.csv")
write.csv2(data.frame(unlist(genus_otu_numbers)), 
           file="./output/01_taxonomy_genus_sums.csv")



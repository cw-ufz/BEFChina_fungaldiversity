############################### Rarefaction curves ############################

# R script (05/13) for analyses in Weissbecker et al. 2018
# version: August 2018

###############################

library(extrafont)  # extrafont_0.17    #make available font "Arial"
library(vegan)      # vegan_2.5-2 
#font_import()      # install fonts if neccessary
#fonts()
#sessionInfo()
#R version 3.5.1 (2018-07-02)
#######################

BEF_Sapro<-readRDS(file="./Cleaned_Data/02_BEF_dom_Sapro_min_250.RDS")
BEF_EcM<-readRDS(file="./Cleaned_Data/02_BEF_dom_EcM_min_210.RDS")
BEF_Patho<-readRDS(file="./Cleaned_Data/02_BEF_dom_Patho_min_40.RDS")
BEF_AM<-readRDS(file="./Cleaned_Data/02_BEF_dom_AM_min_20.RDS")

Sapro_otu<-otu_table(BEF_Sapro)
Sapro_otu_t<-t(Sapro_otu)
Sapro_otu_df<-data.frame(Sapro_otu_t)

Patho_otu<-otu_table(BEF_Patho)
Patho_otu_t<-t(Patho_otu)
Patho_otu_df<-data.frame(Patho_otu_t)

AM_otu<-otu_table(BEF_AM)
AM_otu_t<-t(AM_otu)
AM_otu_df<-data.frame(AM_otu_t)

EcM_otu<-otu_table(BEF_EcM)
EcM_otu_t<-t(EcM_otu)
EcM_otu_df<-data.frame(EcM_otu_t)


############## script output

tiff(filename ="./output/images/05_Rarefaction_curves_minreads.tiff",
     width = 6.8, height = 6.8, units = "in", res=300)
par(mfrow=c(2,2))
par(mar=c(4.5, 5, 1, 0.5))
rarecurve(Sapro_otu_df, step = 10, sample = 3000, label = FALSE, 
          ylab="Saprotroph OTUs")
rarecurve(Patho_otu_df, step = 1, sample = 500, label = FALSE, 
          ylab="Fungal plant pathogen OTUs")
rarecurve(AM_otu_df, step = 1, sample = 300, label = FALSE,
          ylab="AM fungal OTUs")
rarecurve(EcM_otu_df, step = 1, sample = 2500, label = FALSE,
          ylab="EcM fungal OTUs")
dev.off()

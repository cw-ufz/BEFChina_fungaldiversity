####### Procrustes correlation between complete and dominant OTU dataset #######

# R script (04/13) for analyses in Weissbecker et al. 2018
# version: August 2018

######################
library(vegan)            #vegan_2.5-2
library(biomformat)       #biomformat_1.8.0 #read biom file
library(phyloseq)         #phyloseq_1.24.0 
#sessionInfo()
#R version 3.5.1 (2018-07-02)
######################

x = read_biom("./Data/BEF_vipA_454_ITS_1283866816.json_Nov2017.biom")
otu_matrix<-as(biom_data(x), "matrix")                               
BEF_all_otu<- phyloseq(otu_table(otu_matrix, taxa_are_rows = TRUE))
BEF_dom_otu = prune_taxa(taxa_sums(BEF_all_otu) > 3, BEF_all_otu)

all_table<-BEF_all_otu  # 21732 taxa
dom_table<-BEF_dom_otu  # 6693 taxa

### NMDS calculation for the complete dataset 
#(including singletons, doubletons and tripletons)

otu_all<-data.frame(all_table, check.names = FALSE) 
t_all<-t(otu_all)

otu_hell_all<-decostand(t_all,method="hellinger") 

otu_nmds_all<-metaMDS(otu_hell_all, distance = "bray", k = 2, trymax = 30,
                      autotransform =FALSE)
otu_nmds_all$stress #0.1597039 

### NMDS calculation for the dominant OTU dataset 
#(excluding singletons, doubletons and tripletons)

dom_otu<-data.frame(dom_table, check.names = FALSE)
dom_t<-t(dom_otu) 

otu_hell<-decostand(dom_t,method="hellinger") 

otu_nmds<-metaMDS(otu_hell, distance = "bray", k = 2, trymax = 30, 
                  autotransform =FALSE)
otu_nmds$stress # 0.1614698 

### Procrustes correlation analysis of both NMDS ordinations

pro<-protest(otu_nmds_all, otu_nmds)
#Procrustes Sum of Squares (m12 squared):        0.002511 
#Correlation in a symmetric Procrustes rotation: 0.9987   
#Significance:  0.001 
#Number of permutations: 999


############## script output

tiff(filename = "./output/images/04_Procrustes_all_vs_dom_otu.tiff", 
     width = 5, height = 6.8, units = "in", res=300) 
plot(pro, cex=2)
dev.off()





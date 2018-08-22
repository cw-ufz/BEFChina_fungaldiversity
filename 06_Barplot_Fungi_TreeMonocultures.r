############# Barplots fungal functional groups in monocultures ###############

# R script (06/13) for analyses in Weissbecker et al. 2018
# version: August 2018

###############################

library(biomformat) # biomformat_1.8.0
library(phyloseq)   # phyloseq_1.24.0 
library(dplyr)      # dplyr_0.7.6
library(ggplot2)    # ggplot2_3.0.0
library(extrafont)  # extrafont_0.17
  #font_import()
  #fonts()
library(scales)     # scales_0.5.0
library(reshape)    # reshape2_1.4.3
require(grid)       
#sessionInfo()
#R version 3.5.1 (2018-07-02)
################################

x = read_biom("./Data/BEF_vipA_454_ITS_1283866816.json_Nov2017.biom") 
fungal_taxonomy<-read.csv2(file="./Data/BEF_vipA_454_ITS_2011_curated_taxonomy_guild.csv", 
                           row.names = 1)

otu_matrix<-as(biom_data(x), "matrix") 
metadata<-sample_metadata(x)

taxonomy_fungi_pre<-as.matrix(fungal_taxonomy)
row.names(taxonomy_fungi_pre)<-rownames(fungal_taxonomy)
taxonomy_fungi<-taxonomy_fungi_pre


BEF_bar_dom<- phyloseq(otu_table(otu_matrix, taxa_are_rows = TRUE), 
                       tax_table(taxonomy_fungi), 
                       sample_data(metadata))

BEF_mono<- subset_samples(BEF_bar_dom, rownames(sample_data(BEF_bar_dom)) %in% 
                            as.character(seq(1,80,by=1)))
BEF1_mono = prune_taxa(taxa_sums(BEF_mono) > 0, BEF_mono)

# High ECM fungal numbers in samples: 11A, 12A, 48A 
#--> exclude these samples for the plot
sample_list<-rownames(sample_data(BEF_bar_dom))
mylist <- sample_list %in% c("11", "12", "48")
newlist <- sample_list[!mylist]

Dataset<-prune_samples(newlist, BEF1_mono)# 4199 taxa, 77 samples

# Split dataset for each tree species, sum up sequence reads
speciesList <- tapply(sample_names(Dataset), 
                      get_variable(Dataset, "sampling_tree"), c) 
speciesPhyseq <- lapply(speciesList, prune_samples, Dataset)
speciesOTUtable <- lapply(speciesPhyseq,otu_table)

speciesAvg <- lapply(speciesOTUtable,rowSums)
pooledOTUtable = t(do.call(rbind,speciesAvg))
pooledOTUtable = data.frame(OTU=row.names(pooledOTUtable),pooledOTUtable)

TT = tax_table(Dataset)
TT = TT[, which(apply(!apply(TT, 2, is.na), 2, any))]
tdf = data.frame(TT, OTU = taxa_names(Dataset))

pOTUtax = merge(pooledOTUtable, tdf, by.x = "OTU")
pOTU = data.frame(pOTUtax,SeqTotal = rowSums(pOTUtax[,2:17])) #pOTUtax[1,2:17]
pOTU = pOTU[order(-pOTU$SeqTotal),]

pOTU.phylum =pOTU[,c(2:17,19)]
pOTU.phylum$Fungal_guild <- factor(pOTU.phylum$Fungal_guild, 
                                   levels = c("Unknown","Saprotroph","Epiphyte",
                                              "Lichenized", "Endophyte", 
                                              "Ericoid Mycorrhizal","Orchid Mycorrhizal",
                                              "Fungal Parasite", "Lichen Parasite",
                                              "Animal Pathogen","Plant Pathogen",
                                              "Ectomycorrhizal","Arbuscular Mycorrhizal" 
                                              ))

melt.phylum = melt(pOTU.phylum,id.vars="Fungal_guild")
colnames(melt.phylum)[2]="species"
agg.phylum=aggregate(.~Fungal_guild+species,melt.phylum,sum)

# Add the myco tree information and check the correct spelling of the trees!
agg.phylum$species<-gsub("Castanea_henryi", "Castanea henryi (EcM)", agg.phylum$species)
agg.phylum$species<-gsub("Quercus_fabri", "Quercus fabri (EcM)", agg.phylum$species)
agg.phylum$species<-gsub("Cyclobalanopsis_myrsinaefolia", "Cyclobalanopsis myrsinifolia (EcM)", agg.phylum$species)
agg.phylum$species<-gsub("Lithocarpus_glaber", "Lithocarpus glaber (EcM)", agg.phylum$species)
agg.phylum$species<-gsub("Castanopsis_sclerophylla", "Castanopsis sclerophylla (EcM)", agg.phylum$species)
agg.phylum$species<-gsub("Cyclobalanopsis_glauca", "Cyclobalanopsis glauca (EcM)", agg.phylum$species)
agg.phylum$species<-gsub("Quercus_serrata", "Quercus serrata (EcM)", agg.phylum$species)
agg.phylum$species<-gsub("Castanopsis_fargesii", "Castanopsis carlesii (EcM)", agg.phylum$species)
agg.phylum$species<-gsub("Castanopsis_eyrei", "Castanopsis eyrei (EcM)", agg.phylum$species)
agg.phylum$species<-gsub("Nyssa_sinensis", "Nyssa sinensis (AM)", agg.phylum$species)
agg.phylum$species<-gsub("Liquidambar_formosana", "Liquidambar formosana (AM)", agg.phylum$species)
agg.phylum$species<-gsub("Koelreuteria_bipinnata","Koelreuteria bipinnata (AM)", agg.phylum$species)
agg.phylum$species<-gsub("Schima_superba","Schima superba (AM)", agg.phylum$species)
agg.phylum$species<-gsub("Sapindus_mukorossi","Sapindus saponaria (AM)", agg.phylum$species)
agg.phylum$species<-gsub("Sapium_sebiferum", "Triadica sebifera (AM)", agg.phylum$species)
agg.phylum$species<-gsub("Choerospondias_axillaris", "Choerospondias axillaris (AM)", agg.phylum$species)
agg.phylum$species<-gsub("Rhus_chinensis", "Rhus chinensis (AM)", agg.phylum$species)

# Sorted by tree taxonomic relationship were possible (genus/family)
positions <- c("Cyclobalanopsis glauca (EcM)","Cyclobalanopsis myrsinifolia (EcM)",
               "Castanopsis sclerophylla (EcM)",
               "Castanopsis carlesii (EcM)",
               "Quercus serrata (EcM)","Quercus fabri (EcM)",
               "Castanea henryi (EcM)","Lithocarpus glaber (EcM)",
               "Choerospondias axillaris (AM)","Rhus chinensis (AM)",
               "Koelreuteria bipinnata (AM)","Sapindus saponaria (AM)",
               "Nyssa sinensis (AM)","Liquidambar formosana (AM)",
               "Triadica sebifera (AM)","Schima superba (AM)")
     

# Adjust fungal guild names for the plotting
agg.phylum$Fungal_guild<-gsub("Epiphyte", "Other fungal guilds", agg.phylum$Fungal_guild)
agg.phylum$Fungal_guild<-gsub("Lichenized", "Other fungal guilds", agg.phylum$Fungal_guild)
agg.phylum$Fungal_guild<-gsub("Endophyte", "Other fungal guilds", agg.phylum$Fungal_guild)
agg.phylum$Fungal_guild<-gsub("Ericoid Mycorrhizal", "Other fungal guilds", agg.phylum$Fungal_guild)
agg.phylum$Fungal_guild<-gsub("Orchid Mycorrhizal", "Other fungal guilds", agg.phylum$Fungal_guild)
agg.phylum$Fungal_guild<-gsub("Fungal Parasite", "Other fungal guilds", agg.phylum$Fungal_guild)
agg.phylum$Fungal_guild<-gsub("Lichen Parasite", "Other fungal guilds", agg.phylum$Fungal_guild)
agg.phylum$Fungal_guild<-gsub("Animal Pathogen", "Other fungal guilds", agg.phylum$Fungal_guild)
agg.phylum$Fungal_guild<-gsub("Arbuscular Mycorrhizal", "Arbuscular mycorrhizal", agg.phylum$Fungal_guild)
agg.phylum$Fungal_guild<-gsub("Plant Pathogen", "Plant pathogenic", agg.phylum$Fungal_guild)
agg.phylum$Fungal_guild<-gsub("Saprotroph", "Saprotrophic", agg.phylum$Fungal_guild)

############## script output

tiff(filename = "./output/images/06_Barplot_FunGuild_monocultures.tiff", 
     width = 5.5, height = 6, units = "in", res=300) 

legend_title <- "Fungal guild"

ggplot(agg.phylum,aes(x=species,y=value,fill=factor(Fungal_guild, 
      levels=c("Unknown","Saprotrophic","Other fungal guilds","Plant pathogenic",
               "Ectomycorrhizal", "Arbuscular mycorrhizal")))) +
  geom_bar(stat="identity",position="fill") +
  ylab("Relative sequence reads\n") +
  xlab("Tree species monocultures")+
  scale_x_discrete(limits = positions)+
  scale_fill_manual(legend_title,
        values = c("#999999","#F0E442","#d3d3d3","#D5000D", "#009E73","#56B4E9")) +
  theme(axis.title.x = element_text(size=12, face="bold"),
        axis.text.x = element_text(angle=70, colour = "black", vjust=1, hjust = 1,
                                   size=12, family="Arial"),
        axis.text.y = element_text(colour = "black", size=12),
        axis.title.y = element_text(size=12, face="bold"),
        plot.title = element_text(size = 14),
         panel.background = element_rect(fill="white"))+
        theme(text=element_text(size=10, family="Arial"))

dev.off()

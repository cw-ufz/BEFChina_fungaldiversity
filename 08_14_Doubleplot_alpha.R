############################### Fungal alpha diversity plots ##################

# R script (08/13) for analyses in Weissbecker et al. 2018
# version: August 2018

#########################

library("extrafont")#extrafont_0.17
library("gdata")    #gdata_2.18.0
library("vegan")    #vegan_2.5-2
library("ecodist")  #ecodist_2.0.1
library(grid)       #
library(gridExtra)  #gridExtra_2.3
#font_import()
#fonts()
#sessionInfo()
#R version 3.5.1 (2018-07-02)
########################

Sapro_rich_residuals<-readRDS(file="./output/07_Sapro_richness_avreps_residuals.RDS")
EcM_rich_residuals<-readRDS(file="./output/07_EcM_richness_avreps_residuals.RDS")
Patho_rich_residuals<-readRDS(file="./output/07_Patho_richness_avreps_residuals.RDS")
AM_rich_residuals<-readRDS(file="./output/07_AM_richness_avreps_residuals.RDS")

Sapro_metadata<-readRDS(file="./output/07_Sapro_alpha_avreps_metadata.RDS")
EcM_metadata<-readRDS(file="./output/07_EcM_alpha_avreps_metadata.RDS")
Patho_metadata<-readRDS(file="./output/07_Patho_alpha_avreps_metadata.RDS")
AM_metadata<-readRDS(file="./output/07_AM_alpha_avreps_metadata.RDS")

plot(AM_rich_residuals ~  AM_metadata$Plot_factor, ylab="Richness residuals",
     xlab="Forest plot", main="AM fungi")

plot(AM_rich_residuals ~  AM_metadata$Env_var$KAK.eff, ylab="Richness residuals",
     xlab="CEC", main="AM fungi")


### Quick lm Anova tests
EcM_lm<-lm(EcM_rich_residuals~EcM_metadata$richness_vars$local_simpson)
summary(EcM_lm) #R= 0.01983,Radj=0.002926
Anova(EcM_lm) #p=0.2832

Sapro_lm<-lm(Sapro_rich_residuals~Sapro_metadata$richness_vars$local_simpson)
summary(Sapro_lm) #R= 0.005734,Radj=-0.007013 
Anova(Sapro_lm) #p=0.5044

Patho_lm<-lm(Patho_rich_residuals~Patho_metadata$richness_vars$local_simpson)
summary(Patho_lm) #R= 0.007718,Radj=-0.007548  
Anova(Patho_lm) #p=0.4796

AM_lm<-lm(AM_rich_residuals~AM_metadata$richness_vars$local_simpson)
summary(AM_lm) #R= 0.0579,Radj=0.04425 
Anova(AM_lm) #p=0.04324

### Prepare plotting

fungal_residuals<-data.frame()
fungal_residuals<-c(EcM_rich_residuals,Sapro_rich_residuals,Patho_rich_residuals,
                    AM_rich_residuals)

AM_metadata$richness_vars$local_simpson<-data.frame(AM_metadata$richness_vars$local_simpson)
AM_metadata$richness_vars$local_simpson[,2]<-rep("AM", 
                    length(AM_metadata$richness_vars$local_simpson[,1]))

Sapro_metadata$richness_vars$local_simpson<-data.frame(Sapro_metadata$richness_vars$local_simpson)
Sapro_metadata$richness_vars$local_simpson[,2]<-rep("Sapro",
                    length(Sapro_metadata$richness_vars$local_simpson[,1]))

Patho_metadata$richness_vars$local_simpson<-data.frame(Patho_metadata$richness_vars$local_simpson)
Patho_metadata$richness_vars$local_simpson[,2]<-rep("Patho",
                    length(Patho_metadata$richness_vars$local_simpson[,1]))

EcM_metadata$richness_vars$local_simpson<-data.frame(EcM_metadata$richness_vars$local_simpson)
EcM_metadata$richness_vars$local_simpson[,2]<-rep("EcM",
                    length(EcM_metadata$richness_vars$local_simpson[,1]))

names(AM_metadata$richness_vars$local_simpson)<-c("Fungal_div", "Fungal_group")
names(EcM_metadata$richness_vars$local_simpson)<-c("Fungal_div", "Fungal_group")
names(Sapro_metadata$richness_vars$local_simpson)<-c("Fungal_div", "Fungal_group")
names(Patho_metadata$richness_vars$local_simpson)<-c("Fungal_div", "Fungal_group")

tree_simpson<-rbind(EcM_metadata$richness_vars$local_simpson, 
                    Sapro_metadata$richness_vars$local_simpson,
                    Patho_metadata$richness_vars$local_simpson,
                    AM_metadata$richness_vars$local_simpson)

fungal_rich_data_compr<-cbind(fungal_residuals, tree_simpson)

############## script output

#data from script 13:
alpha_model_variances<-read.csv(file="./output/13_varpart_alpha_avreps_modified.csv",
                                sep=";", colClasses=c("factor","factor",
                                                      "numeric","factor"))

tiff(filename = "./output/images/14_alpha_tree_variance.tiff", 
     width = 4 , height =6 , units = "in", res=300) 

p1<-ggplot(fungal_rich_data_compr, aes(x=Fungal_div, y=fungal_residuals, 
                                       shape=fungal_rich_data_compr$Fungal_group))+
 geom_point()+
  scale_shape_manual(values = c(0, 15, 8,16))+
  ylab("fungal richness residuals") +
  xlab("Tree Simpson diversity")+
  theme(axis.title.x = element_text(size=12, face="bold"),
        axis.text.x = element_text(angle=70, colour = "black", 
                                   vjust=1, hjust = 1, size=12, family="Arial"),
        axis.text.y = element_text(colour = "black", size=12),
        axis.title.y = element_text(size=12, face="bold"),
        plot.title = element_text(size = 14),
        legend.title = element_blank(),
        legend.text = element_blank(),
        legend.position="none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
  panel.background = element_rect(colour = "black", size=1, fill=NA))+
  theme(text=element_text(size=10, family="Arial"))#


p2<-ggplot(data=alpha_model_variances, aes(x=group, y=perc, fill=order)) +
  geom_bar(stat="identity")+
  coord_flip()+
  ylab(" % variance of richness residuals") +
  xlab("")+
  scale_fill_manual(values = c("#ffffff","#e0e4e4","#a7b4b4","#436060", "#000000"))+  
  theme(axis.title.x = element_text(size=12, face="bold"),
        axis.text.x = element_text(angle=70, colour = "black", vjust=1,
                                   hjust = 1, size=12, family="Arial"),
        axis.text.y = element_text(colour = "black", size=12),
        axis.title.y = element_text(size=12, face="bold"),
        plot.title = element_text(size = 14),
        legend.title = element_blank(),
        legend.text = element_blank(),
        legend.position="none",
        panel.background = element_rect(fill="white"))+
  theme(text=element_text(size=10, family="Arial"))#

grid.arrange(p1, p2, ncol=1, 
             heights=unit(c(4,1.5), c("in", "in")))

dev.off()







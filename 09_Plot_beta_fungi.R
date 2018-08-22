############################### Plot fungal vs. tree beta-diversity ###########

# R script (09/13) for analyses in Weissbecker et al. 2018
# version: August 2018

#########################

library("extrafont")  #extrafont_0.17
library("gdata")      #gdata_2.18.0
library("vegan")      #vegan_2.5-2
library("ecodist")    #ecodist_2.0.1
library("ggplot2")    #ggplot2_3.0.0
library("Rmisc")      #Rmisc_1.5
  #font_import()
  #fonts()
#sessionInfo()
#R version 3.5.1 (2018-07-02)
########################

source(file="./R_functions/09_Plot_beta_fungi_localtree_function.R")

Sapro_otu<-readRDS(file="./output/03_Sapro_otu_hell_avreps.RDS")
EcM_otu<-readRDS(file="./output/03_EcM_otu_hell_avreps.RDS")
Patho_otu<-readRDS(file="./output/03_Patho_otu_hell_avreps.RDS")
AM_otu<-readRDS(file="./output/03_AM_otu_hell_avreps.RDS")

Sapro_metadata<-readRDS(file="./output/03_Sapro_metadata_trans_avreps.RDS")
EcM_metadata<-readRDS(file="./output/03_EcM_metadata_trans_avreps.RDS")
Patho_metadata<-readRDS(file="./output/03_Patho_metadata_trans_avreps.RDS")
AM_metadata<-readRDS(file="./output/03_AM_metadata_trans_avreps.RDS")

Sapro_tree<-Sapro_metadata$local_tree 
EcM_tree<-EcM_metadata$local_tree 
AM_tree<-AM_metadata$local_tree
Patho_tree<-Patho_metadata$local_tree

# Create beta diversity information
Tree_beta_Sapro<-create_beta(Sapro_tree, "hellinger" , "bray")
Sapro_beta<-create_beta(Sapro_otu, "none" ,"bray")

Tree_beta_EcM<-create_beta(EcM_tree, "hellinger" , "bray")
EcM_beta<-create_beta(EcM_otu, "none" ,"bray")

Tree_beta_AM<-create_beta(AM_tree, "hellinger" , "bray")
AM_beta<-create_beta(AM_otu, "none" ,"bray")

Tree_beta_Patho<-create_beta(Patho_tree, "hellinger" , "bray")
Patho_beta<-create_beta(Patho_otu, "none" ,"bray")


Sapro_Plot_beta<-create_beta(Sapro_metadata$plot, "none" ,"bray")
Patho_Plot_beta<-create_beta(Patho_metadata$plot, "none" ,"bray")
EcM_Plot_beta<-create_beta(EcM_metadata$plot, "none" ,"bray")
AM_Plot_beta<-create_beta(AM_metadata$plot, "none" ,"bray")

Sapro_beta_vector<-as.vector(Sapro_beta)
Sapro_beta_gg<-data.frame(cbind(Sapro_beta,Tree_beta_Sapro, Sapro_Plot_beta ))
color_points_Sapro<-Sapro_beta_gg[,3]
color_points_Sapro[color_points_Sapro=="1"]<-"black"
color_points_Sapro[color_points_Sapro=="0"]<-"blue"

Patho_beta_vector<-as.vector(Patho_beta)
Patho_beta_gg<-data.frame(cbind(Patho_beta,Tree_beta_Patho, Patho_Plot_beta ))
color_points_Patho<-Patho_beta_gg[,3]
color_points_Patho[color_points_Patho=="1"]<-"black"
color_points_Patho[color_points_Patho=="0"]<-"blue"

EcM_beta_vector<-as.vector(EcM_beta)
EcM_beta_gg<-data.frame(cbind(EcM_beta,Tree_beta_EcM, EcM_Plot_beta ))
color_points_EcM<-EcM_beta_gg[,3]
color_points_EcM[color_points_EcM=="1"]<-"black"
color_points_EcM[color_points_EcM=="0"]<-"blue"

AM_beta_vector<-as.vector(AM_beta)
AM_beta_gg<-data.frame(cbind(AM_beta,Tree_beta_AM, AM_Plot_beta ))
color_points_AM<-AM_beta_gg[,3]
color_points_AM[color_points_AM=="1"]<-"black"
color_points_AM[color_points_AM=="0"]<-"blue"

############## script output

p1<-ggplot(Sapro_beta_gg[,], aes(x=Sapro_beta_gg[,2], y=Sapro_beta_gg[,1]))+
  geom_point(shape=16, color = color_points_Sapro, alpha = 1, size=0.5)+
  ylab("Saprotroph beta-diversity")+
  xlab("Tree beta-diversity")+
  theme_bw() +
  scale_y_continuous(limits = c(0.0, 1.0))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

p2<-ggplot(Patho_beta_gg[,], aes(x=Patho_beta_gg[,2], y=Patho_beta_gg[,1]))+
      geom_point(shape=16, color = color_points_Patho, alpha = 1, size=0.5)+
  ylab("Pathogen beta-diversity")+
  xlab("Tree beta-diversity")+
  theme_bw() +
  scale_y_continuous(limits = c(0.0, 1.0))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())
  
p3<-ggplot(AM_beta_gg[,], aes(x=AM_beta_gg[,2], y=AM_beta_gg[,1]))+
  geom_point(shape=16, color = color_points_AM, alpha = 1, size=0.5)+
  ylab("AM fungal beta-diversity")+
  xlab("Tree beta-diversity")+
  theme_bw() +
  scale_y_continuous(limits = c(0.0, 1.0))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

p4<-ggplot(EcM_beta_gg[,], aes(x=EcM_beta_gg[,2], y=EcM_beta_gg[,1]))+
  geom_point(shape=16, color = color_points_EcM, alpha = 1, size=0.5)+
  ylab("EcM fungal beta-diversity")+
  xlab("Tree beta-diversity")+
  theme_bw() +
  scale_y_continuous(limits = c(0.0, 1.0))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


tiff(filename ="./output/images/08_beta_fungi_localtree.tiff",
     width = 6.8, height = 6.8, units = "in", res=300)

multiplot(p1,p2,p3,p4, cols=2)

dev.off()



 








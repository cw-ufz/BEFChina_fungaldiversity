############################### Calculate beta diversity for datasets #########

# R function to script 
# R script (09/13) for analyses in Weissbecker et al. 2018
# version: August 2018

######################
library(vegan)      #vegan_2.5-2
library(gdata)      #gdata_2.18.0
#sessionInfo()
#R version 3.5.1 (2018-07-02)
######################

# The function needs three arguments:
# 1) data, if phyloseq it will extract the OTU table 
# 2) type of transformation ("none" or decostand options)
# 3) type of distance measure ("none" or vergdist options)
# Example: Fungi_beta<-create_beta(fungi_phyloseq,"transform-method", "distance-measure")
#################################

create_beta <- function(mydata, trans, distance_measure) {  
  if (class(mydata)[1]=="phyloseq") {
    print("data recognized as phyloseq object"); mydata<-t(otu_table(mydata))
  } 
  
  if (trans != "none"){
    if(trans == "log(x+1)"){
      mydata<-log(mydata+1)
    } else {
      mydata<-decostand(mydata, method=trans) 
    }
    print("data transformation")
    print(trans)
  }
  
  if (distance_measure != "none"){
    print("dissimilarity calculation")
    print(distance_measure)
    
    mydata<-as.matrix(vegdist(mydata, method=distance_measure))
  }
  
  mydata_dist<-lowerTriangle(mydata, diag=FALSE) 
  
  return(mydata_dist)
}
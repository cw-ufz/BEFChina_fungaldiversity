########################## Variance inflation factor stepwise calculation #####

# R function for
# R script (07/13) for analyses in Weissbecker et al. 2018
# version: August 2018

#author was an anonymous author in one of the r forums:
#https://beckmw.wordpress.com/2013/02/05/collinearity-and-stepwise-vif-selection/

###############################

library(fmsb)  #fmsb_0.6.3
#sessionInfo()
#R version 3.5.1 (2018-07-02)
###############################

vif_func<-function(in_frame,thresh=10,trace=T,...){
  
  require(fmsb)
  if(class(in_frame) != 'data.frame') in_frame<-data.frame(in_frame)
  #get initial vif value for all comparisons of variables
  vif_init<-NULL
  for(val in names(in_frame)){
    form_in<-formula(paste(val,' ~ .'))
    vif_init<-rbind(vif_init,c(val,VIF(lm(form_in,data=in_frame,...))))
  }
  vif_max<-max(as.numeric(vif_init[,2]))
  
  if(vif_max < thresh){
    if(trace==T){ #print output of each iteration
      prmatrix(vif_init,collab=c('var','vif'),rowlab=rep('',nrow(vif_init)),quote=F)
      cat('\n')
      cat(paste('All variables have VIF < ', thresh,', max VIF ',round(vif_max,2), sep=''),'\n\n')
    }
    return(names(in_frame))
  }
  else{
    
    in_dat<-in_frame
    
    while(vif_max >= thresh){
      vif_vals<-NULL
      
      for(val in names(in_dat)){
        form_in<-formula(paste(val,' ~ .'))
        vif_add<-VIF(lm(form_in,data=in_dat,...))
        vif_vals<-rbind(vif_vals,c(val,vif_add))
      }
      max_row<-which(vif_vals[,2] == max(as.numeric(vif_vals[,2])))[1]
      
      vif_max<-as.numeric(vif_vals[max_row,2])
      
      if(vif_max<thresh) break
      if(trace==T){ #print output of each iteration
        prmatrix(vif_vals,collab=c('var','vif'),rowlab=rep('',nrow(vif_vals)),quote=F)
        cat('\n')
        cat('removed: ',vif_vals[max_row,1],vif_max,'\n\n')
        flush.console()
      }
      
      in_dat<-in_dat[,!names(in_dat) %in% vif_vals[max_row,1]]
      
    }
    
    return(names(in_dat))
  }
} 
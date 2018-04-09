## Calculate all possible combinations of variables for best subset approach ##

# R function to script (11/13) for analyses in Weissbecker et al. 2017
# version: April 2018

# The function needs a list with variables as argument
# Example: vector<-list("A", "B", "C"); mrm.subset.list(vector)
#################################

mrm.subset.list<-function(list_of_parameters){
  library(combinat)
  comb_par<-list()
  for (i in seq(from=1, to = length(list_of_parameters), by=1)){
    #print(i)
    comb_par[[length(comb_par)+1]]<-combn(list_of_parameters, i)
  }
  
return(comb_par)
}




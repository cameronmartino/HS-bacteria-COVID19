library(openxlsx)
library(dplyr)
library(tidyverse)

if(!dir.exists('code/glycogene-clusters')){
  stop('glycogene currations are availible upon request')
}

tmp = do.call(rbind,lapply(list.files('data/glycogene-clusters/CAZy/',pattern = 'CAZy_*',full.names = T),read.csv,header=F))
colnames(tmp) = c('protein','CAZy.ID','x1','x2','protein_name','protein_id','x3','x4','status','x5','x6')
write.csv(tmp,file = 'CAZy_all.clean.csv')
ec = read.xlsx('data/glycogene-clusters/CAZy/cazyEC.xlsx',sheet = 1)
ec$CAZy.Type[ec$CAZy.Type=="Glycoside-Hydrolases"] = "GH"
ec$CAZy.Type[ec$CAZy.Type=="Carbohydrate-Esterases"] = "CE"
ec$CAZy.Type[ec$CAZy.Type=="Polysaccharide-Lyases"] = "PL"
ec$CAZy.Type[ec$CAZy.Type=="Auxiliary-Activities"] = "AA"
ec$CAZy.Type[ec$CAZy.Type=="GlycosylTransferases"] = "GT"
ec$CAZy.Type[ec$CAZy.Type=="Carbohydrate-Binding-Modules"] = "CBM"
ec = separate_rows(ec,CAZy.family,sep='\\s')
ec$CAZy.ID = paste0(ec$CAZy.Type,ec$CAZy.family)
write.csv(ec,file='EC2CAZy.csv',row.names = F)
library(vegan)
library(tidyverse)
library(FactoMineR)
library(factoextra)

load('C:/Users/33638/Documents/stage/MICROBIOTA_SPATIAL_FLIGHT/DATA_PROJECT_1.RData')


### Certains OTU ne sont observés pour aucun individu, on les supprime
nb <- colSums(matrix_otu)
matrix_otu <- matrix_otu[,nb!=0]  #10 OTU supprimés

OTU_TSS<-NORMALISATION_MICROBIOTA(matrix_otu, approches = 'TSS')
OTU_CSS<-NORMALISATION_MICROBIOTA(matrix_otu, approches = 'CSS')

res_PCA_tot <- PCA(OTU_TSS, scale.unit = T)
# le premier plan factoriel contient 21% des informations 
OTU_PCA$eig
OTU_PCA$var$coord
OTU_PCA$ind$coord

## Graphe des individus, coloré selon shannon
shannon_tot <- diversity(matrix_otu, index='shannon')
fviz_pca_ind(res_PCA_tot,col.ind = shannon_tot,gradient.cols = c("yellow","red"),pointsize = 2,repel=T)



###############################
#### REGROUPEMENT  ###########
##############################
# on regroupe les OTU pour avoir moins de variables
t_matrix_otu_species <- merge(t(OTU_TSS),species,by="row.names")
t_matrix_otu_species <- tibble::rownames_to_column(t_matrix_otu_species,"OTU")
t_matrix_otu_species <- tibble(t_matrix_otu_species)
summary(t_matrix_otu_species)


### Phylum
res_phylum=levels(t_matrix_otu_species$Phylum)
for(col in 3:30){
  res=aggregate(t_matrix_otu_species[,col],t_matrix_otu_species[,"Phylum"],sum)
  res_phylum=cbind(res_phylum,res[,2])
}
row.names(res_phylum)<-res_phylum[,1]
res_phylum <- apply(as.data.frame(res_phylum[,2:29]),1,FUN=as.double)

#PCA
res_PCA_phylum <- PCA(res_phylum)


#on veut rajouter le shannon sur graphe individu
t_matrix_otu_species_brutes <- merge(t(matrix_otu),species,by="row.names")
t_matrix_otu_species_brutes <- tibble::rownames_to_column(t_matrix_otu_species_brutes,"OTU")
t_matrix_otu_species_brutes <- tibble(t_matrix_otu_species_brutes)

res_Phylum_brut=levels(t_matrix_otu_species_brutes$Phylum)
for(col in 3:30){
  res=aggregate(t_matrix_otu_species_brutes[,col],t_matrix_otu_species_brutes[,"Phylum"],sum)
  res_Phylum_brut=cbind(res_Phylum_brut,res[,2])
}
row.names(res_Phylum_brut)<-res_Phylum_brut[,1]
res_Phylum_brut <- apply(as.data.frame(res_Phylum_brut[,2:29]),1,FUN=as.integer)

shannon_tot_phylum <- diversity(res_Phylum_brut,index='shannon')

#le graphe des individus
fviz_pca_ind(res_PCA_phylum,col.ind = shannon_tot_phylum,gradient.cols = c("yellow","red"),pointsize = 2,repel=T)







library(vegan)
library(tidyverse)
library(FactoMineR)
library(factoextra)


load('C:/Users/33638/Documents/stage/MICROBIOTA_SPATIAL_FLIGHT/DATA_PROJECT_1.RData')


res_PCA_tot <- PCA(matrix_otu, scale.unit = T)
# le premier plan factoriel contient 21% des informations 
OTU_PCA$eig
OTU_PCA$var$coord
OTU_PCA$ind$coord

## Graphe des individus, coloré selon shannon
shannon_tot <- diversity(matrix_otu, index='shannon')
fviz_pca_ind(res_PCA_tot,col.ind = shannon_tot,label="None")



###############################
#### REGROUPEMENT  ###########
##############################
# on regroupe les OTU pour avoir moins de variables
t_matrix_otu_species <- cbind(t(matrix_otu),species)
t_matrix_otu_species <- tibble::rownames_to_column(t_matrix_otu_species,"OTU")
t_matrix_otu_species <- tibble(t_matrix_otu_species)
summary(t_matrix_otu_species)


### Phylum
res_phylum=levels(t_matrix_otu_species$Phylum)
for(col in 2:29){
  res=aggregate(t_matrix_otu_species[,col],t_matrix_otu_species[,"Phylum"],sum)
  res_phylum=cbind(res_phylum,res[,2])
}
row.names(res_phylum)<-res_phylum[,1]
res_phylum <- apply(as.data.frame(res_phylum[,2:29]),1,FUN=as.integer)

shannon_tot_phylum <- diversity(res_phylum,index='shannon')
#PCA
res_PCA_phylum <- PCA(res_phylum)
fviz_pca_ind(OTU_PCA,col.ind = shannon_tot_phylum,label="None")


### Class
res_class=levels(t_matrix_otu_species$Class)
for(col in 2:29){
  res=aggregate(t_matrix_otu_species[,col],t_matrix_otu_species[,"Class"],sum)
  res_class=cbind(res_class,res[,2])
}
row.names(res_class)<-res_class[,1]
res_class <- apply(as.data.frame(res_class[,2:29]),1,FUN=as.integer)

shannon_tot_class <- diversity(res_class,index='shannon')
#PCA
res_PCA_class <- PCA(res_class)
fviz_pca_ind(res_PCA_class,col.ind = shannon_tot_class,label="None")



### Order
res_order=levels(t_matrix_otu_species$Order)
for(col in 2:29){
  res=aggregate(t_matrix_otu_species[,col],t_matrix_otu_species[,"Order"],sum)
  res_order=cbind(res_order,res[,2])
}
row.names(res_order)<-res_order[,1]
res_order <- apply(as.data.frame(res_order[,2:29]),1,FUN=as.integer)

shannon_tot_order <- diversity(res_order,index='shannon')
#PCA
res_PCA_order <- PCA(res_order)
fviz_pca_ind(res_PCA_order,col.ind = shannon_tot_order,label="None")





### Family
res_family=levels(t_matrix_otu_species$Family)
for(col in 2:29){
  res=aggregate(t_matrix_otu_species[,col],t_matrix_otu_species[,"Family"],sum)
  res_family=cbind(res_family,res[,2])
}

row.names(res_family)<-res_family[,1]
res_family <- apply(as.data.frame(res_family[,2:29]),1,FUN=as.integer)

shannon_tot_family <- diversity(res_family,index='shannon')
#PCA
res_PCA_family <- PCA(res_family)
fviz_pca_ind(res_PCA_family,col.ind = shannon_tot_family,label="None")




### Genus
res_genus=levels(t_matrix_otu_species$Genus)
for(col in 2:29){
  res=aggregate(t_matrix_otu_species[,col],t_matrix_otu_species[,"Genus"],sum)
  res_genus=cbind(res_genus,res[,2])
}
row.names(res_genus)<-res_genus[,1]
res_genus <- apply(as.data.frame(res_genus[,2:29]),1,FUN=as.integer)

shannon_tot_genus <- diversity(res_genus,index='shannon')
#PCA
res_PCA_genus <- PCA(res_genus)
fviz_pca_ind(res_PCA_genus,col.ind = shannon_tot_genus,label="None")






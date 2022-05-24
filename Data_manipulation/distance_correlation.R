# Code pour:
# - matrices distances Pearson et Jaccard
# - Dendrogramme => 3 groupes liés à co-occurence OTU
# - Diversité PD sur treefile et dendrogramme
# - Essai Package Magma et distance Chi2
# - PCA globale et PCA en regroupant par phylum
# - Heatmap sur les indices alpha diversity

library(devtools)
install_gitlab("arcgl/rmagma")
library(rMAGMA)
library(tidyverse)
library(corrplot)
library(picante)
library(phylogram)
library(FactoMineR)
library(factoextra)

load('C:/Users/33638/Documents/stage/MICROBIOTA_SPATIAL_FLIGHT/DATA_PROJECT_1.RData')

### Certains OTU ne sont observés pour aucun individu, on les supprime
nb <- colSums(matrix_otu)
matrix_otu <- matrix_otu[,nb!=0]  #10 OTU supprimés


###########################
######## Pearson ########## 
###########################

# normalisation 
OTU_TSS<-NORMALISATION_MICROBIOTA(matrix_otu, approches = 'TSS')
OTU_CSS<-NORMALISATION_MICROBIOTA(matrix_otu, approches = 'CSS')

## matrice de corrélation de Pearson
mcor<- cor(OTU_TSS, method="pearson")

# on la transforme en matrice de distance
dist_pearson <- as.dist(mcor)

# on fait une CAH, avec stratégie d'aggrégation de Ward
cah.pearson <- hclust(dist_pearson,method="ward.D")

#on représente le dendrogramme associé
plot(cah.pearson,labels = F)


##################################### 
######## Présence/abscence ########## 
##################################### 

#on affecte 0 et 1 en focntion de si abscence ou présence de l'otu
matrix_otu_01 <- matrix_otu>0

#On calcule les distances

#distance de Jaccard :
Jaccard = function(i,j){
  return(sum(i!=j)/(sum(i!=j)+sum((i==1 & j==1))))
}
#distance de Dice (double poids accords)
Dice = function(i,j){
  return(sum(i!=j)/(sum(i!=j)+2*sum((i==1 & j==1))))
}

# On calcule la matrice de dissimilarité à l'aide de la fonction `outer` 
# pour calculer la valeur de dissimilarité pour toutes les paires d'OTU possibles.

MatDiss = outer(as.data.frame(matrix_otu_01),as.data.frame(matrix_otu_01),Vectorize(Jaccard))


# on la transforme en matrice de distance
dist_01 <- as.dist(MatDiss)

# on fait une CAH, avec stratégie d'aggrégation de Ward
cah.01 <- hclust(dist_01,method="ward.D")

#on représente le dendrogramme associé
plot(cah.01,labels = F)

# Regardons la courbe de perte d'inertie 
# (on se contente des 20 premières valeurs pour ne pas "noyer" l'information importante)
plot(rev(cah.01$height)[1:20],type="b",main="inertie intra")
   # on peut conserver 3 groupes
plot(cah.01,labels = F)
rect.hclust(cah.01, 3, border ="blue")

## on colore selon phylum
library(dendextend)
dend <- as.dendrogram(cah.01)

# Par défaut les labels n'ont pas de couleurs
labels_colors(dend)

colors_to_use <- as.numeric(species$Phylum)
colors_to_use <- colors_to_use[order.dendrogram(dend)]

# Now we can use them
labels_colors(dend) <- colors_to_use
# Now each state has a color
labels_colors(dend) 

plot(dend, main = "A color for every Phylum")
rect.hclust(cah.01, 3, border ="blue")


##### on récupère les groupes
K=3
gpe = cutree(cah.01,k=K)

matrix_otu_gpe <- rbind(matrix_otu,gpe)

t_matrix_otu_species_brutes <- merge(t(matrix_otu_gpe),species,by="row.names")
t_matrix_otu_species_brutes <- tibble::rownames_to_column(t_matrix_otu_species_brutes,"OTU")

# on va regarder la proportion pour chaque groupe des catégories de Phylum
group1 <- t_matrix_otu_species_brutes[t_matrix_otu_species_brutes[,"gpe"]==1,]
group2 <- t_matrix_otu_species_brutes[t_matrix_otu_species_brutes[,"gpe"]==2,]
group3 <- t_matrix_otu_species_brutes[t_matrix_otu_species_brutes[,"gpe"]==3,]

group1_phylum=levels(group1$Phylum)
for(col in 3:30){
  res=aggregate(group1[,col],list(group1[,"Phylum"]),sum)
  group1_phylum=cbind(group1_phylum,res[,2])
}
group1_phylum <- cbind(group1_phylum[,1],t(as.data.frame(apply(as.data.frame(group1_phylum[,2:29]),1,FUN=as.integer))),rep("group1",4))

group2_phylum=levels(group2$Phylum)
for(col in 3:30){
  res=aggregate(group2[,col],list(group2[,"Phylum"]),sum)
  group2_phylum=cbind(group2_phylum,res[,2])
}
group2_phylum <- cbind(group2_phylum[,1],t(as.data.frame(apply(as.data.frame(group2_phylum[,2:29]),1,FUN=as.integer))),rep("group2",4))

group3_phylum=levels(group3$Phylum)
for(col in 3:30){
  res=aggregate(group3[,col],list(group3[,"Phylum"]),sum)
  group3_phylum=cbind(group3_phylum,res[,2])
}
group3_phylum <- cbind(group3_phylum[,1],t(as.data.frame(apply(as.data.frame(group3_phylum[,2:29]),1,FUN=as.integer))),rep("group3",4))

## on concatene les infos sur les 3 groupes
analyse_phylum_gpe <- rbind(group1_phylum,group2_phylum,group3_phylum) 

#on fait la moyenne de chaque ligne
moy <- apply(analyse_phylum_gpe[,2:29],1,function(x) mean(as.integer(x)))
analyse_phylum_gpe <- data.frame(analyse_phylum_gpe,moy)
analyse_phylum_gpe <- data.frame(analyse_phylum_gpe[,c(1,30)],t(apply(analyse_phylum_gpe[,c(2:29,31)],1,FUN=as.integer)))

#X29 la moyenne (pas trop de sens ?)
ggplot(analyse_phylum_gpe) + aes(x=X30, y=X29, fill=X1)+
  geom_bar(stat="identity", position=position_dodge()) +
  labs(title="Répartition des bactéries dans les groupes \n en moyenne pour les 28 indiv", 
       x="Groupes issus CAH", y = "Nombre d'OTU")


#pour un individu (changer y entre X1.1 et X28)
ggplot(analyse_phylum_gpe) + aes(x=X30, y=X3, fill=X1)+
  geom_bar(stat="identity", position=position_dodge())+
  labs(title="Répartition des bactéries dans les groupes \n pour un individu", 
       x="Groupes issus CAH", y = "Nombre d'OTU")



###########################
######## Diversité PD######
###########################
plot(treefile)
plot.phylo(treefile)
plotTree(treefile,type="fan",fsize=0.7,lwd=1,
         ftype="i")


# On deux arbres :
 # - l'arbre phylogénétique treefile : représente la parenté entre OTU
 # - le dendrogramme qui peut etre vu comme un arbre représentant la co-occurence

## PD de Faith
pd_treefile <- pd(matrix_otu, treefile)

phylo_dend <- as.phylo(dend)
pd_dend <- pd(matrix_otu, phylo_dend)

plot(pd_treefile$PD,pd_dend$PD)


#en prenant en compte l'abondance
# Calculates mean pairwise distance separating taxa in a community
mpd <- mpd(matrix_otu, cophenetic(treefile), abundance.weighted=TRUE)
     #cophenetic : distance entre deux objets calculés à partir d'un dendrogramme
mpd_dend <- mpd(matrix_otu, cophenetic(phylo_dend), abundance.weighted=TRUE)

plot(mpd,mpd_dend)





    


#############################################
## matrice des distances du chi-2
library(ade4)

# l'AFC se base sur les distance du X² (on peut considérer que c'est une ACP avec la métrique du X²)
# on peut donc récupérer la matrice des distances du X² avec les fonctions AFC du package ade4
afc <- dudi.coa(matrix_otu,scan=F)
# calcul de distance
dist_chi2 <- as.matrix(dist.dudi(afc))







### Magma #############
#otu_table_0 <- matrix_otu
#prevalence  <- colMeans(otu_table_0>0)
#sequencing_depth <- rowSums(otu_table_0)

#icol <- prevalence >0.25
#irow <- sequencing_depth >500

#otu_table <- otu_table_0[irow,icol]

#magma_Stool <- magma(data = otu_table)
#plot(magma_Stool)




#######################################
######  PCA ###########################
#######################################

res_PCA_tot <- PCA(OTU_TSS, scale.unit = T)
# le premier plan factoriel contient 21% des informations 
OTU_PCA$eig
OTU_PCA$var$coord
OTU_PCA$ind$coord

## Graphe des individus, coloré selon shannon
shannon_tot <- vegan::diversity(matrix_otu, index='shannon')
fviz_pca_ind(res_PCA_tot,col.ind = shannon_tot,gradient.cols = c("yellow","red"),pointsize = 2,repel=T)



#### REGROUPEMENT  ############
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

shannon_tot_phylum <- vegan::diversity(res_Phylum_brut,index='shannon')

#le graphe des individus
fviz_pca_ind(res_PCA_phylum,col.ind = shannon_tot_phylum,gradient.cols = c("yellow","red"),pointsize = 2,repel=T)



###########################
###### heatmap  ##########
##########################
richness <-apply(matrix_otu, 1, function(x) length(which(x>=1)))
#La richesse varie entre 209 et 316 selon les individus, avec une moyenne de 273 OTU différents détectés
range(richness)

#### Chao1
chao1<-apply(matrix_otu, 1, function(x) chao1(x, taxa.row = T))

simpson <- vegan::diversity(matrix_otu, index='simpson')
shannon <-vegan::diversity(matrix_otu, index='shannon')

indice_alpha_div <- as.matrix(data.frame(richness,chao1,simpson,shannon,pd_treefile$PD,pd_dend$PD,mpd,mpd_dend))

heatmap(indice_alpha_div, scale = "col")







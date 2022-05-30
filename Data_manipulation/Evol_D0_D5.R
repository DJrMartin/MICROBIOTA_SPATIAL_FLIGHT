# 1) RICHNESS
# 2) SIMPSON
# 3) SHANNON
# 4) PD
# 5) PDW
# 6) RANK ABUNDANCE CURVE
# 7) CO-OCCURENCE DIVERSITY

library(phangorn)
library(vegan)  #diversity
library(tidyverse)
library(fossil)
library(tree)
library(ape)  
library(phytools)
library(geiger)
library(picante)
library(stats)
library(RColorBrewer)

#load('C:/Users/33638/Documents/stage/MICROBIOTA_SPATIAL_FLIGHT/DATA_PROJECT_1.RData')

### Certains OTU ne sont observés pour aucun individu, on les supprime
nb <- colSums(matrix_otu)
matrix_otu <- matrix_otu[,nb!=0]  #10 OTU supprimés
species <- species[nb!=0,]


###########################################
### On crée différentes tables ############
###########################################

#matrice avec données brutes et espèces
matrix_otu_species <- merge(t(matrix_otu),species,by="row.names")
matrix_otu_species <- tibble::rownames_to_column(matrix_otu_species,"OTU")


# pour différencier D0 et D5 ---------------------------
dates <- rep(c("D0","D5"),14)
matrix_otu_D <- data.frame(dates,matrix_otu)
matrix_otu_D0 <- matrix_otu_D[matrix_otu_D$dates=="D0",2:ncol(matrix_otu_D)]
matrix_otu_D5 <- matrix_otu_D[matrix_otu_D$dates=="D5",2:ncol(matrix_otu_D)]

#matrice avec données brutes D0 et espèces
matrix_otu_species_D0 <- merge(t(matrix_otu_D0),species,by="row.names")
matrix_otu_species_D0 <- tibble::rownames_to_column(matrix_otu_species_D0,"OTU")

#matrice avec données brutes D5 et espèces
matrix_otu_species_D5 <- merge(t(matrix_otu_D5),species,by="row.names")
matrix_otu_species_D5 <- tibble::rownames_to_column(matrix_otu_species_D5,"OTU")

#matrice avec données nomalisées D0 et espèces
matrix_otu_norm_D0 <-NORMALISATION_MICROBIOTA(matrix_otu_D0, approches = 'TSS')
matrix_otu_norm_species_D0 <- merge(t(matrix_otu_norm_D0),species,by="row.names")
matrix_otu_norm_species_D0 <- tibble::rownames_to_column(matrix_otu_norm_species_D0,"OTU")

#matrice avec données normalisées D5 et espèces
matrix_otu_norm_D5 <-NORMALISATION_MICROBIOTA(matrix_otu_D5, approches = 'TSS')
matrix_otu_norm_species_D5 <- merge(t(matrix_otu_norm_D5),species,by="row.names")
matrix_otu_norm_species_D5 <- tibble::rownames_to_column(matrix_otu_norm_species_D5,"OTU")

###########################################
###    Alpha diversité        ############
###########################################

##### 1) RICHESSE
richness_D0 <-apply(matrix_otu_D0, 1, function(x) length(which(x>=1)))
richness_D5 <-apply(matrix_otu_D5, 1, function(x) length(which(x>=1)))
plot(x=1:14,y=richness_D0, main="Evolution de la richesse entre D0 et D5 pour les 14 individus",ylab="Richesse")
points(x=1:14,y=richness_D5,col="red")
legend(legend=c("D0","D5"), 
       fill=c("black","red"), 'topleft',lty=0)
   #Pas de tendance claire vers une augmentation ou diminution de la richesse entre D0 et D5

##### 2) SIMPSON
simpson_D0 <-vegan::diversity(matrix_otu_D0, index='simpson')
simpson_D5 <-vegan::diversity(matrix_otu_D5, index='simpson')
plot(x=1:14,y=simpson_D0, main="Evolution de simpson entre D0 et D5 pour les 14 individus",ylab="Simpson")
points(x=1:14,y=simpson_D5,col="red")
legend(legend=c("D0","D5"), 
       fill=c("black","red"), 'topleft',lty=0)
#très variable d'un individu à l'autre

##### 3) SHANNON
shannon_D0 <-vegan::diversity(matrix_otu_D0, index='shannon')
shannon_D5 <-vegan::diversity(matrix_otu_D5, index='shannon')
plot(x=1:14,y=shannon_D0, main="Evolution de shannon entre D0 et D5 pour les 14 individus",ylab="Shannon")
points(x=1:14,y=shannon_D5,col="red")
legend(legend=c("D0","D5"), 
       fill=c("black","red"), 'topleft',lty=0)

plot(shannon_D0,shannon_D5)

#### 4) PD
plot(treefile, "f", show.tip.label = FALSE, no.margin = TRUE)
pd_treefile_D0 <- pd(matrix_otu_D0, treefile)
pd_treefile_D5 <- pd(matrix_otu_D5, treefile)
plot(x=1:14,y=pd_treefile_D0$PD, main="Evolution de phylogenetic diversity (PD) \n entre D0 et D5 pour les 14 individus",ylab="PD",ylim=c(12,20))
points(x=1:14,y=pd_treefile_D5$PD,col="red")
legend(legend=c("D0","D5"), 
       fill=c("black","red"), 'topleft',lty=0)

plot(pd_treefile_D0$PD,pd_treefile_D5$PD)
boxplot(data.frame(D0=pd_treefile_D0$PD,D5=pd_treefile_D5$PD))

#### 5) PD Weighted
# Calculates mean pairwise distance separating taxa in a community
mpd_D0 <- mpd(matrix_otu_D0, cophenetic(treefile), abundance.weighted=TRUE)
mpd_D5 <- mpd(matrix_otu_D5, cophenetic(treefile), abundance.weighted=TRUE)
plot(x=1:14,y=mpd_D0, main="Evolution de phylogenetic diversity weighted \n entre D0 et D5 pour les 14 individus",ylab="MPD",ylim=c(0.5,1))
points(x=1:14,y=mpd_D5,col="red")
legend(legend=c("D0","D5"), 
       fill=c("black","red"), 'topleft',lty=0)
#boxplot
boxplot(data.frame(D0=mpd_D0,D5=mpd_D5),main="Boxplot de Weighted PD")
#p value ns.

#### 6) RANK ABUNDANCE CURVE
max_otu <- apply(matrix_otu,1,max)
max_otu_matrice <- matrix(rep(max_otu,617),nrow = 28)
matrix_otu_abund <- matrix_otu/max_otu_matrice

# comparaison D0 D5
indiv1 <- sort(matrix_otu_abund[1,matrix_otu_abund[1,]!=0],decreasing = T)
indiv2 <- sort(matrix_otu_abund[2,matrix_otu_abund[2,]!=0],decreasing = T)
plot(indiv1,type="l",main="Rank abundance curve à D0 et D5 pour indiv A",ylab="Proportion relative OTU")
points(sort(indiv2,decreasing = T),type="l",col="red")
legend(legend=c("D0","D5"), 
       fill=c("black","red"), 'topright',lty=0)

# Comparaisons 3 individus
plot(log(sort(matrix_otu_abund[1,],decreasing = T)),type="l",main="Rank abundance curve pour indiv A,B et C (à D0)",ylab="Proportion relative OTU",xlim=c(0,200))
points(log(sort(matrix_otu_abund[3,],decreasing = T)),type="l",col="red")
points(log(sort(matrix_otu_abund[5,],decreasing = T)),type="l",col="green")
legend(legend=c("A","B","C"), 
       fill=c("black","red","green"), 'topright',lty=0)

# On passe au log
plot(log(sort(matrix_otu_abund[1,matrix_otu_abund[1,]!=0],decreasing = T)),type="l")

# nouvel indice de diversité : coeff pente estimée du log de la courbe
pentes=c()
for(i in 1:28){
  y=log(sort(matrix_otu_abund[i,matrix_otu_abund[i,]!=0],decreasing = T))
  x=1:length(y)
  res <- lm(y~x)
  pentes <- c(pentes,res$coefficients[2])
}
plot(pentes)  #indice pour les 28 individus





#### 7) DIVERSITY OF CO-OCCURENCE
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

# Regardons la courbe de perte d'inertie 
# (on se contente des 20 premières valeurs pour ne pas "noyer" l'information importante)
plot(rev(cah.01$height)[1:20],type="b",main="inertie intra")

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

## UNWEIGHTED
phylo_dend <- as.phylo(dend)
plot(phylo_dend, "f", show.tip.label = FALSE, no.margin = TRUE)
pd_dend <- pd(matrix_otu, phylo_dend)

plot(pd_treefile$PD,pd_dend$PD)

#### WEIGHTED
mpd <- mpd(matrix_otu, cophenetic(treefile), abundance.weighted=TRUE)
mpd_dend <- mpd(matrix_otu, cophenetic(phylo_dend), abundance.weighted=TRUE)
mpd_dend_D0 <- mpd(matrix_otu_D0, cophenetic(phylo_dend), abundance.weighted=TRUE)
mpd_dend_D5 <- mpd(matrix_otu_D5, cophenetic(phylo_dend), abundance.weighted=TRUE)

boxplot(mpd_dend_D0,mpd_dend_D5)

##### on récupère les groupes
K=3
gpe = cutree(cah.01,k=K)
matrix_otu_gpe <- rbind(matrix_otu,gpe)
matrix_otu_gpe_species <- as.data.frame(merge(t(matrix_otu_gpe),species,by="row.names"))

# on va regarder la proportion pour chaque groupe des catégories de Phylum
group1 <- matrix_otu_gpe_species[matrix_otu_gpe_species[,"gpe"]==1,]
group1_prop <- apply(group1[,c(2:29)], 2,function(x) as.numeric(x/sum(x)))

group2 <- matrix_otu_gpe_species[matrix_otu_gpe_species[,"gpe"]==2,]
group2_prop <- apply(group2[,c(2:29)], 2,function(x) as.numeric(x/sum(x)))

group3 <- matrix_otu_gpe_species[matrix_otu_gpe_species[,"gpe"]==3,]
group3_prop <- apply(group3[,c(2:29)], 2,function(x) as.numeric(x/sum(x)))

group1_phylum=levels(group1$Phylum)
for(col in 2:29){
  res=aggregate(group1[,col],list(group1[,"Phylum"]),sum)
  group1_phylum=cbind(group1_phylum,res[,2])
}
group1_phylum <- data.frame(Phylum=group1_phylum[,1],t(as.data.frame(apply(as.data.frame(group1_phylum[,2:29]),1,FUN=as.numeric))),Groupe=rep("group1",4))

group2_phylum=levels(group2$Phylum)
for(col in 2:29){
  res=aggregate(group2[,col],list(group2[,"Phylum"]),sum)
  group2_phylum=cbind(group2_phylum,res[,2])
}
group2_phylum <- data.frame(Phylum=group2_phylum[,1],t(as.data.frame(apply(as.data.frame(group2_phylum[,2:29]),1,FUN=as.numeric))),Groupe=rep("group2",4))

group3_phylum=levels(group3$Phylum)
for(col in 2:29){
  res=aggregate(group3[,col],list(group3[,"Phylum"]),sum)
  group3_phylum=cbind(group3_phylum,res[,2])
}
group3_phylum <- data.frame(Phylum=group3_phylum[,1],t(as.data.frame(apply(as.data.frame(group3_phylum[,2:29]),1,FUN=as.numeric))),Groupe=rep("group3",4))


## on concatene les infos sur les 3 groupes
analyse_phylum_gpe <- rbind(group1_phylum,group2_phylum,group3_phylum) 

#on fait la moyenne de chaque ligne
moy <- apply(analyse_phylum_gpe[,2:29],1,mean)
analyse_phylum_gpe <- data.frame(analyse_phylum_gpe,moy)
#analyse_phylum_gpe <- data.frame(analyse_phylum_gpe[,c(1,30)],t(apply(analyse_phylum_gpe[,c(2:29,31)],1,FUN=as.integer)))

#X29 la moyenne
ggplot(analyse_phylum_gpe, aes(x=Groupe, y=moy, fill=Phylum))+
  geom_bar(stat="identity", position=position_dodge(), color="black")+
  labs(title="Répartition des bactéries dans les groupes \n en moyenne pour les 28 indiv", 
       x="Groupes issus CAH", y = "Nombre d'OTU")+
  scale_fill_brewer('PHYLUMS', palette='RdYlGn')+
  theme_minimal()

##calcul de la distance entre ls individu avec chaque groupe de bacteries.
beta_dist1<-as.matrix(vegdist(t(group1_prop), index='jaccard'))
beta_dist2<-as.matrix(vegdist(t(group2_prop), index='jaccard'))
beta_dist3<-as.matrix(vegdist(t(group3_prop), index='jaccard'))

colnames(beta_dist1)=colnames(beta_dist2)=colnames(beta_dist3)=substr(rownames(beta_dist1), 4,4)
rownames(beta_dist1)=rownames(beta_dist2)=rownames(beta_dist3)=as.factor(substr(rownames(beta_dist1), 6,7))

layout(matrix(c(1,2), nrow = T))
library(RColorBrewer)
#time factor
T_Factor1=T_Factor2=T_Factor3=NULL
for (i in unique(colnames(beta_dist1))){
  T_Factor1=c(T_Factor1,
             sum(beta_dist1[which(colnames(beta_dist1)==i),
                           which(colnames(beta_dist1)==i)])/2)
  T_Factor2=c(T_Factor2,
             sum(beta_dist2[which(colnames(beta_dist2)==i),
                           which(colnames(beta_dist2)==i)])/2)
  T_Factor3=c(T_Factor3,
             sum(beta_dist3[which(colnames(beta_dist3)==i),
                           which(colnames(beta_dist3)==i)])/2)
  
}
boxplot(tibble('Group 1'=T_Factor1,'Group 2'=T_Factor2,'Group 3'=T_Factor3), 
        main='Jaccard distance between \n D0/D5 in OTU groups', 
        col=brewer.pal(3, 'Set3'), cex.main =0.7)

#Inter-difference
T_Factor1=T_Factor2=T_Factor3=NULL
for (i in unique(colnames(beta_dist1))){
  T_Factor1=c(T_Factor1,
              mean(beta_dist1[which(rownames(beta_dist1)=="D0"),
                              which(colnames(beta_dist1)==i)][,1]))
  T_Factor2=c(T_Factor2,
              mean(beta_dist2[which(rownames(beta_dist2)=="D0"),
                              which(colnames(beta_dist2)==i)][,1]))
  T_Factor3=c(T_Factor3,
              mean(beta_dist3[which(rownames(beta_dist3)=="D0"),
                              which(colnames(beta_dist3)==i)][,1]))
}
boxplot(tibble('Group 1'=T_Factor1,'Group 2'=T_Factor2,'Group 3'=T_Factor3), 
        main='Jaccard mean distance between \n individus in OTU groups', 
        col=brewer.pal(3, 'Set3') , cex.main =0.7)



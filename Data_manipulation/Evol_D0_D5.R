# code pour la comparaison entre D0 et D5 pour les 14 individus
# - d'abord sur les indices d'alpha diversity
# - puis en utilisant groupes obtenus par matrice de dissimilarité de Jaccard
#  => évolution des modalités de Phylum sur barplot
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


#########################################
### stats générales sur D0 et D5  #######
#########################################

#OTU les plus fréquents : 
sort(colSums(matrix_otu_D0),decreasing = T)[1:10]
sort(colSums(matrix_otu_D5),decreasing = T)[1:10]
     # pas de changements sur les 4 premiers OTU qui restent majoritaires :
species[names(sort(colSums(matrix_otu_D5),decreasing = T)[1:4]),]

#matrice des différences : on regarde les individus avec les plus grandes évolutions 
matrix_otu_diff <- matrix_otu_D5-matrix_otu_D0
sort(abs(rowSums(matrix_otu_diff)),decreasing = T)
   ## !! très intéressant, certains individus ont un microbiote très stable et 
       # d'autres avec bcp de variations

#on regarde quels OTU sont reponsables d'une différence importante :
sort(matrix_otu_diff["BG-S-D5-CTL",],decreasing = T)[1:10] #cluster_3
sort(matrix_otu_diff["BG-G-D5-CUFF",],decreasing = T)[1:10]  #cluster_2  (++diff)
sort(matrix_otu_diff["BG-Q-D5-CTL",],decreasing = T)[1:10]  #cluster_6
sort(matrix_otu_diff["BG-F-D5-CTL",],decreasing = T)[1:10]  #cluster_25
sort(matrix_otu_diff["BG-A-D5-CUFF",],decreasing = T)[1:10]  #cluster_1
sort(matrix_otu_diff["BG-C-D5-CUFF",],decreasing = T)[1:10]  #cluster_111


###########################################
###    Alpha diversité        ############
###########################################

##### Rank abundance curve
#rankabundance(matrix_otu)

##### RICHESSE

#D0 et D5
richness_D0 <-apply(matrix_otu_D0, 1, function(x) length(which(x>=1)))
richness_D5 <-apply(matrix_otu_D5, 1, function(x) length(which(x>=1)))
plot(x=1:14,y=richness_D0, main="Evolution de la richesse entre D0 et D5 pour les 14 individus",ylab="Richesse")
points(x=1:14,y=richness_D5,col="red")
legend(legend=c("D0","D5"), 
       fill=c("black","red"), 'topleft',lty=0)
   #Pas de tendance claire vers une augmentation ou diminution de la richesse entre D0 et D5

##### SIMPSON
#D0 et D5
simpson_D0 <-vegan::diversity(matrix_otu_D0, index='simpson')
simpson_D5 <-vegan::diversity(matrix_otu_D5, index='simpson')
plot(x=1:14,y=simpson_D0, main="Evolution de simpson entre D0 et D5 pour les 14 individus",ylab="Simpson")
points(x=1:14,y=simpson_D5,col="red")
legend(legend=c("D0","D5"), 
       fill=c("black","red"), 'topleft',lty=0)
#très variable d'un individu à l'autre



##### SHANNON
#D0 et D5
shannon_D0 <-vegan::diversity(matrix_otu_D0, index='shannon')
shannon_D5 <-vegan::diversity(matrix_otu_D5, index='shannon')
plot(x=1:14,y=shannon_D0, main="Evolution de shannon entre D0 et D5 pour les 14 individus",ylab="Shannon")
points(x=1:14,y=shannon_D5,col="red")
legend(legend=c("D0","D5"), 
       fill=c("black","red"), 'topleft',lty=0)




###########################################
###    Beta diversité        ############
###########################################


##################################### 
######## Présence/abscence ########## 

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

matrix_otu_gpe_species <- as.data.frame(merge(t(matrix_otu_gpe),species,by="row.names"))


# on va regarder la proportion pour chaque groupe des catégories de Phylum
group1 <- matrix_otu_gpe_species[matrix_otu_gpe_species[,"gpe"]==1,]
group1_prop <- apply(group1[,c(2:29)], 2,function(x) as.numeric(x/sum(x)))
group1 <- data.frame(group1$Row.names,group1_prop*100,group1[,30:ncol(group1)])
group2 <- matrix_otu_gpe_species[matrix_otu_gpe_species[,"gpe"]==2,]
group2_prop <- apply(group2[,c(2:29)], 2,function(x) as.numeric(x/sum(x)))
group2 <- data.frame(group2$Row.names,group2_prop*100,group2[,30:ncol(group2)])
group3 <- matrix_otu_gpe_species[matrix_otu_gpe_species[,"gpe"]==3,]
group3_prop <- apply(group3[,c(2:29)], 2,function(x) as.numeric(x/sum(x)))
group3 <- data.frame(group3$Row.names,group3_prop*100,group3[,30:ncol(group3)])

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

#X29 la moyenne (pas trop de sens ?)
ggplot(analyse_phylum_gpe) + aes(x=Groupe, y=moy, fill=Phylum)+
  geom_bar(stat="identity", position=position_dodge()) +
  labs(title="Répartition des bactéries dans les groupes \n en moyenne pour les 28 indiv", 
       x="Groupes issus CAH", y = "Nombre d'OTU")

###### EVOLUTION D0 ET D5
ggplot(analyse_phylum_gpe) + aes(x=Groupe, y=X27, fill=Phylum)+
  geom_bar(stat="identity", position=position_dodge())+ylim(0,100)+
  labs(title="Répartition des bactéries dans les groupes \n pour individu S à D0", 
       x="Groupes issus CAH", y = "% d'OTU")+
  theme(plot.title = element_text(hjust = 0.5))
ggplot(analyse_phylum_gpe) + aes(x=Groupe, y=X28, fill=Phylum)+
  geom_bar(stat="identity", position=position_dodge())+ ylim(0,100)+
  labs(title="Répartition des bactéries dans les groupes \n pour individu S à D5", 
       x="Groupes issus CAH", y = "% d'OTU")+
  theme(plot.title = element_text(hjust = 0.5))

##calcul de la distance entre ls individu avec chaque groupe de bacteries.
beta_dist1<-as.matrix(vegdist(t(group1_prop), index='jaccard'))
beta_dist2<-as.matrix(vegdist(t(group2_prop), index='jaccard'))
beta_dist3<-as.matrix(vegdist(t(group3_prop), index='jaccard'))

plot(hclust(vegdist(t(group1_prop), index='jaccard'), method='ward.D'))
plot(hclust(vegdist(t(group2_prop), index='jaccard'), method='ward.D'))
plot(hclust(vegdist(t(group3_prop), index='jaccard'), method='ward.D'))

colnames(beta_dist1)=colnames(beta_dist2)=colnames(beta_dist3)=substr(rownames(beta_dist1), 4,4)
rownames(beta_dist1)=rownames(beta_dist2)=rownames(beta_dist3)=as.factor(substr(rownames(beta_dist1), 6,7))
layout(matrix(c(1,2), nrow = T))
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
        main='Jaccard mean distance between D0/D5 in OTU groups', 
        col=c('gold', 'tomato','green3'), cex.main=0.5)

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
        main='Jaccard mean distance between individu in OTU groups', 
        col=c('gold', 'tomato','green3'), cex.main=0.5)

##ON VA MAINTENANT SE CONCENTRER SUR LE GROUPE 2 ET LE GROUPE 1

morphological_data=data.frame(morphological_data)
experimental_condition$IDD=as.character(paste(experimental_condition$subject, experimental_condition$time))
morphological_data$IDD=as.character(paste(morphological_data$ID, morphological_data$Day))
data_exploration=merge(morphological_data,experimental_condition,by='IDD', x.all=F)

weight=cbind(as.numeric(data_exploration[which(data_exploration$Day=="D0"),10]), 
             as.numeric(data_exploration[which(data_exploration$Day=="D5"),10]))

delta=(weight[,1]-weight[,2])
#PCA
otu_group2=data.frame(t(group2_prop)[which(experimental_condition$time=='D0'),])
otu_group1=data.frame(t(group1_prop)[which(experimental_condition$time=='D0'),])
res.pca=FactoMineR::PCA(otu_group1)
mds.data=res.pca$ind$coord[,c(1,2)]
plot(mds.data)

color=delta
color[which(delta<=2.5)]=brewer.pal(6, 'YlOrRd')[6]
color[which(delta<2.1)]=brewer.pal(6, 'YlOrRd')[5]
color[which(delta<1.6)]=brewer.pal(6, 'YlOrRd')[4]
color[which(delta<1.1)]=brewer.pal(6, 'YlOrRd')[3]
color[which(delta<0.6)]=brewer.pal(6, 'YlOrRd')[2]

plot(mds.data,col=color)
plot(
  delta, 
  as.numeric(diversity(otu_group1,
                       index=c('simpson')))
)
library(randomForest)
otu_group1$class= as.numeric(delta)
model=randomForest::randomForest(class~., otu_group1, mtry=160, ntree=2000)
varImpPlot(model)










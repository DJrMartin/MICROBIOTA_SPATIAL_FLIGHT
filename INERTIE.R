

#### 7) DIVERSITY OF CO-OCCURENCE
#on affecte 0 et 1 en fonction de si abscence ou présence de l'otu
matrix_otu_01 <- matrix_otu>0

#distance de Jaccard :
Jaccard = function(i,j){
  return(sum(i!=j)/(sum(i!=j)+sum((i==1 & j==1))))
}
#distance de Dice (double poids accords)
Dice = function(i,j){
  return(sum(i!=j)/(sum(i!=j)+2*sum((i==1 & j==1))))
}

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
rect.hclust(cah.01, 3, border ="red")


gpe <- cutree(cah.01,k=3)
matrix_otu_groupe <- rbind(matrix_otu_01,as.factor(gpe))
matrix_otu_groupe <- as.data.frame(t(matrix_otu_groupe))

#devrait marcher mais non..
acp <- PCA(matrix_otu_groupe,scale.unit=TRUE,graph=F,quali.sup=29)
plot(acp,choix="var")
plot(acp,choix="ind",habillage = 29, col.hab=c("green","blue","red"))

ind=which(matrix_otu[3,]>0)
length(ind)

pca <- PCA(matrix_otu_groupe[,-29],scale.unit=TRUE,graph=F)
x <- pca$ind$coord[ind, 1L] # Première dimension de la PCA
y <- pca$ind$coord[ind, 2L] # Deuxième dimension de la PCA
z <- pca$ind$coord[ind, 3L] # Troisième dimension de la PCA
#col <- matrix_otu_groupe[,29]
#plot(x, y, col = col,xlim=c(),ylim=c())

gpe_correct=gpe[ind]

library('scatterplot3d')
xlim=range(x)
ylim=range(y)
zlim=range(z)
s3d <- scatterplot3d(x[which((gpe_correct==1)==TRUE)],y[which((gpe_correct==1)==TRUE)],z[which((gpe_correct==1)==TRUE)],
                     xlim=xlim, ylim=ylim, zlim=zlim, xlab="PC1 (32%)", ylab="PC2 (9%)", zlab='PC3', main='INDIVIDU 3')
s3d$points3d(x[which((gpe_correct==2)==TRUE)],y[which((gpe_correct==2)==TRUE)],z[which((gpe_correct==2)==TRUE)],col = "red", pch = 2)
s3d$points3d(x[which((gpe_correct==3)==TRUE)],y[which((gpe_correct==3)==TRUE)],z[which((gpe_correct==3)==TRUE)],col = "green3", pch = 5)

#barycentres
matrix_grp1 <- pca$ind$coord[matrix_otu_groupe[,29]==1,]
bar_grp1 <- colMeans(matrix_grp1[,-29])
matrix_grp2 <- pca$ind$coord[matrix_otu_groupe[,29]==2,]
bar_grp2 <- colMeans(matrix_grp2[,-29])
matrix_grp3 <- pca$ind$coord[matrix_otu_groupe[,29]==3,]
bar_grp3 <- colMeans(matrix_grp3[,-29])

##Calcul de l'inertie intraclasse pour chaque individu sur chaque composante

ind_01=which(matrix_otu[1,]>0)
length(ind_01)
(sum(pca$ind$coord[ind_01, 1][which((gpe_correct==1)==TRUE)]-bar_grp1[1])^2)/length(ind_01)

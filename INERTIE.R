#### 7) DIVERSITY OF CO-OCCURENCE
load('~/Dropbox/Spatial_flight/RData_Microgravity_updata.RData')

#on affecte 0 et 1 en fonction de si abscence ou présence de l'otu
nb <- colSums(matrix_otu)
matrix_otu <- matrix_otu[,nb!=0]  #10 OTU supprimés
species <- species[nb!=0,]
matrix_otu_01 <- matrix_otu>0

#distance de Jaccard :
Jaccard = function(i,j){
  return(sum(i!=j)/(sum(i!=j)+sum((i==1 & j==1))))
}

MatDiss = outer(as.data.frame(matrix_otu_01),as.data.frame(matrix_otu_01),Vectorize(Jaccard))
dist_01 <- as.dist(MatDiss)

# on fait une CAH, avec stratégie d'aggrégation de Ward
cah.01 <- hclust(dist_01,method="ward.D")
plot(cah.01,labels = F)

# Regardons la courbe de perte d'inertie 
plot(rev(cah.01$height)[1:20],type="b",main="inertie intra")
# on peut conserver 3 groupes
plot(cah.01,labels = F)
rect.hclust(cah.01, 3, border ="red")


gpe <- cutree(cah.01,k=3)
matrix_otu_groupe <- rbind(matrix_otu_01,as.factor(gpe))
matrix_otu_groupe <- as.data.frame(t(matrix_otu_groupe))

groupe1=matrix_otu[,which(gpe==3)]
color_visualisation=c('purple', 'tomato', 'lightgreen','lightblue' )
OTU_phylum=NULL
for (i in as.character(unique(species$Phylum))){
  sum_reads=apply(groupe1[, which(species$Phylum[which(gpe==3)]==i)], 1,sum)
  OTU_phylum=cbind(OTU_phylum,sum_reads)
}
rownames(OTU_phylum)=paste(substr(rownames(OTU_phylum), 4,4), rep(c('b','a'), 14))
colnames(OTU_phylum)=as.character(unique(species$Phylum))
barplot(t(OTU_phylum), col=color_visualisation, 
        ylim=c(range(OTU_phylum)[1],range(OTU_phylum)[2]*1.6))
legend(legend=c(as.character(unique(species$Phylum))), 
       fill=color_visualisation, 'topleft', ncol=2,
       cex=0.6, bty='n', lty=0)

ind=which(matrix_otu_01[25,]==TRUE)
length(ind)

library(FactoMineR)
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
                     xlim=xlim, ylim=ylim, zlim=zlim, xlab="PC1 (32%)", ylab="PC2 (9%)", zlab='PC3', main='INDIVIDU Q')
s3d$points3d(x[which((gpe_correct==2)==TRUE)],y[which((gpe_correct==2)==TRUE)],z[which((gpe_correct==2)==TRUE)],col = "red", pch = 2)
s3d$points3d(x[which((gpe_correct==3)==TRUE)],y[which((gpe_correct==3)==TRUE)],z[which((gpe_correct==3)==TRUE)],col = "green3", pch = 5)

#barycentres
bc=NULL
for (grp in 1:3){
  matrix_grp1 <- pca$ind$coord[matrix_otu_groupe[,29]==grp,]
  bc=rbind(bc,colMeans(matrix_grp1[,-29]))
}

##Calcul de l'inertie intra-classe pour chaque individu sur chaque composante WEIGTHED
nk =5
ngrp = 3
PCs =PCs_group=PCs_inertie=INERTIE_sum=INERTIE_local=NULL
for(i in 1:28){
  x=which(matrix_otu[i,]>0)
  for(grp in 1:ngrp){
    for(k in 1:nk){
      gpe_correct=gpe[x]
      PCs=c(PCs,sum((pca$ind$coord[x, k][which((gpe_correct==grp)==TRUE)]-bc[grp,k])^2)/length(x))
    }
    INERTIE_local=c(INERTIE_local,sum(PCs))
    PCs_group=c(PCs_group,PCs)
    PCs=NULL
  }
  INERTIE_sum=rbind(INERTIE_sum,INERTIE_local)
  PCs_inertie=rbind(PCs_inertie,PCs_group)
  PCs_group=INERTIE_local=NULL
}

INERTIE=data.frame('ID' = paste(rep(c(as.character(unique(experimental_condition$subject))),each=2),rep(c('b','a'),2)),
                   PCs_inertie)
rownames(INERTIE)=rownames(INERTIE_sum)=INERTIE$ID
colnames(INERTIE)=c('ID', paste0(rep("CLUSTER_1_", 5),seq(1,5,by=1)),paste0(rep("CLUSTER_2_", 5),seq(1,5,by=1)),
                    paste0(rep("CLUSTER_3_", 5),seq(1,5,by=1)))

colnames(INERTIE_sum)=c('Dim_1','Dim_2','Dim_3')

#Masse maigre (abs et rel)
library(randomForest)
load('~/Dropbox/Spatial_flight/RData_Microgravity_updata.RData')

lM1=lM[-1,]
colnames(lM1)[1]='ID'
lM1$ID=paste(as.character(lM1$ID), rep(c(1,2), each=1))
experimental_condition$ID=paste(as.character(experimental_condition$subject), rep(c(1,2), each=1))
lM2=merge(lM1,experimental_condition, by="ID")


INERTIE_D0=data.frame(INERTIE[which(lM2$time=='D0'),-1])
set.seed(123)
y=as.numeric(as.character(lM2$Whole.Body.Lean.mass[which(lM2$time=='D0')]))>55000
res.pca=PCA(INERTIE_D0)
plot(res.pca$ind$coord[,c(1,2)], col=as.numeric(y)+1)
INERTIE_D0$Y=as.factor(y)

#leaveONEout
x_test=prediction=NULL
cnt=1
while (cnt<20){
  w=sample(1:14, 1)
  rf=randomForest(Y~., INERTIE_D0[-c(w),])
  prediction=c(prediction,predict(rf,INERTIE_D0[w,], type='prob')[,1])
  x_test=c(x_test,INERTIE_D0$Y[c(w)])
  cnt=cnt+1
}

plot(prediction~x_test, xlab='> ou < à 55kg', ylab='Probabilité de prédiction')
pROC::auc(x_test,prediction)
summary(lm((x_test~prediction)))
varImpPlot(rf)

#VO2 max (abs et rel)
morphological_data=data.frame( morphological_data)
morphological_data$ID=paste(as.character(morphological_data$ID), rep(c(1,2), each=1))
experimental_condition$ID=paste(as.character(experimental_condition$subject), rep(c(1,2), each=1))
VO2=merge(morphological_data,experimental_condition, by="ID")

set.seed(123)
INERTIE_D0=data.frame(INERTIE[which(VO2$time=='D0'),])
INERTIE_D0$Y=as.numeric(VO2$VO2max[which(VO2$time=='D0')])/
  as.numeric(VO2$Weight[which(VO2$time=='D0')])>45
INERTIE_D0$Y=as.factor(INERTIE_D0$Y)

#leaveONEout
x_test=prediction=NULL
cnt=1
while (cnt<50){
  w=sample(1:14,2)
  rf=randomForest(Y~., INERTIE_D0[-c(w),-1])
  prediction=c(prediction,predict(rf,INERTIE_D0[w,-c(1,17)], type = 'prob')[,1])
  x_test=c(x_test,INERTIE_D0$Y[c(w)])
  cnt=cnt+1
}

boxplot(prediction~x_test,xlab='< ou > à 45mL/min/kg', ylab='Probabilité de prédiction')
summary(glm((x_test~prediction)))
pROC::auc(x_test,prediction)

###Démontrer la méthode
#Prenons l'exemple d'un cluster de 100 OTUs pour montrer l'intérêt de la méthode
set.seed(1)
abondance=dbinom(x = 1:100, size=100, 0.6)

x <- rnorm(100, mean=0, sd=3)
y <- rnorm(100, mean=0, sd=3)
z <- rnorm(100, mean=0, sd=3)

abundance=data.frame(x,y,z,abondance)

s3d <- scatterplot3d(x,y,z,main='CLUSTERING EXAMPLE')

bc=apply(cbind(x,y,z),2,mean)

s3d$points3d(bc[1],bc[2],
             bc[3],col='red', pch=16, cex=2)

PC1.max=sum((x-bc[1])^2)/length(x)
PC2.max=sum((y-bc[1])^2)/length(x)
PC3.max=sum((z-bc[1])^2)/length(x)
sum(PC1.max,PC2.max,PC3.max)
shannon_total=diversity(abondance, 'shannon')
simpson_total=diversity(abondance, 'simpson')

set.seed(1)
i=1
PC1.fluctuation=PC2.fluctuation=PC3.fluctuation=inertie=sh=sp=rch=NULL

while(i < 2000){
  random=sample(1:100)[1:50]
  new=abundance[random,]
  PC1.fluctuation=c(PC1.fluctuation,sum((new$x-bc[1])^2)/length(new$x))
  PC2.fluctuation=c(PC2.fluctuation,sum((new$y-bc[1])^2)/length(new$y))
  PC3.fluctuation=c(PC3.fluctuation,sum((new$z-bc[1])^2)/length(new$z))
  inertie=c(inertie,sum(sum((new$x-bc[1])^2)/length(new$x),
                        sum((new$y-bc[1])^2)/length(new$y),
                        sum((new$z-bc[1])^2)/length(new$z)))
  sh=c(sh,diversity(new$abondance, 'shannon'))
  sp=c(sp,diversity(new$abondance, 'simpson'))
  rch=c(rch,50)
  i=i+1
}

library(pheatmap)
pheatmap(scale(data.frame(PC1.fluctuation,PC2.fluctuation,PC3.fluctuation,
                          sum_inertie=inertie,
                          'Shannon'=sh,'Simpson'=sp)))


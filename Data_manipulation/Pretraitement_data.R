library(RColorBrewer)
library(vegan)
library(fossil)
library(ggplot2)
###Import#####################################################################
User='David'
if (User=='David'){
  load('~/Dropbox/Spatial_flight/RData_Microgravity_updata.RData')
}
rm(User)

###Courbe de rarefaction######################################################
barplot(apply(matrix_otu, 1, sum), ylab='Profondeur de sequencage', xlab='sujets')
#rarecurve(matrix_otu, step=20)

#On peut voir que chacun des sujets se stabilise voulant dire que la profondeur de séqeunçage
#est suffisante pour observer l'ensemble des espèces présentes dans les échantillons.

###Normalisation##############################################################
#Approche = BRUTE, TSS, CLR, ALR, CSS ou autres (mais il faut les implémenter)
#Geometric Mean of Pairwise Ratio method ??
OTU_CSS<-NORMALISATION_MICROBIOTA(matrix_otu, approches = 'CSS')
OTU_TSS<-NORMALISATION_MICROBIOTA(matrix_otu, approches = 'TSS')
OTU_BRUTES<-NORMALISATION_MICROBIOTA(matrix_otu, approches = 'BRUTES')
OTU_CLR<-NORMALISATION_MICROBIOTA(matrix_otu, approches = 'CLR')

OTU_normalised=OTU_TSS
###Visualisation en fonction des phylums######################################
dev.off()
color_visualisation=c('red', 'grey33', 'green3','lightblue' )
OTU_phylum=NULL
for (i in as.character(unique(species$Phylum))){
  sum_reads=apply(OTU_normalised[, which(species$Phylum==i)], 1,sum)
  OTU_phylum=cbind(OTU_phylum,sum_reads)
}
colnames(OTU_phylum)=as.character(unique(species$Phylum))
barplot(t(OTU_phylum), col=color_visualisation, 
        ylim=c(range(OTU_phylum)[1],range(OTU_phylum)[2]*1.6))
legend(legend=c(as.character(unique(species$Phylum))), 
       fill=color_visualisation, 'topleft', ncol=2,
       cex=0.6, bty='n', lty=0)

###Visualisation en fonction des class#########################################
color_visualisation=c(brewer.pal(11, "Set3"))
OTU_class=NULL
for (i in as.character(unique(species$Class))){
  sum_reads=apply(OTU_normalised[, which(species$Class==i)], 1,sum)
  OTU_class=cbind(OTU_class,sum_reads)
}
OTU_class=cbind(OTU_class,OTU_normalised[, which(species$Class=='Flavobacteriia')])
colnames(OTU_class)=as.character(unique(species$Class))
barplot(t(OTU_class), col=color_visualisation, 
        ylim=c(range(OTU_class)[1],range(OTU_class)[2]*1.5))
legend(legend=c(as.character(unique(species$Class))), 
       fill=color_visualisation, 'topleft', ncol=3,
       cex=0.6, bty='n', lty=0)

###Diversity###################################################################

shannon=diversity(OTU_normalised[which(experimental_condition$time=="D0"), ], index='simpson')

#Time factor
###Beta
#Disimilarity between time factor
beta_total=NULL
for (sample in unique(experimental_condition$subject)){
  beta=vegdist(OTU_normalised[which(experimental_condition$subject==sample),], index='bray')
  beta_total=c(beta_total, beta)
}
time_factor=tibble::tibble('ID'=unique(experimental_condition$subject), 'B-diversity'=beta_total, 
               'Condition'=as.factor(substr(experimental_condition$Sample, 9,12)[-c(2,4,6,8,10,12,14,16,18,20,22,24,26,28)]))
plot(time_factor$`B-diversity`~time_factor$Condition,
     xlab='conditions', ylab='Bray-Curtis diversity') #p=0.37 (t.test)
summary(lm(`B-diversity`~Condition, data=time_factor))

###Visualisation de l'entropy intra-class.
ENTROPY=NULL
otu=NULL
filtre=species[,4]
for (j in as.character(unique(filtre))){
  reads=OTU_normalised[, which(filtre==j)]
  ENTROPY=cbind(ENTROPY, diversity(reads, index='shannon'))
}
colnames(ENTROPY)=colnames(as.character(unique(filtre)))
rownames(ENTROPY)=rownames(matrix_otu)
ENTROPY=cbind(ENTROPY, experimental_condition$time)

Bacterio=cbind(ENTROPY[which(experimental_condition$time=="D0"),1], ENTROPY[which(experimental_condition$time=="D5"),1])
firmicutes=cbind(ENTROPY[which(experimental_condition$time=="D0"),3], ENTROPY[which(experimental_condition$time=="D5"),3])
ylim=range(firmicutes)
plot(firmicutes[1,], type ='l', ylim=ylim, main = 'E of Bacteriodetes',
     ylab='Entropy', xlab='Time')
for(i in 2:length(firmicutes[,1])){
  lines(firmicutes[i,])
}

##Evolution of POWER max during exercice
morphological_data=data.frame(morphological_data)
experimental_condition$IDD=as.character(paste(experimental_condition$subject, experimental_condition$time))
morphological_data$IDD=as.character(paste(morphological_data$ID, morphological_data$Day))
data_exploration=merge(morphological_data,experimental_condition,by='IDD', x.all=F)

weight=cbind(as.numeric(data_exploration[which(data_exploration$Day=="D0"),10]), 
             as.numeric(data_exploration[which(data_exploration$Day=="D5"),10]))

delta=(weight[,1]-weight[,2])

res.pca=FactoMineR::PCA(OTU_normalised[which(experimental_condition$time=='D0'),])
mds.data=res.pca$ind$coord[,c(1,2)]
plot(mds.data)

color=delta
color[which(delta<=2.5)]=brewer.pal(6, 'YlOrRd')[6]
color[which(delta<2.1)]=brewer.pal(6, 'YlOrRd')[5]
color[which(delta<1.6)]=brewer.pal(6, 'YlOrRd')[4]
color[which(delta<1.1)]=brewer.pal(6, 'YlOrRd')[3]
color[which(delta<0.6)]=brewer.pal(6, 'YlOrRd')[2]

plot(mds.data,col=color)
beta_dist<-vegdist(OTU_normalised, index='jaccard')
mds <- metaMDS(beta_dist)
mds_data <- as.data.frame(mds$points)
mds_data$SampleID <- rownames(mds_data)
mds_data$Time=experimental_condition$time
mds_data$Sujets=experimental_condition$subject
plot(mds_data[which(experimental_condition$time=="D5"),c(1,2)],col=color)

plot(firmicutes[,2]-firmicutes[,1], delta)
summary(lm(ENTROPY[which(experimental_condition$time=="D0"),3]~ delta))

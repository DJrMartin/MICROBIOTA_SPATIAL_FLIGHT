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


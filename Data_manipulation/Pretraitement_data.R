library(RColorBrewer)
library(vegan)
###Import#####################################################################
User='David'
if (User=='David'){
  load('C:/Users/David/Desktop/MICROBIOTA_SPATIAL_FLIGHT/DATA_MICROGRAVITY_PROJECT.RData')
}
if (User=='Marion'){
  load()
}
if (User=='Valerie'){
  load()
}
rm(User)

###Courbe de rarefaction######################################################
barplot(apply(matrix_otu, 1, sum), ylab='Profondeur de sequencage', xlab='sujets')
rarecurve(matrix_otu, step=20)
#On peut voir que chacun des sujets se stabilise voulant dire que la profondeur de séqeunçage
#est suffisante pour observer l'ensemble des espèces présentes dans les échantillons.

###Normalisation##############################################################
source('~/')
#Approche = BRUTE, TSS, CLR, ALR, CSS ou autres (mais il faut les implémenter)
OTU_normalised<-NORMALISATION_MICROBIOTA(matrix_otu, approches = 'CSS')
rownames(OTU_normalised)=rownames(matrix_otu)

###Visualisation en fonction des phylums######################################
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
#richesse
ylim=range(apply(OTU_normalised, 1, function(x) length(which(x>=1))), apply(matrix_otu, 1, function(x) length(which(x>=10))))
plot(apply(OTU_normalised, 1, function(x) length(which(x>=1))), ylim=ylim, ylab='Richness', xlab='ind')
points(apply(OTU_normalised, 1, function(x) length(which(x>=5))), col='red')
points(apply(OTU_normalised, 1, function(x) length(which(x>=10))), col='green3')
legend('topleft', legend=c(1,5,10), fill = c('black', "red", "green3"), title='FILTRES', horiz = T)
#Est-ce que l'on filtre les OTUs peu exprimés ? et avec quelle normalisation ?

###Alpha
#Chao1

#Simpson

#Shannon


###Beta
#Disimilarity between time factor

#Disimilarity between 

##log poisson##################################################################

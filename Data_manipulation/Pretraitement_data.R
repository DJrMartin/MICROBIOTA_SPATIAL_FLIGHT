library(RColorBrewer)
library(vegan)
library(fossil)
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
#rarecurve(matrix_otu, step=20)

#On peut voir que chacun des sujets se stabilise voulant dire que la profondeur de séqeunçage
#est suffisante pour observer l'ensemble des espèces présentes dans les échantillons.

###Normalisation##############################################################
source('https://raw.githubusercontent.com/DJrMartin/MICROBIOTA_SPATIAL_FLIGHT/Microgravity_workflow/Data_manipulation/Functions_normalisation.R?token=GHSAT0AAAAAABTX2ZJWQGIXWW542J3QFI4CYTCPHRA')
#Approche = BRUTE, TSS, CLR, ALR, CSS ou autres (mais il faut les implémenter)
OTU_CSS<-NORMALISATION_MICROBIOTA(matrix_otu, approches = 'CSS')
OTU_TSS<-NORMALISATION_MICROBIOTA(matrix_otu, approches = 'TSS')
OTU_BRUTES<-NORMALISATION_MICROBIOTA(matrix_otu, approches = 'BRUTES')
OTU_CLR<-NORMALISATION_MICROBIOTA(matrix_otu, approches = 'CLR')

layout(matrix(c(1,2,
                3,4), nrow=2))
hist(as.numeric(OTU_CSS[2,]), breaks =60, main="CSS", xlab='Distribution OTU')
hist(as.numeric(OTU_TSS[2,]),  breaks =60, main ="TSS", xlab='Distribution OTU')
hist(as.numeric(OTU_BRUTES[2,]),  breaks =60, main ="BRUTES", xlab='Distribution OTU')
hist(as.numeric(OTU_CLR[2,]),  breaks =60, main ="CLR", xlab='Distribution OTU')

###Visualisation en fonction des phylums######################################
color_visualisation=c('red', 'grey33', 'green3','lightblue' )
OTU_phylum=NULL
for (i in as.character(unique(species$Phylum))){
  sum_reads=apply(OTU_CSS[, which(species$Phylum==i)], 1,sum)
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
  sum_reads=apply(OTU_CSS[, which(species$Class==i)], 1,sum)
  OTU_class=cbind(OTU_class,sum_reads)
}
OTU_class=cbind(OTU_class,OTU_CSS[, which(species$Class=='Flavobacteriia')])
colnames(OTU_class)=as.character(unique(species$Class))
barplot(t(OTU_class), col=color_visualisation, 
        ylim=c(range(OTU_class)[1],range(OTU_class)[2]*1.5))
legend(legend=c(as.character(unique(species$Class))), 
       fill=color_visualisation, 'topleft', ncol=3,
       cex=0.6, bty='n', lty=0)

###Diversity###################################################################
#richness on data BRUT
ylim=range(apply(OTU_normalised, 1, function(x) length(which(x>=1))), apply(matrix_otu, 1, function(x) length(which(x>=10))))
plot(apply(OTU_normalised, 1, function(x) length(which(x>=1))), ylim=ylim, ylab='Richness', xlab='ind')
points(apply(OTU_normalised, 1, function(x) length(which(x>=5))), col='red')
points(apply(OTU_normalised, 1, function(x) length(which(x>=10))), col='green3')
legend('topleft', legend=c(1,5,10), fill = c('black', "red", "green3"), title='FILTRES', horiz = T)

RICHNESS<-apply(matrix_otu, 1, function(x) length(which(x>=1)))
#Est-ce que l'on filtre les OTUs peu exprimés ? et avec quelle normalisation ?

###Alpha
#Chao1 on data BRUT, estimation of richness (improve this metrics ??)
#Chao, A. 1984. Non-parametric estimation of the number of classes in a population. 
#Scandinavian Journal of Statistics 11: 265-270
plot(apply(matrix_otu, 1, function(x) length(which(x>=3)))
     ~apply(matrix_otu, 1, function(x) sum(x)))
#On peut se poser la question de savoir s'il est pertinent d'utiliser Chao1 dans notre cas
#On remarque bien que la profondeur de séquençage n'influe pas sur la présence/absence
#de certaines espèces (au contraire)
CHAO1<-apply(matrix_otu, 1, function(x) chao1(x, taxa.row = T))
#Simpson
SIMPSON=diversity(matrix_otu, index='simpson')
#Shannon
SHANNON<-diversity(matrix_otu, index='shannon')
plot(RICHNESS, SHANNON)

#Time factor


###Beta
#Disimilarity between time factor

#Disimilarity between 

##log poisson##################################################################

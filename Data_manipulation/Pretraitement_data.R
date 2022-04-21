###Import
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

###Courbe de rarefaction#############################
library('vegan')
barplot(apply(matrix_otu, 1, sum), ylab='profondeur de sequencage', xlab='sujets')
rarecurve(matrix_otu, step=20)
#On peut voir que chacun des sujets se stabilise voulant dire que la profondeur de séqeunçage
#est suffisante pour observer l'ensemble des espèces présentes dans les échantillons.

###Normalisation#####################################
library('compositions')
#normalisation center logratio
norm_otu_table<-as.data.frame(clr(matrix_otu))
image(as.matrix(matrix_otu))
image(as.matrix(norm_otu_table))

plot(hclust(dist(matrix_otu)))
plot(hclust(dist(norm_otu_table)))

#normalisation pourcentage
pourcentage_otu_table<-as.data.frame(apply(matrix_otu, 2, function(x) x/sum()))

###Diversity##########################################
#richesse
ylim=range(apply(matrix_otu, 1, function(x) length(which(x>=1))), apply(matrix_otu, 1, function(x) length(which(x>=10))))
plot(apply(matrix_otu, 1, function(x) length(which(x>=1))), ylim=ylim, ylab='Richness', xlab='ind')
points(apply(matrix_otu, 1, function(x) length(which(x>=5))), col='red')
points(apply(matrix_otu, 1, function(x) length(which(x>=10))), col='green3')
legend('topleft', legend=c(1,5,10), fill = c('black', "red", "green3"), title='FILTRES', horiz = T)
#Est-ce que l'on filtre les OTUs peu exprimés ? 
#Globalement on garde le même caractètre, ce qui est assez rassurant.

#Alpha

#Beta

##log poisson#########################################

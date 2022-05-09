########## PD diversity
library(tree)
library(ape)  
library(phangorn)
library(phytools)
library(geiger)
library(picante)

load('C:/Users/33638/Documents/stage/MICROBIOTA_SPATIAL_FLIGHT/DATA_PROJECT_1.RData')


class(treefile)  #phylo


plot(treefile)
plot.phylo(treefile)

plotTree(treefile,type="fan",fsize=0.7,lwd=1,
         ftype="i")


# The pd function returns two values for each community, Faith's PD and species richness (SR).
pd <- pd(matrix_otu, treefile)
pd
plot(x=1:28,y=pd$PD,main="PD ")


# positively correlates with species richness across samples.
plot(x=pd$SR,y=pd$PD,main="PD / SR")



# pour prendre en compte l'abondance:
# Calculates mean pairwise distance separating taxa in a community
mpd <- mpd(matrix_otu, cophenetic(treefile), abundance.weighted=TRUE)
             #cophenetic : distance entre deux objets calculés à partir d'un dendrogramme

plot(x=1:28,y=mpd,main="PD ")
plot(x=pd$PD,y=mpd)


##Creation du RData avec l'ensemble des données
setwd(dir='~/Dropbox/Spatial_flight/')
##Tree
metadata<- read.table("MetaCNES_D0_D5paired.tsv", sep = ";", row.names=1, header=TRUE, fill = TRUE, stringsAsFactors = TRUE)
##BIOM
frogs.data <- import_frogs("abundance.biom1")
sample_data(frogs.data) <- metadata
frogs.data

matrix_otu<-t(as.data.frame(otu_table(frogs.data)))
species<-as.data.frame(tax_table(frogs.data))
experimental_condition<-metadata

##Metalome
metalome=read.table("metalome.csv", sep = ";", stringsAsFactors = TRUE, header=TRUE)
metalome=as.matrix(gsub(',', '.', apply(metalome, 2, as.character)))

#Data_anthropometrics
morphological_data=read.table("Anthropo.csv", sep = ";", stringsAsFactors = TRUE, header=TRUE)
morphological_data=as.matrix(gsub(',', '.', apply(morphological_data, 2, as.character)))

##Data IRM
IRM_Fe_content=read.table("IRM_Fe_content.csv", sep = ";", stringsAsFactors = TRUE, header=TRUE)

load("~/Dropbox/Spatial_flight/RData_Microgravity_updata.RData")
##Lean Mass
lM=read.table("LEAN_MASS.csv", sep = ";", stringsAsFactors = TRUE, header=TRUE)

save(lM,metalome,morphological_data,IRM_Fe_content,matrix_otu, species, 
     experimental_condition, metalome, 
     file='~/Dropbox/Spatial_flight/RData_Microgravity_updata.RData')


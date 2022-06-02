
load('C:/Users/33638/Documents/stage/MICROBIOTA_SPATIAL_FLIGHT/DATA_PROJECT_1.RData')

################################################
####### Output : données morphologiques ########
################################################
### Output d'intérêts:
# - VO2max (en absolue)
# - VO2max (en relatif, /weight)
# - Weight
# - Masse maigre (en absolue)
# - Masse maigre (en relatif, /weight)
# - Isocinetique (force max)

#### A L'ETAT BASAL, les 14 individus ##############

morph_data_D0 <- as.data.frame(morphological_data[morphological_data[,2]=='D0',][c(-5,-8,-12,-17),c(1,3,9)])
morph_data_D0$VO2max <- as.numeric(morph_data_D0$VO2max)
morph_data_D0$Weight <- as.numeric(morph_data_D0$Weight)
morph_data_D0 <- data.frame(morph_data_D0,VO2max_rel=morph_data_D0$VO2max/morph_data_D0$Weight)
row.names(morph_data_D0) <- morph_data_D0$ID
morph_data_D0 <- morph_data_D0[,-1]

# Stats desc
summary(morph_data_D0)

# On crée les dataframe sur lesquels faire modèles
nb <- colSums(matrix_otu)
matrix_otu <- matrix_otu[,nb!=0]  #10 OTU supprimés

    # sans regrouper, sur tous les OTU
matrix_otu_D0 <- matrix_otu[seq(1,28,by=2),]
# On merge les deux df
row.names(matrix_otu_D0)<-row.names(morph_data_D0)
otu_output <- merge(matrix_otu_D0,morph_data_D0,by="row.names")
row.names(otu_output)<- otu_output$Row.names
otu_output <- otu_output[,-1]

    # en regroupant selon Phylum
matrix_otu_species <- merge(t(matrix_otu),species,by="row.names")
matrix_phylum=levels(matrix_otu_species$Phylum)
for(col in 2:29){
  res=aggregate(matrix_otu_species[,col],list(matrix_otu_species[,"Phylum"]),sum)
  matrix_phylum=data.frame(matrix_phylum,res[,2])
}
row.names(matrix_phylum)<-matrix_phylum[,1]
matrix_phylum <- matrix_phylum[,-1]
matrix_phylum_D0 <- t(matrix_phylum[,seq(1,28,by=2)])
row.names(matrix_phylum_D0)<-row.names(morph_data_D0)

# On merge les deux df
phylum_output <- merge(matrix_phylum_D0,morph_data_D0,by="row.names")
row.names(phylum_output)<-phylum_output$Row.names
phylum_output <- phylum_output[,-1]


######## Modélisation ##################
library(corrplot)
corr=cor(phylum_output)
corrplot(corr,type = "upper", tl.col = "black", tl.srt = 45)

# Régression linéaire multiple : sur données réduites Phylum
mod_VO2max_rel <- lm(VO2max_rel~.,data=phylum_output[,c(1:4,7)])
summary(mod_VO2max_rel)

mod_VO2max <- lm(VO2max~.,data=phylum_output[,1:5])
summary(mod_VO2max)

mod_weight <- lm(Weight~.,data=phylum_output[,c(1:4,6)])
summary(mod_weight)

# Forêts aléatoires--------------------------
library(randomForest)
ind=sample(1:14,10)
train <- otu_output[ind,]
test <- otu_output[-ind,]

set.seed(123) # permet de fixer les paramètres aléatoires de la rf
rf=randomForest(VO2max~ . , data = train[,c(-619,-620)],importance=T,ntree=500)
rf
plot(rf) #200 arbres suffisent

rf.results <- predict(rf,test)
results <- data.frame(actual = test$VO2max, prediction = rf.results)
head(results)

varImpPlot(rf)
    #%IncMSE : correspond au “poids” de la variable dans la diminution de l’erreur par le modèle.
    #IncNodePurity : correspond à la capacité de la variable à bien discriminer les individus à son noeud

 
# Grande dimension : voir cours Laurent Rouviere





#### DELTA D0-D5, les 14 individus ##############

morph_data_D5 <- as.data.frame(morphological_data[morphological_data[,2]=='D5',][c(-5,-8,-12,-17),c(1,3,9)])
morph_data_D5$VO2max <- as.numeric(morph_data_D5$VO2max)
morph_data_D5$Weight <- as.numeric(morph_data_D5$Weight)
morph_data_D5 <- data.frame(morph_data_D5,VO2max_rel=morph_data_D5$VO2max/morph_data_D5$Weight)
row.names(morph_data_D5) <- morph_data_D5$ID
morph_data_D5 <- morph_data_D5[,-1]

morph_delta <- morph_data_D0-morph_data_D5

## stat desc 
   # corrélation ?


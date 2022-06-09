###Utilisation du package PLNmodels
#https://oliviergimenez.github.io/code_livre_variables_cachees/chiquet.html

library('PLNmodels')
library(corrplot)
library(dplyr)
library(tidyverse)

 ############## PREMIER ESSAI DAVID ############
###Preparation data
experimental_condition$Conditions=as.factor(substr(experimental_condition$Sample, 9,12))
PLN_data<-prepare_data(counts = matrix_otu, covariates = experimental_condition[, c(1,2,4)])

###Model de base

PLN_model=PLN(Abundance~ 0 + time + Conditions, data=PLN_data)

Coefficient<-PLN_model %>% 
  coefficients()
filtre=Coefficient[,2] %>% order(decreasing = T)

Coefficient[filtre[1:10],2]%>% t() %>% round(1) %>% 
  corrplot(is.corr = FALSE, method = 'color', tl.cex = .5, cl.pos = "n")

plot(matrix_otu[,filtre[2]]~experimental_condition$time)
species[filtre[2],]

###Reduction de dimensionalite
PLN_PCA <- PLNPCA(Abundance~ 0 + time + Conditions, 
                  data=PLN_data,ranks = 1:20)







########### ESSAI MARION #########################

##### Preparation data ######

load('C:/Users/33638/Documents/stage/MICROBIOTA_SPATIAL_FLIGHT/DATA_PROJECT_1.RData')

# Matrice d'abondance_________________________
# Certains OTU ne sont observés pour aucun individu, on les supprime
nb <- colSums(matrix_otu)
matrix_otu_D0 <- as.data.frame(matrix_otu[seq(1,28,by=2),nb!=0])  #10 OTU supprimés, on garde juste D0

# Matrice de covariables______________________

## On récupère les variables V02max et vO2max_rel (à D0)
morph_data<- as.data.frame(morphological_data[morphological_data[,2]=='D0',][c(-5,-8,-12,-17),c(1,3,9)])
morph_data$VO2max <- as.numeric(morph_data$VO2max)
morph_data$Weight <- as.numeric(morph_data$Weight)
morph_data <- data.frame(morph_data,VO2max_rel=morph_data$VO2max/morph_data$Weight)
row.names(morph_data) <- morph_data$ID
morph_data <- morph_data[,c(-1,-3)]  

# - 3 : on supprime aussi le poids, on le prends dans tableau lM
# attention on garde donc le poids (+précis) de lM mais pas D0, c'est J-4 avant début expé
# Rq : on a bien divisé la VO2max_rel par le poids au moment de la mesure à D0


# on récupère les autres variables 
lM_OK <- lM[seq(2,37,by=2),]  #on garde J-4 (env D0)
lM_OK <- lM_OK[c(-5,-8,-12,-17),c(4,7)] #on supprime les indiv pour lesquels on a pas les OTU
morph_data <- data.frame(morph_data,Weight=as.numeric(lM_OK$Whole.Body.Total.mass),Lean_mass=as.numeric(lM_OK$Whole.Body.Lean.mass))
morph_data <- data.frame(morph_data,Lean_mass_rel=morph_data$Lean_mass/morph_data$Weight)

#on met sous la bonne forme pour les fct du package PLN
                                                #"A" "B" "C" "D" "F" "G" "I" "J" "K" "N" "O" "P" "Q" "S"
row.names(morph_data)=row.names(matrix_otu_D0)=lM$X[c(2:9,12:15,18:23,26:33,36:37)][seq(1,28,by=2)]
data_PLN <- list(matrix_otu_D0=matrix_otu_D0,morph_data=morph_data)
data_PLN <- prepare_data(counts = data_PLN$matrix_otu_D0, covariates = data_PLN$morph_data)


######### Modélisation PLN ###############

#sans covariable ni offset
myPLN <- PLN(Abundance ~ 1, data_PLN)

fitted   = as.vector(fitted(myPLN))
observed = as.vector(data_PLN$Abundance)

plot(observed,fitted) 
summary(lm(fitted~observed))    #trop bien bizarre..
length(observed)
sum(observed==fitted)

library(corrplot)
corrplot(sigma(myPLN),is.corr=FALSE, tl.cex = .5)


# sans covariable et avec offset : TSS
myPLN_offsets <- 
  PLN(Abundance ~ 1 + offset(log(Offset)),data = data_PLN)
fitted   = as.vector(fitted(myPLN_offsets))
observed = as.vector(data_PLN$Abundance)

plot(observed,fitted) 
summary(lm(fitted~observed))    
length(observed)
sum(observed==fitted)

rbind(
  myPLN$criteria,
  myPLN_offsets$criteria
) %>% knitr::kable()


# Avec covariable weight
myPLN_weight <- 
  PLN(Abundance ~ 1 + Weight, data = data_PLN)
fitted   = as.vector(fitted(myPLN_weight))
observed = as.vector(data_PLN$Abundance)

plot(observed,fitted) 
summary(lm(fitted~observed))    

rbind(
  myPLN$criteria,
  myPLN_weight$criteria
) %>% knitr::kable()


# Avec covariable tree
myPLN_weight <- 
  PLN(Abundance ~ 1 + treefile, data = data_PLN)
fitted   = as.vector(fitted(myPLN_weight))
observed = as.vector(data_PLN$Abundance)

plot(observed,fitted) 
summary(lm(fitted~observed))    

rbind(
  myPLN$criteria,
  myPLN_weight$criteria
) %>% knitr::kable()
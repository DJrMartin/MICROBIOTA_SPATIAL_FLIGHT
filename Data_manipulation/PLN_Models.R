###Utilisation du package PLNmodels
#https://oliviergimenez.github.io/code_livre_variables_cachees/chiquet.html

library('PLNmodels')
library(corrplot)
library(dplyr)

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

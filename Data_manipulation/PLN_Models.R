###Utilisation du package PLNmodels
#https://oliviergimenez.github.io/code_livre_variables_cachees/chiquet.html

library('PLNmodels')

data(trichoptera)
trichoptera <- prepare_data(trichoptera$Abundance, trichoptera$Covariate)

myPLN <- PLN(Abundance ~ 1, data = trichoptera)
coef(myPLN)
myPCA <- PLNPCA(Abundance ~ 1, data = trichoptera, ranks = 1:8)


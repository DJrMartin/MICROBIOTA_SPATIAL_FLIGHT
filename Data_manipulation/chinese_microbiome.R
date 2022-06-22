library('PLNmodels')
library(corrplot)
library(dplyr)
library(tidyverse)

# Individus en colonnes, OTU en lignes
matrix_otu <- read.table("chinese_microbiome/clean_data.txt")

meta <- read.csv("chinese_microbiome/clean_meta.txt",na.strings=c(" "),sep="\t",row.names=1)
taxo <- read.csv("chinese_microbiome/clean_tax.txt",na.strings=c(" "),sep="\t")


# on prépare les données pour PLN
matrix_otu <- as.data.frame(t(matrix_otu))
row.names(meta) <- row.names(matrix_otu)
data_PLN <- list(matrix_otu=matrix_otu,morph_data=meta)
data_PLN <- prepare_data(counts = data_PLN$matrix_otu, covariates = data_PLN$morph_data)


######### Modélisation PLN ###############

#sans covariable ni offset-----------------
myPLN <- PLN(Abundance ~ 1, data_PLN)

fitted   = as.vector(fitted(myPLN))
observed = as.vector(data_PLN$Abundance)

plot(observed,fitted) 
summary(lm(fitted~observed))    #trop de paramètres
myPLN$criteria

# Sigma
heatmap(myPLN$model_par$Sigma,Colv = NA, Rowv = NA)

# sans covariable et avec offset : TSS--------------
myPLN_offsets <- PLN(Abundance ~ 1 + offset(log(Offset)),data = data_PLN)
fitted   = as.vector(fitted(myPLN_offsets))
observed = as.vector(data_PLN$Abundance)

plot(observed,fitted) 
summary(lm(fitted~observed))    


rbind(
  myPLN$criteria,
  myPLN_offsets$criteria
) %>% knitr::kable()


# Avec covariable weight-------------------
myPLN_weight <- PLN(Abundance ~ 1 + Weight, data = data_PLN)

fitted   = as.vector(fitted(myPLN_weight))
observed = as.vector(data_PLN$Abundance)

plot(observed,fitted) 
summary(lm(fitted~observed))    

rbind(
  myPLN$criteria,
  myPLN_weight$criteria
) %>% knitr::kable()


# Theta^
heatmap(myPLN_weight$model_par$Theta,Colv = NA, Rowv = NA)
# Sigma^
heatmap(myPLN_weight$model_par$Sigma,Colv = NA, Rowv = NA)
#MU^
un <- rep(1,14)
X <- matrix(c(un,data_PLN$Weight),nrow=14)
mu_chap <- X%*%t(myPLN_weight$model_par$Theta)+(un%*%t(diag(myPLN_weight$model_par$Sigma)))/2
# on regarde les corrélations entre OTU dans cette matrice mu_chap
heatmap(cor(mu_chap),Colv = NA, Rowv = NA)
heatmap(cor(matrix_otu_D0),Colv = NA, Rowv = NA)




### On force la matrice sigma comme étant diagonale
myPLN_diag <- PLN(Abundance ~ 1, data_PLN,control = list(covariance = "diagonal"))

fitted   = as.vector(fitted(myPLN_diag))
observed = as.vector(data_PLN$Abundance)

plot(observed,fitted) 
summary(lm(fitted~observed))    

rbind(
  myPLN$criteria,
  myPLN_diag$criteria
) %>% knitr::kable()



###### PCA ################
## Modèle sans covariable
myPCA_m0 <- PLNPCA(formula = Abundance ~ 1 + offset(log(Offset)),
                   data = data_PLN, 
                   ranks = 1:3) 
myPCA_m0
plot(myPCA_m0, reverse = TRUE)
#on extrait meilleur modèle (selon BIC)
PCA_m0_BIC <- getBestModel(myPCA_m0, "BIC")
PCA_m0_BIC                       

#le graphe des individus : très intéssant, mais il faudrait comprendre commenton arrive à ca!!!!
factoextra::fviz_pca_ind(PCA_m0_BIC,col.ind = data_PLN$Sex,repel=T)
factoextra::fviz_pca_ind(PCA_m0_BIC,col.ind = data_PLN$Age,repel=T,gradient.cols = c("yellow","red"))
factoextra::fviz_pca_ind(PCA_m0_BIC,col.ind = data_PLN$Group,repel=T)
factoextra::fviz_pca_ind(PCA_m0_BIC,col.ind = as.numeric(data_PLN$Weight_kg),repel=T,gradient.cols = c("yellow","red"))




######################################
###### REG PENALISEE #################
######################################

library(glmnet)

y_weight <- as.numeric(meta$Weight_kg)
ind_NA_weight <- which(is.na(y_weight))
y_weight <- y_weight[-ind_NA_weight]
x_weight <- as.matrix(matrix_otu[-ind_NA_weight,])


###### Weight ------------------------------
#on sépare le jeu de données
train_rows <- sample(1:nrow(x_weight), .66*nrow(x_weight))
x.train <- x_weight[train_rows, ]
x.test <- x_weight[-train_rows, ]
y.train <- y_weight[train_rows]
y.test <- y_weight[-train_rows]

# On utilise une validation croisée 10 blocs pour déterminer la valeur optimale pour lambda

## D'abord avec alpha = 0, Ridge Regression
################################
alpha0.fit <- cv.glmnet(x.train, y.train, type.measure="mse", 
                        alpha=0, family="gaussian")

best_lambda <- alpha0.fit$lambda.min
best_lambda

#produce plot of test MSE by lambda value
plot(alpha0.fit) 

#find coefficients of best model
best_model <- glmnet(x.train, y.train, alpha = 0, lambda = best_lambda)

# Trace plot to visualize how the coefficient estimates changed as a result of increasing lambda:
model <- glmnet(x.train, y.train, alpha = 0)
plot(model, xvar = "lambda")

#use fitted best model to make predictions
alpha0.predicted <- 
  predict(alpha0.fit, s=best_lambda, newx=x.test)

plot(y.test,alpha0.predicted)
#find SST and SSE
sst <- sum((y.test - mean(y.test))^2)
sse <- sum((alpha0.predicted - y.test)^2)

#find R-Squared
rsq <- 1 - sse/sst
rsq


## Puis avec alpha = 1, Lasso Regression
################################
alpha1.fit <- cv.glmnet(x.train, y.train, type.measure="mse", 
                        alpha=1, family="gaussian")

best_lambda <- alpha1.fit$lambda.min
best_lambda

#produce plot of test MSE by lambda value
plot(alpha1.fit) 

#find coefficients of best model
best_model <- glmnet(x.train, y.train, alpha = 1, lambda = best_lambda)
coeff <- coef(best_model)


#use fitted best model to make predictions
alpha1.predicted <- 
  predict(alpha1.fit, s=best_lambda, newx=x.test)

plot(y.test,alpha1.predicted)
#find SST and SSE 
sst <- sum((y.test - mean(y.test))^2)
sse <- sum((alpha1.predicted - y.test)^2)

#find R-Squared
rsq <- 1 - sse/sst
rsq


###### Age ------------------------------
y_age <- as.numeric(meta$Age)
x_age <- as.matrix(matrix_otu)


#on sépare le jeu de données
train_rows <- sample(1:nrow(x_age), .66*nrow(x_age))
x.train <- x_age[train_rows, ]
x.test <- x_age[-train_rows, ]
y.train <- y_age[train_rows]
y.test <- y_age[-train_rows]

# On utilise une validation croisée 10 blocs pour déterminer la valeur optimale pour lambda

## D'abord avec alpha = 0, Ridge Regression
################################
alpha0.fit <- cv.glmnet(x.train, y.train, type.measure="mse", 
                        alpha=0, family="gaussian")

best_lambda <- alpha0.fit$lambda.min
best_lambda

#produce plot of test MSE by lambda value
plot(alpha0.fit) 

#find coefficients of best model
best_model <- glmnet(x.train, y.train, alpha = 0, lambda = best_lambda)

# Trace plot to visualize how the coefficient estimates changed as a result of increasing lambda:
model <- glmnet(x.train, y.train, alpha = 0)
plot(model, xvar = "lambda")

#use fitted best model to make predictions
alpha0.predicted <- 
  predict(alpha0.fit, s=best_lambda, newx=x.test)

plot(y.test,alpha0.predicted)
#find SST and SSE
sst <- sum((y.test - mean(y.test))^2)
sse <- sum((alpha0.predicted - y.test)^2)

#find R-Squared
rsq <- 1 - sse/sst
rsq


## Puis avec alpha = 1, Lasso Regression
################################
alpha1.fit <- cv.glmnet(x.train, y.train, type.measure="mse", 
                        alpha=1, family="gaussian")

best_lambda <- alpha1.fit$lambda.min
best_lambda

#produce plot of test MSE by lambda value
plot(alpha1.fit) 

#find coefficients of best model
best_model <- glmnet(x.train, y.train, alpha = 1, lambda = best_lambda)
coeff <- coef(best_model)


#use fitted best model to make predictions
alpha1.predicted <- 
  predict(alpha1.fit, s=best_lambda, newx=x.test)

plot(y.test,alpha1.predicted)
#find SST and SSE 
sst <- sum((y.test - mean(y.test))^2)
sse <- sum((alpha1.predicted - y.test)^2)

#find R-Squared
rsq <- 1 - sse/sst
rsq


######## GLOBALEMENT LASSO meilleur pour prédire







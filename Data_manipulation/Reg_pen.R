library(glmnet)

############################
##### Données spatiales #####
############################

load('C:/Users/33638/Documents/stage/MICROBIOTA_SPATIAL_FLIGHT/DATA_PROJECT_1.RData')

# Matrice d'abondance_________________________
# Certains OTU ne sont observés pour aucun individu, on les supprime
nb <- colSums(matrix_otu)
x <- as.matrix(matrix_otu[seq(1,28,by=2),nb!=0])  #10 OTU supprimés, on garde juste D0

# Matrice de covariables______________________

## On récupère les variables à prédire
morph_data_D0 <- as.data.frame(morphological_data[morphological_data[,2]=='D0',][c(-5,-8,-12,-17),c(1,3,9)])
morph_data_D5 <- as.data.frame(morphological_data[morphological_data[,2]=='D5',][c(-5,-8,-12,-17),c(1,3,9)])
morph_data_D0$VO2max <- as.numeric(morph_data_D0$VO2max)
morph_data_D0$Weight <- as.numeric(morph_data_D0$Weight)
morph_data_D0 <- data.frame(morph_data_D0,VO2max_rel=morph_data_D0$VO2max/morph_data_D0$Weight)
morph_data_D5$VO2max <- as.numeric(morph_data_D5$VO2max)
morph_data_D5$Weight <- as.numeric(morph_data_D5$Weight)
morph_data_D5 <- data.frame(morph_data_D5,VO2max_rel=morph_data_D5$VO2max/morph_data_D5$Weight)

morph_data_delta <- data.frame(delta_weight=morph_data_D0$Weight-morph_data_D5$Weight,delta_VO2max=morph_data_D0$VO2max-morph_data_D5$VO2max,delta_VO2max_rel=morph_data_D0$VO2max_rel-morph_data_D5$VO2max_rel)

row.names(morph_data_D0)=row.names(morph_data_D5)=row.names(morph_data_delta)=morph_data_D0$ID
morph_data_D0 <- morph_data_D0[,-1]  
morph_data_D5 <- morph_data_D5[,-1] 

# on récupère les autres variables 
lM_OK_D0 <- lM[seq(2,37,by=2),]  #on garde J-4 (env D0)
lM_OK_D0 <- lM_OK_D0[c(-5,-8,-12,-17),c(1,2,4,7)] #on supprime les indiv pour lesquels on a pas les OTU
lM_OK_D5 <- lM[seq(3,37,by=2),]  
lM_OK_D5 <- lM_OK_D5[c(-5,-8,-12,-17),c(1,2,4,7)] #on supprime les indiv pour lesquels on a pas les OTU

## tableau final des y possibles
morph_data <- data.frame(Weight=as.numeric(lM_OK_D0$Whole.Body.Total.mass),
                         delta_Weight=as.numeric(lM_OK_D0$Whole.Body.Total.mass)-as.numeric(lM_OK_D5$Whole.Body.Total.mass),
                         VO2max=morph_data_D0$VO2max,
                         delta_VO2max=morph_data_delta$delta_VO2max,
                         VO2max_rel=morph_data_D0$VO2max_rel,
                         delta_VO2max_rel=morph_data_delta$delta_VO2max_rel,
                         Lean_mass=as.numeric(lM_OK_D0$Whole.Body.Lean.mass),
                         delta_Lean_mass=as.numeric(lM_OK_D0$Whole.Body.Lean.mass)-as.numeric(lM_OK_D5$Whole.Body.Lean.mass))
morph_data <- data.frame(morph_data,Lean_mass_rel=morph_data$Lean_mass/morph_data$Weight,
                            delta_Lean_mass_rel=morph_data$delta_Lean_mass/morph_data$Weight)


# Variable y à prédire:
y <- morph_data$delta_Weight

#on sépare le jeu de données
train_rows <- sample(1:nrow(x), .75*nrow(x))
x.train <- x[train_rows, ]
x.test <- x[-train_rows, ]
y.train <- y[train_rows]
y.test <- y[-train_rows]

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
length(coef(best_model))  #aucun coeff supprimé car Ridge


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
  
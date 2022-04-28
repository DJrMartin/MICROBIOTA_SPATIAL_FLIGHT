##Pretraitement du metalome et gestion du fichier pour le merger aux donnees microbiote.
library(dplyr)
###Import#####################################################################
User='David'
if (User=='David'){
  load('C:/Users/David/Dropbox/EXOMIC 2022/DATA_PROJECT_1.RData')
}
if (User=='Marion'){
  load()
}
rm(User)

metalome=data.frame(metalome)
metalome$ID=paste(substr(as.character(metalome$ID),2,3), rep(c(1,2), each=1))
experimental_condition$ID=paste(as.character(experimental_condition$subject), rep(c(1,2), each=1))

data.merged=merge(metalome, experimental_condition, by='ID')

metals=data.merged%>%
  select(Cuivre,Zinc, Chrome,
         Fer, Se,Vanadium,
         Manganese)
  
colnames(metals)=c('Cu','Zn','Ch','Fe','Se','Va','Mang')
metals=as.matrix(gsub(',', '.', apply(metals, 2, as.character)))
data_normalised=data.frame(scale(apply(metals, 2, as.numeric)))

res.pca=FactoMineR::PCA(scale(data_normalised))
plot(res.pca$ind$coord[which(data.merged$time.y=='D0'), c(1,2)], 
     ylim=c(range(res.pca$ind$coord[,2])), xlim=c(range(res.pca$ind$coord[,1])))
points(res.pca$ind$coord[which(data.merged$time.y=='D5'), c(1,2)], col='red')
legend('topleft', legend=c('Before', 'After'), fill=c('black', 'red'), lty=0, bty='n')
layout(matrix(c(1,2,3,4,
                5,6,7,8), nrow = 2))
for (var in 1:7){
  plot(data_normalised[,var]~data.merged$time.y, 
       ylab=colnames(data_normalised)[var], 
       xlab='Time')
}
pairs(data_normalised)

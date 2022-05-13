library(vegan)
library(fossil)
library(corrplot)
library(ggplot2)
library(GGally)

###Import#####################################################################
User='David'
if (User=='David'){
  load( "/Users/martind/Dropbox/Exomic 2022/RData_Microgravity_updata.RData")
}
if (User=='Marion'){
  load()
}
rm(User)

ind=matrix_otu[2,which(matrix_otu[2,]>0)] #modifier la borne, l'individu..
rare_curve=NULL

for (i in 1:length(ind)){
  new=rep(rownames(data.frame(ind))[i], ind[i])
  rare_curve=c(rare_curve,new)
}

##Expression des indices de diversite en fonction de la profondeur de séquençage.
rc=rshannon=rchao=NULL
for (j in seq(1,length(rare_curve), by=100)){
  u=sample(as.factor(rare_curve), size=j)
  rc=c(rc,length(unique(u)))
  rshannon=c(rshannon,diversity(table(u), index='shannon'))
  rchao=c(rchao,chao1(table(u)))
}
rare_rich=data.frame("deep"=seq(1,length(rare_curve), by=100),"Richesse"=rc)
rare_chao=data.frame("deep"=seq(1,length(rare_curve), by=100),"Richesse"=rchao)

range(rare_rich$Richesse, rare_chao$Richesse)
plot(rare_rich, ylim=c(range(rare_rich$Richesse, rare_chao$Richesse)))
points(rare_chao, col='red')
sd(tail(rc, 10))

plot(x=seq(1,length(rare_curve), by=100) , rshannon, xlab='profondeur de séquençage',
     ylab='Entropy')

##Selection d'variables observé.

entropy_variability=all_matrice=NULL
cnt=1
while(cnt<100){
  entropy=NULL
  for (j in seq(1,length(ind), by=1)){
    u=sample(as.factor(ind), size=j)
    entropy=c(entropy,diversity(as.numeric(as.character(u)), index='shannon'))
  }
  entropy_variability=data.frame(cbind(entropy_variability,entropy))
  cnt=cnt+1
}

boxplot(t(entropy_variability[,-1]), xlab="nb d'èspeces", ylab='Entropy',
        main = "Variabilité de l'entropy (300 tirages)")

#####Visualisation de l'entropy par niveau de phylogénie.
ENTROPY=NULL
otu=NULL
filtre=species[,6]
for (j in as.character(unique(filtre))){
  sum_reads=apply(OTU_normalised[, which(filtre==j)], 1,sum)
  otu=cbind(otu,sum_reads)
}

ENTROPY=cbind(ENTROPY, diversity(matrix_otu, index='shannon'))
colnames(ENTROPY)=colnames(species)[2:7]

###Visualisation de l'entropy intra-class.
ENTROPY=NULL
otu=NULL
filtre=species[,2]
for (j in as.character(unique(filtre))){
  reads=OTU_normalised[, which(filtre==j)]
  ENTROPY=cbind(ENTROPY, diversity(reads, index='shannon'))
}
colnames(ENTROPY)=colnames(as.character(unique(filtre)))
layout(matrix(c(1,2,3,4), nrow = 2))
boxplot(ENTROPY[,1]~experimental_condition$time, xlab='Entropy', ylab='', main=as.character(unique(filtre))[1])
boxplot(ENTROPY[,2]~experimental_condition$time, xlab='Entropy', ylab='', main=as.character(unique(filtre))[2])
plot(ENTROPY[c(27,28),1]~experimental_condition$time[c(1:2)], xlab='Entropy', ylab='', main=as.character(unique(filtre))[3], pch=1)
boxplot(ENTROPY[,4]~experimental_condition$time, xlab='Entropy', ylab='', main=as.character(unique(filtre))[4])

rownames(ENTROPY)=rownames(matrix_otu)
ENTROPY=cbind(ENTROPY, experimental_condition$time)

dotplot(ENTROPY[,1]~experimental_condition$time, col=as.numeric(experimental_condition$subject))
new=cbind(ENTROPY[which(experimental_condition$time=="D0"),1], ENTROPY[which(experimental_condition$time=="D5"),1])
ggparcoord(new, scale = "std", title = "Standard Scale")

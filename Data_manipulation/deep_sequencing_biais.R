library(vegan)
library(fossil)
library(corrplot)
###Import#####################################################################
User='David'
if (User=='David'){
  load( "/Users/martind/Desktop/RData_Microgravity_updata.RData")
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

##Choix au hazard de 200 bactéries 100x pour étudier la variabilité des indices de diversité
set.seed(123)
cnt=1
rsimpson=rshannon=rchao=rc=NULL
while (cnt<=200){
  u=sample(as.factor(rare_curve), size=200)
  rc=c(rc,length(unique(u)))
  rsimpson=c(rsimpson,diversity(table(u), index='simpson'))
  rshannon=c(rshannon,diversity(table(u), index='shannon'))
  rchao=c(rchao,chao1(table(u)))
  cnt=cnt+1
}
layout(matrix(c(1,2,3,4), nrow = 1))
boxplot(rc, main='richesse')
boxplot(rsimpson, main='Simpson') #mean est identique à la valeur réelle (0.95)
boxplot(rshannon, main = "Shannon") #mean est inférieur à la valeur réelle (3.8)
boxplot(rchao, main = "Chao")
corrplot(cor(data.frame('Chao'=rchao, 'richesse'=rc, 'Simspon'=rsimpson, 'Shannon'=rshannon)),type='upper')

##distribution pour les extrémités de shannon
which(rshannon==max(rshannon))
plot(sort(table(u)[which(table(u)>0)], decreasing = T))
which(rc==min(rc))
plot(sort(table(u)[which(table(u)>0)], decreasing = T))



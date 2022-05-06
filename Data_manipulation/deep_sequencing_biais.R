library(vegan)
library(fossil)
###Import#####################################################################
User='David'
if (User=='David'){
  load("/Users/martind/Dropbox/EXOMIC 2022/DATA_PROJECT_1.RData")
}
if (User=='Marion'){
  load()
}
rm(User)
layout(matrix(c(1,2,3), nrow=3))
ind=matrix_otu[2,which(matrix_otu[2,]>0)] #modifier la borne, l'individu..
rare_curve=NULL

for (i in 1:length(ind)){
  new=rep(rownames(data.frame(ind))[i], ind[i])
  rare_curve=c(rare_curve,new)
}

rc=rshannon=rchao=NULL
for (j in seq(1,length(rare_curve), by=100)){
  u=sample(as.factor(rare_curve), size=j)
  rc=c(rc,length(unique(u)))
  rshannon=c(rshannon,diversity(table(u), index='simpson'))
  rchao=c(rchao,chao1(table(u)))
}


rare_rich=data.frame("deep"=seq(1,length(rare_curve), by=100),"Richesse"=rc)
rare_chao=data.frame("deep"=seq(1,length(rare_curve), by=100),"Richesse"=rchao)

range(rare_rich$Richesse, rare_chao$Richesse)
plot(rare_rich, ylim=c(range(rare_rich$Richesse, rare_chao$Richesse)))
points(rare_chao, col='red')
sd(tail(rc, 10))

######déstructuration de la table de comptage pour voir la différence des reads.
new_data=as.numeric(sort(matrix_otu[1,which(matrix_otu[1,]>0)], decreasing = T))
mean(new_data)
plot(new_data)

cnt=1
while (cnt<300){
  inf=new_data<mean(new_data)
  new_data[inf==T]=as.numeric(new_data[inf==T])*2
  new_data[inf==F]=as.numeric(new_data[inf==F])/2
  points(new_data)
  cnt=cnt+1
}


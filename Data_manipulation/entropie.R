library(vegan)

load('C:/Users/33638/Documents/stage/MICROBIOTA_SPATIAL_FLIGHT/DATA_PROJECT_1.RData')

# On regarde pour un individu
vect_otu <- matrix_otu[1,]
vect_otu <- vect_otu[vect_otu!=0]

#vecteur des probabilités:
p <- vect_otu/sum(vect_otu)

#---Fct D'information I(p)---------------
# I(p)=-log(p) est la fonction d'information pour l'entropie/shannon :
# un OTU apparu peu de fois est associé à une information importante, Ip fct décroissante de p
Ip=-log(p)
# on représente Ip en fonction de p:
plot(x=p,y=Ip, main="Fonction d'information pour le premier individu",xlim=c(0,1),ylim=c(0,10))
     #Beaucoup de Proba faibles donc associées à beaucoup d'informations


#----Qt d'information H(p)-----------------
# Quantité d'information : H(p)=somme(p*I(p))
#    H(p) = Shannon si I(p)=-log(p) 
# comme on multiplie par p, qui est un vecteur de proba (p_s entre [0,1])
# on a H(p) qui est :
#    - croissant tant que I(p)=-log(p) est supérieur à 1
#    - décroissant si I(p)=-log(p) est inférieur à 1

# on représente les contributions (absolues) de chaque bactérie 
# à la valeur totale de l'entropie
contrib<-p*Ip
plot(x=p,y=contrib, main="Contribution de chaque espèce à la valeur totale \n de l'entropie pour le premier individu")


#l'indice d'entropie (=shannon) est la somme des termes/contributions -p*log(p)
shannon <- sum(p*Ip)
#ou
shannon <- diversity(vect_otu, index='shannon')


#on peut donc représenter à nouveau les contributions 
# mais en  pourcentage de l'entropie totale (contributions relatives)
contrib_rel<-contrib/shannon
plot(x=p,y=contrib_rel, main="Contribution relative de chaque espèce \n à la valeur totale de l'entropie")



#### On fait varier les choses#####


## ----- Sans les 1:5 ----------
vect_otu <- vect_otu[vect_otu>5]
#vecteur des probabilités:
p <- vect_otu/sum(vect_otu)
# information
Ip=-log(p)
# on représente Ip en fonction de p:
plot(x=p,y=Ip, main="Fonction d'information")
#Beaucoup de Proba faibles donc associées à beaucoup d'informations

# on représente les contributions (absolues) de chaque bactérie 
# à la valeur totale de l'entropie
contrib<-p*Ip
plot(x=p,y=contrib, main="Contribution absolue de chaque espèce \n à la valeur totale de l'entropie")

#l'indice d'entropie (=shannon) est la somme des termes/contributions -p*log(p)
shannon <- sum(p*Ip)

#on peut donc représenter à nouveau les contributions 
# mais en  pourcentage de l'entropie totale (contributions relatives)
contrib_rel<-contrib/shannon
plot(x=p,y=contrib_rel, main="Contribution relative de chaque espèce \n à la valeur totale de l'entropie")





###### Analyse sur le jeu de données des oiseaux ########

data(BCI)
vect_bci <- BCI[1,]
vect_bci <- vect_bci[vect_bci!=0]
#vecteur des probabilités:
p <- vect_bci/sum(vect_bci)
# information
Ip=-log(p)
# on représente Ip en fonction de p:
plot(x=p,y=Ip, main="Fonction d'information")
#Beaucoup de Proba faibles donc associées à beaucoup d'informations

# on représente les contributions (absolues) de chaque bactérie 
# à la valeur totale de l'entropie
contrib<-p*Ip
plot(x=p,y=contrib, main="Contribution absolue de chaque espèce \n à la valeur totale de l'entropie")

#l'indice d'entropie (=shannon) est la somme des termes/contributions -p*log(p)
shannon <- sum(p*Ip)

#on peut donc représenter à nouveau les contributions 
# mais en  pourcentage de l'entropie totale (contributions relatives)
contrib_rel<-contrib/shannon
plot(x=p,y=contrib_rel, main="Contribution relative de chaque espèce \n à la valeur totale de l'entropie")


hist(p,breaks=30)






  ###### Analyse textuelle (texte de Proust)###########
library("tm")
library("tidytext")
library("proustr")
library("tidyverse")
devtools::install_github("ThinkRstat/stopwords")
library("stopwords")
books <- proust_books()

# on récupère les mots avec leur nombre d'apparition
books_tidy <- proust_books() %>%
  mutate(text = stringr::str_replace_all(.$text, "'", " ")) %>% 
  unnest_tokens(word, text) %>%
  filter(!word %in% stopwords_iso$fr) %>%
  count(word, sort = TRUE) %>% 
  head(300)


#vecteur des probabilités:
p <- books_tidy$n/sum(books_tidy$n)
# information
Ip=-log(p)
# on représente Ip en fonction de p:
plot(x=p,y=Ip, main="Fonction d'information")
#Beaucoup de Proba faibles donc associées à beaucoup d'informations

# on représente les contributions (absolues) de chaque bactérie 
# à la valeur totale de l'entropie
contrib<-p*Ip
plot(x=p,y=contrib, main="Contribution absolue de chaque espèce \n à la valeur totale de l'entropie")

#l'indice d'entropie (=shannon) est la somme des termes/contributions -p*log(p)
shannon <- sum(p*Ip)

#on peut donc représenter à nouveau les contributions 
# mais en  pourcentage de l'entropie totale (contributions relatives)
contrib_rel<-contrib/shannon
plot(x=p,y=contrib_rel, main="Contribution relative de chaque espèce \n à la valeur totale de l'entropie")









########### Simulation ##############

# vecteur des entropies
entropies=c()
for(i in 2:200){
  p <- rep(1/i,i)
  Ip <- -log(p)
  contrib <- p*Ip
  entropie_i=sum(p*Ip)
  entropies=c(entropies,entropie_i)
}

plot(entropies,xlab="Nb espèces (toutes équiprobables)",ylab="Indice de Shannon (entropie)",
       main="Evolution de l'indice de Shannon (entropie) \n en fonction du nombre d'espèces")





### Le maximum de contribution à l'entropie est atteint en p tq log(p)=-1 => p=0.368

p=c(0.25,0.05,rep(0.02,3),0.04,0.6)
Ip=-log(p)
contrib<-p*Ip
plot(x=p,y=contrib, main="Simulation de la contribution \n de 5 espèces à la valeur totale de l'entropie")
abline(v=exp(-1),col="red")

###information
p=seq(0.00001,1,by=0.01)
y=-log(p)
plot(x=p,y=y,type="l",ylab="I(ps)=-log(ps)", main="Fonction d'information pour \n l'indice de Shannon")
segments(x0=-2,y0=1,x1=exp(-1),y1=1, col="red")
segments(x0=exp(-1),y0=-1,y1=1, col="red")
text(x=0.4,y=3,"p=0.368 <=> -log(p)=1")


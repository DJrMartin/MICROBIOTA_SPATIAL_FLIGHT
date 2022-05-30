library(vegan)
library(tidyverse)

load('C:/Users/33638/Documents/stage/MICROBIOTA_SPATIAL_FLIGHT/DATA_PROJECT_1.RData')

# On regarde pour un individu 
vect_otu <- matrix_otu["BG-I-D0-CTL",]
vect_otu <- vect_otu[vect_otu!=0]

#vecteur des probabilités:
p <- vect_otu/sum(vect_otu)

#---Fct D'information I(p)---------------
# I(p)=-log(p) est la fonction d'information pour l'entropie/shannon :
# un OTU apparu peu de fois est associé à une information importante, Ip fct décroissante de p
Ip=-log(p)
# on représente Ip en fonction de p:
plot(x=p,y=Ip, main="Fonction d'information pour l'individu \n avec Shannon minimal (3)",xlim=c(0,1),ylim=c(0,10))
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
plot(x=p,y=contrib,xlim=c(0,0.25),ylim=c(0,0.30), main="Contribution de chaque espèce à la valeur totale \n de l'entropie pour individu avec Shannon max (4.23)")
abline(v=mean(p),col="blue")


#l'indice d'entropie (=shannon) est la somme des termes/contributions -p*log(p)
shannon <- sum(p*Ip)
#ou
shannon <- diversity(vect_otu, index='shannon')


#on peut donc représenter à nouveau les contributions 
# mais en  pourcentage de l'entropie totale (contributions relatives)
contrib_rel<-contrib/shannon
plot(x=p,y=contrib_rel, main="Contribution relative de chaque espèce \n à la valeur totale de l'entropie")



###### SHANNON EN FONCTION DE LA SOMME D'INFORMATION
shannon_tot <- diversity(matrix_otu, index='shannon')
I_sum=c()
p_max=c()
p_mean=c()
richness=c()
prop_1=c()
prop_inf5=c()
for(i in 1:28){
  vect_otu <- matrix_otu[i,]
  vect_otu <- vect_otu[vect_otu!=0]
  richness=c(richness,length(vect_otu))
  p=vect_otu/sum(vect_otu)
  p_max=c(p_max,max(p))
  p_mean=c(p_mean,mean(p))
  I_sum=c(I_sum,sum(-log(p)))
  prop_1=c(prop_1,length(vect_otu[vect_otu==1])*100/length(vect_otu))
  prop_inf5=c(prop_inf5,length(vect_otu[vect_otu<=5])*100/length(vect_otu))
}
infos_indiv <- cbind(richness,shannon_tot,p_max,p_mean,prop_1,prop_inf5)



plot(shannon_tot~I_sum, main="SHannon en fonction de la somme d'information")
text(x=2200,y=3.3,"R²=29.7%",col="blue")
###### LA SOMME D'INFORMATION en fonction proba max
plot(I_sum~p_max,main="Somme des informations en fonction de \n la probabilité maximale dans l'échantillon")

###### SHANNON en fonction proba max
plot(shannon_tot~p_max,main="Entropie Shannon en fonction de \n la probabilité maximale dans l'échantillon")


###### LA SOMME D'INFORMATION en fonction proba moyenne
plot(I_sum~p_mean,main="Somme des informations en fonction de \n la probabilité moyenne dans l'échantillon")
text(x=0.0045,y=2200,"R²=96%",col="blue")

###### SHANNON en fonction proba moyenne
plot(shannon_tot~p_mean,main="Entropie Shannon en fonction de \n la probabilité moyenne dans l'échantillon")
text(x=0.0045,y=3.3,"R²=40.6%",col="blue")


###### Shannon en fct de richesse
plot(richness,shannon_tot,main="Corrélation positive entre Shannon et la richesse S",xlab="Richesse S",ylab="Indice de Shannon")
text(x=290,y=3.3,"R²=39%",col="blue")



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

#courbe théorique :
p=seq(0.0001,1,by=0.01)
contrib=-p*log(p)
plot(x=p,y=contrib, main="Contributions à l'entropie \n selon la probabilité p associée",type="l")
abline(v=exp(-1),col="red")


p=c(0.25,0.05,rep(0.02,3),0.04,0.6)
Ip=-log(p)
contrib<-p*Ip
plot(x=p,y=contrib, main="Simulation de la contribution \n de 5 espèces à la valeur totale de l'entropie")
abline(v=exp(-1),col="red")




###information
p=seq(0.00001,1,by=0.01)
y=-log(p)
plot(x=p,y=y,type="l",ylab="I(p)=-log(p)", main="Fonction d'information pour \n l'indice de Shannon")
segments(x0=-2,y0=1,x1=exp(-1),y1=1, col="red")
segments(x0=exp(-1),y0=-1,y1=1, col="red")
text(x=0.4,y=3,"p=0.368 <=> -log(p)=1")










###### Données réduites 
vect_otu <- sort(matrix_otu["BG-I-D0-CTL",])
vect_otu <- vect_otu[vect_otu!=0]
vect_otu <- vect_otu[260:279]

#vecteur des probabilités:
p <- vect_otu/sum(vect_otu)

#---Fct D'information I(p)---------------
# I(p)=-log(p) est la fonction d'information pour l'entropie/shannon :
# un OTU apparu peu de fois est associé à une information importante, Ip fct décroissante de p
Ip=-log(p)
# on représente Ip en fonction de p:
plot(x=p,y=Ip, main="Fonction d'information",xlim=c(0,1),ylim=c(0,10))
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
plot(x=p,y=contrib,xlim=c(0,0.55),ylim=c(0,0.5), main="Contribution de chaque espèce à la valeur totale \n de l'entropie (3)")
abline(v=exp(-1),col="red")


#l'indice d'entropie (=shannon) est la somme des termes/contributions -p*log(p)
shannon <- sum(p*Ip)
#ou
shannon <- diversity(vect_otu, index='shannon')





###############################
#### REGROUPEMENT : pour indiv shannon max ###########
##############################
t_matrix_otu_species <- cbind(t(matrix_otu),species)
t_matrix_otu_species <- tibble::rownames_to_column(t_matrix_otu_species,"OTU")
t_matrix_otu_species <- tibble(t_matrix_otu_species)
summary(t_matrix_otu_species)

#### Phylum-----------
Phylum_groupes <- t_matrix_otu_species %>% 
  select(c("BG-I-D0-CTL" ,31))  %>%   #séléctionne la colonne=l'individu, et la colonne de regroupement (ici family)
  rename(c(Individu="BG-I-D0-CTL",Groupe=Phylum)) %>% 
  group_by(Groupe) %>% 
  summarize(effectif_moyen=sum(Individu))
vect_otu <- ceiling(as.data.frame(Phylum_groupes)$effectif_moyen)
names(vect_otu) <- as.data.frame(Phylum_groupes)$Groupe

#vecteur des probabilités:
p <- vect_otu/sum(vect_otu)
#Fct D'information I(p)
Ip=-log(p)
#Qt d'information H(p)
contrib<-p*Ip
plot(x=p,y=contrib,xlim=c(0,0.55),ylim=c(0,0.5), main="Contribution de chaque groupe de Phylum \n à la valeur totale de l'entropie ")
abline(v=exp(-1),col="red")
#l'indice d'entropie (=shannon) est la somme des termes/contributions -p*log(p)
shannon_Phylum <- sum(p*Ip)

#### Class-----------
Class_groupes <- t_matrix_otu_species %>% 
  select(c("BG-I-D0-CTL" ,32))  %>%   #séléctionne la colonne=l'individu, et la colonne de regroupement (ici family)
  rename(c(Individu="BG-I-D0-CTL",Groupe=Class)) %>% 
  group_by(Groupe) %>% 
  summarize(effectif_moyen=sum(Individu))

vect_otu <- ceiling(as.data.frame(Class_groupes)$effectif_moyen)
names(vect_otu) <- as.data.frame(Class_groupes)$Groupe
vect_otu <- vect_otu[vect_otu!=0]

#vecteur des probabilités:
p <- vect_otu/sum(vect_otu)
#Fct D'information I(p)
Ip=-log(p)
#Qt d'information H(p)
contrib<-p*Ip
plot(x=p,y=contrib,xlim=c(0,0.55),ylim=c(0,0.5), main="Contribution de chaque groupe de Phylum \n à la valeur totale de l'entropie ")
abline(v=exp(-1),col="red")
#l'indice d'entropie (=shannon) est la somme des termes/contributions -p*log(p)
shannon_Class <- sum(p*Ip)

#### Order-----------
Order_groupes <- t_matrix_otu_species %>% 
  select(c("BG-I-D0-CTL" ,33))  %>%   #séléctionne la colonne=l'individu, et la colonne de regroupement (ici order)
  rename(c(Individu="BG-I-D0-CTL",Groupe=Order)) %>% 
  group_by(Groupe) %>% 
  summarize(effectif_moyen=sum(Individu))

vect_otu <- ceiling(as.data.frame(Order_groupes)$effectif_moyen)
names(vect_otu) <- as.data.frame(Order_groupes)$Groupe
vect_otu <- vect_otu[vect_otu!=0]

#vecteur des probabilités:
p <- vect_otu/sum(vect_otu)
#Fct D'information I(p)
Ip=-log(p)
#Qt d'information H(p)
contrib<-p*Ip
plot(x=p,y=contrib,xlim=c(0,0.55),ylim=c(0,0.5), main="Contribution de chaque groupe de Order \n à la valeur totale de l'entropie ")
abline(v=exp(-1),col="red")

#l'indice d'entropie (=shannon) est la somme des termes/contributions -p*log(p)
shannon_Order <- sum(p*Ip)


#### Family-----------
Family_groupes <- t_matrix_otu_species %>% 
  select(c("BG-I-D0-CTL" ,34))  %>%   #séléctionne la colonne=l'individu, et la colonne de regroupement (ici family)
  rename(c(Individu="BG-I-D0-CTL",Groupe=Family)) %>% 
  group_by(Groupe) %>% 
  summarize(effectif_moyen=sum(Individu))

vect_otu <- ceiling(as.data.frame(Family_groupes)$effectif_moyen)
names(vect_otu) <- as.data.frame(Family_groupes)$Groupe
vect_otu <- vect_otu[vect_otu!=0]

#vecteur des probabilités:
p <- vect_otu/sum(vect_otu)
#Fct D'information I(p)
Ip=-log(p)
#Qt d'information H(p)
contrib<-p*Ip
plot(x=p,y=contrib,xlim=c(0,0.55),ylim=c(0,0.5), main="Contribution de chaque groupe de Family \n à la valeur totale de l'entropie ")
abline(v=exp(-1),col="red")
#l'indice d'entropie (=shannon) est la somme des termes/contributions -p*log(p)
shannon_Family <- sum(p*Ip)


#### Genus-----------
Genus_groupes <- t_matrix_otu_species %>% 
  select(c("BG-I-D0-CTL" ,35))  %>%   #séléctionne la colonne=l'individu, et la colonne de regroupement 
  rename(c(Individu="BG-I-D0-CTL",Groupe=Genus)) %>% 
  group_by(Groupe) %>% 
  summarize(effectif_moyen=sum(Individu))

vect_otu <- ceiling(as.data.frame(Genus_groupes)$effectif_moyen)
names(vect_otu) <- as.data.frame(Genus_groupes)$Groupe
vect_otu <- vect_otu[vect_otu!=0]

#vecteur des probabilités:
p <- vect_otu/sum(vect_otu)
#Fct D'information I(p)
Ip=-log(p)
#Qt d'information H(p)
contrib<-p*Ip
plot(x=p,y=contrib,xlim=c(0,0.55),ylim=c(0,0.5), main="Contribution de chaque groupe de Genus \n à la valeur totale de l'entropie")
abline(v=exp(-1),col="red")
#l'indice d'entropie (=shannon) est la somme des termes/contributions -p*log(p)
shannon_Genus <- sum(p*Ip)


########################----------
barplot(c(shannon_Phylum,shannon_Class,shannon_Order,shannon_Family,shannon_Genus,4.2342),names.arg=c("Phylum","Class","Order","Family","Genus","Ensemble"),main="Evolution de l'indice de Shannon (pour individu avec Shannon max.) \n selon les regroupements (du plus au moins fin)")







###############################
#### REGROUPEMENT : pour tous les indiv ###########
##############################
# on regroupe les OTU pour avoir moins de variables
t_matrix_otu_species <- cbind(t(matrix_otu),species)
t_matrix_otu_species <- tibble::rownames_to_column(t_matrix_otu_species,"OTU")
t_matrix_otu_species <- tibble(t_matrix_otu_species)
summary(t_matrix_otu_species)


### Phylum
res_phylum=levels(t_matrix_otu_species$Phylum)
for(col in 2:29){
  res=aggregate(t_matrix_otu_species[,col],t_matrix_otu_species[,"Phylum"],sum)
  res_phylum=cbind(res_phylum,res[,2])
}
row.names(res_phylum)<-res_phylum[,1]
res_phylum <- apply(as.data.frame(res_phylum[,2:29]),1,FUN=as.integer)

shannon_tot_phylum <- diversity(res_phylum)


### Class
res_class=levels(t_matrix_otu_species$Class)
for(col in 2:29){
  res=aggregate(t_matrix_otu_species[,col],t_matrix_otu_species[,"Class"],sum)
  res_class=cbind(res_class,res[,2])
}

row.names(res_class)<-res_class[,1]
res_class <- apply(as.data.frame(res_class[,2:29]),1,FUN=as.integer)

shannon_tot_class <- diversity(res_class)


### Order
res_order=levels(t_matrix_otu_species$Order)
for(col in 2:29){
  res=aggregate(t_matrix_otu_species[,col],t_matrix_otu_species[,"Order"],sum)
  res_order=cbind(res_order,res[,2])
}

row.names(res_order)<-res_order[,1]
res_order <- apply(as.data.frame(res_order[,2:29]),1,FUN=as.integer)

shannon_tot_order <- diversity(res_order)

### Family
res_family=levels(t_matrix_otu_species$Family)
for(col in 2:29){
  res=aggregate(t_matrix_otu_species[,col],t_matrix_otu_species[,"Family"],sum)
  res_family=cbind(res_family,res[,2])
}

row.names(res_family)<-res_family[,1]
res_family <- apply(as.data.frame(res_family[,2:29]),1,FUN=as.integer)

shannon_tot_family <- diversity(res_family)


### Genus
res_genus=levels(t_matrix_otu_species$Genus)
for(col in 2:29){
  res=aggregate(t_matrix_otu_species[,col],t_matrix_otu_species[,"Genus"],sum)
  res_genus=cbind(res_genus,res[,2])
}
row.names(res_genus)<-res_genus[,1]
res_genus <- apply(as.data.frame(res_genus[,2:29]),1,FUN=as.integer)

shannon_tot_genus <- diversity(res_genus)


### On représente tous les Shannon
plot(x=1:28,y=shannon_tot,main="Shannon selon regroupements pour les 28 individus",ylim=c(0,4.5))
points(x=1:28,y=shannon_tot_phylum, col='red')
points(x=1:28,y=shannon_tot_class, col='orange')
points(x=1:28,y=shannon_tot_order, col='blue')
points(x=1:28,y=shannon_tot_family, col='green3')
points(x=1:28,y=shannon_tot_genus, col='purple')
legend('topleft', legend=c("Ensemble","Genus","Family","Order","Class","Phylum"), fill = c('black', "purple", "green3","blue","orange","red"), title='Regroupements')


#Corrélation----------------------
   #Phylum et Order: R²=90%
plot(shannon_tot_phylum,shannon_tot_order)
summary(lm(shannon_tot_phylum~shannon_tot_order))

   #Phylum et Family: R²=29%
plot(shannon_tot_phylum,shannon_tot_family)
summary(lm(shannon_tot_phylum~shannon_tot_family))

   #Phylum et ensemble: R²=6.7%
plot(shannon_tot_phylum,shannon_tot)
summary(lm(shannon_tot_phylum~shannon_tot))


  #Genus et ensemble: R²=67%
plot(shannon_tot_genus,shannon_tot)
summary(lm(shannon_tot_genus~shannon_tot))


# To do...

Run models with "Strain" - may have some impacts

# DR-longevity paper

#basic info

#> summary(factor(Model))
#  0   1 
#190 339 
#> summary(factor(Sex))
#  F   M  MF   N 
#253 197  10  69 

###### Data sort

ZPro<-scale(Prots)
ZCR<-scale(c_nCR)

ZPCR<-data.frame(ZPro=ZPro,ZCR=ZCR)

write.csv(ZPCR,"ZPCR.csv")

###############


# Mac
setwd("/Users/naksh50p/Dropbox/Project1/DR-longevity/Analysis")

rm(list=ls())
detach(Data)

library(MCMCglmm)

Data<-read.csv("DataAll.csv")
attach(Data)
str(Data)

Tree<-read.tree("Tree2.tre")
str(Tree)

#################

# Linax
setwd("/home/labadmin/Dropbox/Project1/DR-longevity/Analysis")
# Linax
library(MCMCglmm)
rm(list=ls())
detach(Data)

Data<-read.csv("Data.csv")
attach(Data)
str(Data)

Tree<-read.tree("Tree.tre")
str(Tree)



#########
##Model 0 - not species effect
#########
prior0<-list(R=list(V=1E-10,nu=-1),G=list(G1=list(V=1E-10,nu=-1)))


m0.1<-MCMCglmm(c_all_logHR~1,random=~Study_ID,family="gaussian",data=Data,verbose=F,mev=c_all_varlogHR,nitt=50000,thin=25,burnin=25000,pr=T,prior=prior0)

m0.2<-MCMCglmm(c_all_logHR~1,random=~Study_ID,family="gaussian",data=Data,verbose=F,mev=c_all_varlogHR,nitt=50000,thin=25,burnin=25000,pr=T,prior=prior0)

m0.3<-MCMCglmm(c_all_logHR~1,random=~Study_ID,family="gaussian",data=Data,verbose=F,mev=c_all_varlogHR,nitt=50000,thin=25,burnin=25000,pr=T,prior=prior0)

gelman.diag(list(m0.1$Sol[,1],m0.2$Sol[,1],m0.3$Sol[,1]))
gelman.diag(list(m0.1$VCV[,c(1,3)],m0.2$VCV[,c(1,3)],m0.3$VCV[,c(1,3)]))
gelman.diag(list(m0.1$Deviance,m0.2$Deviance,m0.3$Deviance))

save(m0.1,m0.2,m0.3,file="model0.RData")

# DIC around 1116 - worse a lot worse - may be mention???


##########
## Model 1 without phylogeny (Model 1)
##########
prior1<-list(R=list(V=1E-10,nu=-1),G=list(G1=list(V=1E-10,nu=-1),G2=list(V=1E-10,nu=-1)))

m1.1<-MCMCglmm(c_all_logHR~1,random=~Species_ID+Study_ID,family="gaussian",data=Data,verbose=F,mev=c_all_varlogHR,nitt=50000,thin=25,burnin=25000,pr=T,prior=prior1)

m1.2<-MCMCglmm(c_all_logHR~1,random=~Species_ID+Study_ID,family="gaussian",data=Data,verbose=F,mev=c_all_varlogHR,nitt=50000,thin=25,burnin=25000,pr=T,prior=prior1)

m1.3<-MCMCglmm(c_all_logHR~1,random=~Species_ID+Study_ID,family="gaussian",data=Data,verbose=F,mev=c_all_varlogHR,nitt=50000,thin=25,burnin=25000,pr=T,prior=prior1)

gelman.diag(list(m1.1$Sol[,1],m1.2$Sol[,1],m1.3$Sol[,1]))
gelman.diag(list(m1.1$VCV[,c(1,2,4)],m1.2$VCV[,c(1,2,4)],m1.3$VCV[,c(1,2,4)]))
gelman.diag(list(m1.1$Deviance,m1.2$Deviance,m1.3$Deviance))


save(m1.1,m1.2,m1.3,file="model1.RData")

# DIC around 1108

##########
## Model 2 with phylogeny (Model 2)
##########
prior2<-list(R=list(V=1E-10,nu=-1),G=list(G1=list(V=1E-10,nu=-1),G2=list(V=1E-10,nu=-1),G3=list(V=1E-10,nu=-1)))

m2.1<-MCMCglmm(c_all_logHR~1,random=~animal+Species_ID+Study_ID,family="gaussian",data=Data,pedigree=Tree,verbose=F,mev=c_all_varlogHR,scale=F,nitt=50000,thin=25,burnin=25000,prior=prior2,pr=T)

m2.2<-MCMCglmm(c_all_logHR~1,random=~animal+Species_ID+Study_ID,family="gaussian",data=Data,pedigree=Tree,verbose=F,mev=c_all_varlogHR,scale=F,nitt=50000,thin=25,burnin=25000,prior=prior2,pr=T)

m2.3<-MCMCglmm(c_all_logHR~1,random=~animal+Species_ID+Study_ID,family="gaussian",data=Data,pedigree=Tree,verbose=F,mev=c_all_varlogHR,scale=F,nitt=50000,thin=25,burnin=25000,prior=prior2,pr=T)

gelman.diag(list(m2.1$Sol[,1],m2.2$Sol[,1],m2.3$Sol[,1]))
gelman.diag(list(m2.1$VCV[,c(1:3,5)],m2.2$VCV[,c(1:3,5)],m2.3$VCV[,c(1:3,5)]))
gelman.diag(list(m2.1$Deviance,m2.2$Deviance,m2.3$Deviance))

save(m2.1,m2.2,m2.3,file="model2.RData")


# note scaling does not chage the results a bit!
# DIC around 1107

#######################
##Model 3 with phylogeny witout species effect
###############################
# ignore this.....

prior3<-list(R=list(V=1E-10,nu=-1),G=list(G1=list(V=1E-10,nu=-1),G2=list(V=1E-10,nu=-1)))

m3.1<-MCMCglmm(c_all_logHR~1,random=~animal+Study_ID,family="gaussian",data=Data,verbose=F,mev=c_all_varlogHR,nitt=50000,thin=25,burnin=25000,pr=T,prior=prior3,scale=F,pedigree=Tree)

m3.2<-MCMCglmm(c_all_logHR~1,random=~animal+Study_ID,family="gaussian",data=Data,verbose=F,mev=c_all_varlogHR,nitt=50000,thin=25,burnin=25000,pr=T,prior=prior3,scale=F,pedigree=Tree)

m3.3<-MCMCglmm(c_all_logHR~1,random=~animal+Study_ID,family="gaussian",data=Data,verbose=F,mev=c_all_varlogHR,nitt=50000,thin=25,burnin=25000,pr=T,prior=prior3,scale=F,pedigree=Tree)

gelman.diag(list(m3.1$Sol[,1],m3.2$Sol[,1],m3.3$Sol[,1]))
gelman.diag(list(m3.1$VCV[,c(1,2,4)],m3.2$VCV[,c(1,2,4)],m3.3$VCV[,c(1,2,4)]))
gelman.diag(list(m3.1$Deviance,m3.2$Deviance,m3.3$Deviance))

# DIC around 1107.5 - not considered

save(m3.1,m3.2,m3.3,file="model3.RData")


################
# Mac again.
setwd("/Users/naksh50p/Dropbox/Project1/DR-longevity/Analysis")

library(MCMCglmm)
load("model1.RData")

gelman.diag(list(m1.1$Sol[,1],m1.2$Sol[,1],m1.3$Sol[,1]))
gelman.diag(list(m1.1$VCV[,c(1,2,4)],m1.2$VCV[,c(1,2,4)],m1.3$VCV[,c(1,2,4)]))
gelman.diag(list(m1.1$Deviance,m1.2$Deviance,m1.3$Deviance))

m1.1$DIC
m1.2$DIC
m1.3$DIC

load("model2.RData")

gelman.diag(list(m2.1$Sol[,1],m2.2$Sol[,1],m2.3$Sol[,1]))
gelman.diag(list(m2.1$VCV[,c(1:3,5)],m2.2$VCV[,c(1:3,5)],m2.3$VCV[,c(1:3,5)]))
gelman.diag(list(m2.1$Deviance,m2.2$Deviance,m2.3$Deviance))

m2.1$DIC
m2.2$DIC
m2.3$DIC

####################################
# paramters estiamtes for Model 1 & 2
#####################################

# meta-analytic mean Model 1
posterior.mode(m1.1$Sol[,1])
summary(m1.1$Sol[,1])
HPDinterval(m1.1$Sol[,1])

# meta-analytic mean Model 2
posterior.mode(m2.1$Sol[,1])
summary(m2.1$Sol[,1])
HPDinterval(m2.1$Sol[,1])



#VCV Model 1

posterior.mode(m1.1$VCV)
summary(m1.1$VCV)
HPDinterval(m1.1$VCV)

# getting "Typical" measurment error variance

W<-1/c_all_varlogHR

s2<-sum(W*(length(W)-1))/(sum(W)^2-sum(W^2)) 

I2<-100*(m1.1$VCV[,"Species_ID"]+m1.1$VCV[,"Study_ID"])/(m1.1$VCV[,"Species_ID"]+m1.1$VCV[,"Study_ID"]+m1.1$VCV[,"units"]+s2)

posterior.mode(I2)
summary(I2)
HPDinterval(I2)


#VCV Model 2
posterior.mode(m2.1$VCV)
summary(m2.1$VCV)
HPDinterval(m2.1$VCV)

H2<-m2.1$VCV[,"animal"]/(m2.1$VCV[,"animal"]+m2.1$VCV[,"Species_ID"]+m2.1$VCV[,"Study_ID"]+m2.1$VCV[,"units"]+s2)

posterior.mode(H2)
summary(H2)
HPDinterval(H2)

I2<-100*(m2.1$VCV[,"animal"]+m2.1$VCV[,"Species_ID"]+m2.1$VCV[,"Study_ID"])/(m2.1$VCV[,"animal"]+m2.1$VCV[,"Species_ID"]+m2.1$VCV[,"Study_ID"]+m2.1$VCV[,"units"]+s2)

posterior.mode(I2)
summary(I2)
HPDinterval(I2)


######################

# Sex and Model effects....

###########
##Model 4 with phylogeny
#############
prior4<-list(R=list(V=1E-10,nu=-1),G=list(G1=list(V=1E-10,nu=-1),G2=list(V=1E-10,nu=-1),G3=list(V=1E-10,nu=-1)))

m4.1<-MCMCglmm(c_all_logHR~Sex+Model,random=~animal+Species_ID+Study_ID,family="gaussian",data=Data,pedigree=Tree,verbose=F,mev=c_all_varlogHR,scale=F,nitt=50000,thin=25,burnin=25000,prior=prior4,pr=T)

m4.2<-MCMCglmm(c_all_logHR~Sex+Model,random=~animal+Species_ID+Study_ID,family="gaussian",data=Data,pedigree=Tree,verbose=F,mev=c_all_varlogHR,scale=F,nitt=50000,thin=25,burnin=25000,prior=prior4,pr=T)

m4.3<-MCMCglmm(c_all_logHR~Sex+Model,random=~animal+Species_ID+Study_ID,family="gaussian",data=Data,pedigree=Tree,verbose=F,mev=c_all_varlogHR,scale=F,nitt=50000,thin=25,burnin=25000,prior=prior4,pr=T)

gelman.diag(list(m4.1$Sol[,1:5],m4.2$Sol[,1:5],m4.3$Sol[,1:5]))
gelman.diag(list(m4.1$VCV[,c(1:3,5)],m4.2$VCV[,c(1:3,5)],m4.3$VCV[,c(1:3,5)]))
gelman.diag(list(m4.1$Deviance,m4.2$Deviance,m4.3$Deviance))

save(m4.1,m4.2,m4.3,file="model4.RData")



###########
##Model 5 without phylogeny
#############
prior5<-list(R=list(V=1E-10,nu=-1),G=list(G1=list(V=1E-10,nu=-1),G2=list(V=1E-10,nu=-1)))

m5.1<-MCMCglmm(c_all_logHR~Sex+Model,random=~Species_ID+Study_ID,family="gaussian",data=Data,verbose=F,mev=c_all_varlogHR,nitt=50000,thin=25,burnin=25000,prior=prior5,pr=T)

m5.2<-MCMCglmm(c_all_logHR~Sex+Model,random=~Species_ID+Study_ID,family="gaussian",data=Data,verbose=F,mev=c_all_varlogHR,nitt=50000,thin=25,burnin=25000,prior=prior5,pr=T)

m5.3<-MCMCglmm(c_all_logHR~Sex+Model,random=~Species_ID+Study_ID,family="gaussian",data=Data,verbose=F,mev=c_all_varlogHR,nitt=50000,thin=25,burnin=25000,prior=prior5,pr=T)

gelman.diag(list(m5.1$Sol[,1:5],m5.2$Sol[,1:5],m5.3$Sol[,1:5]))
gelman.diag(list(m5.1$VCV[,c(1:2,4)],m5.2$VCV[,c(1:2,4)],m5.3$VCV[,c(1:2,4)]))
gelman.diag(list(m5.1$Deviance,m5.2$Deviance,m5.3$Deviance))

save(m5.1,m5.2,m5.3,file="model5.RData")


# just with Study effect is no better around DIC = 1108.5

############
####################################
# paramters estiamtes for Model 4 & 5 
#####################################

library(MCMCglmm)
load("model4.RData")

gelman.diag(list(m4.1$Sol[,1:3],m4.2$Sol[,1:3],m4.3$Sol[,1:3]))
gelman.diag(list(m4.1$VCV[,c(1:3,5)],m4.2$VCV[,c(1:3,5)],m4.3$VCV[,c(1:3,5)]))
gelman.diag(list(m4.1$Deviance,m4.2$Deviance,m4.3$Deviance))

m4.1$DIC
m4.2$DIC
m4.3$DIC

load("model5.RData")

gelman.diag(list(m5.1$Sol[,1:5],m5.2$Sol[,1:5],m5.3$Sol[,1:5]))
gelman.diag(list(m5.1$VCV[,c(1:2,4)],m5.2$VCV[,c(1:2,4)],m5.3$VCV[,c(1:2,4)]))
gelman.diag(list(m5.1$Deviance,m5.2$Deviance,m5.3$Deviance))

m5.1$DIC
m5.2$DIC
m5.3$DIC

# Sol
# Model 4 with phylo
posterior.mode(m4.1$Sol[,1:5])
summary(m4.1$Sol[,1:5])
HPDinterval(m4.1$Sol[,1:100])

# model spp. 
summary(m4.1$Sol[,1:4]+m4.1$Sol[,5])

# Model 5 without phylo
posterior.mode(m5.1$Sol[,1:5])
summary(m5.1$Sol[,1:5])
HPDinterval(m5.1$Sol[,1:5])

#VCV
# Model 4## (Model 4)
posterior.mode(m4.1$VCV)
summary(m4.1$VCV)
HPDinterval(m4.1$VCV)

H2<-m4.1$VCV[,"animal"]/(m4.1$VCV[,"animal"]+m4.1$VCV[,"Species_ID"]+m4.1$VCV[,"Study_ID"]+m4.1$VCV[,"units"]+s2)

posterior.mode(H2)
summary(H2)
HPDinterval(H2)

W<-1/c_all_varlogHR

s2<-sum(W*(length(W)-1))/(sum(W)^2-sum(W^2))

I2<-100*(m4.1$VCV[,"animal"]+m4.1$VCV[,"Species_ID"]+m4.1$VCV[,"Study_ID"])/(m4.1$VCV[,"animal"]+m4.1$VCV[,"Species_ID"]+m4.1$VCV[,"Study_ID"]+m4.1$VCV[,"units"]+s2)

posterior.mode(I2)
summary(I2)
HPDinterval(I2)


# Model 5### (Model 3)
posterior.mode(m5.1$VCV)
summary(m5.1$VCV)
HPDinterval(m5.1$VCV)

I2<-100*(m5.1$VCV[,"Species_ID"]+m5.1$VCV[,"Study_ID"])/(m5.1$VCV[,"Species_ID"]+m5.1$VCV[,"Study_ID"]+m5.1$VCV[,"units"]+s2)

posterior.mode(I2)
summary(I2)
HPDinterval(I2)

########################
##Model 6 with phylogeny
########################
prior6<-list(R=list(V=1E-10,nu=-1),G=list(G1=list(V=1E-10,nu=-1),G2=list(V=1E-10,nu=-1),G3=list(V=1E-10,nu=-1)))

m6.1<-MCMCglmm(c_all_logHR~Sex+Model+c_nCR+I(c_nCR^2)+Prots+I(Prots^2),random=~animal+Species_ID+Study_ID,family="gaussian",data=Data,pedigree=Tree,verbose=F,mev=c_all_varlogHR,scale=F,nitt=50000,thin=25,burnin=25000,prior=prior6,pr=T)

m6.2<-MCMCglmm(c_all_logHR~Sex+Model+c_nCR+I(c_nCR^2)+Prots+I(Prots^2),random=~animal+Species_ID+Study_ID,family="gaussian",data=Data,pedigree=Tree,verbose=F,mev=c_all_varlogHR,scale=F,nitt=50000,thin=25,burnin=25000,prior=prior6,pr=T)

m6.3<-MCMCglmm(c_all_logHR~Sex+Model+c_nCR+I(c_nCR^2)+Prots+I(Prots^2),random=~animal+Species_ID+Study_ID,family="gaussian",data=Data,pedigree=Tree,verbose=F,mev=c_all_varlogHR,scale=F,nitt=50000,thin=25,burnin=25000,prior=prior6,pr=T)

gelman.diag(list(m6.1$Sol[,1:9],m6.2$Sol[,1:9],m6.3$Sol[,1:9]))
gelman.diag(list(m6.1$VCV[,c(1:3,5)],m6.2$VCV[,c(1:3,5)],m6.3$VCV[,c(1:3,5)]))
gelman.diag(list(m6.1$Deviance,m6.2$Deviance,m6.3$Deviance))

m6.1$DIC
m6.2$DIC
m6.3$DIC

save(m6.1,m6.2,m6.3,file="model6.RData")



###########################
##Model 7 without phylogeny
###########################
prior7<-list(R=list(V=1E-10,nu=-1),G=list(G1=list(V=1E-10,nu=-1),G2=list(V=1E-10,nu=-1)))

m7.1<-MCMCglmm(c_all_logHR~Sex+Model+c_nCR+I(c_nCR^2)+Prots+I(Prots^2),random=~Species_ID+Study_ID,family="gaussian",data=Data,verbose=F,mev=c_all_varlogHR,nitt=50000,thin=25,burnin=25000,prior=prior5,pr=T)

m7.2<-MCMCglmm(c_all_logHR~Sex+Model+c_nCR+I(c_nCR^2)+Prots+I(Prots^2),random=~Species_ID+Study_ID,family="gaussian",data=Data,verbose=F,mev=c_all_varlogHR,nitt=50000,thin=25,burnin=25000,prior=prior5,pr=T)

m7.3<-MCMCglmm(c_all_logHR~Sex+Model+c_nCR+I(c_nCR^2)+Prots+I(Prots^2),random=~Species_ID+Study_ID,family="gaussian",data=Data,verbose=F,mev=c_all_varlogHR,nitt=50000,thin=25,burnin=25000,prior=prior5,pr=T)


gelman.diag(list(m7.1$Sol[,1:9],m7.2$Sol[,1:9],m7.3$Sol[,1:9]))
gelman.diag(list(m7.1$VCV[,c(1:2,4)],m7.2$VCV[,c(1:2,4)],m7.3$VCV[,c(1:2,4)]))
gelman.diag(list(m7.1$Deviance,m7.2$Deviance,m7.3$Deviance))

m7.1$DIC
m7.2$DIC
m7.3$DIC

save(m7.1,m7.2,m7.3,file="model7.RData")

#############################

library(MCMCglmm)

load("model6.RData")

gelman.diag(list(m6.1$Sol[,1:9],m6.2$Sol[,1:9],m6.3$Sol[,1:9]))
gelman.diag(list(m6.1$VCV[,c(1:3,5)],m6.2$VCV[,c(1:3,5)],m6.3$VCV[,c(1:3,5)]))
gelman.diag(list(m6.1$Deviance,m6.2$Deviance,m6.3$Deviance))

m6.1$DIC
m6.2$DIC
m6.3$DIC

load("model7.RData")

gelman.diag(list(m7.1$Sol[,1:9],m7.2$Sol[,1:9],m7.3$Sol[,1:9]))
gelman.diag(list(m7.1$VCV[,c(1:2,4)],m7.2$VCV[,c(1:2,4)],m7.3$VCV[,c(1:2,4)]))
gelman.diag(list(m7.1$Deviance,m7.2$Deviance,m7.3$Deviance))

m7.1$DIC
m7.2$DIC
m7.3$DIC


# Sol
# Model 6
posterior.mode(m6.1$Sol[,6:9])
summary(m6.1$Sol[,6:9])
HPDinterval(m6.1$Sol[,6:9])

# weighted intercept

WI<-((253/529)*m6.1$Sol[,1]+(193/529)*m6.1$Sol[,2]+(10/529)*m6.1$Sol[,3]+(69/529)*m6.1$Sol[,4])*0.5+((190/529)*m6.1$Sol[,1]+(339/529)*m6.1$Sol[,5])*0.5

posterior.mode(WI)
summary(WI)
HPDinterval(WI)


# Model 7 (Model 5)
posterior.mode(m7.1$Sol[,6:9])
summary(m7.1$Sol[,6:9])
HPDinterval(m7.1$Sol[,6:9])

WI<-((253/529)*m7.1$Sol[,1]+(193/529)*m7.1$Sol[,2]+(10/529)*m7.1$Sol[,3]+(69/529)*m7.1$Sol[,4])*0.5+((190/529)*m7.1$Sol[,1]+(339/529)*m7.1$Sol[,5])*0.5

posterior.mode(WI)
summary(WI)
HPDinterval(WI)


#VCV
# Model 6##
posterior.mode(m6.1$VCV)
summary(m6.1$VCV)
HPDinterval(m6.1$VCV)

H2<-m6.1$VCV[,"animal"]/(m6.1$VCV[,"animal"]+m6.1$VCV[,"Species_ID"]+m6.1$VCV[,"Study_ID"]+m6.1$VCV[,"units"]+s2)

posterior.mode(H2)
summary(H2)
HPDinterval(H2)

# getting s2

W<-1/c_all_varlogHR

s2<-sum(W*(length(W)-1))/(sum(W)^2-sum(W^2))

I2<-100*(m6.1$VCV[,"animal"]+m6.1$VCV[,"Species_ID"]+m6.1$VCV[,"Study_ID"])/(m6.1$VCV[,"animal"]+m6.1$VCV[,"Species_ID"]+m6.1$VCV[,"Study_ID"]+m6.1$VCV[,"units"]+s2)

posterior.mode(I2)
summary(I2)
HPDinterval(I2)


# Model 7### (model 5)
posterior.mode(m7.1$VCV)
summary(m7.1$VCV)
HPDinterval(m7.1$VCV)

I2<-100*(m7.1$VCV[,"Species_ID"]+m7.1$VCV[,"Study_ID"])/(m7.1$VCV[,"Species_ID"]+m7.1$VCV[,"Study_ID"]+m7.1$VCV[,"units"]+s2)

posterior.mode(I2)
summary(I2)
HPDinterval(I2)

# publication bias

load("model6.RData")

m6.1$Random$formula<-update(m6.1$Random$formula, ~.+leg(mev, -1, FALSE):units)

#Pred1<-predict(m10.1)

Pred<-predict(m6.1, marginal=~leg(mev, -1, FALSE):units)
Prec<-sqrt(1/c_all_varlogHR)
ES<-c_all_logHR-predict(m6.1, marginal=~leg(mev, -1, FALSE):units)

plot(zES,Prec)
plot(Prec,ES)
plot(Prec,ES)


Data2<-data.frame(ES=ES,Prec=Prec)

zES<-ES*Prec

#model<-MCMCglmm(zES~Prec,data=Data2,verbose=F)
#summary(model$Sol)
#HPDinterval(model$Sol)




e0.1<-MCMCglmm(zES~I(1/sqrt(c_all_varlogHR)),family="gaussian",data=Data,verbose=F,nitt=50000,thin=25,burnin=25000)

e0.2<-MCMCglmm(zES~I(1/sqrt(c_all_varlogHR)),family="gaussian",data=Data,verbose=F,nitt=50000,thin=25,burnin=25000)


e0.3<-MCMCglmm(zES~I(1/sqrt(c_all_varlogHR)),family="gaussian",data=Data,verbose=F,nitt=50000,thin=25,burnin=25000)

save(e0.1,e0.2,e0.3,file="egger0.RData")

######
# Orignal data

zES1<-c_all_logHR/sqrt(c_all_varlogHR)


e0.1<-MCMCglmm(zES1~I(1/sqrt(c_all_varlogHR)),family="gaussian",data=Data,verbose=F,nitt=50000,thin=25,burnin=25000)

e0.2<-MCMCglmm(zES1~I(1/sqrt(c_all_varlogHR)),family="gaussian",data=Data,verbose=F,nitt=50000,thin=25,burnin=25000)


e0.3<-MCMCglmm(zES1~I(1/sqrt(c_all_varlogHR)),family="gaussian",data=Data,verbose=F,nitt=50000,thin=25,burnin=25000)

save(e0.1,e0.2,e0.3,file="egger.RData")



#####################
#load("egger0.RData")

load("egger.RData")
library(MCMCglmm)

gelman.diag(list(e0.1$Sol,e0.2$Sol,e0.3$Sol))
gelman.diag(list(e0.1$VCV,e0.2$VCV,e0.3$VCV))
gelman.diag(list(e0.1$Deviance,e0.2$Deviance,e0.3$Deviance))

e0.1$DIC
e0.2$DIC
e0.3$DIC



# fixed effects
posterior.mode(e0.1$Sol[,1:2])
summary(e0.1$Sol[,1:2])
HPDinterval(e0.1$Sol[,1:2])
plot(e0.1$Sol)

# random effects

posterior.mode(e0.1$VCV)
summary(e0.1$VCV)
HPDinterval(e0.1$VCV)
plot(e0.1$VCV)

plot(zES,Prec)


######################################################################################
####################################################################################
##############
####Extra
#############
###############

###############
# preparation...

# Mac
setwd("/Users/naksh50p/Dropbox/Project1/DR-longevity/Analysis")

rm(list=ls())
detach(Data)

library(MCMCglmm)

Data<-read.csv("DataAll.csv")
attach(Data)
str(Data)

Tree<-read.tree("Tree2.tre")
str(Tree)
#########################
#########################



# factors to be considered
# Repro, FoodSched, Type, AL, Intake???, LT25contr (Organism Age):

# probably I should centre the stuff. 

#########
##Model 8 (Model 7)
#########

prior8<-list(R=list(V=1E-10,nu=-1),G=list(G1=list(V=1E-10,nu=-1),G2=list(V=1E-10,nu=-1),G3=list(V=1E-10,nu=-1)))


m8.1<-MCMCglmm(c_all_logHR~Sex+Model+c_nCR+I(c_nCR^2)+Prots+I(Prots^2)+Repro,random=~animal+Species_ID+Study_ID,family="gaussian",data=Data,pedigree=Tree,verbose=F,mev=c_all_varlogHR,scale=F,nitt=100000,thin=50,burnin=50000,prior=prior8,pr=T)

m8.2<-MCMCglmm(c_all_logHR~Sex+Model+c_nCR+I(c_nCR^2)+Prots+I(Prots^2)+Repro,random=~animal+Species_ID+Study_ID,family="gaussian",data=Data,pedigree=Tree,verbose=F,mev=c_all_varlogHR,scale=F,nitt=100000,thin=50,burnin=50000,prior=prior8,pr=T)

m8.3<-MCMCglmm(c_all_logHR~Sex+Model+c_nCR+I(c_nCR^2)+Prots+I(Prots^2)+Repro,random=~animal+Species_ID+Study_ID,family="gaussian",data=Data,pedigree=Tree,verbose=F,mev=c_all_varlogHR,scale=F,nitt=100000,thin=50,burnin=50000,prior=prior8,pr=T)

save(m8.1,m8.2,m8.3,file="model8.RData")

#############################

library(MCMCglmm)

load("model8.RData")

gelman.diag(list(m8.1$Sol[,1:12],m8.2$Sol[,1:12],m8.3$Sol[,1:12]))
gelman.diag(list(m8.1$VCV[,c(1:3,5)],m8.2$VCV[,c(1:3,5)],m8.3$VCV[,c(1:3,5)]))
gelman.diag(list(m8.1$Deviance,m8.2$Deviance,m8.3$Deviance))

m8.1$DIC
m8.2$DIC
m8.3$DIC


# Fixed effects
posterior.mode(m8.1$Sol[,1:12])
summary(m8.1$Sol[,1:12])
HPDinterval(m8.1$Sol[,1:12])
plot(m8.1$Sol[,1:12])

#WI<-((253/529)*m8.1$Sol[,1]+(193/529)*m8.1$Sol[,2]+(10/529)*m8.1$Sol[,3]+(69/529)*m8.1$Sol[,4])*0.5+((190/529)*m8.1$Sol[,1]+(339/529)*m8.1$Sol[,5])*0.5

#posterior.mode(WI)
#summary(WI)
#HPDinterval(WI)

# Random effects
posterior.mode(m8.1$VCV)
summary(m8.1$VCV)
HPDinterval(m8.1$VCV)
plot(m8.1$VCV)

H2<-m8.1$VCV[,"animal"]/(m8.1$VCV[,"animal"]+m8.1$VCV[,"Species_ID"]+m8.1$VCV[,"Study_ID"]+m8.1$VCV[,"units"]+s2)

posterior.mode(H2)
summary(H2)
HPDinterval(H2)

# getting s2

W<-1/c_all_varlogHR

s2<-sum(W*(length(W)-1))/(sum(W)^2-sum(W^2))

I2<-100*(m8.1$VCV[,"animal"]+m8.1$VCV[,"Species_ID"]+m8.1$VCV[,"Study_ID"])/(m8.1$VCV[,"animal"]+m8.1$VCV[,"Species_ID"]+m8.1$VCV[,"Study_ID"]+m8.1$VCV[,"units"]+s2)

posterior.mode(I2)
summary(I2)
HPDinterval(I2)



########################
##Model 9 (Model 8)
###########################


prior9<-list(R=list(V=1E-10,nu=-1),G=list(G1=list(V=1E-10,nu=-1),G2=list(V=1E-10,nu=-1),G3=list(V=1E-10,nu=-1)))

m9.1<-MCMCglmm(c_all_logHR~Sex+Model+c_nCR+I(c_nCR^2)+Prots+I(Prots^2)+FoodSched,random=~animal+Species_ID+Study_ID,family="gaussian",data=Data,pedigree=Tree,verbose=F,mev=c_all_varlogHR,scale=F,nitt=100000,thin=50,burnin=50000,prior=prior9,pr=T)

m9.2<-MCMCglmm(c_all_logHR~Sex+Model+c_nCR+I(c_nCR^2)+Prots+I(Prots^2)+FoodSched,random=~animal+Species_ID+Study_ID,family="gaussian",data=Data,pedigree=Tree,verbose=F,mev=c_all_varlogHR,scale=F,nitt=100000,thin=50,burnin=50000,prior=prior9,pr=T)

m9.3<-MCMCglmm(c_all_logHR~Sex+Model+c_nCR+I(c_nCR^2)+Prots+I(Prots^2)+FoodSched,random=~animal+Species_ID+Study_ID,family="gaussian",data=Data,pedigree=Tree,verbose=F,mev=c_all_varlogHR,scale=F,nitt=100000,thin=50,burnin=50000,prior=prior9,pr=T)

save(m9.1,m9.2,m9.3,file="model9.RData")

#############################

library(MCMCglmm)

load("model9.RData")

gelman.diag(list(m9.1$Sol[,1:12],m9.2$Sol[,1:12],m9.3$Sol[,1:12]))
gelman.diag(list(m9.1$VCV[,c(1:3,5)],m9.2$VCV[,c(1:3,5)],m9.3$VCV[,c(1:3,5)]))
gelman.diag(list(m9.1$Deviance,m9.2$Deviance,m9.3$Deviance))

m9.1$DIC
m9.2$DIC
m9.3$DIC


# Fixed effects
posterior.mode(m9.1$Sol[,1:12])
summary(m9.1$Sol[,1:12])
HPDinterval(m9.1$Sol[,1:12])
plot(m9.1$Sol[,1:12])

# Random effects
posterior.mode(m9.1$VCV)
summary(m9.1$VCV)
HPDinterval(m9.1$VCV)
plot(m9.1$VCV)

# diff 

posterior.mode(m9.1$Sol[,1]+m9.1$Sol[,10]-(m9.1$Sol[,1]+m9.1$Sol[,11]))
summary(m9.1$Sol[,1]+m9.1$Sol[,10]-(m9.1$Sol[,1]+m9.1$Sol[,11]))
HPDinterval(m9.1$Sol[,1]+m9.1$Sol[,10]-(m9.1$Sol[,1]+m9.1$Sol[,11]))



H2<-m9.1$VCV[,"animal"]/(m9.1$VCV[,"animal"]+m9.1$VCV[,"Species_ID"]+m9.1$VCV[,"Study_ID"]+m9.1$VCV[,"units"]+s2)

posterior.mode(H2)
summary(H2)
HPDinterval(H2)

# getting s2

W<-1/c_all_varlogHR

s2<-sum(W*(length(W)-1))/(sum(W)^2-sum(W^2))

I2<-100*(m9.1$VCV[,"animal"]+m9.1$VCV[,"Species_ID"]+m9.1$VCV[,"Study_ID"])/(m9.1$VCV[,"animal"]+m9.1$VCV[,"Species_ID"]+m9.1$VCV[,"Study_ID"]+m9.1$VCV[,"units"]+s2)

posterior.mode(I2)
summary(I2)
HPDinterval(I2)

##########
##Model 10 (Model 9)
##########
##########

prior10<-list(R=list(V=1E-10,nu=-1),G=list(G1=list(V=1E-10,nu=-1),G2=list(V=1E-10,nu=-1),G3=list(V=1E-10,nu=-1)))

m10.1<-MCMCglmm(c_all_logHR~Sex+Model+c_nCR+I(c_nCR^2)+Prots+I(Prots^2)+Type,random=~animal+Species_ID+Study_ID,family="gaussian",data=Data,pedigree=Tree,verbose=F,mev=c_all_varlogHR,scale=F,nitt=100000,thin=50,burnin=50000,prior=prior10,pr=T)

m10.2<-MCMCglmm(c_all_logHR~Sex+Model+c_nCR+I(c_nCR^2)+Prots+I(Prots^2)+Type,random=~animal+Species_ID+Study_ID,family="gaussian",data=Data,pedigree=Tree,verbose=F,mev=c_all_varlogHR,scale=F,nitt=100000,thin=50,burnin=50000,prior=prior10,pr=T)

m10.3<-MCMCglmm(c_all_logHR~Sex+Model+c_nCR+I(c_nCR^2)+Prots+I(Prots^2)+Type,random=~animal+Species_ID+Study_ID,family="gaussian",data=Data,pedigree=Tree,verbose=F,mev=c_all_varlogHR,scale=F,nitt=100000,thin=50,burnin=50000,prior=prior10,pr=T)


save(m10.1,m10.2,m10.3,file="model10.RData")

#############################

library(MCMCglmm)

load("model10.RData")

gelman.diag(list(m10.1$Sol[,1:12],m10.2$Sol[,1:12],m10.3$Sol[,1:12]))
gelman.diag(list(m10.1$VCV[,c(1:3,5)],m10.2$VCV[,c(1:3,5)],m10.3$VCV[,c(1:3,5)]))
gelman.diag(list(m10.1$Deviance,m10.2$Deviance,m10.3$Deviance))

m10.1$DIC
m10.2$DIC
m10.3$DIC



# Fixed effects
posterior.mode(m10.1$Sol[,1:14])
summary(m10.1$Sol[,1:14])
HPDinterval(m10.1$Sol[,1:14])
plot(m10.1$Sol[,1:14])

# Random effects
posterior.mode(m10.1$VCV)
summary(m10.1$VCV)
HPDinterval(m10.1$VCV)
plot(m10.1$VCV)

# diff CNM & FC (11)

posterior.mode(m10.1$Sol[,1]+m10.1$Sol[,10]-(m10.1$Sol[,1]+m10.1$Sol[,11]))
summary(m10.1$Sol[,1]+m10.1$Sol[,10]-(m10.1$Sol[,1]+m10.1$Sol[,11]))
HPDinterval(m10.1$Sol[,1]+m10.1$Sol[,10]-(m10.1$Sol[,1]+m10.1$Sol[,11]))

# diff CNM & FD (12)

posterior.mode(m10.1$Sol[,1]+m10.1$Sol[,10]-(m10.1$Sol[,1]+m10.1$Sol[,12]))
summary(m10.1$Sol[,1]+m10.1$Sol[,10]-(m10.1$Sol[,1]+m10.1$Sol[,12]))
HPDinterval(m10.1$Sol[,1]+m10.1$Sol[,10]-(m10.1$Sol[,1]+m10.1$Sol[,12]))

# diff CNM & FW (13)

posterior.mode(m10.1$Sol[,1]+m10.1$Sol[,10]-(m10.1$Sol[,1]+m10.1$Sol[,13]))
summary(m10.1$Sol[,1]+m10.1$Sol[,10]-(m10.1$Sol[,1]+m10.1$Sol[,13]))
HPDinterval(m10.1$Sol[,1]+m10.1$Sol[,10]-(m10.1$Sol[,1]+m10.1$Sol[,13]))

# diff FC (11) & FD (12)

posterior.mode(m10.1$Sol[,1]+m10.1$Sol[,11]-(m10.1$Sol[,1]+m10.1$Sol[,12]))
summary(m10.1$Sol[,1]+m10.1$Sol[,11]-(m10.1$Sol[,1]+m10.1$Sol[,12]))
HPDinterval(m10.1$Sol[,1]+m10.1$Sol[,11]-(m10.1$Sol[,1]+m10.1$Sol[,12]))

# diff FC (11) & FW (13)

posterior.mode(m10.1$Sol[,1]+m10.1$Sol[,11]-(m10.1$Sol[,1]+m10.1$Sol[,13]))
summary(m10.1$Sol[,1]+m10.1$Sol[,11]-(m10.1$Sol[,1]+m10.1$Sol[,13]))
HPDinterval(m10.1$Sol[,1]+m10.1$Sol[,11]-(m10.1$Sol[,1]+m10.1$Sol[,13]))

# diff FD (12) & FW (13)

posterior.mode(m10.1$Sol[,1]+m10.1$Sol[,12]-(m10.1$Sol[,1]+m10.1$Sol[,13]))
summary(m10.1$Sol[,1]+m10.1$Sol[,12]-(m10.1$Sol[,1]+m10.1$Sol[,13]))
HPDinterval(m10.1$Sol[,1]+m10.1$Sol[,12]-(m10.1$Sol[,1]+m10.1$Sol[,13]))
##################
####################

H2<-m10.1$VCV[,"animal"]/(m10.1$VCV[,"animal"]+m10.1$VCV[,"Species_ID"]+m10.1$VCV[,"Study_ID"]+m10.1$VCV[,"units"]+s2)

posterior.mode(H2)
summary(H2)
HPDinterval(H2)

# getting s2

W<-1/c_all_varlogHR

s2<-sum(W*(length(W)-1))/(sum(W)^2-sum(W^2))

I2<-100*(m10.1$VCV[,"animal"]+m10.1$VCV[,"Species_ID"]+m10.1$VCV[,"Study_ID"])/(m10.1$VCV[,"animal"]+m10.1$VCV[,"Species_ID"]+m10.1$VCV[,"Study_ID"]+m10.1$VCV[,"units"]+s2)

posterior.mode(I2)
summary(I2)
HPDinterval(I2)

### getting residuals and plotting funnel plots....


m10.1$Random$formula<-update(m10.1$Random$formula, ~.+leg(mev, -1, FALSE):units)

#Pred1<-predict(m10.1)

Pred<-predict(m10.1, marginal=~leg(mev, -1, FALSE):units)
Prec<-sqrt(1/c_all_varlogHR)
ES<-c_all_logHR-predict(m10.1, marginal=~leg(mev, -1, FALSE):units)

plot(ES*sqrt(1/c_all_varlogHR),sqrt(1/c_all_varlogHR))
plot(sqrt(1/c_all_varlogHR),ES)

#library(arm)

model<-lm(ES~I(sqrt(1/c_all_varlogHR)))
summary(model)

##########################
##Model 11 (Model 10)
##########################


prior11<-list(R=list(V=1E-10,nu=-1),G=list(G1=list(V=1E-10,nu=-1),G2=list(V=1E-10,nu=-1),G3=list(V=1E-10,nu=-1)))

m11.1<-MCMCglmm(c_all_logHR~Sex+Model+c_nCR+I(c_nCR^2)+Prots+I(Prots^2)+AL,random=~animal+Species_ID+Study_ID,family="gaussian",data=Data,pedigree=Tree,verbose=F,mev=c_all_varlogHR,scale=F,nitt=100000,thin=50,burnin=50000,prior=prior11,pr=T)


m11.2<-MCMCglmm(c_all_logHR~Sex+Model+c_nCR+I(c_nCR^2)+Prots+I(Prots^2)+AL,random=~animal+Species_ID+Study_ID,family="gaussian",data=Data,pedigree=Tree,verbose=F,mev=c_all_varlogHR,scale=F,nitt=100000,thin=50,burnin=50000,prior=prior11,pr=T)


m11.3<-MCMCglmm(c_all_logHR~Sex+Model+c_nCR+I(c_nCR^2)+Prots+I(Prots^2)+AL,random=~animal+Species_ID+Study_ID,family="gaussian",data=Data,pedigree=Tree,verbose=F,mev=c_all_varlogHR,scale=F,nitt=100000,thin=50,burnin=50000,prior=prior11,pr=T)


save(m11.1,m11.2,m11.3,file="model11.RData")

#############################

library(MCMCglmm)

load("model11.RData")

gelman.diag(list(m11.1$Sol[,1:12],m11.2$Sol[,1:12],m11.3$Sol[,1:12]))
gelman.diag(list(m11.1$VCV[,c(1:3,5)],m11.2$VCV[,c(1:3,5)],m11.3$VCV[,c(1:3,5)]))
gelman.diag(list(m11.1$Deviance,m11.2$Deviance,m11.3$Deviance))

m11.1$DIC
m11.2$DIC
m11.3$DIC


# Fixed effects
posterior.mode(m11.1$Sol[,1:12])
summary(m11.1$Sol[,1:12])
HPDinterval(m11.1$Sol[,1:12])
plot(m11.1$Sol[,1:12])

# Random effects
posterior.mode(m11.1$VCV)
summary(m11.1$VCV)
HPDinterval(m11.1$VCV)
plot(m11.1$VCV)

H2<-m11.1$VCV[,"animal"]/(m11.1$VCV[,"animal"]+m11.1$VCV[,"Species_ID"]+m11.1$VCV[,"Study_ID"]+m11.1$VCV[,"units"]+s2)

posterior.mode(H2)
summary(H2)
HPDinterval(H2)

# getting s2

W<-1/c_all_varlogHR

s2<-sum(W*(length(W)-1))/(sum(W)^2-sum(W^2))

I2<-100*(m11.1$VCV[,"animal"]+m11.1$VCV[,"Species_ID"]+m11.1$VCV[,"Study_ID"])/(m11.1$VCV[,"animal"]+m11.1$VCV[,"Species_ID"]+m11.1$VCV[,"Study_ID"]+m11.1$VCV[,"units"]+s2)

posterior.mode(I2)
summary(I2)
HPDinterval(I2)


#########################
##Model 12 (Model 13)
##########################


prior12<-list(R=list(V=1E-10,nu=-1),G=list(G1=list(V=1E-10,nu=-1),G2=list(V=1E-10,nu=-1),G3=list(V=1E-10,nu=-1)))

m12.1<-MCMCglmm(c_all_logHR~Sex+Model+c_nCR+I(c_nCR^2)+Prots+I(Prots^2)+Intake,random=~animal+Species_ID+Study_ID,family="gaussian",data=Data,pedigree=Tree,verbose=F,mev=c_all_varlogHR,scale=F,nitt=100000,thin=50,burnin=50000,prior=prior12,pr=T)


m12.2<-MCMCglmm(c_all_logHR~Sex+Model+c_nCR+I(c_nCR^2)+Prots+I(Prots^2)+Intake,random=~animal+Species_ID+Study_ID,family="gaussian",data=Data,pedigree=Tree,verbose=F,mev=c_all_varlogHR,scale=F,nitt=100000,thin=50,burnin=50000,prior=prior12,pr=T)


m12.3<-MCMCglmm(c_all_logHR~Sex+Model+c_nCR+I(c_nCR^2)+Prots+I(Prots^2)+Intake,random=~animal+Species_ID+Study_ID,family="gaussian",data=Data,pedigree=Tree,verbose=F,mev=c_all_varlogHR,scale=F,nitt=100000,thin=50,burnin=50000,prior=prior12,pr=T)


save(m12.1,m12.2,m12.3,file="model12.RData")

#############################

library(MCMCglmm)

load("model12.RData")

gelman.diag(list(m12.1$Sol[,1:12],m12.2$Sol[,1:12],m12.3$Sol[,1:12]))
gelman.diag(list(m12.1$VCV[,c(1:3,5)],m12.2$VCV[,c(1:3,5)],m12.3$VCV[,c(1:3,5)]))
gelman.diag(list(m12.1$Deviance,m12.2$Deviance,m12.3$Deviance))

m12.1$DIC
m12.2$DIC
m12.3$DIC


# Fixed effects
posterior.mode(m12.1$Sol[,1:12])
summary(m12.1$Sol[,1:12])
HPDinterval(m12.1$Sol[,1:12])
plot(m12.1$Sol[,1:12])

# Random effects
posterior.mode(m12.1$VCV)
summary(m12.1$VCV)
HPDinterval(m12.1$VCV)
plot(m12.1$VCV)

H2<-m12.1$VCV[,"animal"]/(m12.1$VCV[,"animal"]+m12.1$VCV[,"Species_ID"]+m12.1$VCV[,"Study_ID"]+m12.1$VCV[,"units"]+s2)

posterior.mode(H2)
summary(H2)
HPDinterval(H2)

# getting s2

W<-1/c_all_varlogHR

s2<-sum(W*(length(W)-1))/(sum(W)^2-sum(W^2))

I2<-100*(m12.1$VCV[,"animal"]+m12.1$VCV[,"Species_ID"]+m12.1$VCV[,"Study_ID"])/(m12.1$VCV[,"animal"]+m12.1$VCV[,"Species_ID"]+m12.1$VCV[,"Study_ID"]+m12.1$VCV[,"units"]+s2)

posterior.mode(I2)
summary(I2)
HPDinterval(I2)


##########
##Model 13 (Model 11) # species longevity......
##########
prior13<-list(R=list(V=1E-10,nu=-1),G=list(G1=list(V=1E-10,nu=-1),G2=list(V=1E-10,nu=-1),G3=list(V=1E-10,nu=-1)))

m13.1<-MCMCglmm(c_all_logHR~Sex+Model+c_nCR+I(c_nCR^2)+Prots+I(Prots^2)+log(LT25contr),random=~animal+Species_ID+Study_ID,family="gaussian",data=Data,pedigree=Tree,verbose=F,mev=c_all_varlogHR,scale=F,nitt=100000,thin=50,burnin=50000,prior=prior13,pr=T)

m13.2<-MCMCglmm(c_all_logHR~Sex+Model+c_nCR+I(c_nCR^2)+Prots+I(Prots^2)+log(LT25contr),random=~animal+Species_ID+Study_ID,family="gaussian",data=Data,pedigree=Tree,verbose=F,mev=c_all_varlogHR,scale=F,nitt=100000,thin=50,burnin=50000,prior=prior13,pr=T)

m13.3<-MCMCglmm(c_all_logHR~Sex+Model+c_nCR+I(c_nCR^2)+Prots+I(Prots^2)+log(LT25contr),random=~animal+Species_ID+Study_ID,family="gaussian",data=Data,pedigree=Tree,verbose=F,mev=c_all_varlogHR,scale=F,nitt=100000,thin=50,burnin=50000,prior=prior13,pr=T)


save(m13.1,m13.2,m13.3,file="model13.RData")

#############################

library(MCMCglmm)

load("model13.RData")

gelman.diag(list(m13.1$Sol[,1:12],m13.2$Sol[,1:12],m13.3$Sol[,1:12]))
gelman.diag(list(m13.1$VCV[,c(1:3,5)],m13.2$VCV[,c(1:3,5)],m13.3$VCV[,c(1:3,5)]))
gelman.diag(list(m13.1$Deviance,m13.2$Deviance,m13.3$Deviance))

m13.1$DIC
m13.2$DIC
m13.3$DIC


# Fixed effects
posterior.mode(m13.1$Sol[,1:12])
summary(m13.1$Sol[,1:12])
HPDinterval(m13.1$Sol[,1:12])
plot(m13.1$Sol[,1:12])

# Random effects
posterior.mode(m13.1$VCV)
summary(m13.1$VCV)
HPDinterval(m13.1$VCV)
plot(m13.1$VCV)

H2<-m13.1$VCV[,"animal"]/(m13.1$VCV[,"animal"]+m13.1$VCV[,"Species_ID"]+m13.1$VCV[,"Study_ID"]+m13.1$VCV[,"units"]+s2)

posterior.mode(H2)
summary(H2)
HPDinterval(H2)

# getting s2

W<-1/c_all_varlogHR

s2<-sum(W*(length(W)-1))/(sum(W)^2-sum(W^2))

I2<-100*(m13.1$VCV[,"animal"]+m13.1$VCV[,"Species_ID"]+m13.1$VCV[,"Study_ID"])/(m13.1$VCV[,"animal"]+m13.1$VCV[,"Species_ID"]+m13.1$VCV[,"Study_ID"]+m13.1$VCV[,"units"]+s2)

posterior.mode(I2)
summary(I2)
HPDinterval(I2)


#######################
##Model 14 (Model 14)
#######################
# actual intake 

# CR becomes non-sig 

prior14<-list(R=list(V=1E-10,nu=-1),G=list(G1=list(V=1E-10,nu=-1),G2=list(V=1E-10,nu=-1),G3=list(V=1E-10,nu=-1)))

m14.1<-MCMCglmm(c_all_logHR~Sex+Model+c_CR+I(c_CR^2)+Prots+I(Prots^2),random=~animal+Species_ID+Study_ID,family="gaussian",data=Data,pedigree=Tree,verbose=F,mev=c_all_varlogHR,scale=F,nitt=100000,thin=50,burnin=50000,prior=prior14,pr=T)

m14.2<-MCMCglmm(c_all_logHR~Sex+Model+c_CR+I(c_CR^2)+Prots+I(Prots^2),random=~animal+Species_ID+Study_ID,family="gaussian",data=Data,pedigree=Tree,verbose=F,mev=c_all_varlogHR,scale=F,nitt=100000,thin=50,burnin=50000,prior=prior14,pr=T)

m14.3<-MCMCglmm(c_all_logHR~Sex+Model+c_CR+I(c_CR^2)+Prots+I(Prots^2),random=~animal+Species_ID+Study_ID,family="gaussian",data=Data,pedigree=Tree,verbose=F,mev=c_all_varlogHR,scale=F,nitt=100000,thin=50,burnin=50000,prior=prior14,pr=T)


save(m14.1,m14.2,m14.3,file="model14.RData")

#############################

library(MCMCglmm)

load("model14.RData")

gelman.diag(list(m14.1$Sol[,1:12],m14.2$Sol[,1:12],m14.3$Sol[,1:12]))
gelman.diag(list(m14.1$VCV[,c(1:3,5)],m14.2$VCV[,c(1:3,5)],m14.3$VCV[,c(1:3,5)]))
gelman.diag(list(m14.1$Deviance,m14.2$Deviance,m14.3$Deviance))

m14.1$DIC
m14.2$DIC
m14.3$DIC


# Fixed effects
posterior.mode(m14.1$Sol[,1:12])
summary(m14.1$Sol[,1:12])
HPDinterval(m14.1$Sol[,1:12])
plot(m14.1$Sol[,1:12])

# Random effects
posterior.mode(m14.1$VCV)
summary(m14.1$VCV)
HPDinterval(m14.1$VCV)
plot(m14.1$VCV)

H2<-m14.1$VCV[,"animal"]/(m14.1$VCV[,"animal"]+m14.1$VCV[,"Species_ID"]+m14.1$VCV[,"Study_ID"]+m14.1$VCV[,"units"]+s2)

posterior.mode(H2)
summary(H2)
HPDinterval(H2)

# getting s2

W<-1/c_all_varlogHR

s2<-sum(W*(length(W)-1))/(sum(W)^2-sum(W^2))

I2<-100*(m14.1$VCV[,"animal"]+m14.1$VCV[,"Species_ID"]+m14.1$VCV[,"Study_ID"])/(m14.1$VCV[,"animal"]+m14.1$VCV[,"Species_ID"]+m14.1$VCV[,"Study_ID"]+m14.1$VCV[,"units"]+s2)

posterior.mode(I2)
summary(I2)
HPDinterval(I2)

####
#### Model 15 (Model ??) - age on set
#######


prior15<-list(R=list(V=1E-10,nu=-1),G=list(G1=list(V=1E-10,nu=-1),G2=list(V=1E-10,nu=-1),G3=list(V=1E-10,nu=-1)))


m15.1<-MCMCglmm(c_all_logHR~Sex+Model+c_nCR+I(c_nCR^2)+Prots+I(Prots^2)+Page,random=~animal+Species_ID+Study_ID,family="gaussian",data=Data,pedigree=Tree,verbose=F,mev=c_all_varlogHR,scale=F,nitt=100000,thin=50,burnin=50000,prior=prior15,pr=T)

m15.2<-MCMCglmm(c_all_logHR~Sex+Model+c_nCR+I(c_nCR^2)+Prots+I(Prots^2)+Page,random=~animal+Species_ID+Study_ID,family="gaussian",data=Data,pedigree=Tree,verbose=F,mev=c_all_varlogHR,scale=F,nitt=100000,thin=50,burnin=50000,prior=prior15,pr=T)

m15.3<-MCMCglmm(c_all_logHR~Sex+Model+c_nCR+I(c_nCR^2)+Prots+I(Prots^2)+Page,random=~animal+Species_ID+Study_ID,family="gaussian",data=Data,pedigree=Tree,verbose=F,mev=c_all_varlogHR,scale=F,nitt=100000,thin=50,burnin=50000,prior=prior15,pr=T)


save(m15.1,m15.2,m15.3,file="model15.RData")

#############################

library(MCMCglmm)

load("model15.RData")

gelman.diag(list(m15.1$Sol[,1:12],m15.2$Sol[,1:12],m15.3$Sol[,1:12]))
gelman.diag(list(m15.1$VCV[,c(1:3,5)],m15.2$VCV[,c(1:3,5)],m15.3$VCV[,c(1:3,5)]))
gelman.diag(list(m15.1$Deviance,m15.2$Deviance,m15.3$Deviance))

m15.1$DIC
m15.2$DIC
m15.3$DIC


# Fixed effects
posterior.mode(m15.1$Sol[,1:12])
summary(m15.1$Sol[,1:12])
HPDinterval(m15.1$Sol[,1:12])
plot(m15.1$Sol[,1:12])

# Random effects
posterior.mode(m15.1$VCV)
summary(m15.1$VCV)
HPDinterval(m15.1$VCV)
plot(m15.1$VCV)

H2<-m15.1$VCV[,"animal"]/(m15.1$VCV[,"animal"]+m15.1$VCV[,"Species_ID"]+m15.1$VCV[,"Study_ID"]+m15.1$VCV[,"units"]+s2)

posterior.mode(H2)
summary(H2)
HPDinterval(H2)

# getting s2

W<-1/c_all_varlogHR

s2<-sum(W*(length(W)-1))/(sum(W)^2-sum(W^2))

I2<-100*(m15.1$VCV[,"animal"]+m15.1$VCV[,"Species_ID"]+m15.1$VCV[,"Study_ID"])/(m15.1$VCV[,"animal"]+m15.1$VCV[,"Species_ID"]+m15.1$VCV[,"Study_ID"]+m15.1$VCV[,"units"]+s2)

posterior.mode(I2)
summary(I2)
HPDinterval(I2)

##############
# The end....
##############

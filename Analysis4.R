# Analysis 4

# Egger's regression

# Mac
setwd("/Users/naksh50p/Dropbox/Project1/DR-longevity/Analysis")

rm(list=ls())
detach(Data)

library(MCMCglmm)

Data<-read.csv("DataAll.csv")
attach(Data)
str(Data)



# Egger test


prior1<-list(R=list(V=1E-10,nu=-1))


zES<-I(c_all_logHR/sqrt(c_all_varlogHR))




#Prec<-I(1/sqrt(c_all_varlogHR))


e0.1<-MCMCglmm(zES~I(1/sqrt(c_all_varlogHR)),family="gaussian",data=Data,verbose=F,nitt=50000,thin=25,burnin=25000)



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




######
# with all the data
# not to run

##########
prior1<-list(R=list(V=1E-10,nu=-1),G=list(G1=list(V=1E-10,nu=-1),G2=list(V=1E-10,nu=-1)))


zES<-I(c_all_logHR/sqrt(c_all_varlogHR))

Prec<-I(1/sqrt(c_all_varlogHR))

plot(Prec,zES)


e0.1<-MCMCglmm(zES~I(1/sqrt(c_all_varlogHR)),random=~Species_ID+Study_ID,family="gaussian",data=Data,verbose=F,nitt=50000,thin=25,burnin=25000,pr=T,prior=prior1)



###


# fixed effects
posterior.mode(e0.1$Sol[,1:2])
summary(e0.1$Sol[,1:2])
HPDinterval(e0.1$Sol[,1:2])

# random effects

posterior.mode(e0.1$VCV)
summary(e0.1$VCV)
HPDinterval(e0.1$VCV)

################
# makeing parts of datasets....

MaleLab<-which(Model==1 & Sex=="M") # 169
FemaleLab<-which(Model==1 & Sex=="F") # 98
MFLab<-which(Model==1 & Sex=="MF") # 5 - not doing it
NonLab<-which(Model==1 & Sex=="N") # 67 

MaleNon<-which(Model==0 & Sex=="M") # 28
FemaleNon<-which(Model==0 & Sex=="F") # 155
MFNon<-which(Model==0 & Sex=="MF") # 5
NonNon<-which(Model==0 & Sex=="N")  # 2 

# 8 kinds of data - 3 are too small (< 5 data points)

DataML<-Data[MaleLab,]
DataFL<-Data[FemaleLab,]
DataBL<-Data[MFLab,]
DataNL<-Data[NonLab,]

DataMN<-Data[MaleNon,]
DataFN<-Data[FemaleNon,]
DataBN<-Data[MFNon,]
DataNN<-Data[NonNon,]



########
# Female-Lab (e1)

prior1<-list(R=list(V=1E-10,nu=-1),G=list(G1=list(V=1E-10,nu=-1),G2=list(V=1E-10,nu=-1)))


# e1.1
e1.1<-MCMCglmm(c_all_logHR~I(1/sqrt(c_all_varlogHR)),random=~Species_ID+Study_ID,family="gaussian",data=DataFL,verbose=F,nitt=1000000,thin=50,burnin=950000,pr=T,prior=prior1)

# e1.2

e1.2<-MCMCglmm(c_all_logHR~I(1/sqrt(c_all_varlogHR)),random=~Species_ID+Study_ID,family="gaussian",data=DataFL,verbose=F,nitt=1000000,thin=50,burnin=950000,pr=T,prior=prior1)

# e1.3

e1.3<-MCMCglmm(c_all_logHR~I(1/sqrt(c_all_varlogHR)),random=~Species_ID+Study_ID,family="gaussian",data=DataFL,verbose=F,nitt=1000000,thin=50,burnin=950000,pr=T,prior=prior1)


save(e1.1,e1.2,e1.3,file="egger1.RData")

#########

prior2<-list(R=list(V=1E-10,nu=-1),G=list(G1=list(V=1E-10,nu=-1)))


# e1.1
e1.1<-MCMCglmm(c_all_logHR~I(1/sqrt(c_all_varlogHR)),random=~Study_ID,family="gaussian",data=DataFL,verbose=F,nitt=1000000,thin=50,burnin=950000,pr=T,prior=prior2)

# e1.2

e1.2<-MCMCglmm(c_all_logHR~I(1/sqrt(c_all_varlogHR)),random=~Study_ID,family="gaussian",data=DataFL,verbose=F,nitt=1000000,thin=50,burnin=950000,pr=T,prior=prior2)

# e1.3

e1.3<-MCMCglmm(c_all_logHR~I(1/sqrt(c_all_varlogHR)),random=~Study_ID,family="gaussian",data=DataFL,verbose=F,nitt=1000000,thin=50,burnin=950000,pr=T,prior=prior2)


save(e1.1,e1.2,e1.3,file="egger1.RData")


#####################
load("egger1.RData")
library(MCMCglmm)

gelman.diag(list(e1.1$Sol[,1:5],e1.2$Sol[,1:5],e1.3$Sol[,1:5]))
gelman.diag(list(e1.1$VCV[,c(1:3,5)],e1.2$VCV[,c(1:3,5)],e1.3$VCV[,c(1:3,5)]))
gelman.diag(list(e1.1$Deviance,e1.2$Deviance,e1.3$Deviance))

e1.1$DIC
e1.2$DIC
e1.3$DIC



###


# fixed effects
posterior.mode(e1.1$Sol[,1:2])
summary(e1.1$Sol[,1:2])
HPDinterval(e1.1$Sol[,1:2])
plot(e1.1$Sol[,1:2])

# random effects

posterior.mode(e1.1$VCV)
summary(e1.1$VCV)
HPDinterval(e1.1$VCV)
plot(e1.1$VCV)

#######
# Male-Lab (e2)

prior1<-list(R=list(V=1E-10,nu=-1),G=list(G1=list(V=1E-10,nu=-1),G2=list(V=1E-10,nu=-1)))

# e2.1

e2.1<-MCMCglmm(c_all_logHR~I(1/sqrt(c_all_varlogHR)),random=~Species_ID+Study_ID,family="gaussian",data=DataML,verbose=F,nitt=1000000,thin=50,burnin=950000,pr=T,prior=prior1)


# e2.2

e2.2<-MCMCglmm(c_all_logHR~I(1/sqrt(c_all_varlogHR)),random=~Species_ID+Study_ID,family="gaussian",data=DataML,verbose=F,nitt=1000000,thin=50,burnin=950000,pr=T,prior=prior1)

# e2.3

e2.3<-MCMCglmm(c_all_logHR~I(1/sqrt(c_all_varlogHR)),random=~Species_ID+Study_ID,family="gaussian",data=DataML,verbose=F,nitt=1000000,thin=50,burnin=950000,pr=T,prior=prior1)


save(e2.1,e2.2,e2.3,file="egger2.RData")


#########

prior2<-list(R=list(V=1E-10,nu=-1),G=list(G1=list(V=1E-10,nu=-1)))

# e2.1

e2.1<-MCMCglmm(c_all_logHR~I(1/sqrt(c_all_varlogHR)),random=~Study_ID,family="gaussian",data=DataML,verbose=F,nitt=1000000,thin=50,burnin=950000,pr=T,prior=prior2)


# e2.2

e2.2<-MCMCglmm(c_all_logHR~I(1/sqrt(c_all_varlogHR)),random=~Study_ID,family="gaussian",data=DataML,verbose=F,nitt=1000000,thin=50,burnin=950000,pr=T,prior=prior2)

# e2.3

e2.3<-MCMCglmm(c_all_logHR~I(1/sqrt(c_all_varlogHR)),random=~Study_ID,family="gaussian",data=DataML,verbose=F,nitt=1000000,thin=50,burnin=950000,pr=T,prior=prior2)


save(e2.1,e2.2,e2.3,file="egger2.RData")

#####################
load("egger2.RData")
library(MCMCglmm)

gelman.diag(list(e2.1$Sol[,1:5],e2.2$Sol[,1:5],e2.3$Sol[,1:5]))
gelman.diag(list(e2.1$VCV[,c(1:3,5)],e2.2$VCV[,c(1:3,5)],e2.3$VCV[,c(1:3,5)]))
gelman.diag(list(e2.1$Deviance,e2.2$Deviance,e2.3$Deviance))

e2.1$DIC
e2.2$DIC
e2.3$DIC

###


# fixed effects
posterior.mode(e2.1$Sol[,1:2])
summary(e2.1$Sol[,1:2])
HPDinterval(e2.1$Sol[,1:2])
plot(e2.1$Sol[,1:2])

# random effects

posterior.mode(e2.1$VCV)
summary(e2.1$VCV)
HPDinterval(e2.1$VCV)
plot(e2.1$VCV)


#######
# Both-Lab - N = 5 (too samll)

# Non-Lab (e3)

prior1<-list(R=list(V=1E-10,nu=-1),G=list(G1=list(V=1E-10,nu=-1),G2=list(V=1E-10,nu=-1)))

# e3.1

e3.1<-MCMCglmm(c_all_logHR~I(1/sqrt(c_all_varlogHR)),random=~Species_ID+Study_ID,family="gaussian",data=DataNL,verbose=F,nitt=1000000,thin=50,burnin=950000,pr=T,prior=prior1)


# e3.2

e3.2<-MCMCglmm(c_all_logHR~I(1/sqrt(c_all_varlogHR)),random=~Species_ID+Study_ID,family="gaussian",data=DataNL,verbose=F,nitt=1000000,thin=50,burnin=950000,pr=T,prior=prior1)

# e3.3


e3.3<-MCMCglmm(c_all_logHR~I(1/sqrt(c_all_varlogHR)),random=~Species_ID+Study_ID,family="gaussian",data=DataNL,verbose=F,nitt=1000000,thin=50,burnin=950000,pr=T,prior=prior1)


save(e3.1,e3.2,e3.3,file="egger3.RData")

##########

prior2<-list(R=list(V=1E-10,nu=-1),G=list(G1=list(V=1E-10,nu=-1)))
# e3.1

e3.1<-MCMCglmm(c_all_logHR~I(1/sqrt(c_all_varlogHR)),random=~Study_ID,family="gaussian",data=DataNL,verbose=F,nitt=1000000,thin=50,burnin=950000,pr=T,prior=prior2)


# e3.2

e3.2<-MCMCglmm(c_all_logHR~I(1/sqrt(c_all_varlogHR)),random=~Study_ID,family="gaussian",data=DataNL,verbose=F,nitt=1000000,thin=50,burnin=950000,pr=T,prior=prior2)

# e3.3


e3.3<-MCMCglmm(c_all_logHR~I(1/sqrt(c_all_varlogHR)),random=~Study_ID,family="gaussian",data=DataNL,verbose=F,nitt=1000000,thin=50,burnin=950000,pr=T,prior=prior2)


save(e3.1,e3.2,e3.3,file="egger3.RData")


#####################
load("egger3.RData")
library(MCMCglmm)

gelman.diag(list(e3.1$Sol[,1:5],e3.2$Sol[,1:5],e3.3$Sol[,1:5]))
gelman.diag(list(e3.1$VCV[,c(1:3,5)],e3.2$VCV[,c(1:3,5)],e3.3$VCV[,c(1:3,5)]))
gelman.diag(list(e3.1$Deviance,e3.2$Deviance,e3.3$Deviance))

e3.1$DIC
e3.2$DIC
e3.3$DIC

###


# fixed effects
posterior.mode(e3.1$Sol[,1:2])
summary(e3.1$Sol[,1:2])
HPDinterval(e3.1$Sol[,1:2])
plot(e3.1$Sol[,1:2])

# random effects

posterior.mode(e3.1$VCV)
summary(e3.1$VCV)
HPDinterval(e3.1$VCV)
plot(e3.1$VCV)

#############
# Female-Non (e4)

prior1<-list(R=list(V=1E-10,nu=-1),G=list(G1=list(V=1E-10,nu=-1),G2=list(V=1E-10,nu=-1)))


# e4.1

e4.1<-MCMCglmm(c_all_logHR~I(1/sqrt(c_all_varlogHR)),random=~Species_ID+Study_ID,family="gaussian",data=DataFN,verbose=F,nitt=1000000,thin=50,burnin=950000,pr=T,prior=prior1)

# e4.2

e4.2<-MCMCglmm(c_all_logHR~I(1/sqrt(c_all_varlogHR)),random=~Species_ID+Study_ID,family="gaussian",data=DataFN,verbose=F,nitt=1000000,thin=50,burnin=950000,pr=T,prior=prior1)

# e4.3


e4.3<-MCMCglmm(c_all_logHR~I(1/sqrt(c_all_varlogHR)),random=~Species_ID+Study_ID,family="gaussian",data=DataFN,verbose=F,nitt=1000000,thin=50,burnin=950000,pr=T,prior=prior1)


save(e4.1,e4.2,e4.3,file="egger4.RData")

######

prior2<-list(R=list(V=1E-10,nu=-1),G=list(G1=list(V=1E-10,nu=-1)))


# e4.1

e4.1<-MCMCglmm(c_all_logHR~I(1/sqrt(c_all_varlogHR)),random=~Study_ID,family="gaussian",data=DataFN,verbose=F,nitt=1000000,thin=50,burnin=950000,pr=T,prior=prior2)

# e4.2

e4.2<-MCMCglmm(c_all_logHR~I(1/sqrt(c_all_varlogHR)),random=~Study_ID,family="gaussian",data=DataFN,verbose=F,nitt=1000000,thin=50,burnin=950000,pr=T,prior=prior2)

# e4.3


e4.3<-MCMCglmm(c_all_logHR~I(1/sqrt(c_all_varlogHR)),random=~Study_ID,family="gaussian",data=DataFN,verbose=F,nitt=1000000,thin=50,burnin=950000,pr=T,prior=prior2)


save(e4.1,e4.2,e4.3,file="egger4.RData")


#####################
load("egger4.RData")
library(MCMCglmm)

gelman.diag(list(e4.1$Sol[,1:5],e4.2$Sol[,1:5],e4.3$Sol[,1:5]))
gelman.diag(list(e4.1$VCV[,c(1:3,5)],e4.2$VCV[,c(1:3,5)],e4.3$VCV[,c(1:3,5)]))
gelman.diag(list(e4.1$Deviance,e4.2$Deviance,e4.3$Deviance))

e4.1$DIC
e4.2$DIC
e4.3$DIC

###


# fixed effects
posterior.mode(e4.1$Sol[,1:2])
summary(e4.1$Sol[,1:2])
HPDinterval(e4.1$Sol[,1:2])
plot(e4.1$Sol[,1:2])

# random effects

posterior.mode(e4.1$VCV)
summary(e4.1$VCV)
HPDinterval(e4.1$VCV)
plot(e4.1$VCV)

###
plot(1/sqrt(DataFN$c_all_varlogHR),DataFN$c_all_logHR)



#############
# Male-Non (e5)

prior1<-list(R=list(V=1E-10,nu=-1),G=list(G1=list(V=1E-10,nu=-1),G2=list(V=1E-10,nu=-1)))

# e5.1

e5.1<-MCMCglmm(c_all_logHR~I(1/sqrt(c_all_varlogHR)),random=~Species_ID+Study_ID,family="gaussian",data=DataMN,verbose=F,nitt=1000000,thin=50,burnin=950000,pr=T,prior=prior1)

# e5.2

e5.2<-MCMCglmm(c_all_logHR~I(1/sqrt(c_all_varlogHR)),random=~Species_ID+Study_ID,family="gaussian",data=DataMN,verbose=F,nitt=1000000,thin=50,burnin=950000,pr=T,prior=prior1)

# e5.3

e5.3<-MCMCglmm(c_all_logHR~I(1/sqrt(c_all_varlogHR)),random=~Species_ID+Study_ID,family="gaussian",data=DataMN,verbose=F,nitt=1000000,thin=50,burnin=950000,pr=T,prior=prior1)


save(e5.1,e5.2,e5.3,file="egger5.RData")


#######

prior2<-list(R=list(V=1E-10,nu=-1),G=list(G1=list(V=1E-10,nu=-1)))

# e5.1

e5.1<-MCMCglmm(c_all_logHR~I(1/sqrt(c_all_varlogHR)),random=~Study_ID,family="gaussian",data=DataMN,verbose=F,nitt=1000000,thin=50,burnin=950000,pr=T,prior=prior2)

# e5.2

e5.2<-MCMCglmm(c_all_logHR~I(1/sqrt(c_all_varlogHR)),random=~Study_ID,family="gaussian",data=DataMN,verbose=F,nitt=1000000,thin=50,burnin=950000,pr=T,prior=prior2)

# e5.3

e5.3<-MCMCglmm(c_all_logHR~I(1/sqrt(c_all_varlogHR)),random=~Study_ID,family="gaussian",data=DataMN,verbose=F,nitt=1000000,thin=50,burnin=950000,pr=T,prior=prior2)


save(e5.1,e5.2,e5.3,file="egger5.RData")

#####################
load("egger5.RData")
library(MCMCglmm)

gelman.diag(list(e5.1$Sol[,1:5],e5.2$Sol[,1:5],e5.3$Sol[,1:5]))
gelman.diag(list(e5.1$VCV[,c(1:3,5)],e5.2$VCV[,c(1:3,5)],e5.3$VCV[,c(1:3,5)]))
gelman.diag(list(e5.1$Deviance,e5.2$Deviance,e5.3$Deviance))

e5.1$DIC
e5.2$DIC
e5.3$DIC

###


# fixed effects
posterior.mode(e5.1$Sol[,1:2])
summary(e5.1$Sol[,1:2])
HPDinterval(e5.1$Sol[,1:2])
plot(e5.1$Sol[,1:2])

# random effects

posterior.mode(e5.1$VCV)
summary(e5.1$VCV)
HPDinterval(e5.1$VCV)
plot(e5.1$VCV)


#######
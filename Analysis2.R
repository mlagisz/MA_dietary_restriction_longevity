library(MCMCglmm)
#library(ggplots)
#library(arm)

# Mac
setwd("/Users/naksh50p/Dropbox/Project1/DR-longevity/Analysis")

rm(list=ls())
detach(Data)

library(MCMCglmm)


#################

# Linax
setwd("/home/labadmin/Dropbox/Project1/DR-longevity/Analysis")

# Linax
library(MCMCglmm)
rm(list=ls())
detach(Data)



# conclusions: they are not clear patterns
# Mice M and F are sig different but not for Rats and Flies. 
# No sig across-spp interactions meaning that no species specific patterns in sex effects
# this stuff should be in the ESM.....

### Sex effects of separate 3 different species
################

#######
##Rats
#######
rm(list=ls())
detach(Data)
Data<-read.csv("DataR.csv")
attach(Data)
str(Data)

# it is good to have Strain in 

priorR<-list(R=list(V=1E-10,nu=-1),G=list(G1=list(V=1E-10,nu=-1),G2=list(V=1E-10,nu=-1)))

mSR.1<-MCMCglmm(c_all_logHR~Sex,random=~Strain+Study_ID,family="gaussian",data=Data,verbose=F,mev=c_all_varlogHR,scale=F,nitt=75000,thin=25,burnin=50000,prior=priorR)

mSR.2<-MCMCglmm(c_all_logHR~Sex,random=~Strain+Study_ID,family="gaussian",data=Data,verbose=F,mev=c_all_varlogHR,scale=F,nitt=75000,thin=25,burnin=50000,prior=priorR)

mSR.3<-MCMCglmm(c_all_logHR~Sex,random=~Strain+Study_ID,family="gaussian",data=Data,verbose=F,mev=c_all_varlogHR,scale=F,nitt=75000,thin=25,burnin=50000,prior=priorR)

gelman.diag(list(mSR.1$Sol[,1:2],mSR.2$Sol[,1:2],mSR.3$Sol[,1:2]))
gelman.diag(list(mSR.1$VCV[,c(1,2,4)],mSR.2$VCV[,c(1,2,4)],mSR.3$VCV[,c(1,2,4)]))
gelman.diag(list(mSR.1$Deviance,mSR.2$Deviance,mSR.3$Deviance))

mSR.1$DIC
mSR.2$DIC
mSR.3$DIC


save(mSR.1,mSR.2,mSR.3,file="modelSR.RData")

############
# going back to Mac

library(MCMCglmm)
load("modelSR.RData")

# Sol
posterior.mode(mSR.1$Sol)
summary(mSR.1$Sol)
HPDinterval(mSR.1$Sol)

# VCV
posterior.mode(mSR.1$VCV)
summary(mSR.1$VCV)
HPDinterval(mSR.1$VCV)


############
##Mouse
##############
rm(list=ls())
detach(Data)
Data<-read.csv("DataM.csv")
attach(Data)
str(Data)

priorM<-list(R=list(V=1E-10,nu=-1),G=list(G1=list(V=1E-10,nu=-1),G2=list(V=1E-10,nu=-1)))

mSM.1<-MCMCglmm(c_all_logHR~Sex,random=~Strain+Study_ID,family="gaussian",data=Data,verbose=F,mev=c_all_varlogHR,scale=F,nitt=75000,thin=25,burnin=50000,prior=priorM)


mSM.2<-MCMCglmm(c_all_logHR~Sex,random=~Strain+Study_ID,family="gaussian",data=Data,verbose=F,mev=c_all_varlogHR,scale=F,nitt=75000,thin=25,burnin=50000,prior=priorM)


mSM.3<-MCMCglmm(c_all_logHR~Sex,random=~Strain+Study_ID,family="gaussian",data=Data,verbose=F,mev=c_all_varlogHR,scale=F,nitt=75000,thin=25,burnin=50000,prior=priorM)

gelman.diag(list(mSM.1$Sol[,1:2],mSM.2$Sol[,1:2],mSM.3$Sol[,1:2]))
gelman.diag(list(mSM.1$VCV[,c(1,2,4)],mSM.2$VCV[,c(1,2,4)],mSM.3$VCV[,c(1,2,4)]))
gelman.diag(list(mSM.1$Deviance,mSM.2$Deviance,mSM.3$Deviance))

mSM.1$DIC
mSM.2$DIC
mSM.3$DIC


save(mSM.1,mSM.2,mSM.3,file="modelSM.RData")


############
# going back to Mac

library(MCMCglmm)
load("modelSM.RData")

# Sol
posterior.mode(mSM.1$Sol)
summary(mSM.1$Sol)
HPDinterval(mSM.1$Sol)

# VCV
posterior.mode(mSM.1$VCV)
summary(mSM.1$VCV)
HPDinterval(mSM.1$VCV)




#########
##Fly
########

rm(list=ls())
detach(Data)
Data<-read.csv("DataF.csv")
attach(Data)
str(Data)

priorF<-list(R=list(V=1E-10,nu=-1),G=list(G1=list(V=1E-10,nu=-1),G2=list(V=1E-10,nu=-1)))

mSF.1<-MCMCglmm(c_all_logHR~Sex,random=~Strain+Study_ID,family="gaussian",data=Data,verbose=F,mev=c_all_varlogHR,scale=F,nitt=75000,thin=25,burnin=50000,prior=priorF)

mSF.2<-MCMCglmm(c_all_logHR~Sex,random=~Strain+Study_ID,family="gaussian",data=Data,verbose=F,mev=c_all_varlogHR,scale=F,nitt=75000,thin=25,burnin=50000,prior=priorF)

mSF.3<-MCMCglmm(c_all_logHR~Sex,random=~Strain+Study_ID,family="gaussian",data=Data,verbose=F,mev=c_all_varlogHR,scale=F,nitt=75000,thin=25,burnin=50000,prior=priorF)

gelman.diag(list(mSF.1$Sol[,1:2],mSF.2$Sol[,1:2],mSF.3$Sol[,1:2]))
gelman.diag(list(mSF.1$VCV[,c(1,2,4)],mSF.2$VCV[,c(1,2,4)],mSF.3$VCV[,c(1,2,4)]))
gelman.diag(list(mSF.1$Deviance,mSF.2$Deviance,mSF.3$Deviance))

mSF.1$DIC
mSF.2$DIC
mSF.3$DIC


save(mSF.1,mSF.2,mSF.3,file="modelSF.RData")

############
# going back to Mac

library(MCMCglmm)
load("modelSF.RData")

# Sol
posterior.mode(mSF.1$Sol)
summary(mSF.1$Sol)
HPDinterval(mSF.1$Sol)

# VCV
posterior.mode(mSF.1$VCV)
summary(mSF.1$VCV)
HPDinterval(mSF.1$VCV)



#########
##Sex Mix 3 Spp
#########
rm(list=ls())
detach(Data)
Data<-read.csv("DataS.csv")
attach(Data)
str(Data)

priorS<-list(R=list(V=1E-10,nu=-1),G=list(G1=list(V=1E-10,nu=-1),G2=list(V=1E-10,nu=-1)))

mSS.1<-MCMCglmm(c_all_logHR~Sex*Species,random=~Strain+Study_ID,family="gaussian",data=Data,verbose=F,mev=c_all_varlogHR,scale=F,nitt=75000,thin=25,burnin=50000,prior=priorS)


## need to fix the code below
######################

mSS.2<-MCMCglmm(c_all_logHR~Sex,random=~Species_ID+Study_ID,family="gaussian",data=Data,verbose=F,mev=c_all_varlogHR,scale=F,nitt=75000,thin=25,burnin=50000,prior=priorSF)

mSS.3<-MCMCglmm(c_all_logHR~Sex,random=~Species_ID+Study_ID,family="gaussian",data=Data,verbose=F,mev=c_all_varlogHR,scale=F,nitt=75000,thin=25,burnin=50000,prior=priorS)


# to change....
gelman.diag(list(mSS.1$Sol[,1:2],mSS.2$Sol[,1:2],mSS.3$Sol[,1:2]))
gelman.diag(list(mSS.1$VCV[,c(1,2,4)],mSS.2$VCV[,c(1,2,4)],mSS.3$VCV[,c(1,2,4)]))
gelman.diag(list(mSS.1$Deviance,mSS.2$Deviance,mSS.3$Deviance))

mSS.1$DIC
mSS.2$DIC
mSS.3$DIC


save(mSS.1,mSS.2,mSS.3,file="modelSS.RData")

############
# going back to Mac

library(MCMCglmm)
load("modelSS.RData")

# Sol
posterior.mode(mSS.1$Sol)
summary(mSS.1$Sol)
HPDinterval(mSS.1$Sol)

# VCV
plot(mSS.1$VCV)
posterior.mode(mSS.1$VCV)
summary(mSS.1$VCV)
HPDinterval(mSS.1$VCV)

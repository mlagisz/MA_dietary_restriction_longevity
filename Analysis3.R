# 31 Dec 2010


par(mfrow=c(1,2))
plot(SlpES, sqrt(1/SlpV))
plot(IntES, sqrt(1/IntV))

# one data point of fruit fly has large precision
# should I do sensitivity test?????
# probably should - taking #4 out - I do not think it will change the overall result. 

##########################################################################################
# 1 Dec 2010

# need to modfiy to incorprate TreeM5 (not necessary)

# need to use N = 290 rather than N = 329 - with missing data in.....

#

####################
#All the model spp
#####################
# the overall intercept

#install.packages("orthopolynom")
library(MCMCglmm)
library(orthopolynom)
library(ape)

rm(list=ls())
detach(Data)

# Mac
setwd("/Users/naksh50p/Dropbox/Project1/DR-longevity/Analysis")

# Linex
setwd("/home/labadmin/Dropbox/Project1/DR-longevity/Analysis")

#Data<-read.csv("DataM5.csv")
Data<-read.csv("DataM5S.csv")
attach(Data)
str(Data)

Tree<-read.tree("TreeM5.tre")


# Slope

# pirors
prior3<-list(R=list(V=1E-10,nu=-1),G=list(G1=list(V=1E-10,nu=-1),G2=list(V=1E-10,nu=-1),G3=list(V=1E-10,nu=-1),G4=list(V=1E-10,nu=-1)))
prior4<-list(R=list(V=1,nu=0.002),G=list(G1=list(V=1,nu=0.002),G2=list(V=1,nu=0.002),G3=list(V=1,nu=0.002),G4=list(V=1,nu=0.002)))

prior5<-list(R=list(V=1E-10,nu=-1),G=list(G1=list(V=1E-10,nu=-1),G2=list(V=1E-10,nu=-1),G3=list(V=1E-10,nu=-1)))
prior6<-list(R=list(V=1,nu=0.002),G=list(G1=list(V=1,nu=0.002),G2=list(V=1,nu=0.002),G3=list(V=1,nu=0.002)))




# not very good
# with phylo
#mSLM5.1<-MCMCglmm(SlpES~1,random=~animal+Species_ID+Strain+Study_ID,family="gaussian",data=Data,verbose=F,mev=SlpV,nitt=200000,thin=50,burnin=150000,prior=prior3,pedigree=Tree,pr=TRUE)


# good one
# without phylo
mSLM5A.1<-MCMCglmm(SlpES~1,random=~Species_ID+Strain+Study_ID,family="gaussian",data=Data,verbose=F,mev=SlpV,nitt=1000000,thin=50,burnin=950000,prior=prior5,pr=TRUE)


mSLM5A.2<-MCMCglmm(SlpES~1,random=~Species_ID+Strain+Study_ID,family="gaussian",data=Data,verbose=F,mev=SlpV,nitt=1000000,thin=50,burnin=950000,prior=prior5,pr=TRUE)

mSLM5A.3<-MCMCglmm(SlpES~1,random=~Species_ID+Strain+Study_ID,family="gaussian",data=Data,verbose=F,mev=SlpV,nitt=1000000,thin=50,burnin=950000,prior=prior5,pr=TRUE)


gelman.diag(list(mSLM5A.1$Sol[,1:5],mSLM5A.2$Sol[,1:5],mSLM5A.3$Sol[,1:5]))
gelman.diag(list(mSLM5A.1$VCV[,c(1:3,5)],mSLM5A.2$VCV[,c(1:3,5)],mSLM5A.3$VCV[,c(1:3,5)]))
gelman.diag(list(mSLM5A.1$Deviance,mSLM5A.2$Deviance,mSLM5A.3$Deviance))

save(mSLM5A.1,mSLM5A.2,mSLM5A.3,file="modelSLM5A.RData")

#####################
load("modelSLM5A.RData")
library(MCMCglmm)

gelman.diag(list(mSLM5A.1$Sol[,1:5],mSLM5A.2$Sol[,1:5],mSLM5A.3$Sol[,1:5]))
gelman.diag(list(mSLM5A.1$VCV[,c(1:3,5)],mSLM5A.2$VCV[,c(1:3,5)],mSLM5A.3$VCV[,c(1:3,5)]))
gelman.diag(list(mSLM5A.1$Deviance,mSLM5A.2$Deviance,mSLM5A.3$Deviance))

mSLM5A.1$DIC
mSLM5A.2$DIC
mSLM5A.3$DIC

# DIC = 215.959.....
#########################

# Sol
posterior.mode(mSLM5A.1$Sol[,1])
summary(mSLM5A.1$Sol[,1])
HPDinterval(mSLM5A.1$Sol[,1])
plot(mSLM5A.1$Sol[,1])

# VCV
posterior.mode(mSLM5A.1$VCV)
summary(mSLM5A.1$VCV)
HPDinterval(mSLM5A.1$VCV)
plot(mSLM5A.1$VCV)

# getting s2

WS<-1/SlpV

s2S<-sum(WS*(length(WS)-1))/(sum(WS)^2-sum(WS^2))

# s2S added

I2<-100*(mSLM5A.1$VCV[,"Strain"]+mSLM5A.1$VCV[,"Species_ID"]+mSLM5A.1$VCV[,"Study_ID"])/(mSLM5A.1$VCV[,"Strain"]+mSLM5A.1$VCV[,"Species_ID"]+mSLM5A.1$VCV[,"Study_ID"]+mSLM5A.1$VCV[,"units"]+s2S)

posterior.mode(I2)
summary(I2)
HPDinterval(I2)

###################
# Intercept


# with phylo
#mINM5.1<-MCMCglmm(IntES~1,random=~animal+Species_ID+Strain+Study_ID,family="gaussian",data=Data,verbose=F,mev=IntV,nitt=200000,thin=50,burnin=150000,prior=prior3,pedigree=Tree,pr=TRUE)


# withoout phylo
mINM5A.1<-MCMCglmm(IntES~1,random=~Species_ID+Strain+Study_ID,family="gaussian",data=Data,verbose=F,mev=IntV,nitt=1000000,thin=50,burnin=950000, prior=prior5,pr=TRUE)

mINM5A.2<-MCMCglmm(IntES~1,random=~Species_ID+Strain+Study_ID,family="gaussian",data=Data,verbose=F,mev=IntV,nitt=1000000,thin=50,burnin=950000, prior=prior5,pr=TRUE)

mINM5A.3<-MCMCglmm(IntES~1,random=~Species_ID+Strain+Study_ID,family="gaussian",data=Data,verbose=F,mev=IntV,nitt=1000000,thin=50,burnin=950000, prior=prior5,pr=TRUE)


gelman.diag(list(mINM5A.1$Sol[,1:5],mINM5A.2$Sol[,1:5],mINM5A.3$Sol[,1:5]))
gelman.diag(list(mINM5A.1$VCV[,c(1:3,5)],mINM5A.2$VCV[,c(1:3,5)],mINM5A.3$VCV[,c(1:3,5)]))
gelman.diag(list(mINM5A.1$Deviance,mINM5A.2$Deviance,mINM5A.3$Deviance))

save(mINM5A.1,mINM5A.2,mINM5A.3,file="modelINM5A.RData")

##################
load("modelINM5A.RData")

gelman.diag(list(mINM5A.1$Sol[,1:5],mINM5A.2$Sol[,1:5],mINM5A.3$Sol[,1:5]))
gelman.diag(list(mINM5A.1$VCV[,c(1:3,5)],mINM5A.2$VCV[,c(1:3,5)],mINM5A.3$VCV[,c(1:3,5)]))
gelman.diag(list(mINM5A.1$Deviance,mINM5A.2$Deviance,mINM5A.3$Deviance))

# DIC = 718.8883

mINM5A.1$DIC
mINM5A.2$DIC
mINM5A.3$DIC


###########

# Sol
posterior.mode(mINM5A.1$Sol[,1])
summary(mINM5A.1$Sol[,1])
HPDinterval(mINM5A.1$Sol[,1])
plot(mINM5A.1$Sol[,1])

# VCV
posterior.mode(mINM5A.1$VCV)
summary(mINM5A.1$VCV)
HPDinterval(mINM5A.1$VCV)
plot(mINM5A.1$VCV)

# geeting s2I


WI<-1/IntV

s2I<-sum(WI*(length(WI)-1))/(sum(WI)^2-sum(WI^2))

# s2I added

I2<-100*(mINM5A.1$VCV[,"Strain"]+mINM5A.1$VCV[,"Species_ID"]+mINM5A.1$VCV[,"Study_ID"])/(mINM5A.1$VCV[,"Strain"]+mINM5A.1$VCV[,"Species_ID"]+mINM5A.1$VCV[,"Study_ID"]+mINM5A.1$VCV[,"units"]+s2I)

posterior.mode(I2)
summary(I2)
HPDinterval(I2)

#####################################################################

# 30 Nov 2010

####################
#All the model spp
#####################
# 5 intercepts for 5 different model organisms

#install.packages("orthopolynom")
library(MCMCglmm)
library(orthopolynom)

rm(list=ls())
detach(Data)

# Mac
setwd("/Users/naksh50p/Dropbox/Project1/DR-longevity/Analysis")

# Linex
setwd("/home/labadmin/Dropbox/Project1/DR-longevity/Analysis")

#Data<-read.csv("DataM5.csv")
Data<-read.csv("DataM5S.csv")
attach(Data)
str(Data)


# Slope

# pirors
prior3<-list(R=list(V=1E-10,nu=-1),G=list(G1=list(V=1E-10,nu=-1),G2=list(V=1E-10,nu=-1)))
prior4<-list(R=list(V=1,nu=0.002),G=list(G1=list(V=1,nu=0.002),G2=list(V=1,nu=0.002)))


mSLM5.1<-MCMCglmm(SlpES~Species_ID-1,random=~Strain+Study_ID,family="gaussian",data=Data,verbose=F,mev=SlpV,nitt=500000,thin=50,burnin=450000,prior=prior3)

mSLM5.2<-MCMCglmm(SlpES~Species_ID-1,random=~Strain+Study_ID,family="gaussian",data=Data,verbose=F,mev=SlpV,nitt=500000,thin=50,burnin=450000,prior=prior3)

mSLM5.3<-MCMCglmm(SlpES~Species_ID-1,random=~Strain+Study_ID,family="gaussian",data=Data,verbose=F,mev=SlpV,nitt=500000,thin=50,burnin=450000,prior=prior3)


gelman.diag(list(mSLM5.1$Sol[,1:5],mSLM5.2$Sol[,1:5],mSLM5.3$Sol[,1:5]))
gelman.diag(list(mSLM5.1$VCV[,c(1,2,4)],mSLM5.2$VCV[,c(1,2,4)],mSLM5.3$VCV[,c(1,2,4)]))
gelman.diag(list(mSLM5.1$Deviance,mSLM5.2$Deviance,mSLM5.3$Deviance))

save(mSLM5.1,mSLM5.2,mSLM5.3,file="modelSLM5.RData")
# mSLM5.1 is too good compared to the others - need to increase the nitt more....
#####################
load("modelSLM5.RData")
library(MCMCglmm)

gelman.diag(list(mSLM5.1$Sol[,1:5],mSLM5.2$Sol[,1:5],mSLM5.3$Sol[,1:5]))
gelman.diag(list(mSLM5.1$VCV[,c(1,2,4)],mSLM5.2$VCV[,c(1,2,4)],mSLM5.3$VCV[,c(1,2,4)]))
gelman.diag(list(mSLM5.1$Deviance,mSLM5.2$Deviance,mSLM5.3$Deviance))

mSLM5.1$DIC
mSLM5.2$DIC
mSLM5.3$DIC

# DIC = 218.2308
#########################

# Sol
posterior.mode(mSLM5.1$Sol)
summary(mSLM5.1$Sol)
HPDinterval(mSLM5.1$Sol)
plot(mSLM5.1$Sol)

# VCV
posterior.mode(mSLM5.1$VCV)
summary(mSLM5.1$VCV)
HPDinterval(mSLM5.1$VCV)
plot(mSLM5.1$VCV)


WS<-1/SlpV

s2S<-sum(WS*(length(WS)-1))/(sum(WS)^2-sum(WS^2))

# s2S added

I2<-100*(mSLM5.1$VCV[,"Strain"]+mSLM5.1$VCV[,"Study_ID"])/(mSLM5.1$VCV[,"Strain"]+mSLM5.1$VCV[,"Study_ID"]+mSLM5.1$VCV[,"units"]+s2S)

posterior.mode(I2)
summary(I2)
HPDinterval(I2)


###################
# Intercept

mINM5.1<-MCMCglmm(IntES~Species_ID-1,random=~Strain+Study_ID,family="gaussian",data=Data,verbose=F,mev=IntV,nitt=500000,thin=50,burnin=450000, prior=prior3)

mINM5.2<-MCMCglmm(IntES~Species_ID-1,random=~Strain+Study_ID,family="gaussian",data=Data,verbose=F,mev=IntV,nitt=500000,thin=50,burnin=450000, prior=prior3)

mINM5.3<-MCMCglmm(IntES~Species_ID-1,random=~Strain+Study_ID,family="gaussian",data=Data,verbose=F,mev=IntV,nitt=500000,thin=50,burnin=450000, prior=prior3)


gelman.diag(list(mINM5.1$Sol[,1:5],mINM5.2$Sol[,1:5],mINM5.3$Sol[,1:5]))
gelman.diag(list(mINM5.1$VCV[,c(1,2,4)],mINM5.2$VCV[,c(1,2,4)],mINM5.3$VCV[,c(1,2,4)]))
gelman.diag(list(mINM5.1$Deviance,mINM5.2$Deviance,mINM5.3$Deviance))

save(mINM5.1,mINM5.2,mINM5.3,file="modelINM5.RData")

##################
load("modelINM5.RData")

gelman.diag(list(mINM5.1$Sol[,1:5],mINM5.2$Sol[,1:5],mINM5.3$Sol[,1:5]))
gelman.diag(list(mINM5.1$VCV[,c(1,2,4)],mINM5.2$VCV[,c(1,2,4)],mINM5.3$VCV[,c(1,2,4)]))
gelman.diag(list(mINM5.1$Deviance,mINM5.2$Deviance,mINM5.3$Deviance))

# DIC = 718.8883

mINM5.1$DIC
mINM5.2$DIC
mINM5.3$DIC


###########

# Sol
posterior.mode(mINM5.1$Sol)
summary(mINM5.1$Sol)
HPDinterval(mINM5.1$Sol)
plot(mINM5.1$Sol)

# VCV
posterior.mode(mINM5.1$VCV)
summary(mINM5.1$VCV)
HPDinterval(mINM5.1$VCV)
plot(mINM5.1$VCV)



WI<-1/IntV

s2I<-sum(WI*(length(WI)-1))/(sum(WI)^2-sum(WI^2))

# s2I added

I2<-100*(mINM5.1$VCV[,"Strain"]+mINM5.1$VCV[,"Study_ID"])/(mINM5.1$VCV[,"Strain"]+mINM5.1$VCV[,"Study_ID"]+mINM5.1$VCV[,"units"]+s2I)

posterior.mode(I2)
summary(I2)
HPDinterval(I2)
########################################################
##################################################
##################################################
##################################################
##########################################
###below are trials...
###########################################
# 
# this is to look at slopes and intercepts and slopes of model species
# model species: 1) Rat, 2) Mice, 3) Fly, 4) Nematode and 5) Yeast
#

#############################
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

#################


#############
#####Rats
###############
rm(list=ls())
detach(Data)
Data<-read.csv("DataR.csv")
attach(Data)
str(Data)

library(MCMCglmm)

# Slope
# without Strain
# Nees more iteration
prior1<-list(R=list(V=1E-10,nu=-1),G=list(G1=list(V=1E-10,nu=-1)))
prior2<-list(R=list(V=1,nu=0.002),G=list(G1=list(V=1,nu=0.002)))

mSLR.1<-MCMCglmm(SlpES~1,random=~Study_ID,family="gaussian",data=Data,verbose=F,mev=SlpV,nitt=100000,thin=25,burnin=75000,prior=prior1)

# Sol
posterior.mode(mSLR.1$Sol)
summary(mSLR.1$Sol)
HPDinterval(mSLR.1$Sol)

#     Mean             SD       Naive SE Time-series SE 
#	-0.181896       0.056461       0.001785       0.002618

# VCV
posterior.mode(mSLR.1$VCV)
summary(mSLR.1$VCV)
HPDinterval(mSLR.1$VCV)
plot(mSLR.1$VCV)

# with Strain
prior3<-list(R=list(V=1E-10,nu=-1),G=list(G1=list(V=1E-10,nu=-1),G2=list(V=1E-10,nu=-1)))
prior4<-list(R=list(V=1,nu=0.002),G=list(G1=list(V=1,nu=0.002),G2=list(V=1,nu=0.002)))

mSLR.4<-MCMCglmm(SlpES~1,random=~Strain+Study_ID,family="gaussian",data=Data,verbose=F,mev=SlpV,nitt=100000,thin=25,burnin=75000,prior=prior3)

# Sol
posterior.mode(mSLR.4$Sol)
summary(mSLR.4$Sol)
HPDinterval(mSLR.4$Sol)

#     Mean             SD       Naive SE Time-series SE 
# -0.162576       0.079217       0.002505       0.006509

# VCV
posterior.mode(mSLR.4$VCV)
summary(mSLR.4$VCV)
HPDinterval(mSLR.4$VCV)
plot(mSLR.4$VCV)


###################
# Intercept

prior1<-list(R=list(V=1E-10,nu=-1),G=list(G1=list(V=1E-10,nu=-1)))
prior2<-list(R=list(V=1,nu=0.002),G=list(G1=list(V=1,nu=0.002)))

mINR.1<-MCMCglmm(Int25ES~1,random=~Study_ID,family="gaussian",data=Data,verbose=F,mev=Int25V,nitt=100000,thin=25,burnin=75000,prior=prior1)

# Sol
posterior.mode(mINR.1$Sol)
summary(mINR.1$Sol)
HPDinterval(mSLR.1$Sol)

#   Mean             SD       Naive SE Time-series SE 
#  -0.803710       0.072230       0.002284       0.003550
mINR
# VCV
posterior.mode(mINR.1$VCV)
summary(mINR.1$VCV)
HPDinterval(mINR.1$VCV)
plot(mINR.1$VCV)

# with Strain
prior3<-list(R=list(V=1E-10,nu=-1),G=list(G1=list(V=1E-10,nu=-1),G2=list(V=1E-10,nu=-1)))
prior4<-list(R=list(V=1,nu=0.002),G=list(G1=list(V=1,nu=0.002),G2=list(V=1,nu=0.002)))

mINR.4<-MCMCglmm(Int25ES~1,random=~Strain+Study_ID,family="gaussian",data=Data,verbose=F,mev=Int25V,nitt=100000,thin=25,burnin=75000,prior=prior3)

# Sol
posterior.mode(mINR.4$Sol)
summary(mINR.4$Sol)
HPDinterval(mINR.4$Sol)

#     Mean             SD       Naive SE Time-series SE 
#   -0.902498       0.119159       0.003768       0.004943

# VCV
posterior.mode(mINR.4$VCV)
summary(mINR.4$VCV)
HPDinterval(mINR.4$VCV)
plot(mINR.4$VCV)


##########
#####Mice
##########

rm(list=ls())
detach(Data)
Data<-read.csv("DataM.csv")
attach(Data)
str(Data)

# Slope
# without Strain
# Nees more iteration
prior1<-list(R=list(V=1E-10,nu=-1),G=list(G1=list(V=1E-10,nu=-1)))
prior2<-list(R=list(V=1,nu=0.002),G=list(G1=list(V=1,nu=0.002)))

mSLM.1<-MCMCglmm(SlpES~1,random=~Study_ID,family="gaussian",data=Data,verbose=F,mev=SlpV,nitt=100000,thin=25,burnin=75000,prior=prior1)

# Sol
posterior.mode(mSLM.1$Sol)
summary(mSLM.1$Sol)
HPDinterval(mSLM.1$Sol)

# 		Mean             SD       Naive SE Time-series SE 
#     -0.226503       0.049334       0.001560       0.001871

# VCV
posterior.mode(mSLM.1$VCV)
summary(mSLM.1$VCV)
HPDinterval(mSLM.1$VCV)
plot(mSLM.1$VCV)

# with Strain
prior3<-list(R=list(V=1E-10,nu=-1),G=list(G1=list(V=1E-10,nu=-1),G2=list(V=1E-10,nu=-1)))
prior4<-list(R=list(V=1,nu=0.002),G=list(G1=list(V=1,nu=0.002),G2=list(V=1,nu=0.002)))

mSLM.4<-MCMCglmm(SlpES~1,random=~Strain+Study_ID,family="gaussian",data=Data,verbose=F,mev=SlpV,nitt=100000,thin=25,burnin=75000,prior=prior3)

# Sol
posterior.mode(mSLM.4$Sol)
summary(mSLM.4$Sol)
HPDinterval(mSLM.4$Sol)

# 	 Mean             SD       Naive SE Time-series SE 
#	-0.222216       0.055975       0.001770       0.002397

# VCV
posterior.mode(mSLM.4$VCV)
summary(mSLM.4$VCV)
HPDinterval(mSLM.4$VCV)
plot(mSLM.4$VCV)

# Mice
###################
# Intercept

prior1<-list(R=list(V=1E-10,nu=-1),G=list(G1=list(V=1E-10,nu=-1)))
prior2<-list(R=list(V=1,nu=0.002),G=list(G1=list(V=1,nu=0.002)))

mINM.1<-MCMCglmm(Int25ES~1,random=~Study_ID,family="gaussian",data=Data,verbose=F,mev=Int25V,nitt=100000,thin=25,burnin=75000,prior=prior1)

# Sol
posterior.mode(mINM.1$Sol)
summary(mINM.1$Sol)
HPDinterval(mSLM.1$Sol)

#	Mean             SD       Naive SE Time-series SE 
#  -0.668527       0.080620       0.002549       0.002379

# VCV
posterior.mode(mINM.1$VCV)
summary(mINM.1$VCV)
HPDinterval(mINM.1$VCV)
plot(mINM.1$VCV)

# with Strain
prior3<-list(R=list(V=1E-10,nu=-1),G=list(G1=list(V=1E-10,nu=-1),G2=list(V=1E-10,nu=-1)))
prior4<-list(R=list(V=1,nu=0.002),G=list(G1=list(V=1,nu=0.002),G2=list(V=1,nu=0.002)))

mINM.4<-MCMCglmm(Int25ES~1,random=~Strain+Study_ID,family="gaussian",data=Data,verbose=F,mev=Int25V,nitt=100000,thin=25,burnin=75000,prior=prior3)

# Sol
posterior.mode(mINM.4$Sol)
summary(mINM.4$Sol)
HPDinterval(mINM.4$Sol)

#	Mean             SD       Naive SE Time-series SE 
# -0.588495       0.127903       0.004045       0.004708

# VCV
posterior.mode(mINM.4$VCV)
summary(mINM.4$VCV)
HPDinterval(mINM.4$VCV)
plot(mINM.4$VCV)


##########
#Fly
##########

rm(list=ls())
detach(Data)
Data<-read.csv("DataF.csv")
attach(Data)
str(Data)

# Slope
# without Strain
# Nees more iteration
prior1<-list(R=list(V=1E-10,nu=-1),G=list(G1=list(V=1E-10,nu=-1)))
prior2<-list(R=list(V=1,nu=0.002),G=list(G1=list(V=1,nu=0.002)))

mSLF.1<-MCMCglmm(SlpES~1,random=~Study_ID,family="gaussian",data=Data,verbose=F,mev=SlpV,nitt=100000,thin=25,burnin=75000,prior=prior1)

# Sol
posterior.mode(mSLF.1$Sol)
summary(mSLF.1$Sol)
HPDinterval(mSLF.1$Sol)

#     Mean             SD       Naive SE Time-series SE 
#-0.147171       0.070102       0.002217       0.002258

# VCV
posterior.mode(mSLF.1$VCV)
summary(mSLF.1$VCV)
HPDinterval(mSLF.1$VCV)
plot(mSLF.1$VCV)


# with Strain
prior3<-list(R=list(V=1E-10,nu=-1),G=list(G1=list(V=1E-10,nu=-1),G2=list(V=1E-10,nu=-1)))
prior4<-list(R=list(V=1,nu=0.002),G=list(G1=list(V=1,nu=0.002),G2=list(V=1,nu=0.002)))

mSLF.4<-MCMCglmm(SlpES~1,random=~Strain+Study_ID,family="gaussian",data=Data,verbose=F,mev=SlpV,nitt=100000,thin=25,burnin=75000,prior=prior3)

# Sol
posterior.mode(mSLF.4$Sol)
summary(mSLF.4$Sol)
HPDinterval(mSLF.4$Sol)

#     Mean             SD       Naive SE Time-series SE 
#   -0.153735       0.092610       0.002929       0.003193

# VCV
posterior.mode(mSLF.4$VCV)
summary(mSLF.4$VCV)
HPDinterval(mSLF.4$VCV)
plot(mSLF.4$VCV)


###################
# Intercept

prior1<-list(R=list(V=1E-10,nu=-1),G=list(G1=list(V=1E-10,nu=-1)))
prior2<-list(R=list(V=1,nu=-0.002),G=list(G1=list(V=1,nu=-0.002)))

mINF.1<-MCMCglmm(Int25ES~1,random=~Study_ID,family="gaussian",data=Data,verbose=F,mev=Int25V,nitt=100000,thin=25,burnin=75000,prior=prior1)

# Sol
posterior.mode(mINF.1$Sol)
summary(mINF.1$Sol)
HPDinterval(mINF.1$Sol)

#	Mean             SD       Naive SE Time-series SE 
#   -0.55700        0.14959        0.00473        0.00559



# VCV
posterior.mode(mINF.1$VCV)
summary(mINF.1$VCV)
HPDinterval(mINF.1$VCV)
plot(mINF.1$VCV)

# with Strain
prior3<-list(R=list(V=1E-10,nu=-1),G=list(G1=list(V=1E-10,nu=-1),G2=list(V=1E-10,nu=-1)))
prior4<-list(R=list(V=1,nu=0.002),G=list(G1=list(V=1,nu=0.002),G2=list(V=1,nu=0.002)))

mINF.4<-MCMCglmm(Int25ES~1,random=~Strain+Study_ID,family="gaussian",data=Data,verbose=F,mev=Int25V,nitt=100000,thin=25,burnin=75000,prior=prior3)

# Sol
posterior.mode(mINF.4$Sol)
summary(mINF.4$Sol)
HPDinterval(mINF.4$Sol)

#    Mean             SD       Naive SE Time-series SE 
#   -0.543744       0.197312       0.006240       0.005858

# VCV
posterior.mode(mINF.4$VCV)
summary(mINF.4$VCV)
HPDinterval(mINF.4$VCV)
plot(mINF.4$VCV)


#################

##########
#Nematode
##########
rm(list=ls())
detach(Data)
Data<-read.csv("DataN.csv")
attach(Data)
str(Data)


# Slope
# without Strain
# Nees more iteration
prior1<-list(R=list(V=1E-10,nu=-1),G=list(G1=list(V=1E-10,nu=-1)))
prior2<-list(R=list(V=1,nu=0.002),G=list(G1=list(V=1,nu=0.002)))

mSLN.1<-MCMCglmm(SlpES~1,random=~Study_ID,family="gaussian",data=Data,verbose=F,mev=SlpV,nitt=100000,thin=25,burnin=75000,prior=prior1)

# Sol
posterior.mode(mSLN.1$Sol)
summary(mSLN.1$Sol)
HPDinterval(mSLN.1$Sol)

#	Mean             SD       Naive SE Time-series SE 
#    -0.260637       0.077656       0.002456       0.002337 


# VCV
posterior.mode(mSLN.1$VCV)
summary(mSLN.1$VCV)
HPDinterval(mSLN.1$VCV)
plot(mSLN.1$VCV)

# with Strain
prior3<-list(R=list(V=1E-10,nu=-1),G=list(G1=list(V=1E-10,nu=-1),G2=list(V=1E-10,nu=-1)))
prior4<-list(R=list(V=1,nu=0.002),G=list(G1=list(V=1,nu=0.002),G2=list(V=1,nu=0.002)))

mSLN.4<-MCMCglmm(SlpES~1,random=~Strain+Study_ID,family="gaussian",data=Data,verbose=F,mev=SlpV,nitt=100000,thin=25,burnin=75000,prior=prior3)

# Sol
posterior.mode(mSLN.4$Sol)
summary(mSLN.4$Sol)
HPDinterval(mSLN.4$Sol)

#     Mean             SD       Naive SE Time-series SE 
#	-0.273260       0.166336       0.005260       0.005513

# VCV
posterior.mode(mSLN.4$VCV)
summary(mSLN.4$VCV)
HPDinterval(mSLM.4$VCV)
plot(mSLN.4$VCV)


# Nematode
###################
# Intercept

prior1<-list(R=list(V=1E-10,nu=-1),G=list(G1=list(V=1E-10,nu=-1)))
prior2<-list(R=list(V=1,nu=0.002),G=list(G1=list(V=1,nu=0.002)))

mINN.1<-MCMCglmm(Int25ES~1,random=~Study_ID,family="gaussian",data=Data,verbose=F,mev=Int25V,nitt=100000,thin=25,burnin=75000,prior=prior1)

# Sol
posterior.mode(mINN.1$Sol)
summary(mINN.1$Sol)
HPDinterval(mINN.1$Sol)

#	Mean             SD       Naive SE Time-series SE 
#     -1.022261       0.177653       0.005618       0.005951



# VCV
posterior.mode(mINN.1$VCV)
summary(mINN.1$VCV)
HPDinterval(mINN.1$VCV)
plot(mINN.1$VCV)

# with Strain
prior3<-list(R=list(V=1E-10,nu=-1),G=list(G1=list(V=1E-10,nu=-1),G2=list(V=1E-10,nu=-1)))
prior4<-list(R=list(V=1,nu=0.002),G=list(G1=list(V=1,nu=0.002),G2=list(V=1,nu=0.002)))

mINN.4<-MCMCglmm(Int25ES~1,random=~Strain+Study_ID,family="gaussian",data=Data,verbose=F,mev=Int25V,nitt=100000,thin=25,burnin=75000,prior=prior3)

# Sol
posterior.mode(mINN.4$Sol)
summary(mINN.4$Sol)
HPDinterval(mINN.4$Sol)

#Mean             SD       Naive SE Time-series SE 
#     -1.027620       0.264674       0.008370       0.008846

# VCV
posterior.mode(mINN.4$VCV)
summary(mINN.4$VCV)
HPDinterval(mINN.4$VCV)
plot(mINN.4$VCV)



###############
########Yeast
################

rm(list=ls())
detach(Data)
Data<-read.csv("DataY.csv")
attach(Data)
str(Data)


# Slope
# without Strain
# Nees more iteration
prior1<-list(R=list(V=1E-10,nu=-1),G=list(G1=list(V=1E-10,nu=-1)))
prior2<-list(R=list(V=1,nu=0.002),G=list(G1=list(V=1,nu=0.002)))

mSLY.1<-MCMCglmm(SlpES~1,random=~Study_ID,family="gaussian",data=Data,verbose=F,mev=SlpV,nitt=100000,thin=25,burnin=75000,prior=prior1)

# Sol
posterior.mode(mSLY.1$Sol)
summary(mSLY.1$Sol)
HPDinterval(mSLY.1$Sol)

#	Mean             SD       Naive SE Time-series SE 
#    0.018644       0.149387       0.004724       0.004767

# VCV
posterior.mode(mSLY.1$VCV)
summary(mSLY.1$VCV)
HPDinterval(mSLY.1$VCV)
plot(mSLY.1$VCV)

# with Strain
prior3<-list(R=list(V=1E-10,nu=-1),G=list(G1=list(V=1E-10,nu=-1),G2=list(V=1E-10,nu=-1)))
prior4<-list(R=list(V=1,nu=-0.002),G=list(G1=list(V=1,nu=0.002),G2=list(V=1,nu=0.002)))

mSLY.4<-MCMCglmm(SlpES~1,random=~Strain+Study_ID,family="gaussian",data=Data,verbose=F,mev=SlpV,nitt=100000,thin=25,burnin=75000,prior=prior3)

# Sol
posterior.mode(mSLY.4$Sol)
summary(mSLY.4$Sol)
HPDinterval(mSLY.4$Sol)

# 	Mean             SD       Naive SE Time-series SE 
#      0.020023       0.170795       0.005401       0.005261

# VCV
posterior.mode(mSLY.4$VCV)
summary(mSLY.4$VCV)
HPDinterval(mSLY.4$VCV)
plot(mSLY.4$VCV)

# Yeast
###################
# Intercept

prior1<-list(R=list(V=1E-10,nu=-1),G=list(G1=list(V=1E-10,nu=-1)))
prior2<-list(R=list(V=1,nu=0.002),G=list(G1=list(V=1,nu=0.002)))

mINY.1<-MCMCglmm(Int25ES~1,random=~Study_ID,family="gaussian",data=Data,verbose=F,mev=Int25V,nitt=100000,thin=25,burnin=75000,prior=prior1)

# Sol
posterior.mode(mINY.1$Sol)
summary(mINY.1$Sol)
HPDinterval(mINY.1$Sol)

#	Mean             SD       Naive SE Time-series SE 
#     -0.741013       0.245476       0.007763       0.007639 

# VCV
posterior.mode(mINY.1$VCV)
summary(mINY.1$VCV)
HPDinterval(mINY.1$VCV)
plot(mINY.1$VCV)

# with Strain
prior3<-list(R=list(V=1E-10,nu=-1),G=list(G1=list(V=1E-10,nu=-1),G2=list(V=1E-10,nu=-1)))
prior4<-list(R=list(V=1,nu=0.002),G=list(G1=list(V=1,nu=0.002),G2=list(V=1,nu=0.002)))

mINY.4<-MCMCglmm(Int25ES~1,random=~Strain+Study_ID,family="gaussian",data=Data,verbose=F,mev=Int25V,nitt=100000,thin=25,burnin=75000,prior=prior3)

# Sol
posterior.mode(mINY.4$Sol)
summary(mINY.4$Sol)
HPDinterval(mINY.4$Sol)

#   Mean             SD       Naive SE Time-series SE 
#	-0.79026        0.34351        0.01086        0.00953

# VCV
posterior.mode(mINY.4$VCV)
summary(mINY.4$VCV)
HPDinterval(mINY.4$VCV)
plot(mINY.4$VCV)


####################
#All the model spp
#####################

### here it is essential to fix "Strain"!!!! - fixed....
library(MCMCglmm)

rm(list=ls())
detach(Data)
Data<-read.csv("DataM5.csv")
attach(Data)
str(Data)


# Slope
# without Strain
# Nees more iteration
prior3<-list(R=list(V=1E-10,nu=-1),G=list(G1=list(V=1E-10,nu=-1),G2=list(V=1E-10,nu=-1)))
prior4<-list(R=list(V=1,nu=0.002),G=list(G1=list(V=1,nu=0.002),G2=list(V=1,nu=0.002)))

mSLM5.1<-MCMCglmm(SlpES~1,random=~Species_ID+Study_ID,family="gaussian",data=Data,verbose=F,mev=SlpV,nitt=100000,thin=25,burnin=75000,prior=prior3)

# Sol
posterior.mode(mSLM5.1$Sol)
summary(mSLM5.1$Sol)
HPDinterval(mSLM5.1$Sol)

#    Mean             SD       Naive SE Time-series SE 
#    -0.183827       0.077600       0.002454       0.002707

# VCV
posterior.mode(mSLM5.1$VCV)
summary(mSLM5.1$VCV)
HPDinterval(mSLM5.1$VCV)
plot(mSLM5.1$VCV)

# with Strain
prior5<-list(R=list(V=1E-10,nu=-1),G=list(G1=list(V=1E-10,nu=-1),G2=list(V=1E-10,nu=-1),G3=list(V=1E-10,nu=-1)))
prior6<-list(R=list(V=1,nu=0.002),G=list(G1=list(V=1,nu=0.002),G2=list(V=1,nu=0.002),G3=list(V=1,nu=0.002)))

#!!!!!!!!!!! strains have to be fixed!!!
# It has been fixed

#mSLM5.4<-MCMCglmm(SlpES~1,random=~Species_ID+Strain+Study_ID,family="gaussian",data=Data,verbose=F,mev=SlpV,nitt=100000,thin=25,burnin=75000,prior=prior5)

mSLM5.4<-MCMCglmm(SlpES~Species_ID-1,random=~Strain+Study_ID,family="gaussian",data=Data,verbose=F,mev=SlpV,nitt=100000,thin=25,burnin=75000,prior=prior3)

# Sol
posterior.mode(mSLM5.4$Sol)
summary(mSLM5.4$Sol)
HPDinterval(mSLM5.4$Sol)

#Species_IDSpp001 Species_IDSpp002 Species_IDSpp004 Species_IDSpp006 Species_IDSpp008 
#     -0.16415516      -0.25676596      -0.04524156      -0.21701118      -0.16341740

#Species_IDSpp001 -0.14665 0.07689 0.002432       0.002978
#Species_IDSpp002 -0.24813 0.06762 0.002138       0.002048
#Species_IDSpp004  0.00596 0.12140 0.003839       0.004126
#Species_IDSpp006 -0.26010 0.08950 0.002830       0.002782
#Species_IDSpp008 -0.16003 0.06385 0.002019       0.002095

# VCV
posterior.mode(mSLM5.4$VCV)
summary(mSLM5.4$VCV)
HPDinterval(mSLM5.4$VCV)
plot(mSLM5.4$VCV)

###################
# Intercept

prior3<-list(R=list(V=1E-10,nu=-1),G=list(G1=list(V=1E-10,nu=-1),G2=list(V=1E-10,nu=-1)))
prior4<-list(R=list(V=1,nu=0.002),G=list(G1=list(V=1,nu=0.002),G2=list(V=1,nu=0.002)))

mINM5.1<-MCMCglmm(Int25ES~1,random=~Species_ID+Study_ID,family="gaussian",data=Data,verbose=F,mev=Int25V,nitt=100000,thin=25,burnin=75000,prior=prior1)

# Sol
posterior.mode(mINM5.1$Sol)
summary(mINM5.1$Sol)
HPDinterval(mINM5.1$Sol)

#Mean             SD       Naive SE Time-series SE 
#-0.806267       0.073708       0.002331       0.003037
mINR
# VCV
posterior.mode(mINM5.1$VCV)
summary(mINM5.1$VCV)
HPDinterval(mINM5.1$VCV)
plot(mINM5.1$VCV)

# with Strain
prior5<-list(R=list(V=1E-10,nu=-1),G=list(G1=list(V=1E-10,nu=-1),G2=list(V=1E-10,nu=-1),G3=list(V=1E-10,nu=-1)))
prior6<-list(R=list(V=1,nu=0.002),G=list(G1=list(V=1,nu=-0.002),G2=list(V=1,nu=0.002),G3=list(V=1,nu=0.002)))

#mINM5.4<-MCMCglmm(Int25ES~1,random=~Strain+Study_ID,family="gaussian",data=Data,verbose=F,mev=Int25V,nitt=100000,thin=25,burnin=75000,prior=prior5)

mINM5.4<-MCMCglmm(IntES~Species_ID-1,random=~Strain+Study_ID,family="gaussian",data=Data,verbose=F,mev=IntV,nitt=100000,thin=25,burnin=75000,prior=prior3)

# Sol
posterior.mode(mINM5.4$Sol)
summary(mINM5.4$Sol)
HPDinterval(mINM5.4$Sol)

#Species_IDSpp001 Species_IDSpp002 Species_IDSpp004 Species_IDSpp006 Species_IDSpp008 
#      -0.6858567       -0.3929344       -0.6143786       -0.6965876       -0.3415701
                    Mean     SD Naive SE Time-series SE
#Species_IDSpp001 -0.7276 0.1892 0.005984       0.007682
#Species_IDSpp002 -0.4016 0.1563 0.004942       0.005231
#Species_IDSpp004 -0.6789 0.2852 0.009017       0.009006
#Species_IDSpp006 -0.6769 0.2221 0.007023       0.006605
#Species_IDSpp008 -0.3483 0.1658 0.005242       0.006214

# VCV
posterior.mode(mINM5.4$VCV)
summary(mINM5.4$VCV)
HPDinterval(mINM5.4$VCV)
plot(mINM5.4$VCV)


##########

#Extra

# quadratic effect +ve - I do not understand what this means.....

# trying to labtime - not sig plus the opposite direction from expected. 

Sub<-which(Labtime!="NA")

DataSub<-Data[Sub,]

model<-MCMCglmm(c_all_logHR~I(log(Labtime))+I(log(Labtime)^2)+Species_ID,random=~Strain+Study_ID,family="gaussian",data=DataSub,verbose=F,mev=DataSub$c_all_varlogHR,nitt=50000,thin=25,burnin=25000,prior=prior3,scale=F)

# trying random slope 
# bascally nothing is happening - more complex I guess....

prior7<-list(R=list(V=1,nu=0.002),G=list(G1=list(V=diag(2),nu=2)))

model<-MCMCglmm(c_all_logHR~I(log(Labtime)),random=~us(1+I(log(Labtime))):Species_ID,family="gaussian",data=DataSub,verbose=F,mev=DataSub$c_all_varlogHR,nitt=50000,thin=25,burnin=25000,prior=prior7,scale=F)


# Sol
posterior.mode(model$Sol)
summary(model$Sol)
HPDinterval(model$Sol)


# VCV
posterior.mode(model$VCV)
summary(model$VCV)
HPDinterval(model$VCV)
plot(model$VCV)

plot(log(Labtime),c_all_logHR,cex=1/(c_all_varlogHR)^(1/5))


#############

# moved from Anaalysis1.R
# this should be in the different file
# and it is in Analysis 2 - do not use below.......
##############################################################################
##############################################################################
##############################################################################
##############################################################################

## The five spp. 

# Mac
setwd("/Users/naksh50p/Dropbox/Projects/DR-longevity/Analysis")

rm(list=ls())
detach(Data)

library(MCMCglmm)

Data<-read.csv("Data.csv")
attach(Data)
str(Data)

# need to do sex effect interaction with 3 species rats and mice. 
# it is better to make separte data frames.....



###############
###########
##Rats
###########
SubR<-which(Species_ID=="Spp001")

DataR<-Data[SubR,]


priorR<-list(R=list(V=1E-10,nu=-1),G=list(G1=list(V=1E-10,nu=-1),G2=list(V=1E-10,nu=-1)))

mR.1<-MCMCglmm(c_all_logHR~Sex+c_nCR+I(c_nCR^2)+Prots+I(Prots^2),random=~Strain+Study_ID,family="gaussian",data=DataR,verbose=F,mev=c_all_varlogHR[SubR],scale=F,nitt=50000,thin=25,burnin=25000,prior=priorR)


summary(mR.1$Sol)

mR.1$DIC


mR.2<-MCMCglmm(c_all_logHR~Sex+c_nCR+Prots,random=~Strain+Study_ID,family="gaussian",data=DataR,verbose=F,mev=c_all_varlogHR[SubR],scale=F,nitt=50000,thin=25,burnin=25000,prior=priorR)

summary(mR.2$Sol)

mR.2$DIC

priorR1<-list(R=list(V=1E-10,nu=-1),G=list(G1=list(V=1E-10,nu=-1)))
priorR2<-list(R=list(V=0.002,nu=1),G=list(G1=list(V=0.002,nu=1)))

mR.3<-MCMCglmm(c_all_logHR~Sex+c_nCR,random=~Study_ID,family="gaussian",data=DataR,verbose=F,mev=c_all_varlogHR[SubR],scale=F,nitt=50000,thin=25,burnin=25000,prior=priorR1,singular.ok=T)

summary(mR.3$Sol)

mR.3$DIC

###########
##Mice
###########

SubM<-which(Species_ID=="Spp002")


#### Sex effects but better version in Analysis2.R

#############################################
###############################################################
# Sex Effect --- Mice Rats and Flies separate analysis
###############################################################################
################################################################
# Rats - using DataR
library(MCMCglmm)
rm(list=ls())
detach(Data)
Data<-read.csv("DataR.csv")
attach(Data)
str(Data)

priorR<-list(R=list(V=1E-10,nu=-1),G=list(G1=list(V=1E-10,nu=-1)))

mR.1<-MCMCglmm(c_all_logHR~Sex,random=~Study_ID,family="gaussian",data=Data,verbose=F,mev=c_all_varlogHR,scale=F,nitt=50000,thin=25,burnin=25000,prior=priorR,pr=T) # could take out pr = T

# Sex effect is not sig

plot(mR.1$Sol[,1:2])
summary(mR.1$Sol[,1:2])
posterior.mode(mR.1$Sol[,1:2])
HPDinterval(mR.1$Sol[,1:2])

plot(mR.1$VCV)
summary(mR.1$VCV)
posterior.mode(mR.1$VCV)
HPDinterval(mR.1$VCV)

###################################
# Mice - using DataM
library(MCMCglmm)
rm(list=ls())
detach(Data)
Data<-read.csv("DataM.csv")
attach(Data)
str(Data)

priorM<-list(R=list(V=1E-10,nu=-1),G=list(G1=list(V=1E-10,nu=-1)))

mM.1<-MCMCglmm(c_all_logHR~Sex,random=~Study_ID,family="gaussian",data=Data,verbose=F,mev=c_all_varlogHR,scale=F,nitt=50000,thin=25,burnin=25000,prior=priorM,pr=T) # could take out pr = T

# Sex effect is sig. 

plot(mM.1$Sol[,1:2])
summary(mM.1$Sol[,1:2])
posterior.mode(mM.1$Sol[,1:2])
HPDinterval(mM.1$Sol[,1:2])

plot(mM.1$VCV)
summary(mM.1$VCV)
posterior.mode(mM.1$VCV)
HPDinterval(mM.1$VCV)

################################
# Flies - using DataF
library(MCMCglmm)
rm(list=ls())
detach(Data)
Data<-read.csv("DataF.csv")
attach(Data)
str(Data)

priorF<-list(R=list(V=1E-10,nu=-1),G=list(G1=list(V=1E-10,nu=-1)))

mF.1<-MCMCglmm(c_all_logHR~Sex,random=~Study_ID,family="gaussian",data=Data,verbose=F,mev=c_all_varlogHR,scale=F,nitt=50000,thin=25,burnin=25000,prior=priorF,pr=T) # could take out pr = T

# Sex effect not sig. 

plot(mF.1$Sol[,1:2])
summary(mF.1$Sol[,1:2])
posterior.mode(mF.1$Sol[,1:2])
HPDinterval(mF.1$Sol[,1:2])

plot(mF.1$VCV)
summary(mF.1$VCV)
posterior.mode(mF.1$VCV)
HPDinterval(mF.1$VCV)


#################################
#################################
# Use DataS - sexual spp
library(MCMCglmm)
rm(list=ls())
detach(Data)
Data<-read.csv("DataS.csv")
attach(Data)
str(Data)

priorS<-list(R=list(V=1E-10,nu=-1),G=list(G1=list(V=1E-10,nu=-1)))

mS.1<-MCMCglmm(c_all_logHR~Sex*Species_ID,random=~Study_ID,family="gaussian",data=Data,verbose=F,mev=c_all_varlogHR,scale=F,nitt=50000,thin=25,burnin=25000,prior=priorS,pr=T) # could take out pr = T

#plot(mS.1$Sol[,1:6])
summary(mS.1$Sol[,1:6])
posterior.mode(mS.1$Sol[,1:6])
HPDinterval(mS.1$Sol[,1:6])

plot(mS.1$VCV)
summary(mS.1$VCV)
posterior.mode(mS.1$VCV)
HPDinterval(mS.1$VCV)


################################################
################
###################





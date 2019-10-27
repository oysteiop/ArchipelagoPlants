##################################################################################################
# ARCHIPELAGO PLANTS: OUTPUT OF COLONISATION MODEL
##################################################################################################

rm(list=ls())
gc()
memory.size()
memory.limit(64000)

#### - Read packages - ####
library(devtools)
library(BayesLogit)
library(Hmsc)
library(abind)
library(corrplot)
library(ggplot2)

#### - Load the fitted model - ####
setwd("F:/HelsinkiData23102019/archipelago/hmsc/Rcode/colonisation")
load("modelNS.RData")

# Extract posterior distribution
post=convertToCodaObject(m)

# Compute effective sample sizes and PSRFs
esBeta=effectiveSize(post$Beta)
esGamma=effectiveSize(post$Gamma)
esRho=effectiveSize(post$Rho)
esAlpha=effectiveSize(post$Alpha[[1]])
esEta=effectiveSize(post$Eta[[1]])
esLambda=effectiveSize(post$Lambda[[1]])
str(post$Omega[[1]])
set.seed(1)
ss=sample(1:229441, 1000, replace=F)
esOmega=effectiveSize(post$Omega[[1]][,ss])

psBeta=gelman.diag(post$Beta, multivariate=F)$psrf
psGamma=gelman.diag(post$Gamma, multivariate=F)$psrf
psRho=gelman.diag(post$Rho, multivariate=F)$psrf
psAlpha=gelman.diag(post$Alpha[[1]], multivariate=F)$psrf
psEta=gelman.diag(post$Eta[[1]], multivariate=F)$psrf
psLambda=gelman.diag(post$Lambda[[1]], multivariate=F)$psrf
psOmega=gelman.diag(post$Omega[[1]][,ss], multivariate=F)$psrf

mixing = list(esBeta, esGamma, esRho, esAlpha, esEta, esLambda, esOmega,
              psBeta, psGamma, psRho, psAlpha, psEta, psLambda, psOmega)
#save(mixing, file="mixing.RData")

load(file="mixing.RData")

summary(mixing[[1]]) # Beta
summary(mixing[[2]]) # Gamma
summary(mixing[[3]]) # Rho
summary(mixing[[7]]) # Omega

summary(mixing[[8]]) # Beta
summary(mixing[[9]]) # Gamma
summary(mixing[[10]]) # Rho
summary(mixing[[14]]) # Omega

# Produce posterior trace plots
plot(post$Rho)
summary(post$Rho[[1]]) #Phylogenetic signal

pdf("posteriorplots/betapost.pdf") #Regression coefficients
plot(post$Beta[,1:200])
dev.off()

pdf("posteriorplots/gammapost.pdf") #Trait effects
plot(post$Gamma)
dev.off()

pdf("posteriorplots/omegapost.pdf") #Species associations
plot(post$Omega[[1]][,1:200])
dev.off()

#### - Evaluate model fit - ####
#load("CrossVal5PredY.RData")
predY = computePredictedValues(m)
predYm=apply(simplify2array(predY), 1:2, mean) #Posterior mean
#save(predYm, file="predYm.RData")

MF=evaluateModelFit(m, predY)
#save(MF, file="MF.RData")

load(file="MF.RData")

AUC = MF$AUC
R2 = MF$TjurR2

mean(R2,na.rm=T) #0.17
range(R2,na.rm=T) #0 - 0.63
mean(AUC,na.rm=T) #0.89
range(AUC,na.rm=T) #.60 - 1.00

# Explanatory power at island level
load(file="predYm.RData")

tmp=(m$Y>(-Inf))*1
predYm2=predYm*tmp
plot(rowSums(m$Y, na.rm=T), rowSums(predYm2, na.rm=T))
lines(0:100,0:100)
cor(rowSums(m$Y, na.rm=T), rowSums(predYm2, na.rm=T))^2

# Plot Tjur r^2 vs species prevalence
pdf("tjur_vs_prev.pdf", height=4, width=4, family="Times")
par(mfrow=c(1,1),mar=c(4,5,2,1))
plot(prev,R2,las=1,pch=16,col="grey",cex=.8,main=paste("Repeated measures model: Mean = ",signif(mean(R2,na.rm=T),2),".", sep=""),ylim = c(0,1),xlab = "",ylab=expression(paste("Coefficient of discrimination (Tjur's",r^2,")")))
mtext("Species prevalence",1,line=2.5)
dev.off()

#### - Compute and plot variance partitioning - ####
group = c(1,1,2,3)
groupnames = c(m$covNames[-1])
groupnames
VP = computeVariancePartitioning(m, group = group, groupnames = groupnames)
#save(VP, file="VP.RData")

load(file="VP.RData")

str(VP)

sums=apply(VP$vals[1:3,],2,sum)
VP$vals=VP$vals[,rev(order(sums))]

pdf("plots/varpartCol.pdf",height=5,width=60)
plotVariancePartitioning(m, VP = VP)
dev.off()

### - Association networks ordered by taxonomy - ####
OmegaCorCol = computeAssociations(m)
#save(OmegaCorCol, file="OmegaCorCol.RData")

load(file="OmegaCorCol.RData")

tree=m$phyloTree
tree=untangle(tree,"read.tree")

orderC=m$C[tree$tip.label,tree$tip.label]
orderC[1:5,1:5]
orderOmega=OmegaCorCol[[1]]$mean[tree$tip.label,tree$tip.label]
orderOmega[1:5,1:5]
cor(c(orderC),c(orderOmega)) #0.014

plotOrderCol=match(tree$tip.label,colnames(OmegaCorCol[[1]]$mean))

#### - Cross-validation - ####
a=Sys.time()
partition=createPartition(m, nfolds=2, column=1)
predY_CV2 = computePredictedValues(m, partition=partition, updater=list(GammaEta=FALSE))
MF_CV2=evaluateModelFit(m, predY_CV2)
Sys.time()-a

save(MF_CV2,file="MF_CV2.Rdata") #2fold
save(predY_CV2,file="crossvalY.Rdata") #2fold

load("crossvalY.Rdata") #5fold

R2 = MF_CV2$TjurR2
AUC = MF_CV2$AUC

mean(R2,na.rm=T)
mean(AUC,na.rm=T)
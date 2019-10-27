##################################################################################################
# ARCHIPELAGO PLANTS: OUTPUT OF REPEATED MEASURES MODEL
##################################################################################################

rm(list=ls())
gc()
memory.size()
memory.limit(64000)

#### - Read packages - ####
#library(devtools)
library(BayesLogit)
library(Hmsc)
#library(abind)
#library(corrplot)
#library(ggplot2)
library(parallel)

#### - Load the fitted model - ####
setwd("F:/HelsinkiData23102019/archipelago/hmsc/Rcode/repeatedmeasures")
load("modelNS.RData")

# Extract posterior distribution
post=convertToCodaObject(m)

# Compute effective sample sizes and PSRFs
esBeta = effectiveSize(post$Beta)
esGamma = effectiveSize(post$Gamma)
esAlpha = effectiveSize(post$Alpha[[1]])
esEta = effectiveSize(post$Eta[[1]])
esEta2 = effectiveSize(post$Eta[[2]])
esLambda = effectiveSize(post$Lambda[[1]])
esLambda2 = effectiveSize(post$Lambda[[2]])
str(post$Omega[[1]])
set.seed(1)
ss = sample(1:342225, 1000, replace=F)
esOmega = effectiveSize(post$Omega[[1]][,ss])
esOmega2 = effectiveSize(post$Omega[[2]][,ss])

psBeta = gelman.diag(post$Beta, multivariate=F)$psrf
psGamma = gelman.diag(post$Gamma, multivariate=F)$psrf
psAlpha = gelman.diag(post$Alpha[[1]], multivariate=F)$psrf
psEta = gelman.diag(post$Eta[[1]], multivariate=F)$psrf
psEta2 = gelman.diag(post$Eta[[2]], multivariate=F)$psrf
psLambda = gelman.diag(post$Lambda[[1]], multivariate=F)$psrf
psLambda2 = gelman.diag(post$Lambda[[2]], multivariate=F)$psrf
psOmega = gelman.diag(post$Omega[[1]][,ss], multivariate=F)$psrf
psOmega2 = gelman.diag(post$Omega[[2]][,ss], multivariate=F)$psrf

mixing = list(esBeta, esGamma, esAlpha, esEta, esEta2, esLambda, esLambda2, esOmega, esOmega2,
              psBeta, psGamma, psAlpha, psEta, psEta2, psLambda, psLambda2, psOmega, psOmega2)
#save(mixing, file="mixing.RData")

load(file="mixing.RData")

summary(mixing[[1]]) # Beta
summary(mixing[[2]]) # Gamma
summary(mixing[[8]]) # Omega 1 (Island)
summary(mixing[[9]]) # Omega 2 (SU)

summary(mixing[[10]]) # Beta
summary(mixing[[11]]) # Gamma
summary(mixing[[17]]) # Omega 1 (Island)
summary(mixing[[18]]) # Omega 2 (SU)

# Produce posterior trace plots
pdf("posteriorplots/betapost.pdf") #Regression parameters
plot(post$Beta[,1:200])
dev.off()

pdf("posteriorplots/gammapost.pdf") #Trait effects
plot(post$Gamma)
dev.off()

pdf("posteriorplots/omegapost.pdf") #Species associations
plot(post$Omega[[1]][,1:200])
dev.off()

# Evaluate model fit
predY = computePredictedValues(m)
predYm = apply(simplify2array(predY), 1:2, mean) #Posterior mean
#save(predYm, file="predYm.RData")

MF = evaluateModelFit(m, predY)
#save(MF, file="MF.RData")

load(file="MF.RData")

AUC = MF$AUC
R2 = MF$TjurR2

mean(R2,na.rm=T) #0.30
range(R2,na.rm=T) #0 - 0.71
mean(AUC,na.rm=T) #0.94
range(AUC,na.rm=T) #0.74 - 1.00

# Plot Tjur r^2 vs species prevalence
prev=colSums(m$Y)/(m$ny)

pdf("plots/tjur_vs_prev.pdf", height=4, width=4, family="Times")
par(mfrow=c(1,1),mar=c(4,5,2,1))
plot(prev,MF$TjurR2,las=1,pch=16,col="grey",cex=.8,main=paste("Repeated measures model: Mean = ",signif(mean(MF$TjurR2,na.rm=T),2),".", sep=""),ylim = c(0,1),xlab = "",ylab=expression(paste("Coefficient of discrimination (Tjur's",r^2,")")))
mtext("Species prevalence",1,line=2.5)
dev.off()

# Explanatory power at island level
load(file="predYm.RData")

plot(rowSums(m$Y),rowSums(predYm)) #Species richness per island per sampling time
round(cor(rowSums(m$Y),rowSums(predYm))^2,2)
#full: 1.00

plot(colSums(m$Y),colSums(predYm)) #Occurrences per species
round(cor(colSums(m$Y),colSums(predYm))^2,2)
#full: 1.0

#### - Compute and plot variance partitioning - ####
m$covNames
group = c(1,1,2,3,4)
groupnames = c(m$covNames[-1])
groupnames

VP = computeVariancePartitioning(m, group = group, groupnames = groupnames)
#save(VP, file="VP.RData")

load(file="VP.RData")

str(VP) #Trait r^2 = 9.5%

pdf("plots/varpart.pdf", height=5, width=60)
plotVariancePartitioning(m, VP = VP)
dev.off()

# Trait effects
gammaPost = getPostEstimate(m, "Gamma")
gammaPost
plotGamma(m, post = gammaPost, param = "Mean", covOrder = "Vector", covVector = c(2:5))
summary(post$Gamma)

#### - Plot association networks - ####
OmegaCor = computeAssociations(m)
#save(OmegaCor, file="OmegaCor.RData")

load(file="OmegaCor.RData")

pdf("plots/rmCor.pdf",width=4,height=4)
supportLevel = 0.75
for (r in 1:m$nr){
  plotOrder = corrMatOrder(OmegaCor[[r]]$mean,order="AOE")
  toPlot = ((OmegaCor[[r]]$support>supportLevel) + (OmegaCor[[r]]$support<(1-supportLevel))>0)*OmegaCor[[r]]$mean
  corrplot(toPlot[plotOrder,plotOrder], type="lower",tl.pos="n",method = "color", col=colorRampPalette(c("blue","white","red"))(200),title=paste("random effect level:",m$levelNames[r]), mar=c(0,0,1,0))
}
dev.off()

#### - Plot map of change in species richness - ####
load(file="predYm.RData")

S=rowSums(predYm)
predT = (predYm%*%m$Tr)/matrix(rep(S,m$nt),ncol=m$nt)
RCP = kmeans(predYm, 5)
RCP$cluster = as.factor(RCP$cluster)
mapData=data.frame(m$rL[[1]]$s,S,predT,RCP$cluster)
head(mapData)

obsdelta=NULL
for(i in 1:471){
  obsdelta[i]=sum(m$Y[i+471,])-sum(m$Y[i,])
}

delta=NULL
for(i in 1:471){
  delta[i]=log(S[i+471])-log(S[i])
  #delta[i]=(S[i+471]-S[i])/S[i]
  
}

hist(delta)
mean(delta)

mapData$delta=delta
mean(S[1:471])
mean(S[472:942])
mean(S[472:942])-mean(S[1:471]) #Change in species richness between inventories

# PLOT PREDICTED SPECIES RICHNESS
sp <- ggplot(data = mapData, aes(x=X_manif, y=Y_manif, color=S))+geom_point(size=3)
sp + ggtitle("Predicted species richness") + scale_color_gradient(low="blue", high="red")

# PLOT PREDICTED DELTA SPECIES RICHNESS
#ENVIRONMENTAL COVARIATES
X = as.matrix(read.csv("Z:/data/archipelago/hmsc/data/X.csv"))
X=X[1:471,-c(1,3)]
X=data.frame(rbind(X,X))
head(X)

dim(X)
dim(mapData)

library(ggsn)
mapData$long=mapData$X_manif
mapData$lat=mapData$Y_manif

pdf("plots/delta_map.pdf",width=6,height=5,family="Times")
sp <- ggplot(data = mapData, aes(x=X_manif, y=Y_manif, color=delta*100))+geom_point(alpha=.7,size=3*sqrt(X$area)/200) + scalebar(mapData,location="bottomleft", dist=10,dd2km=NULL,model="WGS84",st.size=4) +north(mapData, symbol=1)
sp + ggtitle("Percent change in species richness between inventories") + scale_colour_gradient2(low="blue",mid="lightgrey",high="red", name="Delta SR")+theme_bw() + xlab("") + ylab("")
dev.off()

#PLOT PREDICTED REGIONS OF COMMON PROFILE
sp <- ggplot(data = mapData, aes(x=X_manif, y=Y_manif, color=RCP.cluster))+geom_point(size=3)
sp + ggtitle("Regions of common profile") 

#### - Cross-validation - ####
a=Sys.time()
partition=createPartition(m, nfolds=4, column=1)
predY_CV4 = computePredictedValues(m, partition=partition, nParallel = 2, updater=list(GammaEta=FALSE))
MF_CV4=evaluateModelFit(m, predY_CV4)
Sys.time()-a

#save(MF_CV2,file="MF_CV2.Rdata") #2fold
save(MF_CV4,file="MF_CV4.Rdata") #4fold

predYm_CV = apply(simplify2array(predY_CV2), 1:2, mean) #Posterior mean
#save(predYm_CV, file="predYm_CV.RData")

#save(predY_CV2,file="crossvalY.Rdata") #2fold
#save(predY_CV4,file="crossval4Y.Rdata") #4fold

load("crossvalY.Rdata") #2fold

load(file="MF_CV2.Rdata") #2fold

R2 = MF_CV2$TjurR2
AUC = MF_CV2$AUC

mean(R2,na.rm=T)
mean(AUC,na.rm=T)


load(file="predYm_CV.RData")

plot(rowSums(m$Y),rowSums(predYm_CV)) #Species richness per island per sampling time
round(cor(rowSums(m$Y),rowSums(predYm_CV))^2,2)
#full: 1.00
#2fold = 0.53

plot(colSums(m$Y),colSums(predYm_CV)) #Occurrences per species
round(cor(colSums(m$Y),colSums(predYm_CV))^2,2)
#full: 1.0
#2fold = 1.0

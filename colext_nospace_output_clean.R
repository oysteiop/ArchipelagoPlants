##################################################################################################
# ARCHIPELAGO PLANTS: OUTPUT OF COL-EXT NO SPACE MODEL
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
setwd("D:/HelsinkiData23102019/archipelago/hmsc/Rcode/colext_nospace")
setwd("C:/data/archipelago/hmsc/Rcode/colext_nospace")

load("model.RData")

#### - Evaluate mixing - ####
post = convertToCodaObject(m)

esBeta = effectiveSize(post$Beta)
esGamma = effectiveSize(post$Gamma)
esEta = effectiveSize(post$Eta[[1]])
esLambda = effectiveSize(post$Lambda[[1]])
str(post$Omega[[1]])
ss = sample(1:3080025, 1000, replace=F)
esOmega = effectiveSize(post$Omega[[1]][,ss])

psBeta = gelman.diag(post$Beta, multivariate=F)$psrf
psGamma = gelman.diag(post$Gamma, multivariate=F)$psrf
psEta = gelman.diag(post$Eta[[1]], multivariate=F)$psrf
psLambda = gelman.diag(post$Lambda[[1]], multivariate=F)$psrf
psOmega = gelman.diag(post$Omega[[1]][,ss], multivariate=F)$psrf

mixing=list(esBeta, esGamma, esEta, esLambda, esOmega,
            psBeta, psGamma, psEta, psLambda, psOmega)
#save(mixing, file="mixing.RData")

load("mixing.RData")

summary(mixing[[1]]) #Beta
summary(mixing[[2]]) #Gamma
summary(mixing[[5]]) #Omega

summary(mixing[[6]]) #Beta
summary(mixing[[7]]) #Gamma
summary(mixing[[10]]) #Omega

#### - Produce posterior trace plots - ####

pdf("posteriorplots/betapost.pdf") #Regression coefficients
plot(post$Beta[,1:200])
dev.off()

pdf("posteriorplots/gammapost.pdf") #Trait effects
plot(post$Gamma)
dev.off()

pdf("posteriorplots/omegapost.pdf") #Species associations
plot(post$Omega[[1]][,1:200])
dev.off()

#### - Compute and plot variance partitioning - ####
group = c(1,1,2,3)
groupnames = c(m$covNames[-1])
VP = computeVariancePartitioning(m, group=group, groupnames=groupnames)
#save(VP, file="VP.RData")

load("VP.RData")
str(VP) #Trait r^2 = 39.6%

# Variance partitioning plot
pdf("plots/varpartColExt.pdf", width=40, height=15)
par(mfrow=c(3,1), mar=c(5,5,4,8))

VPold = VP
VPold$vals = VP$vals[,seq(1,1753,3)]
ng = dim(VPold$vals)[1]
leg = VPold$groupnames
for (r in 1:m$nr) {
  leg = c(leg, paste("random: ", m$levelNames[r], sep = ""))
}
means = round(100 * rowMeans(VPold$vals), 1)
for (i in 1:ng) {
  leg[i] = paste(leg[i], " (mean = ", toString(means[i]), 
                 ")", sep = "")
}
barplot(VPold$vals, main = "Historial inventory", xlab = "Species", 
        ylab = "Variance proportion", legend = leg, col = heat.colors(ng, alpha = 1))

VPcol = VP
VPcol$vals = VP$vals[,seq(2,1754,3)]
ng = dim(VPcol$vals)[1]
leg = VPcol$groupnames
for (r in 1:m$nr) {
  leg = c(leg, paste("random: ", m$levelNames[r], sep = ""))
}
means = round(100 * rowMeans(VPcol$vals), 1)
for (i in 1:ng) {
  leg[i] = paste(leg[i], " (mean = ", toString(means[i]), 
                 ")", sep = "")
}
barplot(VPcol$vals, main = "Colonisation probability", xlab = "Species", 
        ylab = "Variance proportion", legend = leg, col = heat.colors(ng, alpha = 1))

VPext = VP
VPext$vals = VP$vals[,seq(3,1755,3)]
ng = dim(VPext$vals)[1]
leg = VPext$groupnames
for (r in 1:m$nr) {
  leg = c(leg, paste("random: ", m$levelNames[r], sep = ""))
}
means = round(100 * rowMeans(VPext$vals), 1)
for (i in 1:ng) {
  leg[i] = paste(leg[i], " (mean = ", toString(means[i]), 
                 ")", sep = "")
}
barplot(VPext$vals, main = "Extinction probability", xlab = "Species", 
        ylab = "Variance proportion", legend = leg, col = heat.colors(ng, alpha = 1))

dev.off()

#### - Evaluate explanatory power - ####
predY = computePredictedValues(m)
MF = evaluateModelFit(m, predY)
#save(MF,file="MF.RData")

load("MF.RData")

R2 = MF$TjurR2
AUC = MF$AUC

signif(mean(R2[seq(1,1753,3)], na.rm=T), 2) #.31
signif(mean(R2[seq(2,1754,3)], na.rm=T), 2) #.20
signif(mean(R2[seq(3,1755,3)], na.rm=T), 2) #.24
signif(range(R2[seq(1,1753,3)], na.rm=T), 2) #0 - 0.73
signif(range(R2[seq(2,1754,3)], na.rm=T), 2) #0 - 0.63
signif(range(R2[seq(3,1755,3)], na.rm=T), 2) #0 - 0.64

signif(mean(AUC[seq(1,1753,3)], na.rm=T), 2) #0.94
signif(mean(AUC[seq(2,1754,3)], na.rm=T), 2) #0.92
signif(mean(AUC[seq(3,1755,3)], na.rm=T), 2) #0.92
signif(range(AUC[seq(1,1753,3)], na.rm=T), 2) #0.78 - 1.00
signif(range(AUC[seq(2,1754,3)], na.rm=T), 2) #0.68 - 1.00
signif(range(AUC[seq(3,1755,3)], na.rm=T), 2) #0.50 - 1.00

pdf("plots/tjur_vs_prev.pdf",height=3,width=10,family="Times")
xx = colSums(m$Y[,seq(1,1753,3)])/(m$ny)
par(mfrow=c(1,3), mar=c(4,5,2,1))
plot(xx,R2[seq(1,1753,3)], las=1, main=paste("Historical occurrence: Mean = ", round(mean(R2[seq(1,1753,3)], na.rm=T), 2), sep=""), 
     ylim = c(0,1), xlab = "Historical prevalence", ylab=expression(paste("Tjur's ", r^2)))
plot(xx,R2[seq(2,1754,3)], las=1, main=paste("Colonisation probability: Mean = ", round(mean(R2[seq(2,1754,3)], na.rm=T), 2), sep=""), 
     ylim = c(0,1), xlab = "Historical prevalence", ylab=expression(paste("Tjur's ", r^2)))
plot(xx,R2[seq(3,1755,3)], las=1, main=paste("Extinction probability: Mean = ", round(mean(R2[seq(3,1755,3)], na.rm=T), 2), sep=""), 
     ylim = c(0,1), xlab = "Historical prevalence", ylab=expression(paste("Tjur's ", r^2)))
dev.off()

# Explanatory power at island level
predYm = apply(simplify2array(predY), 1:2, mean)
#save(predYm, file = "predYm.RData")

load(file = "predYm.RData")
str(predYm)

tmp = (m$Y>-Inf)*1
pred2 = predYm*tmp
str(pred2)

S = rowSums(pred2[,seq(1,1753,3)], na.rm=T)
Col = rowMeans(pred2[,seq(2,1754,3)], na.rm=T)
Ext = rowMeans(pred2[,seq(3,1755,3)], na.rm=T)
nCol = rowSums(pred2[,seq(2,1754,3)], na.rm=T)
nExt = rowSums(pred2[,seq(3,1755,3)], na.rm=T)

obsS = rowSums(m$Y[,seq(1,1753,3)], na.rm=T)
plot(obsS, S)
lines(0:400, 0:400)
cor(obsS, S)^2 #0.99

obsExt = rowMeans(m$Y[,seq(3,1755,3)], na.rm=T)
plot(obsExt, Ext)
lines(0:1, 0:1)
cor(obsExt, Ext)^2 #0.82

obsnExt = rowSums(m$Y[,seq(3,1755,3)], na.rm=T)
plot(obsnExt, nExt)
lines(0:200, 0:200)
cor(obsnExt, nExt)^2 #0.92

obsCol = rowMeans(m$Y[,seq(2,1754,3)], na.rm=T)
plot(obsCol, Col)
lines(0:1,0:1)
cor(obsCol, Col)^2 #0.90

obsnCol = rowSums(m$Y[,seq(2,1754,3)], na.rm=T)
plot(obsnCol, nCol)
lines(0:100, 0:100)
cor(obsnCol, nCol)^2 #0.90

#### - Trait effects - ####
post = convertToCodaObject(m)
mgamma = getPostEstimate(m, "Gamma")
mgamma
summary(post$Gamma)

plotGamma(m, post=mgamma, param="Mean", cex=c(.7,.7,.7))
plotGamma(m, post=mgamma, param="Support", supportLevel=.95, covOrder="Vector",
          covVector=2:4, trOrder="Vector", trVector=c(1,4,2,5,3,6), trNamesNumbers=c(T,T), cex=c(.6,.8,.8))

#### - Plot association networks - ####
#OmegaCor=computeAssociations(m)
#save(OmegaCor,file="OmegaCor_Full.RData")

load(file="OmegaCor_Full.RData")
str(OmegaCor)

oldCor = OmegaCor
oldCor[[1]]$mean = OmegaCor[[1]]$mean[seq(1,1753,3), seq(1,1753,3)]
oldCor[[1]]$support = OmegaCor[[1]]$support[seq(1,1753,3), seq(1,1753,3)]
str(oldCor)

colCor = OmegaCor
colCor[[1]]$mean = OmegaCor[[1]]$mean[seq(2,1754,3), seq(2,1754,3)]
colCor[[1]]$support = OmegaCor[[1]]$support[seq(2,1754,3), seq(2,1754,3)]
str(colCor)

extCor = OmegaCor
extCor[[1]]$mean = OmegaCor[[1]]$mean[seq(3,1755,3), seq(3,1755,3)]
extCor[[1]]$support = OmegaCor[[1]]$support[seq(3,1755,3), seq(3,1755,3)]
str(extCor)

signif(cor(c(oldCor[[1]]$mean), c(colCor[[1]]$mean))^2, 4)*100 #23.12%
signif(cor(c(oldCor[[1]]$mean), c(extCor[[1]]$mean))^2, 3)*100 #6.51%

signif(cor(c(colCor[[1]]$mean), c(extCor[[1]]$mean))^2, 3)*100 #9.39%

col = colorRampPalette(c("blue3","white","red3"))(200)

pdf("plots/colextCor.pdf",width=9,height=3)
par(mfrow=c(1,3))
supportLevel = 0.75
for (r in 1:m$nr){
  plotOrder = corrMatOrder(oldCor[[r]]$mean, order="AOE")
  toPlot = ((oldCor[[r]]$support>supportLevel) + (oldCor[[r]]$support<(1-supportLevel))>0)*oldCor[[r]]$mean
  corrplot(toPlot[plotOrder,plotOrder], tl.pos="n", method = "color", type="lower", col=colorRampPalette(c("navyblue","white","firebrick"))(200), title=expression(paste("Historical occurrence")), mar=c(0,0,1,0))
}

for (r in 1:m$nr){
  #plotOrder = corrMatOrder(colCor[[r]]$mean,order="AOE")
  toPlot = ((colCor[[r]]$support>supportLevel) + (colCor[[r]]$support<(1-supportLevel))>0)*colCor[[r]]$mean
  corrplot(toPlot[plotOrder,plotOrder], tl.pos="n", method = "color", type="lower", col=colorRampPalette(c("navyblue","white","firebrick"))(200), title=expression(paste("Colonisation probability")), mar=c(0,0,1,0))
}

for (r in 1:m$nr){
  #plotOrder = corrMatOrder(extCor[[r]]$mean,order="AOE")
  toPlot = ((extCor[[r]]$support>supportLevel) + (extCor[[r]]$support<(1-supportLevel))>0)*extCor[[r]]$mean
  corrplot(toPlot[plotOrder,plotOrder], tl.pos="n",method = "color", type="lower", col=colorRampPalette(c("navyblue","white","firebrick"))(200), title=expression(paste("Extinction probability")), mar=c(0,0,1,0))
}
dev.off()

#### - Correlation summary - ####
cmat = OmegaCor[[1]]$mean
cmat[1:6, 1:6]

# (Set weakly supported coefs to NA)
for(i in 1:nrow(cmat)){
  for(j in 1:ncol(cmat)){
    if(OmegaCor[[1]]$support[i,j]>0.25 & OmegaCor[[1]]$support[i,j]<0.75){
      cmat[i,j]=NA
    }
  }
}

# Set diagonal to NA
diag(cmat) = NA

# Historical associations
old = cmat[seq(1,1753,3), seq(1,1753,3)]
n = sum(old>-Inf, na.rm=T)/2 # Number of coefs

oldpold = sum(old>0, na.rm=T)/2 # Positive historical associations
(sum(old>0, na.rm=T)/2)/n*100 # Proportion
oldnold = sum(old<0, na.rm=T)/2 # Negative historical associations
(sum(old<0, na.rm=T)/2)/n*100 # Proportion

mean(old[old>0], na.rm=T) # Mean strength of positive coefs
mean(old[old<0], na.rm=T) # Mean strength of negative coefs

oldpold = c(old[old>0])
oldnold = c(old[old<0])

# Historical-colonisation associations
oldcol = cmat[seq(1,1753,3), seq(2,1754,3)]
oldcol[1:5,1:5]
oldpcol = na.omit(c(oldcol[oldcol>0])) # Positive
oldncol = na.omit(c(oldcol[oldcol<0])) # Negative

# Historical-extinction associations
oldext = cmat[seq(1,1753,3), seq(3,1755,3)]
oldext[1:5,1:5]
oldpext = na.omit(c(oldext[oldext>0])) # Positive
oldnext = na.omit(c(oldext[oldext<0])) # Negative

# Colonisation-colonisation associations
col = cmat[seq(2,1754,3), seq(2,1754,3)]
rcc = rowMeans(col)
n = sum(col>-Inf, na.rm=T)/2 # Number of coefs
col[1:10,1:5]

colpcol = c(col[col>0]) #Positive colonisation-colonisation associations
colncol = c(col[col<0]) #Negative colonisation-colonisation associations

# Extinction-extinction associations
ext = cmat[seq(3,1755,3), seq(3,1755,3)]
ree = rowMeans(ext)
n = sum(ext>-10, na.rm=T)/2 #Number of coefs
ext[1:5,1:5]

extpext = c(ext[ext>0]) # Positive extinction-extinction associations
extnext = c(ext[ext<0]) # Negative extinction-extinction associations

# Colonisation-extinction associations
colext = cmat[seq(2,1754,3), seq(3,1755,3)]
rce = rowMeans(colext)
rec = colMeans(colext)

diag(colext)=NA #Set intraspecific associations to NA
n = sum(colext>-Inf, na.rm=T)/2 # Number of coefs
colext[1:5,1:5]

colpext = na.omit(c(colext[colext>0])) # Positive colonisation-extinction associations
length(colpcol)
colnext = na.omit(c(colext[colext<0])) # Negative colonisation-extinction associations
length(colncol)

# Extinction-colonisation associations
extcol = cmat[seq(3,1755,3), seq(2,1754,3)]
#mean(diag(extcol),na.rm=T)
diag(extcol) = NA
n = (dim(extcol)[1]*dim(extcol)[1]-1)/2
n = sum(extcol>-Inf, na.rm=T)/2

extcol[1:5,1:5]
sum(extcol>0, na.rm=T)/2         # Positive
(sum(extcol>0, na.rm=T)/2)/n*100 # Proportion
sum(extcol<0, na.rm=T)/2         # Negative
(sum(extcol<0, na.rm=T)/2)/n*100 # Proportion

mean(extcol[extcol>0], na.rm=T) # Mean strength of positive extinction-colonisation associations
mean(extcol[extcol<0], na.rm=T) # Mean strength of negative extinction-colonisation associations

#plot(old,col)
#plot(old,ext)
#plot(col,ext)

cor(c(old), c(col), use="pairwise")^2 # 22.3%, well-supported 46.1%
cor(c(old), c(ext), use="pairwise")^2 # 5.6%, well-supported 19.0%
cor(c(col), c(ext), use="pairwise")^2 # 7.7%, well-supported 29.9%

#### -  Boxplot figure - ####
pdf("plots/Figure4.pdf", height=3.5, width=9, family="Times")
par(mfrow=c(1,3), mar=c(8,5,2,1))

# Historical
oldold = c(oldpold,-oldnold)
dat = data.frame(corr=oldold, group=c(rep("old+old",length(oldpold)), rep("old-old",length(oldnold))))
dat$group = relevel(dat$group, ref="old+old")
boxplot(dat$corr~dat$group, varwidth=T, col=c("firebrick","navyblue"), at=c(1,2),
        xlim=c(0,17), xaxt="n", xlab="", ylab="", las=1, cex=0, main=expression(paste("(a) Residual correlations")))
mtext("| r |", 2, line=2.5, cex=1)
# Colonisation
par(new=T)
colcol = c(colpcol,-colncol)
dat = data.frame(corr=colcol, group=c(rep("col+col",length(colpcol)), rep("col-col",length(colncol))))
dat$group = relevel(dat$group, ref="col+col")
boxplot(dat$corr~dat$group, varwidth=T, col=c("firebrick","navyblue"), at=c(4,5),
        xlab="", ylab="", xlim=c(0,17), xaxt="n", yaxt="n", las=1, cex=0)

# Extinction
par(new=T)
extext = c(extpext,-extnext)
dat = data.frame(corr=extext, group=c(rep("ext+ext",length(extpext)), rep("ext-ext",length(extnext))))
dat$group = relevel(dat$group, ref="ext+ext")
boxplot(dat$corr~dat$group, varwidth=T, col=c("firebrick","navyblue"), at=c(7,8),
        xlab="", ylab="", xlim=c(0,17), xaxt="n", yaxt="n", cex=0)

# Historical-colonisation
par(new=T)
oldcol = c(oldpcol,-oldncol)
dat = data.frame(corr=oldcol, group=c(rep("old+col",length(oldpcol)), rep("old-col",length(oldncol))))
dat$group = relevel(dat$group, ref="old+col")
boxplot(dat$corr~dat$group, varwidth=T, col=c("firebrick","navyblue"), at=c(10,11),
        xlab="", ylab="", xlim=c(0,17), xaxt="n", yaxt="n", cex=0)

# Historial-extinction
par(new=T)
oldext = c(oldpext,-oldnext)
dat = data.frame(corr=oldext, group=c(rep("old+ext",length(oldpext)), rep("old-ext",length(oldnext))))
dat$group = relevel(dat$group, ref="old+ext")
boxplot(dat$corr~dat$group, varwidth=T, col=c("firebrick","navyblue"), at=c(13,14),
        xlab="", ylab="", xlim=c(0,17), xaxt="n", yaxt="n", cex=0)

# Colonisation-extinction
par(new=T)
colext = c(colpext,-colnext)
dat = data.frame(corr=colext, group=c(rep("col+ext",length(colpext)), rep("col-ext",length(colnext))))
dat$group = relevel(dat$group, ref="col+ext")
boxplot(dat$corr~dat$group, varwidth=T, col=c("firebrick","navyblue"), at=c(16,17),
        xlab="", ylab="", xlim=c(0,17), xaxt="n", yaxt="n", cex=0)

axis(1, at=c(1.5,4.5,7.5,10.5,13.5,16.5), labels=rep("", 6))
labs = c("Hist - Hist","Col - Col","Ext - Ext","Hist - Col","Hist - Ext","Col - Ext")
text(c(1.5,4.5,7.5,10.5,13.5,16.5), par("usr")[3]-.05, srt=45, adj=1, labels=labs, xpd=T)
legend("topright", pch=15, cex=1.25, col=c("firebrick","navyblue"), c("Positive", "Negative"), bty="n")
#dev.off()

#### - Correlations Shoreline species only - ####

oldCor = OmegaCor
oldCor[[1]]$mean = OmegaCor[[1]]$mean[seq(1,1753,3), seq(1,1753,3)]
oldCor[[1]]$support = OmegaCor[[1]]$support[seq(1,1753,3), seq(1,1753,3)]
oldCor[[1]]$mean = oldCor[[1]]$mean[which(m$TrData[seq(1,1753,3),1]=="yes"), which(m$TrData[seq(1,1753,3),1]=="yes")]
oldCor[[1]]$support = oldCor[[1]]$support[which(m$TrData[seq(1,1753,3),1]=="yes"), which(m$TrData[seq(1,1753,3),1]=="yes")]
str(oldCor)

colCor = OmegaCor
colCor[[1]]$mean = OmegaCor[[1]]$mean[seq(2,1754,3), seq(2,1754,3)]
colCor[[1]]$support = OmegaCor[[1]]$support[seq(2,1754,3), seq(2,1754,3)]
colCor[[1]]$mean = colCor[[1]]$mean[which(m$TrData[seq(2,1754,3),1]=="yes"), which(m$TrData[seq(2,1754,3),1]=="yes")]
colCor[[1]]$support = colCor[[1]]$support[which(m$TrData[seq(2,1754,3),1]=="yes"), which(m$TrData[seq(2,1754,3),1]=="yes")]
str(colCor)

extCor = OmegaCor
extCor[[1]]$mean = OmegaCor[[1]]$mean[seq(3,1755,3), seq(3,1755,3)]
extCor[[1]]$support = OmegaCor[[1]]$support[seq(3,1755,3), seq(3,1755,3)]
extCor[[1]]$mean = extCor[[1]]$mean[which(m$TrData[seq(3,1755,3),1]=="yes"), which(m$TrData[seq(3,1755,3),1]=="yes")]
extCor[[1]]$support = extCor[[1]]$support[which(m$TrData[seq(3,1755,3),1]=="yes"), which(m$TrData[seq(3,1755,3),1]=="yes")]
str(extCor)

pdf("plots/colextCorShoreline.pdf",width=9,height=3)
par(mfrow=c(1,3))
supportLevel = 0.75
for (r in 1:m$nr){
  plotOrder = corrMatOrder(oldCor[[r]]$mean, order="AOE")
  toPlot = ((oldCor[[r]]$support>supportLevel) + (oldCor[[r]]$support<(1-supportLevel))>0)*oldCor[[r]]$mean
  corrplot(toPlot[plotOrder,plotOrder], tl.pos="n", method = "color", type="lower", col=colorRampPalette(c("navyblue","white","firebrick"))(200), title=expression(paste("Historical occurrence")), mar=c(0,0,1,0))
}

for (r in 1:m$nr){
  #plotOrder = corrMatOrder(colCor[[r]]$mean,order="AOE")
  toPlot = ((colCor[[r]]$support>supportLevel) + (colCor[[r]]$support<(1-supportLevel))>0)*colCor[[r]]$mean
  corrplot(toPlot[plotOrder,plotOrder], tl.pos="n", method = "color", type="lower", col=colorRampPalette(c("navyblue","white","firebrick"))(200), title=expression(paste("Colonisation probability")), mar=c(0,0,1,0))
}

for (r in 1:m$nr){
  #plotOrder = corrMatOrder(extCor[[r]]$mean,order="AOE")
  toPlot = ((extCor[[r]]$support>supportLevel) + (extCor[[r]]$support<(1-supportLevel))>0)*extCor[[r]]$mean
  corrplot(toPlot[plotOrder,plotOrder], tl.pos="n", method = "color", type="lower", col=colorRampPalette(c("navyblue","white","firebrick"))(200), title=expression(paste("Extinction probability")), mar=c(0,0,1,0))
}
dev.off()


cor(c(oldCor[[1]]$mean), c(colCor[[1]]$mean))^2 #34.5%
cor(c(oldCor[[1]]$mean), c(extCor[[1]]$mean))^2 #22.2%
cor(c(colCor[[1]]$mean), c(extCor[[1]]$mean))^2 #20.4%

#### - Cross-validation - ####
a = Sys.time()
partition = createPartition(m, nfolds=2, column=1)
predY_CV2 = computePredictedValues(m, partition=partition, updater=list(GammaEta=FALSE))
MF_CV2 = evaluateModelFit(m, predY_CV2)
Sys.time()-a

#save(MF_CV2,file="MF_CV2.Rdata") #2fold

#save(predY_CV2,file="crossvalY.Rdata") #2fold

load("MF_CV2.Rdata") #2fold
load("crossvalY.Rdata") #5fold

R2 = MF_CV2$TjurR2
AUC = MF_CV2$AUC

signif(mean(R2[seq(1,1753,3)], na.rm=T), 2) #.033
signif(mean(R2[seq(2,1754,3)], na.rm=T), 2) #.029
signif(mean(R2[seq(3,1755,3)], na.rm=T), 2) # 0
signif(range(R2[seq(1,1753,3)], na.rm=T), 2) #0 - 0.26
signif(range(R2[seq(2,1754,3)], na.rm=T), 2) #0 - 0.19
signif(range(R2[seq(3,1755,3)], na.rm=T), 2) #0 - 0.17

signif(mean(AUC[seq(1,1753,3)], na.rm=T), 2) #0.74
signif(mean(AUC[seq(2,1754,3)], na.rm=T), 2) #0.72
signif(mean(AUC[seq(3,1755,3)], na.rm=T), 2) #0.66
signif(range(AUC[seq(1,1753,3)], na.rm=T), 2) #0.44 - 0.99
signif(range(AUC[seq(2,1754,3)], na.rm=T), 2) #0.41 - 0.98
signif(range(AUC[seq(3,1755,3)], na.rm=T), 2) #0.38 - 1.00

predYm_CV = apply(simplify2array(predY_CV2), 1:2, mean) #Posterior mean
#save(predYm_CV, file="predYm_CV.RData")

load("predYm_CV.RData")

tmp = (m$Y>-Inf)*1
pred2 = predYm_CV*tmp
str(pred2)

S = rowSums(pred2[,seq(1,1753,3)], na.rm=T)
Col = rowMeans(pred2[,seq(2,1754,3)], na.rm=T)
Ext = rowMeans(pred2[,seq(3,1755,3)], na.rm=T)
nCol = rowSums(pred2[,seq(2,1754,3)], na.rm=T)
nExt = rowSums(pred2[,seq(3,1755,3)], na.rm=T)

obsS = rowSums(m$Y[,seq(1,1753,3)], na.rm=T)
plot(obsS, S)
lines(0:400, 0:400)
cor(obsS, S)^2
#2fold = 0.45

obsnExt = rowSums(m$Y[,seq(3,1755,3)], na.rm=T)
plot(obsnExt, nExt)
lines(0:200, 0:200)
cor(obsnExt, nExt)^2 #0.92
#2fold = 0.65

obsnCol = rowSums(m$Y[,seq(2,1754,3)], na.rm=T)
plot(obsnCol, nCol)
lines(0:100, 0:100)
cor(obsnCol, nCol)^2 #0.89
#2fold = 0.39

#### - Gradient plot - ####

# type = 1 fixes to the most likely value (defined as expected value for covarites, mode for factors)
# type = 2 fixes to most likely value, given the value of focal variable,based on linear relationship
# type = 3 fixes to the value given.

GradientA2 = constructGradient(m, focalVariable = "area", non.focalVariables = list(
  buff5=list(2),
  sd_height=list(2)))

GradientT2 = constructGradient(m, focalVariable = "sd_height", non.focalVariables = list(
  buff5=list(2),
  area=list(2)))

GradientB2 = constructGradient(m, focalVariable = "buff5", non.focalVariables = list(
  sd_height=list(2),
  area=list(2)))

GradientA1 = constructGradient(m, focalVariable = "area", non.focalVariables = list(
  buff5=list(1),
  sd_height=list(1)))

GradientT1 = constructGradient(m, focalVariable = "sd_height", non.focalVariables = list(
  buff5=list(1),
  area=list(1)))

GradientB1 = constructGradient(m, focalVariable = "buff5", non.focalVariables = list(
  sd_height=list(1),
  area=list(1)))

predYA2 = predict(m, Gradient = GradientA2, expected = TRUE)
predYT2 = predict(m, Gradient = GradientT2, expected = TRUE)
predYB2 = predict(m, Gradient = GradientB2, expected = TRUE)

predYA1 = predict(m, Gradient = GradientA1, expected = TRUE)
predYT1 = predict(m, Gradient = GradientT1, expected = TRUE)
predYB1 = predict(m, Gradient = GradientB1, expected = TRUE)

# Plotting function
plotpred = function(predY=NA , Gradient =NA , xlabel=NA, ylabel=NA){
  predShist = lapply(predY, function(x) rowMeans(x[,seq(1,1753,3)], na.rm=T))
  predScol = lapply(predY, function(x) rowMeans(x[,seq(2,1754,3)], na.rm=T))
  predSext = lapply(predY, function(x) rowMeans(x[,seq(3,1755,3)], na.rm=T))

  qpredhist = apply(abind(predShist, along = 2), c(1), quantile, prob = c(.025,0.5,.975))
  qpredcol = apply(abind(predScol, along = 2), c(1), quantile, prob = c(.025,0.5,.975))
  qpredext = apply(abind(predSext, along = 2), c(1), quantile, prob = c(.025,0.5,.975))
  colCond = qpredcol*qpredhist[c(2,2,2),]
  extCond = qpredext*qpredhist[c(2,2,2),]

  xx = Gradient$XDataNew[,1]
  plot(xx, qpredhist[2,], col="white", ylim=c(0,.75), bty="l", las=1, xlab="", ylab="")
  mtext(xlabel, 1, line=2.5, cex=.8)
  mtext(ylabel, 2, line=2.5, cex=.8)
  polygon(c(xx, rev(xx)), c(qpredhist[1, ], rev(qpredhist[3, ])), 
        col = "grey75", border = FALSE)
  lines(xx, qpredhist[2, ], lwd = 2)

  polygon(c(xx, rev(xx)), c(qpredext[1, ], rev(qpredext[3, ])), 
        col = rgb(.698,.133,.133,.5), border = FALSE)
  lines(xx, qpredext[2, ], lwd = 2)

  polygon(c(xx, rev(xx)), c(qpredcol[1, ], rev(qpredcol[3, ])), 
        col = rgb(0,0,0.502,.5), border = FALSE)
  lines(xx, qpredcol[2, ], lwd = 2)
}

pdf("plots/gradient_plot.pdf", width=7, height=4.3, family="Times")
par(mfrow=c(2,3), mar=c(2,4,3,1))
plotpred(predY=predYA2, Gradient=GradientA2, xlabel="", ylabel="Probability")
plotpred(predY=predYT2, Gradient=GradientT2, xlabel="", ylabel="")
plotpred(predY=predYB2, Gradient=GradientB2, xlabel="", ylabel="")
par(xpd=T)
legend(x=12, y=.8, pch=15, col=c(rgb(0,0,0.502,.5), rgb(.698,.133,.133,.5)), 
       c("Colonisation", "Extinction"), cex=1.25, bty="n")
par(mar=c(4,4,1,1), xpd=F)
plotpred(predY=predYA1, Gradient=GradientA1, xlabel="log (Island area)", ylabel="Probability")
plotpred(predY=predYT1, Gradient=GradientT1, xlabel="Topographic complexity", ylabel="")
plotpred(predY=predYB1, Gradient=GradientB1, xlabel="log (Neighbourhood size)", ylabel="")
dev.off()


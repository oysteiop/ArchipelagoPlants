##################################################################################################
# ARCHIPELAGO PLANTS: CONDITIONAL PREDICTIONS
##################################################################################################
rm(list=ls())
gc()
memory.size()
memory.limit(64000)

##################################################################################################
# READ PACKAGES
##################################################################################################
library(devtools)
library(BayesLogit)
#install_github("hmsc-r/HMSC", build_opts = c("--no-resave-data", "--no-manual"))
library(Hmsc)
library(abind)
library(reshape2)
library(plyr)

##################################################################################################
# SET DIRECTORIES AND READ THE FITTED MODEL
##################################################################################################
setwd("F:/HelsinkiData23102019/archipelago/hmsc/Rcode/colext_nospace")
load("model.RData")

Y = as.matrix(read.csv("F:/HelsinkiData23102019/archipelago/hmsc/data/oldnewY.csv"))
Y=Y[472:942,]
Y=Y[,order(colnames(Y))]

#### - Conditional predictions - ####
predList = list()
qpredList = list()
cpredList = list()
epredList = list()
tjurMat = matrix(NA, nrow=585, ncol=12)

post = poolMcmcChains(m$postList, thin=10)
mcmc = 10

a = Sys.time()
for(i in 1:1){
  print(paste("Now computing species", i))

sp = seq(((3*i)-2),(3*i),by=1)

# True occurrences
obs = Y[,i]
qobs = m$Y[,sp[1]]
cobs = m$Y[,sp[2]]
eobs = m$Y[,sp[3]]

# Full model predictions: m0
Yc = m$Y
Yc[,sp] = NA #Set focal species to NA
predC = predict(m, post=post, Yc=Yc, mcmcStep=mcmc, expected=TRUE)
m0m = apply(simplify2array(predC), 1:2, mean)
Yc[,seq(2,1754,3)] = NA #Set all colonisation to NA
Yc[,seq(3,1755,3)] = NA #Set all extinction to NA
predC = predict(m, post=post, Yc=Yc, mcmcStep=mcmc, expected=TRUE)
m0bm = apply(simplify2array(predC), 1:2, mean)

q0 = m0bm[,sp[1]]
c0 = m0m[,sp[2]]
e0 = m0m[,sp[3]]
p0 = q0*(1-e0)+(1-q0)*c0

# Effect of covariates: m1
Yc = m$Y
Yc[,sp] = NA #Set focal species to NA
Xc = m$X
Xc[,2] = mean(m$X[,2]) #Set covariates to mean
Xc[,3] = mean(m$X[,3])
Xc[,4] = mean(m$X[,4])

predC = predict(m, post=post, Yc=Yc, X=Xc, mcmcStep=mcmc, expected=TRUE)
m1m = apply(simplify2array(predC), 1:2, mean)
Yc[,seq(2,1754,3)] = NA #Set all colonisation to NA
Yc[,seq(3,1755,3)] = NA #Set all extinction to NA
predC = predict(m, post=post, Yc=Yc, X=Xc, mcmcStep=mcmc, expected=TRUE)
m1bm = apply(simplify2array(predC), 1:2, mean)

# Effect of covariates on past occurrence
q1 = m1bm[,sp[1]]
c1 = m0m[,sp[2]]
e1 = m0m[,sp[3]]
p1 = q1*(1-e1)+(1-q1)*c1

# Effect of covariates on colonization
q2 = m0bm[,sp[1]]
c2 = m1m[,sp[2]]
e2 = m0m[,sp[3]]
p2 = q2*(1-e2)+(1-q2)*c2

# Effect of covariates on extinction
q3 = m0bm[,sp[1]]
c3 = m0m[,sp[2]]
e3 = m1m[,sp[3]]
p3 = q3*(1-e3)+(1-q3)*c3

# Effect of historical community: m2
Yc = m$Y
Yc[,sp] = NA #Set focal species to NA
Yc[,seq(1,1753,3)] = NA #Set all historical occurrences to NA
predC = predict(m, post=post, Yc=Yc, mcmcStep=mcmc, expected=TRUE)
m2m = apply(simplify2array(predC), 1:2, mean)
Yc[,seq(2,1754,3)] = NA #Set all colonisation to NA
Yc[,seq(3,1755,3)] = NA #Set all extinction to NA
predC = predict(m, post=post, Yc=Yc, mcmcStep=mcmc, expected=TRUE)
m2bm = apply(simplify2array(predC), 1:2, mean)

# Effect of past occurrence on past occurrence
q4 = m2bm[,sp[1]]
c4 = m0m[,sp[2]]
e4 = m0m[,sp[3]]
p4 = q4*(1-e4)+(1-q4)*c4

# Effect of past occurrence on colonization
q5 = m0bm[,sp[1]]
c5 = m2m[,sp[2]]
e5 = m0m[,sp[3]]
p5 = q5*(1-e5)+(1-q5)*c5

# Effect of past occurrence on extinction
q6 = m0bm[,sp[1]]
c6 = m0m[,sp[2]]
e6 = m2m[,sp[3]]
p6 = q6*(1-e6)+(1-q6)*c6

# Effect of species colonisations: m3
Yc = m$Y
Yc[,sp] = NA #Set focal species to NA
Yc[,seq(2,1754,3)] = NA #Set all colonisations to NA
predC = predict(m, post=post, Yc=Yc, mcmcStep=mcmc, expected=TRUE)
m3m = apply(simplify2array(predC), 1:2, mean)

# Effect of colonisation on colonization
q7 = m0bm[,sp[1]]
c7 = m3m[,sp[2]]
e7 = m0m[,sp[3]]
p7 = q7*(1-e7)+(1-q7)*c7

# Effect of colonisation on extinction
q8 = m0bm[,sp[1]]
c8 = m0m[,sp[2]]
e8 = m3m[,sp[3]]
p8 = q8*(1-e8)+(1-q8)*c8

# Effect of species extinctions: m4
Yc = m$Y
Yc[,sp] = NA #Set focal species to NA
Yc[,seq(3,1755,3)] = NA #Set all extinctions to NA
predC = predict(m, post=post, Yc=Yc, mcmcStep=mcmc, expected=TRUE)
m4m = apply(simplify2array(predC), 1:2, mean)

# Effect of extinction on colonization
q9 = m0bm[,sp[1]]
c9 = m4m[,sp[2]]
e9 = m0m[,sp[3]]
p9 = q9*(1-e9)+(1-q9)*c9

# Effect of extinction on extinction
q10 = m0bm[,sp[1]]
c10 = m0m[,sp[2]]
e10 = m4m[,sp[3]]
p10 = q10*(1-e10)+(1-q10)*c10

# Compile
condPreds = data.frame(island=1:471, p0=p0, p1=p1, p2=p2, p3=p3, p4=p4, p5=p5, p6=p6, p7=p7, p8=p8, p9=p9, p10=p10, obs=obs)
qPreds = data.frame(island=1:471, q0=q0, q1=q1, q2=q2, q3=q3, q4=q4, q5=q5, q6=q6, q7=q7, q8=q8, q9=q9, q10=q10, obs=qobs)
cPreds = data.frame(island=1:471, c0=c0, c1=c1, c2=c2, c3=c3, c4=c4, c5=c5, c6=c6, c7=c7, c8=c8, c9=c9, c10=c10, obs=cobs)
ePreds = data.frame(island=1:471, e0=e0, e1=e1, e2=e2, e3=e3, e4=e4, e5=e5, e6=e6, e7=e7, e8=e8, e9=e9, e10=e10, obs=eobs)

tjur_vals = apply(condPreds[,-1], 2, function(x){mean(x[which(obs==1)])-mean(x[which(obs==0)])})

predList[[i]] = condPreds
qpredList[[i]] = qPreds
cpredList[[i]] = cPreds
epredList[[i]] = ePreds

tjurMat[i,] = tjur_vals
}

b = Sys.time()
b-a
#5 species thin 10 mcmc 10 = 14.15 hours, 2.83 hours per species

save(predList, file="cond_preds/predList.RData")

save(tjurMat, file="cond_preds/tjurMat.RData")

load(file="predList.RData")

#Combine from remote computer
list.files()
load(file = "cond_preds/predList1-100.RData")
pred1_100 = predList
load(file = "cond_preds/predList101-200.RData")
pred101_200 = predList
load(file = "cond_preds/predList201-300.RData")
pred201_300 = predList
load(file = "cond_preds/predList301-400.RData")
pred301_400 = predList
load(file = "cond_preds/predList401-500.RData")
pred401_500 = predList
load(file = "cond_preds/predList501-585.RData")
pred501_585 = predList

pred1_100 = lapply(c(1:100), function(x){pred1_100[[x]]})
pred101_200 = lapply(c(101:200), function(x){pred101_200[[x]]})
pred201_300 = lapply(c(201:300), function(x){pred201_300[[x]]})
pred301_400 = lapply(c(301:400), function(x){pred301_400[[x]]})
pred401_500 = lapply(c(401:500), function(x){pred401_500[[x]]})
pred501_585 = lapply(c(501:585), function(x){pred501_585[[x]]})

predList=c(pred1_100, pred101_200, pred201_300, pred301_400, pred401_500, pred501_585)
length(predList)

# RMSD
mse = function(xx){
  se = NULL
  for(i in 1:length(xx)){
    se[i] = (xx[i]-x[i,2])^2
  }
  mean(sqrt(se))
  }

out = list()
for(i in 1:length(predList)){
  x = predList[[i]]
  out[[i]] = apply(x[,2:12], 2, function(xx){mse(xx)})
}
names(out) = paste("sp", 1:length(predList), sep="_")

out2 = lapply(out, as.data.frame)
vals = rbind.fill(out2)
vals = vals[,1]
mat = matrix(vals, nrow=11, ncol=length(predList), byrow=F)
means = apply(mat,1,mean)
ses = apply(mat,1,sd)

#rmsd
plot(1:11,means,ylim=c(-0.01,0.04),pch=c(1,16,17,15,16,17,15,17,15,17,15),
     xaxt="n",ylab="Root mean square deviation",xlab="",main=expression(paste("(b) Influence on model predictions")))
segments(1:11,means-ses,1:11,means+ses)
abline(h=0,lty=2)
axis(1,c(1,3,6,8.5,10.5),labels=F,xlab="")
text(c(1,3,6,8.5,10.5),par("usr")[3]-.003,srt=45,adj=1,labels=c("Full model", "Environment","Historical occurrence", "Colonisation", "Extinction"),xpd=T)

abline(v=1.5,lty=2)
abline(v=4.5,lty=2)
abline(v=7.5,lty=2)
abline(v=9.5,lty=2)
legend("topright",pch=c(1,16,17,15),c("Full model", "Historical occurrence", "Colonisation", "Extinction"),bg="white")

# Tjur####
tjur = function(x){
tval = NULL
for(i in 2:12){
tval[i-1] = mean(x[which(x[,13]==1),i])-mean(x[which(x[,13]==0),i])
}
return(tval)
}

out = lapply(predList,function(xx){tjur(xx)})
out2 = lapply(out,as.data.frame)
vals = rbind.fill(out2)
vals = vals[,1]

mat = matrix(vals,nrow=11,ncol=length(predList),byrow=F)
means = apply(mat,1,mean,na.rm=T)
ses = apply(mat,1,sd,na.rm=T)/sqrt(ncol(mat))

plot(1:11,means,ylim=c(0.2,.35),xlab="",xaxt="n",pch=c(1,16,17,15,16,17,15,17,15,17,15),las=1,
     main=expression(paste("(c) Influence on Tjur ",r^2)),ylab=expression(paste("Tjur ",r^2)))
abline(h=means[1])
segments(1:11,means-ses,1:11,means+ses)
axis(1,c(1,3,6,8.5,10.5),labels=F,xlab="")
text(c(1,3,6,8.5,10.5),par("usr")[3]-.007,srt=45,adj=1,labels=c("Full model", "Environment","Historical occurrence", "Colonisation", "Extinction"),xpd=T)

abline(v=1.5,lty=2)
abline(v=4.5,lty=2)
abline(v=7.5,lty=2)
abline(v=9.5,lty=2)
legend("topright",pch=c(1,16,17,15),c("Full model", "Historical occurrence", "Colonisation", "Extinction"),bg="white")

dev.off()


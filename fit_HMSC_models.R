##################################################################################################
### ARCHIPELAGO PLANTS: HMSC MODEL FITTING
##################################################################################################
rm(list=ls())
gc()

# READ PACKAGES
library(devtools)
#install_url('https://cran.r-project.org/src/contrib/Archive/BayesLogit/BayesLogit_0.6.tar.gz')
library(BayesLogit)
library(mvtnorm)
#install_github("hmsc-r/HMSC", build_opts = c("--no-resave-data", "--no-manual"))
library(Hmsc)
library(ape)

setwd("D:/HelsinkiData23102019/archipelago/hmsc/")

################################################################################
# MODEL 1: BOTH DATASETS WITH TIME AS COVARIATE ('REPEATED MEASURES' MODEL) ####
################################################################################

# Load datafiles
Y = read.csv("Rcode/repeatedmeasures/Y.csv")
Y = as.matrix(Y)
X = read.csv("Rcode/repeatedmeasures/X.csv")
Tr = read.csv("Rcode/repeatedmeasures/Tr.csv")
dfPi = read.csv("Rcode/repeatedmeasures/dfPi.csv")
colnames(dfPi) = c("island","su")
xy = as.matrix(read.csv("Rcode/repeatedmeasures/xy.csv"))

# RANDOM EFFECT STRUCTURE
rL1 = HmscRandomLevel(units=unique(as.character(dfPi[,1]))) #Non-spatial random effect for islands
rL1$nfMax=4
rL2 = HmscRandomLevel(units=unique(as.character(dfPi[,2]))) #Non-spatial random effects for observation (sampling unit)
rL2$nfMax=4

#Set model formulae
XFormula = as.formula(paste("~", paste0(colnames(X), collapse=" + ")))
TrFormula = ~-1+Shoreplant

m = Hmsc(Y = Y, XData = X, XFormula = XFormula, TrData = Tr, TrFormula = TrFormula, dist = "probit", 
         studyDesign = dfPi, ranLevels=list(island=rL1, su=rL2))

# RUN MCMC
samples = 1000
thin = 200
transient = .5*(thin*samples)
adaptNf = 0.4*(thin*samples)
nChains = 2

a1 = Sys.time()
m = sampleMcmc(m, samples = samples, thin = thin, adaptNf=rep(adaptNf, m$nr), 
               transient = transient, nChains = nChains, updater=list(GammaEta=FALSE))
b1 = Sys.time()
#2*300k in 7.4 days

# SAVE RESULTS TO FILE
save(m, file="Rcode/repeatedmeasures/modelNS.RData")

#############################################################################
# MODEL 2A: HISTORICAL OCCURRENCE WITH PHYLO ('TAXONOMIC-SIGNAL MODEL A")####
#############################################################################

# Load datafiles
Y = read.csv("Rcode/histocc/Y.csv")
Y = as.matrix(Y)
X = read.csv("Rcode/histocc/X.csv")
Tr = read.csv("Rcode/histocc/Tr.csv")
dfPi = read.csv("Rcode/histocc/dfPi.csv")
colnames(dfPi) = c("island")
xy = as.matrix(read.csv("Rcode/histocc/xy.csv"))
tree2 = read.tree("Rcode/histocc/tree2.txt")

# RANDOM EFFECT STRUCTURE
rL1 = HmscRandomLevel(units=unique(as.character(dfPi[,1]))) #Spatial random effect for islands
rL1$nfMax=4

XFormula = as.formula(paste("~", paste0(colnames(X)[-1], collapse=" + ")))
TrFormula = ~Shoreplant

m = Hmsc(Y = Y, XData = X, XFormula = XFormula, TrData = Tr, TrFormula = TrFormula, phyloTree = tree2, 
         dist = "probit", studyDesign = dfPi, ranLevels = list(island=rL1))

# RUN MCMC
samples = 1000
thin = 40
transient = 0.5*(thin*samples)
adaptNf = 0.4*(thin*samples)
nChains = 2

a2 = Sys.time()
m = sampleMcmc(m, samples = samples, thin = thin, adaptNf=rep(adaptNf, m$nr), 
               transient = transient, nChains = nChains, updater=list(GammaEta=FALSE))
b2 = Sys.time()

# SAVE RESULTS TO FILE
save(m, file="Rcode/histocc/modelNS.Rdata")

####################################################################
# MODEL 2B: COLONISATION WITH PHYLO ('TAXONOMIC-SIGNAL MODEL B")####
####################################################################

# Load datafiles
Y = read.csv("Rcode/colonisation/Y.csv")
Y = as.matrix(Y)
X = read.csv("Rcode/colonisation/X.csv")
Tr = read.csv("Rcode/colonisation/Tr.csv")
dfPi = read.csv("Rcode/colonisation/dfPi.csv")
colnames(dfPi) = c("island")
xy = as.matrix(read.csv("Rcode/colonisation/xy.csv"))
tree2 = read.tree("Rcode/colonisation/tree2.txt")

# RANDOM EFFECT STRUCTURE
rL1 = HmscRandomLevel(units=unique(as.character(dfPi[,1]))) #Spatial random effect for islands
rL1$nfMax=4

XFormula = as.formula(paste("~", paste0(colnames(X)[-1], collapse=" + ")))
TrFormula = ~Shoreplant

m = Hmsc(Y = Y, XData = X, XFormula = XFormula, TrData = Tr, TrFormula = TrFormula, phyloTree = tree2, 
         dist = "probit", studyDesign = dfPi, ranLevels=list(island=rL1))

# RUN MCMC
samples = 1000
thin = 40
transient = 0.5*(thin*samples)
adaptNf = 0.4*(thin*samples)
nChains = 2

a3 = Sys.time()
m = sampleMcmc(m, samples = samples, thin = thin, adaptNf=rep(adaptNf, m$nr), 
               transient = transient, nChains = nChains, updater=list(GammaEta=FALSE))
b3 = Sys.time()

# SAVE RESULTS TO FILE
save(m, file="Rcode/colonisation/modelNS.Rdata")

##################################################################################################
# MODEL 2C: EXTINCTION WITH PHYLO ('TAXONOMIC-SIGNAL MODEL C")####
##################################################################################################

# Load datafiles
Y = read.csv("Rcode/extinction/Y.csv")
Y = as.matrix(Y)
X = read.csv("Rcode/extinction/X.csv")
Tr = read.csv("Rcode/extinction/Tr.csv")
dfPi = read.csv("Rcode/extinction/dfPi.csv")
colnames(dfPi) = c("island")
xy = as.matrix(read.csv("Rcode/extinction/xy.csv"))
tree2 = read.tree("Rcode/extinction/tree2.txt")

# RANDOM EFFECT STRUCTURE
rL1 = HmscRandomLevel(units=unique(as.character(dfPi[,1]))) #Non-spatial random effect for islands
rL1$nfMax=4

XFormula = as.formula(paste("~", paste0(colnames(X)[-1], collapse=" + ")))
TrFormula = ~Shoreplant

m = Hmsc(Y = Y, XData = X, XFormula = XFormula, TrData = Tr, TrFormula = TrFormula, phyloTree = tree2, 
         dist = "probit", studyDesign = dfPi, ranLevels = list(island=rL1))

# RUN MCMC
samples = 1000
thin = 40
transient = 0.5*(thin*samples)
adaptNf = 0.4*(thin*samples)
nChains = 2

a4 = Sys.time()
m = sampleMcmc(m, samples = samples, thin = thin, adaptNf=rep(adaptNf, m$nr), 
               transient = transient, nChains = nChains, updater=list(GammaEta=FALSE))
b4 = Sys.time()

# SAVE RESULTS TO FILE
save(m, file="Rcode/extinction/modelNS.Rdata")

##################################################################################################
# MOD 3: COL/EXT MODEL####
##################################################################################################

# Load datafiles
Y = read.csv("Rcode/colext_nospace/Y.csv")
Y = as.matrix(Y)
X = read.csv("Rcode/colext_nospace/X.csv")
Tr2 = read.csv("Rcode/colext_nospace/Tr2.csv")
dfPi = read.csv("Rcode/colext_nospace/dfPi.csv")
colnames(dfPi) = c("island")
xy = as.matrix(read.csv("Rcode/colext_nospace/xy.csv"))

# RANDOM EFFECT STRUCTURE
sRL = xy
rownames(sRL) = as.character(dfPi[,1])
rL = HmscRandomLevel(units=unique(dfPi[,1]))

# CONSTRUCT THE MODEL
XFormula = as.formula(paste("~", paste0(colnames(X), collapse=" + ")))
TrFormula = ~-1+column_type:Shoreplant

m = Hmsc(Y = Y, XData = X, XFormula = XFormula, TrData = Tr2, TrFormula = TrFormula, 
             dist = "probit", studyDesign = dfPi, ranLevels = list(island=rL))

# RUN MCMC
samples = 1000
thin = 200
transient = .5*(thin*samples)
adaptNf = 0.4*(thin*samples)
nChains = 2

a6 = Sys.time()
m = sampleMcmc(m, samples = samples, thin = thin, adaptNf=rep(adaptNf, m$nr), 
               transient = transient, nChains = nChains, updater=list(GammaEta=FALSE))
b6 = Sys.time()

# SAVE RESULTS TO FILE
save(m, file="Rcode/colext_nospace/model.Rdata")

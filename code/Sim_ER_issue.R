#=============================================================
#Run samSim simulations 
#using the harrison chinook data as example
#samSim updates
#Catarina Wor
#October 2022 
#=============================================================

#install samsim 
#remotes::install_github("Pacific-salmon-assess/samSim", ref="timevar", force=TRUE)


library(samSim)
library(ggplot2)
library(devtools)
library(gridExtra)
library(dplyr)
library(here)
## Load relevant input data
# Simulation run parameters describing different scenarios
simPar <- read.csv("data/harckER/harcnkSimPars_ER.csv")
# CU-specific parameters
cuPar <- read.csv("data/harckER/harcnkCUPars.csv")

## Store relevant object names to help run simulation 
scenNames <- unique(simPar$scenario)
dirNames <- sapply(scenNames, function(x) paste(x, unique(simPar$species),sep = "_"))

## First check to ensure that a single scenario can be run (only a small number
# of trials necessary)

#Run and save simulated data

a<-4
  
   genericRecoverySim(simPar=simPar[a,], 
    cuPar=cuPar, 
    catchDat=NULL, 
    srDat=NULL,
    variableCU=FALSE, 
    ricPars=NULL, 
    larkPars=NULL,
    cuCustomCorrMat= NULL,
    outDir="outs", 
    nTrials=100, 
    makeSubDirs=TRUE, 
    random=FALSE, 
    uniqueProd=TRUE,
    uniqueSurv=FALSE)


dat <- readRDS(paste0("outs/SamSimOutputs/simData/", simPar$nameOM[a],"/",simPar$scenario[a],"/",
                         paste(simPar$nameOM[a],"_", simPar$nameMP[a], "_", "CUsrDat.RData",sep="")))$srDatout
  

dat[dat$ER>.999,]
 
  
  

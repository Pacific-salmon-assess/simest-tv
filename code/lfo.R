#=============================================================
#Example use of samsim for time varying simulation evaluation
#using the coho data as it is compliant with the most recent
#samSim updates
#Catarina Wor
# March 2022 
#=============================================================

#TODO
#add estimation routines
#

#install samsim 
#remotes::install_github("Pacific-salmon-assess/samSim", ref="timevar", force=TRUE)

#install samest
#remotes::install_git('https://github.com/Pacific-salmon-assess/samEst')

library(samEst)
library(samSim)
library(ggplot2)
library(devtools)
library(gridExtra)
library(dplyr)
library(here)
#source("sgen_functions.R")
source("utils.R")
#TODO estimate only for 40 yrs of data.

#here::here()
## Load relevant input data
# Simulation run parameters describing different scenarios
#simPar <- read.csv("../data/samsimIFcoho/cohoSimPars_test.csv")
simPar <- read.csv("../data/samsimHarCk/harcnkSimPars.csv")
mat <- read.csv("../data/samsimIFcoho/cohoCorrMat.csv", header=F)

## Store relevant object names to help run simulation 
scenNames <- unique(simPar$scenario)
dirNames <- sapply(scenNames, function(x) paste(x, unique(simPar$species),sep = "_"))

## First check to ensure that a single scenario can be run (only a small number
# of trials necessary)
plotscn <- TRUE
p <- list()
lfoTMB <- list()

for(a in seq_len(nrow(simPar))){
   #a=1
 
  simData <- readRDS(paste0("../test/SamSimOutputs/simData/", simPar$nameOM[a],"/",simPar$scenario[a],"/",
                         paste(simPar$nameOM[a],"_", simPar$nameMP[a], "_", "CUsrDat.RData",sep="")))$srDatout

  lfodf<-matrix(NA,nrow=length(unique(simData$iteration)),ncol=14,
  dimnames = list(unique(simData$iteration),
    c("simple", "autocorr", "rwa_lastparam","rwa_last3paramavg","rwa_last5paramavg",
      "rwb_lastparam","rwb_last3paramavg","rwb_last5paramavg",
      "hmm_regime_pick","hmm_regime_average","hmma_regime_pick",
      "hmma_regime_average","hmmb_regime_pick","hmmb_regime_average",
      )))
 
  for(u in unique(simData$iteration)){
    dat<-simData[simData$iteration==u,]
    dat<-dat[dat$year>(max(dat$year)-46),]

    dat <- dat[!is.na(dat$obsRecruits),]
    length(unique(dat$year))
    df <- data.frame(by=dat$year,
                  S=dat$obsSpawners,
                  R=dat$obsRecruits,
                  logRS=log(dat$obsRecruits/dat$obsSpawners))

    #=======================
    #lfo comparison
    lfostatic<-tmb_mod_lfo_cv(data=df,model='static')
    lfoac <- tmb_mod_lfo_cv(data=df,model='staticAC')
    lfoalpha <- tmb_mod_lfo_cv(data=df,model='alpha', siglfo="obs")
    lfobeta <- tmb_mod_lfo_cv(data=df,model='beta', siglfo="obs")
    lfohmm <- tmb_mod_lfo_cv(data=df,model='HMM')
    lfohmma <- tmb_mod_lfo_cv(data=df,model='HMM_a')
    lfohmmb <- tmb_mod_lfo_cv(data=df,model='HMM_b')

    lfodf[u,] <- c(sum(lfostatic), sum(lfoac), 
      sum(lfoalpha$lastparam), sum(lfoalpha$last3paramavg), sum(lfoalpha$last5paramavg), 
      sum(lfobeta$lastparam), sum(lfobeta$last3paramavg), sum(lfobeta$last5paramavg),
      sum(lfohmm$regime_pick),sum(lfohmm$regime_average),
      sum(lfohmma$regime_pick),sum(lfohmma$regime_average),
      sum(lfohmmb$regime_pick),sum(lfohmmb$regime_average))
  }

lfoTMB[[a]] <-lfodf
}







head(lfodf)
dimnames(lfodf)[[2]]
lfo<-apply(lfodf,1,which.max)

lfochoice<-data.frame(
  chsnmod=dimnames(lfodf)[[2]][apply(lfodf,1,which.max)])

lfochoice$chsnmod<-factor(lfochoice$chsnmod, levels=c("simple", "autocorr", "rwa_lastparam", "rwa_last3paramavg", "rwa_last5paramavg",
 "rwb_lastparam", "rwb_last3paramavg", "rwb_last5paramavg", "hmm_regime_pick", "hmm_regime_average", 
 "hmma_regime_pick", "hmma_regime_average", "hmmb_regime_pick", "hmmb_regime_average"))



ggplot(lfochoice) +  
 geom_bar(aes(chsnmod))+
theme_bw(14)




#==========================
#todo
#test
#add bayes



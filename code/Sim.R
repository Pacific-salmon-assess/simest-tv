#=============================================================
#Run samSim simulations 
#using the harrison chinook data as example
#samSim updates
#Catarina Wor
#October 2022 
#=============================================================

#install samsim 
#remotes::install_github("Pacific-salmon-assess/samSim", ref="timevar", force=TRUE)

#install samest
#remotes::install_git('https://github.com/Pacific-salmon-assess/samEst')

here::here()
library(samSim)
library(ggplot2)
library(devtools)
library(gridExtra)
library(dplyr)
library(here)
source("code/utils.R")

## Load relevant input data
# Simulation run parameters describing different scenarios
simPar <- read.csv(("data/harck/harcnkSimPars.csv"))
# CU-specific parameters
cuPar <- read.csv(("data/harck/harcnkCUPars.csv"))

## Store relevant object names to help run simulation 
scenNames <- unique(simPar$scenario)
dirNames <- sapply(scenNames, function(x) paste(x, unique(simPar$species),sep = "_"))

## First check to ensure that a single scenario can be run (only a small number
# of trials necessary)

#Run and save simulated data

for(a in seq_len(nrow(simPar))){

  
   genericRecoverySim(simPar=simPar[a,], cuPar=cuPar, catchDat=NULL, srDat=NULL,
            variableCU=FALSE, ricPars=NULL, larkPars=NULL, cuCustomCorrMat= NULL,
            outDir="outs", nTrials=100, makeSubDirs=TRUE, random=FALSE, uniqueProd=TRUE,
                               uniqueSurv=FALSE)

}


simData <- readRDS(paste0("../test/SamSimOutputs/simData/", simPar$nameOM[a],"/",simPar$scenario[1],"/",
                         paste(simPar$nameOM[1],"_", simPar$nameMP[1], "_", "CUsrDat.RData",sep="")))$srDatout
  

simData<-simData[simData$CU==1,]

lfodf<-matrix(NA,nrow=length(unique(simData$iteration)),ncol=14,
  dimnames = list(unique(simData$iteration),
    c("simple", "autocorr", "rwa_lastparam","rwa_last3paramavg","rwa_last5paramavg",
      "rwb_lastparam","rwb_last3paramavg","rwb_last5paramavg",
      "hmm_regime_pick","hmm_regime_average","hmma_regime_pick",
      "hmma_regime_average","hmmb_regime_pick","hmmb_regime_average")))

simest<-list()
rmse<-list()

#compiled Bayesian models
simple_mod <- sr_mod(type='static', ac=FALSE, par='n', loglik=FALSE, modelcode=TRUE)
simpleac_mod <- sr_mod(type='static', ac=TRUE, par='n', loglik=FALSE, modelcode=TRUE)
rwa_mod <- sr_mod(type='rw',ac=FALSE,par="a",loglik=FALSE, modelcode=TRUE)
rwb_mod <- sr_mod(type='rw',ac=FALSE,par="b",loglik=FALSE, modelcode=TRUE)
rwab_mod <- sr_mod(type='rw',ac=FALSE,par="both",loglik=FALSE, modelcode=TRUE)
hmma_mod<-sr_mod(type='hmm',ac=FALSE,par="a",loglik=FALSE, modelcode=TRUE)
hmmb_mod<-sr_mod(type='hmm',ac=FALSE,par="b",loglik=FALSE, modelcode=TRUE)
hmmab_mod<-sr_mod(type='hmm',ac=FALSE,par="both",loglik=FALSE, modelcode=TRUE)
hmmabcaphi_mod<-sr_mod(type='hmm',ac=FALSE,par="both",loglik=FALSE, modelcode=TRUE,caphigh=TRUE)


  

  if(plotscn ==TRUE){
    df<-simData[[i]] 
    df<-df[df$iteration==1,]
    stackcu1<-cbind(df[,-c(9,10,11)],stack(df[,9:11]))

    p[[i]] <- ggplot(stackcu1) +
      geom_line(aes(x=year,y=values, col=as.factor(CU)))+
      geom_point(aes(x=year,y=values, col=as.factor(CU)))+
      theme_bw(14)+ theme(legend.position="none")+
      facet_wrap(~ind, scales="free_y")+
      scale_colour_viridis_d() +
      labs(title = simPar$nameOM[i])

    
  }

  simData[[i]] <- simData[[i]] %>% 
    filter(CU == 1, year %in% (nyr-50+1):nyr) %>% 
    filter(!is.na(obsRecruits))%>%
    mutate() %>% 
    rename(byr=year, spwn=spawners, rec=recruits, alpha_true=alpha, beta_true=beta) %>% 
    mutate(scenario = scenNames[i], alpha=99., beta=99., alpha_se=99., beta_se=99.) %>% # cols for output
    select(scenario, everything())  #reorder cols


#ggsave(
#      filename = "../plots/scenarios.pdf", 
#      plot = marrangeGrob(p, nrow=1, ncol=1), 
#      width = 12, height = 5
#    )
  



#==========================
#todo





# add sceatios for changes of age at spawners - later



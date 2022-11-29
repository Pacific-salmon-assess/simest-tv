#=============================================================
#Example use of samsim for time varying simulation evaluation 
#of lfo model selection in TMB
#using the coho data as it is compliant with the most recent
#samSim updates
#Catarina Wor
# October 2022 
#=============================================================


#install samest
#remotes::install_git('https://github.com/Pacific-salmon-assess/samEst')

library(samEst)
library(samSim)
library(ggplot2)
library(devtools)
library(gridExtra)
library(dplyr)

#source("code/utils.R")
#TODO estimate only for 40 yrs of data.

#here::here()
## Load relevant input data
# Simulation run parameters describing different scenario
simPar <- read.csv(("data/harck/harcnkSimPars.csv"))
# CU-specific parameters
cuPar <- read.csv(("data/harck/harcnkCUPars.csv"))
## Store relevant object names to help run simulation 

scenNames <- unique(simPar$scenario)
dirNames <- sapply(scenNames, function(x) paste(x, unique(simPar$species),sep = "_"))

## First check to ensure that a single scenario can be run (only a small number
# of trials necessary)
plotscn <- TRUE
p <- list()
lfoTMB <- list()
lfomwTMB <- list()
aicTMB <- list()
bicTMB <- list()

for(a in seq_len(nrow(simPar))){
   #a=2
 
  simData <- readRDS(paste0("outs/SamSimOutputs/simData/", simPar$nameOM[a],"/",simPar$scenario[a],"/",
                         paste(simPar$nameOM[a],"_", simPar$nameMP[a], "_", "CUsrDat.RData",sep="")))$srDatout

  lfodf<-matrix(NA,nrow=length(unique(simData$iteration)),ncol=20,
  dimnames = list(unique(simData$iteration),
    c("simple", "autocorr", 
      "rwa_last","rwa_last3","rwa_last5",
      "rwb_last","rwb_last3","rwb_last5",
      "rwab_last","rwab_last3","rwab_last5",
      "hmma_last","hmma_last3","hmma_last5",
      "hmmb_last","hmmb_last3","hmmb_last5",
      "hmm_last", "hmm_last3", "hmm_last5"
      )))
  lfomwdf<-matrix(NA,nrow=length(unique(simData$iteration)),ncol=20,
  dimnames = list(unique(simData$iteration),
    c("simple", "autocorr", 
      "rwa_last","rwa_last3","rwa_last5",
      "rwb_last","rwb_last3","rwb_last5",
      "rwab_last","rwab_last3","rwab_last5",
      "hmma_last","hmma_last3","hmma_last5",
      "hmmb_last","hmmb_last3","hmmb_last5",
      "hmm_last", "hmm_last3", "hmm_last5"
      )))
 
   aicdf<-matrix(NA,nrow=length(unique(simData$iteration)),ncol=8,
  dimnames = list(unique(simData$iteration),
    c("simple", "autocorr", 
      "rwa",
      "rwb",
      "rwab",
      "hmma",
      "hmmb",
      "hmm"
      )))

   bicdf<-matrix(NA,nrow=length(unique(simData$iteration)),ncol=8,
  dimnames = list(unique(simData$iteration),
    c("simple", "autocorr", 
      "rwa",
      "rwb",
      "rwab",
      "hmma",
      "hmmb",
      "hmm"
      )))
  
  for(u in unique(simData$iteration)){
    #u=30
    dat<-simData[simData$iteration==u,]
    dat<-dat[dat$year>(max(dat$year)-46),]
    dat <- dat[!is.na(dat$obsRecruits),]
   
    df <- data.frame(by=dat$year,
                  S=dat$obsSpawners,
                  R=dat$obsRecruits,
                  logRS=log(dat$obsRecruits/dat$obsSpawners))

    #=======================
    #lfo comparison
    lfostatic <- tmb_mod_lfo_cv(data=df,model='static', L=10)
    lfoac <- tmb_mod_lfo_cv(data=df,model='staticAC', L=10)
    lfoalpha <- tmb_mod_lfo_cv(data=df,model='rw_a', siglfo="obs", L=10)
    lfobeta <- tmb_mod_lfo_cv(data=df,model='rw_b', siglfo="obs", L=10)
    lfoalphabeta <- tmb_mod_lfo_cv(data=df,model='rw_both', siglfo="obs", L=10)
    lfohmma <- tmb_mod_lfo_cv(data=df,model='HMM_a', L=10)
    lfohmmb <- tmb_mod_lfo_cv(data=df,model='HMM_b', L=10)
    lfohmm <- tmb_mod_lfo_cv(data=df,model='HMM', L=10)
    
    TMBstatic <- ricker_TMB(data=df, priors=1)
    TMBac <- ricker_TMB(data=df, AC=TRUE,priors=1)
    TMBtva <- ricker_rw_TMB(data=df,tv.par='a',priors=1)
    TMBtvb <- ricker_rw_TMB(data=df, tv.par='b',priors=1)
    TMBtvab <- ricker_rw_TMB(data=df, tv.par='both',priors=1)
    TMBhmma <- ricker_hmm_TMB(data=df, tv.par='a',priors=1)
    TMBhmmb <- ricker_hmm_TMB(data=df, tv.par='b',priors=1)
    TMBhmm  <- ricker_hmm_TMB(data=df, tv.par='both',priors=1)

  
    
    LLdf<-rbind(ifelse(is.na(lfostatic$lastparam),-999,lfostatic$lastparam),
      ifelse(is.na(lfoac$lastparam),-999,lfoac$lastparam),
      ifelse(is.na(lfoalpha$lastparam),-999,lfoalpha$lastparam),
      ifelse(is.na(lfoalpha$last3param),-999,lfoalpha$last3param),
      ifelse(is.na(lfoalpha$last5param),-999,lfoalpha$last5param),
      ifelse(is.na(lfobeta$lastparam),-999,lfobeta$lastparam),
      ifelse(is.na(lfobeta$last3param),-999,lfobeta$last3param),
      ifelse(is.na(lfobeta$last5param),-999,lfobeta$last5param),
      ifelse(is.na(lfoalphabeta$lastparam),-999,lfoalphabeta$lastparam),
      ifelse(is.na(lfoalphabeta$last3param),-999,lfoalphabeta$last3param),
      ifelse(is.na(lfoalphabeta$last5param),-999,lfoalphabeta$last5param),
      ifelse(is.na(lfohmma$lastregime_pick),-999,lfohmma$lastregime_pick),
      ifelse(is.na(lfohmma$last3regime_pick),-999,lfohmma$last3regime_pick),
      ifelse(is.na(lfohmma$last5regime_pick),-999,lfohmma$last5regime_pick),
      ifelse(is.na(lfohmmb$lastregime_pick),-999,lfohmmb$lastregime_pick),
      ifelse(is.na(lfohmmb$last3regime_pick),-999,lfohmmb$last3regime_pick),
      ifelse(is.na(lfohmmb$last5regime_pick),-999,lfohmmb$last5regime_pick),
      ifelse(is.na(lfohmm$lastregime_pick),-999,lfohmm$lastregime_pick),
      ifelse(is.na(lfohmm$last3regime_pick),-999,lfohmm$last3regime_pick),
      ifelse(is.na(lfohmm$last5regime_pick),-999,lfohmm$last5regime_pick)
      )
    rownames(LLdf)<-colnames(lfodf)
    
    
    convdf<-rbind(lfostatic$conv_problem,
      lfoac$conv_problem,
      lfoalpha$conv_problem,
      lfoalpha$conv_problem,
      lfoalpha$conv_problem,
      lfobeta$conv_problem,
      lfobeta$conv_problem,
      lfobeta$conv_problem,
      lfoalphabeta$conv_problem,
      lfoalphabeta$conv_problem,
      lfoalphabeta$conv_problem,
      lfohmma$conv_problem,
      lfohmma$conv_problem,
      lfohmma$conv_problem,
      lfohmmb$conv_problem,
      lfohmmb$conv_problem,
      lfohmmb$conv_problem,
      lfohmm$conv_problem,
      lfohmm$conv_problem,
      lfohmm$conv_problem
      )

    
    mw <- model_weights(LLdf, form='PBMA',type='full')
    mw[as.logical(apply(convdf,1,sum))]<-0
    lfomwdf[u,] <- mw
    lfodf[u,] <- c(ifelse(sum(lfostatic$conv_problem)>0,999,sum(lfostatic$lastparam)), 
      ifelse(sum(lfoac$conv_problem)>0,999,sum(lfoac$lastparam)), 
      ifelse(sum(lfoalpha$conv_problem)>0,999,sum(lfoalpha$lastparam)), 
      ifelse(sum(lfoalpha$conv_problem)>0,999,sum(lfoalpha$last3paramavg)), 
      ifelse(sum(lfoalpha$conv_problem)>0,999,sum(lfoalpha$last5paramavg)), 
      ifelse(sum(lfobeta$conv_problem)>0,999,sum(lfobeta$lastparam)), 
      ifelse(sum(lfobeta$conv_problem)>0,999,sum(lfobeta$last3paramavg)), 
      ifelse(sum(lfobeta$conv_problem)>0,999,sum(lfobeta$last5paramavg)), 
      ifelse(sum(lfoalphabeta$conv_problem)>0,999,sum(lfoalphabeta$lastparam)), 
      ifelse(sum(lfoalphabeta$conv_problem)>0,999,sum(lfoalphabeta$last3paramavg)), 
      ifelse(sum(lfoalphabeta$conv_problem)>0,999,sum(lfoalphabeta$last5paramavg)),    
      ifelse(sum(lfohmma$conv_problem)>0,999,sum(lfohmma$lastregime_pick)),
      ifelse(sum(lfohmma$conv_problem)>0,999,sum(lfohmma$last3regime_pick)),
      ifelse(sum(lfohmma$conv_problem)>0,999,sum(lfohmma$last5regime_pick)),
      ifelse(sum(lfohmmb$conv_problem)>0,999,sum(lfohmmb$lastregime_pick)),
      ifelse(sum(lfohmmb$conv_problem)>0,999,sum(lfohmmb$last3regime_pick)),
      ifelse(sum(lfohmmb$conv_problem)>0,999,sum(lfohmmb$last5regime_pick)),
      ifelse(sum(lfohmm$conv_problem)>0,999,sum(lfohmm$lastregime_pick)),
      ifelse(sum(lfohmm$conv_problem)>0,999,sum(lfohmm$last3regime_pick)),
      ifelse(sum(lfohmm$conv_problem)>0,999,sum(lfohmm$last5regime_pick))
      )

    aicdf[u,]<-c(ifelse(TMBstatic$conv_problem,999,TMBstatic$AICc),
      ifelse(TMBac$conv_problem,999,TMBac$AICc),
      ifelse(TMBtva$conv_problem,999,TMBtva$AICc),
      ifelse(TMBtvb$conv_problem,999,TMBtvb$AICc),
      ifelse(TMBtvab$conv_problem,999,TMBtvab$AICc),
      ifelse(TMBhmma$conv_problem,999,TMBhmma$AICc),
      ifelse(TMBhmmb$conv_problem,999,TMBhmmb$AICc),
      ifelse(TMBhmm$conv_problem,999,TMBhmm$AICc))

    bicdf[u,]<-c(ifelse(TMBstatic$conv_problem,999,TMBstatic$BIC),
      ifelse(TMBac$conv_problem,999,TMBac$BIC),
      ifelse(TMBtva$conv_problem,999,TMBtva$BIC),
      ifelse(TMBtvb$conv_problem,999,TMBtvb$BIC),
      ifelse(TMBtvab$conv_problem,999,TMBtvab$BIC),
      ifelse(TMBhmma$conv_problem,999,TMBhmma$BIC),
      ifelse(TMBhmmb$conv_problem,999,TMBhmmb$BIC),
      ifelse(TMBhmm$conv_problem,999,TMBhmm$BIC))
    

  }

  lfoTMB[[a]] <- lfodf
  lfomwTMB[[a]] <- lfomwdf
  aicTMB[[a]] <- aicdf
  bicTMB[[a]] <- bicdf
}

  
if(!file.exists("outs/simestlfo")){
  dir.create("outs/simestlfo") 
}

save(lfoTMB, lfomwTMB,aicTMB,bicTMB,file="outs/simest/simestlfo_prodcapscenarios.Rdata")




#todo
#processing of lfo output
lfochoicel<-list()
aicchoicel<-list()
bicchoicel<-list()


for(a in seq_len(nrow(simPar))){

  lfomwdf<- lfomwTMB[[a]] 
  
  
  #fix a bug in the naming
  dimnames(lfomwdf)[[2]]<-c("simple", "autocorr", 
      "rwa_last","rwa_last3","rwa_last5",
      "rwb_last","rwb_last3","rwb_last5",
      "rwab_last","rwab_last3","rwab_last5",
      "hmma_last_pick","hmma_last3_pick","hmma_last5_pick",
      "hmma_last_average","hmma_last3_average","hmma_last5_average",
      "hmmb_last_pick","hmmb_last3_pick","hmmb_last5_pick",
      "hmmb_last_average","hmmb_last3_average","hmmb_last5_average",
      "hmm_last_pick", "hmm_last3_pick", "hmm_last5_pick",
      "hmm_last_average","hmm_last3_average","hmm_last5_average"
      )

  lfo<-apply(lfomwdf,1,which.max)

  
  lfochoice<-data.frame(
    chsnmodorig=dimnames(lfomwdf)[[2]][apply(lfomwdf,1,which.max)])
  
  lfochoice$chsnmod <- dplyr::recode(lfochoice$chsnmodorig, 
      "rwa_last"="rwa","rwa_last3" ="rwa","rwa_last5"="rwa",
      "rwb_last"="rwb","rwb_last3"="rwb","rwb_last5"="rwb",
      "rwab_last"="rwab","rwab_last3"="rwab","rwab_last5"="rwab",
      "hmma_last_pick"="hmma","hmma_last3_pick"="hmma","hmma_last5_pick"="hmma",
      "hmma_last_average"="hmma","hmma_last3_average"="hmma","hmma_last5_average"="hmma",
      "hmmb_last_pick"="hmmb","hmmb_last3_pick"="hmmb","hmmb_last5_pick"="hmmb",
      "hmmb_last_average"="hmmb","hmmb_last3_average"="hmmb","hmmb_last5_average"="hmmb",
      "hmm_last_pick"="hmm", "hmm_last3_pick"="hmm", "hmm_last5_pick"="hmm",
      "hmm_last_average"="hmm","hmm_last3_average"="hmm","hmm_last5_average"="hmm")    

  

  lfochoice$chsnmod<-factor(lfochoice$chsnmod, levels=c("simple", "autocorr", 
      "rwa",
      "rwb",
      "rwab",
      "hmma",
      "hmmb",
      "hmm"
      ))
  lfochoice$scenario<- simPar$nameOM[a]
  lfochoice$method <-"LFO"
  lfochoicel[[a]]<-lfochoice

  aicdf<-aicTMB[[a]]
  
  dimnames(aicdf)[[2]]<-c("simple", "autocorr", 
      "rwa",
      "rwb",
      "rwab",
      "hmma",
      "hmmb",
      "hmm"
      )

  
  aicchoice<-data.frame(
    chsnmod=dimnames(aicdf)[[2]][apply(aicdf,1,which.min)])

  aicchoice$chsnmod<-factor(aicchoice$chsnmod, levels=c("simple", "autocorr", 
      "rwa",
      "rwb",
      "rwab",
      "hmma",
      "hmmb",
      "hmm"
      ))
  aicchoice$scenario<- simPar$nameOM[a]
  aicchoice$method<- "AICc"
   aicchoice$chsnmodorig<-NA
  aicchoicel[[a]]<-aicchoice

  bicdf<-bicTMB[[a]]
  
  dimnames(bicdf)[[2]]<-c("simple", "autocorr", 
      "rwa",
      "rwb",
      "rwab",
      "hmma",
      "hmmb",
      "hmm"
      )

  
  bicchoice<-data.frame(
    chsnmod=dimnames(bicdf)[[2]][apply(bicdf,1,which.min)])
  bicchoice$chsnmodorig<-NA
  bicchoice$chsnmod<-factor(bicchoice$chsnmod, levels=c("simple", "autocorr", 
      "rwa",
      "rwb",
      "rwab",
      "hmma",
      "hmmb",
      "hmm"
      ))
  bicchoice$scenario<- simPar$nameOM[a]
  bicchoice$method<- "BIC"

  bicchoicel[[a]]<-bicchoice

}

dflfo <- do.call("rbind", lfochoicel)
dfaic <- do.call("rbind", aicchoicel)
dfbic <- do.call("rbind", bicchoicel)

df<-rbind(dflfo,dfaic,dfbic)


modsel<-ggplot(df) +  
 geom_bar(aes(chsnmod,fill=method), 
    position = position_dodge(width = 0.9, preserve = "single"))+
 facet_wrap(~scenario)+ 
 theme_bw(14) +
 theme(axis.text.x = element_text(angle = 90))+
 scale_fill_viridis_d(begin=.3, end=.9) 



ggsave(
      filename = "outs/SamSimOutputs/plotcheck/model_selectionLFOallopt.pdf", 
      plot = modsel, 
      width = 12, height = 5
    )


#==========================
#todo
#test
#add bayes



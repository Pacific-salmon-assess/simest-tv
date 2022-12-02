# Choose model with model selection criteria
#save pbias for that model
# boxplot of mean pbias by model.

#pseudo code
#1. Simulate
#2.fit all models
#3. run selection criteria
#4. save parameter estimates only for selected model
#5.do boxplots of %bias across all models, sum of numbers above bars = nsim

#install samest
#remotes::install_git('https://github.com/Pacific-salmon-assess/samEst', force=TRUE)



library(samEst)
library(samSim)
library(ggplot2)
library(devtools)
library(gridExtra)
library(dplyr)

source("code/lfopbias_func.R")
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

selestall<-list()



for(a in seq_len(nrow(simPar))){
   #a=2
 
  simData <- readRDS(paste0("outs/SamSimOutputs/simData/", simPar$nameOM[a],"/",simPar$scenario[a],"/",
                         paste(simPar$nameOM[a],"_", simPar$nameMP[a], "_", "CUsrDat.RData",sep="")))$srDatout

  
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

    selest<-list()

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
    
    mw <- model_weights(LLdf, form='PBMA',type='full')
    
    lfomwdf[u,] <- mw

    lfoch <- dimnames(lfomwdf)[[2]][which.max(mw)]
    
    paramslfo<-modselecfit(ch=lfoch,df=df)

    TMBstatic <- ricker_TMB(data=df, priors=1)
    TMBac <- ricker_TMB(data=df, AC=TRUE,priors=1)
    TMBtva <- ricker_rw_TMB(data=df,tv.par='a',priors=1)
    TMBtvb <- ricker_rw_TMB(data=df, tv.par='b',priors=1)
    TMBtvab <- ricker_rw_TMB(data=df, tv.par='both',priors=1)
    TMBhmma <- ricker_hmm_TMB(data=df, tv.par='a',priors=1)
    TMBhmmb <- ricker_hmm_TMB(data=df, tv.par='b',priors=1)
    TMBhmm  <- ricker_hmm_TMB(data=df, tv.par='both',priors=1)


    aicdf[u,]<-c(ifelse(TMBstatic$conv_problem,999,TMBstatic$AICc),
      ifelse(TMBac$conv_problem,999,TMBac$AICc),
      ifelse(TMBtva$conv_problem,999,TMBtva$AICc),
      ifelse(TMBtvb$conv_problem,999,TMBtvb$AICc),
      ifelse(TMBtvab$conv_problem,999,TMBtvab$AICc),
      ifelse(TMBhmma$conv_problem,999,TMBhmma$AICc),
      ifelse(TMBhmmb$conv_problem,999,TMBhmmb$AICc),
      ifelse(TMBhmm$conv_problem,999,TMBhmm$AICc))
    
    aicch <-dimnames(aicdf)[[2]][which.min(aicdf[u,])]
    paramsaic<-modselecfit(ch=aicch,df=df)


    bicdf[u,]<-c(ifelse(TMBstatic$conv_problem,999,TMBstatic$BIC),
      ifelse(TMBac$conv_problem,999,TMBac$BIC),
      ifelse(TMBtva$conv_problem,999,TMBtva$BIC),
      ifelse(TMBtvb$conv_problem,999,TMBtvb$BIC),
      ifelse(TMBtvab$conv_problem,999,TMBtvab$BIC),
      ifelse(TMBhmma$conv_problem,999,TMBhmma$BIC),
      ifelse(TMBhmmb$conv_problem,999,TMBhmmb$BIC),
      ifelse(TMBhmm$conv_problem,999,TMBhmm$BIC))

    bicch <-dimnames(bicdf)[[2]][which.min(bicdf[u,])]
    paramsbic<-modselecfit(ch=bicch,df=df)


    dfa<- data.frame(parameter="alpha",
      iteration=u,
      method=rep(c("lfo","aic","bic"),each=length(dat$year)),
      model=rep(c(lfoch,aicch,bicch),each=length(dat$year)),
      by=dat$year,
      sim=dat$alpha,
      est=c(paramslfo$alpha,paramsaic$alpha,paramsbic$alpha))
    dfa$pbias<- ((dfa$est-dfa$sim)/dfa$sim)*100

    dfsmax<-data.frame(parameter="Smax",
      iteration=u,
      method=rep(c("lfo","aic","bic"),each=length(dat$year)),
      model=rep(c(lfoch,aicch,bicch),each=length(dat$year)),
      by=dat$year,
      sim=1/dat$beta,
      est=c(paramslfo$smax,paramsaic$smax,paramsbic$smax))
    dfsmax$pbias<- ((dfa$est-dfa$sim)/dfa$sim)*100
  
    dfsig<-data.frame(parameter="sig",
      iteration=u,
      method=rep(c("lfo","aic","bic"),each=length(dat$year)),
      model=rep(c(lfoch,aicch,bicch),each=length(dat$year)),
      by=dat$year,
      sim=dat$sigma,
      est=c(paramslfo$sig,paramsaic$sig,paramsbic$sig))
    dfsig$pbias<- ((dfa$est-dfa$sim)/dfa$sim)*100
    

    smsysim<-smsyCalc(dat$alpha,dat$beta)
    dfsmsy <- data.frame(parameter="Smsy",
      iteration=u,
      method=rep(c("lfo","aic","bic"),each=length(dat$year)),
      model=rep(c(lfoch,aicch,bicch),each=length(dat$year)),
      by=dat$year,
      sim=smsysim,
      est=c(paramslfo$smsy,paramsaic$smsy,paramsbic$smsy))
      
    dfsmsy$pbias<- ((dfa$est-dfa$sim)/dfa$sim)*100
    

    dfsgen<- data.frame(parameter="Smsy",
      iteration=u,
      method=rep(c("lfo","aic","bic"),each=length(dat$year)),
      model=rep(c(lfoch,aicch,bicch),each=length(dat$year)),
      by=dat$year,
      sim= unlist(mapply(sGenCalc,a=dat$alpha,Smsy=smsysim, b=dat$beta)),
      est=c(paramslfo$sgen,paramsaic$sgen,paramsbic$sgen))
    dfsig$pbias<- ((dfa$est-dfa$sim)/dfa$sim)*100
    
   
    dfumsy<-data.frame(parameter="Smsy",
      iteration=u,
      method=rep(c("lfo","aic","bic"),each=length(dat$year)),
      model=rep(c(lfoch,aicch,bicch),each=length(dat$year)),
      by=dat$year,
      sim= umsyCalc(dat$alpha),
      est=c(paramslfo$umsy,paramsaic$umsy,paramsbic$umsy))
     
    dfsig$pbias<- ((dfa$est-dfa$sim)/dfa$sim)*100
    
    selest[[u]]<-rbind(dfa,dfsmax,dfsig,dfsmsy,dfsgen,dfumsy)

  }

  selestall[[a]] <- selest
}
       



      
      
      
      
 
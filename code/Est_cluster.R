

packages <- c("ggplot2", "remotes", "dplyr", "here", "rstan", "TMB")
install.packages(setdiff(packages, rownames(installed.packages()))) 

#install samsim
#load samsim
#setwd("/home/caw001/Documents/samSim/R")
#devtools::document()
#devtools::load_all()
#remotes::install_github("Pacific-salmon-assess/samSim", ref="timevar", force=TRUE, upgrade="always" )

#install samest
#setwd("/home/caw001/Documents/samEst/R")
#devtools::document()
#devtools::load_all()
#remotes::install_git('https://github.com/Pacific-salmon-assess/samEst', force=TRUE)

setwd("/home/caw001/Documents/simest-tv")


library(ggplot2)
library(dplyr)
library(rstan)
library(TMB)
library(parallel)
source("code/utils.R")

rstan_options(auto_write = TRUE)
ncores<-detectCores()-1

options(mc.cores = parallel::detectCores())



# Simulation run parameters describing different scenarios
simPar <- read.csv("data/HarCk/harcnkSimPars.csv")

## Store relevant object names to help run simulation 
scenNames <- unique(simPar$scenario)
dirNames <- sapply(scenNames, function(x) paste(x, unique(simPar$species),sep = "_"))

## First check to ensure that a single scenario can be run (only a small number
# of trials necessary)
simData <- list()




#compiled Bayesian models
simple_mod <- compile_code(type='static', ac=FALSE, par='n')
simpleac_mod <- compile_code(type='static', ac=TRUE, par='n')
rwa_mod <- compile_code(type='rw',ac=FALSE,par="a")
rwb_mod <- compile_code(type='rw',ac=FALSE,par="b")
rwab_mod <- compile_code(type='rw',ac=FALSE,par="both")
hmma_mod <- compile_code(type='hmm',ac=FALSE,par="a")
hmmb_mod <- compile_code(type='hmm',ac=FALSE,par="b")
hmmab_mod <- compile_code(type='hmm',ac=FALSE,par="both")
hmmabcaphi_mod <- compile_code(type='hmm',ac=FALSE,par="both",caphigh=TRUE)



simresult<-mclapply(mc.cores=ncores,
  X=seq_len(nrow(simPar)),
  FUN=function(a){
    simData <- readRDS(paste0("outs/SamSimOutputs/simData/", simPar$nameOM[a],"/",simPar$scenario[a],"/",
                         paste(simPar$nameOM[a],"_", simPar$nameMP[a], "_", "CUsrDat.RData",sep="")))$srDatout
     
    out<-mclapply(mc.cores=ncores,
      X=unique(simData$iteration),
      FUN=function(u){
        dat<-simData[simData$iteration==u,]
        dat<-dat[dat$year>(max(dat$year)-46),]
        dat <- dat[!is.na(dat$obsRecruits),]
        df <- data.frame(by=dat$year,
                        S=dat$obsSpawners,
                        R=dat$obsRecruits,
                        logRS=log(dat$obsRecruits/dat$obsSpawners))
 

        p <- ricker_TMB(data=df)
        pac <- ricker_TMB(data=df, AC=TRUE)
        ptva <- ricker_rw_TMB(data=df,tv.par='a')
        ptvb <- ricker_rw_TMB(data=df, tv.par='b')
        ptvab <- ricker_rw_TMB(data=df, tv.par='both')
        phmma <- ricker_hmm_TMB(data=df, tv.par='a')
        phmmb <- ricker_hmm_TMB(data=df, tv.par='b')
        phmm  <- ricker_hmm_TMB(data=df, tv.par='both')
        
        b <- ricker_stan(data=df,iter = 800, mod=simple_mod)
        #ricker autocorr
        bac <- ricker_stan(data=df,iter = 800, AC=TRUE, mod=simpleac_mod)
        #ricker tva
        btva <- ricker_rw_stan(data=df, par="a",iter = 800, mod=rwa_mod)
        #ricker tvb
        btvb <- ricker_rw_stan(data=df, par="b",iter = 800, mod=rwb_mod)
        #ricker tvab
        btvab <- ricker_rw_stan(data=df, par="both",iter = 800, mod=rwab_mod) 
        #ricker tvhmma
        bhmma <- ricker_hmm_stan(data=df, par="a",iter = 800, mod=hmma_mod)
        #ricker tvhmmb
        bhmmb <- ricker_hmm_stan(data=df, par="b",iter = 800, mod=hmmb_mod)
        #ricker tvhmmab
        bhmmab <- ricker_hmm_stan(data=df, par="both",iter = 800, mod=hmmab_mod) 
        #ricker tvhmmab capacity high
        bhmmabcaphi <- ricker_hmm_stan(data=df, par="both",iter = 800, mod=hmmabcaphi_mod)
  
        #a estimates
        dfa<- data.frame(parameter="alpha",
          iteration=u,
          method=rep(c(rep("MLE",10),rep("MCMC",12)),each=nrow(df)),
          model=rep(c("simple","autocorr","rwa","rwb","rwab","hmma_regime","hmma_average",
              "hmmb_regime","hmmab_regime","hmmab_average","simple","autocorr",
              "rwa","rwb","rwab","hmma_regime","hmma_average",
              "hmmb_regime","hmmab_regime","hmmab_average", "hmmabhc_regime","hmmabhc_average"),each=nrow(df)),
          by=rep(dat$year,22),
          sim=rep(dat$alpha,22),
          est=c(rep(p$alpha,nrow(df)),
              rep(pac$alpha,nrow(df)),
              ptva$alpha,
              rep(ptvb$alpha,nrow(df)),
              ptvab$alpha,
              phmma$alpha[phmma$regime],
              c(phmma$alpha%*%phmma$probregime),
              rep(phmmb$alpha,nrow(df)),
              phmm$alpha[phmm$regime],
              phmm$alpha%*%phmm$probregime,
              rep(b$alpha,nrow(df)),
              rep(bac$alpha,nrow(df)),
              btva$alpha[-1],
              rep(btvb$alpha,nrow(df)),
              btvab$alpha[-1],
              bhmma$alpha_regime,
              bhmma$alpha_wgt,
              rep(bhmmb$alpha,nrow(df)),
              bhmmab$alpha_regime,
              bhmmab$alpha_wgt,
              bhmmabcaphi$alpha_regime,
              bhmmabcaphi$alpha_wgt
              ),
          convergence=c(rep(c(p$model$convergence + p$conv_problem,
              pac$model$convergence + pac$conv_problem,
              ptva$model$convergence + ptva$conv_problem,
              ptvb$model$convergence + ptvb$conv_problem,
              ptvab$model$convergence + ptvab$conv_problem,
              phmma$model$convergence + phmma$conv_problem,
              phmma$model$convergence + phmma$conv_problem,
              phmmb$model$convergence + phmmb$conv_problem,
              phmm$model$convergence + phmm$conv_problem,
              phmm$model$convergence + phmm$conv_problem,
              as.numeric(abs(b$mcmcsummary["log_a","Rhat"]-1)>.1),
              as.numeric(abs(bac$mcmcsummary["log_a","Rhat"]-1)>.1)
              ),each=nrow(df)),
              as.numeric(abs(btva$mcmcsummary[grep("log_a\\[",rownames(btva$mcmcsummary)),"Rhat"]-1)>.1),
              rep(as.numeric(abs(btvb$mcmcsummary[grep("log_a",rownames(btvb$mcmcsummary)),"Rhat"]-1)>.1),nrow(df)),
              as.numeric(abs(btvab$mcmcsummary[grep("log_a\\[",rownames(btvab$mcmcsummary)),"Rhat"]-1)>.1),
              #hmma pick
              c(as.numeric(abs(bhmma$mcmcsummary[grep("log_a\\[",rownames(bhmma$mcmcsummary)),
               "Rhat"]-1)>.1)[bhmma$mcmcsummary[grep("zstar\\[",rownames(bhmma$mcmcsummary)),"50%"]]+
               as.numeric(abs(bhmma$mcmcsummary[grep("zstar\\[",rownames(bhmma$mcmcsummary)),"Rhat"]-1)>.1)),
              #hmma avg
              sum(as.numeric(abs(bhmma$mcmcsummary[grep("log_a\\[",rownames(bhmma$mcmcsummary)),"Rhat"]-1)>.1))+
               sumpair(as.numeric(abs(bhmma$mcmcsummary[grep("^gamma\\[",rownames(bhmma$mcmcsummary)),"Rhat"]-1)>.1)),
              #hmmb 
              rep(as.numeric(abs(bhmmb$mcmcsummary[grep("log_a",rownames(bhmmb$mcmcsummary)),"Rhat"]-1)>.1),nrow(df)),
              #hmmab pick
              c(as.numeric(abs(bhmmab$mcmcsummary[grep("log_a\\[",rownames(bhmmab$mcmcsummary)),
                "Rhat"]-1)>.1)[bhmmab$mcmcsummary[grep("zstar\\[",rownames(bhmmab$mcmcsummary)),"50%"]]+
                as.numeric(abs(bhmmab$mcmcsummary[grep("zstar\\[",rownames(bhmmab$mcmcsummary)),"Rhat"]-1)>.1)),
              #hmmab avg
              sum(as.numeric(abs(bhmmab$mcmcsummary[grep("log_a\\[",rownames(bhmmab$mcmcsummary)),"Rhat"]-1)>.1))+
               sumpair(as.numeric(abs(bhmmab$mcmcsummary[grep("^gamma\\[",rownames(bhmmab$mcmcsummary)),"Rhat"]-1)>.1)),
              #hmmab caphig pick
              c(as.numeric(abs(bhmmabcaphi$mcmcsummary[grep("log_a\\[",rownames(bhmmabcaphi$mcmcsummary)),
                "Rhat"]-1)>.1)[bhmmabcaphi$mcmcsummary[grep("zstar\\[",rownames(bhmmabcaphi$mcmcsummary)),"50%"]]+
                as.numeric(abs(bhmmabcaphi$mcmcsummary[grep("zstar\\[",rownames(bhmmabcaphi$mcmcsummary)),"Rhat"]-1)>.1)),
              #hmmab caphig avg
             sum(as.numeric(abs(bhmmabcaphi$mcmcsummary[grep("log_a\\[",rownames(bhmmabcaphi$mcmcsummary)),"Rhat"]-1)>.1))+
              sumpair(as.numeric(abs(bhmmabcaphi$mcmcsummary[grep("^gamma\\[",rownames(bhmmabcaphi$mcmcsummary)),"Rhat"]-1)>.1))))
        dfa$pbias<- ((dfa$est-dfa$sim)/dfa$sim)*100

        rmsea<-aggregate((dfa$est-dfa$sim)^2, list(model=dfa$model,method=dfa$method), function(x)sqrt(mean(x)))
        rmsea$iteration<-u
        rmsea$parameter<-"alpha"
        rmsea$convergence <- aggregate(dfa$convergence, list(model=dfa$model,method=dfa$method), function(x)sum(x,na.rm=T))$x
        
        #Smax
        dfsmax<- data.frame(parameter="Smax",
          iteration=u,
          method=rep(c(rep("MLE",10),rep("MCMC",12)),each=nrow(df)),
          model=rep(c("simple","autocorr","rwa","rwb","rwab","hmma_regime",
              "hmmb_regime","hmmb_average","hmmab_regime","hmmab_average",
              "simple","autocorr","rwa","rwb","rwab","hmma_regime",
              "hmmb_regime","hmmb_average","hmmab_regime","hmmab_average",
              "hmmabhc_regime","hmmabhc_average"),each=nrow(df)),
              by=rep(dat$year,22),
          sim=rep(1/dat$beta,22),
          est=c(rep(p$Smax,nrow(df)),
              rep(pac$Smax,nrow(df)),
              rep(ptva$Smax,nrow(df)),
              ptvb$Smax,
              ptvab$Smax,
              rep(phmma$Smax,nrow(df)),
              phmmb$Smax[phmmb$regime],
              phmmb$Smax%*%phmmb$probregime,
              phmm$Smax[phmm$regime],
              phmm$Smax%*%phmm$probregime,
              rep(b$Smax,nrow(df)),
              rep(bac$Smax,nrow(df)),
              rep(btva$Smax,nrow(df)),
              btvb$Smax,
              btvab$Smax, 
              rep(bhmma$Smax,nrow(df)),
              bhmmb$Smax_regime,
              bhmmb$Smax_wgt,
              bhmmab$Smax_regime,
              bhmmab$Smax_wgt,
              bhmmabcaphi$Smax_regime,
              bhmmabcaphi$Smax_wgt),
          convergence=c(rep(c(p$model$convergence + p$conv_problem,
              pac$model$convergence + pac$conv_problem,
              ptva$model$convergence + ptva$conv_problem,
              ptvb$model$convergence + ptvb$conv_problem,
              ptvab$model$convergence + ptvab$conv_problem,
              phmma$model$convergence + phmma$conv_problem,
              phmmb$model$convergence + phmmb$conv_problem,
              phmmb$model$convergence + phmmb$conv_problem,
              phmm$model$convergence + phmm$conv_problem,
              phmm$model$convergence + phmm$conv_problem,
              as.numeric(abs(b$mcmcsummary["S_max","Rhat"]-1)>.1),
              as.numeric(abs(bac$mcmcsummary["S_max","Rhat"]-1)>.1),
              as.numeric(abs(btva$mcmcsummary["S_max","Rhat"]-1)>.1)
              ),each=nrow(df)),
              as.numeric(abs(btvb$mcmcsummary[grep("S_max",rownames(btvb$mcmcsummary)),"Rhat"]-1)>.1),
              as.numeric(abs(btvab$mcmcsummary[grep("S_max",rownames(btvab$mcmcsummary)),"Rhat"]-1)>.1),     
              rep(as.numeric(abs(bhmma$mcmcsummary[grep("S_max",rownames(bhmma$mcmcsummary)),"Rhat"]-1)>.1),nrow(df)),
              #hmmb pick
              c(as.numeric(abs(bhmmb$mcmcsummary[grep("S_max\\[",rownames(bhmmb$mcmcsummary)),
              "Rhat"]-1)>.1)[bhmmb$mcmcsummary[grep("zstar\\[",rownames(bhmmb$mcmcsummary)),"50%"]]+
              as.numeric(abs(bhmmb$mcmcsummary[grep("zstar\\[",rownames(bhmmb$mcmcsummary)),"Rhat"]-1)>.1)),
              #hmmb avg
              sum(as.numeric(abs(bhmmb$mcmcsummary[grep("S_max\\[",rownames(bhmmb$mcmcsummary)),"Rhat"]-1)>.1))+
               sumpair(as.numeric(abs(bhmmb$mcmcsummary[grep("^gamma\\[",rownames(bhmmb$mcmcsummary)),"Rhat"]-1)>.1)),
              #hmmab pick 
              c(as.numeric(abs(bhmmab$mcmcsummary[grep("S_max\\[",rownames(bhmmab$mcmcsummary)),
                "Rhat"]-1)>.1)[bhmmab$mcmcsummary[grep("zstar\\[",rownames(bhmmab$mcmcsummary)),"50%"]]+
                as.numeric(abs(bhmmab$mcmcsummary[grep("zstar\\[",rownames(bhmmab$mcmcsummary)),"Rhat"]-1)>.1)),
              #hmmab avg
              sum(as.numeric(abs(bhmmab$mcmcsummary[grep("S_max\\[",rownames(bhmmab$mcmcsummary)),"Rhat"]-1)>.1))+
               sumpair(as.numeric(abs(bhmmab$mcmcsummary[grep("^gamma\\[",rownames(bhmmab$mcmcsummary)),"Rhat"]-1)>.1)),
              #hmmab pick caphi
              c(as.numeric(abs(bhmmabcaphi$mcmcsummary[grep("S_max\\[",rownames(bhmmabcaphi$mcmcsummary)),
                "Rhat"]-1)>.1)[bhmmabcaphi$mcmcsummary[grep("zstar\\[",rownames(bhmmabcaphi$mcmcsummary)),"50%"]]+
                as.numeric(abs(bhmmabcaphi$mcmcsummary[grep("zstar\\[",rownames(bhmmabcaphi$mcmcsummary)),"Rhat"]-1)>.1)),
              #hmmab avg caphi
              sum(as.numeric(abs(bhmmabcaphi$mcmcsummary[grep("S_max\\[",rownames(bhmmabcaphi$mcmcsummary)),"Rhat"]-1)>.1))+
              sumpair(as.numeric(abs(bhmmabcaphi$mcmcsummary[grep("^gamma\\[",rownames(bhmmabcaphi$mcmcsummary)),"Rhat"]-1)>.1)))
            )
        dfsmax$pbias<- ((dfsmax$est-dfsmax$sim)/dfsmax$sim)*100

        rmsesmax<-aggregate((dfsmax$est-dfsmax$sim)^2, 
            list(model=dfsmax$model,method=dfsmax$method), function(x)sqrt(mean(x)))
        rmsesmax$iteration<-u
        rmsesmax$parameter<-"Smax"
        rmsesmax$convergence <- aggregate(dfsmax$convergence, list(model=dfsmax$model,method=dfsmax$method), function(x)sum(x,na.rm=T))$x
         
        #sigma
        dfsig<- data.frame(parameter="sigma",
          iteration=u,
          method=rep(c(rep("MLE",8),rep("MCMC",9)),each=nrow(df)),
          model=rep(c("simple","autocorr","rwa","rwb","rwab","hmma_regime",
              "hmmb_regime","hmmab_regime","simple","autocorr","rwa","rwb","rwab",
              "hmma_regime","hmmb_regime","hmmab_regime","hmmabhc_regime"),each=nrow(df)),
              by=rep(dat$year,17),
          sim=rep(dat$sigma,17),
          est=c(rep(p$sig,nrow(df)),
              rep(pac$sig,nrow(df)),
              rep(ptva$sig,nrow(df)),
              rep(ptvb$sig,nrow(df)),
              rep(ptvab$sig,nrow(df)),
              rep(phmma$sigma,nrow(df)),
              rep(phmmb$sigma,nrow(df)),
              rep(phmm$sigma,nrow(df)),
              rep(b$sigobs,nrow(df)),
              rep(bac$sigobs,nrow(df)),
              rep(btva$sigobs,nrow(df)),
              rep(btvb$sigobs,nrow(df)),
              rep(btvab$sigobs,nrow(df)),
              rep(bhmma$sigobs,nrow(df)),
              rep(bhmmb$sigobs,nrow(df)),
              rep(bhmmab$sigobs,nrow(df)),
              rep(bhmmabcaphi$sigobs,nrow(df))),
          convergence=rep(c(p$model$convergence + p$conv_problem,
              pac$model$convergence + pac$conv_problem,
              ptva$model$convergence + ptva$conv_problem,
              ptvb$model$convergence + ptvb$conv_problem,
              ptvab$model$convergence + ptvab$conv_problem,
              phmma$model$convergence + phmma$conv_problem,
              phmmb$model$convergence + phmmb$conv_problem,
              phmm$model$convergence + phmm$conv_problem,
              as.numeric(abs(b$mcmcsummary["sigma_e","Rhat"]-1)>.1),
              as.numeric(abs(bac$mcmcsummary["sigma_e","Rhat"]-1)>.1),
              as.numeric(abs(btva$mcmcsummary["sigma_e","Rhat"]-1)>.1),
              as.numeric(abs(btvb$mcmcsummary["sigma_e","Rhat"]-1)>.1),
              as.numeric(abs(btvab$mcmcsummary["sigma_e","Rhat"]-1)>.1),
              as.numeric(abs(bhmma$mcmcsummary["sigma","Rhat"]-1)>.1),
              as.numeric(abs(bhmmb$mcmcsummary["sigma","Rhat"]-1)>.1),
              as.numeric(abs(bhmmab$mcmcsummary["sigma","Rhat"]-1)>.1),
              as.numeric(abs(bhmmabcaphi$mcmcsummary["sigma","Rhat"]-1)>.1)
              ),each=nrow(df)))
        dfsig$pbias<- ((dfsig$est-dfsig$sim)/dfsig$sim)*100
  
        rmsesig<-aggregate((dfsig$est-dfsig$sim)^2, list(model=dfsig$model,method=dfsig$method), function(x)sqrt(mean(x)))
        rmsesig$iteration<-u
        rmsesig$parameter<-"sigma"
        rmsesig$convergence <- aggregate(dfsig$convergence, list(model=dfsig$model,method=dfsig$method), function(x)sum(x,na.rm=T))$x
      
        #Smsy
        smsysim<-smsyCalc(dat$alpha,dat$beta)
        
        dfsmsy<- data.frame(parameter="smsy",
          iteration=u,
          method=rep(c(rep("MLE",11),rep("MCMC",13)),each=nrow(df)),
          model=rep(c("simple","autocorr","rwa","rwb","rwab","hmma_regime","hmma_average",
              "hmmb_regime","hmmb_average","hmmab_regime","hmmab_average",
              "simple", "autocorr", "rwa","rwb","rwab","hmma_regime","hmma_average",
              "hmmb_regime","hmmb_average","hmmab_regime","hmmab_average",
              "hmmabhc_regime","hmmabhc_average"),each=nrow(df)),
              by=rep(dat$year,24),
          sim=rep(smsysim,24),
          est=c(rep(p$Smsy,nrow(df)),
                rep(pac$Smsy,nrow(df)),
                ptva$Smsy,
                ptvb$Smsy,
                ptvab$Smsy,
                phmma$Smsy[phmma$regime],
                smsyCalc(c(phmma$alpha%*%phmma$probregime),phmma$beta),
                phmmb$Smsy[phmmb$regime],
                smsyCalc(phmmb$alpha,phmmb$beta%*%phmmb$probregime),
                phmm$Smsy[phmm$regime],
                smsyCalc(phmm$alpha%*%phmm$probregime,phmm$beta%*%phmm$probregime),
                rep(b$Smsy,nrow(df)),
                rep(bac$Smsy,nrow(df)),
                btva$Smsy,
                btvb$Smsy,
                btvab$Smsy,
                bhmma$Smsy_regime,
                bhmma$Smsy_wgt,
                bhmmb$Smsy_regime,
                bhmmb$Smsy_wgt,
                bhmmab$Smsy_regime,
                bhmmab$Smsy_wgt,
                bhmmabcaphi$Smsy_regime,
                bhmmabcaphi$Smsy_wgt
            ),    
          convergence=c(rep(c(p$model$convergence + p$conv_problem,
              pac$model$convergence + pac$conv_problem,
              ptva$model$convergence + ptva$conv_problem,
              ptvb$model$convergence + ptvb$conv_problem,
              ptvab$model$convergence + ptvab$conv_problem,
              phmma$model$convergence + phmma$conv_problem,
              phmma$model$convergence + phmma$conv_problem,
              phmmb$model$convergence + phmmb$conv_problem,
              phmmb$model$convergence + phmmb$conv_problem,
              phmm$model$convergence + phmm$conv_problem,
              phmm$model$convergence + phmm$conv_problem,
              as.numeric(abs(b$mcmcsummary["S_msy","Rhat"]-1)>.1),
              as.numeric(abs(bac$mcmcsummary["S_msy","Rhat"]-1)>.1)),each=nrow(df)),
              as.numeric(abs(btva$mcmcsummary[grep("S_msy",rownames(btva$mcmcsummary)),"Rhat"]-1)>.1),
              as.numeric(abs(btvb$mcmcsummary[grep("S_msy",rownames(btvb$mcmcsummary)),"Rhat"]-1)>.1),
              as.numeric(abs(btvab$mcmcsummary[grep("S_msy",rownames(btvab$mcmcsummary)),"Rhat"]-1)>.1),
              #
              c(as.numeric(abs(bhmma$mcmcsummary[grep("S_msy\\[",rownames(bhmma$mcmcsummary)),
              "Rhat"]-1)>.1)[bhmma$mcmcsummary[grep("zstar\\[",rownames(bhmma$mcmcsummary)),"50%"]]+
              as.numeric(abs(bhmma$mcmcsummary[grep("zstar\\[",rownames(bhmma$mcmcsummary)),"Rhat"]-1)>.1)),
              
              sum(as.numeric(abs(bhmma$mcmcsummary[grep("S_msy\\[",rownames(bhmma$mcmcsummary)),"Rhat"]-1)>.1))+
              sumpair(as.numeric(abs(bhmma$mcmcsummary[grep("^gamma\\[",rownames(bhmma$mcmcsummary)),"Rhat"]-1)>.1)),
              
              c(as.numeric(abs(bhmmb$mcmcsummary[grep("S_msy\\[",rownames(bhmmb$mcmcsummary)),
              "Rhat"]-1)>.1)[bhmmb$mcmcsummary[grep("zstar\\[",rownames(bhmmb$mcmcsummary)),"50%"]]+
              as.numeric(abs(bhmmb$mcmcsummary[grep("zstar\\[",rownames(bhmmb$mcmcsummary)),"Rhat"]-1)>.1)),
              
              sum(as.numeric(abs(bhmmb$mcmcsummary[grep("S_msy\\[",rownames(bhmmb$mcmcsummary)),"Rhat"]-1)>.1))+
              sumpair(as.numeric(abs(bhmmb$mcmcsummary[grep("^gamma\\[",rownames(bhmmb$mcmcsummary)),"Rhat"]-1)>.1)),
              
              c(as.numeric(abs(bhmmab$mcmcsummary[grep("S_msy\\[",rownames(bhmmab$mcmcsummary)),
              "Rhat"]-1)>.1)[bhmmab$mcmcsummary[grep("zstar\\[",rownames(bhmmab$mcmcsummary)),"50%"]]+
              as.numeric(abs(bhmmab$mcmcsummary[grep("zstar\\[",rownames(bhmmab$mcmcsummary)),"Rhat"]-1)>.1)),
              
              sum(as.numeric(abs(bhmmab$mcmcsummary[grep("S_msy\\[",rownames(bhmmab$mcmcsummary)),"Rhat"]-1)>.1))+
              sumpair(as.numeric(abs(bhmmab$mcmcsummary[grep("^gamma\\[",rownames(bhmmab$mcmcsummary)),"Rhat"]-1)>.1)),
              
              c(as.numeric(abs(bhmmabcaphi$mcmcsummary[grep("S_msy\\[",rownames(bhmmabcaphi$mcmcsummary)),
              "Rhat"]-1)>.1)[bhmmabcaphi$mcmcsummary[grep("^zstar\\[",rownames(bhmmabcaphi$mcmcsummary)),"50%"]]+
              as.numeric(abs(bhmmabcaphi$mcmcsummary[grep("^zstar\\[",rownames(bhmmabcaphi$mcmcsummary)),"Rhat"]-1)>.1)),
              
              sum(as.numeric(abs(bhmmabcaphi$mcmcsummary[grep("S_msy\\[",rownames(bhmmabcaphi$mcmcsummary)),"Rhat"]-1)>.1))+
              sumpair(as.numeric(abs(bhmmabcaphi$mcmcsummary[grep("^gamma\\[",rownames(bhmmabcaphi$mcmcsummary)),"Rhat"]-1)>.1)))
        ) 
        dfsmsy$pbias<- ((dfsmsy$est-dfsmsy$sim)/dfsmsy$sim)*100

        rmsesmsy<-aggregate((dfsmsy$est-dfsmsy$sim)^2, list(model=dfsmsy$model,method=dfsmsy$method), function(x)sqrt(mean(x)))
        rmsesmsy$iteration<-u
        rmsesmsy$parameter<-"smsy"
        rmsesmsy$convergence <- aggregate(dfsmsy$convergence, list(model=dfsmsy$model,method=dfsmsy$method), function(x)sum(x,na.rm=T))$x
         
        #Sgen
        
        sgen_b<-median(unlist(mapply(sGenCalc,a=c(b$samples[,,"log_a"]),
            b=c(b$samples[,,"b"]),
            Smsy=c(b$samples[,,"S_msy"]))))
        sgen_bac<-median(unlist(mapply(sGenCalc,a=c(bac$samples[,,"log_a"]),
            b=c(bac$samples[,,"b"]),
            Smsy=bac$samples[,,"S_msy"])))
        
        sgen_tva<-NULL
        sgen_tvb<-NULL
        sgen_tvab<-NULL  

        for(j in seq_len(nrow(dat))){
        
          #tva
          sgen_tva[j]<-median(unlist(mapply(sGenCalc,a=c(btva$samples[,,paste0("log_a[",j,"]")]),
            b=c(btva$samples[,,paste0("b")]),
            Smsy= c(btva$samples[,,paste0("S_msy[",j,"]")]))),na.rm=T)
          
          #tvb    
          sgen_tvb[j]<-median(unlist(mapply(sGenCalc,a=c(btvb$samples[,,paste0("log_a")]),
            b=c(btvb$samples[,,paste0("b[",j,"]")]),
            Smsy= c(btvb$samples[,,paste0("S_msy[",j,"]")]))),na.rm=T)
          
          #tvab        
          sgen_tvab[j]<-median(unlist(mapply(sGenCalc,a=c(btvab$samples[,,paste0("log_a[",j,"]")]),
            b=c(btvab$samples[,,paste0("b[",j,"]")]),
            Smsy=c(btvab$samples[,,paste0("S_msy[",j,"]")]))),na.rm=T)
        }
        
        sgen_hmma <- NULL
        sgen_hmmb <- NULL
        sgen_hmmab <- NULL
        sgen_hmmabhc <- NULL
        kregime <- length(bhmma$alpha) 
  
        for(k in seq_len(kregime)){
          sgen_hmma[k] <- median(unlist(mapply(sGenCalc,a=c(bhmma$samples[,,paste0("log_a[",k,"]")]),
            b=c(bhmma$samples[,,paste0("b")]),
            Smsy=c(bhmma$samples[,,paste0("S_msy[",k,"]")]))),na.rm=T) 
        
          sgen_hmmb[k] <- median(unlist(mapply(sGenCalc,a=c(bhmmb$samples[,,paste0("log_a")]),
            b=c(bhmmb$samples[,,paste0("b[",k,"]")]),
            Smsy=c(bhmmb$samples[,,paste0("S_msy[",k,"]")]))),na.rm=T)
          
          sgen_hmmab[k] <- median(unlist(mapply(sGenCalc,a=c(bhmmab$samples[,,paste0("log_a[",k,"]")]),
            b=c(bhmmab$samples[,,paste0("b[",k,"]")]),
            Smsy=c(bhmmab$samples[,,paste0("S_msy[",k,"]")]))),na.rm=T)
          
          sgen_hmmabhc[k] <- median(unlist(mapply(sGenCalc,a=c(bhmmabcaphi$samples[,,paste0("log_a[",k,"]")]),
            b=c(bhmmabcaphi$samples[,,paste0("b[",k,"]")]),
            Smsy=c(bhmmabcaphi$samples[,,paste0("S_msy[",k,"]")]))),na.rm=T)
        }

        sgen_hmma_regime <- sgen_hmma[bhmma$mcmcsummary[grep("zstar\\[",rownames(bhmma$mcmcsummary)),"50%"]]
        sgen_hmma_average<-matrix(bhmma$mcmcsummary[grep("^gamma\\[",rownames(bhmma$mcmcsummary)),"50%"],ncol=2,byrow=T)%*%sgen_hmma
        
        sgen_hmmb_regime<- sgen_hmmb[bhmmb$mcmcsummary[grep("zstar\\[",rownames(bhmmb$mcmcsummary)),"50%"]]
        sgen_hmmb_average<-matrix(bhmmb$mcmcsummary[grep("^gamma\\[",rownames(bhmmb$mcmcsummary)),"50%"],ncol=2,byrow=T)%*%sgen_hmmb
        
        sgen_hmmab_regime<-sgen_hmmab[bhmmab$mcmcsummary[grep("zstar\\[",rownames(bhmmab$mcmcsummary)),"50%"]]
        sgen_hmmab_average<-matrix(bhmmab$mcmcsummary[grep("^gamma\\[",rownames(bhmmab$mcmcsummary)),"50%"],ncol=2,byrow=T)%*%sgen_hmmab
        
        sgen_hmmabhc_regime<-sgen_hmmabhc[bhmmabcaphi$mcmcsummary[grep("zstar\\[",rownames(bhmmabcaphi$mcmcsummary)),"50%"]]
        sgen_hmmabhc_average<-matrix(bhmmabcaphi$mcmcsummary[grep("^gamma\\[",rownames(bhmmabcaphi$mcmcsummary)),"50%"],ncol=2,byrow=T)%*%sgen_hmmabhc 

        dfsgen <- data.frame(parameter="sgen",
          iteration=u,
          method=rep(c(rep("MLE",11),rep("MCMC",13)),each=nrow(df)),
          model=rep(c("simple","autocorr","rwa","rwb","rwab","hmma_regime","hmma_average",
            "hmmb_regime","hmmb_average","hmmab_regime","hmmab_average",
            "simple","autocorr","rwa","rwb","rwab","hmma_regime","hmma_average",
            "hmmb_regime","hmmb_average","hmmab_regime","hmmab_average","hmmabhc_regime","hmmabhc_average"),each=nrow(df)),
            by=rep(dat$year,24),
          sim=rep(unlist(mapply(sGenCalc,a=dat$alpha,Smsy=smsysim, b=dat$beta)),24),
          est=c(unlist(mapply(sGenCalc,a=dfa$est[dfa$model=="simple"&dfa$method=="MLE"],
               Smsy=dfsmsy$est[dfsmsy$model=="simple"&dfsmsy$method=="MLE"], 
               b=1/dfsmax$est[dfsmax$model=="simple"&dfsmax$method=="MLE"])),
              unlist(mapply(sGenCalc,a=dfa$est[dfa$model=="autocorr"&dfa$method=="MLE"],
                Smsy=dfsmsy$est[dfsmsy$model=="autocorr"&dfsmsy$method=="MLE"], 
                b=1/dfsmax$est[dfsmax$model=="autocorr"&dfsmax$method=="MLE"])),
              unlist(mapply(sGenCalc,a=dfa$est[dfa$model=="rwa"&dfa$method=="MLE"],
                Smsy=dfsmsy$est[dfsmsy$model=="rwa"&dfsmsy$method=="MLE"], 
                b=1/dfsmax$est[dfsmax$model=="rwa"&dfsmax$method=="MLE"])),
              unlist(mapply(sGenCalc,a=dfa$est[dfa$model=="rwb"&dfa$method=="MLE"],
                Smsy=dfsmsy$est[dfsmsy$model=="rwb"&dfsmsy$method=="MLE"], 
                b=1/dfsmax$est[dfsmax$model=="rwb"&dfsmax$method=="MLE"])),
              unlist(mapply(sGenCalc,a=dfa$est[dfa$model=="rwab"&dfa$method=="MLE"],
                Smsy=dfsmsy$est[dfsmsy$model=="rwab"&dfsmsy$method=="MLE"], 
                b=1/dfsmax$est[dfsmax$model=="rwab"&dfsmax$method=="MLE"])),
              unlist(mapply(sGenCalc,a=dfa$est[dfa$model=="hmma_regime"&dfa$method=="MLE"],
                Smsy=dfsmsy$est[dfsmsy$model=="hmma_regime"&dfsmsy$method=="MLE"], 
                b=1/dfsmax$est[dfsmax$model=="hmma_regime"&dfsmax$method=="MLE"])),
              unlist(mapply(sGenCalc,a=dfa$est[dfa$model=="hmma_average"&dfa$method=="MLE"],
                Smsy=dfsmsy$est[dfsmsy$model=="hmma_average"&dfsmsy$method=="MLE"], 
                b=1/dfsmax$est[dfsmax$model=="hmma_regime"&dfsmax$method=="MLE"])),
              unlist(mapply(sGenCalc,a=dfa$est[dfa$model=="hmmb_regime"&dfa$method=="MLE"],
                Smsy=dfsmsy$est[dfsmsy$model=="hmmb_regime"&dfsmsy$method=="MLE"],
                 b=1/dfsmax$est[dfsmax$model=="hmmb_regime"&dfsmax$method=="MLE"])),
              unlist(mapply(sGenCalc,a=dfa$est[dfa$model=="hmmb_regime"&dfa$method=="MLE"],
                Smsy=dfsmsy$est[dfsmsy$model=="hmmb_average"&dfsmsy$method=="MLE"], 
                b=1/dfsmax$est[dfsmax$model=="hmmb_average"&dfsmax$method=="MLE"])),
              unlist(mapply(sGenCalc,a=dfa$est[dfa$model=="hmmab_regime"&dfa$method=="MLE"],
                Smsy=dfsmsy$est[dfsmsy$model=="hmmab_regime"&dfsmsy$method=="MLE"], 
                b=1/dfsmax$est[dfsmax$model=="hmmab_regime"&dfsmax$method=="MLE"])),
              unlist(mapply(sGenCalc,a=dfa$est[dfa$model=="hmmab_average"&dfa$method=="MLE"],
                Smsy=dfsmsy$est[dfsmsy$model=="hmmab_average"&dfsmsy$method=="MLE"],
                 b=1/dfsmax$est[dfsmax$model=="hmmab_average"&dfsmax$method=="MLE"])),
              rep(sgen_b,nrow(df)),
              rep(sgen_bac,nrow(df)),
              sgen_tva,
              sgen_tvb,
              sgen_tvab,
              sgen_hmma_regime,
              sgen_hmma_average,
              sgen_hmmb_regime,
              sgen_hmmb_average,
              sgen_hmmab_regime,
              sgen_hmmab_average,
              sgen_hmmabhc_regime,
              sgen_hmmabhc_average
              ),
          convergence=c(rep(c(p$model$convergence + p$conv_problem,
            pac$model$convergence + pac$conv_problem,
            ptva$model$convergence + ptvab$conv_problem,
            ptvb$model$convergence + ptvab$conv_problem,
            ptvab$model$convergence + ptvab$conv_problem,
            phmma$model$convergence + phmma$conv_problem,
            phmma$model$convergence + phmma$conv_problem,
            phmmb$model$convergence + phmmb$conv_problem,
            phmmb$model$convergence + phmmb$conv_problem,
            phmm$model$convergence + phmm$conv_problem,
            phmm$model$convergence + phmm$conv_problem,
            sum(as.numeric(abs(b$mcmcsummary["log_a","Rhat"]-1)>.1),
                as.numeric(abs(b$mcmcsummary["b","Rhat"]-1)>.1)),
            sum(as.numeric(abs(bac$mcmcsummary["log_a","Rhat"]-1)>.1),
                as.numeric(abs(bac$mcmcsummary["b","Rhat"]-1)>.1))),
            each=nrow(df)),
            (as.numeric(abs(btva$mcmcsummary[grep("log_a\\[",rownames(btva$mcmcsummary)),"Rhat"]-1)>.1)+
             as.numeric(abs(btva$mcmcsummary[grep("S_msy\\[",rownames(btva$mcmcsummary)),"Rhat"]-1)>.1)+
                as.numeric(abs(btva$mcmcsummary["b","Rhat"]-1)>.1)),
            (as.numeric(abs(btvb$mcmcsummary["log_a","Rhat"]-1)>.1)+
             as.numeric(abs(btvb$mcmcsummary[grep("S_msy\\[",rownames(btvb$mcmcsummary)),"Rhat"]-1)>.1)+
                as.numeric(abs(btvb$mcmcsummary[grep("^b\\[",rownames(btvb$mcmcsummary)),"Rhat"]-1)>.1)),
            (as.numeric(abs(btvab$mcmcsummary[grep("log_a\\[",rownames(btvab$mcmcsummary)),"Rhat"]-1)>.1)+
             as.numeric(abs(btvab$mcmcsummary[grep("S_msy\\[",rownames(btvab$mcmcsummary)),"Rhat"]-1)>.1)+
                as.numeric(abs(btvab$mcmcsummary[grep("^b\\[",rownames(btvab$mcmcsummary)),"Rhat"]-1)>.1)),
            
             c(as.numeric(abs(bhmma$mcmcsummary[grep("S_msy\\[",rownames(bhmma$mcmcsummary)),
            "Rhat"]-1)>.1)[bhmma$mcmcsummary[grep("zstar\\[",rownames(bhmma$mcmcsummary)),"50%"]]+
             as.numeric(abs(bhmma$mcmcsummary[grep("log_a\\[",rownames(bhmma$mcmcsummary)),
            "Rhat"]-1)>.1)[bhmma$mcmcsummary[grep("zstar\\[",rownames(bhmma$mcmcsummary)),"50%"]]+
            as.numeric(abs(bhmma$mcmcsummary[grep("zstar\\[",rownames(bhmma$mcmcsummary)),"Rhat"]-1)>.1)+
            as.numeric(abs(bhmma$mcmcsummary["b","Rhat"]-1)>.1)),
          
            sum(as.numeric(abs(bhmma$mcmcsummary[grep("S_msy\\[",rownames(bhmma$mcmcsummary)),"Rhat"]-1)>.1))+
            sum(as.numeric(abs(bhmma$mcmcsummary[grep("log_a\\[",rownames(bhmma$mcmcsummary)),"Rhat"]-1)>.1))+
            as.numeric(abs(bhmma$mcmcsummary["b","Rhat"]-1)>.1)+
            sumpair(as.numeric(abs(bhmma$mcmcsummary[grep("^gamma\\[",rownames(bhmma$mcmcsummary)),"Rhat"]-1)>.1)),
          
            c(as.numeric(abs(bhmmb$mcmcsummary[grep("S_msy\\[",rownames(bhmmb$mcmcsummary)),
            "Rhat"]-1)>.1)[bhmmb$mcmcsummary[grep("zstar\\[",rownames(bhmmb$mcmcsummary)),"50%"]]+
             as.numeric(abs(bhmmb$mcmcsummary[grep("^b\\[",rownames(bhmmb$mcmcsummary)),
            "Rhat"]-1)>.1)[bhmmb$mcmcsummary[grep("zstar\\[",rownames(bhmmb$mcmcsummary)),"50%"]]+
            as.numeric(abs(bhmmb$mcmcsummary[grep("zstar\\[",rownames(bhmmb$mcmcsummary)),"Rhat"]-1)>.1)+
            as.numeric(abs(bhmmb$mcmcsummary["log_a","Rhat"]-1)>.1)),
          
            sum(as.numeric(abs(bhmmb$mcmcsummary[grep("S_msy\\[",rownames(bhmmb$mcmcsummary)),"Rhat"]-1)>.1))+
            sum(as.numeric(abs(bhmmb$mcmcsummary[grep("b\\[",rownames(bhmmb$mcmcsummary)),"Rhat"]-1)>.1))+
            as.numeric(abs(bhmmb$mcmcsummary["log_a","Rhat"]-1)>.1)+
            sumpair(as.numeric(abs(bhmmb$mcmcsummary[grep("^gamma\\[",rownames(bhmmb$mcmcsummary)),"Rhat"]-1)>.1)),
          
            c(as.numeric(abs(bhmmab$mcmcsummary[grep("S_msy\\[",rownames(bhmmab$mcmcsummary)),
            "Rhat"]-1)>.1)[bhmmab$mcmcsummary[grep("zstar\\[",rownames(bhmmab$mcmcsummary)),"50%"]]+
             as.numeric(abs(bhmmab$mcmcsummary[grep("log_a\\[",rownames(bhmmab$mcmcsummary)),
            "Rhat"]-1)>.1)[bhmmab$mcmcsummary[grep("zstar\\[",rownames(bhmmab$mcmcsummary)),"50%"]]+
            as.numeric(abs(bhmmab$mcmcsummary[grep("zstar\\[",rownames(bhmmab$mcmcsummary)),"Rhat"]-1)>.1)+
            as.numeric(abs(bhmmab$mcmcsummary[grep("^b\\[",rownames(bhmmab$mcmcsummary)),
            "Rhat"]-1)>.1)[bhmmab$mcmcsummary[grep("zstar\\[",rownames(bhmmab$mcmcsummary)),"50%"]]),
          
            sum(as.numeric(abs(bhmmab$mcmcsummary[grep("S_msy\\[",rownames(bhmmab$mcmcsummary)),"Rhat"]-1)>.1))+
            sum(as.numeric(abs(bhmmab$mcmcsummary[grep("log_a\\[",rownames(bhmmab$mcmcsummary)),"Rhat"]-1)>.1))+
            sum(as.numeric(abs(bhmmab$mcmcsummary[grep("^b\\[",rownames(bhmmab$mcmcsummary)),"Rhat"]-1)>.1))+
            sumpair(as.numeric(abs(bhmmab$mcmcsummary[grep("^gamma\\[",rownames(bhmmab$mcmcsummary)),"Rhat"]-1)>.1)),
        
            c(as.numeric(abs(bhmmabcaphi$mcmcsummary[grep("S_msy\\[",rownames(bhmmabcaphi$mcmcsummary)),
            "Rhat"]-1)>.1)[bhmmabcaphi$mcmcsummary[grep("zstar\\[",rownames(bhmmabcaphi$mcmcsummary)),"50%"]]+
             as.numeric(abs(bhmmabcaphi$mcmcsummary[grep("log_a\\[",rownames(bhmmabcaphi$mcmcsummary)),
            "Rhat"]-1)>.1)[bhmmabcaphi$mcmcsummary[grep("zstar\\[",rownames(bhmmabcaphi$mcmcsummary)),"50%"]]+
            as.numeric(abs(bhmmabcaphi$mcmcsummary[grep("zstar\\[",rownames(bhmmabcaphi$mcmcsummary)),"Rhat"]-1)>.1)+
            as.numeric(abs(bhmmabcaphi$mcmcsummary[grep("^b\\[",rownames(bhmmabcaphi$mcmcsummary)),
            "Rhat"]-1)>.1)[bhmmabcaphi$mcmcsummary[grep("zstar\\[",rownames(bhmmabcaphi$mcmcsummary)),"50%"]]),
         
            sum(as.numeric(abs(bhmmabcaphi$mcmcsummary[grep("S_msy\\[",rownames(bhmmabcaphi$mcmcsummary)),"Rhat"]-1)>.1))+
            sum(as.numeric(abs(bhmmabcaphi$mcmcsummary[grep("log_a\\[",rownames(bhmmabcaphi$mcmcsummary)),"Rhat"]-1)>.1))+
            sum(as.numeric(abs(bhmmabcaphi$mcmcsummary[grep("^b\\[",rownames(bhmmabcaphi$mcmcsummary)),"Rhat"]-1)>.1))+
            sumpair(as.numeric(abs(bhmmabcaphi$mcmcsummary[grep("^gamma\\[",rownames(bhmmabcaphi$mcmcsummary)),"Rhat"]-1)>.1))
            ))
  
        dfsgen$pbias<- ((dfsgen$est-dfsgen$sim)/dfsgen$sim)*100
        
        rmsesgen<-aggregate((dfsgen$est-dfsgen$sim)^2, list(model=dfsgen$model,method=dfsgen$method),
            function(x)sqrt(mean(x)))
        rmsesgen$iteration<-u
        rmsesgen$parameter<-"sgen"
        rmsesgen$convergence <- aggregate(dfsgen$convergence, list(model=dfsgen$model,method=dfsgen$method),
          function(x)sum(x,na.rm=T))$x
         
  
        dfumsy<- data.frame(parameter="umsy",
          iteration=u,
          method=rep(c(rep("MLE",10),rep("MCMC",12)),each=nrow(df)),
          model=rep(c("simple","autocorr","rwa","rwb","rwab","hmma_regime","hmma_average",
            "hmmb_regime","hmmab_regime","hmmab_average","simple","autocorr","rwa","rwb",
            "rwab","hmma_regime","hmma_average","hmmb_regime","hmmab_regime",
            "hmmab_average","hmmabhc_regime","hmmabhc_average"),each=nrow(df)),
             by=rep(dat$year,22),
          sim=rep(umsyCalc(dat$alpha),22),
          est=c(rep(p$umsy, nrow(df)),
                rep(pac$umsy, nrow(df)),
                 ptva$umsy,
                rep(ptvb$umsy, nrow(df)),
                ptvab$umsy,
                phmma$umsy[phmma$regime],
                phmma$umsy%*%phmma$probregime,
                rep(phmmb$umsy,nrow(df)),
                phmm$umsy[phmm$regime],
                phmm$umsy%*%phmm$probregime,
                rep(b$umsy,nrow(df)),
                rep(bac$umsy,nrow(df)),
                btva$umsy,
                rep(btvb$umsy,nrow(df)),
                btvab$umsy,
                bhmma$umsy_regime,
                bhmma$umsy_wgt,
                rep(bhmmb$umsy,nrow(df)),
                bhmmab$umsy_regime,
                bhmmab$umsy_wgt,
                bhmmabcaphi$umsy_regime,
                bhmmabcaphi$umsy_wgt
                ),
          convergence=c(rep(c(p$model$convergence + p$conv_problem,
            pac$model$convergence + pac$conv_problem,
            ptva$model$convergence+ ptva$conv_problem,
            ptvb$model$convergence+ ptvb$conv_problem,
            ptvab$model$convergence + ptvab$conv_problem,
            phmma$model$convergence + phmma$conv_problem,
            phmma$model$convergence + phmma$conv_problem,
            phmmb$model$convergence + phmmb$conv_problem,
            phmm$model$convergence + phmm$conv_problem,
            phmm$model$convergence + phmm$conv_problem,
            as.numeric(abs(b$mcmcsummary["U_msy","Rhat"]-1)>.1),
            as.numeric(abs(bac$mcmcsummary["U_msy","Rhat"]-1)>.1)
            ),each=nrow(df)),
            as.numeric(abs(btva$mcmcsummary[grep("U_msy\\[",rownames(btva$mcmcsummary)),"Rhat"]-1)>.1),
            rep(as.numeric(abs(btvb$mcmcsummary[grep("U_msy",rownames(btvb$mcmcsummary)),"Rhat"]-1)>.1),each=nrow(df)),
            as.numeric(abs(btvab$mcmcsummary[grep("U_msy\\[",rownames(btvab$mcmcsummary)),"Rhat"]-1)>.1),
            
            c(as.numeric(abs(bhmma$mcmcsummary[grep("U_msy\\[",rownames(bhmma$mcmcsummary)),
            "Rhat"]-1)>.1)[bhmma$mcmcsummary[grep("zstar\\[",rownames(bhmma$mcmcsummary)),"50%"]]+
            as.numeric(abs(bhmma$mcmcsummary[grep("zstar\\[",rownames(bhmma$mcmcsummary)),"Rhat"]-1)>.1)),
            
            sum(as.numeric(abs(bhmma$mcmcsummary[grep("U_msy\\[",rownames(bhmma$mcmcsummary)),"Rhat"]-1)>.1))+
            sumpair(as.numeric(abs(bhmma$mcmcsummary[grep("^gamma\\[",rownames(bhmma$mcmcsummary)),"Rhat"]-1)>.1)),
            
            rep(as.numeric(abs(bhmmb$mcmcsummary["U_msy","Rhat"]-1)>.1),nrow(df)),
              
            c(as.numeric(abs(bhmmab$mcmcsummary[grep("U_msy\\[",rownames(bhmmab$mcmcsummary)),
            "Rhat"]-1)>.1)[bhmmab$mcmcsummary[grep("zstar\\[",rownames(bhmmab$mcmcsummary)),"50%"]]+
            as.numeric(abs(bhmmab$mcmcsummary[grep("zstar\\[",rownames(bhmmab$mcmcsummary)),"Rhat"]-1)>.1)),
            
            sum(as.numeric(abs(bhmmab$mcmcsummary[grep("U_msy\\[",rownames(bhmmab$mcmcsummary)),"Rhat"]-1)>.1))+
            sumpair(as.numeric(abs(bhmmab$mcmcsummary[grep("^gamma\\[",rownames(bhmmab$mcmcsummary)),"Rhat"]-1)>.1)),
            
            c(as.numeric(abs(bhmmabcaphi$mcmcsummary[grep("U_msy\\[",rownames(bhmmabcaphi$mcmcsummary)),
            "Rhat"]-1)>.1)[bhmmabcaphi$mcmcsummary[grep("^zstar\\[",rownames(bhmmabcaphi$mcmcsummary)),"50%"]]+
            as.numeric(abs(bhmmabcaphi$mcmcsummary[grep("^zstar\\[",rownames(bhmmabcaphi$mcmcsummary)),"Rhat"]-1)>.1)),
            
            sum(as.numeric(abs(bhmmabcaphi$mcmcsummary[grep("U_msy\\[",rownames(bhmmabcaphi$mcmcsummary)),"Rhat"]-1)>.1))+
            sumpair(as.numeric(abs(bhmmabcaphi$mcmcsummary[grep("^gamma\\[",rownames(bhmmabcaphi$mcmcsummary)),"Rhat"]-1)>.1)))       
        )
        dfumsy$pbias<- ((dfumsy$est-dfumsy$sim)/dfumsy$sim)*100

        rmseumsy<-aggregate((dfumsy$est-dfumsy$sim)^2, list(model=dfumsy$model,method=dfumsy$method), function(x)sqrt(mean(x)))
        rmseumsy$iteration<-u
        rmseumsy$parameter<-"umsy"
        rmseumsy$convergence <- aggregate(dfumsy$convergence, list(model=dfumsy$model,method=dfumsy$method),
          function(x)sum(x,na.rm=T))$x
   
        rmse<-rbind(rmsea,rmsesmax,rmsesig,rmsesmsy,rmsesgen,rmseumsy)
        simest<-rbind(dfa,dfsmax,dfsig,dfsmsy,dfsgen,dfumsy)
        return(list(rmse,simest))
    })
    return(out)

  })
  

save(simresult, file=paste0("outs/simest/simest_",simPar$scenario[a],".Rdata"))

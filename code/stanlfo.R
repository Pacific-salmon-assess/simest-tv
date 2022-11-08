#=============================================================
#Testing routine 
#=============================================================


#
# run these if first time running script or if updates were implemented. 
#NOTE if samsim is already installed, this should override it, note we are isntalling
#from timevar branch
#install samsim 
remotes::install_github("Pacific-salmon-assess/samSim", ref="timevar", force=TRUE)

#install samest
remotes::install_git('https://github.com/Pacific-salmon-assess/samEst', force=TRUE)


#======================================================
library(samEst)
library(samSim)
library(ggplot2)
library(devtools)
library(gridExtra)
library(dplyr)
library(here)
#source("sgen_functions.R")
source("code/utils.R")
#TODO estimate only for 40 yrs of data.
options(mc.cores = parallel::detectCores())
#here::here()
## Load relevant input data
# Simulation run parameters describing different scenarios
simPar <- read.csv("data/HarCk/harcnkSimPars.csv")

## Store relevant object names to help run simulation 
scenNames <- unique(simPar$scenario)
dirNames <- sapply(scenNames, function(x) paste(x, unique(simPar$species),sep = "_"))

## First check to ensure that a single scenario can be run (only a small number
# of trials necessary)
plotscn <- TRUE
p <- list()
simData <- list()

#Define models (helps prevent crashing)
m1=sr_mod(type='static',ac = FALSE,par='n',loglik=T)
m2=sr_mod(type='static',ac = TRUE,par='n',loglik=T)
m3=sr_mod(type='rw',par='a',loglik=T)
m4=sr_mod(type='rw',par='b',loglik=T)
m5=sr_mod(type='rw',par='both',loglik=T)
m6=sr_mod(type='hmm',par='a',loglik=T)
m7=sr_mod(type='hmm',par='b',loglik=T)
m8_1=sr_mod(type='hmm',par='both',loglik=T,caphigh=F)
m8_2=sr_mod(type='hmm',par='both',loglik=T,caphigh=T)

#Create folders to store lfo outputs from simulations
if(!file.exists("outs/simest")){
  dir.create("outs/simest")
  dir.create("outs/simest/lfostan")
  dir.create("outs/simest/lfostan/decLinearCap")
  dir.create("outs/simest/lfostan/decLinearProd")
  dir.create("outs/simest/lfostan/decLinearProdshiftCap")
  dir.create("outs/simest/lfostan/randomwalkProd")
  dir.create("outs/simest/lfostan/regimeCap")
  dir.create("outs/simest/lfostan/regimeProd")
  dir.create("outs/simest/lfostan/regimeProdCap")
  dir.create("outs/simest/lfostan/shiftCap")
  dir.create("outs/simest/lfostan/shiftSigma")
  dir.create("outs/simest/lfostan/sineProd")
  dir.create("outs/simest/lfostan/stationary")
}


for(a in 1:nrow(simPar)){
  simData[[a]] <- readRDS(paste0("outs/SamSimOutputs/simData/", simPar$nameOM[a],"/",simPar$scenario[a],"/",
                                 paste(simPar$nameOM[a],"_", simPar$nameMP[a], "_", "CUsrDat.RData",sep="")))$srDatout
  
  for(u in unique(simData[[a]]$iteration)){
    dat<-simData[[a]][simData[[a]]$iteration==u,]
    dat<-dat[dat$year>(max(dat$year)-40),]
    
    dat <- dat[!is.na(dat$obsRecruits),]
    df <- data.frame(by=dat$year,
                     S=dat$obsSpawners,
                     R=dat$obsRecruits,
                     logRS=log(dat$obsRecruits/dat$obsSpawners))

    #Assess model fits for each model type
    #model 1 - static Ricker
    ll1<- stan_lfo_cv(mod=m1,type='static',df=df,L=12)
    #model 2 - static autocorrelated Ricker
    ll2<- stan_lfo_cv(mod=m2,type='tv',df=df,L=12)
    #model 3 - dynamic productivity Ricker
    ll3<- stan_lfo_cv(mod=m3,type='tv',df=df,L=12)
    #model 4 - dynamic capacity Ricker
    ll4<- stan_lfo_cv(mod=m4,type='tv',df=df,L=12)
    #model 5 - dynamic productivity & capacity Ricker
    ll5<- stan_lfo_cv(mod=m5,type='tv',df=df,L=12)
    #model 6 - productivity regime shift - 2 regimes
    ll6<- stan_lfo_cv(mod=m6,type='regime',df=df,L=12,K=2)
    #model 7 - capacity regime shift
    ll7<- stan_lfo_cv(mod=m7,type='regime',df=df,L=12,K=2)
    #model 8 - productivity and capacity regime shift
    ll8_1<- stan_lfo_cv(mod=m8_1,type='regime',df=df,L=12,K=2)
    #model 8 - productivity and capacity regime shift
    ll8_2<- stan_lfo_cv(mod=m8_2,type='regime',df=df,L=12,K=2)
    
    full_loglik[[i]]=rbind(ll1,ll2,ll3,ll4,ll5,ll6,ll7,ll8_1,ll8_2)
    rownames(full_loglik[[i]])=c('m1',paste('m2',c(1,3,5),sep="_"),paste('m3',c(1,3,5),sep="_"),paste('m4',c(1,3,5),sep="_"),paste('m5',c(1,3,5),sep="_"),paste('m6',c(1,3,5,'1w','3w','5w'),sep="_"),paste('m7',c(1,3,5,'1w','3w','5w'),sep="_"),paste('m8_1',c(1,3,5,'1w','3w','5w'),sep="_"),paste('m8_2',c(1,3,5,'1w','3w','5w'),sep="_"))
    wm2=which.max(apply(ll2,1,sum)) #select best likelihood from different timeframes (1-y back, 3-y back, 5-y back)
    wm3=which.max(apply(ll3,1,sum))#best fit for model 3
    wm4=which.max(apply(ll4,1,sum))#best fit for model 4
    wm5=which.max(apply(ll5,1,sum)) #best fit for model 5
    wm6=which.max(apply(ll6,1,sum)) #best fit for model 6
    wm7=which.max(apply(ll7,1,sum)) #best fit for model 7
    wm8=which.max(apply(rbind(ll8_1,ll8_2),1,sum)) #best fit for model 8 (two combinations of higher alpha, lower cap; higher alpha, higher cap)
    if(wm8<=6){
      mod_loglik[[i]]=rbind(ll1,ll2[wm2,],ll3[wm3,],ll4[wm4,],ll5[wm5,],ll6[wm6,],ll7[wm7,],ll8_1[wm8,])
      rownames(mod_loglik[[i]])[1:8]=paste('m',seq(1:8),sep='');rownames(mod_loglik[[i]])[8]='m8_1'
    }
    if(wm8>6){
      mod_loglik[[i]]=rbind(ll1,ll2[wm2,],ll3[wm3,],ll4[wm4,],ll5[wm5,],ll6[wm6,],ll7[wm7,],ll8_2[wm8-6,])
      rownames(mod_loglik[[i]])[1:8]=paste('m',seq(1:8),sep='');rownames(mod_loglik[[i]])[8]='m8_2'
    }
    write.csv(full_loglik[[i]],here('outs','simest','lfostan',simPar$scenario[a],paste(u,'full_loglik.csv',sep='')))
    write.csv(mod_loglik[[i]],here('outs','simest','lfostan',simPar$scenario[a],paste(u,'mod_loglik.csv',sep='')))
    
  }
}

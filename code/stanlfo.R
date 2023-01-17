
#NOTE if samsim is already installed, this should override it, note we are isntalling
#from timevar branch
#install samsim 
#remotes::install_github("Pacific-salmon-assess/samSim", ref="timevar", force=TRUE)

#install samest
#remotes::install_git('https://github.com/Pacific-salmon-assess/samEst', force=TRUE)
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

#LFO====
#Define models for LFO sets
m1=samEst::sr_mod(type='static',ac = FALSE,par='n',lfo=T)
m2=samEst::sr_mod(type='static',ac = TRUE,par='n',lfo=T)
m3=samEst::sr_mod(type='rw',par='a',lfo=T)
m4=samEst::sr_mod(type='rw',par='b',lfo=T)
m5=samEst::sr_mod(type='rw',par='both',lfo=T)
m6=samEst::sr_mod(type='hmm',par='a',lfo=T)
m7=samEst::sr_mod(type='hmm',par='b',lfo=T)
m8=samEst::sr_mod(type='hmm',par='both',lfo=T)

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

full_loglik=list()
mod_loglik=list()
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
    
    full_loglik[[u]]=rbind(ll1,ll2,ll3,ll4,ll5,ll6,ll7,ll8_1,ll8_2)
    rownames(full_loglik[[u]])=c('m1',paste('m2',c(1,3,5),sep="_"),paste('m3',c(1,3,5),sep="_"),paste('m4',c(1,3,5),sep="_"),paste('m5',c(1,3,5),sep="_"),paste('m6',c(1,3,5,'1w','3w','5w'),sep="_"),paste('m7',c(1,3,5,'1w','3w','5w'),sep="_"),paste('m8_1',c(1,3,5,'1w','3w','5w'),sep="_"),paste('m8_2',c(1,3,5,'1w','3w','5w'),sep="_"))
    wm2=which.max(apply(ll2,1,sum)) #select best likelihood from different timeframes (1-y back, 3-y back, 5-y back)
    wm3=which.max(apply(ll3,1,sum))#best fit for model 3
    wm4=which.max(apply(ll4,1,sum))#best fit for model 4
    wm5=which.max(apply(ll5,1,sum)) #best fit for model 5
    wm6=which.max(apply(ll6,1,sum)) #best fit for model 6
    wm7=which.max(apply(ll7,1,sum)) #best fit for model 7
    wm8=which.max(apply(rbind(ll8_1,ll8_2),1,sum)) #best fit for model 8 (two combinations of higher alpha, lower cap; higher alpha, higher cap)
    if(wm8<=6){
      mod_loglik[[u]]=rbind(ll1,ll2[wm2,],ll3[wm3,],ll4[wm4,],ll5[wm5,],ll6[wm6,],ll7[wm7,],ll8_1[wm8,])
      rownames(mod_loglik[[u]])[1:8]=paste('m',seq(1:8),sep='');rownames(mod_loglik[[u]])[8]='m8_1'
    }
    if(wm8>6){
      mod_loglik[[u]]=rbind(ll1,ll2[wm2,],ll3[wm3,],ll4[wm4,],ll5[wm5,],ll6[wm6,],ll7[wm7,],ll8_2[wm8-6,])
      rownames(mod_loglik[[u]])[1:8]=paste('m',seq(1:8),sep='');rownames(mod_loglik[[u]])[8]='m8_2'
    }
    write.csv(full_loglik[[u]],here('outs','simest','lfostan',simPar$scenario[a],paste(u,'full_loglik.csv',sep='')))
    write.csv(mod_loglik[[u]],here('outs','simest','lfostan',simPar$scenario[a],paste(u,'mod_loglik.csv',sep='')))
    
  }
}

#LOO====
#Define models for LOO sets
m1f=samEst::sr_mod(type='static',ac = FALSE,par='n',lfo=F)
m2f=samEst::sr_mod(type='static',ac = TRUE,par='n',lfo=F)
m3f=samEst::sr_mod(type='rw',par='a',lfo=F)
m4f=samEst::sr_mod(type='rw',par='b',lfo=F)
m5f=samEst::sr_mod(type='rw',par='both',lfo=F)
m6f=samEst::sr_mod(type='hmm',par='a',lfo=F)
m7f=samEst::sr_mod(type='hmm',par='b',lfo=F)
m8f=samEst::sr_mod(type='hmm',par='both',lfo=F)

loo_elpd=list()
mw_loo_pbma=list()
mw_loo_pbma_l30=list()
mw_loo_sw=list()
mw_loo_sw_l30=list()
for(a in 1:nrow(simPar)){
  simData[[a]] <- readRDS(paste0("outs/SamSimOutputs/simData/", simPar$nameOM[a],"/",simPar$scenario[a],"/",
                                 paste(simPar$nameOM[a],"_", simPar$nameMP[a], "_", "CUsrDat.RData",sep="")))$srDatout
  
  loodf<-matrix(NA,nrow=length(unique(simData[[a]]$iteration)),ncol=8,
                 dimnames = list(unique(simData$iteration),
                                 c("simple", "autocorr", 
                                   "rwa",
                                   "rwb",
                                   "rwab",
                                   "hmma",
                                   "hmmb",
                                   "hmmab")))
  loomw1df<-matrix(NA,nrow=length(unique(simData[[a]]$iteration)),ncol=8,
                   dimnames = list(unique(simData[[a]]$iteration),
                                   c("simple", "autocorr", 
                                     "rwa",
                                     "rwb",
                                     "rwab",
                                     "hmma",
                                     "hmmb",
                                     "hmmab")))
  
  loomw2df<-matrix(NA,nrow=length(unique(simData[[a]]$iteration)),ncol=8,
                   dimnames = list(unique(simData[[a]]$iteration),
                                   c("simple", "autocorr", 
                                     "rwa",
                                     "rwb",
                                     "rwab",
                                     "hmma",
                                     "hmmb",
                                     "hmmab")))
  
  loomw1dfl30<-matrix(NA,nrow=length(unique(simData[[a]]$iteration)),ncol=8,
                   dimnames = list(unique(simData[[a]]$iteration),
                                   c("simple", "autocorr", 
                                     "rwa",
                                     "rwb",
                                     "rwab",
                                     "hmma",
                                     "hmmb",
                                     "hmmab")))
  
  loomw2dfl30<-matrix(NA,nrow=length(unique(simData[[a]]$iteration)),ncol=8,
                      dimnames = list(unique(simData[[a]]$iteration),
                                      c("simple", "autocorr", 
                                        "rwa",
                                        "rwb",
                                        "rwab",
                                        "hmma",
                                        "hmmb",
                                        "hmmab")))
  
  for(u in unique(simData[[a]]$iteration)){
    dat<-simData[[a]][simData[[a]]$iteration==u,]
    dat<-dat[dat$year>(max(dat$year)-40),]
    
    dat <- dat[!is.na(dat$obsRecruits),]
    df <- data.frame(by=dat$year,
                     S=dat$obsSpawners,
                     R=dat$obsRecruits,
                     logRS=log(dat$obsRecruits/dat$obsSpawners))
    
    f1 = rstan::sampling(m1f, 
                         data = list(N=nrow(df),
                                     L=max(df$by)-min(df$by)+1,
                                     ii=df$by-min(df$by)+1,
                                     R_S=df$logRS,
                                     S=df$S),
                         control = list(adapt_delta = 0.99,max_treedepth=15), warmup = 200, chains = 6, iter = 700)
    #model 2 - static autocorrelated Ricker
    f2 = rstan::sampling(m2f, 
                         data = list(N=nrow(df),
                                     L=max(df$by)-min(df$by)+1,
                                     ii=df$by-min(df$by)+1,
                                     R_S=df$logRS,
                                     S=df$S),
                         control = list(adapt_delta = 0.99,max_treedepth=15), warmup = 200, chains = 6, iter = 700)
    
    #model 3 - dynamic productivity Ricker
    f3 = rstan::sampling(m3f, 
                         data = list(N=nrow(df),
                                     L=max(df$by)-min(df$by)+1,
                                     ii=df$by-min(df$by)+1,
                                     R_S=df$logRS,
                                     S=df$S),
                         control = list(adapt_delta = 0.99,max_treedepth=15), warmup = 200, chains = 6, iter = 700)
    
    #model 4 - dynamic capacity Ricker
    f4 = rstan::sampling(m4f, 
                         data = list(N=nrow(df),
                                     L=max(df$by)-min(df$by)+1,
                                     ii=df$by-min(df$by)+1,
                                     R_S=df$logRS,
                                     S=df$S),
                         control = list(adapt_delta = 0.99,max_treedepth=15), warmup = 200, chains = 6, iter = 700)
    
    #model 5 - dynamic productivity & capacity Ricker
    f5 = rstan::sampling(m5f, 
                         data = list(N=nrow(df),
                                     L=max(df$by)-min(df$by)+1,
                                     ii=df$by-min(df$by)+1,
                                     R_S=df$logRS,
                                     S=df$S),
                         control = list(adapt_delta = 0.99,max_treedepth=15), warmup = 200, chains = 6, iter = 700)
    
    #model 6 - productivity regime shift - 2 regimes
    f6 = rstan::sampling(m6f, 
                         data = list(N=nrow(df),
                                     R_S=df$logRS,
                                     S=df$S,
                                     K=2,
                                     alpha_dirichlet=rep(1,2)), #prior for state transition probabilities (this makes them equal)
                         control = list(adapt_delta = 0.99,max_treedepth=15), warmup = 200, chains = 6, iter = 700)
    
    #model 7 - capacity regime shift
    f7 = rstan::sampling(m7f, 
                         data = list(N=nrow(df),
                                     R_S=df$logRS,
                                     S=df$S,
                                     K=2,
                                     alpha_dirichlet=rep(1,2)), #prior for state transition probabilities (this makes them equal)
                         control = list(adapt_delta = 0.99,max_treedepth=15), warmup = 200, chains = 6, iter = 700)
    
    #model 8 - productivity and capacity regime shift
    f8= rstan::sampling(m8f, 
                           data = list(N=nrow(df),
                                       R_S=df$logRS,
                                       S=df$S,
                                       K=2,
                                       alpha_dirichlet=rep(1,2)), #prior for state transition probabilities (this makes them equal)
                           control = list(adapt_delta = 0.99,max_treedepth=15), warmup = 200, chains = 6, iter = 700)
    
    
    elpd.m1=loo::loo(f1,cores=4)  
    elpd.m2=loo::loo(f2,cores=4)
    elpd.m3=loo::loo(f3,cores=4)
    elpd.m4=loo::loo(f4,cores=4)
    elpd.m5=loo::loo(f5,cores=4)
    elpd.m6=loo::loo(f6,cores=4)
    elpd.m7=loo::loo(f7,cores=4)
    elpd.m8=loo::loo(f8,cores=4)
    
    loodf[u,]=c(elpd.m1$estimates[3,1],
                elpd.m2$estimates[3,1],
                elpd.m3$estimates[3,1],
                elpd.m4$estimates[3,1],
                elpd.m5$estimates[3,1],
                elpd.m6$estimates[3,1],
                elpd.m7$estimates[3,1],
                elpd.m8$estimates[3,1])
    
    lpd_point <- rbind(
        elpd.m1$pointwise[,"elpd_loo"], 
        elpd.m2$pointwise[,"elpd_loo"],
        elpd.m3$pointwise[,"elpd_loo"], 
        elpd.m4$pointwise[,"elpd_loo"],
        elpd.m5$pointwise[,"elpd_loo"], 
        elpd.m6$pointwise[,"elpd_loo"],
        elpd.m7$pointwise[,"elpd_loo"], 
        elpd.m8$pointwise[,"elpd_loo"]
      )
    
    loomw1df[u,]=samEst::model_weights(lpd_point,form='PBMA',type='full')
    
    lpd_point30 <- rbind(
        elpd.m1$pointwise[,"elpd_loo"][(round((2/3)*nrow(elpd.m1$pointwise))+1):nrow(elpd.m1$pointwise)], 
        elpd.m2$pointwise[,"elpd_loo"][(round((2/3)*nrow(elpd.m1$pointwise))+1):nrow(elpd.m1$pointwise)],
        elpd.m3$pointwise[,"elpd_loo"][(round((2/3)*nrow(elpd.m1$pointwise))+1):nrow(elpd.m1$pointwise)], 
        elpd.m4$pointwise[,"elpd_loo"][(round((2/3)*nrow(elpd.m1$pointwise))+1):nrow(elpd.m1$pointwise)],
        elpd.m5$pointwise[,"elpd_loo"][(round((2/3)*nrow(elpd.m1$pointwise))+1):nrow(elpd.m1$pointwise)], 
        elpd.m6$pointwise[,"elpd_loo"][(round((2/3)*nrow(elpd.m1$pointwise))+1):nrow(elpd.m1$pointwise)],
        elpd.m7$pointwise[,"elpd_loo"][(round((2/3)*nrow(elpd.m1$pointwise))+1):nrow(elpd.m1$pointwise)], 
        elpd.m8$pointwise[,"elpd_loo"][(round((2/3)*nrow(elpd.m1$pointwise))+1):nrow(elpd.m1$pointwise)]
      )
    loomw1dfl30[u,]=samEst::model_weights(lpd_point30,form='PBMA',type='full')
    
    
    loomw2df[u,]=loo::stacking_weights(cbind(
      elpd.m1$pointwise[,"elpd_loo"], 
      elpd.m2$pointwise[,"elpd_loo"],
      elpd.m3$pointwise[,"elpd_loo"], 
      elpd.m4$pointwise[,"elpd_loo"],
      elpd.m5$pointwise[,"elpd_loo"], 
      elpd.m6$pointwise[,"elpd_loo"],
      elpd.m7$pointwise[,"elpd_loo"], 
      elpd.m8$pointwise[,"elpd_loo"]
    ))
    
    loomw2dfl30[u,]=loo::stacking_weights(cbind(
      elpd.m1$pointwise[,"elpd_loo"][(round((2/3)*nrow(elpd.m1$pointwise))+1):nrow(elpd.m1$pointwise)], 
      elpd.m2$pointwise[,"elpd_loo"][(round((2/3)*nrow(elpd.m1$pointwise))+1):nrow(elpd.m1$pointwise)],
      elpd.m3$pointwise[,"elpd_loo"][(round((2/3)*nrow(elpd.m1$pointwise))+1):nrow(elpd.m1$pointwise)], 
      elpd.m4$pointwise[,"elpd_loo"][(round((2/3)*nrow(elpd.m1$pointwise))+1):nrow(elpd.m1$pointwise)],
      elpd.m5$pointwise[,"elpd_loo"][(round((2/3)*nrow(elpd.m1$pointwise))+1):nrow(elpd.m1$pointwise)], 
      elpd.m6$pointwise[,"elpd_loo"][(round((2/3)*nrow(elpd.m1$pointwise))+1):nrow(elpd.m1$pointwise)],
      elpd.m7$pointwise[,"elpd_loo"][(round((2/3)*nrow(elpd.m1$pointwise))+1):nrow(elpd.m1$pointwise)], 
      elpd.m8$pointwise[,"elpd_loo"][(round((2/3)*nrow(elpd.m1$pointwise))+1):nrow(elpd.m1$pointwise)]
    ))
  }
  
  loo_elpd[[a]] <- loodf
  mw_loo_pbma[[a]] <- loomw1df
  mw_loo_pbma_l30[[a]] <- loomw1dfl30
  mw_loo_sw[[a]] <- loomw2df
  mw_loo_sw_l30[[a]] <- loomw2dfl30
}
save(loo_elpd, mw_loo_pbma,mw_loo_pbma_l30,mw_loo_sw,mw_loo_sw_l30,file="outs/simest/simestloo_stan_prodcapscenarios.Rdata")


#Process LOO MS data####
load("outs/simest/simestloo_stan_prodcapscenarios.Rdata") 

#todo
#processing of lfo output
loochoicel<-list()
for(a in seq_len(nrow(simPar))){
  
  loomwdf<- mw_loo_pbma[[a]] 
  
  #fix a bug in the naming
 loo<-apply(loomwdf,1,which.max)
  
  
 loochoice<-data.frame(
    chsnmod=dimnames(loomwdf)[[2]][apply(loomwdf,1,which.max)])
  
  
 loochoice$chsnmod<-factor(loochoice$chsnmod, levels=c("simple", "autocorr", 
                                                        "rwa",
                                                        "rwb",
                                                        "rwab",
                                                        "hmma",
                                                        "hmmb",
                                                        "hmm"))
  
  loochoice$scenario<- simPar$nameOM[a]
  
  loochoice$method <-"LOO-CV"
  loochoicel[[a]]<-loochoice
}

dflfo <- do.call("rbind", lfochoicel)
dfaic <- do.call("rbind", aicchoicel)
dfbic <- do.call("rbind", bicchoicel)

df<-rbind(dflfo,dfaic,dfbic)
summary(df)

df$simulated <- dplyr::recode(df$scenario, 
                              "stationary"="simple",
                              "decLinearProd"="rwa",
                              "regimeProd"="hmma",
                              "sineProd"="rwa",
                              "regimeCap"="hmmb",
                              "decLinearCap"="rwb",
                              "sigmaShift"="simple",
                              "regimeProdCap"="hmm",
                              "shiftCap"="hmmb",
                              "decLinearProdshiftCap"="rwab")   

df$simulated_f<-factor(df$simulated, levels=c("simple",
                                              "autocorr",
                                              "rwa",
                                              "rwb",
                                              "rwab", 
                                              "hmma",
                                              "hmmb",
                                              "hmm"))

df$scenario_f <- factor(df$scenario,levels=c("stationary",
                                             "autocorr",
                                             "decLinearProd",        
                                             "regimeProd",
                                             "sineProd",
                                             "regimeCap",            
                                             "decLinearCap",
                                             "sigmaShift",            
                                             "regimeProdCap",
                                             "shiftCap",
                                             "decLinearProdshiftCap")
)



mytheme = list(
  theme_bw(16)+
    theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
          legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="bold", size=12),
          axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
)


modsel<-ggplot(df) +  
  geom_bar(aes(chsnmod,fill=method), 
           position = position_dodge(width = 0.9, preserve = "single"))+
  geom_rect(aes(xmin=as.numeric(simulated_f)-.5,
                xmax=as.numeric(simulated_f)+.5,
                ymin=-Inf,ymax=Inf),
            color="gray90",alpha=0.002)+
  facet_wrap(~scenario_f)+ 
  # theme_bw(14) +
  mytheme+
  theme(axis.text.x = element_text(angle = 90))+
  xlab("chosen estimation model")+
  scale_fill_viridis_d(begin=.3, end=.9) 
modsel


ggsave(
  filename = "outs/SamSimOutputs/plotcheck/model_selectionLFOallopt.pdf", 
  plot = modsel, 
  width = 14, height = 8
)


unique(dflfocm$chsnmod)

dt<-df |> count(chsnmod,simulated,method)
dt$simulated_f <- factor(dt$simulated, levels=
                           c("simple","autocorr", "rwa","rwb","rwab","hmma", "hmmb", "hmm"))
dt$estimated_f <- factor(dt$chsnmod, levels=
                           c("simple","autocorr", "rwa","rwb","rwab","hmma", "hmmb", "hmm"))
dt$nst<-0
summary(dt)

for(j in seq_along(unique(dt$simulated))){
  for(i in unique(dt$method)){
    dp<-dt[dt$simulated==unique(dt$simulated)[j]&dt$method==i,]
    dp$nst<-dp$n/sum(dp$n)*100
    dt$nst[dt$simulated==unique(dt$simulated)[j]&dt$method==i]<-dp$nst
  }
  
}


dt$diag<-dt$simulated_f==dt$estimated_f

confmat<-ggplot(data =  dt, mapping = aes(x = simulated_f, y = estimated_f)) +
  geom_tile(aes(fill = nst)) +
  geom_segment(data=transform(subset(dt, !!diag), 
                              simulated=as.numeric(simulated_f), 
                              estimated=as.numeric(estimated_f)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="white", size=2)+
  
  #scale_color_manual(guide = FALSE, values = c(`TRUE` = "black", `FALSE`=NA))+
  geom_text(aes(label =  nst), vjust = 1) +
  scale_fill_gradient(low="white", high="#009194") +
  #scale_colour_manual(values = c("white", "black"))+
  theme_bw() + theme(legend.position = "none")+
  facet_wrap(~method)+
  mytheme+
  theme(axis.text.x = element_text(angle = 90),legend.position="none")+
  ylab("estimated")+xlab("simulated")


ggsave(
  filename = "outs/SamSimOutputs/plotcheck/confmatMLE.png", 
  plot = confmat, 
  width = 12, height = 5
)


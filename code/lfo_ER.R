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
simPar <- read.csv("data/HarCkER/harcnkSimPars_ER.csv")
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
    
    if(sum(df$S==0)){
      next
    }
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

  
if(!file.exists("outs/simest")){
  dir.create("outs/simest") 
}

save(lfoTMB, lfomwTMB,aicTMB,bicTMB,file="outs/simest/lfo_erscenarios.Rdata")


load("outs/simest/lfo_erscenarios.Rdata") 

#todo
#processing of lfo output
lfochoicel<-list()
aicchoicel<-list()
bicchoicel<-list()


for(a in seq_len(nrow(simPar))){

  lfomwdf<- lfomwTMB[[a]] 
   
  lfo<-apply(lfomwdf,1,which.max)
  lfo[sapply( lfo,  function(x) length(x)==0L)] <- NA
  
  lfochoice<-data.frame(
    chsnmodorig=dimnames(lfomwdf)[[2]][unlist(lfo)])
  
  lfochoice$chsnmod <- dplyr::recode(lfochoice$chsnmodorig, 
      "rwa_last"="rwa","rwa_last3" ="rwa","rwa_last5"="rwa",
      "rwb_last"="rwb","rwb_last3"="rwb","rwb_last5"="rwb",
      "rwab_last"="rwab","rwab_last3"="rwab","rwab_last5"="rwab",
      "hmma_last_pick"="hmma","hmma_last3_pick"="hmma","hmma_last5_pick"="hmma",    
      "hmmb_last_pick"="hmmb","hmmb_last3_pick"="hmmb","hmmb_last5_pick"="hmmb",
      "hmm_last_pick"="hmm", "hmm_last3_pick"="hmm", "hmm_last5_pick"="hmm")    

  
  lfochoice$chsnmod<-factor(lfochoice$chsnmod, levels=c("simple", "autocorr", 
      "rwa",
      "rwb",
      "rwab",
      "hmma",
      "hmmb",
      "hmm"))

  lfochoice$scenario<- simPar$scenario[a]

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
  aictmp<-apply(aicdf,1,which.min)
  aictmp[sapply(aictmp,  function(x) length(x)==0L)] <- NA
  
  
  aicchoice<-data.frame(
    chsnmod=dimnames(aicdf)[[2]][unlist(aictmp)])

  aicchoice$chsnmod<-factor(aicchoice$chsnmod, levels=c("simple", "autocorr", 
      "rwa",
      "rwb",
      "rwab",
      "hmma",
      "hmmb",
      "hmm"
      ))
  aicchoice$scenario<- simPar$scenario[a]
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
  bictmp<-apply(bicdf,1,which.min)
  bictmp[sapply(bictmp,  function(x) length(x)==0L)] <- NA
  
  
  bicchoice<-data.frame(
    chsnmod=dimnames(bicdf)[[2]][unlist(bictmp)])
  bicchoice$chsnmodorig<-NA
  bicchoice$chsnmod<-factor(bicchoice$chsnmod, levels=c("simple", "autocorr", 
      "rwa",
      "rwb",
      "rwab",
      "hmma",
      "hmmb",
      "hmm"
      ))
  bicchoice$scenario<- simPar$scenario[a]
  bicchoice$method<- "BIC"

  bicchoicel[[a]]<-bicchoice

}

dflfo <- do.call("rbind", lfochoicel)
dfaic <- do.call("rbind", aicchoicel)
dfbic <- do.call("rbind", bicchoicel)

df<-rbind(dflfo,dfaic,dfbic)
unique(df$scenario)

df$simulated <- "simple"


df$simulated_f<-factor(df$simulated, levels=c("simple",
                                            "autocorr",
                                            "rwa",
                                            "rwb",
                                            "rwab", 
                                            "hmma",
                                            "hmmb",
                                            "hmm"))

df$scenario_f <- factor(df$scenario,levels=c("lowERLowError",
                                             "lowERHighError",   
                                             "highERLowError",  
                                             "highERHighError", 
                                             "ShiftERLowError", 
                                            "ShiftERHighError",
                                            "trendERLowError",  
                                            "trendERHighError")
)



mytheme = list(
    theme_bw(16)+
        theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
              legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="bold", size=12),
              axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
)

summary(df)
dff<-df[!is.na(df$chsnmod),]


modsel<-ggplot(dff) +  
 geom_bar(aes(chsnmod,fill=method), 
    position = position_dodge(width = 0.9, preserve = "single"))+
 geom_rect(aes(xmin=as.numeric(simulated_f)-.5,
                         xmax=as.numeric(simulated_f)+.5,
                         ymin=-Inf,ymax=Inf),
                    color="gray90",alpha=0.002)+
 facet_wrap(~scenario_f)+ 
 mytheme+
 theme(axis.text.x = element_text(angle = 90))+
 xlab("chosen estimation model")+
 scale_fill_viridis_d(begin=.3, end=.9) 
modsel


ggsave(
      filename = "outs/SamSimOutputs/plotcheck/model_selectionERscen.png", 
      plot = modsel, 
      width = 14, height = 8
    )




dt<-dff |> count(chsnmod,simulated,method,scenario)
dt$simulated_f <- factor(dt$simulated, levels=
  c("simple","autocorr", "rwa","rwb","rwab","hmma", "hmmb", "hmm"))
dt$estimated_f <- factor(dt$chsnmod, levels=
  c("simple","autocorr", "rwa","rwb","rwab","hmma", "hmmb", "hmm"))

dt$scenario_f <- factor(dt$scenario,levels=c("lowERLowError",
                                             "lowERHighError",   
                                             "highERLowError",  
                                             "highERHighError", 
                                             "ShiftERLowError", 
                                            "ShiftERHighError",
                                            "trendERLowError",  
                                            "trendERHighError")
)
 dt$nst<-0

for(j in seq_along(unique(dt$scenario))){
  for(i in unique(dt$method)){
    dp<-dt[dt$scenario==unique(dt$scenario)[j]&dt$method==i,]
    dp$nst<-dp$n/sum(dp$n)*100
    dt$nst[dt$scenario==unique(dt$scenario)[j]&dt$method==i]<-dp$nst
  }

}


confmat<-ggplot(data =  dt, mapping = aes(x = scenario_f, y = chsnmod)) +
  geom_tile(aes(fill = nst), colour = "white") +
  geom_text(aes(label = sprintf("%1.0f", nst)), vjust = 1) +
  scale_fill_gradient(low="white", high="#009194") +
  theme_bw() + theme(legend.position = "none")+
  facet_wrap(~method)+
  mytheme+
  theme(axis.text.x = element_text(angle = 90),legend.position="none")+
  ylab("estimated")+xlab("simulated")
confmat

ggsave(
      filename = "outs/SamSimOutputs/plotcheck/confmatMLE_er.png", 
      plot = confmat, 
      width = 10, height = 5
    )


#==========================
#todo




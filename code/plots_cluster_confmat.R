#routines to visualize smulation estimation Results
#Catarina Wor
#November 2022
#============================================

library(gridExtra)
library(ggplot2)
library(ggpubr)
source("code/utils.R")


mytheme = list(
    theme_classic(16)+
        theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
              legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="bold", size=12),
              axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
)



#========================================================================================================
#base case
#read in data
simPar <- read.csv("data/generic/SimPars.csv")

## Store relevant object names to help run simulation 


res<-readRDS(file = "outs/simest/generic/resbase.rds")
res<-res[res$convergence==0,]

aic=subset(res,parameter=='AIC'&method=='MLE')
bic=subset(res,parameter=='BIC'&method=='MLE')
lfo=subset(res,parameter=='LFO'&method=='MLE')

lfo<-lfo[lfo$model %in% c("simple","autocorr","rwa","rwb","rwab",
    "hmma","hmmb","hmmab"),]

scn<-factor(unique(aic$scenario), levels=c(
  "stationary", 
  "sigmaShift",
  "autocorr",
  "decLinearProd",  
  "sineProd",
  "regimeProd",
  "shiftProd", 
  "decLinearCap",
  "regimeCap", 
  "shiftCap",
  "decLinearProdshiftCap",
  "regimeProdCap" ) )


EM=c("stationary",
     "autocorr",
     "dynamic.a","regime.a",
     "dynamic.b","regime.b",
     "dynamic.ab","regime.ab")

##Confusion matrices
conf_matrix<-expand.grid(EM=EM,OM=scn)
conf_matrix$w_AIC=NA
conf_matrix$BIC=NA
conf_matrix$LFO=NA

cn1<-list()
cn2<-list()
cn3<-list()

aic_set<-list()
bic_set<-list()
lfo_set<-list()
o=0
for(a in seq_along(scn)){

  #AIC
  aica<-subset(aic,scenario==scn[a])
  aic_set[[a]]=tidyr::spread(aica[,-c(10,11,12,13)],key=model,value=mode)
  naicsim<-length(unique(aic_set[[a]]$iteration))
  aic_set[[a]]=aic_set[[a]][c(15,8,12,9,14,11,13,10)] #reorder estimation models
  #head(aic_set[[a]])
  sc1=apply(aic_set[[a]],1,which.min)
  
  cn1[[a]]=summary(factor(sc1,levels=seq(1:ncol(aic_set[[a]]))))/naicsim

  bica=subset(bic,scenario==scn[a])
  bic_set[[a]]=tidyr::spread(bica[,-c(10,11,12,13)],key=model,value=mode)
  nbicsim<-length(unique( bic_set[[a]]$iteration))
  bic_set[[a]]=bic_set[[a]][c(15,8,12,9,14,11,13,10)] #reorder estimation models

  sc2=apply(bic_set[[a]],1,which.min)
 
  cn2[[a]]=summary(factor(sc2,levels=seq(1:ncol(bic_set[[a]]))))/nbicsim

  lfoa=subset(lfo,scenario==scn[a])
  lfo_set[[a]]=tidyr::spread(lfoa[,-c(10,11,12,13)],key=model,value=mode)
  nlfosim <- length(unique( lfo_set[[a]]$iteration))
  lfo_set[[a]]=lfo_set[[a]][c(15,8,12,9,14,11,13,10)] #reorder estimation models
  #head(lfo_set[[a]])

  sc3=apply(lfo_set[[a]],1,which.max)
  cn3[[a]]=summary(factor(sc3,levels=seq(1:ncol(lfo_set[[a]]))))/ nlfosim


  myseq<-seq(from=o+1, length.out=length(EM))
  conf_matrix$w_AIC[myseq]<-cn1[[a]]
  conf_matrix$BIC[myseq]<-cn2[[a]]
  conf_matrix$LFO[myseq]<-cn3[[a]]
  o=max(myseq)

}

conf_matrix$eqem_om <- dplyr::recode(conf_matrix$OM, 
      "stationary"="stationary",
      "autocorr"="autocorr",
      "sigmaShift"="stationary", 
      "decLinearProd"="dynamic.a",
      "sineProd"="dynamic.a",
      "decLinearCap"="dynamic.b",
      "decLinearProdshiftCap"="dynamic.ab",
      "regimeProd"="regime.a",
      "shiftProd"="regime.a",
      "regimeCap"="regime.b",
      "shiftCap"="regime.b", 
      "regimeProdCap"="regime.ab",
      )   
conf_matrix$diag<-conf_matrix$eqem_om==conf_matrix$EM


conf_matrix$scentype<-dplyr::case_match(conf_matrix$OM, 
      "stationary"~"stationary",
      "autocorr"~"autocorrelation",
      "sigmaShift"~"sigma", 
      "decLinearProd"~"log(alpha)",
      "sineProd"~"log(alpha)",
      "regimeProd"~"log(alpha)",
      "shiftProd"~"log(alpha)",
      "decLinearCap"~"Smax",
      "regimeCap"~"Smax",
      "shiftCap"~"Smax", 
      "regimeProdCap"~"both",
      "decLinearProdshiftCap"~"both"
      )   

conf_matrix$scendesc<-dplyr::case_match(conf_matrix$OM, 
      "stationary"~"",
      "autocorr"~"",
      "sigmaShift"~"shift up", 
      "decLinearProd"~"linear decline",
      "sineProd"~"sine trend",
      "regimeProd"~"shift up + down",
      "shiftProd"~"shift down + up",
      "decLinearCap"~"linear decline",
      "regimeCap"~"shift up + down",
      "shiftCap"~"shift down", 
      "regimeProdCap"~"regime both",
      "decLinearProdshiftCap"~"trend & shift"
      )   
 regimeProd shiftProd decLinearCap regimeCap shiftCap decLinearProdshiftCap regimeProdCap

conf_matrix$fullname <- apply( conf_matrix[ , c("scendesc", "scentype") ] , 1 , paste , collapse = " " )


conf_matrix$fullname<-factor(conf_matrix$fullname, levels=c(
  " stationary",
  "shift up sigma",
  " autocorrelation",           
  "linear decline log(alpha)",
  "sine trend log(alpha)",  
  "shift up + down log(alpha)", 
  "shift down + up log(alpha)",       
  "linear decline Smax",
  "shift up + down Smax",              
  "shift down Smax",           
  "regime both both",           
 "trend & shift both"         
 ) )



paic=ggplot(data = conf_matrix, mapping = aes(x = fullname, y = EM)) +
  geom_tile(aes(fill = w_AIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_AIC,2)), vjust = 1,size=6) +
  ggtitle("AIC")+
  scale_x_discrete(labels = ~ stringr::str_wrap(as.character(.x), 15))+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray90", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1)) +
  xlab("Simulation Scenario")+ylab("Estimation Model")
paic


ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/confusion_matrices/base/AIC_MLE.png", 
  plot=paic, width = 10, height = 8)



pbic=ggplot(data =  conf_matrix, mapping = aes(x = fullname, y = EM)) +
  geom_tile(aes(fill = BIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(BIC,2)), vjust = 1,size=6) +
  ggtitle("BIC")+
  scale_x_discrete(labels = ~ stringr::str_wrap(as.character(.x), 15))+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray90", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
pbic

ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/confusion_matrices/base/BIC_MLE.png",
  plot=pbic, width = 10, height = 8)



plfo=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = LFO), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(LFO,2)), vjust = 1, size=6) +
  ggtitle("LFO MLE")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray90", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
plfo

ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/confusion_matrices/base/LFO_MLE.png", 
  plot=plfo, width = 10, height = 8)


#---------------------------------------------------------------------------------------------
#stan


resstan1<-readRDS(file = "outs/simest/generic/resstan1.rds")
resstan2<-readRDS(file = "outs/simest/generic/resstan2.rds")
resstan<-rbind(resstan1,resstan2)


aic=subset(resstan,parameter=='AIC')
bic=subset(resstan,parameter=='BIC')


scn<-factor(unique(aic$scenario), levels=c(
  "stationary", 
  "sigmaShift",
  "autocorr",
  "decLinearProd",  
  "sineProd",
  "decLinearCap",
  "decLinearProdshiftCap",      
 "regimeProd",
 "shiftProd", 
 "regimeCap", 
 "shiftCap",
 "regimeProdCap" ) )


EM=c("stationary",
     "autocorr",
     "dynamic.a","dynamic.b","dynamic.ab",
     "regime.a","regime.b","regime.ab")

##Confusion matrices
conf_matrix<-expand.grid(EM=EM,OM=scn)
conf_matrix$AIC=NA
conf_matrix$BIC=NA


cn1<-list()
cn2<-list()


aic_set<-list()
bic_set<-list()

o=0
for(a in seq_along(scn)){

  #AIC
  aica<-subset(aic,scenario==scn[a])
  aic_set[[a]]=tidyr::spread(aica[,-c(10,11,12,13)],key=model,value=median)
  aic_set[[a]]=aic_set[[a]][c(15,8,12,14,13,9,11,10)] #reorder estimation models

  sc1=apply(aic_set[[a]],1,which.max)
  cn1[[a]]=summary(factor(sc1,levels=seq(1:ncol(aic_set[[a]]))))/1000

  bica=subset(bic,scenario==scn[a])
  bic_set[[a]]=tidyr::spread(bica[,-c(10,11,12,13)],key=model,value=median)
  bic_set[[a]]=bic_set[[a]][c(15,8,12,14,13,9,11,10)] #reorder estimation models

  sc2=apply(bic_set[[a]],1,which.max)
  cn2[[a]]=summary(factor(sc2,levels=seq(1:ncol(bic_set[[a]]))))/1000



  myseq<-seq(from=o+1, length.out=length(EM))
  conf_matrix$w_AIC[myseq]<-cn1[[a]]
  conf_matrix$BIC[myseq]<-cn2[[a]]

  o=max(myseq)

}


conf_matrix$eqem_om <- dplyr::recode(conf_matrix$OM, 
      "stationary"="stationary",
      "autocorr"="autocorr",
      "sigmaShift"="stationary", 
      "decLinearProd"="dynamic.a",
      "sineProd"="dynamic.a",
      "decLinearCap"="dynamic.b",
      "decLinearProdshiftCap"="dynamic.ab",
      "regimeProd"="regime.a",
      "shiftProd"="regime.a",
      "regimeCap"="regime.b",
      "shiftCap"="regime.b", 
      "regimeProdCap"="regime.ab",
      )   
conf_matrix$diag<-conf_matrix$eqem_om==conf_matrix$EM

conf_matrix$EM





pmcaic=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = w_AIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_AIC,2)), vjust = 1,size=6) +
  ggtitle(" MCMC AIC")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray90", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1)) +
  xlab("Simulation Scenario")+ylab("Estimation Model")
pmcaic

ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/confusion_matrices/base/AIC_MCMC.png", 
  plot=pmcaic, width = 10, height = 8)


pmcbic=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = BIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(BIC,2)), vjust = 1,size=6) +
  ggtitle("MCMC BIC")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray90", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1)) +
  xlab("Simulation Scenario")+ylab("Estimation Model")
pmcbic
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/confusion_matrices/base/BIC_HMC.png",
 plot=pmcbic, width = 10, height = 8)


#===================================================================================
reslfo<-readRDS(file = "outs/simest/generic/resstanloo.rds")
head(reslfo)


scn<-factor(unique(reslfo$scenario), levels=c(
  "stationary", 
  "sigmaShift",
  "autocorr",
  "decLinearProd",  
  "sineProd",
  "regimeProd",
  "shiftProd", 
  "decLinearCap",
  "regimeCap", 
  "shiftCap",
  "decLinearProdshiftCap",
  "regimeProdCap"  ) )


EM=c("stationary",
     "autocorr",
     "dynamic.a","regime.a",
     "dynamic.b","regime.b",
     "dynamic.ab","regime.ab")
##Confusion matrices
conf_matrix_lfo<-expand.grid(EM=EM,OM=scn)
conf_matrix_lfo$LFOmcmc=NA


cn3<-list()
lfomcmc_set<-list()

o=0
for(a in seq_along(scn)){


  lfoa=subset(reslfo,scenario==scn[a])
  lfomcmc_set[[a]]=tidyr::spread(lfoa[,-9],key=model,value=est)
  nsim<-length(unique( lfomcmc_set[[a]]$iteration))
  lfomcmc_set[[a]]=lfomcmc_set[[a]][c(15,8,12,14,13,9,11,10)] #reorder estimation models

  sc3=apply(lfomcmc_set[[a]],1,which.max)
  cn3[[a]]=summary(factor(sc3,levels=seq(1:ncol(lfomcmc_set[[a]]))))/nsim

  myseq<-seq(from=o+1, length.out=length(EM))
  conf_matrix_lfo$LFOmcmc[myseq]<-cn3[[a]]
  o=max(myseq)

}


conf_matrix_lfo$eqem_om <- dplyr::recode(conf_matrix_lfo$OM, 
      "stationary"="stationary",
      "autocorr"="autocorr",
      "sigmaShift"="stationary", 
      "decLinearProd"="dynamic.a",
      "sineProd"="dynamic.a",
      "decLinearCap"="dynamic.b",
      "decLinearProdshiftCap"="dynamic.ab",
      "regimeProd"="regime.a",
      "shiftProd"="regime.a",
      "regimeCap"="regime.b",
      "shiftCap"="regime.b", 
      "regimeProdCap"="regime.ab",
      )   
conf_matrix_lfo$diag<-conf_matrix_lfo$eqem_om==conf_matrix_lfo$EM


conf_matrix_lfo$scentype<-dplyr::case_match(conf_matrix_lfo$OM, 
      "stationary"~"stationary",
      "autocorr"~"autocorrelation",
      "sigmaShift"~"sigma", 
      "decLinearProd"~"log(alpha)",
      "sineProd"~"log(alpha)",
      "regimeProd"~"log(alpha)",
      "shiftProd"~"log(alpha)",
      "decLinearCap"~"Smax",
      "regimeCap"~"Smax",
      "shiftCap"~"Smax", 
      "regimeProdCap"~"both",
      "decLinearProdshiftCap"~"both"
      )   

conf_matrix_lfo$scendesc<-dplyr::case_match(conf_matrix_lfo$OM, 
      "stationary"~"",
      "autocorr"~"",
      "sigmaShift"~"shift up", 
      "decLinearProd"~"linear decline",
      "sineProd"~"sine trend",
      "regimeProd"~"shift up + down",
      "shiftProd"~"shift down + up",
      "decLinearCap"~"linear decline",
      "regimeCap"~"shift up + down",
      "shiftCap"~"shift down", 
      "regimeProdCap"~"regime both",
      "decLinearProdshiftCap"~"trend & shift"
      )   
 
conf_matrix_lfo$fullname <- apply( conf_matrix_lfo[ , c("scendesc", "scentype") ] ,
                                   1 , paste , collapse = " " )


conf_matrix_lfo$fullname<-factor(conf_matrix_lfo$fullname, levels=c(
  " stationary",
  "shift up sigma",
  " autocorrelation",           
  "linear decline log(alpha)",
  "sine trend log(alpha)",  
  "shift up + down log(alpha)", 
  "shift down + up log(alpha)",       
  "linear decline Smax",
  "shift up + down Smax",              
  "shift down Smax",           
  "regime both both",           
 "trend & shift both"         
 ) )


pmclfo=ggplot(data =  conf_matrix_lfo, mapping = aes(x = fullname, y = EM)) +
  geom_tile(aes(fill = LFOmcmc), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(LFOmcmc,2)), vjust = 1, size=6) +
  ggtitle("LFO MCMC")+
  scale_fill_gradient(low="white", high="#009194")  +
  scale_x_discrete(labels = ~ stringr::str_wrap(as.character(.x), 15))+
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray70", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
pmclfo
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/confusion_matrices/base/LFO_MCMC.png", 
  plot=pmclfo, width = 10, height = 8)



multi.confusion <- ggarrange(paic,pbic, pmclfo,
                        nrow = 2, ncol = 2)
multi.confusion

#lfo with 3 last yrs average


reslfo3<-reslfo[reslfo$model %in% c("simple","autocorr","rwa_last3","rwb_last3","rwab_last3",
    "hmma_last3","hmmb_last3","hmmab_last3"),]


reslfo3$model<-dplyr::recode(reslfo3$model, 
      "simple"="simple",
      "autocorr"="autocorr",
      "rwa_last3"="rwa",
      "rwb_last3"="rwb",
      "rwab_last3"="rwab",
    "hmma_last3"="hmma",
    "hmmb_last3"="hmmb",
    "hmmab_last3"="hmmab"
      ) 

conf_matrix$LFO3mcmc<-NA

cn4<-list()
lfomcmc_set3<-list()

o=0
for(a in seq_along(scn)){

  lfoa=subset(reslfo3,scenario==scn[a])
  lfomcmc_set3[[a]]=tidyr::spread(lfoa[,-9],key=model,value=est)
  nsim<-length(unique( lfomcmc_set3[[a]]$iteration))
  lfomcmc_set3[[a]]=lfomcmc_set3[[a]][c(15,8,12,14,13,9,11,10)] #reorder estimation models

  sc4=apply(lfomcmc_set3[[a]],1,which.max)
  cn3[[a]]=summary(factor(sc4,levels=seq(1:ncol(lfomcmc_set3[[a]]))))/nsim


  myseq<-seq(from=o+1, length.out=length(EM))
  conf_matrix$LFO3mcmc[myseq]<-cn3[[a]]
  o=max(myseq)

}

conf_matrix$diag<-conf_matrix$eqem_om==conf_matrix$EM

pmclfo3=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = LFO3mcmc), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(LFO3mcmc,2)), vjust = 1, size=6) +
  ggtitle("LFO MCMC avg 3 yrs")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray70", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
pmclfo3
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/confusion_matrices/base/LFO_MCMC_3yrs.png",
 plot=pmclfo3, width = 10, height = 8)



#=============================================================================================
#avg5 years

reslfo5<-reslfo[reslfo$model %in% c("simple","autocorr","rwa_last5","rwb_last5","rwab_last5",
    "hmma_last5","hmmb_last5","hmmab_last5"),]

reslfo5$model<-dplyr::recode(reslfo5$model, 
      "simple"="simple",
      "autocorr"="autocorr",
      "rwa_last5"="rwa",
      "rwb_last5"="rwb",
      "rwab_last5"="rwab",
    "hmma_last5"="hmma",
    "hmmb_last5"="hmmb",
    "hmmab_last5"="hmmab"
      ) 


conf_matrix$LFO5mcmc<-NA


cn5<-list()
lfomcmc_set5<-list()


o=0
for(a in seq_along(scn)){

  lfoa=subset(reslfo5,scenario==scn[a])
  lfomcmc_set5[[a]]=tidyr::spread(lfoa[,-9],key=model,value=est)
  nsim<-length(unique( lfomcmc_set5[[a]]$iteration))
  lfomcmc_set5[[a]]=lfomcmc_set5[[a]][c(15,8,12,14,13,9,11,10)] #reorder estimation models

  sc5=apply(lfomcmc_set5[[a]],1,which.max)
  cn5[[a]]=summary(factor(sc5,levels=seq(1:ncol(lfomcmc_set5[[a]]))))/nsim


  myseq<-seq(from=o+1, length.out=length(EM))
  conf_matrix$LFO5mcmc[myseq]<-cn5[[a]]
  o=max(myseq)

}


conf_matrix$diag<-conf_matrix$eqem_om==conf_matrix$EM

pmclfo5=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = LFO5mcmc), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(LFO5mcmc,2)), vjust = 1, size=6) +
  ggtitle("LFO MCMC avg 5 yrs")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray70", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
pmclfo5
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/confusion_matrices/base/LFO_MCMC_5yrs.png", 
  plot=pmclfo5, width = 10, height = 8)



#========================================================================================================
#sensitivity a scenario
#read in data
simPar <- read.csv("data/sensitivity/SimPars.csv")

## Store relevant object names to help run simulation 
scenNames <- unique(simPar$scenario)

resa1<-readRDS(file = "outs/simest/sensitivity/res_a1.rds")
resa2<-readRDS(file = "outs/simest/sensitivity/res_a2.rds")


res_a<-rbind(resa1,resa2)


#res_a <- res_a[res_a$convergence==0,]

aic=subset(res_a, parameter=='AIC'&method=='MLE')
bic=subset(res_a, parameter=='BIC'&method=='MLE')
lfo=subset(res_a, parameter=='LFO'&method=='MLE')


lfo<-lfo[lfo$model %in% c("simple","autocorr","rwa","rwb","rwab","hmma","hmmb","hmmab"),]


aic$mode[aic$convergence>0]<-Inf
bic$mode[aic$convergence>0]<-Inf

scn<-factor(unique(aic$scenario), levels=c(
  "trendLinearProd1",
  "trendLinearProd1.3", 
  "trendLinearProd2", 
  "trendLinearProd5",
  "trendLinearProd7",
  "trendLinearProd10",
  "regimeProd1",
  "regimeProd1.3",
  "regimeProd2",      
  "regimeProd5",      
  "regimeProd7",
  "regimeProd10" ) )


EM=c("stationary",
     "autocorr",
     "dynamic.a","regime.a",
     "dynamic.b","regime.b",
     "dynamic.ab",
     "regime.ab")
##Confusion matrices
conf_matrix<-expand.grid(EM=EM,OM=scn)
conf_matrix$w_AIC=NA
conf_matrix$BIC=NA
conf_matrix$LFO=NA

cn1<-list()
cn2<-list()
cn3<-list()

#summarize model selection
aic_set=list()
bic_set=list()
lfo_set=list()


o=0
for(a in seq_along(scn)){

   #AIC
  aica<-subset(aic,scenario==scn[a])
  aic_set[[a]]=tidyr::spread(aica[,-c(10,11,12,13)],key=model,value=mode)
  aic_set[[a]]=aic_set[[a]][c(15,8,12,9,14,11,13,10)] #reorder estimation models

  sc1=apply(aic_set[[a]],1,which.min)
  cn1[[a]]=summary(factor(sc1,levels=seq(1:ncol(aic_set[[a]]))))/1000

  bica=subset(bic,scenario==scn[a])
  bic_set[[a]]=tidyr::spread(bica[,-c(10,11,12,13)],key=model,value=mode)
  bic_set[[a]]=bic_set[[a]][c(15,8,12,9,14,11,13,10)] #reorder estimation models
    head( bic_set[[a]])

  sc2=apply(bic_set[[a]],1,which.min)
  cn2[[a]]=summary(factor(sc2,levels=seq(1:ncol(bic_set[[a]]))))/1000

  lfoa=subset(lfo,scenario==scn[a])
  lfo_set[[a]]=tidyr::spread(lfoa[,-c(10,11,12,13)],key=model,value=mode)
  lfo_set[[a]]=lfo_set[[a]][c(15,8,12,9,14,11,13,10)] #reorder estimation models

  sc3=apply(lfo_set[[a]],1,which.max)
  cn3[[a]]=summary(factor(sc3,levels=seq(1:ncol(lfo_set[[a]]))))/1000


  myseq<-seq(from=o+1, length.out=length(EM))
  conf_matrix$w_AIC[myseq]<-cn1[[a]]
  conf_matrix$BIC[myseq]<-cn2[[a]]
  conf_matrix$LFO[myseq]<-cn3[[a]]
  o=max(myseq)
  
}


conf_matrix$eqem_om <- dplyr::recode(conf_matrix$OM, 
       "trendLinearProd1"="dynamic.a", 
       "trendLinearProd1.3"="dynamic.a", 
       "trendLinearProd2"="dynamic.a", 
       "trendLinearProd5"="dynamic.a",
       "trendLinearProd7"="dynamic.a",
       "trendLinearProd10"="dynamic.a",
       "regimeProd1"="regime.a",
       "regimeProd1.3"="regime.a",
       "regimeProd2"="regime.a",      
       "regimeProd5"="regime.a",      
       "regimeProd7"="regime.a",
       "regimeProd10"="regime.a"
      )   

conf_matrix$eqem_om<-factor(conf_matrix$eqem_om, 
    levels=c("stationary", 
        "autocorr", 
        "dynamic.a", "regime.a",
        "dynamic.b", "regime.b", 
        "dynamic.ab", "regime.ab" ))
conf_matrix$diag<-conf_matrix$eqem_om==conf_matrix$EM

conf_matrix$EM




pa_aic=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = w_AIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_AIC,2)), vjust = 1,size=6) +
  ggtitle("AIC sens a")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray90", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1)) +
  xlab("Simulation Scenario")+ylab("Estimation Model")
pa_aic
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/confusion_matrices/sens_a/AIC_MLE_sensa.png", 
  plot=pa_aic,width = 10, height = 8)



pa_bic=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = BIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(BIC,2)), vjust = 1,size=6) +
  ggtitle("BIC sens a")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray90", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
pa_bic
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/confusion_matrices/sens_a/BIC_MLE_sensa.png",
 plot=pa_bic, width = 10, height = 8)



p=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = LFO), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(LFO,2)), vjust = 1, size=6) +
  ggtitle("LFO")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray90", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
p
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/confusion_matrices/sens_a/LFO_MLE_sensa.png", 
  plot=p, width = 10, height = 8)


#LFO
reslfo<-readRDS(file = "outs/simest/sensitivity/resstanloo_a.rds")
head(reslfo)

#reslfo<-reslfo[reslfo$model %in% c("simple","autocorr","rwa_last5","rwb_last5","rwab_last5",
#    "hmma_last5","hmmb_last5","hmmab_last5"),]

#reslfo$model<-dplyr::recode(reslfo$model, 
#      "simple"="simple",
#      "autocorr"="autocorr",
#      "rwa_last5"="rwa",
#      "rwb_last5"="rwb",
#      "rwab_last5"="rwab",
#    "hmma_last5"="hmma",
#    "hmmb_last5"="hmmb",
#    "hmmab_last5"="hmmab"
#      ) 


#reslfo<-reslfo[reslfo$model %in% c("simple","autocorr","rwa_last3","rwb_last3","rwab_last3",
#    "hmma_last3","hmmb_last3","hmmab_last3"),]


#reslfo$model<-dplyr::recode(reslfo$model, 
#      "simple"="simple",
#      "autocorr"="autocorr",
#      "rwa_last3"="rwa",
#      "rwb_last3"="rwb",
#      "rwab_last3"="rwab",
#    "hmma_last3"="hmma",
#    "hmmb_last3"="hmmb",
#    "hmmab_last3"="hmmab"
#      ) 


##Confusion matrices
conf_matrix<-expand.grid(EM=EM,OM=scn)
conf_matrix$LFOmcmc=NA


cn3<-list()
lfomcmc_set_a<-list()


o=0
for(a in seq_along(scn)){


  lfoa=subset(reslfo,scenario==scn[a])
  lfomcmc_set_a[[a]]=tidyr::spread(lfoa[,-9],key=model,value=est)
  head(lfomcmc_set_a[[a]])
  nsim<-length(unique( lfomcmc_set_a[[a]]$iteration))
  lfomcmc_set_a[[a]]=lfomcmc_set_a[[a]][c(27,8,18,9,24,15,21,12)] #reorder estimation models

  sc3=apply(lfomcmc_set_a[[a]],1,which.max)
  cn3[[a]]=summary(factor(sc3,levels=seq(1:ncol(lfomcmc_set_a[[a]]))))/nsim


  myseq<-seq(from=o+1, length.out=length(EM))
  conf_matrix$LFOmcmc[myseq]<-cn3[[a]]
  o=max(myseq)

}


conf_matrix$eqem_om <- dplyr::recode(conf_matrix$OM, 
       "trendLinearProd1"="dynamic.a", 
       "trendLinearProd2"="dynamic.a", 
       "trendLinearProd5"="dynamic.a",
       "trendLinearProd7"="dynamic.a",
       "trendLinearProd10"="dynamic.a",
       "regimeProd1"="regime.a",
       "regimeProd2"="regime.a",      
       "regimeProd5"="regime.a",      
       "regimeProd7"="regime.a",
       "regimeProd10"="regime.a"
      ) 
conf_matrix$eqem_om<-factor(conf_matrix$eqem_om, 
    levels=c("stationary", 
        "autocorr", 
        "dynamic.a", "regime.a",
        "dynamic.b", "regime.b", 
        "dynamic.ab", "regime.ab" 
         
        
       ))
conf_matrix$diag<-conf_matrix$eqem_om==conf_matrix$EM







pmclfo_sensa=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = LFOmcmc), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(LFOmcmc,2)), vjust = 1, size=6) +
  ggtitle("sens a LFO MCMC")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray70", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
pmclfo_sensa
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/confusion_matrices/sens_a/LFO_MCMC_sensa.png", 
  plot=pmclfo_sensa, width = 10, height = 8)


#lfo with 3 last yrs average

reslfo3<-reslfo[reslfo$model %in% c("simple","autocorr","rwa_last3","rwb_last3","rwab_last3",
    "hmma_last3","hmmb_last3","hmmab_last3"),]


reslfo3$model<-dplyr::recode(reslfo3$model, 
      "simple"="simple",
      "autocorr"="autocorr",
      "rwa_last3"="rwa",
      "rwb_last3"="rwb",
      "rwab_last3"="rwab",
    "hmma_last3"="hmma",
    "hmmb_last3"="hmmb",
    "hmmab_last3"="hmmab"
      ) 

conf_matrix$LFO3mcmc<-NA


cn4<-list()
lfomcmc_set3_a<-list()


o=0
for(a in seq_along(scn)){
 
  lfoa=subset(reslfo3,scenario==scn[a])
  lfomcmc_set3_a[[a]]=tidyr::spread(lfoa[,-9],key=model,value=est)
  nsim<-length(unique( lfomcmc_set3_a[[a]]$iteration))
  lfomcmc_set3_a[[a]]=lfomcmc_set3_a[[a]][c(15,8,12,9,14,11,13,10)] #reorder estimation models

  #head(lfomcmc_set3_a[[a]])

  sc4=apply(lfomcmc_set3_a[[a]],1,which.max)
  cn4[[a]]=summary(factor(sc4,levels=seq(1:ncol(lfomcmc_set3_a[[a]]))))/nsim


  myseq<-seq(from=o+1, length.out=length(EM))
  conf_matrix$LFO3mcmc[myseq]<-cn4[[a]]
  o=max(myseq)

}


conf_matrix$diag<-conf_matrix$eqem_om==conf_matrix$EM





pmclfo3_sensa=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = LFO3mcmc), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(LFO3mcmc,2)), vjust = 1, size=6) +
  ggtitle("sens a LFO MCMC avg 3 yrs")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray70", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
pmclfo3_sensa
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/confusion_matrices/sens_a/LFO_MCMC_3yrs.png", plot=pmclfo3_sensa)




#=============================================================================================
#avg5 years

reslfo5<-reslfo[reslfo$model %in% c("simple","autocorr","rwa_last5","rwb_last5","rwab_last5",
    "hmma_last5","hmmb_last5","hmmab_last5"),]

reslfo5$model<-dplyr::recode(reslfo5$model, 
      "simple"="simple",
      "autocorr"="autocorr",
      "rwa_last5"="rwa",
      "rwb_last5"="rwb",
      "rwab_last5"="rwab",
    "hmma_last5"="hmma",
    "hmmb_last5"="hmmb",
    "hmmab_last5"="hmmab"
      ) 


conf_matrix$LFO5mcmc<-NA


cn5<-list()
lfomcmc_set5_a<-list()


o=0
for(a in seq_along(scn)){

  lfoa=subset(reslfo5,scenario==scn[a])
  lfomcmc_set5_a[[a]]=tidyr::spread(lfoa[,-9],key=model,value=est)
  nsim<-length(unique( lfomcmc_set5_a[[a]]$iteration))
  lfomcmc_set5_a[[a]]=lfomcmc_set5_a[[a]][c(15,8,12,14,13,9,11,10)] #reorder estimation models

  sc5=apply(lfomcmc_set5_a[[a]],1,which.max)
  cn5[[a]]=summary(factor(sc5,levels=seq(1:ncol(lfomcmc_set5_a[[a]]))))/nsim


  myseq<-seq(from=o+1, length.out=length(EM))
  conf_matrix$LFO5mcmc[myseq]<-cn5[[a]]
  o=max(myseq)

}


conf_matrix$diag<-conf_matrix$eqem_om==conf_matrix$EM





pmclfo5_sensa=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = LFO5mcmc), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(LFO5mcmc,2)), vjust = 1, size=6) +
  ggtitle("sens a LFO MCMC avg 5 yrs")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray70", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
pmclfo5_sensa
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/confusion_matrices/sens_a/LFO_MCMC_5yrs.png", plot=pmclfo5_sensa)






#========================================================================================================
#sensitivity a scenario -half Smax
#read in data
simPar <- read.csv("data/sensitivity_halfSmax/SimPars.csv")

## Store relevant object names to help run simulation 
scenNames <- unique(simPar$scenario)

res_asmax<-readRDS(file = "outs/simest/sensitivity_halfSmax/res_asmax.rds")

aic_asmax=subset(res_asmax, parameter=='AIC'&method=='MLE')
bic_asmax=subset(res_asmax, parameter=='BIC'&method=='MLE')
lfo_asmax=subset(res_asmax, parameter=='LFO'&method=='MLE')


lfo_asmax<-lfo_asmax[lfo_asmax$model %in% c("simple","autocorr","rwa","rwb","rwab","hmma","hmmb","hmmab"),]
lfo_asmax[is.na(lfo_asmax$est),]<--Inf

lfo_asmax<-lfo_asmax[lfo_asmax$model %in% c("simple","autocorr","rwa","rwb","rwab","hmma","hmmb","hmmab"),]

lfo_asmax[is.na(lfo_asmax$est),]<--Inf
aic_asmax$est[aic_asmax$convergence>0]<-Inf
bic_asmax$est[aic_asmax$convergence>0]<-Inf

scn<-factor(unique(aic_asmax$scenario), levels=c(
  "trendLinearProd1_halfSmax", 
  "trendLinearProd2_halfSmax", 
  "trendLinearProd5_halfSmax",
  "trendLinearProd7_halfSmax",
  "trendLinearProd10_halfSmax",
  "regimeProd1_halfSmax",
  "regimeProd2_halfSmax",      
  "regimeProd5_halfSmax",      
  "regimeProd7_halfSmax",
  "regimeProd10_halfSmax" ) )


EM=c("stationary",
     "autocorr",
     "dynamic.a","dynamic.b","dynamic.ab",
     "regime.a","regime.b","regime.ab")
##Confusion matrices
conf_matrix<-expand.grid(EM=EM,OM=scn)
conf_matrix$w_AIC=NA
conf_matrix$BIC=NA
conf_matrix$LFO=NA

cn1<-list()
cn2<-list()
cn3<-list()


#summarize model selection
aic_set=list()
bic_set=list()
lfo_set=list()


o=0
for(a in seq_along(scn)){

  #AIC
  aica<-subset(aic,scenario==scn[a])
  dim(aica)
  unique(aica$scenario)


   #AIC
  aica<-subset(aic_asmax,scenario==scn[a])
  aic_set[[a]]=tidyr::spread(aica[,-9],key=model,value=est)
  aic_set[[a]]=aic_set[[a]][c(15,8,12,14,13,9,11,10)] #reorder estimation models

  sc1=apply(aic_set[[a]],1,which.min)
  cn1[[a]]=summary(factor(sc1,levels=seq(1:ncol(aic_set[[a]]))))/1000

  bica=subset(bic_asmax,scenario==scn[a])
  bic_set[[a]]=tidyr::spread(bica[,-9],key=model,value=est)
  bic_set[[a]]=bic_set[[a]][c(15,8,12,14,13,9,11,10)] #reorder estimation models

  sc2=apply(bic_set[[a]],1,which.min)
  cn2[[a]]=summary(factor(sc2,levels=seq(1:ncol(bic_set[[a]]))))/1000

  lfoa=subset(lfo_asmax,scenario==scn[a])
  lfo_set[[a]]=tidyr::spread(lfoa[,-9],key=model,value=est)
  lfo_set[[a]]=lfo_set[[a]][c(15,8,12,14,13,9,11,10)] #reorder estimation models

  sc3=apply(lfo_set[[a]],1,which.max)
  cn3[[a]]=summary(factor(sc3,levels=seq(1:ncol(lfo_set[[a]]))))/1000


  myseq<-seq(from=o+1, length.out=length(EM))
  conf_matrix$w_AIC[myseq]<-cn1[[a]]
  conf_matrix$BIC[myseq]<-cn2[[a]]
  conf_matrix$LFO[myseq]<-cn3[[a]]
  o=max(myseq)

  
}


conf_matrix$eqem_om <- dplyr::recode(conf_matrix$OM, 
       "trendLinearProd1_halfSmax"="dynamic.a", 
       "trendLinearProd2_halfSmax"="dynamic.a", 
       "trendLinearProd5_halfSmax"="dynamic.a",
       "trendLinearProd7_halfSmax"="dynamic.a",
       "trendLinearProd10_halfSmax"="dynamic.a",
       "regimeProd1_halfSmax"="regime.a",
       "regimeProd2_halfSmax"="regime.a",      
       "regimeProd5_halfSmax"="regime.a",      
       "regimeProd7_halfSmax"="regime.a",
       "regimeProd10_halfSmax"="regime.a",
      )   

conf_matrix$eqem_om<-factor(conf_matrix$eqem_om, 
    levels=c("stationary", 
        "autocorr", 
        "dynamic.a", 
        "dynamic.b", 
        "dynamic.ab", 
        "regime.a", 
        "regime.b", 
        "regime.ab"))
conf_matrix$diag<-conf_matrix$eqem_om==conf_matrix$EM





p=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = w_AIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_AIC,2)), vjust = 1,size=6) +
  ggtitle("AIC sens a half smax")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray90", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1)) +
  xlab("Simulation Scenario")+ylab("Estimation Model")
p



p=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = BIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(BIC,2)), vjust = 1,size=6) +
  ggtitle("BIC sens a half smax")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray90", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
p



p=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = LFO), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(LFO,2)), vjust = 1, size=6) +
  ggtitle("LFO sens a half smax")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray90", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
p

reslfo<-readRDS(file = "outs/simest/generic/resstanloo.rds")
head(reslfo)

#reslfo<-reslfo[reslfo$model %in% c("simple","autocorr","rwa_last5","rwb_last5","rwab_last5",
#    "hmma_last5","hmmb_last5","hmmab_last5"),]

#reslfo$model<-dplyr::recode(reslfo$model, 
#      "simple"="simple",
#      "autocorr"="autocorr",
#      "rwa_last5"="rwa",
#      "rwb_last5"="rwb",
#      "rwab_last5"="rwab",
#    "hmma_last5"="hmma",
#    "hmmb_last5"="hmmb",
#    "hmmab_last5"="hmmab"
#      ) 


#reslfo<-reslfo[reslfo$model %in% c("simple","autocorr","rwa_last3","rwb_last3","rwab_last3",
#    "hmma_last3","hmmb_last3","hmmab_last3"),]


#reslfo$model<-dplyr::recode(reslfo$model, 
#      "simple"="simple",
#      "autocorr"="autocorr",
#      "rwa_last3"="rwa",
#      "rwb_last3"="rwb",
#      "rwab_last3"="rwab",
#    "hmma_last3"="hmma",
#    "hmmb_last3"="hmmb",
#    "hmmab_last3"="hmmab"
#      ) 

scn<-factor(unique(reslfo$scenario), levels=c(
  "stationary", 
  "sigmaShift",
  "autocorr",
  "decLinearProd",  
  "sineProd",
  "regimeProd",
  "shiftProd", 
  "decLinearCap",
  "regimeCap", 
  "shiftCap",
  "decLinearProdshiftCap",
  "regimeProdCap"  ) )


EM=c("stationary",
     "autocorr",
     "dynamic.a","regime.a",
     "dynamic.b","regime.b",
     "dynamic.ab","regime.ab")
##Confusion matrices
conf_matrix<-expand.grid(EM=EM,OM=scn)
conf_matrix$LFOmcmc=NA


cn3<-list()
lfomcmc_set<-list()


o=0
for(a in seq_along(scn)){


  lfoa=subset(reslfo,scenario==scn[a])
  lfomcmc_set[[a]]=tidyr::spread(lfoa[,-9],key=model,value=est)
  nsim<-length(unique( lfomcmc_set[[a]]$iteration))
  lfomcmc_set[[a]]=lfomcmc_set[[a]][c(15,8,12,14,13,9,11,10)] #reorder estimation models

  sc3=apply(lfomcmc_set[[a]],1,which.max)
  cn3[[a]]=summary(factor(sc3,levels=seq(1:ncol(lfomcmc_set[[a]]))))/nsim


  myseq<-seq(from=o+1, length.out=length(EM))
  conf_matrix$LFOmcmc[myseq]<-cn3[[a]]
  o=max(myseq)

}


conf_matrix$eqem_om <- dplyr::recode(conf_matrix$OM, 
      "stationary"="stationary",
      "autocorr"="autocorr",
      "sigmaShift"="stationary", 
      "decLinearProd"="dynamic.a",
      "sineProd"="dynamic.a",
      "decLinearCap"="dynamic.b",
      "decLinearProdshiftCap"="dynamic.ab",
      "regimeProd"="regime.a",
      "shiftProd"="regime.a",
      "regimeCap"="regime.b",
      "shiftCap"="regime.b", 
      "regimeProdCap"="regime.ab",
      )   
conf_matrix$diag<-conf_matrix$eqem_om==conf_matrix$EM

conf_matrix$EM






pmclfo=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = LFOmcmc), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(LFOmcmc,2)), vjust = 1, size=6) +
  ggtitle("LFO MCMC")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray70", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
pmclfo
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/confusion_matrices/base/LFO_MCMC.png", plot=pmclfo)


#lfo with 3 last yrs average


reslfo3<-reslfo[reslfo$model %in% c("simple","autocorr","rwa_last3","rwb_last3","rwab_last3",
    "hmma_last3","hmmb_last3","hmmab_last3"),]


reslfo3$model<-dplyr::recode(reslfo3$model, 
      "simple"="simple",
      "autocorr"="autocorr",
      "rwa_last3"="rwa",
      "rwb_last3"="rwb",
      "rwab_last3"="rwab",
    "hmma_last3"="hmma",
    "hmmb_last3"="hmmb",
    "hmmab_last3"="hmmab"
      ) 

conf_matrix$LFO3mcmc<-NA


cn4<-list()
lfomcmc_set3<-list()


o=0
for(a in seq_along(scn)){

  lfoa=subset(reslfo3,scenario==scn[a])
  lfomcmc_set3[[a]]=tidyr::spread(lfoa[,-9],key=model,value=est)
  nsim<-length(unique( lfomcmc_set3[[a]]$iteration))
  lfomcmc_set3[[a]]=lfomcmc_set3[[a]][c(15,8,12,14,13,9,11,10)] #reorder estimation models

  sc4=apply(lfomcmc_set3[[a]],1,which.max)
  cn3[[a]]=summary(factor(sc4,levels=seq(1:ncol(lfomcmc_set3[[a]]))))/nsim


  myseq<-seq(from=o+1, length.out=length(EM))
  conf_matrix$LFO3mcmc[myseq]<-cn3[[a]]
  o=max(myseq)

}


conf_matrix$diag<-conf_matrix$eqem_om==conf_matrix$EM





pmclfo3=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = LFO3mcmc), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(LFO3mcmc,2)), vjust = 1, size=6) +
  ggtitle("LFO MCMC avg 3 yrs")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray70", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
pmclfo3
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/confusion_matrices/base/LFO_MCMC_3yrs.png", plot=pmclfo3)




#=============================================================================================
#avg5 years

reslfo5<-reslfo[reslfo$model %in% c("simple","autocorr","rwa_last5","rwb_last5","rwab_last5",
    "hmma_last5","hmmb_last5","hmmab_last5"),]

reslfo5$model<-dplyr::recode(reslfo5$model, 
      "simple"="simple",
      "autocorr"="autocorr",
      "rwa_last5"="rwa",
      "rwb_last5"="rwb",
      "rwab_last5"="rwab",
    "hmma_last5"="hmma",
    "hmmb_last5"="hmmb",
    "hmmab_last5"="hmmab"
      ) 


conf_matrix$LFO5mcmc<-NA


cn5<-list()
lfomcmc_set5<-list()


o=0
for(a in seq_along(scn)){

  lfoa=subset(reslfo5,scenario==scn[a])
  lfomcmc_set5[[a]]=tidyr::spread(lfoa[,-9],key=model,value=est)
  nsim<-length(unique( lfomcmc_set5[[a]]$iteration))
  lfomcmc_set5[[a]]=lfomcmc_set5[[a]][c(15,8,12,14,13,9,11,10)] #reorder estimation models

  sc5=apply(lfomcmc_set5[[a]],1,which.max)
  cn5[[a]]=summary(factor(sc5,levels=seq(1:ncol(lfomcmc_set5[[a]]))))/nsim


  myseq<-seq(from=o+1, length.out=length(EM))
  conf_matrix$LFO5mcmc[myseq]<-cn5[[a]]
  o=max(myseq)

}


conf_matrix$diag<-conf_matrix$eqem_om==conf_matrix$EM





pmclfo5=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = LFO5mcmc), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(LFO5mcmc,2)), vjust = 1, size=6) +
  ggtitle("LFO MCMC avg 5 yrs")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray70", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
pmclfo5
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/confusion_matrices/base/LFO_MCMC_5yrs.png", plot=pmclfo5)





#========================================================================================================
#sensitivity smax scenario 
#read in data
simPar <- read.csv("data/Smax_sensitivity/SimPars.csv")

## Store relevant object names to help run simulation 
scenNames <- unique(simPar$scenario)

res_smax<-readRDS(file = "outs/simest/Smax_sensitivity/res_smax.rds")
#res_smax95<-readRDS(file = "outs/simest/res_gamma_alpha/res_smax_95.rds")


ressmax<-rbind(res_smax)

#res_a <- res_a[res_a$convergence==0,]

aic_smax=subset(ressmax, parameter=='AIC'&method=='MLE')
bic_smax=subset(ressmax, parameter=='BIC'&method=='MLE')
lfo_smax=subset(ressmax, parameter=='LFO'&method=='MLE')


lfo_smax<-lfo_smax[lfo_smax$model %in% c("simple","autocorr","rwa","rwb","rwab","hmma","hmmb","hmmab"),]
lfo_smax[is.na(lfo_smax$est),]<--Inf

lfo_smax<-lfo_smax[lfo_smax$model %in% c("simple","autocorr","rwa","rwb","rwab","hmma","hmmb","hmmab"),]

lfo_smax[is.na(lfo_smax$est),]<--Inf
aic_smax$est[aic_smax$convergence>0]<-Inf
bic_smax$est[aic_smax$convergence>0]<-Inf

scn<-factor(unique(aic_smax$scenario), levels=c(
    "trendLinearSmax025",
    "trendLinearSmax050", 
    "trendLinearSmax150", 
    "trendLinearSmax200", 
    "trendLinearSmax300", 
    "regimeSmax025",     
    "regimeSmax050",      
    "regimeSmax150",      
    "regimeSmax200",      
    "regimeSmax300"  ) )


EM=c("stationary",
     "autocorr",
     "dynamic.a","dynamic.b","dynamic.ab",
     "regime.a","regime.b","regime.ab")
##Confusion matrices
conf_matrix<-expand.grid(EM=EM,OM=scn)
conf_matrix$w_AIC=NA
conf_matrix$BIC=NA
conf_matrix$LFO=NA

cn1<-list()
cn2<-list()
cn3<-list()


#summarize model selection
aic_set=list()
bic_set=list()
lfo_set=list()


o=0
for(a in seq_along(scn)){

  #AIC
  aica<-subset(aic,scenario==scn[a])
  dim(aica)
  unique(aica$scenario)


   #AIC
  aica<-subset(aic_smax,scenario==scn[a])
  aic_set[[a]]=tidyr::spread(aica[,-9],key=model,value=est)
  aic_set[[a]]=aic_set[[a]][c(15,8,12,14,13,9,11,10)] #reorder estimation models

  sc1=apply(aic_set[[a]],1,which.min)
  cn1[[a]]=summary(factor(sc1,levels=seq(1:ncol(aic_set[[a]]))))/1000

  bica=subset(bic_smax,scenario==scn[a])
  bic_set[[a]]=tidyr::spread(bica[,-9],key=model,value=est)
  bic_set[[a]]=bic_set[[a]][c(15,8,12,14,13,9,11,10)] #reorder estimation models

  sc2=apply(bic_set[[a]],1,which.min)
  cn2[[a]]=summary(factor(sc2,levels=seq(1:ncol(bic_set[[a]]))))/1000

  lfoa=subset(lfo_smax,scenario==scn[a])
  lfo_set[[a]]=tidyr::spread(lfoa[,-9],key=model,value=est)
  lfo_set[[a]]=lfo_set[[a]][c(15,8,12,14,13,9,11,10)] #reorder estimation models

  sc3=apply(lfo_set[[a]],1,which.max)
  cn3[[a]]=summary(factor(sc3,levels=seq(1:ncol(lfo_set[[a]]))))/1000


  myseq<-seq(from=o+1, length.out=length(EM))
  conf_matrix$w_AIC[myseq]<-cn1[[a]]
  conf_matrix$BIC[myseq]<-cn2[[a]]
  conf_matrix$LFO[myseq]<-cn3[[a]]
  o=max(myseq)

  
}


conf_matrix$eqem_om <- dplyr::recode(conf_matrix$OM, 
    "trendLinearSmax025"="dynamic.b",
    "trendLinearSmax050"="dynamic.b", 
    "trendLinearSmax150"="dynamic.b", 
    "trendLinearSmax200"="dynamic.b", 
    "trendLinearSmax300"="dynamic.b", 
    "regimeSmax025"="regime.b",     
    "regimeSmax050"="regime.b",      
    "regimeSmax150"="regime.b",      
    "regimeSmax200"="regime.b",      
    "regimeSmax300"="regime.b"
      )   

conf_matrix$eqem_om<-factor(conf_matrix$eqem_om, 
    levels=c("stationary", 
        "autocorr", 
        "dynamic.a", 
        "dynamic.b", 
        "dynamic.ab", 
        "regime.a", 
        "regime.b", 
        "regime.ab"))
conf_matrix$diag<-conf_matrix$eqem_om==conf_matrix$EM





p=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = w_AIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_AIC,2)), vjust = 1,size=6) +
  ggtitle("AIC sens smax")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray90", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1)) +
  xlab("Simulation Scenario")+ylab("Estimation Model")
p



p=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = BIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(BIC,2)), vjust = 1,size=6) +
  ggtitle("BIC sens smax")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray90", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
p



p=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = LFO), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(LFO,2)), vjust = 1, size=6) +
  ggtitle("LFO sens smax")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray90", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
p


reslfo<-readRDS(file = "outs/simest/generic/resstanloo.rds")
head(reslfo)

#reslfo<-reslfo[reslfo$model %in% c("simple","autocorr","rwa_last5","rwb_last5","rwab_last5",
#    "hmma_last5","hmmb_last5","hmmab_last5"),]

#reslfo$model<-dplyr::recode(reslfo$model, 
#      "simple"="simple",
#      "autocorr"="autocorr",
#      "rwa_last5"="rwa",
#      "rwb_last5"="rwb",
#      "rwab_last5"="rwab",
#    "hmma_last5"="hmma",
#    "hmmb_last5"="hmmb",
#    "hmmab_last5"="hmmab"
#      ) 


#reslfo<-reslfo[reslfo$model %in% c("simple","autocorr","rwa_last3","rwb_last3","rwab_last3",
#    "hmma_last3","hmmb_last3","hmmab_last3"),]


#reslfo$model<-dplyr::recode(reslfo$model, 
#      "simple"="simple",
#      "autocorr"="autocorr",
#      "rwa_last3"="rwa",
#      "rwb_last3"="rwb",
#      "rwab_last3"="rwab",
#    "hmma_last3"="hmma",
#    "hmmb_last3"="hmmb",
#    "hmmab_last3"="hmmab"
#      ) 

scn<-factor(unique(reslfo$scenario), levels=c(
  "stationary", 
  "sigmaShift",
  "autocorr",
  "decLinearProd",  
  "sineProd",
  "regimeProd",
  "shiftProd", 
  "decLinearCap",
  "regimeCap", 
  "shiftCap",
  "decLinearProdshiftCap",
  "regimeProdCap"  ) )


EM=c("stationary",
     "autocorr",
     "dynamic.a","regime.a",
     "dynamic.b","regime.b",
     "dynamic.ab","regime.ab")
##Confusion matrices
conf_matrix<-expand.grid(EM=EM,OM=scn)
conf_matrix$LFOmcmc=NA


cn3<-list()
lfomcmc_set<-list()


o=0
for(a in seq_along(scn)){


  lfoa=subset(reslfo,scenario==scn[a])
  lfomcmc_set[[a]]=tidyr::spread(lfoa[,-9],key=model,value=est)
  nsim<-length(unique( lfomcmc_set[[a]]$iteration))
  lfomcmc_set[[a]]=lfomcmc_set[[a]][c(15,8,12,14,13,9,11,10)] #reorder estimation models

  sc3=apply(lfomcmc_set[[a]],1,which.max)
  cn3[[a]]=summary(factor(sc3,levels=seq(1:ncol(lfomcmc_set[[a]]))))/nsim


  myseq<-seq(from=o+1, length.out=length(EM))
  conf_matrix$LFOmcmc[myseq]<-cn3[[a]]
  o=max(myseq)

}


conf_matrix$eqem_om <- dplyr::recode(conf_matrix$OM, 
      "stationary"="stationary",
      "autocorr"="autocorr",
      "sigmaShift"="stationary", 
      "decLinearProd"="dynamic.a",
      "sineProd"="dynamic.a",
      "decLinearCap"="dynamic.b",
      "decLinearProdshiftCap"="dynamic.ab",
      "regimeProd"="regime.a",
      "shiftProd"="regime.a",
      "regimeCap"="regime.b",
      "shiftCap"="regime.b", 
      "regimeProdCap"="regime.ab",
      )   
conf_matrix$diag<-conf_matrix$eqem_om==conf_matrix$EM

conf_matrix$EM






pmclfo=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = LFOmcmc), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(LFOmcmc,2)), vjust = 1, size=6) +
  ggtitle("LFO MCMC")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray70", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
pmclfo
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/confusion_matrices/base/LFO_MCMC.png", plot=pmclfo)


#lfo with 3 last yrs average


reslfo3<-reslfo[reslfo$model %in% c("simple","autocorr","rwa_last3","rwb_last3","rwab_last3",
    "hmma_last3","hmmb_last3","hmmab_last3"),]


reslfo3$model<-dplyr::recode(reslfo3$model, 
      "simple"="simple",
      "autocorr"="autocorr",
      "rwa_last3"="rwa",
      "rwb_last3"="rwb",
      "rwab_last3"="rwab",
    "hmma_last3"="hmma",
    "hmmb_last3"="hmmb",
    "hmmab_last3"="hmmab"
      ) 

conf_matrix$LFO3mcmc<-NA


cn4<-list()
lfomcmc_set3<-list()


o=0
for(a in seq_along(scn)){

  lfoa=subset(reslfo3,scenario==scn[a])
  lfomcmc_set3[[a]]=tidyr::spread(lfoa[,-9],key=model,value=est)
  nsim<-length(unique( lfomcmc_set3[[a]]$iteration))
  lfomcmc_set3[[a]]=lfomcmc_set3[[a]][c(15,8,12,14,13,9,11,10)] #reorder estimation models

  sc4=apply(lfomcmc_set3[[a]],1,which.max)
  cn3[[a]]=summary(factor(sc4,levels=seq(1:ncol(lfomcmc_set3[[a]]))))/nsim


  myseq<-seq(from=o+1, length.out=length(EM))
  conf_matrix$LFO3mcmc[myseq]<-cn3[[a]]
  o=max(myseq)

}


conf_matrix$diag<-conf_matrix$eqem_om==conf_matrix$EM





pmclfo3=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = LFO3mcmc), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(LFO3mcmc,2)), vjust = 1, size=6) +
  ggtitle("LFO MCMC avg 3 yrs")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray70", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
pmclfo3
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/confusion_matrices/base/LFO_MCMC_3yrs.png", plot=pmclfo3)




#=============================================================================================
#avg5 years

reslfo5<-reslfo[reslfo$model %in% c("simple","autocorr","rwa_last5","rwb_last5","rwab_last5",
    "hmma_last5","hmmb_last5","hmmab_last5"),]

reslfo5$model<-dplyr::recode(reslfo5$model, 
      "simple"="simple",
      "autocorr"="autocorr",
      "rwa_last5"="rwa",
      "rwb_last5"="rwb",
      "rwab_last5"="rwab",
    "hmma_last5"="hmma",
    "hmmb_last5"="hmmb",
    "hmmab_last5"="hmmab"
      ) 


conf_matrix$LFO5mcmc<-NA


cn5<-list()
lfomcmc_set5<-list()


o=0
for(a in seq_along(scn)){

  lfoa=subset(reslfo5,scenario==scn[a])
  lfomcmc_set5[[a]]=tidyr::spread(lfoa[,-9],key=model,value=est)
  nsim<-length(unique( lfomcmc_set5[[a]]$iteration))
  lfomcmc_set5[[a]]=lfomcmc_set5[[a]][c(15,8,12,14,13,9,11,10)] #reorder estimation models

  sc5=apply(lfomcmc_set5[[a]],1,which.max)
  cn5[[a]]=summary(factor(sc5,levels=seq(1:ncol(lfomcmc_set5[[a]]))))/nsim


  myseq<-seq(from=o+1, length.out=length(EM))
  conf_matrix$LFO5mcmc[myseq]<-cn5[[a]]
  o=max(myseq)

}


conf_matrix$diag<-conf_matrix$eqem_om==conf_matrix$EM





pmclfo5=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = LFO5mcmc), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(LFO5mcmc,2)), vjust = 1, size=6) +
  ggtitle("LFO MCMC avg 5 yrs")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray70", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
pmclfo5
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/confusion_matrices/sens_a/sensa_LFO_MCMC_5yrs.png", plot=pmclfo5)











#========================================================================================================
#sensitivity smax scenario double alpha
#read in data
simPar <- read.csv("data/Smax_sensitivity_doublealpha/SimPars.csv")

## Store relevant object names to help run simulation 
scenNames <- unique(simPar$scenario)

res_smaxda<-readRDS(file = "outs/simest/Smax_sensitivity_doublealpha/res_smaxda.rds")
#res_smaxda56<-readRDS(file = "outs/simest/res_gamma_alpha/res_smaxda_56.rds")


ressmaxda<-rbind(res_smaxda)

#res_a <- res_a[res_a$convergence==0,]

aic_smaxda=subset(ressmaxda, parameter=='AIC'&method=='MLE')
bic_smaxda=subset(ressmaxda, parameter=='BIC'&method=='MLE')
lfo_smaxda=subset(ressmaxda, parameter=='LFO'&method=='MLE')


lfo_smaxda<-lfo_smaxda[lfo_smaxda$model %in% c("simple","autocorr","rwa","rwb","rwab","hmma","hmmb","hmmab"),]
lfo_smaxda[is.na(lfo_smaxda$est),]<--Inf

lfo_smaxda<-lfo_smaxda[lfo_smaxda$model %in% c("simple","autocorr","rwa","rwb","rwab","hmma","hmmb","hmmab"),]

lfo_smaxda[is.na(lfo_smaxda$est),]<--Inf
aic_smaxda$est[aic_smaxda$convergence>0]<-Inf
bic_smaxda$est[aic_smaxda$convergence>0]<-Inf

scn<-factor(unique(aic_smaxda$scenario), levels=c(
    "trendLinearSmax025_da",
    "trendLinearSmax050_da", 
    "trendLinearSmax150_da", 
    "trendLinearSmax200_da", 
    "trendLinearSmax300_da", 
    "regimeSmax025_da",     
    "regimeSmax050_da",      
    "regimeSmax150_da",      
    "regimeSmax200_da",      
    "regimeSmax300_da"  ) )


EM=c("stationary",
     "autocorr",
     "dynamic.a","dynamic.b","dynamic.ab",
     "regime.a","regime.b","regime.ab")
##Confusion matrices
conf_matrix<-expand.grid(EM=EM,OM=scn)
conf_matrix$w_AIC=NA
conf_matrix$BIC=NA
conf_matrix$LFO=NA

cn1<-list()
cn2<-list()
cn3<-list()


#summarize model selection
aic_set=list()
bic_set=list()
lfo_set=list()


o=0
for(a in seq_along(scn)){

  #AIC
  aica<-subset(aic,scenario==scn[a])
  dim(aica)
  unique(aica$scenario)


   #AIC
  aica<-subset(aic_smaxda,scenario==scn[a])
  aic_set[[a]]=tidyr::spread(aica[,-9],key=model,value=est)
  aic_set[[a]]=aic_set[[a]][c(15,8,12,14,13,9,11,10)] #reorder estimation models

  sc1=apply(aic_set[[a]],1,which.min)
  cn1[[a]]=summary(factor(sc1,levels=seq(1:ncol(aic_set[[a]]))))/1000

  bica=subset(bic_smaxda,scenario==scn[a])
  bic_set[[a]]=tidyr::spread(bica[,-9],key=model,value=est)
  bic_set[[a]]=bic_set[[a]][c(15,8,12,14,13,9,11,10)] #reorder estimation models

  sc2=apply(bic_set[[a]],1,which.min)
  cn2[[a]]=summary(factor(sc2,levels=seq(1:ncol(bic_set[[a]]))))/1000

  lfoa=subset(lfo_smaxda,scenario==scn[a])
  lfo_set[[a]]=tidyr::spread(lfoa[,-9],key=model,value=est)
  lfo_set[[a]]=lfo_set[[a]][c(15,8,12,14,13,9,11,10)] #reorder estimation models

  sc3=apply(lfo_set[[a]],1,which.max)
  cn3[[a]]=summary(factor(sc3,levels=seq(1:ncol(lfo_set[[a]]))))/1000


  myseq<-seq(from=o+1, length.out=length(EM))
  conf_matrix$w_AIC[myseq]<-cn1[[a]]
  conf_matrix$BIC[myseq]<-cn2[[a]]
  conf_matrix$LFO[myseq]<-cn3[[a]]
  o=max(myseq)

  
}


conf_matrix$eqem_om <- dplyr::recode(conf_matrix$OM, 
    "trendLinearSmax025_da"="dynamic.b",
    "trendLinearSmax050_da"="dynamic.b", 
    "trendLinearSmax150_da"="dynamic.b", 
    "trendLinearSmax200_da"="dynamic.b", 
    "trendLinearSmax300_da"="dynamic.b", 
    "regimeSmax025_da"="regime.b",     
    "regimeSmax050_da"="regime.b",      
    "regimeSmax150_da"="regime.b",      
    "regimeSmax200_da"="regime.b",      
    "regimeSmax300_da"="regime.b"
      )   

conf_matrix$eqem_om<-factor(conf_matrix$eqem_om, 
    levels=c("stationary", 
        "autocorr", 
        "dynamic.a", 
        "dynamic.b", 
        "dynamic.ab", 
        "regime.a", 
        "regime.b", 
        "regime.ab"))
conf_matrix$diag<-conf_matrix$eqem_om==conf_matrix$EM





p=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = w_AIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_AIC,2)), vjust = 1,size=6) +
  ggtitle("AIC sens smax double alpha")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray90", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1)) +
  xlab("Simulation Scenario")+ylab("Estimation Model")
p



p=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = BIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(BIC,2)), vjust = 1,size=6) +
  ggtitle("BIC sens smax  double alpha")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray90", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
p



p=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = LFO), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(LFO,2)), vjust = 1, size=6) +
  ggtitle("LFO sens smax  double alpha")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray90", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
p




#========================================================================================================
#sigma low sensitivity 
#read in data
simPar <- read.csv("data/sigmalow_sensitivity/SimPars.csv")

## Store relevant object names to help run simulation 
scenNames <- unique(simPar$scenario)

res_siglow<-readRDS(file = "outs/simest/sigmalow_sensitivity/res_siglow.rds")

res_siglow36<-readRDS(file = "outs/simest/sigmalow_sensitivity/res_siglow_36.rds")
res_siglow84<-readRDS(file = "outs/simest/sigmalow_sensitivity/res_siglow_84.rds")

ressiglow<-rbind(res_siglow,res_siglow36,res_siglow84)

#res_a <- res_a[res_a$convergence==0,]

aic_siglow=subset(ressiglow, parameter=='AIC'&method=='MLE')
bic_siglow=subset(ressiglow, parameter=='BIC'&method=='MLE')
lfo_siglow=subset(ressiglow, parameter=='LFO'&method=='MLE')


lfo_siglow<-lfo_siglow[lfo_siglow$model %in% c("simple","autocorr","rwa","rwb","rwab","hmma","hmmb","hmmab"),]
lfo_siglow[is.na(lfo_siglow$est),]<--Inf

lfo_siglow<-lfo_siglow[lfo_siglow$model %in% c("simple","autocorr","rwa","rwb","rwab","hmma","hmmb","hmmab"),]

lfo_siglow[is.na(lfo_siglow$est),]<--Inf
aic_siglow$est[aic_siglow$convergence>0]<-Inf
bic_siglow$est[aic_siglow$convergence>0]<-Inf

scn<-factor(unique(aic_siglow$scenario), levels=c(
     "sigmalow_stationary",
     "sigmalow_decLinearProd",        
     "sigmalow_regimeProd",            
     "sigmalow_sineProd",             
     "sigmalow_regimeCap",             
     "sigmalow_decLinearCap",         
     "sigmalow_regimeProdCap",         
     "sigmalow_shiftCap",             
     "sigmalow_decLinearProdshiftCap"
  ) )


EM=c("stationary",
     "autocorr",
     "dynamic.a","dynamic.b","dynamic.ab",
     "regime.a","regime.b","regime.ab")
##Confusion matrices
conf_matrix<-expand.grid(EM=EM,OM=scn)
conf_matrix$w_AIC=NA
conf_matrix$BIC=NA
conf_matrix$LFO=NA

cn1<-list()
cn2<-list()
cn3<-list()


#summarize model selection
aic_set=list()
bic_set=list()
lfo_set=list()

unique(aica$model)
o=0

aica[7282, ]
aica[8002, ]

for(a in seq_along(scn)){

   #AIC
  aica<-subset(aic_siglow,scenario==scn[a])
  aic_set[[a]]=tidyr::spread(aica[,-9],key=model,value=est)
  aic_set[[a]]=aic_set[[a]][c(15,8,12,14,13,9,11,10)] #reorder estimation models

  sc1=apply(aic_set[[a]],1,which.min)
  cn1[[a]]=summary(factor(sc1,levels=seq(1:ncol(aic_set[[a]]))))/1000

  bica=subset(bic_siglow,scenario==scn[a])
  bic_set[[a]]=tidyr::spread(bica[,-9],key=model,value=est)
  bic_set[[a]]=bic_set[[a]][c(15,8,12,14,13,9,11,10)] #reorder estimation models

  sc2=apply(bic_set[[a]],1,which.min)
  cn2[[a]]=summary(factor(sc2,levels=seq(1:ncol(bic_set[[a]]))))/1000

  lfoa=subset(lfo_siglow,scenario==scn[a])
  lfo_set[[a]]=tidyr::spread(lfoa[,-9],key=model,value=est)
  lfo_set[[a]]=lfo_set[[a]][c(15,8,12,14,13,9,11,10)] #reorder estimation models

  sc3=apply(lfo_set[[a]],1,which.max)
  cn3[[a]]=summary(factor(sc3,levels=seq(1:ncol(lfo_set[[a]]))))/1000


  myseq<-seq(from=o+1, length.out=length(EM))
  conf_matrix$w_AIC[myseq]<-cn1[[a]]
  conf_matrix$BIC[myseq]<-cn2[[a]]
  conf_matrix$LFO[myseq]<-cn3[[a]]
  o=max(myseq)

  
}


conf_matrix$eqem_om <- dplyr::recode(conf_matrix$OM, 
    "sigmalow_stationary"="stationary",
     "sigmalow_decLinearProd"="dynamic.a",        
     "sigmalow_regimeProd"="regime.a",            
     "sigmalow_sineProd"="dynamic.a",             
     "sigmalow_regimeCap"="regime.b",             
     "sigmalow_decLinearCap"="dynamic.b",         
     "sigmalow_regimeProdCap"="regime.ab",         
     "sigmalow_shiftCap"="regime.b",             
     "sigmalow_decLinearProdshiftCap"="dynamic.ab")
    

conf_matrix$eqem_om<-factor(conf_matrix$eqem_om, 
    levels=c("stationary", 
        "autocorr", 
        "dynamic.a", 
        "dynamic.b", 
        "dynamic.ab", 
        "regime.a", 
        "regime.b", 
        "regime.ab"))
conf_matrix$diag<-conf_matrix$eqem_om==conf_matrix$EM





p=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = w_AIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_AIC,2)), vjust = 1,size=6) +
  ggtitle("AIC sens sigma low")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray90", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1)) +
  xlab("Simulation Scenario")+ylab("Estimation Model")
p



p=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = BIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(BIC,2)), vjust = 1,size=6) +
  ggtitle("BIC sens sigma low")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray90", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
p



p=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = LFO), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(LFO,2)), vjust = 1, size=6) +
  ggtitle("LFO sens smax sigma low")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray90", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
p



#========================================================================================================
#sigma med sensitivity 
#read in data
simPar <- read.csv("data/sigmamed_sensitivity/SimPars.csv")

## Store relevant object names to help run simulation 
scenNames <- unique(simPar$scenario)

res_sigmed<-readRDS(file = "outs/simest/sigmamed_sensitivity/res_sigmed.rds")

#res_sigmed73<-readRDS(file = "outs/simest/res_gamma_alpha/res_sigmed_73.rds")


ressigmed<-rbind(res_sigmed)

#res_a <- res_a[res_a$convergence==0,]

aic_sigmed=subset(ressigmed, parameter=='AIC'&method=='MLE')
bic_sigmed=subset(ressigmed, parameter=='BIC'&method=='MLE')
lfo_sigmed=subset(ressigmed, parameter=='LFO'&method=='MLE')


lfo_sigmed<-lfo_sigmed[lfo_sigmed$model %in% c("simple","autocorr","rwa","rwb","rwab","hmma","hmmb","hmmab"),]
lfo_sigmed[is.na(lfo_sigmed$est),]<--Inf

lfo_sigmed<-lfo_sigmed[lfo_sigmed$model %in% c("simple","autocorr","rwa","rwb","rwab","hmma","hmmb","hmmab"),]

lfo_sigmed[is.na(lfo_sigmed$est),]<--Inf
aic_sigmed$est[aic_sigmed$convergence>0]<-Inf
bic_sigmed$est[aic_sigmed$convergence>0]<-Inf

scn<-factor(unique(aic_sigmed$scenario), levels=c(
     "sigmamed_stationary",
     "sigmamed_decLinearProd",        
     "sigmamed_regimeProd",            
     "sigmamed_sineProd",             
     "sigmamed_regimeCap",             
     "sigmamed_decLinearCap",         
     "sigmamed_regimeProdCap",         
     "sigmamed_shiftCap",             
     "sigmamed_decLinearProdshiftCap"
  ) )


EM=c("stationary",
     "autocorr",
     "dynamic.a","dynamic.b","dynamic.ab",
     "regime.a","regime.b","regime.ab")
##Confusion matrices
conf_matrix<-expand.grid(EM=EM,OM=scn)
conf_matrix$w_AIC=NA
conf_matrix$BIC=NA
conf_matrix$LFO=NA

cn1<-list()
cn2<-list()
cn3<-list()


#summarize model selection
aic_set=list()
bic_set=list()
lfo_set=list()


o=0
for(a in seq_along(scn)){

   #AIC
  aica<-subset(aic_sigmed,scenario==scn[a])
  aic_set[[a]]=tidyr::spread(aica[,-9],key=model,value=est)
  aic_set[[a]]=aic_set[[a]][c(15,8,12,14,13,9,11,10)] #reorder estimation models

  sc1=apply(aic_set[[a]],1,which.min)
  cn1[[a]]=summary(factor(sc1,levels=seq(1:ncol(aic_set[[a]]))))/1000

  bica=subset(bic_sigmed,scenario==scn[a])
  bic_set[[a]]=tidyr::spread(bica[,-9],key=model,value=est)
  bic_set[[a]]=bic_set[[a]][c(15,8,12,14,13,9,11,10)] #reorder estimation models

  sc2=apply(bic_set[[a]],1,which.min)
  cn2[[a]]=summary(factor(sc2,levels=seq(1:ncol(bic_set[[a]]))))/1000

  lfoa=subset(lfo_sigmed,scenario==scn[a])
  lfo_set[[a]]=tidyr::spread(lfoa[,-9],key=model,value=est)
  lfo_set[[a]]=lfo_set[[a]][c(15,8,12,14,13,9,11,10)] #reorder estimation models

  sc3=apply(lfo_set[[a]],1,which.max)
  cn3[[a]]=summary(factor(sc3,levels=seq(1:ncol(lfo_set[[a]]))))/1000


  myseq<-seq(from=o+1, length.out=length(EM))
  conf_matrix$w_AIC[myseq]<-cn1[[a]]
  conf_matrix$BIC[myseq]<-cn2[[a]]
  conf_matrix$LFO[myseq]<-cn3[[a]]
  o=max(myseq)

  
}


conf_matrix$eqem_om <- dplyr::recode(conf_matrix$OM, 
    "sigmamed_stationary"="stationary",
     "sigmamed_decLinearProd"="dynamic.a",        
     "sigmamed_regimeProd"="regime.a",            
     "sigmamed_sineProd"="dynamic.a",             
     "sigmamed_regimeCap"="regime.b",             
     "sigmamed_decLinearCap"="dynamic.b",         
     "sigmamed_regimeProdCap"="regime.ab",         
     "sigmamed_shiftCap"="regime.b",             
     "sigmamed_decLinearProdshiftCap"="dynamic.ab")
    

conf_matrix$eqem_om<-factor(conf_matrix$eqem_om, 
    levels=c("stationary", 
        "autocorr", 
        "dynamic.a", 
        "dynamic.b", 
        "dynamic.ab", 
        "regime.a", 
        "regime.b", 
        "regime.ab"))
conf_matrix$diag<-conf_matrix$eqem_om==conf_matrix$EM





p=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = w_AIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_AIC,2)), vjust = 1,size=6) +
  ggtitle("AIC sens sigma med")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray90", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1)) +
  xlab("Simulation Scenario")+ylab("Estimation Model")
p



p=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = BIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(BIC,2)), vjust = 1,size=6) +
  ggtitle("BIC sens sigma med")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray90", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
p



p=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = LFO), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(LFO,2)), vjust = 1, size=6) +
  ggtitle("LFO sens smax sigma med")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray90", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
p





#routines to visualize smulation estimation Results
#Catarina Wor
#November 2022
#============================================



library(gridExtra)
library(ggplot2)
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
simPar <- read.csv("data/sigmalow_sensitivity/SimPars.csv")

## Store relevant object names to help run simulation 
scenNames <- unique(simPar$scenario)

res1<-readRDS(file = "outs/simest/sigmalow_sensitivity/res_siglow1.rds")
res2<-readRDS(file = "outs/simest/sigmalow_sensitivity/res_siglow2.rds")
res105<-readRDS(file = "outs/simest/sigmalow_sensitivity/res_siglow_105.rds")
res<-rbind(res1,res2,res105)#,resstan16,resstan712)

unique(res$scenario)
res<-res[res$convergence==0,]



res$scenario<-case_match(
  res$scenario,
  "sigmalow_stationary"~"stationary",
  "sigmalow_decLinearProd"~"decLinearProd",
  "sigmalow_regimeProd"~ "regimeProd",  
 "sigmalow_sineProd"~ "sineProd",     
 "sigmalow_regimeCap"~ "regimeCap",    
 "sigmalow_decLinearCap"~ "decLinearCap", 
 "sigmalow_regimeProdCap"~ "regimeProdCap",
 "sigmalow_shiftCap"~"shiftCap")

 
unique(res$scenario)

aic=subset(res,parameter=='AIC'&method=='MLE')
bic=subset(res,parameter=='BIC'&method=='MLE')
lfo=subset(res,parameter=='LFO'&method=='MLE')

lfo<-lfo[lfo$model %in% c("simple","autocorr","rwa","rwb","rwab",
    "hmma","hmmb","hmmab"),]
unique(lfo$model)


#lfo$model<-dplyr::recode(lfo$model, 
#      "simple"="simple",
#      "autocorr"="autocorr",
#      "rwa_last5"="rwa",
#      "rwb_last5"="rwb",
#      "rwab_last5"="rwab",
#    "hmma_last5"="hmma",
#    "hmmb_last5"="hmmb",
#    "hmmab_last5"="hmmab"
#      )   



scn<-factor(unique(aic$scenario), levels=c(
  "stationary", 
  #"autocorr",
  "decLinearProd",  
  "sineProd",
  "regimeProd", 
  "decLinearCap",
  "regimeCap", 
  "shiftCap",
  #"decLinearProdshiftCap",
  "regimeProdCap" ) )


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

aic_set<-list()
bic_set<-list()
lfo_set<-list()
o=0
for(a in seq_along(scn)){


  #head(aic_set[[a]][c(15,8,12,9,14,11,13,10)])
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
      #"autocorr"="autocorr",
      "decLinearProd"="dynamic.a",
      "sineProd"="dynamic.a",
      "decLinearCap"="dynamic.b",
      #"decLinearProdshiftCap"="dynamic.ab",
      "regimeProd"="regime.a",
      "regimeCap"="regime.b",
      "shiftCap"="regime.b", 
      "regimeProdCap"="regime.ab",
      ) 

      conf_matrix$eqem_om<-factor(conf_matrix$eqem_om,levels=c("stationary",
     "autocorr",
     "dynamic.a","regime.a",
     "dynamic.b","regime.b",
     "dynamic.ab",
     "regime.ab"))  
conf_matrix$diag<-conf_matrix$eqem_om==conf_matrix$EM

unique(conf_matrix$EM)





paic=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = w_AIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_AIC,2)), vjust = 1,size=6) +
  ggtitle("AIC MLE siglow")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray90", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1)) +
  xlab("Simulation Scenario")+ylab("Estimation Model")
paic

ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/confusion_matrices/sig_low/AIC_MLE_siglow.png",
 plot=paic)



pbic=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = BIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(BIC,2)), vjust = 1,size=6) +
  ggtitle("BIC MLE siglow")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray90", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
pbic

ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/confusion_matrices/sig_low/BIC_MLE_siglow.png", plot=pbic)



plfo=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = LFO), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(LFO,2)), vjust = 1, size=6) +
  ggtitle("LFO MLE siglow")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray90", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
plfo

ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/confusion_matrices/sig_low/LFO_MLE_siglow.png", plot=plfo)


#---------------------------------------------------------------------------------------------
#stan


resstan<-readRDS(file = "outs/simest/sigmalow_sensitivity/resstan_siglow.rds")
#resstan2<-readRDS(file = "outs/simest/sigmalow_sensitivity/resstan2.rds")
#resstan<-rbind(resstan1,resstan2)
head(resstan)


unique(resstan$parameter)

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
      "decLinearProd"="dynamic.a",
      "sineProd"="dynamic.a",
      "decLinearCap"="dynamic.b",
      "decLinearProdshiftCap"="dynamic.ab",
      "regimeProd"="regime.a",
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

ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/confusion_matrices/base/AIC_MCMC.png", plot=pmcaic)






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
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/confusion_matrices/base/BIC_MCMC.png", plot=pmcbic)


#===================================================================================
reslfo<-readRDS(file = "outs/simest/sigmalow_sensitivity/resstanloo_siglow.rds")
head(reslfo)



reslfo$scenario<-case_match(
 reslfo$scenario,
   "sigmalow_stationary"~"stationary",
   "sigmalow_decLinearProd"~"decLinearProd",        
 "sigmalow_regimeProd"~ "regimeProd",            
 "sigmalow_sineProd"~ "sineProd",             
 "sigmalow_regimeCap"~ "regimeCap",             
 "sigmalow_decLinearCap"~ "decLinearCap",        
 "sigmalow_regimeProdCap"~ "regimeProdCap",
 "sigmalow_shiftCap"~"shiftCap",             
 "sigmalow_decLinearProdshiftCap"~"decLinearProdshiftCap")

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
  "decLinearProd",  
  "sineProd",
  "regimeProd",
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
  head(lfomcmc_set[[a]])

  sc3=apply(lfomcmc_set[[a]],1,which.max)
  cn3[[a]]=summary(factor(sc3,levels=seq(1:ncol(lfomcmc_set[[a]]))))/nsim


  myseq<-seq(from=o+1, length.out=length(EM))
  conf_matrix$LFOmcmc[myseq]<-cn3[[a]]
  o=max(myseq)

}



conf_matrix$eqem_om <- dplyr::recode(conf_matrix$OM, 
      "stationary"="stationary",
      "decLinearProd"="dynamic.a",
      "sineProd"="dynamic.a",
      "decLinearCap"="dynamic.b",
      "decLinearProdshiftCap"="dynamic.ab",
      "regimeProd"="regime.a",
      "regimeCap"="regime.b",
      "shiftCap"="regime.b", 
      "regimeProdCap"="regime.ab",
      )   

   conf_matrix$eqem_om<-factor(conf_matrix$eqem_om,levels=c("stationary",
     "autocorr",
     "dynamic.a","regime.a",
     "dynamic.b","regime.b",
     "dynamic.ab",
     "regime.ab"))  
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
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/confusion_matrices/sig_low/LFO_MCMC_siglow.png", plot=pmclfo)


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
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/confusion_matrices/sig_low/LFO_MCMC_3yrs_siglow.png", plot=pmclfo3)




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
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/confusion_matrices/sig_low/LFO_MCMC_5yrs_siglow.png", plot=pmclfo5)






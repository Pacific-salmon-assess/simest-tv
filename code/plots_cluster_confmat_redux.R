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
#alpha scenarios
#read in data
simPar <- read.csv("data/generic/SimPars.csv")

## Store relevant object names to help run simulation 
scenNames <- unique(simPar$scenario)

res1<-readRDS(file = "outs/simest/generic/resbase1.rds")
res2<-readRDS(file = "outs/simest/generic/resbase2.rds")

res<-rbind(res1,res2)#,resstan16,resstan712)
unique(res$model)
res<-res[res$convergence==0,]


res_tva<-res[res$model%in%c("simple", "autocorr", "rwa","hmma")&
res$scenario%in%c("stationary","autocorr","decLinearProd","regimeProd"),]

#head(res_tva)


aic_tva=subset(res_tva,parameter=='AIC'&method=='MLE')
bic_tva=subset(res_tva,parameter=='BIC'&method=='MLE')
lfo_tva=subset(res_tva,parameter=='LFO'&method=='MLE')



scn<-factor(unique(aic_tva$scenario), levels=c(
  "stationary", 
  "autocorr",
  "decLinearProd",  
  "regimeProd" ) )


EM=c("stationary",
     "autocorr",
     "dynamic.a","regime.a")
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
  aica<-subset(aic_tva,scenario==scn[a])
  aic_set[[a]]=tidyr::spread(aica[,-c(10,11,12,13)],key=model,value=mode)
  naicsim<-length(unique(aic_set[[a]]$iteration))
  aic_set[[a]]=aic_set[[a]][c(11,8,10,9)] #reorder estimation models
  #head(aic_set[[a]])
  sc1=apply(aic_set[[a]],1,which.min)
  
  cn1[[a]]=summary(factor(sc1,levels=seq(1:ncol(aic_set[[a]]))))/naicsim

  bica=subset(bic_tva,scenario==scn[a])
  bic_set[[a]]=tidyr::spread(bica[,-c(10,11,12,13)],key=model,value=mode)
  nbicsim<-length(unique( bic_set[[a]]$iteration))
  bic_set[[a]]=bic_set[[a]][c(11,8,10,9)] #reorder estimation models

  sc2=apply(bic_set[[a]],1,which.min)
 
  cn2[[a]]=summary(factor(sc2,levels=seq(1:ncol(bic_set[[a]]))))/nbicsim

  lfoa=subset(lfo_tva,scenario==scn[a])
  lfo_set[[a]]=tidyr::spread(lfoa[,-c(10,11,12,13)],key=model,value=mode)
  nlfosim <- length(unique( lfo_set[[a]]$iteration))
  lfo_set[[a]]=lfo_set[[a]][c(11,8,10,9)] #reorder estimation models
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
      "decLinearProd"="dynamic.a",
      "regimeProd"="regime.a"
      )   
conf_matrix$diag<-conf_matrix$eqem_om==conf_matrix$EM



head(conf_matrix)

conf_matrix$OM2<-dplyr::case_match(conf_matrix$OM,
  "stationary"~"stationary",
      "autocorr"~"autocorr",
  "decLinearProd"~"linear decline",
  "regimeProd"~"regime increase")


conf_matrix$OM2<-factor(conf_matrix$OM2, levels=c(
  "stationary", 
  "autocorr",
  "linear decline",  
  "regime increase" ) )



paic_alphascn=ggplot(data =  conf_matrix, mapping = aes(x = OM2, y = EM)) +
  geom_tile(aes(fill = w_AIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_AIC,2)), vjust = 1,size=6) +
  ggtitle("AIC")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray90", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1)) +
  xlab("simple and time-varying alpha scenarios")+ylab("Estimation Model")
paic_alphascn

ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/confusion_matrices/base_redux/AIC_alphascn_MLE.png", 
  plot=paic_alphascn)



pbic_alphascn=ggplot(data =  conf_matrix, mapping = aes(x = OM2, y = EM)) +
  geom_tile(aes(fill = BIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(BIC,2)), vjust = 1,size=6) +
  ggtitle("BIC")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray90", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("simple and time-varying alpha scenarios")+ylab("Estimation Model")
pbic_alphascn

ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/confusion_matrices/base_redux/BIC_alphascn_MLE.png", 
  plot=pbic_alphascn)



plfo_alphascn=ggplot(data =  conf_matrix, mapping = aes(x = OM2, y = EM)) +
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
  xlab("simple and time-varying alpha scenarios")+ylab("Estimation Model")
plfo_alphascn

ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/confusion_matrices/base_redux/LFO_alphascn_MLE.png",
 plot=plfo_alphascn)


#---------------------------------------------------------------------------------------------
#stan


#===================================================================================
reslfo<-readRDS(file = "outs/simest/generic/resstanloo.rds")
head(reslfo)


reslfo1<-reslfo[reslfo$model %in% c("simple","autocorr","rwa","rwb","rwab",
    "hmma","hmmb","hmmab"),]
##Confusion matrices
conf_matrix_lfo<-expand.grid(EM=EM,OM=scn)
conf_matrix_lfo$LFOmcmc=NA


cn3<-list()
lfomcmc_set<-list()


o=0
for(a in seq_along(scn)){


  lfoa=subset(reslfo1,scenario==scn[a])
  lfomcmc_set[[a]]=tidyr::spread(lfoa[,-9],key=model,value=est)
  nsim<-length(unique( lfomcmc_set[[a]]$iteration))
  lfomcmc_set[[a]]=lfomcmc_set[[a]][c(15,8,12,9)] #reorder estimation models
  head(lfomcmc_set[[a]])

  sc3=apply(lfomcmc_set[[a]],1,which.max)
  cn3[[a]]=summary(factor(sc3,levels=seq(1:ncol(lfomcmc_set[[a]]))))/nsim


  myseq<-seq(from=o+1, length.out=length(EM))
  conf_matrix_lfo$LFOmcmc[myseq]<-cn3[[a]]
  o=max(myseq)

}


conf_matrix_lfo$eqem_om <- dplyr::recode(conf_matrix_lfo$OM, 
      "stationary"="stationary",
      "autocorr"="autocorr",
      
      "decLinearProd"="dynamic.a",
      
      "regimeProd"="regime.a"
      )   
conf_matrix_lfo$diag<-conf_matrix_lfo$eqem_om==conf_matrix$EM







pmclfo_alphascn=ggplot(data =  conf_matrix_lfo, mapping = aes(x = OM, y = EM)) +
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
pmclfo_alphascn
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/confusion_matrices/base/LFO_alphascn_MCMC.png", plot=pmclfo_alphascn)

#=====================================================================================
#=====================================================================================
#=====================================================================================
#time varying beta scenarios. 


#alpha scenarios
#read in data
simPar <- read.csv("data/generic/SimPars.csv")

## Store relevant object names to help run simulation 
scenNames <- unique(simPar$scenario)

res1<-readRDS(file = "outs/simest/generic/resbase1.rds")
res2<-readRDS(file = "outs/simest/generic/resbase2.rds")

res<-rbind(res1,res2)#,resstan16,resstan712)
unique(res$scenario)
res<-res[res$convergence==0,]


res_tvb<-res[res$model%in%c("simple", "autocorr", "rwb","rwab")&
res$scenario%in%c("stationary","autocorr","decLinearCap","regimeProdCap"),]

#head(res_tva)


aic_tvb=subset(res_tvb,parameter=='AIC'&method=='MLE')
bic_tvb=subset(res_tvb,parameter=='BIC'&method=='MLE')
lfo_tvb=subset(res_tvb,parameter=='LFO'&method=='MLE')




scn_b<-factor(unique(aic_tvb$scenario), levels=c(
  "stationary", 
  "autocorr",
  "decLinearCap",  
  "regimeProdCap" ) )

EM=c("stationary",
     "autocorr",
     "dynamic.b",
     "dynamic.ab")
##Confusion matrices
conf_matrix_tvb<-expand.grid(EM=EM,OM=scn_b)
conf_matrix_tvb$w_AIC=NA
conf_matrix_tvb$BIC=NA
conf_matrix_tvb$LFO=NA

cn1<-list()
cn2<-list()
cn3<-list()

aic_set<-list()
bic_set<-list()
lfo_set<-list()
o=0
for(a in seq_along(scn_b)){

  #AIC 
  aica<-subset(aic_tvb,scenario==scn_b[a])
  aic_set[[a]]=tidyr::spread(aica[,-c(10,11,12,13)],key=model,value=mode)
  naicsim<-length(unique(aic_set[[a]]$iteration))
  aic_set[[a]]=aic_set[[a]][c(11,8,10,9)] #reorder estimation models
  #head(aic_set[[a]])
  sc1=apply(aic_set[[a]],1,which.min)
  
  cn1[[a]]=summary(factor(sc1,levels=seq(1:ncol(aic_set[[a]]))))/naicsim

  bica=subset(bic_tvb,scenario==scn_b[a])
  bic_set[[a]]=tidyr::spread(bica[,-c(10,11,12,13)],key=model,value=mode)
  nbicsim<-length(unique( bic_set[[a]]$iteration))
  bic_set[[a]]=bic_set[[a]][c(11,8,10,9)] #reorder estimation models

  sc2=apply(bic_set[[a]],1,which.min)
 
  cn2[[a]]=summary(factor(sc2,levels=seq(1:ncol(bic_set[[a]]))))/nbicsim

  lfoa=subset(lfo_tvb,scenario==scn_b[a])
  lfo_set[[a]]=tidyr::spread(lfoa[,-c(10,11,12,13)],key=model,value=mode)
  nlfosim <- length(unique( lfo_set[[a]]$iteration))
  lfo_set[[a]]=lfo_set[[a]][c(11,8,10,9)] #reorder estimation models
  #head(lfo_set[[a]])

  sc3=apply(lfo_set[[a]],1,which.max)
  cn3[[a]]=summary(factor(sc3,levels=seq(1:ncol(lfo_set[[a]]))))/ nlfosim


  myseq<-seq(from=o+1, length.out=length(EM))
  conf_matrix_tvb$w_AIC[myseq]<-cn1[[a]]
  conf_matrix_tvb$BIC[myseq]<-cn2[[a]]
  conf_matrix_tvb$LFO[myseq]<-cn3[[a]]
  o=max(myseq)

}


conf_matrix_tvb$eqem_om <- dplyr::recode(conf_matrix_tvb$OM, 
      "stationary"="stationary",
      "autocorr"="autocorr",
      "decLinearCap"="dynamic.b",  
      "regimeProdCap"="dynamic.ab"
      )   
conf_matrix_tvb$diag<-conf_matrix_tvb$eqem_om==conf_matrix_tvb$EM



head(conf_matrix)

conf_matrix_tvb$OM2<-dplyr::case_match(conf_matrix_tvb$OM,
  "stationary"~"stationary",
      "autocorr"~"autocorr",
  "decLinearCap"~"linear decline",
  "regimeProdCap"~"regime both")


conf_matrix_tvb$OM2<-factor(conf_matrix_tvb$OM2, levels=c(
  "stationary", 
  "autocorr",
  "linear decline",  
  "regime both" ) )



paic_smaxscn=ggplot(data =  conf_matrix_tvb, mapping = aes(x = OM2, y = EM)) +
  geom_tile(aes(fill = w_AIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_AIC,2)), vjust = 1,size=6) +
  ggtitle("AIC")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix_tvb, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray90", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1)) +
  xlab("simple and time-varying smax scenarios")+ylab("Estimation Model")
paic_smaxscn

ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/confusion_matrices/base_redux/AIC_smaxscn_MLE.png", 
  plot=paic_smaxscn)



pbic_smaxscn=ggplot(data =  conf_matrix_tvb, mapping = aes(x = OM2, y = EM)) +
  geom_tile(aes(fill = BIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(BIC,2)), vjust = 1,size=6) +
  ggtitle("BIC")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix_tvb, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray90", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("simple and time-varying smax scenarios")+ylab("Estimation Model")
pbic_smaxscn

ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/confusion_matrices/base_redux/BIC_smaxscn_MLE.png", 
  plot=pbic_smaxscn)



plfo_smaxscn=ggplot(data =  conf_matrix_tvb, mapping = aes(x = OM2, y = EM)) +
  geom_tile(aes(fill = LFO), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(LFO,2)), vjust = 1, size=6) +
  ggtitle("LFO MLE")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix_tvb, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray90", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("simple and time-varying smax scenarios")+ylab("Estimation Model")
plfo_smaxscn

ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/confusion_matrices/base_redux/LFO_smaxscn_MLE.png",
 plot=plfo_smaxscn)


#---------------------------------------------------------------------------------------------
#stan


#===================================================================================
reslfo<-readRDS(file = "outs/simest/generic/resstanloo.rds")
head(reslfo)


reslfo1<-reslfo[reslfo$model %in% c("simple","autocorr","rwb","rwab"),]
##Confusion matrices
conf_matrix_lfo<-expand.grid(EM=EM,OM=scn_b)
conf_matrix_lfo$LFOmcmc=NA


cn3<-list()
lfomcmc_set<-list()


o=0
for(a in seq_along(scn_b)){


  lfoa=subset(reslfo1,scenario==scn_b[a])
  lfomcmc_set[[a]]=tidyr::spread(lfoa[,-9],key=model,value=est)
  nsim<-length(unique( lfomcmc_set[[a]]$iteration))
  lfomcmc_set[[a]]=lfomcmc_set[[a]][c(11,8,10,9)] #reorder estimation models
  #head(lfomcmc_set[[a]])

  sc3=apply(lfomcmc_set[[a]],1,which.max)
  cn3[[a]]=summary(factor(sc3,levels=seq(1:ncol(lfomcmc_set[[a]]))))/nsim


  myseq<-seq(from=o+1, length.out=length(EM))
  conf_matrix_lfo$LFOmcmc[myseq]<-cn3[[a]]
  o=max(myseq)

}


conf_matrix_lfo$eqem_om <- dplyr::recode(conf_matrix_lfo$OM, 
      "stationary"="stationary",
      "autocorr"="autocorr",
      
      "decLinearCap"="dynamic.b",
      
      "regimeProdCap"="dynamic.ab"
      )   
conf_matrix_lfo$diag<-conf_matrix_lfo$eqem_om==conf_matrix_lfo$EM







pmclfo_smaxscn=ggplot(data =  conf_matrix_lfo, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = LFOmcmc), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(LFOmcmc,2)), vjust = 1, size=6) +
  ggtitle("LFO MCMC")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix_lfo, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray70", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
pmclfo_smaxscn
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/confusion_matrices/base/LFO_smaxscn_MCMC.png", plot=pmclfo_smaxscn)






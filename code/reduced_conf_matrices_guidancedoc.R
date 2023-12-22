#guideline document plots - reduced confusion matrices#
library(here)
library(gridExtra)
library(ggplot2)
source("code/utils.R")


mytheme = list(
  theme_classic(16)+
    theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
          legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="bold", size=12),
          axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
)


#1: Q2. I2####
#BIC - base scenarios (-om: sineprod, sigmashift, regimeProdCap; -em: dynamic ab, hmm ab)

#base case
#read in data
simPar <- read.csv("data/generic/SimPars.csv")

## Store relevant object names to help run simulation 
scenNames <- unique(simPar$scenario)

res1<-readRDS(file = "outs/cluster/generic/resbase1.rds")
res2<-readRDS(file = "outs/cluster/generic/resbase2.rds")
res<-rbind(res1,res2)

bic=subset(res,parameter=='BIC'&method=='MLE')

bic<-bic[bic$model %in% c("simple","autocorr","rwa","rwb",
                          "hmma","hmmb"),]
bic<-bic[bic$scenario %in% c("stationary", 
                        "autocorr",
                        "decLinearProd",  
                        "regimeProd",
                        "decLinearCap",
#                        "decLinearProdshiftCap",      
                
                         
                        "regimeCap" 
                       ),]

bic$est=bic$mode
unique(bic$model)
unique(bic$scenario)


scn<-factor(unique(bic$scenario), levels=c(
  "stationary", 
  "autocorr",
  "decLinearProd",  
  "regimeProd",
  "decLinearCap",
  "regimeCap") )


EM=c("stationary",
     "autocorr",
     "random walk-Prod",
     "regime-Prod","random walk-Cap","regime-Cap")

conf_matrix<-expand.grid(EM=EM,OM=scn)
conf_matrix$BIC=NA

cn2<-list()

bic_set<-list()
o=0
for(a in seq_along(scn)){

  bica=subset(bic,scenario==scn[a])
  bic_set[[a]]=tidyr::spread(bica[,-c(9:11)],key=model,value=est)
  bic_set[[a]]=bic_set[[a]][c(15,10,13,14,11,12)] #reorder estimation models
  
  sc2=apply(bic_set[[a]],1,which.min)
  cn2[[a]]=summary(factor(sc2,levels=seq(1:ncol(bic_set[[a]]))))/1000

  myseq<-seq(from=o+1, length.out=length(EM))
  conf_matrix$BIC[myseq]<-cn2[[a]]
  o=max(myseq)
  
}

conf_matrix$eqem_om <- dplyr::recode(conf_matrix$OM, 
                                     "stationary"="stationary",
                                     "autocorr"="autocorr",
                                     "decLinearProd"="random walk-Prod",
                                     "decLinearCap"="random walk-Cap",
 #                                    "decLinearProdshiftCap"="random walk-Prod",
                                     "regimeProd"="regime-Prod",
                                     "regimeCap"="regime-Cap",
                                     
)   
conf_matrix$diag<-conf_matrix$eqem_om==conf_matrix$EM

conf_matrix$EM

p=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = BIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(BIC,2)), vjust = 1,size=6) +
  ggtitle("MLE BIC")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                              simulated=as.numeric(OM), 
                              estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray90", linewidth=1)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
p


png(filename=here('outs','figures','reduced matrices','Q2I2_reduced_base_BIC_MLE.png'),width=8,height=6,units='in',res=300)
p
dev.off()



#base case + ab models
simPar <- read.csv("data/generic/SimPars.csv")

## Store relevant object names to help run simulation 
scenNames <- unique(simPar$scenario)

res1<-readRDS(file = "outs/cluster/generic/resbase1.rds")
res2<-readRDS(file = "outs/cluster/generic/resbase2.rds")
res<-rbind(res1,res2)

bic=subset(res,parameter=='BIC'&method=='MLE')

bic<-bic[bic$model %in% c("simple","autocorr","rwa","rwb",
                          "hmma","hmmb","rwab","hmmab"),]
bic<-bic[bic$scenario %in% c("stationary", 
                             "autocorr",
                             "decLinearProd",  
                             "regimeProd",
                             "decLinearCap",
                             "regimeCap",
                             "decLinearProdshiftCap",
                             "regimeProdCap"
                             ),]

bic$est=bic$mode
unique(bic$model)
unique(bic$scenario)


scn<-factor(unique(bic$scenario), levels=c(
  "stationary", 
  "autocorr",
  "decLinearProd",  
  "regimeProd",
  "decLinearCap",
  "regimeCap",
  "decLinearProdshiftCap",
  "regimeProdCap") )


EM=c("stationary",
     "autocorr",
     "random walk-Prod",
     "regime-Prod","random walk-Cap","regime-Cap",'random walk-ProdCap','regime-ProdCap')

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
  
  bica=subset(bic,scenario==scn[a])
  bic_set[[a]]=tidyr::spread(bica[,-c(9:11)],key=model,value=est)
  bic_set[[a]]=bic_set[[a]][c(17,10,14,11,16,13,15,12)] #reorder estimation models
  
  sc2=apply(bic_set[[a]],1,which.min)
  cn2[[a]]=summary(factor(sc2,levels=seq(1:ncol(bic_set[[a]]))))/1000
  
  myseq<-seq(from=o+1, length.out=length(EM))
  conf_matrix$BIC[myseq]<-cn2[[a]]
  o=max(myseq)
  
}

conf_matrix$eqem_om <- dplyr::recode(conf_matrix$OM, 
                                     "stationary"="stationary",
                                     "autocorr"="autocorr",
                                     "decLinearProd"="random walk-Prod",
                                     "decLinearCap"="random walk-Cap",
                                     #                                    "decLinearProdshiftCap"="dynamic.a",
                                     "regimeProd"="regime-Prod",
                                     "regimeCap"="regime-Cap",
                                     "decLinearProdshiftCap"="random walk-ProdCap",
                                     "regimeProdCap"="regime-ProdCap"
                                     
)   
conf_matrix$diag<-conf_matrix$eqem_om==conf_matrix$EM

conf_matrix$EM

p=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = BIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(BIC,2)), vjust = 1,size=6) +
  ggtitle("MLE BIC")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                              simulated=as.numeric(OM), 
                              estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray90", linewidth=1)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
p


png(filename=here('outs','figures','reduced matrices','Q2I2_rwab_hmmab_BIC_MLE.png'),width=8,height=6,units='in',res=300)
p
dev.off()


#2: Q2. I8####
#LFO MCMC - base scenarios (-om: sineprod, sigmashift, regimeProdCap; -em: dynamic ab, hmm ab)

resstan1<-readRDS(file = "outs/cluster/generic/resstanloo1.rds")
resstan2<-readRDS(file = "outs/cluster/generic/resstanloo2.rds")
resstan=rbind(resstan1,resstan2)

head(resstan)

unique(resstan$parameter)

lfo=subset(resstan,parameter=='LFO')
lfo<-lfo[lfo$model %in% c("simple","autocorr","rwa","rwb",
                          "hmma","hmmb"),]
lfo<-lfo[lfo$scenario %in% c("stationary", 
                             "autocorr",
                             "decLinearProd",  
                          
                             #"decLinearProdshiftCap",      
                             "regimeProd",
                             "decLinearCap",
                             "regimeCap"),]

lfo$by=NA

scn<-factor(unique(lfo$scenario), levels=c(
  "stationary", 
  "autocorr",
  "decLinearProd",  
  "regimeProd",
  "decLinearCap",
  "regimeCap") )


EM=c("stationary",
     "autocorr",
     "random walk-Prod",
     "regime-Prod","random walk-Cap","regime-Cap")

##Confusion matrices
conf_matrix<-expand.grid(EM=EM,OM=scn)
conf_matrix$lfo=NA



cn1<-list()


lfo_set<-list()

o=0
for(a in seq_along(scn)){
  
  #lfo
  lfoa<-subset(lfo,scenario==scn[a])
  lfo_set[[a]]=tidyr::spread(lfoa,key=model,value=est)
  lfo_set[[a]]=lfo_set[[a]][c(14,9,12,10,13,11)] #reorder estimation models
  
  sc1=apply(lfo_set[[a]],1,which.max)
  cn1[[a]]=summary(factor(sc1,levels=seq(1:ncol(lfo_set[[a]]))))/1000
  
  
  myseq<-seq(from=o+1, length.out=length(EM))
  conf_matrix$lfo[myseq]<-cn1[[a]]
  
  o=max(myseq)
  
}


conf_matrix$eqem_om <- dplyr::recode(conf_matrix$OM, 
                                     "stationary"="stationary",
                                     "autocorr"="autocorr",
                                     "decLinearProd"="random walk-Prod",
                                     "regimeProd"="regime-Prod",
                                     "decLinearCap"="random walk-Cap",
                                     "regimeCap"="regime-Cap"
)   
conf_matrix$diag<-conf_matrix$eqem_om==conf_matrix$EM

conf_matrix$EM


pmcaic=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = lfo), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(lfo,2)), vjust = 1,size=6) +
  ggtitle(" MCMC LFO-CV")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                              simulated=as.numeric(OM), 
                              estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray90", linewidth=1)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1)) +
  xlab("Simulation Scenario")+ylab("Estimation Model")
pmcaic

png(filename=here('outs','figures','reduced matrices','Q2I8_reduced_base_LFO_MCMC.png'),width=8,height=6,units='in',res=300)
pmcaic
dev.off()

#Q3 I1####
#low, med, high sigmas for base case scenarios (regimeProdCap; -em: dynamic ab, dynamic b, hmm b, hmm ab)

##high sigma, base case####
simPar <- read.csv("data/generic/SimPars.csv")

## Store relevant object names to help run simulation 
scenNames <- unique(simPar$scenario)

res1<-readRDS(file = "outs/cluster/generic/resbase1.rds")
res2<-readRDS(file = "outs/cluster/generic/resbase2.rds")
res<-rbind(res1,res2)

bic=subset(res,parameter=='BIC'&method=='MLE')

bic<-bic[bic$model %in% c("simple","autocorr","rwa","rwb",
                          "hmma","hmmb"),]
bic<-bic[bic$scenario %in% c("stationary", 
                             "autocorr",
                             "decLinearProd",  
                             "decLinearCap",
                             #                        "decLinearProdshiftCap",      
                             "regimeProd",
                             "regimeCap" 
                             ),]

bic$est=bic$mode
unique(bic$model)
unique(bic$scenario)


scn<-factor(unique(bic$scenario), levels=c(
  "stationary", 
  "autocorr",
  "decLinearProd",  
  #  "decLinearProdshiftCap",      
  "regimeProd",
  "decLinearCap",
  "regimeCap") )


EM=c("stationary",
     "autocorr",
     "random walk-Prod",
     "regime-Prod","random walk-Cap","regime-Cap")

conf_matrix<-expand.grid(EM=EM,OM=scn)
conf_matrix$BIC=NA

cn2<-list()

bic_set<-list()
o=0
for(a in seq_along(scn)){
  
  bica=subset(bic,scenario==scn[a])
  bic_set[[a]]=tidyr::spread(bica[,-c(9:11)],key=model,value=est)
  bic_set[[a]]=bic_set[[a]][c(15,10,13,11,14,12)] #reorder estimation models
  
  sc2=apply(bic_set[[a]],1,which.min)
  cn2[[a]]=summary(factor(sc2,levels=seq(1:ncol(bic_set[[a]]))))/1000
  
  myseq<-seq(from=o+1, length.out=length(EM))
  conf_matrix$BIC[myseq]<-cn2[[a]]
  o=max(myseq)
  
}

conf_matrix$eqem_om <- dplyr::recode(conf_matrix$OM, 
                                     "stationary"="stationary",
                                     "autocorr"="autocorr",
                                     "decLinearProd"="random walk-Prod",
                                     "regimeProd"="regime-Prod",
                                     "decLinearCap"="random walk-Cap",
                                     "regimeCap"="regime-Cap",

)   
conf_matrix$diag<-conf_matrix$eqem_om==conf_matrix$EM

conf_matrix$EM



p_sigh=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = BIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(BIC,2)), vjust = 1,size=6) +
  ggtitle("MLE BIC - high sigma (0.7)")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                              simulated=as.numeric(OM), 
                              estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray90", linewidth=1)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")

png(filename=here('outs','figures','reduced matrices','Q3I1_highsigma_MLE_BIC.png'),width=8,height=6,units='in',res=300)
p_sigh
dev.off()


##med sigma, base case####
simPar <- read.csv("data/sigmamed_sensitivity/SimPars.csv")

## Store relevant object names to help run simulation 
scenNames <- unique(simPar$scenario)

res_sigmed1<-readRDS(file = "outs/cluster/sigmamed_sensitivity/res_sigmed1.rds")
res_sigmed2<-readRDS(file = "outs/cluster/sigmamed_sensitivity/res_sigmed2.rds")
res_sigmed12<-readRDS(file = "outs/cluster/sigmamed_sensitivity/res_sigmed_12.rds")

ressigmed<-rbind(res_sigmed1,res_sigmed2,res_sigmed12)

bic=subset(ressigmed,parameter=='BIC'&method=='MLE')

bic<-bic[bic$model %in% c("simple","autocorr","rwa","rwb",
                          "hmma","hmmb"),]
bic$scenario=gsub('sigmamed_','',bic$scenario)

bic<-bic[bic$scenario %in% c("stationary", 
                             "autocorr",
                             "decLinearProd",  
                             "decLinearCap",
                             #                        "decLinearProdshiftCap",      
                             "regimeProd",
                             "regimeCap" 
),]

unique(bic$model)
unique(bic$scenario)

bic$scenario=gsub('sigmamed_','',bic$scenario)

scn<-factor(unique(bic$scenario), levels=c(
  "stationary", 
  "autocorr",
  "decLinearProd",  
  #  "decLinearProdshiftCap",      
  "regimeProd",
  "decLinearCap",
  "regimeCap") )


EM=c("stationary",
     "autocorr",
     "random walk-Prod",
     "regime-Prod","random walk-Cap","regime-Cap")

conf_matrix<-expand.grid(EM=EM,OM=scn)
conf_matrix$BIC=NA

cn2<-list()

bic_set<-list()
o=0
for(a in seq_along(scn)){
  
  bica=subset(bic,scenario==scn[a])
  bica=subset(bica,duplicated(est)==F) #some doubled up results
  bic_set[[a]]=tidyr::spread(bica[,-c(9:10)],key=model,value=est)
  bic_set[[a]]=bic_set[[a]][c(12,7,10,8,11,9)] #reorder estimation models
  sc2=apply(bic_set[[a]],1,which.min)
  cn2[[a]]=summary(factor(sc2,levels=seq(1:ncol(bic_set[[a]]))))/nrow(bic_set[[a]])
  
  myseq<-seq(from=o+1, length.out=length(EM))
  conf_matrix$BIC[myseq]<-cn2[[a]]
  o=max(myseq)
  
}

conf_matrix$eqem_om <- dplyr::recode(conf_matrix$OM, 
                                     "stationary"="stationary",
                                     "autocorr"="autocorr",
                                     "decLinearProd"="random walk-Prod",
                                     "regimeProd"="regime-Prod",
                                     "decLinearCap"="random walk-Cap",
                                     "regimeCap"="regime-Cap",
                                     
)   
conf_matrix$diag<-conf_matrix$eqem_om==conf_matrix$EM

conf_matrix$EM


p_med=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = BIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(BIC,2)), vjust = 1,size=6) +
  ggtitle("MLE BIC - med sigma (0.3)")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                              simulated=as.numeric(OM), 
                              estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray90", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")

png(filename=here('outs','figures','reduced matrices','Q3I1_medsigma_MLE_BIC.png'),width=8,height=6,units='in',res=300)
p_med
dev.off()


#low sigma, base case####
#read in data
simPar <- read.csv("data/sigmalow_sensitivity/SimPars.csv")

## Store relevant object names to help run simulation 
scenNames <- unique(simPar$scenario)

res_siglow1<-readRDS(file = "outs/cluster/sigmalow_sensitivity/res_siglow1.rds")
res_siglow2<-readRDS(file = "outs/cluster/sigmalow_sensitivity/res_siglow2.rds")
res_siglow15<-readRDS(file = "outs/cluster/sigmalow_sensitivity/res_siglow_15.rds")

ressiglow<-rbind(res_siglow1,res_siglow2,res_siglow15)

bic=subset(ressiglow,parameter=='BIC'&method=='MLE')

bic<-bic[bic$model %in% c("simple","autocorr","rwa","rwb",
                          "hmma",'hmmb'),]
bic$scenario=gsub('sigmalow_','',bic$scenario)
bic<-bic[bic$scenario %in% c("stationary", 
                             "autocorr",
                             "decLinearProd",  
                             "regimeProd",
                             "decLinearCap",
                             #                        "decLinearProdshiftCap",      
                             "regimeCap" 
                         
                             
                        ),]


unique(bic$model)
unique(bic$scenario)


scn<-factor(unique(bic$scenario), levels=c(
  "stationary", 
  "autocorr",
  "decLinearProd",  
  "regimeProd",
  "decLinearCap",
  #  "decLinearProdshiftCap",      

  "regimeCap") )


EM=c("stationary",
     "autocorr",
     "random walk-Prod",
     "regime-Prod","random walk-Cap","regime-Cap")

conf_matrix<-expand.grid(EM=EM,OM=scn)
conf_matrix$BIC=NA

cn2<-list()

bic_set<-list()
o=0
for(a in seq_along(scn)){
  
  bica=subset(bic,scenario==scn[a])
  bica=subset(bica,duplicated(est)==F) #some doubled up results
  bic_set[[a]]=tidyr::spread(bica[,-c(9:10)],key=model,value=est)
  bic_set[[a]]=bic_set[[a]][c(12,7,10,8,11,9)] #reorder estimation models
  
  mwbic=apply(wbic,2,mean)
  sc2=apply(bic_set[[a]],1,which.min)
  
  cn2[[a]]=summary(factor(sc2,levels=seq(1:ncol(bic_set[[a]]))))/nrow(bic_set[[a]])
  
  myseq<-seq(from=o+1, length.out=length(EM))
  conf_matrix$BIC[myseq]<-cn2[[a]]
  o=max(myseq)
  
}

conf_matrix$eqem_om <- dplyr::recode(conf_matrix$OM, 
                                     "stationary"="stationary",
                                     "autocorr"="autocorr",
                                     "decLinearProd"="random walk-Prod",
                                     "regimeProd"="regime-Prod",
                                     "decLinearCap"="random walk-Cap",
                                     "regimeCap"="regime-Cap",
                                     
)   
conf_matrix$diag<-conf_matrix$eqem_om==conf_matrix$EM

conf_matrix$EM



p_low=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = BIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(BIC,2)), vjust = 1,size=6) +
  ggtitle("MLE BIC - low sigma (0.1)")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                              simulated=as.numeric(OM), 
                              estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray90", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")

png(filename=here('outs','figures','reduced matrices','Q3I1_lowsigma_MLE_BIC.png'),width=8,height=6,units='in',res=300)
p_low
dev.off()




#ER scenarios ####
#base case
#read in data
simPar <- read.csv("data/genericER/SimPars_ER.csv")

## Store relevant object names to help run simulation 
scenNames <- unique(simPar$scenario)

res1<-readRDS(file = "outs/cluster/genericER/res_er1.rds")
res2<-readRDS(file = "outs/cluster/genericER/res_er2.rds")
res<-rbind(res1,res2)

bic=subset(res,parameter=='BIC'&method=='MLE')

bic<-bic[bic$model %in% c("simple","autocorr","rwa","rwb",
                          "hmma","hmmb"),]

bic$est=bic$mode
unique(bic$model)
unique(bic$scenario)


scn<-factor(unique(bic$scenario))


EM=c("stationary",
     "autocorr",
     "random walk-Prod",
     "regime-Prod","random walk-Cap","regime-Cap")

conf_matrix<-expand.grid(EM=EM,OM=scn)
conf_matrix$BIC=NA

cn2<-list()

bic_set<-list()
o=0
for(a in seq_along(scn)){
  
  bica=subset(bic,scenario==scn[a])
  bica=subset(bica,duplicated(est)==F)
  bic_set[[a]]=tidyr::spread(bica[,-c(9:11)],key=model,value=est)
  bic_set[[a]]=bic_set[[a]][c(15,10,13,11,14,12)] #reorder estimation models
  
  sc2=apply(bic_set[[a]],1,which.min)
  cn2[[a]]=summary(factor(sc2,levels=seq(1:ncol(bic_set[[a]]))))/1000
  
  myseq<-seq(from=o+1, length.out=length(EM))
  conf_matrix$BIC[myseq]<-cn2[[a]]
  o=max(myseq)
  
}

conf_matrix$EM

p=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = BIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(BIC,2)), vjust = 1,size=6) +
  ggtitle("MLE BIC")+
  scale_fill_gradient(low="white", high="#009194")  +
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
p


png(filename=here('outs','figures','reduced matrices','ERscn_reduced_base_BIC_MLE.png'),width=12,height=6,units='in',res=300)
p
dev.off()


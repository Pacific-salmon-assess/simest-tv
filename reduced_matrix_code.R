#routines to visualize smulation estimation Results
#Catarina Wor
#November 2022




library(gridExtra)
library(ggplot2)
source("code/utils.R")


mytheme = list(
  theme_classic(16)+
    theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
          legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="bold", size=12),
          axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
)



#base case####
#read in data
simPar <- read.csv("data/generic/SimPars.csv")

## Store relevant object names to help run simulation 
scenNames <- unique(simPar$scenario)

res1<-readRDS(file = "outs/simest/generic/resbase1.rds")
res2<-readRDS(file = "outs/simest/generic/resbase2.rds")

res<-rbind(res1,res2)#,resstan16,resstan712)

res<-res[res$convergence==0,]

aic=subset(res,parameter=='AIC'&method=='MLE')
bic=subset(res,parameter=='BIC'&method=='MLE')
lfo=subset(res,parameter=='LFO'&method=='MLE')

lfo<-lfo[lfo$model %in% c("simple","autocorr","rwa","rwb"),]
unique(lfo$model)


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
     "dynamic.a",
     "dynamic.b")
##Confusion matrices
conf_matrix<-expand.grid(EM=EM,OM=scn)
conf_matrix$w_AIC=NA
conf_matrix$BIC=NA
conf_matrix$LFO=NA

cn1<-list()
cn1.2<-list()
cn1.5<-list()
cn2<-list()
cn2.2<-list()
cn2.5<-list()
cn3<-list()

aic_set<-list()
bic_set<-list()
lfo_set<-list()
o=0

#rwa
for(a in seq_along(scn)){
  
  #AIC
  aica<-subset(aic,scenario==scn[a])
  aic_set[[a]]=tidyr::spread(aica[,-c(10,11,12,13)],key=model,value=mode)
  naicsim<-length(unique(aic_set[[a]]$iteration))
  aic_set[[a]]=aic_set[[a]][c(15,8,12,14)] #reorder estimation models
  #head(aic_set[[a]])
  
  #set1 = min. AIC
  sc1=apply(aic_set[[a]],1,which.min)

  #set2 = static parsimony - delta 2 AIC
  #delta AIC
  wAIC=t(apply(aic_set[[a]],1,samEst::model_weights,form='AIC'))
  dAIC=  t(apply(aic_set[[a]],1,function(x) x-min(x)))
  sc1.2=sc1
  sc1.2=ifelse(sc1.2==2&dAIC[,1]<2,1,sc1.2)
  sc1.2=ifelse(sc1.2==3&dAIC[,1]<2,1,sc1.2)
  
  #set3 = static parsimony - delta 5 AIC
  sc1.5=sc1
  sc1.5=ifelse(sc1.5==2&dAIC[,1]<5,1,sc1.5)
  sc1.5=ifelse(sc1.5==3&dAIC[,1]<5,1,sc1.5)
  
  cn1[[a]]=summary(factor(sc1,levels=seq(1:ncol(aic_set[[a]]))))/naicsim
  cn1.2[[a]]=summary(factor(sc1.2,levels=seq(1:ncol(aic_set[[a]]))))/naicsim
  cn1.5[[a]]=summary(factor(sc1.5,levels=seq(1:ncol(aic_set[[a]]))))/naicsim
  
  bica=subset(bic,scenario==scn[a])
  bic_set[[a]]=tidyr::spread(bica[,-c(10,11,12,13)],key=model,value=mode)
  nbicsim<-length(unique( bic_set[[a]]$iteration))
  bic_set[[a]]=bic_set[[a]][c(15,8,12,14)] #reorder estimation models
  
  sc2=apply(bic_set[[a]],1,which.min)
  
  #set2 = static parsimony - delta 2 BIC
  #delta AIC
  wBIC=t(apply(bic_set[[a]],1,samEst::model_weights,form='AIC'))
  dBIC=  t(apply(bic_set[[a]],1,function(x) x-min(x)))
  sc2.2=sc2
  sc2.2=ifelse(sc1.2==2&dBIC[,1]<2,1,sc1.2)
  sc2.2=ifelse(sc1.2==3&dBIC[,1]<2,1,sc1.2)
  
  #set3 = static parsimony - delta 5 BIC
  sc2.5=sc2
  sc2.5=ifelse(sc2.5==2&dBIC[,1]<5,1,sc2.5)
  sc2.5=ifelse(sc2.5==3&dBIC[,1]<5,1,sc2.5)
  
  cn2[[a]]=summary(factor(sc2,levels=seq(1:ncol(bic_set[[a]]))))/nbicsim
  cn2.2[[a]]=summary(factor(sc2.2,levels=seq(1:ncol(bic_set[[a]]))))/nbicsim
  cn2.5[[a]]=summary(factor(sc2.5,levels=seq(1:ncol(bic_set[[a]]))))/nbicsim
  
  lfoa=subset(lfo,scenario==scn[a])
  lfo_set[[a]]=tidyr::spread(lfoa[,-c(10,11,12,13)],key=model,value=mode)
  nlfosim <- length(unique( lfo_set[[a]]$iteration))
  lfo_set[[a]]=lfo_set[[a]][c(10,8,9)] #reorder estimation models
  #head(lfo_set[[a]])
  
  sc3=apply(lfo_set[[a]],1,which.max)
  cn3[[a]]=summary(factor(sc3,levels=seq(1:ncol(lfo_set[[a]]))))/ nlfosim
  wLFO=t(apply(lfo_set[[a]],1,samEst::model_weights,form='AIC'))
  
  #use weight thresholds - 0.55 ~ 2 dAIC
 # sc3.2=sc3
#  sc3.2=ifelse(na.omit(wLFO[,2])<0.55,1,sc3.2)
#  sc3.2=ifelse(sc3.2==3&wLFO[,3]<0.55,1,sc3.2)
  
  
  myseq<-seq(from=o+1, length.out=length(EM))
  conf_matrix$w_AIC[myseq]<-cn1[[a]]
  conf_matrix$w_AIC.d2[myseq]<-cn1.2[[a]]
  conf_matrix$w_AIC.d5[myseq]<-cn1.5[[a]]
  conf_matrix$BIC[myseq]<-cn2[[a]]
  conf_matrix$BIC.d2[myseq]<-cn2.2[[a]]
  conf_matrix$BIC.d5[myseq]<-cn2.5[[a]]
  conf_matrix$LFO[myseq]<-cn3[[a]]
  o=max(myseq)
}

paic=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = w_AIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_AIC,2)), vjust = 1,size=6) +
  ggtitle("AIC")+
  scale_fill_gradient(low="white", high="#009194")+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1)) +
  xlab("Simulation Scenario")+ylab("Estimation Model")
paic

paic2=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = w_AIC.d2), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_AIC.d2,2)), vjust = 1,size=6) +
  ggtitle("AIC - min. d2")+
  scale_fill_gradient(low="white", high="#009194")+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1)) +
  xlab("Simulation Scenario")+ylab("Estimation Model")
paic2

paic3=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = w_AIC.d5), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_AIC.d5,2)), vjust = 1,size=6) +
  ggtitle("AIC - min. d5")+
  scale_fill_gradient(low="white", high="#009194")+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1)) +
  xlab("Simulation Scenario")+ylab("Estimation Model")
paic3

pbic=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = BIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(BIC,2)), vjust = 1,size=6) +
  ggtitle("BIC")+
  scale_fill_gradient(low="white", high="#009194")  +

  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
pbic

pbic2=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = BIC.d2), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(BIC.d2,2)), vjust = 1,size=6) +
  ggtitle("BIC - min. d2")+
  scale_fill_gradient(low="white", high="#009194")  +
  
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
pbic2

pbic3=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = BIC.d5), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(BIC.d5,2)), vjust = 1,size=6) +
  ggtitle("BIC - min. d5")+
  scale_fill_gradient(low="white", high="#009194")  +
  
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
pbic3


plfo=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = LFO), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(LFO,2)), vjust = 1, size=6) +
  ggtitle("LFO MLE")+
  scale_fill_gradient(low="white", high="#009194")  +
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
plfo

plfo=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = LFO), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(LFO,2)), vjust = 1, size=6) +
  ggtitle("LFO MLE")+
  scale_fill_gradient(low="white", high="#009194")  +
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
plfo

ggsave(here::here('outs','figures','aic_base.jpg')
, plot=paic,width=10)
ggsave(here::here('outs','figures','aic_base_d2.jpg')
       , plot=paic2,width=10)
ggsave(here::here('outs','figures','aic_base_d5.jpg')
       , plot=paic3,width=10)
ggsave(here::here('outs','figures','bic_base.jpg')
       , plot=pbic,width=10)
ggsave(here::here('outs','figures','bic_base_d2.jpg')
       , plot=pbic2,width=10)
ggsave(here::here('outs','figures','bic_base_d5.jpg')
       , plot=pbic3,width=10)
ggsave(here::here('outs','figures','lfo_base.jpg')
       , plot=plfo,width=10)

#mcmc lfo####
simPar <- read.csv("data/generic/SimPars.csv")

## Store relevant object names to help run simulation 
scenNames <- unique(simPar$scenario)

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

##Confusion matrices
conf_matrix<-expand.grid(EM=EM,OM=scn)
conf_matrix$w_AIC=NA
conf_matrix$BIC=NA
conf_matrix$LFO=NA


cn3<-list()

lfo_set<-list()
o=0
#rwa
for(a in seq_along(scn)){

  
  lfoa=subset(reslfo,scenario==scn[a])
  lfo_set[[a]]=tidyr::spread(lfoa[,-c(10,11,12,13)],key=model,value=est)
  nlfosim <- length(unique( lfo_set[[a]]$iteration))
  lfo_set[[a]]=lfo_set[[a]][c(10,8,9)] #reorder estimation models
  #head(lfo_set[[a]])
  
  sc3=apply(lfo_set[[a]],1,which.max)
  cn3[[a]]=summary(factor(sc3,levels=seq(1:ncol(lfo_set[[a]]))))/ nlfosim
  wLFO=t(apply(lfo_set[[a]],1,samEst::model_weights,form='AIC'))
  
  #use weight thresholds - 0.55 ~ 2 dAIC
  # sc3.2=sc3
  #  sc3.2=ifelse(na.omit(wLFO[,2])<0.55,1,sc3.2)
  #  sc3.2=ifelse(sc3.2==3&wLFO[,3]<0.55,1,sc3.2)
  myseq<-seq(from=o+1, length.out=length(EM))
  conf_matrix$LFO[myseq]<-cn3[[a]]
  o=max(myseq)
}

plfo=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = LFO), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(LFO,2)), vjust = 1, size=6) +
  ggtitle("LFO MLE")+
  scale_fill_gradient(low="white", high="#009194")  +
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
plfo

ggsave(here::here('outs','figures','lfo_base_mcmc.jpg')
       , plot=plfo,width=10)

#sensitivty a####
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


lfo<-lfo[lfo$model %in% c("simple","autocorr","rwa"),]
unique(lfo$model)

scn<-factor(unique(aic$scenario), levels=c(
  "trendLinearProd1", 
  "trendLinearProd2", 
  "trendLinearProd5",
  "trendLinearProd7",
  "trendLinearProd10",
  "regimeProd1",
  "regimeProd2",      
  "regimeProd5",      
  "regimeProd7",
  "regimeProd10" ) )

conf_matrix<-expand.grid(EM=EM,OM=scn)
conf_matrix$w_AIC=NA
conf_matrix$BIC=NA
conf_matrix$LFO=NA

cn1<-list()
cn1.2<-list()
cn1.5<-list()
cn2<-list()
cn2.2<-list()
cn2.5<-list()
cn3<-list()

aic_set<-list()
bic_set<-list()
lfo_set<-list()
o=0

#rwa
for(a in seq_along(scn)){
  
  #AIC
  aica<-subset(aic,scenario==scn[a])
  aic_set[[a]]=tidyr::spread(aica[,-c(10,11,12,13)],key=model,value=mode)
  naicsim<-length(unique(aic_set[[a]]$iteration))
  aic_set[[a]]=aic_set[[a]][c(15,8,12)] #reorder estimation models
  #head(aic_set[[a]])
  
  #set1 = min. AIC
  sc1=apply(aic_set[[a]],1,which.min)
  
  #set2 = static parsimony - delta 2 AIC
  #delta AIC
  wAIC=t(apply(aic_set[[a]],1,samEst::model_weights,form='AIC'))
  dAIC=  t(apply(aic_set[[a]],1,function(x) x-min(x)))
  sc1.2=sc1
  sc1.2=ifelse(sc1.2==2&dAIC[,1]<2,1,sc1.2)
  sc1.2=ifelse(sc1.2==3&dAIC[,1]<2,1,sc1.2)
  
  #set3 = static parsimony - delta 5 AIC
  sc1.5=sc1
  sc1.5=ifelse(sc1.5==2&dAIC[,1]<5,1,sc1.5)
  sc1.5=ifelse(sc1.5==3&dAIC[,1]<5,1,sc1.5)
  
  cn1[[a]]=summary(factor(sc1,levels=seq(1:ncol(aic_set[[a]]))))/naicsim
  cn1.2[[a]]=summary(factor(sc1.2,levels=seq(1:ncol(aic_set[[a]]))))/naicsim
  cn1.5[[a]]=summary(factor(sc1.5,levels=seq(1:ncol(aic_set[[a]]))))/naicsim
  
  bica=subset(bic,scenario==scn[a])
  bic_set[[a]]=tidyr::spread(bica[,-c(10,11,12,13)],key=model,value=mode)
  nbicsim<-length(unique( bic_set[[a]]$iteration))
  bic_set[[a]]=bic_set[[a]][c(15,8,12)] #reorder estimation models
  
  sc2=apply(bic_set[[a]],1,which.min)
  
  #set2 = static parsimony - delta 2 BIC
  #delta AIC
  wBIC=t(apply(bic_set[[a]],1,samEst::model_weights,form='AIC'))
  dBIC=  t(apply(bic_set[[a]],1,function(x) x-min(x)))
  sc2.2=sc2
  sc2.2=ifelse(sc1.2==2&dBIC[,1]<2,1,sc1.2)
  sc2.2=ifelse(sc1.2==3&dBIC[,1]<2,1,sc1.2)
  
  #set3 = static parsimony - delta 5 BIC
  sc2.5=sc2
  sc2.5=ifelse(sc2.5==2&dBIC[,1]<5,1,sc2.5)
  sc2.5=ifelse(sc2.5==3&dBIC[,1]<5,1,sc2.5)
  
  cn2[[a]]=summary(factor(sc2,levels=seq(1:ncol(bic_set[[a]]))))/nbicsim
  cn2.2[[a]]=summary(factor(sc2.2,levels=seq(1:ncol(bic_set[[a]]))))/nbicsim
  cn2.5[[a]]=summary(factor(sc2.5,levels=seq(1:ncol(bic_set[[a]]))))/nbicsim
  
  lfoa=subset(lfo,scenario==scn[a])
  lfo_set[[a]]=tidyr::spread(lfoa[,-c(10,11,12,13)],key=model,value=mode)
  nlfosim <- length(unique( lfo_set[[a]]$iteration))
  lfo_set[[a]]=lfo_set[[a]][c(10,8,9)] #reorder estimation models
  #head(lfo_set[[a]])
  
  sc3=apply(lfo_set[[a]],1,which.max)
  cn3[[a]]=summary(factor(sc3,levels=seq(1:ncol(lfo_set[[a]]))))/ nlfosim
  wLFO=t(apply(lfo_set[[a]],1,samEst::model_weights,form='AIC'))
  
  #use weight thresholds - 0.55 ~ 2 dAIC
  # sc3.2=sc3
  #  sc3.2=ifelse(na.omit(wLFO[,2])<0.55,1,sc3.2)
  #  sc3.2=ifelse(sc3.2==3&wLFO[,3]<0.55,1,sc3.2)
  
  
  myseq<-seq(from=o+1, length.out=length(EM))
  conf_matrix$w_AIC[myseq]<-cn1[[a]]
  conf_matrix$w_AIC.d2[myseq]<-cn1.2[[a]]
  conf_matrix$w_AIC.d5[myseq]<-cn1.5[[a]]
  conf_matrix$BIC[myseq]<-cn2[[a]]
  conf_matrix$BIC.d2[myseq]<-cn2.2[[a]]
  conf_matrix$BIC.d5[myseq]<-cn2.5[[a]]
  conf_matrix$LFO[myseq]<-cn3[[a]]
  o=max(myseq)
}


paic=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = w_AIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_AIC,2)), vjust = 1,size=6) +
  ggtitle("AIC")+
  scale_fill_gradient(low="white", high="#009194")+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1)) +
  xlab("Simulation Scenario")+ylab("Estimation Model")
paic

paic2=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = w_AIC.d2), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_AIC.d2,2)), vjust = 1,size=6) +
  ggtitle("AIC - min. d2")+
  scale_fill_gradient(low="white", high="#009194")+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1)) +
  xlab("Simulation Scenario")+ylab("Estimation Model")
paic2

paic5=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = w_AIC.d5), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_AIC.d5,2)), vjust = 1,size=6) +
  ggtitle("AIC - min. d5")+
  scale_fill_gradient(low="white", high="#009194")+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1)) +
  xlab("Simulation Scenario")+ylab("Estimation Model")
paic5

pbic=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = BIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(BIC,2)), vjust = 1,size=6) +
  ggtitle("BIC")+
  scale_fill_gradient(low="white", high="#009194")  +
  
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
pbic

pbic2=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = BIC.d2), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(BIC.d2,2)), vjust = 1,size=6) +
  ggtitle("BIC - min. d2")+
  scale_fill_gradient(low="white", high="#009194")  +
  
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
pbic2

pbic3=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = BIC.d5), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(BIC.d5,2)), vjust = 1,size=6) +
  ggtitle("BIC - min. d5")+
  scale_fill_gradient(low="white", high="#009194")  +
  
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
pbic3


plfo=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = LFO), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(LFO,2)), vjust = 1, size=6) +
  ggtitle("LFO MLE")+
  scale_fill_gradient(low="white", high="#009194")  +
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
plfo

ggsave(here::here('outs','figures','aic_sensa.jpg')
       , plot=paic,width=10)
ggsave(here::here('outs','figures','aic_sensa_d2.jpg')
       , plot=paic2,width=10)
ggsave(here::here('outs','figures','aic_sensa_d5.jpg')
       , plot=paic3,width=10)
ggsave(here::here('outs','figures','bic_sensa.jpg')
       , plot=pbic,width=10)
ggsave(here::here('outs','figures','bic_sensa_d2.jpg')
       , plot=pbic2,width=10)
ggsave(here::here('outs','figures','bic_sensa_d5.jpg')
       , plot=pbic3,width=10)
ggsave(here::here('outs','figures','lfo_sensa.jpg')
       , plot=plfo,width=10)

#sensitivity b####
simPar <- read.csv("data/Smax_sensitivity/SimPars.csv")

## Store relevant object names to help run simulation 
scenNames <- unique(simPar$scenario)

res_smax<-readRDS(file = "outs/results/Smax_sensitivity/res_smax.rds")
#res_smax95<-readRDS(file = "outs/simest/res_gamma_alpha/res_smax_95.rds")


ressmax<-rbind(res_smax)

#res_a <- res_a[res_a$convergence==0,]

aic=subset(res_a, parameter=='AIC'&method=='MLE')
bic=subset(res_a, parameter=='BIC'&method=='MLE')
lfo=subset(res_a, parameter=='LFO'&method=='MLE')


lfo<-lfo[lfo$model %in% c("simple","autocorr","rwa"),]
unique(lfo$model)

#lfo[is.na(lfo$est),]<--Inf
#aic$est[aic_s$convergence>0]<-Inf
#bic$est[aic_smax$convergence>0]<-Inf

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

conf_matrix<-expand.grid(EM=EM,OM=scn)
conf_matrix$w_AIC=NA
conf_matrix$BIC=NA
conf_matrix$LFO=NA

cn1<-list()
cn1.2<-list()
cn1.5<-list()
cn2<-list()
cn2.2<-list()
cn2.5<-list()
cn3<-list()

aic_set<-list()
bic_set<-list()
lfo_set<-list()
o=0

#rwa
for(a in seq_along(scn)){
  
  #AIC
  aica<-subset(aic,scenario==scn[a])
  aic_set[[a]]=tidyr::spread(aica[,-c(10,11,12,13)],key=model,value=mode)
  naicsim<-length(unique(aic_set[[a]]$iteration))
  aic_set[[a]]=aic_set[[a]][c(15,8,12)] #reorder estimation models
  #head(aic_set[[a]])
  
  #set1 = min. AIC
  sc1=apply(aic_set[[a]],1,which.min)
  
  #set2 = static parsimony - delta 2 AIC
  #delta AIC
  wAIC=t(apply(aic_set[[a]],1,samEst::model_weights,form='AIC'))
  dAIC=  t(apply(aic_set[[a]],1,function(x) x-min(x)))
  sc1.2=sc1
  sc1.2=ifelse(sc1.2==2&dAIC[,1]<2,1,sc1.2)
  sc1.2=ifelse(sc1.2==3&dAIC[,1]<2,1,sc1.2)
  
  #set3 = static parsimony - delta 5 AIC
  sc1.5=sc1
  sc1.5=ifelse(sc1.5==2&dAIC[,1]<5,1,sc1.5)
  sc1.5=ifelse(sc1.5==3&dAIC[,1]<5,1,sc1.5)
  
  cn1[[a]]=summary(factor(sc1,levels=seq(1:ncol(aic_set[[a]]))))/naicsim
  cn1.2[[a]]=summary(factor(sc1.2,levels=seq(1:ncol(aic_set[[a]]))))/naicsim
  cn1.5[[a]]=summary(factor(sc1.5,levels=seq(1:ncol(aic_set[[a]]))))/naicsim
  
  bica=subset(bic,scenario==scn[a])
  bic_set[[a]]=tidyr::spread(bica[,-c(10,11,12,13)],key=model,value=mode)
  nbicsim<-length(unique( bic_set[[a]]$iteration))
  bic_set[[a]]=bic_set[[a]][c(15,8,12)] #reorder estimation models
  
  sc2=apply(bic_set[[a]],1,which.min)
  
  #set2 = static parsimony - delta 2 BIC
  #delta AIC
  wBIC=t(apply(bic_set[[a]],1,samEst::model_weights,form='AIC'))
  dBIC=  t(apply(bic_set[[a]],1,function(x) x-min(x)))
  sc2.2=sc2
  sc2.2=ifelse(sc1.2==2&dBIC[,1]<2,1,sc1.2)
  sc2.2=ifelse(sc1.2==3&dBIC[,1]<2,1,sc1.2)
  
  #set3 = static parsimony - delta 5 BIC
  sc2.5=sc2
  sc2.5=ifelse(sc2.5==2&dBIC[,1]<5,1,sc2.5)
  sc2.5=ifelse(sc2.5==3&dBIC[,1]<5,1,sc2.5)
  
  cn2[[a]]=summary(factor(sc2,levels=seq(1:ncol(bic_set[[a]]))))/nbicsim
  cn2.2[[a]]=summary(factor(sc2.2,levels=seq(1:ncol(bic_set[[a]]))))/nbicsim
  cn2.5[[a]]=summary(factor(sc2.5,levels=seq(1:ncol(bic_set[[a]]))))/nbicsim
  
  lfoa=subset(lfo,scenario==scn[a])
  lfo_set[[a]]=tidyr::spread(lfoa[,-c(10,11,12,13)],key=model,value=mode)
  nlfosim <- length(unique( lfo_set[[a]]$iteration))
  lfo_set[[a]]=lfo_set[[a]][c(10,8,9)] #reorder estimation models
  #head(lfo_set[[a]])
  
  sc3=apply(lfo_set[[a]],1,which.max)
  cn3[[a]]=summary(factor(sc3,levels=seq(1:ncol(lfo_set[[a]]))))/ nlfosim
  wLFO=t(apply(lfo_set[[a]],1,samEst::model_weights,form='AIC'))
  
  #use weight thresholds - 0.55 ~ 2 dAIC
  # sc3.2=sc3
  #  sc3.2=ifelse(na.omit(wLFO[,2])<0.55,1,sc3.2)
  #  sc3.2=ifelse(sc3.2==3&wLFO[,3]<0.55,1,sc3.2)
  
  
  myseq<-seq(from=o+1, length.out=length(EM))
  conf_matrix$w_AIC[myseq]<-cn1[[a]]
  conf_matrix$w_AIC.d2[myseq]<-cn1.2[[a]]
  conf_matrix$w_AIC.d5[myseq]<-cn1.5[[a]]
  conf_matrix$BIC[myseq]<-cn2[[a]]
  conf_matrix$BIC.d2[myseq]<-cn2.2[[a]]
  conf_matrix$BIC.d5[myseq]<-cn2.5[[a]]
  conf_matrix$LFO[myseq]<-cn3[[a]]
  o=max(myseq)
}


paic=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = w_AIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_AIC,2)), vjust = 1,size=6) +
  ggtitle("AIC")+
  scale_fill_gradient(low="white", high="#009194")+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1)) +
  xlab("Simulation Scenario")+ylab("Estimation Model")
paic

paic2=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = w_AIC.d2), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_AIC.d2,2)), vjust = 1,size=6) +
  ggtitle("AIC - min. d2")+
  scale_fill_gradient(low="white", high="#009194")+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1)) +
  xlab("Simulation Scenario")+ylab("Estimation Model")
paic2

paic5=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = w_AIC.d5), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_AIC.d5,2)), vjust = 1,size=6) +
  ggtitle("AIC - min. d5")+
  scale_fill_gradient(low="white", high="#009194")+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1)) +
  xlab("Simulation Scenario")+ylab("Estimation Model")
paic5

pbic=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = BIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(BIC,2)), vjust = 1,size=6) +
  ggtitle("BIC")+
  scale_fill_gradient(low="white", high="#009194")  +
  
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
pbic

pbic2=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = BIC.d2), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(BIC.d2,2)), vjust = 1,size=6) +
  ggtitle("BIC - min. d2")+
  scale_fill_gradient(low="white", high="#009194")  +
  
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
pbic2

pbic3=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = BIC.d5), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(BIC.d5,2)), vjust = 1,size=6) +
  ggtitle("BIC - min. d5")+
  scale_fill_gradient(low="white", high="#009194")  +
  
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
pbic3


plfo=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = LFO), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(LFO,2)), vjust = 1, size=6) +
  ggtitle("LFO MLE")+
  scale_fill_gradient(low="white", high="#009194")  +
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
plfo

ggsave(here::here('outs','figures','aic_sensb.jpg')
       , plot=paic,width=10)
ggsave(here::here('outs','figures','aic_sensb_d2.jpg')
       , plot=paic2,width=10)
ggsave(here::here('outs','figures','aic_sensb_d5.jpg')
       , plot=paic3,width=10)
ggsave(here::here('outs','figures','bic_sensb.jpg')
       , plot=pbic,width=10)
ggsave(here::here('outs','figures','bic_sensb_d2.jpg')
       , plot=pbic2,width=10)
ggsave(here::here('outs','figures','bic_sensb_d5.jpg')
       , plot=pbic3,width=10)
ggsave(here::here('outs','figures','lfo_sensb.jpg')
       , plot=plfo,width=10)


#sigma low####
simPar <- read.csv("data/sigmalow_sensitivity/SimPars.csv")

## Store relevant object names to help run simulation 
scenNames <- unique(simPar$scenario)

res_siglow1<-readRDS(file = "outs/results/sigmalow_sensitivity/res_siglow1.rds")
res_siglow2<-readRDS(file = "outs/results/sigmalow_sensitivity/res_siglow2.rds")

ressiglow<-rbind(res_siglow1,res_siglow2)

#res_a <- res_a[res_a$convergence==0,]

aic=subset(ressiglow, parameter=='AIC'&method=='MLE')
bic=subset(ressiglow, parameter=='BIC'&method=='MLE')
lfo=subset(ressiglow, parameter=='LFO'&method=='MLE')


lfo<-lfo[lfo$model %in% c("simple","autocorr","rwa","rwb","rwab","hmma","hmmb","hmmab"),]
lfo[is.na(lfo$est),]<--Inf

lfo<-lfo[lfo$model %in% c("simple","autocorr","rwa","rwb","rwab","hmma","hmmb","hmmab"),]

lfo[is.na(lfo$est),]<--Inf
aic$est[aic$convergence>0]<-Inf
bic$est[aic$convergence>0]<-Inf

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

conf_matrix<-expand.grid(EM=EM,OM=scn)
conf_matrix$w_AIC=NA
conf_matrix$BIC=NA
conf_matrix$LFO=NA

cn1<-list()
cn1.2<-list()
cn1.5<-list()
cn2<-list()
cn2.2<-list()
cn2.5<-list()
cn3<-list()

aic_set<-list()
bic_set<-list()
lfo_set<-list()
o=0

#rwa
for(a in seq_along(scn)){
  
  #AIC
  aica<-subset(aic,scenario==scn[a])
  aic_set[[a]]=tidyr::spread(aica[,-c(10,11,12,13)],key=model,value=mode)
  naicsim<-length(unique(aic_set[[a]]$iteration))
  aic_set[[a]]=aic_set[[a]][c(15,8,12)] #reorder estimation models
  #head(aic_set[[a]])
  
  #set1 = min. AIC
  sc1=apply(aic_set[[a]],1,which.min)
  
  #set2 = static parsimony - delta 2 AIC
  #delta AIC
  wAIC=t(apply(aic_set[[a]],1,samEst::model_weights,form='AIC'))
  dAIC=  t(apply(aic_set[[a]],1,function(x) x-min(x)))
  sc1.2=sc1
  sc1.2=ifelse(sc1.2==2&dAIC[,1]<2,1,sc1.2)
  sc1.2=ifelse(sc1.2==3&dAIC[,1]<2,1,sc1.2)
  
  #set3 = static parsimony - delta 5 AIC
  sc1.5=sc1
  sc1.5=ifelse(sc1.5==2&dAIC[,1]<5,1,sc1.5)
  sc1.5=ifelse(sc1.5==3&dAIC[,1]<5,1,sc1.5)
  
  cn1[[a]]=summary(factor(sc1,levels=seq(1:ncol(aic_set[[a]]))))/naicsim
  cn1.2[[a]]=summary(factor(sc1.2,levels=seq(1:ncol(aic_set[[a]]))))/naicsim
  cn1.5[[a]]=summary(factor(sc1.5,levels=seq(1:ncol(aic_set[[a]]))))/naicsim
  
  bica=subset(bic,scenario==scn[a])
  bic_set[[a]]=tidyr::spread(bica[,-c(10,11,12,13)],key=model,value=mode)
  nbicsim<-length(unique( bic_set[[a]]$iteration))
  bic_set[[a]]=bic_set[[a]][c(15,8,12)] #reorder estimation models
  
  sc2=apply(bic_set[[a]],1,which.min)
  
  #set2 = static parsimony - delta 2 BIC
  #delta AIC
  wBIC=t(apply(bic_set[[a]],1,samEst::model_weights,form='AIC'))
  dBIC=  t(apply(bic_set[[a]],1,function(x) x-min(x)))
  sc2.2=sc2
  sc2.2=ifelse(sc1.2==2&dBIC[,1]<2,1,sc1.2)
  sc2.2=ifelse(sc1.2==3&dBIC[,1]<2,1,sc1.2)
  
  #set3 = static parsimony - delta 5 BIC
  sc2.5=sc2
  sc2.5=ifelse(sc2.5==2&dBIC[,1]<5,1,sc2.5)
  sc2.5=ifelse(sc2.5==3&dBIC[,1]<5,1,sc2.5)
  
  cn2[[a]]=summary(factor(sc2,levels=seq(1:ncol(bic_set[[a]]))))/nbicsim
  cn2.2[[a]]=summary(factor(sc2.2,levels=seq(1:ncol(bic_set[[a]]))))/nbicsim
  cn2.5[[a]]=summary(factor(sc2.5,levels=seq(1:ncol(bic_set[[a]]))))/nbicsim
  
  lfoa=subset(lfo,scenario==scn[a])
  lfo_set[[a]]=tidyr::spread(lfoa[,-c(10,11,12,13)],key=model,value=mode)
  nlfosim <- length(unique( lfo_set[[a]]$iteration))
  lfo_set[[a]]=lfo_set[[a]][c(10,8,9)] #reorder estimation models
  #head(lfo_set[[a]])
  
  sc3=apply(lfo_set[[a]],1,which.max)
  cn3[[a]]=summary(factor(sc3,levels=seq(1:ncol(lfo_set[[a]]))))/ nlfosim
  wLFO=t(apply(lfo_set[[a]],1,samEst::model_weights,form='AIC'))
  
  #use weight thresholds - 0.55 ~ 2 dAIC
  # sc3.2=sc3
  #  sc3.2=ifelse(na.omit(wLFO[,2])<0.55,1,sc3.2)
  #  sc3.2=ifelse(sc3.2==3&wLFO[,3]<0.55,1,sc3.2)
  
  
  myseq<-seq(from=o+1, length.out=length(EM))
  conf_matrix$w_AIC[myseq]<-cn1[[a]]
  conf_matrix$w_AIC.d2[myseq]<-cn1.2[[a]]
  conf_matrix$w_AIC.d5[myseq]<-cn1.5[[a]]
  conf_matrix$BIC[myseq]<-cn2[[a]]
  conf_matrix$BIC.d2[myseq]<-cn2.2[[a]]
  conf_matrix$BIC.d5[myseq]<-cn2.5[[a]]
  conf_matrix$LFO[myseq]<-cn3[[a]]
  o=max(myseq)
}


paic=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = w_AIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_AIC,2)), vjust = 1,size=6) +
  ggtitle("AIC")+
  scale_fill_gradient(low="white", high="#009194")+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1)) +
  xlab("Simulation Scenario")+ylab("Estimation Model")
paic

paic2=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = w_AIC.d2), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_AIC.d2,2)), vjust = 1,size=6) +
  ggtitle("AIC - min. d2")+
  scale_fill_gradient(low="white", high="#009194")+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1)) +
  xlab("Simulation Scenario")+ylab("Estimation Model")
paic2

paic5=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = w_AIC.d5), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_AIC.d5,2)), vjust = 1,size=6) +
  ggtitle("AIC - min. d5")+
  scale_fill_gradient(low="white", high="#009194")+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1)) +
  xlab("Simulation Scenario")+ylab("Estimation Model")
paic5

pbic=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = BIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(BIC,2)), vjust = 1,size=6) +
  ggtitle("BIC")+
  scale_fill_gradient(low="white", high="#009194")  +
  
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
pbic

pbic2=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = BIC.d2), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(BIC.d2,2)), vjust = 1,size=6) +
  ggtitle("BIC - min. d2")+
  scale_fill_gradient(low="white", high="#009194")  +
  
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
pbic2

pbic3=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = BIC.d5), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(BIC.d5,2)), vjust = 1,size=6) +
  ggtitle("BIC - min. d5")+
  scale_fill_gradient(low="white", high="#009194")  +
  
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
pbic3


plfo=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = LFO), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(LFO,2)), vjust = 1, size=6) +
  ggtitle("LFO MLE")+
  scale_fill_gradient(low="white", high="#009194")  +
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
plfo

ggsave(here::here('outs','figures','aic_siglow.jpg')
       , plot=paic,width=10)
ggsave(here::here('outs','figures','aic_siglow_d2.jpg')
       , plot=paic2,width=10)
ggsave(here::here('outs','figures','aic_siglow_d5.jpg')
       , plot=paic3,width=10)
ggsave(here::here('outs','figures','bic_siglow.jpg')
       , plot=pbic,width=10)
ggsave(here::here('outs','figures','bic_siglow_d2.jpg')
       , plot=pbic2,width=10)
ggsave(here::here('outs','figures','bic_siglow_d5.jpg')
       , plot=pbic3,width=10)
ggsave(here::here('outs','figures','lfo_siglow.jpg')
       , plot=plfo,width=10)

#sigma med####
simPar <- read.csv("data/sigmamed_sensitivity/SimPars.csv")

## Store relevant object names to help run simulation 
scenNames <- unique(simPar$scenario)

res_sigmed1<-readRDS(file = "outs/results/sigmamed_sensitivity/res_sigmed1.rds")
res_sigmed2<-readRDS(file = "outs/results/sigmamed_sensitivity/res_sigmed2.rds")

ressigmed<-rbind(res_sigmed1,res_sigmed2)

#res_a <- res_a[res_a$convergence==0,]

aic=subset(ressigmed, parameter=='AIC'&method=='MLE')
bic=subset(ressigmed, parameter=='BIC'&method=='MLE')
lfo=subset(ressigmed, parameter=='LFO'&method=='MLE')


lfo<-lfo[lfo$model %in% c("simple","autocorr","rwa","rwb","rwab","hmma","hmmb","hmmab"),]
lfo[is.na(lfo$est),]<--Inf

lfo<-lfo[lfo$model %in% c("simple","autocorr","rwa","rwb","rwab","hmma","hmmb","hmmab"),]

lfo[is.na(lfo$est),]<--Inf
aic$est[aic$convergence>0]<-Inf
bic$est[aic$convergence>0]<-Inf

scn<-factor(unique(aic$scenario), levels=c(
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

conf_matrix<-expand.grid(EM=EM,OM=scn)
conf_matrix$w_AIC=NA
conf_matrix$BIC=NA
conf_matrix$LFO=NA

cn1<-list()
cn1.2<-list()
cn1.5<-list()
cn2<-list()
cn2.2<-list()
cn2.5<-list()
cn3<-list()

aic_set<-list()
bic_set<-list()
lfo_set<-list()
o=0

#rwa
for(a in seq_along(scn)){
  
  #AIC
  aica<-subset(aic,scenario==scn[a])
  aic_set[[a]]=tidyr::spread(aica[,-c(10,11,12,13)],key=model,value=mode)
  naicsim<-length(unique(aic_set[[a]]$iteration))
  aic_set[[a]]=aic_set[[a]][c(15,8,12)] #reorder estimation models
  #head(aic_set[[a]])
  
  #set1 = min. AIC
  sc1=apply(aic_set[[a]],1,which.min)
  
  #set2 = static parsimony - delta 2 AIC
  #delta AIC
  wAIC=t(apply(aic_set[[a]],1,samEst::model_weights,form='AIC'))
  dAIC=  t(apply(aic_set[[a]],1,function(x) x-min(x)))
  sc1.2=sc1
  sc1.2=ifelse(sc1.2==2&dAIC[,1]<2,1,sc1.2)
  sc1.2=ifelse(sc1.2==3&dAIC[,1]<2,1,sc1.2)
  
  #set3 = static parsimony - delta 5 AIC
  sc1.5=sc1
  sc1.5=ifelse(sc1.5==2&dAIC[,1]<5,1,sc1.5)
  sc1.5=ifelse(sc1.5==3&dAIC[,1]<5,1,sc1.5)
  
  cn1[[a]]=summary(factor(sc1,levels=seq(1:ncol(aic_set[[a]]))))/naicsim
  cn1.2[[a]]=summary(factor(sc1.2,levels=seq(1:ncol(aic_set[[a]]))))/naicsim
  cn1.5[[a]]=summary(factor(sc1.5,levels=seq(1:ncol(aic_set[[a]]))))/naicsim
  
  bica=subset(bic,scenario==scn[a])
  bic_set[[a]]=tidyr::spread(bica[,-c(10,11,12,13)],key=model,value=mode)
  nbicsim<-length(unique( bic_set[[a]]$iteration))
  bic_set[[a]]=bic_set[[a]][c(15,8,12)] #reorder estimation models
  
  sc2=apply(bic_set[[a]],1,which.min)
  
  #set2 = static parsimony - delta 2 BIC
  #delta AIC
  wBIC=t(apply(bic_set[[a]],1,samEst::model_weights,form='AIC'))
  dBIC=  t(apply(bic_set[[a]],1,function(x) x-min(x)))
  sc2.2=sc2
  sc2.2=ifelse(sc1.2==2&dBIC[,1]<2,1,sc1.2)
  sc2.2=ifelse(sc1.2==3&dBIC[,1]<2,1,sc1.2)
  
  #set3 = static parsimony - delta 5 BIC
  sc2.5=sc2
  sc2.5=ifelse(sc2.5==2&dBIC[,1]<5,1,sc2.5)
  sc2.5=ifelse(sc2.5==3&dBIC[,1]<5,1,sc2.5)
  
  cn2[[a]]=summary(factor(sc2,levels=seq(1:ncol(bic_set[[a]]))))/nbicsim
  cn2.2[[a]]=summary(factor(sc2.2,levels=seq(1:ncol(bic_set[[a]]))))/nbicsim
  cn2.5[[a]]=summary(factor(sc2.5,levels=seq(1:ncol(bic_set[[a]]))))/nbicsim
  
  lfoa=subset(lfo,scenario==scn[a])
  lfo_set[[a]]=tidyr::spread(lfoa[,-c(10,11,12,13)],key=model,value=mode)
  nlfosim <- length(unique( lfo_set[[a]]$iteration))
  lfo_set[[a]]=lfo_set[[a]][c(10,8,9)] #reorder estimation models
  #head(lfo_set[[a]])
  
  sc3=apply(lfo_set[[a]],1,which.max)
  cn3[[a]]=summary(factor(sc3,levels=seq(1:ncol(lfo_set[[a]]))))/ nlfosim
  wLFO=t(apply(lfo_set[[a]],1,samEst::model_weights,form='AIC'))
  
  #use weight thresholds - 0.55 ~ 2 dAIC
  # sc3.2=sc3
  #  sc3.2=ifelse(na.omit(wLFO[,2])<0.55,1,sc3.2)
  #  sc3.2=ifelse(sc3.2==3&wLFO[,3]<0.55,1,sc3.2)
  
  
  myseq<-seq(from=o+1, length.out=length(EM))
  conf_matrix$w_AIC[myseq]<-cn1[[a]]
  conf_matrix$w_AIC.d2[myseq]<-cn1.2[[a]]
  conf_matrix$w_AIC.d5[myseq]<-cn1.5[[a]]
  conf_matrix$BIC[myseq]<-cn2[[a]]
  conf_matrix$BIC.d2[myseq]<-cn2.2[[a]]
  conf_matrix$BIC.d5[myseq]<-cn2.5[[a]]
  conf_matrix$LFO[myseq]<-cn3[[a]]
  o=max(myseq)
}


paic=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = w_AIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_AIC,2)), vjust = 1,size=6) +
  ggtitle("AIC")+
  scale_fill_gradient(low="white", high="#009194")+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1)) +
  xlab("Simulation Scenario")+ylab("Estimation Model")
paic

paic2=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = w_AIC.d2), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_AIC.d2,2)), vjust = 1,size=6) +
  ggtitle("AIC - min. d2")+
  scale_fill_gradient(low="white", high="#009194")+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1)) +
  xlab("Simulation Scenario")+ylab("Estimation Model")
paic2

paic5=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = w_AIC.d5), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_AIC.d5,2)), vjust = 1,size=6) +
  ggtitle("AIC - min. d5")+
  scale_fill_gradient(low="white", high="#009194")+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1)) +
  xlab("Simulation Scenario")+ylab("Estimation Model")
paic5

pbic=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = BIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(BIC,2)), vjust = 1,size=6) +
  ggtitle("BIC")+
  scale_fill_gradient(low="white", high="#009194")  +
  
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
pbic

pbic2=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = BIC.d2), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(BIC.d2,2)), vjust = 1,size=6) +
  ggtitle("BIC - min. d2")+
  scale_fill_gradient(low="white", high="#009194")  +
  
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
pbic2

pbic3=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = BIC.d5), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(BIC.d5,2)), vjust = 1,size=6) +
  ggtitle("BIC - min. d5")+
  scale_fill_gradient(low="white", high="#009194")  +
  
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
pbic3


plfo=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = LFO), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(LFO,2)), vjust = 1, size=6) +
  ggtitle("LFO MLE")+
  scale_fill_gradient(low="white", high="#009194")  +
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
plfo

ggsave(here::here('outs','figures','aic_sigmed.jpg')
       , plot=paic,width=10)
ggsave(here::here('outs','figures','aic_sigmed_d2.jpg')
       , plot=paic2,width=10)
ggsave(here::here('outs','figures','aic_sigmed_d5.jpg')
       , plot=paic3,width=10)
ggsave(here::here('outs','figures','bic_sigmed.jpg')
       , plot=pbic,width=10)
ggsave(here::here('outs','figures','bic_sigmed_d2.jpg')
       , plot=pbic2,width=10)
ggsave(here::here('outs','figures','bic_sigmed_d5.jpg')
       , plot=pbic3,width=10)
ggsave(here::here('outs','figures','lfo_sigmed.jpg')
       , plot=plfo,width=10)
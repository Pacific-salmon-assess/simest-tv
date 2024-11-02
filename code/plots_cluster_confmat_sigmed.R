#routines to visualize smulation estimation Results
#Catarina Wor
#November 2022
#============================================



library(gridExtra)
library(ggplot2)
library(dplyr)
library(cowplot)
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
simPar <- read.csv("data/sigmamed_sensitivity/SimPars.csv")

## Store relevant object names to help run simulation 
scenNames <- unique(simPar$scenario)

res<-readRDS(file = "outs/simest/sigmamed_sensitivity/res_sigmed.rds")


res<-res[res$convergence==0,]

res$sigscenario<-"medium"
res$scenario<-case_match(
res$scenario,
 "sigmamed_stationary"~"stationary",
 "sigmamed_decLinearProd"~"decLinearProd",        
 "sigmamed_regimeProd"~ "regimeProd",            
 "sigmamed_sineProd"~ "sineProd",             
 "sigmamed_regimeCap"~ "regimeCap",             
 "sigmamed_decLinearCap"~ "decLinearCap",        
 "sigmamed_regimeProdCap"~ "regimeProdCap",
 "sigmamed_shiftCap"~"shiftCap",
 "sigmamed_decLinearProdshiftCap"~"decLinearProdshiftCap")

aic=subset(res,parameter=='AIC'&method=='MLE')
bic=subset(res,parameter=='BIC'&method=='MLE')
lfo=subset(res,parameter=='LFO'&method=='MLE')

lfo<-lfo[lfo$model %in% c("simple","autocorr","rwa","rwb","rwab",
    "hmma","hmmb","hmmab"),]


scn<-factor(unique(aic$scenario), levels=c(
   "stationary",
   "decLinearProd",
   "sineProd",
   "regimeProd",
   "decLinearCap",         
   "regimeCap",
   "shiftCap",
   "decLinearProdshiftCap",
   "regimeProdCap"
 ) )


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
  aic_set[[a]]=aic_set[[a]][c(16,9,13,10,15,12,14,11)] #reorder estimation models
  #head(aic_set[[a]])
  sc1=apply(aic_set[[a]],1,which.min)
  
  cn1[[a]]=summary(factor(sc1,levels=seq(1:ncol(aic_set[[a]]))))/naicsim

  bica=subset(bic,scenario==scn[a])
  bic_set[[a]]=tidyr::spread(bica[,-c(10,11,12,13)],key=model,value=mode)
  nbicsim<-length(unique( bic_set[[a]]$iteration))
  bic_set[[a]]=bic_set[[a]][c(16,9,13,10,15,12,14,11)] #reorder estimation models

  sc2=apply(bic_set[[a]],1,which.min)
 
  cn2[[a]]=summary(factor(sc2,levels=seq(1:ncol(bic_set[[a]]))))/nbicsim

  lfoa=subset(lfo,scenario==scn[a])
  lfo_set[[a]]=tidyr::spread(lfoa[,-c(10,11,12,13)],key=model,value=mode)
  nlfosim <- length(unique( lfo_set[[a]]$iteration))
  lfo_set[[a]]=lfo_set[[a]][c(16,9,13,10,15,12,14,11)] #reorder estimation models
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
 conf_matrix$eqem_om<-factor( conf_matrix$eqem_om,levels=c( "stationary",
 "autocorr", "dynamic.a", "regime.a", "dynamic.b", "regime.b", "dynamic.ab", "regime.ab"
  ))      

conf_matrix$diag<-conf_matrix$eqem_om==conf_matrix$EM


conf_matrix$scentype<-dplyr::case_match(as.character(conf_matrix$OM), 
      "stationary"~"stationary",
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

conf_matrix$scendesc<-dplyr::case_match(as.character(conf_matrix$OM), 
      "stationary"~"",
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
 
conf_matrix$fullname <- apply( conf_matrix[ , c("scendesc", "scentype") ] ,
                                   1 , paste , collapse = " " )


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

paic_sigmed=ggplot(data =  conf_matrix, mapping = aes(x = fullname, y = EM)) +
  geom_tile(aes(fill = w_AIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_AIC,2)), vjust = 1,size=6) +
  ggtitle(expression(paste("AIC", ~sigma, "=0.3")))+
  scale_fill_gradient(low="white", high="#009194")  +
  scale_x_discrete(labels = ~ stringr::str_wrap(as.character(.x), 15))+
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray90", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1)) +
  xlab("Simulation Scenario")+ylab("Estimation Model")
paic_sigmed

ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/confusion_matrices/sig_med/AIC_MLE_sigmed.png",
 plot=paic_sigmed, width = 9,height = 7)



pbic_sigmed=ggplot(data =  conf_matrix, mapping = aes(x = fullname, y = EM)) +
  geom_tile(aes(fill = BIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(BIC,2)), vjust = 1,size=6) +
  ggtitle(expression(paste("BIC", ~sigma, "=0.3"))) +
  scale_fill_gradient(low="white", high="#009194") +
  scale_x_discrete(labels = ~ stringr::str_wrap(as.character(.x), 15))+
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray90", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
pbic_sigmed

ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/confusion_matrices/sig_med/BIC_MLE_sigmed.png", 
  plot=pbic_sigmed, width = 9,height = 7)



plfo=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = LFO), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(LFO,2)), vjust = 1, size=6) +
   ggtitle(expression(paste("LFO MLE", ~sigma, "=0.3")))+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray90", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
plfo

ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/confusion_matrices/sig_med/LFO_MLE_sigmed.png", 
  plot=plfo, width = 8,height = 7)







#===============================================================================

ressiglow<-readRDS(file = "outs/simest/sigmalow_sensitivity/res_siglow.rds")
#res2siglow<-readRDS(file = "outs/simest/sigmalow_sensitivity/res_siglow2.rds")
#res105siglow<-readRDS(file = "outs/simest/sigmalow_sensitivity/res_siglow_105.rds")
#ressiglow<-rbind(res1siglow,res2siglow,res105siglow)#,resstan16,resstan712)

unique(ressiglow$scenario)

ressiglow<-ressiglow[ressiglow$convergence==0,]


ressiglow$sigscenario<-"low"
ressiglow$scenario<-case_match(
  ressiglow$scenario,
  "sigmalow_stationary"~"stationary",
  "sigmalow_decLinearProd"~"decLinearProd",
  "sigmalow_regimeProd"~ "regimeProd",  
 "sigmalow_sineProd"~ "sineProd",     
 "sigmalow_regimeCap"~ "regimeCap",    
 "sigmalow_decLinearCap"~ "decLinearCap", 
 "sigmalow_regimeProdCap"~ "regimeProdCap",
 "sigmalow_shiftCap"~"shiftCap",
 "sigmalow_decLinearProdshiftCap"~"decLinearProdshiftCap")


aic_siglow=subset(ressiglow,parameter=='AIC'&method=='MLE')
bic_siglow=subset(ressiglow,parameter=='BIC'&method=='MLE')
lfo_siglow=subset(ressiglow,parameter=='LFO'&method=='MLE')

lfo_siglow<-lfo_siglow[lfo_siglow$model %in% c("simple","autocorr","rwa","rwb","rwab",
    "hmma","hmmb","hmmab"),]


scn_siglow<-factor(unique(aic_siglow$scenario), levels=c(
  "stationary", 
  "decLinearProd",  
  "sineProd",
  "regimeProd", 
  "decLinearCap",
  "regimeCap", 
  "shiftCap",
  "regimeProdCap",
  "decLinearProdshiftCap" ) )


##Confusion matrices
conf_matrix_siglow<-expand.grid(EM=EM,OM=scn_siglow)
conf_matrix_siglow$w_AIC=NA
conf_matrix_siglow$BIC=NA
conf_matrix_siglow$LFO=NA

cn1_siglow<-list()
cn2_siglow<-list()
cn3_siglow<-list()

aic_set_siglow<-list()
bic_set_siglow<-list()
lfo_set_siglow<-list()
o=0
for(a in seq_along(scn_siglow)){

  #head(aic_set[[a]][c(15,8,12,9,14,11,13,10)])
  #AIC
  aica_siglow<-subset(aic_siglow,scenario==scn_siglow[a])
  aic_set_siglow[[a]]=tidyr::spread(aica_siglow[,-c(10,11,12,13)],key=model,value=mode)
  naicsim<-length(unique(aic_set_siglow[[a]]$iteration))
  aic_set_siglow[[a]]=aic_set_siglow[[a]][c(16,9,13,10,15,12,14,11)] #reorder estimation models
  #head(aic_set[[a]])
  sc1_siglow=apply(aic_set_siglow[[a]],1,which.min)
  
  cn1_siglow[[a]]=summary(factor(sc1_siglow,levels=seq(1:ncol(aic_set_siglow[[a]]))))/naicsim

  bica_siglow=subset(bic_siglow,scenario==scn_siglow[a])
  bic_set_siglow[[a]]=tidyr::spread(bica_siglow[,-c(10,11,12,13)],key=model,value=mode)
  nbicsim<-length(unique( bic_set_siglow[[a]]$iteration))
  bic_set_siglow[[a]]=bic_set_siglow[[a]][c(16,9,13,10,15,12,14,11)] #reorder estimation models

  sc2_siglow=apply(bic_set_siglow[[a]],1,which.min)
 
  cn2_siglow[[a]]=summary(factor(sc2_siglow,levels=seq(1:ncol(bic_set_siglow[[a]]))))/nbicsim

  lfoa_siglow=subset(lfo_siglow,scenario==scn_siglow[a])
  lfo_set_siglow[[a]]=tidyr::spread(lfoa_siglow[,-c(10,11,12,13)],key=model,value=mode)
  nlfosim <- length(unique(lfo_set_siglow[[a]]$iteration))
  lfo_set_siglow[[a]]=lfo_set_siglow[[a]][c(16,9,13,10,15,12,14,11)] #reorder estimation models
  #head(lfo_set[[a]])

  sc3_siglow=apply(lfo_set_siglow[[a]],1,which.max)
  cn3_siglow[[a]]=summary(factor(sc3_siglow,levels=seq(1:ncol(lfo_set_siglow[[a]]))))/ nlfosim


  myseq<-seq(from=o+1, length.out=length(EM))
  conf_matrix_siglow$w_AIC[myseq]<-cn1_siglow[[a]]
  conf_matrix_siglow$BIC[myseq]<-cn2_siglow[[a]]
  conf_matrix_siglow$LFO[myseq]<-cn3_siglow[[a]]
  o=max(myseq)

}


conf_matrix_siglow$eqem_om <- dplyr::recode(conf_matrix_siglow$OM, 
      "stationary"="stationary",  
      "decLinearProd"="dynamic.a",
      "sineProd"="dynamic.a",
      "decLinearCap"="dynamic.b",
      "regimeProd"="regime.a",
      "regimeCap"="regime.b",
      "shiftCap"="regime.b", 
      "regimeProdCap"="regime.ab",
      "decLinearProdshiftCap"="dynamic.ab"
      ) 

      conf_matrix_siglow$eqem_om<-factor(conf_matrix_siglow$eqem_om,levels=c("stationary",
     "autocorr",
     "dynamic.a","regime.a",
     "dynamic.b","regime.b",
     "dynamic.ab",
     "regime.ab"))  
conf_matrix_siglow$diag<-conf_matrix_siglow$eqem_om==conf_matrix_siglow$EM





#========================================================================================================
#base case - sigma=0.6

resbase<-readRDS(file = "outs/simest/generic/resbase.rds")
#res2base<-readRDS(file = "outs/simest/generic/resbase2.rds")

#resbase<-rbind(res1base,res2base)#,resstan16,resstan712)

resbase<-resbase[resbase$convergence==0,]

resbase$sigscenario<-"high"

unique(resbase$scenario)[unique(resbase$scenario)%in%c("stationary",
                                       "decLinearProd",
                                       "sineProd",
                                       "regimeProd",
                                       "decLinearCap",         
                                       "regimeCap",
                                       "shiftCap",
                                       "decLinearProdshiftCap",
                                       "regimeProdCap")]

resbase<-resbase[resbase$scenario%in%c("stationary",
                                       "decLinearProd",
                                       "sineProd",
                                       "regimeProd",
                                       "decLinearCap",         
                                       "regimeCap",
                                       "shiftCap",
                                       "decLinearProdshiftCap",
                                       "regimeProdCap"),]


aicbase=subset(resbase,parameter=='AIC'&method=='MLE')
bicbase=subset(resbase,parameter=='BIC'&method=='MLE')
lfobase=subset(resbase,parameter=='LFO'&method=='MLE')

lfobase<-lfobase[lfobase$model %in% c("simple","autocorr","rwa","rwb","rwab",
    "hmma","hmmb","hmmab"),]


##Confusion matrices
conf_matrix_base<-expand.grid(EM=EM,OM=scn)
conf_matrix_base$w_AIC=NA
conf_matrix_base$BIC=NA
conf_matrix_base$LFO=NA

cn1_base<-list()
cn2_base<-list()
cn3_base<-list()

aic_set_base<-list()
bic_set_base<-list()
lfo_set_base<-list()
o=0
for(a in seq_along(scn)){

  #AIC
  aica_base<-subset(aicbase,scenario==scn[a])
  aic_set_base[[a]]=tidyr::spread(aica_base[,-c(10,11,12,13)],key=model,value=mode)
  naicsim<-length(unique(aic_set_base[[a]]$iteration))
  aic_set_base[[a]]=aic_set_base[[a]][c(16,9,13,10,15,12,14,11)] #reorder estimation models
  #head(aic_set[[a]])
  sc1_base=apply(aic_set_base[[a]],1,which.min)
  
  cn1_base[[a]]=summary(factor(sc1_base,levels=seq(1:ncol(aic_set_base[[a]]))))/naicsim

  bica_base=subset(bicbase,scenario==scn[a])
  bic_set_base[[a]]=tidyr::spread(bica_base[,-c(10,11,12,13)],key=model,value=mode)
  nbicsim<-length(unique( bic_set_base[[a]]$iteration))
  bic_set_base[[a]]=bic_set_base[[a]][c(16,9,13,10,15,12,14,11)] #reorder estimation models

  sc2_base=apply(bic_set_base[[a]],1,which.min)
 
  cn2_base[[a]]=summary(factor(sc2_base,levels=seq(1:ncol(bic_set_base[[a]]))))/nbicsim

  lfoa_base=subset(lfobase,scenario==scn[a])
  lfo_set_base[[a]]=tidyr::spread(lfoa_base[,-c(10,11,12,13)],key=model,value=mode)
  nlfosim <- length(unique( lfo_set_base[[a]]$iteration))
  lfo_set_base[[a]]=lfo_set_base[[a]][c(16,9,13,10,15,12,14,11)] #reorder estimation models
  #head(lfo_set[[a]])

  sc3_base=apply(lfo_set_base[[a]],1,which.max)
  cn3_base[[a]]=summary(factor(sc3_base,levels=seq(1:ncol(lfo_set_base[[a]]))))/ nlfosim


  myseq<-seq(from=o+1, length.out=length(EM))
  conf_matrix_base$w_AIC[myseq]<-cn1_base[[a]]
  conf_matrix_base$BIC[myseq]<-cn2_base[[a]]
  conf_matrix_base$LFO[myseq]<-cn3_base[[a]]
  o=max(myseq)

}


conf_matrix_base$eqem_om <- dplyr::recode(conf_matrix_base$OM, 
      "stationary"="stationary", 
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
 conf_matrix_base$eqem_om<-factor( conf_matrix_base$eqem_om,levels=c( "stationary",
 "autocorr", "dynamic.a", "regime.a", "dynamic.b", "regime.b", "dynamic.ab", "regime.ab"
  ))      

conf_matrix_base$diag<-conf_matrix_base$eqem_om==conf_matrix_base$EM




#======================================================================================================================
#aggregate plot

conf_matrix$sigma<-0.3
conf_matrix_base$sigma<-0.6
conf_matrix_siglow$sigma<-0.1


conf_matrix_sigcomp<-rbind(conf_matrix,conf_matrix_base)#,conf_matrix_siglow)
conf_matrix_sigcomp$OM

conf_matrix_sigcomp<-conf_matrix_sigcomp[conf_matrix_sigcomp$OM%in%c("decLinearProd","regimeProd","decLinearCap","shiftCap"), ]

conf_matrix_sigcomp$OM<-as.character(conf_matrix_sigcomp$OM)
conf_matrix_a$EM<-as.character(conf_matrix_a$EM)

conf_matrix_siglow[conf_matrix_siglow$OM=="regimeProd",]
conf_matrix[conf_matrix$OM=="regimeProd",]
conf_matrix_base[conf_matrix_base$OM=="regimeProd",]

conf_matrix_sigcomp_sc1<-conf_matrix_sigcomp[conf_matrix_sigcomp$OM=="decLinearProd"&conf_matrix_sigcomp$EM=="dynamic.a",]
conf_matrix_sigcomp_sc2<-conf_matrix_sigcomp[conf_matrix_sigcomp$OM=="regimeProd"&conf_matrix_sigcomp$EM%in%c("dynamic.a"),]
conf_matrix_sigcomp_sc3<-conf_matrix_sigcomp[conf_matrix_sigcomp$OM=="decLinearCap"&conf_matrix_sigcomp$EM=="dynamic.b",]
conf_matrix_sigcomp_sc4<-conf_matrix_sigcomp[conf_matrix_sigcomp$OM=="shiftCap"&conf_matrix_sigcomp$EM%in%c("dynamic.b"),]
#conf_matrix_sigcomp_sc5<-conf_matrix_sigcomp[conf_matrix_sigcomp$OM=="decLinearProdshiftCap"&conf_matrix_sigcomp$EM=="dynamic.ab",]



conf_matrix_sigcomp_right<-rbind(conf_matrix_sigcomp_sc1,
  conf_matrix_sigcomp_sc2,
  conf_matrix_sigcomp_sc3,
  conf_matrix_sigcomp_sc4)

unique(conf_matrix_sigcomp$OM)

conf_matrix_sigcomp_right$OM2<-case_match(
  conf_matrix_sigcomp_right$OM,
  "decLinearProd" ~"decline log(a)",
  "regimeProd" ~ "shift up log(a)",   
  "decLinearCap" ~"decline Smax",  
  "shiftCap"  ~"shift down Smax")
conf_matrix_sigcomp_right$OM2<-factor(conf_matrix_sigcomp_right$OM2,levels=c("shift up log(a)","decline log(a)",
  "shift down Smax", "decline Smax"))

head(conf_matrix_sigcomp)

lineAIC_sig<-ggplot(conf_matrix_sigcomp_right)+
geom_point(aes(x=sigma,y=w_AIC,color=OM2),size=4,show.legend = F)+
geom_line(aes(x=sigma,y=w_AIC,color=OM2),linewidth=1.2)+
#scale_colour_manual(values = c("#95D840FF","#482677FF", "#95D840FF", "#482677FF"))+
#scale_linetype_manual(values = c(1,1,2,2))+
ylab("proportion of correct model assignment with AICc")+
xlab(expression(paste("value of", ~sigma)))+
scale_color_viridis_d(begin=.1, end=.8,option = "A") +
#coord_cartesian(ylim=c(0.1,0.88))+
mytheme+
theme(axis.title=element_text(size=14,face="bold"),legend.key.width = unit(3, "line"))+
 guides(color=guide_legend(nrow=2, byrow=TRUE)) 
lineAIC_sig

unique( conf_matrix_sigcomp_right$OM2) 

barAIC_sig<-ggplot(conf_matrix_sigcomp_right)+
geom_bar(stat="identity",aes(x=as.factor(sigma),y=w_AIC,fill=OM2), position=position_dodge())+
#geom_line(aes(x=sigma,y=w_AIC,color=OM2),linewidth=1.2)+
#scale_colour_manual(values = c("#95D840FF","#482677FF", "#95D840FF", "#482677FF"))+
#scale_linetype_manual(values = c(1,1,2,2))+
ylab("proportion of correct model assignment with AICc")+
xlab(expression(paste("value of", ~sigma)))+
scale_fill_viridis_d(begin=.1, end=.8,option = "B", labels = c(expression(shift~up~log(alpha)) , 
                                                              expression(decline~log(alpha)),
                                                              expression(shift~down~S[max]),
                                                              expression(decline~S[max]))) +
#coord_cartesian(ylim=c(0.1,0.88))+
mytheme+
theme(axis.title=element_text(size=14,face="bold"),legend.key.width = unit(1.5, "line"))+
guides(fill=guide_legend(nrow=2, byrow=TRUE)) 
barAIC_sig



lineBIC_sig<-ggplot(conf_matrix_sigcomp_right)+
geom_point(aes(x=sigma,y=BIC,color=OM2),size=4,show.legend = F)+
geom_line(aes(x=sigma,y=BIC,color=OM2),linewidth=1.2)+
#scale_colour_manual(values = c("#95D840FF","#482677FF", "#95D840FF", "#482677FF"))+
#scale_linetype_manual(values = c(1,1,2,2))+
ylab("proportion of correct model assignment with BIC")+
xlab(expression(paste("value of", ~sigma)))+
scale_color_viridis_d(begin=.1, end=.8,option = "B") +
#coord_cartesian(ylim=c(0.1,0.88))+
mytheme+
theme(axis.title=element_text(size=14,face="bold"),legend.key.width = unit(1.5, "line"))+
guides(color=guide_legend(nrow=2, byrow=TRUE)) 
lineBIC_sig


barBIC_sig<-ggplot(conf_matrix_sigcomp_right)+
geom_bar(stat="identity",aes(x=as.factor(sigma),y=BIC,fill=OM2), position=position_dodge())+
#scale_colour_manual(values = c("#95D840FF","#482677FF", "#95D840FF", "#482677FF"))+
#scale_linetype_manual(values = c(1,1,2,2))+
ylab("proportion of correct model assignment with BIC")+
xlab(expression(paste("value of", ~sigma)))+
scale_fill_viridis_d(begin=.1, end=.9,option = "B", labels = c(expression(shift~up~log(alpha)) , 
                                                              expression(decline~log(alpha)),
                                                              expression(shift~down~S[max]),
                                                              expression(decline~S[max]))) +
mytheme +
theme(axis.title=element_text(size=14,face="bold"),legend.key.width = unit(1.5, "line")) +
guides(fill=guide_legend(nrow=2, byrow=TRUE)) 
barBIC_sig


#----------------------------------
#need to run plots_cluster_confmat_sensa before this works
#require that you run confmat_sens_a first
linesAIC_alpha_smax <- ggarrange(lineAIC_a, lineAICsmax, #lineAIC_sig,
                        nrow = 1, ncol = 2,
                        common.legend = TRUE,
                        legend="bottom")
linesAIC_alpha_smax




plot_grid(linesAIC_alpha_smax,lineAIC_sig,  rel_widths = c(2, 1) )


plot_grid(linesAIC_alpha_smax,barAIC_sig,  rel_widths = c(2, 1) )


ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/confusion_matrices/lineAICsens_alpha_smax.png",
 plot=linesAIC_alpha_smax, width = 13,height = 6)




linesBIC_alpha_smax <- ggarrange(lineBIC_a, lineBICsmax, 
                        nrow = 1, ncol = 2,
                        common.legend = TRUE,
                        legend="bottom")
linesBIC_alpha_smax


linesBIC_alpha_smax_sig<-plot_grid(linesBIC_alpha_smax,lineBIC_sig,  rel_widths = c(2, 1) )
linesBIC_alpha_smax_sig
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/confusion_matrices/lineBICsens_alpha_smax_sig.png",
 plot=linesBIC_alpha_smax_sig, width = 18,height = 6)



linesBIC_alpha_smax_sigbar<-plot_grid(linesBIC_alpha_smax,barBIC_sig,  rel_widths = c(2, 1) )
linesBIC_alpha_smax_sigbar
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/confusion_matrices/lineBICsens_alpha_smax_sigbar.png",
 plot=linesBIC_alpha_smax_sigbar, width = 18,height = 6)


#blue and red 

linesBIC_alpha_smax_br <- ggarrange(lineBIC_a_br, lineBICsmax_br, 
                        nrow = 1, ncol = 2,
                        common.legend = TRUE,
                        legend="bottom")
linesBIC_alpha_smax_br
linesBIC_alpha_smax_sigbar_br<-plot_grid(linesBIC_alpha_smax_br,barBIC_sig,  rel_widths = c(2, 1) )
linesBIC_alpha_smax_sigbar_br


ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/confusion_matrices/lineBICsens_alpha_smax_sigbar_br.png",
 plot=linesBIC_alpha_smax_sigbar_br, width = 18,height = 6)






#---------------------------------------------------------------------------------------------
#stan



#===================================================================================
reslfo<-readRDS(file = "outs/simest/sigmamed_sensitivity/resstanloo_sigmed.rds")
head(reslfo)


reslfo$scenario<-case_match(
 reslfo$scenario,
 "sigmamed_stationary"~"stationary",
 "sigmamed_decLinearProd"~"decLinearProd",        
 "sigmamed_regimeProd"~ "regimeProd",            
 "sigmamed_sineProd"~ "sineProd",             
 "sigmamed_regimeCap"~ "regimeCap",             
 "sigmamed_decLinearCap"~ "decLinearCap",        
 "sigmamed_regimeProdCap"~ "regimeProdCap",
 "sigmamed_shiftCap"~"shiftCap",             
 "sigmamed_decLinearProdshiftCap"~"decLinearProdshiftCap")



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
conf_matrixsigmedlfo<-expand.grid(EM=EM,OM=scn)
conf_matrixsigmedlfo$LFOmcmc=NA

cn3<-list()
lfomcmc_set<-list()

o=0
for(a in seq_along(scn)){

  lfoa=subset(reslfo,scenario==scn[a])
  lfomcmc_set[[a]]=tidyr::spread(lfoa[,-9],key=model,value=est)
  nsim<-length(unique( lfomcmc_set[[a]]$iteration))
  lfomcmc_set[[a]]=lfomcmc_set[[a]][c(27,8,18,9,24,15,21,12)] #reorder estimation models
  head(lfomcmc_set[[a]])

  sc3=apply(lfomcmc_set[[a]],1,which.max)
  cn3[[a]]=summary(factor(sc3,levels=seq(1:ncol(lfomcmc_set[[a]]))))/nsim


  myseq<-seq(from=o+1, length.out=length(EM))
  conf_matrixsigmedlfo$LFOmcmc[myseq]<-cn3[[a]]
  o=max(myseq)

}


conf_matrixsigmedlfo$eqem_om <- dplyr::recode(conf_matrixsigmedlfo$OM, 
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


   conf_matrixsigmedlfo$eqem_om<-factor(conf_matrixsigmedlfo$eqem_om,levels=c("stationary",
     "autocorr",
     "dynamic.a","regime.a",
     "dynamic.b","regime.b",
     "dynamic.ab",
     "regime.ab"))  
conf_matrixsigmedlfo$diag<-conf_matrixsigmedlfo$eqem_om==conf_matrixsigmedlfo$EM

conf_matrixsigmedlfo$scentype<-dplyr::case_match(as.character(conf_matrixsigmedlfo$OM), 
      "stationary"~"stationary",
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

conf_matrixsigmedlfo$scendesc<-dplyr::case_match(as.character(conf_matrixsigmedlfo$OM), 
      "stationary"~"",
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
 
conf_matrixsigmedlfo$fullname <- apply( conf_matrixsigmedlfo[ , c("scendesc", "scentype") ] ,
                                   1 , paste , collapse = " " )


conf_matrixsigmedlfo$fullname<-factor(conf_matrixsigmedlfo$fullname, levels=c(
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


pmclfo=ggplot(data =  conf_matrixsigmedlfo, mapping = aes(x = fullname, y = EM)) +
  geom_tile(aes(fill = LFOmcmc), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(LFOmcmc,2)), vjust = 1, size=6) +
  ggtitle(expression(paste("LFO", ~sigma, "=0.3")))+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrixsigmedlfo, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray70", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
pmclfo
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/confusion_matrices/sig_med/LFO_MCMC_sigmed.png",
 plot=pmclfo, width = 9,height = 7)


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







#ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/confusion_matrices/LFOstanbase.png", plot=pmclfo)


#routines to visualize smulation estimation Results
#Catarina Wor
#November 2022
#============================================


library("ggpubr")
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

#========================================================================================================
#sensitivity smax scenario
#read in data

restmb_smax<-readRDS(file = "outs/simest/Smax_sensitivity/res_smax.rds")

aic_smax=subset(restmb_smax, parameter=='AIC'&method=='MLE')
bic_smax=subset(restmb_smax, parameter=='BIC'&method=='MLE')
lfo_smax=subset(restmb_smax, parameter=='LFO'&method=='MLE')


lfo_smax<-lfo_smax[lfo_smax$model %in% c("simple","autocorr","rwa","rwb","rwab","hmma","hmmb","hmmab"),]


aic_smax$mode[aic_smax$convergence>0]<-Inf
bic_smax$mode[aic_smax$convergence>0]<-Inf

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
  "regimeSmax300" 
 ) )


EM=c("stationary",
     "autocorr",
     "dynamic.a","regime.a",
     "dynamic.b","regime.b",
     "dynamic.ab",
     "regime.ab")
##Confusion matrices
conf_matrix_smax<-expand.grid(EM=EM,OM=scn)
conf_matrix_smax$w_AIC=NA
conf_matrix_smax$BIC=NA
conf_matrix_smax$LFO=NA

cn1<-list()
cn2<-list()
cn3<-list()


#summarize model selection
aic_set=list()
bic_set=list()
lfo_set=list()


o=0
for(a in seq_along(scn)){

  #head(aic_set[[a]])

   #AIC
  aica<-subset(aic_smax,scenario==scn[a])
  aic_set[[a]]=tidyr::spread(aica[,-c(10,11,12,13)],key=model,value=mode)
  aic_set[[a]]=aic_set[[a]][c(15,8,12,9,14,11,13,10)] #reorder estimation models


  sc1=apply(aic_set[[a]],1,which.min)
  cn1[[a]]=summary(factor(sc1,levels=seq(1:ncol(aic_set[[a]]))))/1000

  bica=subset(bic_smax,scenario==scn[a])
  bic_set[[a]]=tidyr::spread(bica[,-c(10,11,12,13)],key=model,value=mode)
  bic_set[[a]]=bic_set[[a]][c(15,8,12,9,14,11,13,10)] #reorder estimation models
    head( bic_set[[a]])

  sc2=apply(bic_set[[a]],1,which.min)
  cn2[[a]]=summary(factor(sc2,levels=seq(1:ncol(bic_set[[a]]))))/1000

  lfoa=subset(lfo_smax,scenario==scn[a])
  lfo_set[[a]]=tidyr::spread(lfoa[,-c(10,11,12,13)],key=model,value=mode)
  lfo_set[[a]]=lfo_set[[a]][c(15,8,12,9,14,11,13,10)] #reorder estimation models

  sc3=apply(lfo_set[[a]],1,which.max)
  cn3[[a]]=summary(factor(sc3,levels=seq(1:ncol(lfo_set[[a]]))))/1000


  myseq<-seq(from=o+1, length.out=length(EM))
  conf_matrix_smax$w_AIC[myseq]<-cn1[[a]]
  conf_matrix_smax$BIC[myseq]<-cn2[[a]]
  conf_matrix_smax$LFO[myseq]<-cn3[[a]]
  o=max(myseq)

  
}


conf_matrix_smax$eqem_om <- dplyr::recode(conf_matrix_smax$OM, 
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

conf_matrix_smax$eqem_om<-factor(conf_matrix_smax$eqem_om, 
    levels=c("stationary", 
        "autocorr", 
        "dynamic.a", "regime.a",
        "dynamic.b", "regime.b", 
        "dynamic.ab", "regime.ab"        
       ))
conf_matrix_smax$diag<-conf_matrix_smax$eqem_om==conf_matrix_smax$EM

conf_matrix_smax$EM


conf_matrix_smax$OM2<-dplyr::case_match(conf_matrix_smax$OM,
  "trendLinearSmax025" ~ "trend 25% Smax",
  "trendLinearSmax050" ~ "trend 50% Smax",
  "trendLinearSmax150" ~ "trend 150% Smax",
  "trendLinearSmax200" ~ "trend 200% Smax",
  "trendLinearSmax300" ~ "trend 300% Smax",
  "regimeSmax025" ~ "regime 25% Smax",
  "regimeSmax050" ~ "regime 50% Smax",
  "regimeSmax150" ~ "regime 150% Smax",
  "regimeSmax200" ~ "regime 200% Smax",
  "regimeSmax300"~ "regime 300% Smax")

conf_matrix_a$OM2<-factor(conf_matrix_a$OM2, levels=c( "trend 25% Smax", 
 "trend 50% Smax",
  "trend 150% Smax", 
  "trend 200% Smax",
  "trend 300% Smax",
  "regime 25% Smax",
  "regime 50% Smax", 
  "regime 150% Smax",
  "regime 200% Smax",      
  "regime 300% Smax"))


pa_aic_smax=ggplot(data =  conf_matrix_smax, mapping = aes(x = OM2, y = EM)) +
  geom_tile(aes(fill = w_AIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_AIC,2)), vjust = 1,size=6) +
  ggtitle(expression("AIC Sensitivity"~S[max]))+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix_smax, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray90", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1)) +
  xlab("Simulation Scenario")+ylab("Estimation Model")
pa_aic_smax
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/confusion_matrices/sens_smax/AIC_MLE_senssmax.png",
 plot=pa_aic_smax,  width = 9,height = 7)


pa_bic_smax=ggplot(data =  conf_matrix_smax, mapping = aes(x = OM2, y = EM)) +
  geom_tile(aes(fill = BIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(BIC,2)), vjust = 1,size=6) +
  ggtitle(expression("BIC Sensitivity"~S[max]))+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix_smax, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray90", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
pa_bic_smax
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/confusion_matrices/sens_smax/BIC_MLE_senssmax.png",
 plot=pa_bic_smax,  width = 9,height = 7)



p_lfo_smax=ggplot(data =  conf_matrix_smax, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = LFO), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(LFO,2)), vjust = 1, size=6) +
  ggtitle("LFO")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix_smax, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray90", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
p_lfo_smax
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/confusion_matrices/sens_smax/LFO_MLE_senssmax.png",
 plot=p_lfo)




#calculate these 
conf_matrix_smax$diffsmax_ch<- dplyr::case_match(conf_matrix_smax$OM,
  "trendLinearSmax025"~"decrease x0.25",
   "trendLinearSmax050"~"decrease x0.5",
   "trendLinearSmax150"~"increase x1.5",
 "trendLinearSmax200"~"increase x2.0", 
 "trendLinearSmax300"~"increase x3.0", 
 "regimeSmax025"~"decrease x0.25",     
 "regimeSmax050"~"decrease x0.5",      
 "regimeSmax150"~"increase x1.5",      
 "regimeSmax200"~"increase x2.0",      
 "regimeSmax300" ~"increase x3.0")   




conf_matrix_smax$diffsmax<- dplyr::case_match(conf_matrix_smax$OM,
  "trendLinearSmax025"~4,
   "trendLinearSmax050"~2,
   "trendLinearSmax150"~1.5,
 "trendLinearSmax200"~2.0, 
 "trendLinearSmax300"~3.0, 
 "regimeSmax025"~4,     
 "regimeSmax050"~2,      
 "regimeSmax150"~1.5,      
 "regimeSmax200"~2.0,      
 "regimeSmax300" ~3.0)



conf_matrix_smax$type<- dplyr::case_match(conf_matrix_smax$OM,
 "trendLinearSmax025"~"decrease trend",
   "trendLinearSmax050"~"decrease trend",
   "trendLinearSmax150"~"increase trend",
 "trendLinearSmax200"~"increase trend", 
 "trendLinearSmax300"~"increase trend", 
 "regimeSmax025"~"decrease regime",     
 "regimeSmax050"~"decrease regime",      
 "regimeSmax150"~"increase regime",      
 "regimeSmax200"~"increase regime",      
 "regimeSmax300" ~"increase regime")   

 
conf_matrix_right_smax<-conf_matrix_smax[conf_matrix_smax$EM=="dynamic.b",]

lineAICsmax<-ggplot(conf_matrix_right_smax)+
geom_point(aes(x=diffsmax,y=w_AIC,color=type),size=3)+
geom_line(aes(x=diffsmax,y=w_AIC,color=type,linetype=type),linewidth=1.2)+
scale_colour_manual(values = c("#95D840FF","#482677FF", "#95D840FF", "#482677FF"))+
scale_linetype_manual(values = c(1,1,2,2))+
ylab("     ")+
xlab("multiplication factor for Smax")+
coord_cartesian(ylim=c(0.1,0.88))+
mytheme +
theme(axis.title=element_text(size=14,face="bold"))
lineAICsmax

ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/confusion_matrices/sens_smax/lineAICsenssmax.png",
 plot=lineAICsmax,width = 10,height = 6)



lineBICsmax<-ggplot(conf_matrix_right_smax)+
geom_point(aes(x=diffsmax,y=BIC,color=type),size=3)+
geom_line(aes(x=diffsmax,y=BIC,color=type,linetype=type),linewidth=1.2)+
scale_colour_manual(values = c("#95D840FF","#482677FF", "#95D840FF", "#482677FF"))+
scale_linetype_manual(values = c(1,1,3,3))+
ylab("% of correct model assignment with BIC")+
xlab(expression("multiplication factor for"~ S[max]))+
#scale_color_viridis_d(begin=.1, end=.8) +
mytheme
lineBICsmax

ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/confusion_matrices/sens_a/lineBICsenssmax.png",
 plot=lineBICsmax)


#different colous schemes
lineAICsmax_br<-ggplot(conf_matrix_right_smax)+
geom_point(aes(x=diffsmax,y=w_AIC,color=type),size=3)+
geom_line(aes(x=diffsmax,y=w_AIC,color=type,linetype=type),linewidth=1.2)+
scale_colour_manual(values = c("#B2182B","#B2182B","#2166AC",  "#2166AC"))+
scale_linetype_manual(values = c(2,1,2,1))+
#scale_colour_manual(values = c("#95D840FF","#482677FF", "#95D840FF", "#482677FF"))+
#scale_linetype_manual(values = c(1,1,2,2))+
ylab("     ")+
xlab(expression("multiplication factor for"~ S[max]))+
coord_cartesian(ylim=c(0.1,0.88))+
mytheme +
theme(axis.title=element_text(size=14,face="bold"))
lineAICsmax_br




lineBICsmax_br<-ggplot(conf_matrix_right_smax)+
geom_point(aes(x=diffsmax,y=BIC,color=type),size=4)+
geom_line(aes(x=diffsmax,y=BIC,color=type,linetype=type),linewidth=1.2)+
scale_colour_manual(values = c("#B2182B","#B2182B","#2166AC",  "#2166AC"))+
scale_linetype_manual(values = c(2,1,2,1))+
ylab("     ")+
xlab(expression("multiplication factor for"~ S[max]))+
#scale_color_viridis_d(begin=.1, end=.8) +
mytheme
lineBICsmax_br


#----------------------------------
#need to run plots_cluster_confmat_sensa before this works

linesAIC_alpha_smax <- ggarrange(lineAIC_a, lineAICsmax,
                        nrow = 1, ncol = 2,
                        common.legend = TRUE,
                        legend="bottom")
linesAIC_alpha_smax

ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/confusion_matrices/lineAICsens_alpha_smax.png",
 plot=linesAIC_alpha_smax)


#blue and red colors




linesBIC_alpha_smax_br <- ggarrange(lineBIC_a_br, lineBICsmax_br,
                        nrow = 1, ncol = 2,
                        common.legend = TRUE,
                        legend="bottom")
linesBIC_alpha_smax_br

ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/confusion_matrices/lineBICsens_alpha_smax_br.png",
 plot=linesBIC_alpha_smax_br)















#LFO
reslfo<-readRDS(file = "outs/simest/Smax_sensitivity/resstanloo_smax.rds")
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
        "dynamic.a", "regime.a",
        "dynamic.b", "regime.b", 
        "dynamic.ab", "regime.ab" 
                 
       ))
conf_matrix$diag<-conf_matrix$eqem_om==conf_matrix$EM



conf_matrix$OM2<-dplyr::case_match(conf_matrix$OM,
  "trendLinearSmax025" ~ "trend 25% Smax",
  "trendLinearSmax050" ~ "trend 50% Smax",
  "trendLinearSmax150" ~ "trend 150% Smax",
  "trendLinearSmax200" ~ "trend 200% Smax",
  "trendLinearSmax300" ~ "trend 300% Smax",
  "regimeSmax025" ~ "regime 25% Smax",
  "regimeSmax050" ~ "regime 50% Smax",
  "regimeSmax150" ~ "regime 150% Smax",
  "regimeSmax200" ~ "regime 200% Smax",
  "regimeSmax300"~ "regime 300% Smax")

conf_matrix$OM2<-factor(conf_matrix$OM2, levels=c( "trend 25% Smax", 
 "trend 50% Smax",
  "trend 150% Smax", 
  "trend 200% Smax",
  "trend 300% Smax",
  "regime 25% Smax",
  "regime 50% Smax", 
  "regime 150% Smax",
  "regime 200% Smax",      
  "regime 300% Smax"))


pmclfo_senssmax=ggplot(data =  conf_matrix, mapping = aes(x = OM2, y = EM)) +
  geom_tile(aes(fill = LFOmcmc), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(LFOmcmc,2)), vjust = 1, size=6) +
  ggtitle(expression("LFO Sensitivity"~S[max]))+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray70", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
pmclfo_senssmax
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/confusion_matrices/sens_smax/LFO_MCMC_senssmax.png",
 plot=pmclfo_senssmax,  width = 9,height = 7)


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





pmclfo3_senssmax=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = LFO3mcmc), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(LFO3mcmc,2)), vjust = 1, size=6) +
  ggtitle("sens smax LFO MCMC avg 3 yrs")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray70", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
pmclfo3_senssmax
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/confusion_matrices/sens_smax/LFO_MCMC_3yrs_smax.png",
 plot=pmclfo3_senssmax)




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





pmclfo5_senssmax=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = LFO5mcmc), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(LFO5mcmc,2)), vjust = 1, size=6) +
  ggtitle("sens smax LFO MCMC avg 5 yrs")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray70", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
pmclfo5_senssmax
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/confusion_matrices/sens_a/LFO_MCMC_5yrs_smax.png",
 plot=pmclfo5_senssmax)




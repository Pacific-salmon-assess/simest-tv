#routines to visualize smulation estimation Results
#Catarina Wor
#November 2022
#============================================




library(ggplot2)
library(gridExtra)
library(dplyr)
source("code/utils.R")
source("code/cluster_func_plots.R")

library("ggpubr")

mytheme = list(
    theme_classic(14)+
        theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),
            panel.border = element_rect(fill = NA, color = "black"), legend.title = element_blank(),
            legend.position="bottom", strip.text = element_text( size=11),
            axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),
            plot.title = element_text(face = "bold", hjust = 0.5,size=12))
)


#========================================================================================================
#base case
#read in data
source("code/read_base_data.R")
#exclude simple model
resparam<-resparam[resparam$model!="simple",]
resparam<-resparam[resparam$method!="HMC",]

resparam$model<-factor(resparam$model, levels=c( "autocorr", "rwa", "rwb", "rwab", "hmma","hmmb",  "hmmab"   ))

resparam$modeltype<-dplyr::case_match(resparam$model,
     "autocorr"~"simple",
      "rwa"~"RW",
       "hmma"~"HMM", 
       "rwb"~"RW", 
       "hmmb"~"HMM", 
       "rwab"~"RW",  
       "hmmab"~"HMM")



resparam$scenario<-factor(resparam$scenario, levels=c("stationary",
                                        "autocorr",
                                        "sigmaShift",
                                        "decLinearProd",
                                        "sineProd", 
                                        "regimeProd", 
                                        "shiftProd",                       
                                        "decLinearCap",
                                        "regimeCap",
                                        "shiftCap",                      
                                        "regimeProdCap",         
                                        "decLinearProdshiftCap"  ))





resparam$scentype<-dplyr::case_match(resparam$scenario, 
      "stationary"~"stationary",
      "autocorr"~"stationary",
      "sigmaShift"~"stationary", 
      "decLinearProd"~"trend-a",
      "sineProd"~"trend-a",
      "regimeProd"~"shift-a",
      "shiftProd"~"shift-a",
      "decLinearCap"~"trend-smax",
      "regimeCap"~"shift-smax",
      "shiftCap"~"shift-smax", 
      "regimeProdCap"~"shift-both",
      "decLinearProdshiftCap"~"tv-both"
      )   

resparam$scentype<-factor(resparam$scentype,levels=c( "stationary", "trend-a", "shift-a",
   "trend-smax", "shift-smax", "tv-both",  "shift-both"  ))
unique(resparam$scentype)
#mean percent absolute bias
mapbias<-aggregate(resparam$pbias, list(parameter=resparam$parameter,
                               scenario=resparam$scenario,
                               method=resparam$method,
                               model=resparam$model,
                               modeltype=resparam$modeltype,
                               scentype=resparam$scentype), function(x){median(abs(x),na.rm=T)})


head(mapbias)


#alpha
mapbias_p_alpha<-ggplot(mapbias[mapbias$parameter%in%c("alpha"),])+
geom_bar(aes(x=model,y=x,fill=modeltype),stat="identity",position="dodge")+
facet_wrap(scentype~scenario)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
mytheme+
ylab(expression(paste("median absolute % bias on log(",~alpha,")")))+
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
mapbias_p_alpha 
ggsave("C:/Users/worc/OneDrive - DFO-MPO/timevar_simeval/summary_stats/median_abs_pct_bias_alpha.png",
    plot=mapbias_p_alpha, width = 10,height = 7)
                        



#Smax


head(mapbias[mapbias$parameter%in%c("smax"),])

mapbias_p_smax<-ggplot(mapbias[mapbias$parameter%in%c("smax"),])+
geom_bar(aes(x=model,y=x,fill=modeltype),stat="identity",position="dodge")+
facet_wrap(scentype~scenario)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
mytheme+
ylab(expression(paste("median absolute % bias on ",~S[max])))+
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
mapbias_p_smax
ggsave("C:/Users/worc/OneDrive - DFO-MPO/timevar_simeval/summary_stats/median_abs_pct_bias_smax.png",
    plot=mapbias_p_smax, width = 10,height = 7)


unique(mapbias$parameter)

mapbias_p_smsy<-ggplot(mapbias[mapbias$parameter%in%c("smsy"),])+
geom_bar(aes(x=model,y=x,fill=modeltype),stat="identity",position="dodge")+
facet_wrap(scentype~scenario)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
mytheme+
ylab(expression(paste("median absolute % bias on ",~S[MSY],"")))+
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
mapbias_p_smsy
ggsave("C:/Users/worc/OneDrive - DFO-MPO/timevar_simeval/summary_stats/median_abs_pct_bias_smsy.png",
    plot=mapbias_p_smsy, width = 10,height = 7)




mapbias_p_sgen<-ggplot(mapbias[mapbias$parameter%in%c("sgen"),])+
geom_bar(aes(x=model,y=x,fill=modeltype),stat="identity",position="dodge")+
facet_wrap(scentype~scenario)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
mytheme+
ylab(expression(paste("median absolute % bias on ",~S[gen],"")))+
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
mapbias_p_sgen
ggsave("C:/Users/worc/OneDrive - DFO-MPO/timevar_simeval/summary_stats/median_abs_pct_bias_sgen.png",
    plot=mapbias_p_sgen, width = 10,height = 7)


#there is an error in UMSY calcs, SMY being reported instead for some models
#mapbias_p_umsy<-ggplot(mapbias[mapbias$parameter%in%c("umsy"),])+
#geom_bar(aes(x=model,y=x,fill=method),stat="identity",position="dodge")+
#facet_wrap(~scenario)+
#scale_color_viridis_d(begin=.1, end=.8) +
#scale_fill_viridis_d(begin=.1, end=.8) +
#mytheme+
#coord_cartesian(ylim = c(0,50))+
#ylab(expression(paste("median absolute % bias on ",U[MSY],"")))+
#theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
#mapbias_p_umsy
#ggsave("C:/Users/worc/OneDrive - DFO-MPO/timevar_simeval/summary_stats/mean_abs_pct_bias_smsy.png",
#    plot=mapbias_p_smsy)



mpbias<-aggregate(resparam$pbias, list(parameter=resparam$parameter,
                               scenario=resparam$scenario,
                               method=resparam$method,
                               model=resparam$model), function(x){median((x),na.rm=T)})






#year based cv
cvdf<-aggregate(resparam$mode, list(parameter=resparam$parameter,
                               scenario=resparam$scenario,
                               method=resparam$method,
                               model=resparam$model,
                               by=resparam$by,
                               modeltype=resparam$modeltype,
                               scentype=resparam$scentype ), function(x){sd(x,na.rm=T)/abs(mean(x,na.rm=T))})

sddf<-aggregate(resparam$mode, list(parameter=resparam$parameter,
                               scenario=resparam$scenario,
                               method=resparam$method,
                               model=resparam$model,
                               by=resparam$by,
                               modeltype=resparam$modeltype,
                               scentype=resparam$scentype ), function(x){sd(x,na.rm=T)})


avgdf<-aggregate(resparam$mode, list(parameter=resparam$parameter,
                               scenario=resparam$scenario,
                               method=resparam$method,
                               model=resparam$model,
                               by=resparam$by,
                               modeltype=resparam$modeltype,
                               scentype=resparam$scentype  ), function(x){mean(x,na.rm=T)})


meancvdf<-aggregate(cvdf$x, list(parameter=cvdf$parameter,
                               scenario=cvdf$scenario,
                               method=cvdf$method,
                               model=cvdf$model,
                                modeltype=cvdf$modeltype,
                               scentype=cvdf$scentype 
                              ), function(x){mean(x,na.rm=T)})



meancv_alpha_p<-ggplot(meancvdf[meancvdf$parameter%in%c("alpha"),])+
geom_bar(aes(x=model,y=x,fill=modeltype),stat="identity",position="dodge")+
facet_wrap(scentype~scenario)+
scale_color_viridis_d(begin=.1, end=.8, option = "A") +
scale_fill_viridis_d(begin=.1, end=.8, option = "A") +
mytheme+
ylab(expression(paste("mean cv for by estimates of log(",~alpha,")")))+
#coord_cartesian(ylim = c(0,.5))+ 
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
meancv_alpha_p
ggsave("C:/Users/worc/OneDrive - DFO-MPO/timevar_simeval/summary_stats/mean_by_cv_alpha.png",
    plot=meancv_alpha_p, width = 10,height = 7)
  


meancv_smax_p<-ggplot(meancvdf[meancvdf$parameter%in%c("smax"),])+
geom_bar(aes(x=model,y=x,fill=modeltype),stat="identity",position="dodge")+
facet_wrap(scentype~scenario)+
scale_color_viridis_d(begin=.1, end=.8, option = "A") +
scale_fill_viridis_d(begin=.1, end=.8, option = "A") +
mytheme+
ylab(expression(paste("mean cv for by estimates of ",~S[max])))+
coord_cartesian(ylim = c(0,.5))+ 
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
meancv_smax_p
ggsave("C:/Users/worc/OneDrive - DFO-MPO/timevar_simeval/summary_stats/mean_by_cv_smax.png",
    plot=meancv_smax_p, width = 10,height = 7)




meancv_smsy_p<-ggplot(meancvdf[meancvdf$parameter%in%c("smsy"),])+
geom_bar(aes(x=model,y=x,fill=modeltype),stat="identity",position="dodge")+
facet_wrap(scentype~scenario)+
scale_color_viridis_d(begin=.1, end=.8, option = "A") +
scale_fill_viridis_d(begin=.1, end=.8, option = "A") +
mytheme+
ylab(expression(paste("mean cv for by estimates of ",S[MSY])))+
coord_cartesian(ylim = c(0,1.5))+ 
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
meancv_smsy_p
ggsave("C:/Users/worc/OneDrive - DFO-MPO/timevar_simeval/summary_stats/mean_by_cv_smsy.png",
    plot=meancv_smsy_p, width = 10,height = 7)


meancv_sgen_p<-ggplot(meancvdf[meancvdf$parameter%in%c("sgen"),])+
geom_bar(aes(x=model,y=x,fill=modeltype),stat="identity",position="dodge")+
facet_wrap(scentype~scenario)+
scale_color_viridis_d(begin=.1, end=.8, option = "A") +
scale_fill_viridis_d(begin=.1, end=.8, option = "A") +
mytheme+
ylab(expression(paste("mean cv for by estimates of ",S[gen])))+
coord_cartesian(ylim = c(0,.7))+ 
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
meancv_sgen_p
ggsave("C:/Users/worc/OneDrive - DFO-MPO/timevar_simeval/summary_stats/mean_by_cv_sgen.png",
    plot=meancv_sgen_p, width = 10,height = 7)

#write out the data
mapbias$mean.abs.pct.bias<-mapbias$x
mapbias$mean.pct.bias<-mpbias$x
mapbias$mean.cv<-meancvdf$x


mapbias<-select(mapbias,!c("x"))

mapbias$scentype<-case_match(mapbias$scenario,
      "stationary"~"stationary",
      "autocorr"~"stationary",
      "sigmaShift"~"stationary", 
      "decLinearProd"~"tv-a",
      "sineProd"~"tv-a",
      "regimeProd"~"tv-a",
      "shiftProd"~"tv-a",
      "decLinearCap"~"tv-smax",
      "regimeCap"~"tv-smax",
      "shiftCap"~"tv-smax", 
      "regimeProdCap"~"tv-both",
      "decLinearProdshiftCap"~"tv-both"
      )   


write.csv(mapbias,file = "C:/Users/worc/OneDrive - DFO-MPO/timevar_simeval/summary_stats/summary_stats_base.csv", row.names = FALSE)





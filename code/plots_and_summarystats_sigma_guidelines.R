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
    theme_classic(11)+
        theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),
            panel.border = element_rect(fill = NA, color = "black"), legend.title = element_blank(),
            legend.position="bottom", strip.text = element_text( size=12),
            axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),
            plot.title = element_text(face = "bold", hjust = 0.5,size=12))
)



#===================================================================================================
#--------------------------------------------------------------------------------------------------
# sensitivity sig low scenarios. 

#sensitivity a
simPar <- read.csv("data/sensitivity/SimPars.csv")

## Store relevant object names to help run simulation 
scenNames <- unique(simPar$scenario)

ressiglow1<-readRDS(file = "outs/simest/sigmalow_sensitivity2/res_siglow1.rds")
ressiglow2<-readRDS(file = "outs/simest/sigmalow_sensitivity2/res_siglow2.rds")
ressiglow15<-readRDS(file = "outs/simest/sigmalow_sensitivity2/res_siglow_15.rds")


resertmb<-rbind(ressiglow1,ressiglow2,ressiglow15)

#resstaner1<-readRDS(file = "outs/simest/genericER/resstan_erq1.rds")
#resstaner2<-readRDS(file = "outs/simest/genericER/resstan_erq2.rds")
#resstaner3<-readRDS(file = "outs/simest/genericER/resstan_erq3.rds")
#resstaner4<-readRDS(file = "outs/simest/genericER/resstan_erq4.rds")

#reserstan<-rbind(resstaner1,resstaner2,resstaner3,resstaner4)


ressiglow<-rbind(resertmb)#,reserstan)

#res<-resstan
ressiglow$parameter[ressiglow$parameter=="Smax"]<-"smax"
#reser$method[reser$method=="MCMC"]<-"HMC"
resaparam<-ressiglow[ressiglow$parameter%in%c("alpha","smax","smsy","sgen","umsy"),]

#exclude outliers

resaparam$convergence[resaparam$parameter=="alpha"&resaparam$mode>40]<-1
resaparam$convergence[resaparam$parameter=="smax"&resaparam$mode>1e8]<-1

#unique(resaparam$iteration[resaparam$parameter=="smax"&resaparam$mode>1e8])


#convstat<-aggregate(resaparam$convergence,
#    list(scenario=resaparam$scenario,
#        model=resaparam$model,
#        method=resaparam$method,
#        iteration=resaparam$iteration),
#    function(x){sum(x)})
#convstatMLE<-convstat[convstat$x==0&convstat$method=="MLE",]
#convstatMCMC<-convstat[convstat$x==0&convstat$method=="HMC",]


#allconv<-inner_join(convstatMLE[,-3], convstatMCMC[,-3])

#convsum<-aggregate(allconv$iteration,
#    list(model=allconv$model,scenario=allconv$scenario),
#    function(x){length(unique(x))})



#conv_iter<-aggregate(allconv$iteration,
#    list(model=allconv$model,scenario=allconv$scenario),
#    function(x){(unique(x))})

#convsnc<-as.numeric(rownames(convsum))

#resl<-list()
#for(i in seq_along(convsnc)){
#    sel<-conv_iter[convsnc[i],]
#    resl[[i]]<-resaparam %>% filter(model==sel$model&
#                            scenario==sel$scenario&
#                            iteration%in%sel$x[[1]])   
#}

#resaparam<-as.data.frame(data.table::rbindlist(resl))

resaparam$model<-factor(resaparam$model, levels=c("simple", "autocorr", "rwa", "hmma", "rwb", "hmmb", "rwab",  "hmmab"   ))

unique(resaparam$scenario)


resaparam$scens<-case_match(resaparam$scenario,"sigmalow_stationary"~"stationary",
            "sigmalow_autocorr"~"autocorr", 
            "sigmalow_decLinearProd"~"decLinearProd",        
            "sigmalow_regimeProd"~"regimeProd", 
            "sigmalow_sineProd"~"sineProd" , 
            "sigmalow_regimeCap"~"regimeCap" ,           
            "sigmalow_decLinearCap"~"decLinearCap",        
            "sigmalow_regimeProdCap"~"regimeProdCap",
            "sigmalow_shiftCap"~"shiftCap" ,            
            "sigmalow_decLinearProdshiftCap"~"decLinearProdshiftCap")

unique(resaparam$scens)

resaparam$scens<-factor(resaparam$scens, levels=c("stationary",
            "autocorr",
            "decLinearProd",         
            "regimeProd",           
            "sineProd" ,             
            "regimeCap",             
            "decLinearCap",          
            "regimeProdCap",        
            "shiftCap",
            "decLinearProdshiftCap"))

head(resaparam)



 ggplot(resparam[resparam$method=="MLE"&
resparam$model=="hmmb"&
resparam$parameter %in% c("smax")&
resparam$scenario== "stationary" ,]) + 
  geom_histogram(aes(x=mode))


  mapbias<-aggregate(resaparam$pbias, list(parameter=resaparam$parameter,
                               scenario=resaparam$scenario,
                               scens=resaparam$scens,
                               method=resaparam$method,
                               model=resaparam$model), function(x){median(abs(x),na.rm=T)})

head( mapbias)

#alpha
mapbias_p_alpha<-ggplot(mapbias[mapbias$parameter%in%c("alpha"),])+
geom_bar(aes(x=model,y=x,fill=method),stat="identity",position="dodge")+
facet_wrap(~scens, ncol=3)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
mytheme+
ylab(expression(paste("median absolute % bias on log(",alpha,")")))+
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
mapbias_p_alpha 
ggsave("C:/Users/worc/OneDrive - DFO-MPO/timevar_simeval/summary_stats/sens_siglow_median_abs_pct_bias_alpha.png",
    plot=mapbias_p_alpha)
                        

#Smax


head(mapbias[mapbias$parameter%in%c("smax"),])

mapbias_p_smax<-ggplot(mapbias[mapbias$parameter%in%c("smax"),])+
geom_bar(aes(x=model,y=x,fill=method),stat="identity",position="dodge")+
facet_wrap(~scens, ncol=3)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
mytheme+
ylab(expression(paste("median absolute % bias on ",S[max],"")))+
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
mapbias_p_smax
ggsave("C:/Users/worc/OneDrive - DFO-MPO/timevar_simeval/summary_stats/sens_siglow_median_abs_pct_bias_smax.png",
    plot=mapbias_p_smax)


unique(mapbias$parameter)

mapbias_p_smsy<-ggplot(mapbias[mapbias$parameter%in%c("smsy"),])+
geom_bar(aes(x=model,y=x,fill=method),stat="identity",position="dodge")+
facet_wrap(~scens, ncol=3)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
mytheme+
ylab(expression(paste("median absolute % bias on ",S[MSY],"")))+
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
mapbias_p_smsy
ggsave("C:/Users/worc/OneDrive - DFO-MPO/timevar_simeval/summary_stats/sens_siglow_median_abs_pct_bias_smsy.png",
    plot=mapbias_p_smsy)




mapbias_p_sgen<-ggplot(mapbias[mapbias$parameter%in%c("sgen"),])+
geom_bar(aes(x=model,y=x,fill=method),stat="identity",position="dodge")+
facet_wrap(~scens, ncol=3)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
mytheme+
ylab(expression(paste("median absolute % bias on ",S[gen],"")))+
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
mapbias_p_sgen
ggsave("C:/Users/worc/OneDrive - DFO-MPO/timevar_simeval/summary_stats/sens_siglow_median_abs_pct_bias_sgen.png",
    plot=mapbias_p_sgen)


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



mpbias<-aggregate(resaparam$pbias, list(parameter=resaparam$parameter,
                               scenario=resaparam$scenario,
                               method=resaparam$method,
                               model=resaparam$model), function(x){median((x),na.rm=T)})




mpbias_p<-ggplot(mpbias[mpbias$parameter%in%c("alpha"),])+
geom_bar(aes(x=model,y=x,fill=method),stat="identity",position="dodge")+
facet_wrap(~scenario)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
mytheme+
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
mpbias_p


head(resaparam)

#year based cv
cvdf<-aggregate(resaparam$est, list(parameter=resaparam$parameter,
                               scenario=resaparam$scenario,
                               scens=resaparam$scens,
                               method=resaparam$method,
                               model=resaparam$model,
                               by=resaparam$by ), function(x){sd(x,na.rm=T)/abs(mean(x,na.rm=T))})



sddf<-aggregate(resaparam$est, list(parameter=resaparam$parameter,
                               scenario=resaparam$scenario,
                               scens=resaparam$scens,
                               method=resaparam$method,
                               model=resaparam$model,
                               by=resaparam$by ), function(x){sd(x,na.rm=T)})


avgdf<-aggregate(resaparam$est, list(parameter=resaparam$parameter,
                               scenario=resaparam$scenario,
                               scens=resaparam$scens,
                               method=resaparam$method,
                               model=resaparam$model,
                               by=resaparam$by ), function(x){mean(x,na.rm=T)})



meancvdf<-aggregate(cvdf$x, list(parameter=cvdf$parameter,
                               scenario=cvdf$scenario,
                               scens=cvdf$scens,
                               method=cvdf$method,
                               model=cvdf$model
                              ), function(x){mean(x,na.rm=T)})



meancv_alpha_p<-ggplot(meancvdf[meancvdf$parameter%in%c("alpha"),])+
geom_bar(aes(x=model,y=x,fill=method),stat="identity",position="dodge")+
facet_wrap(~scens,ncol=3)+
scale_color_viridis_d(begin=.1, end=.8, option = "A") +
scale_fill_viridis_d(begin=.1, end=.8, option = "A") +
mytheme+
ylab(expression(paste("mean cv for by estimates of log(",alpha,")")))+
#coord_cartesian(ylim = c(0,2.5))+ 
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
meancv_alpha_p
ggsave("C:/Users/worc/OneDrive - DFO-MPO/timevar_simeval/summary_stats/sens_siglow_mean_by_cv_alpha.png",
    plot=meancv_alpha_p)
  


meancv_smax_p<-ggplot(meancvdf[meancvdf$parameter%in%c("smax"),])+
geom_bar(aes(x=model,y=x,fill=method),stat="identity",position="dodge")+
facet_wrap(~scens,ncol=3)+
scale_color_viridis_d(begin=.1, end=.8, option = "A") +
scale_fill_viridis_d(begin=.1, end=.8, option = "A") +
mytheme+
ylab(expression(paste("mean cv for by estimates of ",S[max])))+
coord_cartesian(ylim = c(0,.8))+ 
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
meancv_smax_p
ggsave("C:/Users/worc/OneDrive - DFO-MPO/timevar_simeval/summary_stats/sens_siglow_mean_by_cv_smax.png",
    plot=meancv_smax_p)




meancv_smsy_p<-ggplot(meancvdf[meancvdf$parameter%in%c("smsy"),])+
geom_bar(aes(x=model,y=x,fill=method),stat="identity",position="dodge")+
facet_wrap(~scens,ncol=3)+
scale_color_viridis_d(begin=.1, end=.8, option = "A") +
scale_fill_viridis_d(begin=.1, end=.8, option = "A") +
mytheme+
ylab(expression(paste("mean cv for by estimates of ",S[MSY])))+
coord_cartesian(ylim = c(0,5))+ 
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
meancv_smsy_p
ggsave("C:/Users/worc/OneDrive - DFO-MPO/timevar_simeval/summary_stats/sens_siglow_mean_by_cv_smsy.png",
    plot=meancv_smsy_p)


meancv_sgen_p<-ggplot(meancvdf[meancvdf$parameter%in%c("sgen"),])+
geom_bar(aes(x=model,y=x,fill=method),stat="identity",position="dodge")+
facet_wrap(~scens,ncol=3)+
scale_color_viridis_d(begin=.1, end=.8, option = "A") +
scale_fill_viridis_d(begin=.1, end=.8, option = "A") +
mytheme+
ylab(expression(paste("mean cv for by estimates of ",S[gen])))+
coord_cartesian(ylim = c(0,1))+ 
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
meancv_sgen_p
ggsave("C:/Users/worc/OneDrive - DFO-MPO/timevar_simeval/summary_stats/sens_siglow_mean_by_cv_sgen.png",
    plot=meancv_sgen_p)

#write out the data
mapbias$mean.abs.pct.bias<-mapbias$x
mapbias$mean.pct.bias<-mpbias$x
mapbias$mean.cv<-meancvdf$x


mapbias<-select(mapbias,!c("x"))




write.csv(mapbias,file = "C:/Users/worc/OneDrive - DFO-MPO/timevar_simeval/summary_stats/summary_stats_sens_siglow.csv", row.names = FALSE)



#========================================================================================================================


#===================================================================================================
#--------------------------------------------------------------------------------------------------
# sensitivity sig medium scenarios. 



## Store relevant object names to help run simulation 

ressigmed1<-readRDS(file = "outs/simest/sigmamed_sensitivity2/res_sigmed1.rds")
ressigmed2<-readRDS(file = "outs/simest/sigmamed_sensitivity2/res_sigmed2.rds")
ressigmed15<-readRDS(file = "outs/simest/sigmamed_sensitivity2/res_sigmed_12.rds")


ressigmedtmb<-rbind(ressigmed1,ressigmed2,ressigmed15)

#resstaner1<-readRDS(file = "outs/simest/genericER/resstan_erq1.rds")
#resstaner2<-readRDS(file = "outs/simest/genericER/resstan_erq2.rds")
#resstaner3<-readRDS(file = "outs/simest/genericER/resstan_erq3.rds")
#resstaner4<-readRDS(file = "outs/simest/genericER/resstan_erq4.rds")

#reserstan<-rbind(resstaner1,resstaner2,resstaner3,resstaner4)


ressigmed<-rbind(ressigmedtmb)#,reserstan)

#res<-resstan
ressigmed$parameter[ressigmed$parameter=="Smax"]<-"smax"
#reser$method[reser$method=="MCMC"]<-"HMC"
resaparam<-ressigmed[ressigmed$parameter%in%c("alpha","smax","smsy","sgen","umsy"),]

#exclude outliers

resaparam$convergence[resaparam$parameter=="alpha"&resaparam$mode>40]<-1
resaparam$convergence[resaparam$parameter=="smax"&resaparam$mode>1e8]<-1

#unique(resaparam$iteration[resaparam$parameter=="smax"&resaparam$mode>1e8])


#convstat<-aggregate(resaparam$convergence,
#    list(scenario=resaparam$scenario,
#        model=resaparam$model,
#        method=resaparam$method,
#        iteration=resaparam$iteration),
#    function(x){sum(x)})
#convstatMLE<-convstat[convstat$x==0&convstat$method=="MLE",]
#convstatMCMC<-convstat[convstat$x==0&convstat$method=="HMC",]


#allconv<-inner_join(convstatMLE[,-3], convstatMCMC[,-3])

#convsum<-aggregate(allconv$iteration,
#    list(model=allconv$model,scenario=allconv$scenario),
#    function(x){length(unique(x))})



#conv_iter<-aggregate(allconv$iteration,
#    list(model=allconv$model,scenario=allconv$scenario),
#    function(x){(unique(x))})

#convsnc<-as.numeric(rownames(convsum))

#resl<-list()
#for(i in seq_along(convsnc)){
#    sel<-conv_iter[convsnc[i],]
#    resl[[i]]<-resaparam %>% filter(model==sel$model&
#                            scenario==sel$scenario&
#                            iteration%in%sel$x[[1]])   
#}

#resaparam<-as.data.frame(data.table::rbindlist(resl))

resaparam$model<-factor(resaparam$model, levels=c("simple", "autocorr", "rwa", "hmma", "rwb", "hmmb", "rwab",  "hmmab"   ))

unique(resaparam$scenario)


resaparam$scens<-case_match(resaparam$scenario,"sigmamed_stationary"~"stationary",
            "sigmamed_autocorr"~"autocorr", 
            "sigmamed_decLinearProd"~"decLinearProd",        
            "sigmamed_regimeProd"~"regimeProd", 
            "sigmamed_sineProd"~"sineProd" , 
            "sigmamed_regimeCap"~"regimeCap" ,           
            "sigmamed_decLinearCap"~"decLinearCap",        
            "sigmamed_regimeProdCap"~"regimeProdCap",
            "sigmamed_shiftCap"~"shiftCap" ,            
            "sigmamed_decLinearProdshiftCap"~"decLinearProdshiftCap")

unique(resaparam$scens)

resaparam$scens<-factor(resaparam$scens, levels=c("stationary",
            "autocorr",
            "decLinearProd",         
            "regimeProd",           
            "sineProd" ,             
            "regimeCap",             
            "decLinearCap",          
            "regimeProdCap",        
            "shiftCap",
            "decLinearProdshiftCap"))

head(resaparam)



 ggplot(resparam[resparam$method=="MLE"&
resparam$model=="hmmb"&
resparam$parameter %in% c("smax")&
resparam$scenario== "stationary" ,]) + 
  geom_histogram(aes(x=mode))


  mapbias<-aggregate(resaparam$pbias, list(parameter=resaparam$parameter,
                               scenario=resaparam$scenario,
                               scens=resaparam$scens,
                               method=resaparam$method,
                               model=resaparam$model), function(x){median(abs(x),na.rm=T)})

head( mapbias)

#alpha
mapbias_p_alpha<-ggplot(mapbias[mapbias$parameter%in%c("alpha"),])+
geom_bar(aes(x=model,y=x,fill=method),stat="identity",position="dodge")+
facet_wrap(~scens, ncol=3)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
mytheme+
ylab(expression(paste("median absolute % bias on log(",alpha,")")))+
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
mapbias_p_alpha 
ggsave("C:/Users/worc/OneDrive - DFO-MPO/timevar_simeval/summary_stats/sens_sigmed_median_abs_pct_bias_alpha.png",
    plot=mapbias_p_alpha)
                        

#Smax


head(mapbias[mapbias$parameter%in%c("smax"),])

mapbias_p_smax<-ggplot(mapbias[mapbias$parameter%in%c("smax"),])+
geom_bar(aes(x=model,y=x,fill=method),stat="identity",position="dodge")+
facet_wrap(~scens, ncol=3)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
mytheme+
ylab(expression(paste("median absolute % bias on ",S[max],"")))+
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
mapbias_p_smax
ggsave("C:/Users/worc/OneDrive - DFO-MPO/timevar_simeval/summary_stats/sens_sigmed_median_abs_pct_bias_smax.png",
    plot=mapbias_p_smax)


unique(mapbias$parameter)

mapbias_p_smsy<-ggplot(mapbias[mapbias$parameter%in%c("smsy"),])+
geom_bar(aes(x=model,y=x,fill=method),stat="identity",position="dodge")+
facet_wrap(~scens, ncol=3)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
mytheme+
ylab(expression(paste("median absolute % bias on ",S[MSY],"")))+
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
mapbias_p_smsy
ggsave("C:/Users/worc/OneDrive - DFO-MPO/timevar_simeval/summary_stats/sens_sigsigmed_median_abs_pct_bias_smsy.png",
    plot=mapbias_p_smsy)




mapbias_p_sgen<-ggplot(mapbias[mapbias$parameter%in%c("sgen"),])+
geom_bar(aes(x=model,y=x,fill=method),stat="identity",position="dodge")+
facet_wrap(~scens, ncol=3)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
mytheme+
ylab(expression(paste("median absolute % bias on ",S[gen],"")))+
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
mapbias_p_sgen
ggsave("C:/Users/worc/OneDrive - DFO-MPO/timevar_simeval/summary_stats/sens_sigmed_median_abs_pct_bias_sgen.png",
    plot=mapbias_p_sgen)


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



mpbias<-aggregate(resaparam$pbias, list(parameter=resaparam$parameter,
                               scenario=resaparam$scenario,
                               method=resaparam$method,
                               model=resaparam$model), function(x){median((x),na.rm=T)})




mpbias_p<-ggplot(mpbias[mpbias$parameter%in%c("alpha"),])+
geom_bar(aes(x=model,y=x,fill=method),stat="identity",position="dodge")+
facet_wrap(~scenario)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
mytheme+
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
mpbias_p


head(resaparam)

#year based cv
cvdf<-aggregate(resaparam$est, list(parameter=resaparam$parameter,
                               scenario=resaparam$scenario,
                               scens=resaparam$scens,
                               method=resaparam$method,
                               model=resaparam$model,
                               by=resaparam$by ), function(x){sd(x,na.rm=T)/abs(mean(x,na.rm=T))})



sddf<-aggregate(resaparam$est, list(parameter=resaparam$parameter,
                               scenario=resaparam$scenario,
                               scens=resaparam$scens,
                               method=resaparam$method,
                               model=resaparam$model,
                               by=resaparam$by ), function(x){sd(x,na.rm=T)})


avgdf<-aggregate(resaparam$est, list(parameter=resaparam$parameter,
                               scenario=resaparam$scenario,
                               scens=resaparam$scens,
                               method=resaparam$method,
                               model=resaparam$model,
                               by=resaparam$by ), function(x){mean(x,na.rm=T)})



meancvdf<-aggregate(cvdf$x, list(parameter=cvdf$parameter,
                               scenario=cvdf$scenario,
                               scens=cvdf$scens,
                               method=cvdf$method,
                               model=cvdf$model
                              ), function(x){mean(x,na.rm=T)})



meancv_alpha_p<-ggplot(meancvdf[meancvdf$parameter%in%c("alpha"),])+
geom_bar(aes(x=model,y=x,fill=method),stat="identity",position="dodge")+
facet_wrap(~scens,ncol=3)+
scale_color_viridis_d(begin=.1, end=.8, option = "A") +
scale_fill_viridis_d(begin=.1, end=.8, option = "A") +
mytheme+
ylab(expression(paste("mean cv for by estimates of log(",alpha,")")))+
#coord_cartesian(ylim = c(0,2.5))+ 
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
meancv_alpha_p
ggsave("C:/Users/worc/OneDrive - DFO-MPO/timevar_simeval/summary_stats/sens_sigmed_mean_by_cv_alpha.png",
    plot=meancv_alpha_p)
  


meancv_smax_p<-ggplot(meancvdf[meancvdf$parameter%in%c("smax"),])+
geom_bar(aes(x=model,y=x,fill=method),stat="identity",position="dodge")+
facet_wrap(~scens,ncol=3)+
scale_color_viridis_d(begin=.1, end=.8, option = "A") +
scale_fill_viridis_d(begin=.1, end=.8, option = "A") +
mytheme+
ylab(expression(paste("mean cv for by estimates of ",S[max])))+
coord_cartesian(ylim = c(0,3))+ 
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
meancv_smax_p
ggsave("C:/Users/worc/OneDrive - DFO-MPO/timevar_simeval/summary_stats/sens_sigmed_mean_by_cv_smax.png",
    plot=meancv_smax_p)




meancv_smsy_p<-ggplot(meancvdf[meancvdf$parameter%in%c("smsy"),])+
geom_bar(aes(x=model,y=x,fill=method),stat="identity",position="dodge")+
facet_wrap(~scens,ncol=3)+
scale_color_viridis_d(begin=.1, end=.8, option = "A") +
scale_fill_viridis_d(begin=.1, end=.8, option = "A") +
mytheme+
ylab(expression(paste("mean cv for by estimates of ",S[MSY])))+
coord_cartesian(ylim = c(0,5))+ 
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
meancv_smsy_p
ggsave("C:/Users/worc/OneDrive - DFO-MPO/timevar_simeval/summary_stats/sens_sigmed_mean_by_cv_smsy.png",
    plot=meancv_smsy_p)


meancv_sgen_p<-ggplot(meancvdf[meancvdf$parameter%in%c("sgen"),])+
geom_bar(aes(x=model,y=x,fill=method),stat="identity",position="dodge")+
facet_wrap(~scens,ncol=3)+
scale_color_viridis_d(begin=.1, end=.8, option = "A") +
scale_fill_viridis_d(begin=.1, end=.8, option = "A") +
mytheme+
ylab(expression(paste("mean cv for by estimates of ",S[gen])))+
coord_cartesian(ylim = c(0,1))+ 
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
meancv_sgen_p
ggsave("C:/Users/worc/OneDrive - DFO-MPO/timevar_simeval/summary_stats/sens_sigmed_mean_by_cv_sgen.png",
    plot=meancv_sgen_p)

#write out the data
mapbias$mean.abs.pct.bias<-mapbias$x
mapbias$mean.pct.bias<-mpbias$x
mapbias$mean.cv<-meancvdf$x


mapbias<-select(mapbias,!c("x"))




write.csv(mapbias,file = "C:/Users/worc/OneDrive - DFO-MPO/timevar_simeval/summary_stats/summary_stats_sens_sigmed.csv", row.names = FALSE)


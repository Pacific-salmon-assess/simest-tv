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




#========================================================================================================
#base case
#read in data
simPar <- read.csv("data/generic/SimPars.csv")

## Store relevant object names to help run simulation 
scenNames <- unique(simPar$scenario)

res1<-readRDS(file = "outs/simest/generic/resbase1.rds")
res2<-readRDS(file = "outs/simest/generic/resbase2.rds")


restmb<-rbind(res1,res2)


resstan1<-readRDS(file = "outs/simest/generic/resstan1.rds")
resstan2<-readRDS(file = "outs/simest/generic/resstan2.rds")
resstan<-rbind(resstan1,resstan2)
#resstan<-readRDS(file = "outs/simest/generic/resstan.rds")

res<-rbind(restmb,resstan)
head(res)
#res<-resstan
res$parameter[res$parameter=="Smax"]<-"smax"
res$method[res$method=="MCMC"]<-"HMC"
resparam<-res[res$parameter%in%c("alpha","smax","smsy","sgen","umsy"),]

#exclude outliers

resparam$convergence[resparam$parameter=="alpha"&resparam$mode>40]<-1
resparam$convergence[resparam$parameter=="smax"&resparam$mode>1e8]<-1



convstat<-aggregate(resparam$convergence,
    list(scenario=resparam$scenario,
        model=resparam$model,
        method=resparam$method,
        iteration=resparam$iteration),
    function(x){sum(x)})
convstatMLE<-convstat[convstat$x==0&convstat$method=="MLE",]
convstatMCMC<-convstat[convstat$x==0&convstat$method=="HMC",]


allconv<-inner_join(convstatMLE[,-3], convstatMCMC[,-3])

convsum<-aggregate(allconv$iteration,
    list(model=allconv$model,scenario=allconv$scenario),
    function(x){length(unique(x))})



conv_iter<-aggregate(allconv$iteration,
    list(model=allconv$model,scenario=allconv$scenario),
    function(x){(unique(x))})

convsnc<-as.numeric(rownames(convsum))

resl<-list()
for(i in seq_along(convsnc)){

    sel<-conv_iter[convsnc[i],]
    resl[[i]]<-resparam %>% filter(model==sel$model&
                            scenario==sel$scenario&
                            iteration%in%sel$x[[1]])
    
}

resparam<-as.data.frame(data.table::rbindlist(resl))

resparam$model<-factor(resparam$model, levels=c("simple", "autocorr", "rwa", "hmma", "rwb", "hmmb", "rwab",  "hmmab"   ))

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







range(resparam[resparam$method=="MLE"&
resparam$parameter %in% c("smax")&
resparam$scenario== "stationary" ,"mode"])


range(resparam[resparam$method=="MLE"&
resparam$parameter %in% c("smsy")&
resparam$scenario== "stationary" ,"mode"])

head(resparam[resparam$method=="HMC"&
resparam$model=="hmma"&
resparam$parameter %in% c("smsy"),])

 ggplot(resparam[resparam$method=="MLE"&
resparam$model=="hmmb"&
resparam$parameter %in% c("smax")&
resparam$scenario== "stationary" ,]) + 
  geom_histogram(aes(x=mode))


#mean percent absolute bias
mapbias<-aggregate(resparam$pbias, list(parameter=resparam$parameter,
                               scenario=resparam$scenario,
                               method=resparam$method,
                               model=resparam$model), function(x){median(abs(x),na.rm=T)})

#alpha
mapbias_p_alpha<-ggplot(mapbias[mapbias$parameter%in%c("alpha"),])+
geom_bar(aes(x=model,y=x,fill=method),stat="identity",position="dodge")+
facet_wrap(~scenario)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
mytheme+
ylab(expression(paste("median absolute % bias on log(",alpha,")")))+
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
mapbias_p_alpha 
ggsave("C:/Users/worc/OneDrive - DFO-MPO/timevar_simeval/summary_stats/median_abs_pct_bias_alpha.png",
    plot=mapbias_p_alpha)
                        

#Smax


head(mapbias[mapbias$parameter%in%c("smax"),])

mapbias_p_smax<-ggplot(mapbias[mapbias$parameter%in%c("smax"),])+
geom_bar(aes(x=model,y=x,fill=method),stat="identity",position="dodge")+
facet_wrap(~scenario)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
mytheme+
ylab(expression(paste("median absolute % bias on ",S[max],"")))+
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
mapbias_p_smax
ggsave("C:/Users/worc/OneDrive - DFO-MPO/timevar_simeval/summary_stats/median_abs_pct_bias_smax.png",
    plot=mapbias_p_smax)


unique(mapbias$parameter)

mapbias_p_smsy<-ggplot(mapbias[mapbias$parameter%in%c("smsy"),])+
geom_bar(aes(x=model,y=x,fill=method),stat="identity",position="dodge")+
facet_wrap(~scenario)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
mytheme+
ylab(expression(paste("median absolute % bias on ",S[MSY],"")))+
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
mapbias_p_smsy
ggsave("C:/Users/worc/OneDrive - DFO-MPO/timevar_simeval/summary_stats/median_abs_pct_bias_smsy.png",
    plot=mapbias_p_smsy)




mapbias_p_sgen<-ggplot(mapbias[mapbias$parameter%in%c("sgen"),])+
geom_bar(aes(x=model,y=x,fill=method),stat="identity",position="dodge")+
facet_wrap(~scenario)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
mytheme+
ylab(expression(paste("median absolute % bias on ",S[gen],"")))+
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
mapbias_p_sgen
ggsave("C:/Users/worc/OneDrive - DFO-MPO/timevar_simeval/summary_stats/median_abs_pct_bias_sgen.png",
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



mpbias<-aggregate(resparam$pbias, list(parameter=resparam$parameter,
                               scenario=resparam$scenario,
                               method=resparam$method,
                               model=resparam$model), function(x){median((x),na.rm=T)})




mpbias_p<-ggplot(mpbias[mpbias$parameter%in%c("alpha"),])+
geom_bar(aes(x=model,y=x,fill=method),stat="identity",position="dodge")+
facet_wrap(~scenario)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
mytheme+
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
mpbias_p




#year based cv
cvdf<-aggregate(resparam$mode, list(parameter=resparam$parameter,
                               scenario=resparam$scenario,
                               method=resparam$method,
                               model=resparam$model,
                               by=resparam$by ), function(x){sd(x,na.rm=T)/abs(mean(x,na.rm=T))})

unique(cvdf$scenario)
cvdf[cvdf$scenario=="decLinearProdshiftCap"&
cvdf$model=="simple"&
cvdf$parameter=="alpha",]

sddf<-aggregate(resparam$mode, list(parameter=resparam$parameter,
                               scenario=resparam$scenario,
                               method=resparam$method,
                               model=resparam$model,
                               by=resparam$by ), function(x){sd(x,na.rm=T)})

sddf[sddf$scenario=="decLinearProdshiftCap"&
sddf$model=="simple"&
sddf$parameter=="alpha",]


avgdf<-aggregate(resparam$mode, list(parameter=resparam$parameter,
                               scenario=resparam$scenario,
                               method=resparam$method,
                               model=resparam$model,
                               by=resparam$by ), function(x){mean(x,na.rm=T)})

avgdf[avgdf$scenario=="decLinearProdshiftCap"&
avgdf$model=="simple"&
avgdf$parameter=="alpha",]

meancvdf<-aggregate(cvdf$x, list(parameter=cvdf$parameter,
                               scenario=cvdf$scenario,
                               method=cvdf$method,
                               model=cvdf$model
                              ), function(x){mean(x,na.rm=T)})



meancv_alpha_p<-ggplot(meancvdf[meancvdf$parameter%in%c("alpha"),])+
geom_bar(aes(x=model,y=x,fill=method),stat="identity",position="dodge")+
facet_wrap(~scenario)+
scale_color_viridis_d(begin=.1, end=.8, option = "A") +
scale_fill_viridis_d(begin=.1, end=.8, option = "A") +
mytheme+
ylab(expression(paste("mean cv for by estimates of log(",alpha,")")))+
#coord_cartesian(ylim = c(0,.5))+ 
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
meancv_alpha_p
ggsave("C:/Users/worc/OneDrive - DFO-MPO/timevar_simeval/summary_stats/mean_by_cv_alpha.png",
    plot=meancv_alpha_p)
  


meancv_smax_p<-ggplot(meancvdf[meancvdf$parameter%in%c("smax"),])+
geom_bar(aes(x=model,y=x,fill=method),stat="identity",position="dodge")+
facet_wrap(~scenario)+
scale_color_viridis_d(begin=.1, end=.8, option = "A") +
scale_fill_viridis_d(begin=.1, end=.8, option = "A") +
mytheme+
ylab(expression(paste("mean cv for by estimates of ",S[max])))+
coord_cartesian(ylim = c(0,.5))+ 
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
meancv_smax_p
ggsave("C:/Users/worc/OneDrive - DFO-MPO/timevar_simeval/summary_stats/mean_by_cv_smax.png",
    plot=meancv_smax_p)




meancv_smsy_p<-ggplot(meancvdf[meancvdf$parameter%in%c("smsy"),])+
geom_bar(aes(x=model,y=x,fill=method),stat="identity",position="dodge")+
facet_wrap(~scenario)+
scale_color_viridis_d(begin=.1, end=.8, option = "A") +
scale_fill_viridis_d(begin=.1, end=.8, option = "A") +
mytheme+
ylab(expression(paste("mean cv for by estimates of ",S[MSY])))+
coord_cartesian(ylim = c(0,1.5))+ 
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
meancv_smsy_p
ggsave("C:/Users/worc/OneDrive - DFO-MPO/timevar_simeval/summary_stats/mean_by_cv_smsy.png",
    plot=meancv_smsy_p)


meancv_sgen_p<-ggplot(meancvdf[meancvdf$parameter%in%c("sgen"),])+
geom_bar(aes(x=model,y=x,fill=method),stat="identity",position="dodge")+
facet_wrap(~scenario)+
scale_color_viridis_d(begin=.1, end=.8, option = "A") +
scale_fill_viridis_d(begin=.1, end=.8, option = "A") +
mytheme+
ylab(expression(paste("mean cv for by estimates of ",S[gen])))+
coord_cartesian(ylim = c(0,.7))+ 
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
meancv_sgen_p
ggsave("C:/Users/worc/OneDrive - DFO-MPO/timevar_simeval/summary_stats/mean_by_cv_sgen.png",
    plot=meancv_sgen_p)

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





#routines to visualize smulation estimation Results
#Catarina Wor
#November 2022
#============================================




library(ggplot2)
library(gridExtra)
source("code/utils.R")
source("code/cluster_func_plots.R")

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
scenNames <- unique(simPar$scenario)

res1<-readRDS(file = "outs/simest/generic/resbase1.rds")
res2<-readRDS(file = "outs/simest/generic/resbase2.rds")


restmb<-rbind(res1,res2)

dim(restmb)
unique(restmb$iteration)

resstan1<-readRDS(file = "outs/simest/generic/resstan1.rds")
resstan2<-readRDS(file = "outs/simest/generic/resstan2.rds")
resstan<-rbind(resstan1,resstan2)
#resstan<-readRDS(file = "outs/simest/generic/resstan.rds")
head(resstan)
unique(resstan$iteration)
resstan
head(res)
res<-rbind(restmb,resstan)

#res<-resstan
res$parameter[res$parameter=="Smax"]<-"smax"

resparam<-res[res$parameter%in%c("alpha","smax","smsy","sgen","umsy"),]


convstat<-aggregate(resparam$convergence,
    list(scenario=resparam$scenario,
        model=resparam$model,
        method=resparam$method,
        iteration=resparam$iteration),
    function(x){sum(x)})
convstatMLE<-convstat[convstat$x==0&convstat$method=="MLE",]
convstatMCMC<-convstat[convstat$x==0&convstat$method=="MCMC",]


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



df<-reshape2::melt(resparam, id.vars=c("parameter","iteration","scenario","method","model","by", "convergence","pbias","bias"))

#df_alpha<-df[df$parameter%in%c("alpha"),]
df$col<-factor(df$variable,levels=c("median","mode", "sim"))


df_alpha_sim<- df[df$parameter=="alpha"&df$variable=="sim",]

df_alpha_est<- df[df$parameter=="alpha"&df$variable=="mode",]

summarydf_alpha<-aggregate(df_alpha_est$value,by=list(scenario=df_alpha_est$scenario, 
    method=df_alpha_est$method, 
    model=df_alpha_est$model,
    by=df_alpha_est$by ),
    function(x) {quantile(x,probs = c(0.025, .5, 0.975))})
summarydf_alpha<-do.call(data.frame, summarydf_alpha)

summarydf_alpha_sim<-aggregate(df_alpha_sim$value,by=list(scenario=df_alpha_est$scenario, 
    method=df_alpha_est$method, 
    model=df_alpha_est$model,
    by=df_alpha_est$by ),
    function(x) {unique(x)})



summarydf_alpha$scenario<-factor(summarydf_alpha$scenario,levels=c("stationary",
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


summarydf_alpha_sim$scenario<-factor(summarydf_alpha_sim$scenario,levels=c("stationary",
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


summarydf_alpha_sim1<-summarydf_alpha_sim[summarydf_alpha_sim$scenario%in%c( "decLinearProd",        
"regimeProd",   "shiftProd",  "sineProd",   "stationary" ),]#&summarydf_alpha_sim$model%in%c( "hmma","hmmb" "hmmab", "rwa",  "rwb", "rwab", "simple","autocorr" ),]

summarydf_alpha1<-summarydf_alpha[summarydf_alpha$scenario%in%c(       "decLinearProd",        
"regimeProd",   "shiftProd",  "sineProd",   "stationary" ),]#&summarydf_alpha$model%in%c( "hmma", "hmmab", "rwa", "rwab", "simple" ),]

ggplot() + 
geom_pointrange(data=summarydf_alpha1,aes(x=by,y= x.50.,ymin = x.2.5., ymax = x.97.5., col=method),alpha=.6)+
geom_line(data=summarydf_alpha_sim1,aes(x=by,y= x),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
coord_cartesian(ylim = c(0.2,2.7))+ 
mytheme+ 
ylab("alpha") +
xlab("year") +
facet_grid(scenario~model, scales="free_y")
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/compareMCMC_MLE_alpha1.png")
#MLE estimates are less biased and higher than MCMC


summarydf_alpha_sim2<-summarydf_alpha_sim[summarydf_alpha_sim$scenario%in%c( "autocorr", "decLinearCap",        
"decLinearProdshiftCap", "regimeCap",  "regimeProdCap", "shiftCap" ),]#&summarydf_alpha_sim$model%in%c( "hmma","hmmb" "hmmab", "rwa",  "rwb", "rwab", "simple","autocorr" ),]

summarydf_alpha2<-summarydf_alpha[summarydf_alpha$scenario%in%c( "autocorr", "decLinearCap",        
"decLinearProdshiftCap", "regimeCap",  "regimeProdCap", "shiftCap" ),]#&summarydf_alpha$model%in%c( "hmma", "hmmab", "rwa", "rwab", "simple" ),]

ggplot() + 
geom_pointrange(data=summarydf_alpha2,aes(x=by,y= x.50.,ymin = x.2.5., ymax = x.97.5., col=method),alpha=.6)+
geom_line(data=summarydf_alpha_sim2,aes(x=by,y= x),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
coord_cartesian(ylim = c(0.2,2.7))+ 
mytheme + 
ylab("alpha") +
xlab("year") +
facet_grid(scenario~model, scales="free_y")
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/compareMCMC_MLE_alpha2.png")
#MLE estimates are less biased and higher than MCMC




#=======================================================
#b estimates


df_smax_sim<- df[df$parameter=="smax"&df$variable=="sim",]
df_smax_est<- df[df$parameter=="smax"&df$variable=="mode",]


summarydf_smax<-aggregate(df_smax_est$value,by=list(scenario=df_smax_est$scenario, 
    method=df_smax_est$method, 
    model=df_smax_est$model,
    by=df_smax_est$by ),
    function(x) {quantile(x,probs = c(0.025, .5, 0.975))})
summarydf_smax<-do.call(data.frame, summarydf_smax)

summarydf_smax_sim<-aggregate(df_smax_sim$value,by=list(scenario=df_smax_est$scenario, 
    method=df_smax_est$method, 
    model=df_smax_est$model,
    by=df_smax_est$by ),
    function(x) {unique(x)})


summarydf_smax$scenario<-factor(summarydf_smax$scenario,levels=c("stationary",
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


summarydf_smax_sim$scenario<-factor(summarydf_smax_sim$scenario,levels=c("stationary",
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




#head(summarydf_smax)
#head(summarydf_alpha_sim)
#unique(summarydf_smax_sim$scenario)


summarydf_smax_sim1<-summarydf_smax_sim[summarydf_smax_sim$scenario%in%c( "decLinearProd",        
"regimeProd",   "shiftProd",  "sineProd",   "stationary"  ),]#&summarydf_smax_sim$model%in%c( "hmma","hmmb", "hmmab", "rwa", "rwb", "rwab", "simple" ),]

summarydf_smax1<-summarydf_smax[summarydf_smax$scenario%in%c("decLinearProd",        
"regimeProd",   "shiftProd",  "sineProd",   "stationary"  ),]#&summarydf_smax$model%in%c( "hmma","hmmb", "hmmab", "rwa", "rwb", "rwab", "simple" ),]
head(summarydf_smax)

ggplot() + 
geom_pointrange(data=summarydf_smax1,aes(x=by,y= x.50.,ymin = x.2.5., ymax = x.97.5., col=method),alpha=.6)+
geom_line(data=summarydf_smax_sim1,aes(x=by,y= x),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
mytheme + 
ylab("Smax") +
xlab("year") +
coord_cartesian(ylim = c(60000,400000))+ 
facet_grid(scenario~model, scales="free_y")
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/compareMCMC_MLE_smax1.png")



summarydf_smax_sim2<-summarydf_smax_sim[summarydf_smax_sim$scenario%in%c( "autocorr", "decLinearCap",        
"decLinearProdshiftCap", "regimeCap",  "regimeProdCap", "shiftCap" ),]#&summarydf_smax_sim$model%in%c( "hmma","hmmb", "hmmab", "rwa", "rwb", "rwab", "simple" ),]

summarydf_smax2<-summarydf_smax[summarydf_smax$scenario%in%c("autocorr", "decLinearCap",        
"decLinearProdshiftCap", "regimeCap",  "regimeProdCap", "shiftCap"),]#&summarydf_smax$model%in%c( "hmma","hmmb", "hmmab", "rwa", "rwb", "rwab", "simple" ),]
head(summarydf_smax2)

ggplot() + 
geom_pointrange(data=summarydf_smax2,aes(x=by,y= x.50.,ymin = x.2.5., ymax = x.97.5., col=method),alpha=.6)+
geom_line(data=summarydf_smax_sim2,aes(x=by,y= x),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
mytheme + 
ylab("Smax") +
xlab("year") +
coord_cartesian(ylim = c(60000,400000))+ 
facet_grid(scenario~model, scales="free_y")
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/compareMCMC_MLE_smax2.png")

dim(summarydf_smax2)



#=======================================================
#smsy estimates


df_smsy_sim<- df[df$parameter=="smsy"&df$variable=="sim",]

df_smsy_est<- df[df$parameter=="smsy"&df$variable=="mode",]

tdf<-df_smsy_est[df_smsy_est$model=="hmmab"&df_smsy_est$scenario=="sineProd"&df_smsy_est$iteration==40,]

sdf<-df_smsy_sim[df_smsy_sim$method=="MLE"&df_smsy_sim$model=="hmmab"&
df_smsy_sim$scenario=="sineProd"&df_smsy_sim$iteration==40,]

#ggplot() + 
#geom_line(data=tdf,aes(x=by,y= value,color=method, group=interaction(iteration,method)), alpha=.6,linewidth=1.2)+
#geom_line(data=sdf,aes(x=by,y= value),color="black", alpha=.6,linewidth=1.2)+
#scale_color_viridis_d(begin=.1, end=.8) +
#scale_fill_viridis_d(begin=.1, end=.8) +
#mytheme


summarydf_smsy<-aggregate(df_smsy_est$value,by=list(scenario=df_smsy_est$scenario, 
    method=df_smsy_est$method, 
    model=df_smsy_est$model,
    by=df_smsy_est$by ),
    function(x) {quantile(x,probs = c(0.025, .5, 0.975))})
summarydf_smsy<-do.call(data.frame, summarydf_smsy)

summary(df_smsy_sim)

summarydf_smsy_sim<-aggregate(df_smsy_sim$value,by=list(scenario=df_smsy_est$scenario, 
    method=df_smsy_est$method, 
    model=df_smsy_est$model,
    by=df_smsy_est$by ),
    function(x) {unique(x)})

summarydf_smsy_sim<-summarydf_smsy_sim[summarydf_smsy_sim$method=="MLE",]

summarydf_smsy<-aggregate(df_smsy_est$value,by=list(scenario=df_smsy_est$scenario, 
    method=df_smsy_est$method, 
    model=df_smsy_est$model,
    by=df_smsy_est$by ),
    function(x) {quantile(x,probs = c(0.025, .5, 0.975))})
summarydf_smsy<-do.call(data.frame, summarydf_smsy)


summarydf_smsy$scenario<-factor(summarydf_smsy$scenario,levels=c("stationary",
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


summarydf_smsy_sim$scenario<-factor(summarydf_smsy_sim$scenario,levels=c("stationary",
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



summarydf_smsy_sim1<-summarydf_smsy_sim[summarydf_smsy_sim$scenario%in%c( "decLinearProd",        
"regimeProd",   "shiftProd",  "sineProd",   "stationary"  ),]#&summarydf_smax_sim$model%in%c( "hmma","hmmb", "hmmab", "rwa", "rwb", "rwab", "simple" ),]


summarydf_smsy1<-summarydf_smsy[summarydf_smsy$scenario%in%c("decLinearProd",        
"regimeProd",   "shiftProd",  "sineProd",   "stationary"  ),]#&summarydf_smax$model%in%c( "hmma","hmmb", "hmmab", "rwa", "rwb", "rwab", "simple" ),]


summarydf_smsy1[summarydf_smsy1$model=="hmmab"&summarydf_smsy1$scenario=="sineProd",]


ggplot() + 
geom_pointrange(data=summarydf_smsy1,aes(x=by,y= x.50.,ymin = x.2.5., ymax = x.97.5., col=method),alpha=.6)+
geom_line(data=summarydf_smsy_sim1,aes(x=by,y= x),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
mytheme + 
ylab("Smsy") +
xlab("year") +
coord_cartesian(ylim = c(20000,150000))+ 
facet_grid(scenario~model, scales="free_y")
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/compareMCMC_MLE_smsy1.png")



summarydf_smsy_sim2<-summarydf_smsy_sim[summarydf_smsy_sim$scenario%in%c( "autocorr", "decLinearCap",        
"decLinearProdshiftCap", "regimeCap",  "regimeProdCap", "shiftCap" ),]#&summarydf_smax_sim$model%in%c( "hmma","hmmb", "hmmab", "rwa", "rwb", "rwab", "simple" ),]

summarydf_smsy2<-summarydf_smsy[summarydf_smsy$scenario%in%c("autocorr", "decLinearCap",        
"decLinearProdshiftCap", "regimeCap",  "regimeProdCap", "shiftCap"),]#&summarydf_smax$model%in%c( "hmma","hmmb", "hmmab", "rwa", "rwb", "rwab", "simple" ),]
head(summarydf_smsy)





ggplot() + 
geom_pointrange(data=summarydf_smsy2,aes(x=by,y= x.50.,ymin = x.2.5., ymax = x.97.5., col=method),alpha=.6)+
geom_line(data=summarydf_smsy_sim2,aes(x=by,y= x),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
mytheme + 
ylab("Smax") +
xlab("year") +
coord_cartesian(ylim = c(20000,150000))+ 
facet_grid(scenario~model, scales="free_y")
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/compareMCMC_MLE_smsy2.png")





#=================================================================================================================
#sensitivity a
simPar <- read.csv("data/sensitivity/SimPars.csv")

## Store relevant object names to help run simulation 
scenNames <- unique(simPar$scenario)

resa1<-readRDS(file = "outs/simest/sensitivity/res_a1.rds")
resa2<-readRDS(file = "outs/simest/sensitivity/res_a2.rds")


restmb<-rbind(resa1,resa2)

head(restmb)

#res<-readRDS(file = "outs/simest/generic/res.rds")
#restmb<-restmb[restmb$convergence==0,]

#

resstana1<-readRDS(file = "outs/simest/sensitivity/resstan_a1.rds")
resstana2<-readRDS(file = "outs/simest/sensitivity/resstan_a2.rds")
resstan<-rbind(resstana1,resstana2)

#resstan<-readRDS(file = "outs/simest/generic/resstan.rds")

res<-rbind(restmb,resstan)

res$parameter[res$parameter=="Smax"]<-"smax"

resparam<-res[res$parameter%in%c("alpha","smax","sigma","smsy","sgen","umsy"),]


convstat<-aggregate(resparam$convergence,
    list(scenario=resparam$scenario,
        model=resparam$model,
        method=resparam$method,
        iteration=resparam$iteration),
    function(x){sum(x)})

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



df<-reshape2::melt(resparam, id.vars=c("parameter","iteration","scenario","method","model","by", "convergence","pbias","bias"))

#df_alpha<-df[df$parameter%in%c("alpha"),]
df$col<-factor(df$variable,levels=c("median","mode", "sim"))

head(df)
df[df$parameter=="umsy",][1:50,]

summarydf<-aggregate(df$value,by=list(scenario=df$scenario,
    parameter=df$parameter,  
    method=df$method, 
    model=df$model,
    by=df$by, 
    variable=df$variable),
    function(x) {quantile(x,probs = c(0.025, .5, 0.975),na.rm=T)})
summarydf<-do.call(data.frame, summarydf)
head(summarydf)

summarydf$scenario<-factor(summarydf$scenario,levels=c( "trendLinearProd1" ,
                                                        "trendLinearProd2", 
                                                        "trendLinearProd5", 
                                                        "trendLinearProd7",
                                                        "trendLinearProd10" ,    
                                                        "regimeProd1",  
                                                        "regimeProd2" ,     
                                                        "regimeProd5",     
                                                        "regimeProd7",
                                                        "regimeProd10"
  ))



summarydf$magnitude_ch<-"a regime 3.7 -> 1.03"
summarydf$magnitude_ch[summarydf$scenario%in%c("trendLinearProd1")]<-"a trend 3.7 -> 1.03"    
summarydf$magnitude_ch[summarydf$scenario%in%c("trendLinearProd10")]<-"a trend 3.7 -> 10"  
summarydf$magnitude_ch[summarydf$scenario%in%c( "regimeProd10")]<-"a regime 3.7 -> 10"  
summarydf$magnitude_ch[summarydf$scenario%in%c("trendLinearProd2" )]<-"a trend 3.7 -> 2"   
summarydf$magnitude_ch[summarydf$scenario%in%c("regimeProd2")]<-"a regime 3.7 -> 2"   
summarydf$magnitude_ch[summarydf$scenario%in%c("trendLinearProd5")]<-"a trend 3.7 -> 5"  
summarydf$magnitude_ch[summarydf$scenario%in%c("regimeProd5")]<-"a regime 3.7 -> 5" 
summarydf$magnitude_ch[summarydf$scenario%in%c("regimeProd7")]<-"a regime 3.7 -> 7" 
summarydf$magnitude_ch[summarydf$scenario%in%c("trendLinearProd7")]<-"a trend 3.7 -> 7"               
summarydf$magnitude_ch<-factor(summarydf$magnitude_ch,levels=c("a trend 3.7 -> 1.03", 
       "a trend 3.7 -> 2",
       "a trend 3.7 -> 5",
       "a trend 3.7 -> 7", 
       "a trend 3.7 -> 10",
       "a regime 3.7 -> 1.03", 
       "a regime 3.7 -> 2",
       "a regime 3.7 -> 5",
       "a regime 3.7 -> 7", 
       "a regime 3.7 -> 10"))

unique(summarydf$magnitude_ch)
unique(summarydf$scenario)
summarydf_alpha_sim<- summarydf[summarydf$parameter=="alpha"&summarydf$variable=="sim",]

summarydf_alpha<- summarydf[summarydf$parameter=="alpha"&summarydf$variable=="mode",]




summarydf_alpha_sim1<-summarydf_alpha_sim[summarydf_alpha_sim$scenario%in%c( "regimeProd1",  
                                                                    "regimeProd2" ,     
                                                                    "regimeProd5",     
                                                                    "regimeProd7",
                                                                    "regimeProd10" ),]
summarydf_alpha1<-summarydf_alpha[summarydf_alpha$scenario%in%c( "regimeProd1",  
                                                                    "regimeProd2" ,     
                                                                    "regimeProd5",     
                                                                    "regimeProd7",
                                                                    "regimeProd10" ),]
head(summarydf_alpha1)
head(summarydf_alpha_sim1)

unique(summarydf_alpha1$magnitude_ch)

ggplot() + 
geom_pointrange(data=summarydf_alpha1,aes(x=by,y= x.50.,ymin = x.2.5., ymax = x.97.5., col=method),alpha=.6)+
geom_line(data=summarydf_alpha_sim1,aes(x=by,y= x.50.),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
coord_cartesian(ylim = c(-3.2,2.5))+ 
mytheme+ 
ylab("alpha") +
xlab("year") +
facet_grid(magnitude_ch~model, scales="free_y")
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/sens_a/compareMCMC_MLE_sensa_alpha1.png")





summarydf_alpha_sim2<-summarydf_alpha_sim[summarydf_alpha_sim$scenario%in%c( "trendLinearProd1" ,
                                                                    "trendLinearProd2", 
                                                                    "trendLinearProd5", 
                                                                    "trendLinearProd7",
                                                                    "trendLinearProd10" ),]
summarydf_alpha2<-summarydf_alpha[summarydf_alpha$scenario%in%c( "trendLinearProd1" ,
                                                                    "trendLinearProd2", 
                                                                    "trendLinearProd5", 
                                                                    "trendLinearProd7",
                                                                    "trendLinearProd10" ),]
  

head(summarydf_alpha2)
head(summarydf_alpha_sim2)  
                                                                
ggplot() + 
geom_pointrange(data=summarydf_alpha2,aes(x=by,y= x.50.,ymin = x.2.5., ymax = x.97.5., col=method),alpha=.6)+
geom_line(data=summarydf_alpha_sim2,aes(x=by,y= x.50.),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
coord_cartesian(ylim = c(-0.4,2.7))+ 
mytheme + 
ylab("alpha") +
xlab("year") +
facet_grid(magnitude_ch~model, scales="free_y")
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/sens_a/compareMCMC_MLE_sensa_alpha2.png")
#MLE estimates are less biased and higher than MCMC


#smax



summarydf_smax_sim<- summarydf[summarydf$parameter=="smax"&summarydf$variable=="sim",]

summarydf_smax<- summarydf[summarydf$parameter=="smax"&summarydf$variable=="mode",]




summarydf_smax_sim1<-summarydf_smax_sim[summarydf_smax_sim$scenario%in%c( "regimeProd1",  
                                                                    "regimeProd2" ,     
                                                                    "regimeProd5",     
                                                                    "regimeProd7",
                                                                    "regimeProd10" ),]
summarydf_smax1<-summarydf_smax[summarydf_smax$scenario%in%c( "regimeProd1",  
                                                                    "regimeProd2" ,     
                                                                    "regimeProd5",     
                                                                    "regimeProd7",
                                                                    "regimeProd10" ),]
head(summarydf_smax1)

ggplot() + 
geom_pointrange(data=summarydf_smax1,aes(x=by,y= x.50.,ymin = x.2.5., ymax = x.97.5., col=method),alpha=.6)+
geom_line(data=summarydf_smax_sim1,aes(x=by,y= x.50.),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
coord_cartesian(ylim = c(10000,370000))+ 
mytheme+ 
ylab("smax") +
xlab("year") +
facet_grid(magnitude_ch~model, scales="free_y")
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/sens_a/compareMCMC_MLE_sensa_smax1.png")



summarydf_smax_sim2<-summarydf_smax_sim[summarydf_smax_sim$scenario%in%c( "trendLinearProd1" ,
                                                                    "trendLinearProd2", 
                                                                    "trendLinearProd5", 
                                                                    "trendLinearProd7",
                                                                    "trendLinearProd10" ),]
summarydf_smax2<-summarydf_smax[summarydf_smax$scenario%in%c( "trendLinearProd1" ,
                                                                    "trendLinearProd2", 
                                                                    "trendLinearProd5", 
                                                                    "trendLinearProd7",
                                                                    "trendLinearProd10" ),]

head(summarydf_smax2)
head(summarydf_smax_sim2)
ggplot() + 
geom_pointrange(data=summarydf_smax2,aes(x=by,y= x.50.,ymin = x.2.5., ymax = x.97.5., col=method),alpha=.6)+
geom_line(data=summarydf_smax_sim2,aes(x=by,y= x.50.),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
coord_cartesian(ylim = c(10000,470000))+ 
mytheme + 
ylab("smax") +
xlab("year") +
facet_grid(magnitude_ch~model, scales="free_y")
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/sens_a/compareMCMC_MLE_sensa_smax2.png")
#MLE estimates are less biased and higher than MCMC



#=======================================================================
#smsy


summarydf_smsy_sim<- summarydf[summarydf$parameter=="smsy"&summarydf$variable=="sim",]

summarydf_smsy<- summarydf[summarydf$parameter=="smsy"&summarydf$variable=="mode",]




summarydf_smsy_sim1<-summarydf_smsy_sim[summarydf_smsy_sim$scenario%in%c( "regimeProd1",  
                                                                    "regimeProd2" ,     
                                                                    "regimeProd5",     
                                                                    "regimeProd7",
                                                                    "regimeProd10" ),]
summarydf_smsy1<-summarydf_smsy[summarydf_smsy$scenario%in%c( "regimeProd1",  
                                                                    "regimeProd2" ,     
                                                                    "regimeProd5",     
                                                                    "regimeProd7",
                                                                    "regimeProd10" ),]
head(summarydf_smsy1)

ggplot() + 
geom_pointrange(data=summarydf_smsy1,aes(x=by,y= x.50.,ymin = x.2.5., ymax = x.97.5., col=method),alpha=.6)+
geom_line(data=summarydf_smsy_sim1,aes(x=by,y= x.50.),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
coord_cartesian(ylim = c(-100000,170000))+ 
mytheme+ 
ylab("smsy") +
xlab("year") +
facet_grid(magnitude_ch~model, scales="free_y")
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/sens_a/compareMCMC_MLE_sensa_smsy1.png")



summarydf_smsy_sim2<-summarydf_smax_sim[summarydf_smsy_sim$scenario%in%c( "trendLinearProd1" ,
                                                                    "trendLinearProd2", 
                                                                    "trendLinearProd5", 
                                                                    "trendLinearProd7",
                                                                    "trendLinearProd10" ),]
summarydf_smsy2<-summarydf_smax[summarydf_smsy$scenario%in%c( "trendLinearProd1" ,
                                                                    "trendLinearProd2", 
                                                                    "trendLinearProd5", 
                                                                    "trendLinearProd7",
                                                                    "trendLinearProd10" ),]

head(summarydf_smsy2)
head(summarydf_smsy_sim2)
ggplot() + 
geom_pointrange(data=summarydf_smsy2,aes(x=by,y= x.50.,ymin = x.2.5., ymax = x.97.5., col=method),alpha=.6)+
geom_line(data=summarydf_smsy_sim2,aes(x=by,y= x.50.),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
coord_cartesian(ylim = c(10000,470000))+ 
mytheme + 
ylab("smsy") +
xlab("year") +
facet_grid(scenario~model, scales="free_y")
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/sens_a/compareMCMC_MLE_sensa_smsy2.png")
#MLE estimates are less biased and higher than MCMC





#=======================================================================
#umsy -- these are wrong need to be checked



summarydf_umsy_sim<- summarydf[summarydf$parameter=="umsy"& 
                              summarydf$variable=="sim"&
                              summarydf$method=="MCMC",]

summarydf_umsy<- summarydf[summarydf$parameter=="umsy"&summarydf$variable=="mode",]


summarydf_umsy_sim[summarydf_umsy_sim$method=="MCMC"&
                   summarydf_umsy_sim$scenario=="regimeProd1"&
                   summarydf_umsy_sim$model=="rwa",]


summarydf_umsy_sim1<-summarydf_umsy_sim[summarydf_umsy_sim$scenario%in%c( "regimeProd1",  
                                                                    "regimeProd2" ,     
                                                                    "regimeProd5",     
                                                                    "regimeProd7",
                                                                    "regimeProd10" ),]
summarydf_umsy1<-summarydf_umsy[summarydf_umsy$scenario%in%c( "regimeProd1",  
                                                                    "regimeProd2" ,     
                                                                    "regimeProd5",     
                                                                    "regimeProd7",
                                                                    "regimeProd10" ),]
head(summarydf_umsy1)

head(summarydf_umsy_sim)
summarydf_umsy_sim1[300:330,]
summarydf_umsy1[summarydf_umsy1$scenario=="regimeProd1"&
                summarydf_umsy1$method=="MCMC"&
                summarydf_umsy1$model=="hmma",]

summarydf_alpha1[summarydf_alpha1$scenario=="regimeProd1"&
                summarydf_alpha1$method=="MCMC"&
                summarydf_alpha1$model=="hmma",]



summarydf_umsy_sim1[summarydf_umsy_sim1$scenario=="regimeProd1"&
                summarydf_umsy_sim1$method=="MCMC"&
                summarydf_umsy_sim1$model=="hmma",]

ggplot() + 
geom_pointrange(data=summarydf_umsy1,aes(x=by,y= x.50.,ymin = x.2.5., ymax = x.97.5., col=method),alpha=.6)+
geom_line(data=summarydf_umsy_sim1,aes(x=by,y= x.50.),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
coord_cartesian(ylim = c(-.3,1))+ 
mytheme+ 
ylab("umsy") +
xlab("year") +
facet_grid(scenario~model, scales="free_y")
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/sens_a/compareMCMC_MLE_sensa_smsy1.png")



summarydf_smsy_sim2<-summarydf_smax_sim[summarydf_smsy_sim$scenario%in%c( "trendLinearProd1" ,
                                                                    "trendLinearProd2", 
                                                                    "trendLinearProd5", 
                                                                    "trendLinearProd7",
                                                                    "trendLinearProd10" ),]
summarydf_smsy2<-summarydf_smax[summarydf_smsy$scenario%in%c( "trendLinearProd1" ,
                                                                    "trendLinearProd2", 
                                                                    "trendLinearProd5", 
                                                                    "trendLinearProd7",
                                                                    "trendLinearProd10" ),]

head(summarydf_smsy2)
head(summarydf_smsy_sim2)
ggplot() + 
geom_pointrange(data=summarydf_smsy2,aes(x=by,y= x.50.,ymin = x.2.5., ymax = x.97.5., col=method),alpha=.6)+
geom_line(data=summarydf_smsy_sim2,aes(x=by,y= x.50.),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
coord_cartesian(ylim = c(10000,470000))+ 
mytheme + 
ylab("alpha") +
xlab("year") +
facet_grid(magnitude_ch~model, scales="free_y")
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/sens_a/compareMCMC_MLE_sensa_smsy2.png")
#MLE estimates are less biased and higher than MCMC



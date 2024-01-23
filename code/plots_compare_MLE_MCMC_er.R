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
    theme_classic(16)+
        theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
              legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="bold", size=12),
              axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
)




#========================================================================================================
#base case
#read in data
simPar <- read.csv("data/genericER/SimPars_ER.csv")

## Store relevant object names to help run simulation 
scenNames <- unique(simPar$scenario)

res1<-readRDS(file = "outs/simest/genericER/res_erq1.rds")
res2<-readRDS(file = "outs/simest/genericER/res_erq2.rds")
res3<-readRDS(file = "outs/simest/genericER/res_erq3.rds")
res4<-readRDS(file = "outs/simest/genericER/res_erq4.rds")



restmb<-rbind(res1,res2,res3,res4)

head(restmb)
unique(restmb$iteration)

resstan1<-readRDS(file = "outs/simest/genericER/resstan_baseER1.rds")
resstan2<-readRDS(file = "outs/simest/genericER/resstan_baseER2.rds")
resstan<-rbind(resstan1,resstan2)
#resstan<-readRDS(file = "outs/simest/generic/resstan.rds")

res<-rbind(restmb,resstan)

#res<-resstan
res$parameter[res$parameter=="Smax"]<-"smax"
res$method[res$method=="MCMC"]<-"HMC"
resparam<-res[res$parameter%in%c("alpha","smax","smsy","sgen","umsy"),]


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



df<-reshape2::melt(resparam, id.vars=c("parameter","iteration","scenario","method","model","by", "convergence","pbias","bias"))

#df_alpha<-df[df$parameter%in%c("alpha"),]
df$col<-factor(df$variable,levels=c("median","mode", "sim"))

unique(df$scenario)
df$scenario<-factor(df$scenario,levels=c("decLinearProd_highERLowError",
                                         "decLinearProd_ShiftERLowError", 
                                         "incLinearProd_highERLowError",  
                                         "incLinearProd_ShiftERLowError", 
                                         "highERHighError",               
                                         "highERLowError",               
                                         "lowERHighError",                
                                         "lowERLowError",                
                                          "ShiftERHighError",              
                                          "ShiftERLowError",               
                                          "trendERHighError",              
                                          "trendERLowError" )   


df$ERtrend<-case_match(df$scenario,"decLinearProd_highERLowError"~"highER",
                                   "decLinearProd_ShiftERLowError"~"ShiftER", 
                                   "incLinearProd_highERLowError"~"highER",  
                                   "incLinearProd_ShiftERLowError"~"ShiftER", 
                                   "highERHighError"~"highER",               
                                   "highERLowError"~"highER",               
                                   "lowERHighError"~"lowER",                
                                   "lowERLowError"~"lowER",                
                                    "ShiftERHighError"~"ShiftER",              
                                    "ShiftERLowError"~"ShiftER",               
                                    "trendERHighError"~"trendER",              
                                    "trendERLowError"~"trendER")



df$ERerror<-case_match(df$scenario,"decLinearProd_highERLowError"~"LowError",
                                   "decLinearProd_ShiftERLowError"~"LowError", 
                                   "incLinearProd_highERLowError"~"LowError",  
                                   "incLinearProd_ShiftERLowError"~"LowError", 
                                   "highERHighError"~"HighError",               
                                   "highERLowError"~"LowError",               
                                   "lowERHighError"~"HighError",                
                                   "lowERLowError"~"LowError",                
                                    "ShiftERHighError"~"HighError",              
                                    "ShiftERLowError"~"LowError",               
                                    "trendERHighError"~"HighError",              
                                    "trendERLowError"~"LowError")


df$dynamics<-case_match(df$scenario,"decLinearProd_highERLowError"~"decLinear",
                                   "decLinearProd_ShiftERLowError"~"decLinear", 
                                   "incLinearProd_highERLowError"~"incLinear",  
                                   "incLinearProd_ShiftERLowError"~"incLinear", 
                                   "highERHighError"~"stationary",               
                                   "highERLowError"~"stationary",               
                                   "lowERHighError"~"stationary",                
                                   "lowERLowError"~"stationary",                
                                    "ShiftERHighError"~"stationary",              
                                    "ShiftERLowError"~"stationary",               
                                    "trendERHighError"~"stationary",              
                                    "trendERLowError"~"stationary")
df$model<-factor(df$model,levels=c("simple",
                                   "autocorr", 
                                   "rwa",
                                   "hmma",
                                   "rwb",
                                   "hmmb",
                                   "rwab",
                                   "hmmab"  ))


df_alpha_sim<- df[df$parameter=="alpha"&df$variable=="sim",]

df_alpha_est<- df[df$parameter=="alpha"&df$variable=="mode",]

head(df_alpha_est)
summarydf_alpha<-aggregate(df_alpha_est$value,by=list(scenario=df_alpha_est$scenario, 
    method=df_alpha_est$method, 
    model=df_alpha_est$model,
    by=df_alpha_est$by ,
    dynamics=df_alpha_est$dynamics,
    ERerror=df_alpha_est$ERerror,
    ERtrend=df_alpha_est$ERtrend),
    function(x) {quantile(x,probs = c(0.025, .5, 0.975))})
summarydf_alpha<-do.call(data.frame, summarydf_alpha)

summarydf_alpha_sim<-aggregate(df_alpha_sim$value,by=list(scenario=df_alpha_sim$scenario, 
    method=df_alpha_sim$method, 
    model=df_alpha_sim$model,
    by=df_alpha_sim$by ,
    dynamics=df_alpha_sim$dynamics,
    ERerror=df_alpha_sim$ERerror,
    ERtrend=df_alpha_sim$ERtrend),
    function(x) {unique(x)})

head(df_alpha_sim)
unique(df_alpha_sim$scenario)

summarydf_alpha_sim1<-summarydf_alpha_sim[summarydf_alpha_sim$scenario%in%c("decLinearProd_highERLowError",
                                         "decLinearProd_ShiftERLowError", 
                                         "incLinearProd_highERLowError",  
                                         "incLinearProd_ShiftERLowError"
                                            ),]

summarydf_alpha1<-summarydf_alpha[summarydf_alpha$scenario%in%c("decLinearProd_highERLowError",
                                         "decLinearProd_ShiftERLowError", 
                                         "incLinearProd_highERLowError",  
                                         "incLinearProd_ShiftERLowError"
                                          ),]
scenlab1<-c("highER","ShiftER", "highER","ShiftER")
names(scenlab1) <- c("decLinearProd_highERLowError",
                                         "decLinearProd_ShiftERLowError", 
                                         "incLinearProd_highERLowError",  
                                         "incLinearProd_ShiftERLowError")

head(summarydf_alpha)

er_alpha1<-ggplot() + 
geom_pointrange(data=summarydf_alpha1,aes(x=by,y= x.50.,ymin = x.2.5., ymax = x.97.5., col=method),alpha=.6)+
geom_line(data=summarydf_alpha_sim1,aes(x=by,y= x),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
coord_cartesian(ylim = c(0.2,2.7))+ 
mytheme+ 
ylab("alpha") +
xlab("year") +
facet_grid(dynamics+scenario~model, scales="free_y",labeller = labeller(scenario= scenlab1))
er_alpha1
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/ER/compareMCMC_MLE_alpha1.png")
#MLE estimates are less biased and higher than MCMC


summarydf_alpha_sim2<-summarydf_alpha_sim[summarydf_alpha_sim$scenario%in%c(    "highERHighError",                              
                                         "lowERHighError",                           
                                          "ShiftERHighError",                                                                   
                                          "trendERHighError"),]#&summarydf_alpha_sim$model%in%c( "hmma","hmmb" "hmmab", "rwa",  "rwb", "rwab", "simple","autocorr" ),]

summarydf_alpha2<-summarydf_alpha[summarydf_alpha$scenario%in%c( "highERHighError",                              
                                         "lowERHighError",                           
                                          "ShiftERHighError",                                                                   
                                          "trendERHighError"),]#&summarydf_alpha$model%in%c( "hmma", "hmmab", "rwa", "rwab", "simple" ),]

                                       ),]
scenlab2<-c("highER","lowER", "ShiftER","trendER")
names(scenlab2) <- c( "highERHighError",                              
                                         "lowERHighError",                           
                                          "ShiftERHighError",                                                                   
                                          "trendERHighError")

er_alpha2<-ggplot() + 
geom_pointrange(data=summarydf_alpha2,aes(x=by,y= x.50.,ymin = x.2.5., ymax = x.97.5., col=method),alpha=.6)+
geom_line(data=summarydf_alpha_sim2,aes(x=by,y= x),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
coord_cartesian(ylim = c(0.2,2.7))+ 
mytheme + 
ylab("alpha") +
xlab("year") +
facet_grid(dynamics+scenario~model, scales="free_y",labeller = labeller(scenario= scenlab2))
er_alpha2
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/ER/compareMCMC_MLE_alpha2.png",
    plot=er_alpha2)
#MLE estimates are less biased and higher than MCMC


summarydf_alpha_sim3<-summarydf_alpha_sim[summarydf_alpha_sim$scenario%in%c(               
                                         "highERLowError", 
                                         "lowERLowError",  
                                          "ShiftERLowError",      
                                          "trendERLowError"),]#&summarydf_alpha_sim$model%in%c( "hmma","hmmb" "hmmab", "rwa",  "rwb", "rwab", "simple","autocorr" ),]

summarydf_alpha3<-summarydf_alpha[summarydf_alpha$scenario%in%c(  "highERLowError", 
                                         "lowERLowError",  
                                          "ShiftERLowError",      
                                          "trendERLowError" ),]#&summarydf_alpha$model%in%c( "hmma", "hmmab", "rwa", "rwab", "simple" ),]


                                     
scenlab3<-c("highER","lowER", "ShiftER","trendER")
names(scenlab3) <- c( "highERLowError", 
                                         "lowERLowError",  
                                          "ShiftERLowError",      
                                          "trendERLowError" )



er_alpha3<-ggplot() + 
geom_pointrange(data=summarydf_alpha3,aes(x=by,y= x.50.,ymin = x.2.5., ymax = x.97.5., col=method),alpha=.6)+
geom_line(data=summarydf_alpha_sim3,aes(x=by,y= x),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
coord_cartesian(ylim = c(0.2,2.7))+ 
mytheme + 
ylab("alpha") +
xlab("year") +
facet_grid(dynamics+scenario~model, scales="free_y",labeller = labeller(scenario= scenlab3))
er_alpha3
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/ER/compareMCMC_MLE_alpha3.png",
    plot=er_alpha3)


#=======================================================
#b estimates


df_smax_sim<- df[df$parameter=="smax"&df$variable=="sim",]
df_smax_est<- df[df$parameter=="smax"&df$variable=="mode",]


summarydf_smax<-aggregate(df_smax_est$value,by=list(scenario=df_smax_est$scenario, 
    method=df_smax_est$method, 
    model=df_smax_est$model,
    by=df_smax_est$by,
    dynamics=df_smax_est$dynamics,
    ERerror=df_smax_est$ERerror,
    ERtrend=df_smax_est$ERtrend ),
    function(x) {quantile(x,probs = c(0.025, .5, 0.975))})
summarydf_smax<-do.call(data.frame, summarydf_smax)

summarydf_smax_sim<-aggregate(df_smax_sim$value,by=list(scenario=df_smax_sim$scenario, 
    method=df_smax_sim$method, 
    model=df_smax_sim$model,
    by=df_smax_sim$by,
     dynamics=df_smax_sim$dynamics,
    ERerror=df_smax_sim$ERerror,
    ERtrend=df_smax_sim$ERtrend ),
    function(x) {unique(x)})




summarydf_smax_sim1<-summarydf_smax_sim[summarydf_smax_sim$scenario%in%c( "decLinearProd_highERLowError",
                                         "decLinearProd_ShiftERLowError", 
                                         "incLinearProd_highERLowError",  
                                         "incLinearProd_ShiftERLowError" ),]

summarydf_smax1<-summarydf_smax[summarydf_smax$scenario%in%c("decLinearProd_highERLowError",
                                         "decLinearProd_ShiftERLowError", 
                                         "incLinearProd_highERLowError",  
                                         "incLinearProd_ShiftERLowError" ),]


er_smax1<-ggplot() + 
geom_pointrange(data=summarydf_smax1,aes(x=by,y= x.50.,ymin = x.2.5., ymax = x.97.5., col=method),alpha=.6)+
geom_line(data=summarydf_smax_sim1,aes(x=by,y= x),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
mytheme + 
ylab("Smax") +
xlab("year") +
coord_cartesian(ylim = c(60000,400000))+ 
facet_grid(dynamics+scenario~model, scales="free_y",labeller = labeller(scenario= scenlab1))
er_smax1
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/ER/compareMCMC_MLE_smax1.png",
    plot=er_smax1)



summarydf_smax_sim2<-summarydf_smax_sim[summarydf_smax_sim$scenario%in%c(  "highERHighError",                              
                                         "lowERHighError",                           
                                          "ShiftERHighError",                                                                   
                                          "trendERHighError" ),]

summarydf_smax2<-summarydf_smax[summarydf_smax$scenario%in%c( "highERHighError",                              
                                         "lowERHighError",                           
                                          "ShiftERHighError",                                                                   
                                          "trendERHighError" ),]#&summarydf_smax$model%in%c( "hmma","hmmb", "hmmab", "rwa", "rwb", "rwab", "simple" ),]
head(summarydf_smax2)

er_smax2<-ggplot() + 
geom_pointrange(data=summarydf_smax2,aes(x=by,y= x.50.,ymin = x.2.5., ymax = x.97.5., col=method),alpha=.6)+
geom_line(data=summarydf_smax_sim2,aes(x=by,y= x),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
mytheme + 
ylab("Smax") +
xlab("year") +
coord_cartesian(ylim = c(60000,400000))+ 
facet_grid(dynamics+scenario~model, scales="free_y",labeller = labeller(scenario= scenlab2))
er_smax2
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/ER/compareMCMC_MLE_smax2.png")




summarydf_smax_sim3<-summarydf_smax_sim[summarydf_smax_sim$scenario%in%c(  "highERHighError",                              
                                         "lowERHighError",                           
                                          "ShiftERHighError",                                                                   
                                          "trendERHighError" ),]

summarydf_smax3<-summarydf_smax[summarydf_smax$scenario%in%c( "highERHighError",                              
                                         "lowERHighError",                           
                                          "ShiftERHighError",                                                                   
                                          "trendERHighError" ),]#&summarydf_smax$model%in%c( "hmma","hmmb", "hmmab", "rwa", "rwb", "rwab", "simple" ),]


er_smax3<-ggplot() + 
geom_pointrange(data=summarydf_smax3,aes(x=by,y= x.50.,ymin = x.2.5., ymax = x.97.5., col=method),alpha=.6)+
geom_line(data=summarydf_smax_sim3,aes(x=by,y= x),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
mytheme + 
ylab("Smax") +
xlab("year") +
coord_cartesian(ylim = c(60000,400000))+ 
facet_grid(dynamics+scenario~model, scales="free_y",labeller = labeller(scenario= scenlab2))
er_smax3
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/ER/compareMCMC_MLE_smax2.png")


#=======================================================
#smsy estimates


df_smsy_sim<- df[df$parameter=="smsy"&df$variable=="sim",]

df_smsy_est<- df[df$parameter=="smsy"&df$variable=="mode",]

df<-df_smsy_est[df_smsy_est$model=="hmmab"&df_smsy_est$scenario=="sineProd"&df_smsy_est$iteration==40,]

df<-df_smsy_sim[df_smsy_sim$method=="MLE"&df_smsy_sim$model=="hmmab"&
df_smsy_sim$scenario=="sineProd"&df_smsy_sim$iteration==40,]


summarydf_smsy<-aggregate(df_smsy_est$value,by=list(scenario=df_smsy_est$scenario, 
    method=df_smsy_est$method, 
    model=df_smsy_est$model,
    by=df_smsy_est$by,
     dynamics=df_smsy_est$dynamics,
    ERerror=df_smsy_est$ERerror,
    ERtrend=df_smsy_est$ERtrend ),
    function(x) {quantile(x,probs = c(0.025, .5, 0.975))})
summarydf_smsy<-do.call(data.frame, summarydf_smsy)

summary(df_smsy_sim)

summarydf_smsy_sim<-aggregate(df_smsy_sim$value,by=list(scenario=df_smsy_sim$scenario, 
    method=df_smsy_sim$method, 
    model=df_smsy_sim$model,
    by=df_smsy_sim$by,
    dynamics=df_smsy_sim$dynamics,
    ERerror=df_smsy_sim$ERerror,
    ERtrend=df_smsy_sim$ERtrend ),
    function(x) {unique(x)})




summarydf_smsy_sim1<-summarydf_smsy_sim[summarydf_smsy_sim$scenario%in%c( "decLinearProd_highERLowError",
                                         "decLinearProd_ShiftERLowError", 
                                         "incLinearProd_highERLowError",  
                                         "incLinearProd_ShiftERLowError" ),]

summarydf_smsy1<-summarydf_smsy[summarydf_smsy$scenario%in%c("decLinearProd_highERLowError",
                                         "decLinearProd_ShiftERLowError", 
                                         "incLinearProd_highERLowError",  
                                         "incLinearProd_ShiftERLowError" ),]#&summarydf_smax$model%in%c( "hmma","hmmb", "hmmab", "rwa", "rwb", "rwab", "simple" ),]



er_smsy1<-ggplot() + 
geom_pointrange(data=summarydf_smsy1,aes(x=by,y= x.50.,ymin = x.2.5., ymax = x.97.5., col=method),alpha=.6)+
geom_line(data=summarydf_smsy_sim1,aes(x=by,y= x),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
mytheme + 
ylab("Smsy") +
xlab("year") +
coord_cartesian(ylim = c(20000,150000))+ 
facet_grid(dynamics+scenario~model, scales="free_y",labeller = labeller(scenario= scenlab1))
er_smsy1
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/base/compareMCMC_MLE_smsy1.png",
    plot=er_smsy1)



summarydf_smsy_sim2<-summarydf_smsy_sim[summarydf_smsy_sim$scenario%in%c( "highERHighError",                              
                                         "lowERHighError",                           
                                          "ShiftERHighError",                                                                   
                                          "trendERHighError"  ),]

summarydf_smsy2<-summarydf_smsy[summarydf_smsy$scenario%in%c( "highERHighError",                              
                                         "lowERHighError",                           
                                          "ShiftERHighError",                                                                   
                                          "trendERHighError" ),]



head(summarydf_smsy2)

er_smsy2 <- ggplot() + 
geom_pointrange(data=summarydf_smsy2,aes(x=by,y= x.50.,ymin = x.2.5., ymax = x.97.5., col=method),alpha=.6)+
geom_line(data=summarydf_smsy_sim2,aes(x=by,y= x),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
mytheme + 
ylab("Smax") +
xlab("year") +
coord_cartesian(ylim = c(20000,150000))+ 
facet_grid(dynamics+scenario~model, scales="free_y",labeller = labeller(scenario= scenlab2))
er_smsy2
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/ER/compareMCMC_MLE_smsy2.png")



summarydf_smsy_sim3<-summarydf_smsy_sim[summarydf_smsy_sim$scenario%in%c( "highERHighError",                              
                                         "lowERHighError",                           
                                          "ShiftERHighError",                                                                   
                                          "trendERHighError"  ),]

summarydf_smsy3<-summarydf_smsy[summarydf_smsy$scenario%in%c( "highERHighError",                              
                                         "lowERHighError",                           
                                          "ShiftERHighError",                                                                   
                                          "trendERHighError" ),]



head(summarydf_smsy2)

er_smsy3 <- ggplot() + 
geom_pointrange(data=summarydf_smsy3,aes(x=by,y= x.50.,ymin = x.2.5., ymax = x.97.5., col=method),alpha=.6)+
geom_line(data=summarydf_smsy_sim3,aes(x=by,y= x),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
mytheme + 
ylab("Smax") +
xlab("year") +
coord_cartesian(ylim = c(20000,150000))+ 
facet_grid(dynamics+scenario~model, scales="free_y",labeller = labeller(scenario= scenlab2))
er_smsy3
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/ER/compareMCMC_MLE_smsy3.png")




#-----------------------------------------------------------------------------------
#summary plots 
unique(summarydf_smsy_sim$scenario)


summarydf_smsy_sim_redux<-summarydf_smsy_sim[summarydf_smsy_sim$scenario%in%c("highERLowError",  "ShiftERLowError", 
    "decLinearProd_highERLowError",
  "incLinearProd_highERLowError", "decLinearProd_ShiftERLowError", "incLinearProd_ShiftERLowError")&
summarydf_smsy_sim$model%in%c("autocorr", "rwa")&summarydf_smsy_sim$method=="MLE",]


summarydf_smsy_redux<-summarydf_smsy[summarydf_smsy$scenario%in%c( "highERLowError",  "ShiftERLowError", 
    "decLinearProd_highERLowError",
  "incLinearProd_highERLowError","decLinearProd_ShiftERLowError", "incLinearProd_ShiftERLowError" )&
  summarydf_smsy$model%in%c("autocorr", "rwa") &summarydf_smsy$method=="MLE",]



head(summarydf_smsy)
head(summarydf_smsy_redux_alpha)

psmsy_highERscn_line<-ggplot() + 
geom_pointrange(data=summarydf_smsy_redux,aes(x=by,y= x.50.,ymin = x.2.5., ymax = x.97.5., color=model),alpha=.9)+
geom_line(data=summarydf_smsy_sim_redux,aes(x=by,y= x),color="black", alpha=.8,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8,option = "E") +
scale_fill_viridis_d(begin=.1, end=.8,option = "E") +
coord_cartesian(ylim = c(20000,180000))+ 
mytheme+ 
ylab("Smsy") +
xlab("year") +
#theme( strip.text.y = element_blank(),strip.text.x = element_blank())+
facet_grid(dynamics+ERtrend ~model, scales="free_y")
psmsy_highERscn_line


head(df$scenario)
df_smsy_est_redux<- df[df$parameter=="smsy"&df$variable=="mode"&
df$scenario%in%c("highERLowError",  "ShiftERLowError", 
    "decLinearProd_highERLowError",
  "incLinearProd_highERLowError","decLinearProd_ShiftERLowError", "incLinearProd_ShiftERLowError")&
df$model%in%c("autocorr", "rwa")&
df$method=="MLE",]



head(df_smsy_est_redux)

psmsy_erscn_violin_abs<-ggplot(df_smsy_est_redux) + 
geom_violin(aes(x=model,y=abs(bias), fill=model), scale="width", trim=TRUE, alpha=.7,adjust = 1.8)+
geom_boxplot(aes(x=model,y=abs(bias), fill=model),outlier.shape = NA, width=0.1)+
geom_hline(yintercept=0,color="black", alpha=.6,linewidth=1.2)+
 scale_fill_viridis_d(begin=.1, end=.8,,option = "E") +
 facet_grid(dynamics ~ERtrend, scales="free_y")+
 ylab("absolute bias in Smsy") +
 coord_cartesian(ylim = c(0,75000))+ 
mytheme
psmsy_erscn_violin_abs

multi.page.abs.smsy_erscn <- ggarrange(psmsy_highERscn_line, psmsy_erscn_violin_abs,
                        nrow = 1, ncol = 2,
                        common.legend = TRUE,
                        legend="bottom")
multi.page.abs.smsy_erscn

ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/summary_figs/bias_trends_violins_smsy_ERscn.png",
    plot=multi.page.abs.smsy_erscn)


#=======================================================================================================
#summary plots 
#with ow ER scenarios
unique(summarydf_smsy_sim$scenario)


summarydf_smsy_sim_redux2<-summarydf_smsy_sim[summarydf_smsy_sim$scenario%in%c("highERLowError",  "ShiftERLowError", "lowERLowError",
    "decLinearProd_highERLowError",
  "incLinearProd_highERLowError", "decLinearProd_ShiftERLowError", "incLinearProd_ShiftERLowError")&
summarydf_smsy_sim$model%in%c("autocorr", "rwa")&summarydf_smsy_sim$method=="MLE",]


summarydf_smsy_redux2<-summarydf_smsy[summarydf_smsy$scenario%in%c( "highERLowError",  "ShiftERLowError", "lowERLowError",
    "decLinearProd_highERLowError","incLinearProd_highERLowError","decLinearProd_ShiftERLowError", "incLinearProd_ShiftERLowError" )&
  summarydf_smsy$model%in%c("autocorr", "rwa") &summarydf_smsy$method=="MLE",]



psmsy_highERscn_line2<-ggplot() + 
geom_pointrange(data=summarydf_smsy_redux2,aes(x=by,y= x.50.,ymin = x.2.5., ymax = x.97.5., color=model),alpha=.9)+
geom_line(data=summarydf_smsy_sim_redux2,aes(x=by,y= x),color="black", alpha=.8,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8,option = "E") +
scale_fill_viridis_d(begin=.1, end=.8,option = "E") +
coord_cartesian(ylim = c(20000,180000))+ 
mytheme+ 
ylab("Smsy") +
xlab("year") +
#theme( strip.text.y = element_blank(),strip.text.x = element_blank())+
facet_grid(dynamics+ERtrend ~model, scales="free_y")
psmsy_highERscn_line2


head(df$scenario)
df_smsy_est_redux2<- df[df$parameter=="smsy"&df$variable=="mode"&
df$scenario%in%c("highERLowError",  "ShiftERLowError", "lowERLowError","decLinearProd_highERLowError", 
    "incLinearProd_highERLowError","decLinearProd_ShiftERLowError", "incLinearProd_ShiftERLowError")&
df$model%in%c("autocorr", "rwa")&
df$method=="MLE",]



head(df_smsy_est_redux)

psmsy_erscn_violin_abs2<-ggplot(df_smsy_est_redux2) + 
geom_violin(aes(x=model,y=abs(bias), fill=model), scale="width", trim=TRUE, alpha=.7,adjust = 1.8)+
geom_boxplot(aes(x=model,y=abs(bias), fill=model),outlier.shape = NA, width=0.1)+
geom_hline(yintercept=0,color="black", alpha=.6,linewidth=1.2)+
 scale_fill_viridis_d(begin=.1, end=.8,,option = "E") +
 facet_grid(dynamics ~ERtrend, scales="free_y")+
 ylab("absolute bias in Smsy") +
 coord_cartesian(ylim = c(0,75000))+ 
mytheme
psmsy_erscn_violin_abs2

multi.page.abs.smsy_erscn2 <- ggarrange(psmsy_highERscn_line2, psmsy_erscn_violin_abs2,
                        nrow = 1, ncol = 2,
                        common.legend = TRUE,
                        legend="bottom")
multi.page.abs.smsy_erscn2
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/summary_figs/bias_trends_violins_smsy_ERscn_wlowER.png",
    plot=multi.page.abs.smsy_erscn2)


#=============================
#Smax


summarydf_smax_sim_redux<-summarydf_smax_sim[summarydf_smax_sim$scenario%in%c("highERLowError",  "ShiftERLowError", "lowERLowError",
    "decLinearProd_highERLowError",
  "incLinearProd_highERLowError", "decLinearProd_ShiftERLowError", "incLinearProd_ShiftERLowError")&
summarydf_smax_sim$model%in%c("autocorr", "rwa")&summarydf_smax_sim$method=="MLE",]


summarydf_smax_redux<-summarydf_smax[summarydf_smax$scenario%in%c( "highERLowError",  "ShiftERLowError", "lowERLowError",
    "decLinearProd_highERLowError","incLinearProd_highERLowError","decLinearProd_ShiftERLowError", "incLinearProd_ShiftERLowError" )&
  summarydf_smax$model%in%c("autocorr", "rwa") &summarydf_smax$method=="MLE",]




psmax_highERscn_line<-ggplot() + 
geom_pointrange(data=summarydf_smax_redux,aes(x=by,y= x.50.,ymin = x.2.5., ymax = x.97.5., color=model),alpha=.9)+
geom_line(data=summarydf_smax_sim_redux,aes(x=by,y= x),color="black", alpha=.8,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8,option = "E") +
scale_fill_viridis_d(begin=.1, end=.8,option = "E") +
coord_cartesian(ylim = c(50000,400000))+ 
mytheme+ 
ylab("Smax") +
xlab("year") +
#theme( strip.text.y = element_blank(),strip.text.x = element_blank())+
facet_grid(dynamics+ERtrend ~model, scales="free_y")
psmax_highERscn_line


head(df$scenario)
df_smax_est_redux<- df[df$parameter=="smax"&df$variable=="mode"&
df$scenario%in%c("highERLowError",  "ShiftERLowError", "lowERLowError","decLinearProd_highERLowError", 
    "incLinearProd_highERLowError","decLinearProd_ShiftERLowError", "incLinearProd_ShiftERLowError")&
df$model%in%c("autocorr", "rwa")&
df$method=="MLE",]



head(df_smsy_est_redux)

psmax_erscn_violin_abs<-ggplot(df_smax_est_redux) + 
geom_violin(aes(x=model,y=abs(bias), fill=model), scale="width", trim=TRUE, alpha=.7,adjust = 1.8)+
geom_boxplot(aes(x=model,y=abs(bias), fill=model),outlier.shape = NA, width=0.1)+
geom_hline(yintercept=0,color="black", alpha=.6,linewidth=1.2)+
 scale_fill_viridis_d(begin=.1, end=.8,,option = "E") +
 facet_grid(dynamics ~ERtrend, scales="free_y")+
 ylab("absolute bias in Smax") +
 coord_cartesian(ylim = c(0,150000))+ 
mytheme
psmax_erscn_violin_abs

multi.page.abs.smax_erscn <- ggarrange(psmax_highERscn_line, psmax_erscn_violin_abs,
                        nrow = 1, ncol = 2,
                        common.legend = TRUE,
                        legend="bottom")
multi.page.abs.smax_erscn
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/summary_figs/bias_trends_violins_smax_ERscn_wlowER.png",
    plot=multi.page.abs.smax_erscn)



#=============================
#alpha


summarydf_alpha_sim_redux<-summarydf_alpha_sim[summarydf_alpha_sim$scenario%in%c("highERLowError",  "ShiftERLowError", 
    "lowERLowError","decLinearProd_highERLowError",
  "incLinearProd_highERLowError", "decLinearProd_ShiftERLowError", "incLinearProd_ShiftERLowError")&
summarydf_alpha_sim$model%in%c("autocorr", "rwa")&summarydf_alpha_sim$method=="MLE",]


summarydf_alpha_redux<-summarydf_alpha[summarydf_alpha$scenario%in%c( "highERLowError",  "ShiftERLowError", 
    "lowERLowError", "decLinearProd_highERLowError","incLinearProd_highERLowError",
    "decLinearProd_ShiftERLowError", "incLinearProd_ShiftERLowError" )&
  summarydf_alpha$model%in%c("autocorr", "rwa") &summarydf_alpha$method=="MLE",]




palpha_highERscn_line<-ggplot() + 
geom_pointrange(data=summarydf_alpha_redux,aes(x=by,y= x.50.,ymin = x.2.5., ymax = x.97.5., color=model),alpha=.9)+
geom_line(data=summarydf_alpha_sim_redux,aes(x=by,y= x),color="black", alpha=.8,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8,option = "E") +
scale_fill_viridis_d(begin=.1, end=.8,option = "E") +
coord_cartesian(ylim = c(0,2.3))+ 
mytheme+ 
ylab("log(alpha)") +
xlab("year") +
#theme( strip.text.y = element_blank(),strip.text.x = element_blank())+
facet_grid(dynamics+ERtrend ~model, scales="free_y")
palpha_highERscn_line


df_alpha_est_redux<- df[df$parameter=="alpha"&df$variable=="mode"&
df$scenario%in%c("highERLowError",  "ShiftERLowError", "lowERLowError","decLinearProd_highERLowError", 
    "incLinearProd_highERLowError","decLinearProd_ShiftERLowError", "incLinearProd_ShiftERLowError")&
df$model%in%c("autocorr", "rwa")&
df$method=="MLE",]



head(df_smsy_est_redux)

palpha_erscn_violin_abs<-ggplot(df_alpha_est_redux) + 
geom_violin(aes(x=model,y=abs(bias), fill=model), scale="width", trim=TRUE, alpha=.7,adjust = 1.8)+
geom_boxplot(aes(x=model,y=abs(bias), fill=model),outlier.shape = NA, width=0.1)+
geom_hline(yintercept=0,color="black", alpha=.6,linewidth=1.2)+
 scale_fill_viridis_d(begin=.1, end=.8,,option = "E") +
 facet_grid(dynamics ~ERtrend, scales="free_y")+
 ylab("absolute bias in log(alpha)") +
 coord_cartesian(ylim = c(0,2))+ 
mytheme
psmax_erscn_violin_abs

multi.page.abs.smax_erscn <- ggarrange(psmax_highERscn_line, palpha_erscn_violin_abs,
                        nrow = 1, ncol = 2,
                        common.legend = TRUE,
                        legend="bottom")
multi.page.abs.smax_erscn
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/summary_figs/bias_trends_violins_smax_ERscn_wlowER.png",
    plot=multi.page.abs.smax_erscn)

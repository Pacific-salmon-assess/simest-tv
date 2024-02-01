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


source("code/read_base_data.R")


#========================================================================================================
#base case
#read in data
simPar <- read.csv("data/generic/SimPars.csv")

## Store relevant object names to help run simulation 
scenNames <- unique(simPar$scenario)

res1<-readRDS(file = "outs/simest/generic/resbase1.rds")
res2<-readRDS(file = "outs/simest/generic/resbase2.rds")


restmb<-rbind(res1,res2)



df<-reshape2::melt(resparam, id.vars=c("parameter","iteration","scenario","method","model","by", 
                                      "convergence","conv_warning","pbias","bias"))


df$scenario<-factor(df$scenario,levels=c("stationary",
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

df$scencode<-dplyr::case_match(df$scenario, 
      "stationary"~"Base1",
      "autocorr"~"Base2",
      "sigmaShift"~"Base3", 
      "decLinearProd"~"Base4",
      "sineProd"~"Base5",
      "regimeProd"~"Base6",
      "shiftProd"~"Base7",
      "decLinearCap"~"Base8",
      "regimeCap"~"Base9",
      "shiftCap"~"Base10", 
      "regimeProdCap"~"Base11",
      "decLinearProdshiftCap"~"Base12"
      )   




df$scencode <-factor(df$scencode, levels=c("Base1","Base2","Base3",
             "Base4","Base5","Base6",
              "Base7","Base8","Base9",
               "Base10","Base11","Base12"))



df$scentype<-dplyr::case_match(df$scenario, 
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


df$model<-factor(df$model,levels=c("simple",
                                   "autocorr", 
                                   "rwa",
                                   "hmma",
                                   "rwb",
                                   "hmmb",
                                   "rwab",
                                   "hmmab"  ))


head(df)

summarydf  <- df %>%
   group_by(scenario,parameter,
    method,model,by,variable,scencode,scentype) %>%
   reframe(qs = quantile(value, c(0.025, .5, 0.975),na.rm=T), prob = c("lower","median", "upper"))



summarydf <- reshape2::dcast(data=summarydf,  
    scenario + parameter + method + model + by + variable + scencode + scentype ~ prob, 
    value.var= "qs",fun.aggregate=mean)


head(summarydf)
summarydf_alpha<-summarydf[summarydf$parameter=="alpha"&
                            summarydf$variable=="mode",
                            ]

summarydf_alpha_sim<-summarydf[summarydf$parameter=="alpha"&
                                summarydf$variable=="sim",
                                ]



alphabase<-ggplot() + 
geom_pointrange(data=summarydf_alpha,aes(x=by,y= median,ymin = lower, ymax = upper, col=method),alpha=.6)+
geom_line(data=summarydf_alpha_sim,aes(x=by,y= median),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
coord_cartesian(ylim = c(0.2,2.7))+ 
mytheme+ 
ylab("alpha") +
xlab("year") +
facet_grid(scencode+scentype~model, scales="free_y")
alphabase
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/base/compareMCMC_MLE_alpha.png",
    plot=alpha1base)


summarydf_alpha_sim1<-summarydf_alpha_sim[summarydf_alpha_sim$scenario%in%c( "decLinearProd",        
"regimeProd",   "shiftProd",  "sineProd",   "stationary" ,"autocorr"),]#&summarydf_alpha_sim$model%in%c( "hmma","hmmb" "hmmab", "rwa",  "rwb", "rwab", "simple","autocorr" ),]

summarydf_alpha1<-summarydf_alpha[summarydf_alpha$scenario%in%c(       "decLinearProd",        
"regimeProd",   "shiftProd",  "sineProd",   "stationary","autocorr" ),]#&summarydf_alpha$model%in%c( "hmma", "hmmab", "rwa", "rwab", "simple" ),]

alpha1base<-ggplot() + 
geom_pointrange(data=summarydf_alpha1,aes(x=by,y=median,ymin =lower, ymax = upper, col=method),alpha=.6)+
geom_line(data=summarydf_alpha_sim1,aes(x=by,y=median),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
coord_cartesian(ylim = c(0.2,2.7))+ 
mytheme+ 
ylab(expression(log(alpha))) +
xlab("year") +
facet_grid(scencode+scentype~model, scales="free_y")
alpha1base
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/base/compareMCMC_MLE_alpha1.png",
    plot=alpha1base)
#MLE estimates are less biased and higher than MCMC


summarydf_alpha_sim2<-summarydf_alpha_sim[summarydf_alpha_sim$scenario%in%c(  "decLinearCap",  "sigmaShift",          
"decLinearProdshiftCap", "regimeCap",  "regimeProdCap", "shiftCap" ),]

summarydf_alpha2<-summarydf_alpha[summarydf_alpha$scenario%in%c(  "decLinearCap",  "sigmaShift",      
"decLinearProdshiftCap", "regimeCap",  "regimeProdCap", "shiftCap" ),]



alpha2base<-ggplot() + 
geom_pointrange(data=summarydf_alpha2,aes(x=by,y=median,ymin =lower, ymax =upper, col=method),alpha=.6)+
geom_line(data=summarydf_alpha_sim2,aes(x=by,y= median),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
coord_cartesian(ylim = c(0.2,2.7))+ 
mytheme + 
ylab(expression(log(alpha))) +
xlab("year") +
facet_grid(scencode+scentype~model, scales="free_y")
alpha2base
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/base/compareMCMC_MLE_alpha2.png",
    plot=alpha2base)
#MLE estimates are less biased and higher than MCMC




#=======================================================
#b estimates

summarydf_smax<-summarydf[summarydf$parameter=="smax"&
                            summarydf$variable=="mode",
                            ]

summarydf_smax_sim<-summarydf[summarydf$parameter=="smax"&
                                summarydf$variable=="sim",
                                ]


summarydf_smax_sim1<-summarydf_smax_sim[summarydf_smax_sim$scenario%in%c( "decLinearProd",        
"regimeProd",   "shiftProd",  "sineProd",   "stationary","autocorr"  ),]

summarydf_smax1<-summarydf_smax[summarydf_smax$scenario%in%c("decLinearProd",        
"regimeProd",   "shiftProd",  "sineProd",   "stationary","autocorr"  ),]



smax1base<-ggplot() + 
geom_pointrange(data=summarydf_smax1,aes(x=by,y=median,ymin =lower, ymax = upper, col=method),alpha=.6)+
geom_line(data=summarydf_smax_sim1,aes(x=by,y= median),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
mytheme + 
ylab(expression(S[max])) +
xlab("year") +
coord_cartesian(ylim = c(60000,400000))+ 
facet_grid(scencode+scentype~model, scales="free_y")
smax1base
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/base/compareMCMC_MLE_smax1.png",
    plot=smax1base)



summarydf_smax_sim2<-summarydf_smax_sim[summarydf_smax_sim$scenario%in%c(  "decLinearCap",   "sigmaShift",      
"decLinearProdshiftCap", "regimeCap",  "regimeProdCap", "shiftCap" ),]

summarydf_smax2<-summarydf_smax[summarydf_smax$scenario%in%c( "decLinearCap",    "sigmaShift",     
"decLinearProdshiftCap", "regimeCap",  "regimeProdCap", "shiftCap"),]

smax2base<-ggplot() + 
geom_pointrange(data=summarydf_smax2,aes(x=by,y=median,ymin =lower, ymax = upper, col=method),alpha=.6)+
geom_line(data=summarydf_smax_sim2,aes(x=by,y= median),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
mytheme + 
ylab(expression(S[max])) +
xlab("year") +
coord_cartesian(ylim = c(60000,400000))+ 
facet_grid(scencode+scentype~model, scales="free_y")
smax2base
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/base/compareMCMC_MLE_smax2.png",
    plot=smax2base)





#=======================================================
#smsy estimates


summarydf_smsy<-summarydf[summarydf$parameter=="smsy"&
                            summarydf$variable=="mode",
                            ]

summarydf_smsy_sim<-summarydf[summarydf$parameter=="smsy"&
                                summarydf$variable=="sim",
                                ]



summarydf_smsy_sim1<-summarydf_smsy_sim[summarydf_smsy_sim$scenario%in%c( "decLinearProd",        
"regimeProd",   "shiftProd",  "sineProd",   "stationary","autocorr"  ),]


summarydf_smsy1<-summarydf_smsy[summarydf_smsy$scenario%in%c("decLinearProd",        
"regimeProd",   "shiftProd",  "sineProd",   "stationary","autocorr"  ),]

summarydf_smsy1[summarydf_smsy1$model=="hmmab"&summarydf_smsy1$scenario=="sineProd",]


smsy1base<-ggplot() + 
geom_pointrange(data=summarydf_smsy1,aes(x=by,y=median,ymin =lower, ymax = upper, col=method),alpha=.6)+
geom_line(data=summarydf_smsy_sim1,aes(x=by,y= median),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
mytheme + 
ylab(expression(S[MSY])) +
xlab("year") +
coord_cartesian(ylim = c(20000,150000))+ 
facet_grid(scencode+scentype~model, scales="free_y")
smsy1base
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/base/compareMCMC_MLE_smsy1.png",
    plot=smsy1base)



summarydf_smsy_sim2<-summarydf_smsy_sim[summarydf_smsy_sim$scenario%in%c(  "decLinearCap",   "sigmaShift",      
"decLinearProdshiftCap", "regimeCap",  "regimeProdCap", "shiftCap"  ),]#&summarydf_smax_sim$model%in%c( "hmma","hmmb", "hmmab", "rwa", "rwb", "rwab", "simple" ),]

summarydf_smsy2<-summarydf_smsy[summarydf_smsy$scenario%in%c( "decLinearCap",   "sigmaShift",      
"decLinearProdshiftCap", "regimeCap",  "regimeProdCap", "shiftCap" ),]#&summarydf_smax$model%in%c( "hmma","hmmb", "hmmab", "rwa", "rwb", "rwab", "simple" ),]



smsy2base<-ggplot() + 
geom_pointrange(data=summarydf_smsy2,aes(x=by,y=median,ymin = lower, ymax = upper, col=method),alpha=.6)+
geom_line(data=summarydf_smsy_sim2,aes(x=by,y=median),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
mytheme + 
ylab(expression(S[MSY])) +
xlab("year") +
coord_cartesian(ylim = c(20000,150000))+ 
facet_grid(scencode+scentype~model, scales="free_y")
smsy2base
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/base/compareMCMC_MLE_smsy2.png",
    plot=smsy2base)



#=======================================================
#smsy estimates


summarydf_umsy<-summarydf[summarydf$parameter=="umsy"&
                            summarydf$variable=="mode",
                            ]

summarydf_umsy_sim<-summarydf[summarydf$parameter=="umsy"&
                                summarydf$variable=="sim",
                                ]



summarydf_umsy_sim1<-summarydf_umsy_sim[summarydf_umsy_sim$scenario%in%c( "decLinearProd",        
"regimeProd",   "shiftProd",  "sineProd",   "stationary","autocorr"  ),]


summarydf_umsy1<-summarydf_umsy[summarydf_umsy$scenario%in%c("decLinearProd",        
"regimeProd",   "shiftProd",  "sineProd",   "stationary","autocorr"  ),]




umsy1base<-ggplot() + 
geom_pointrange(data=summarydf_umsy1,aes(x=by,y=median,ymin =lower, ymax = upper, col=method),alpha=.6)+
geom_line(data=summarydf_umsy_sim1,aes(x=by,y= median),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
mytheme + 
ylab(expression(U[MSY])) +
xlab("year") +
coord_cartesian(ylim = c(0,1))+ 
facet_grid(scencode+scentype~model, scales="free_y")
umsy1base
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/base/compareMCMC_MLE_umsy1.png",
    plot=umsy1base)



summarydf_umsy_sim2<-summarydf_umsy_sim[summarydf_umsy_sim$scenario%in%c(  "decLinearCap",   "sigmaShift",      
"decLinearProdshiftCap", "regimeCap",  "regimeProdCap", "shiftCap"  ),]

summarydf_umsy2<-summarydf_umsy[summarydf_umsy$scenario%in%c( "decLinearCap",   "sigmaShift",      
"decLinearProdshiftCap", "regimeCap",  "regimeProdCap", "shiftCap" ),]



umsy2base<-ggplot() + 
geom_pointrange(data=summarydf_umsy2,aes(x=by,y=median,ymin = lower, ymax = upper, col=method),alpha=.6)+
geom_line(data=summarydf_umsy_sim2,aes(x=by,y=median),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
mytheme + 
ylab(expression(S[MSY])) +
xlab("year") +
coord_cartesian(ylim = c(20000,150000))+ 
facet_grid(scencode+scentype~model, scales="free_y")
smsy2base
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/base/compareMCMC_MLE_umsy2.png",
    plot=smsy2base)




#-----------------------------------------------------------------------------------
#new plots 

unique(summarydf_alpha_sim$scenario)
head(summarydf_alpha)
# Bias in alpha
summarydf_alpha_sim_redux<-summarydf_alpha_sim[summarydf_alpha_sim$scenario%in%c("decLinearProd", "regimeProd", "sineProd")&
summarydf_alpha$model%in%c("autocorr", "rwa", "hmma")&summarydf_alpha$method=="MLE",]

summarydf_alpha_sim_redux$scenario2<-case_match(summarydf_alpha_sim_redux$scenario,
    "decLinearProd"~ "linear decline",
     "regimeProd" ~ "shift increase" , 
     "sineProd" ~ "sine fluctuation")

summarydf_alpha_redux<-summarydf_alpha[summarydf_alpha$scenario%in%c( "decLinearProd", "regimeProd", "sineProd" )&
  summarydf_alpha$model%in%c("autocorr", "rwa", "hmma") &summarydf_alpha$method=="MLE",]


summarydf_alpha_redux$scenario2<-case_match(summarydf_alpha_redux$scenario,
    "decLinearProd"~ "linear decline",
     "regimeProd" ~ "shift increase" , 
     "sineProd" ~ "sine fluctuation")

head(summarydf_alpha_redux)

palpha_line<-ggplot() + 
geom_pointrange(data=summarydf_alpha_redux,aes(x=by,y= median,ymin =lower, ymax = upper, color=model),alpha=.9)+
geom_line(data=summarydf_alpha_sim_redux,aes(x=by,y=median),color="black", alpha=.8,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8,option = "E") +
scale_fill_viridis_d(begin=.1, end=.8,option = "E") +
#coord_cartesian(ylim = c(0.2,3.0))+ 
mytheme+ 
ylab("log(alpha)") +
xlab("year") +
theme( strip.text.y = element_blank(),strip.text.x = element_blank())+
facet_grid(scenario2~model, scales="free_y")
palpha_line


df_alpha_est_redux<- df[df$parameter=="alpha"&df$variable=="mode"&
df$scenario%in%c("decLinearProd", "shiftProd", "sineProd")&
df$model%in%c("autocorr", "rwa", "hmma")&
df$method=="MLE",]


df_alpha_est_redux$scenario2<-case_match(df_alpha_est_redux$scenario,
    "decLinearProd"~ "linear decline",
     "shiftProd" ~ "shift decline" , 
     "sineProd" ~ "sine fluctuation")

head(df)
palpha_violin<-ggplot(df_alpha_est_redux) + 
geom_violin(aes(x=model,y=bias, fill=model), scale="width", trim=TRUE, alpha=.7,adjust = 1.8)+
geom_boxplot(aes(x=model,y=bias, fill=model),outlier.shape = NA, width=0.1)+
geom_hline(yintercept=0,color="black", alpha=.6,linewidth=1.2)+
 scale_fill_viridis_d(begin=.1, end=.8,,option = "E") +
 facet_grid(scenario2~., scales="free_y")+
 mytheme
palpha_violin

multi.page <- ggarrange(palpha_line, palpha_violin,
                        nrow = 1, ncol = 2,
                        common.legend = TRUE,
                        legend="bottom")
multi.page




palpha_violin_abs<-ggplot(df_alpha_est_redux) + 
geom_violin(aes(x=model,y=abs(bias), fill=model), scale="width", trim=TRUE, alpha=.7,adjust = 1.8)+
geom_boxplot(aes(x=model,y=abs(bias), fill=model),outlier.shape = NA, width=0.1)+
geom_hline(yintercept=0,color="black", alpha=.6,linewidth=1.2)+
 scale_fill_viridis_d(begin=.1, end=.8,,option = "E") +
 facet_grid(scenario2~., scales="free_y")+
 ylab("absolute bias in log(alpha)") +
mytheme
palpha_violin_abs

multi.page.abs <- ggarrange(palpha_line, palpha_violin_abs,
                        nrow = 1, ncol = 2,
                        common.legend = TRUE,
                        legend="bottom")
multi.page.abs

ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/summary_figs/bias_trends_violins_log_alpha.png",
    plot=multi.page.abs)


#------
#GUIDELINES

#cv for line plot


#year based cv

cvdf<-aggregate(resparam$mode, list(parameter=resparam$parameter,
                               scenario=resparam$scenario,
                               method=resparam$method,
                               model=resparam$model,
                               by=resparam$by ), function(x){sd(x,na.rm=T)/abs(mean(x,na.rm=T))})




meancvdf<-aggregate(cvdf$x, list(parameter=cvdf$parameter,
                               scenario=cvdf$scenario,
                               method=cvdf$method,
                               model=cvdf$model
                              ), function(x){mean(x,na.rm=T)})

meancvdf_redux<-meancvdf[meancvdf$parameter=="alpha"&
meancvdf$scenario%in%c("decLinearProd", "regimeProd", "sineProd")&
meancvdf$model%in%c("autocorr", "rwa", "hmma")&
meancvdf$method=="MLE",]

meancvdf_redux$meancv<-round(meancvdf_redux$x,2)
meancvdf_redux$scenario2<-case_match(meancvdf_redux$scenario,
    "decLinearProd"~ "linear decline",
     "regimeProd" ~ "shift increase" , 
     "sineProd" ~ "sine fluctuation")
meancvdf_redux$model<-factor(meancvdf_redux$model, levels=c("autocorr", "rwa", "hmma"))

head(summarydf_alpha_redux)

palpha_line_cv<-ggplot() + 
geom_pointrange(data=summarydf_alpha_redux,aes(x=by,y= median,ymin = lower, ymax = upper, color=model),alpha=.9)+
geom_line(data=summarydf_alpha_sim_redux,aes(x=by,y= median),color="black", alpha=.8,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8,option = "E") +
scale_fill_viridis_d(begin=.1, end=.8,option = "E") +
#coord_cartesian(ylim = c(0.2,3.0))+ 
mytheme+ 
ylab("log(alpha)") +
xlab("year") +
theme( strip.text.y = element_blank(),strip.text.x = element_blank())+
geom_text(data = meancvdf_redux, aes(x=-Inf,y=Inf,hjust=0,
                vjust=1.0,label=meancv), size=6)+
facet_grid(scenario2~model, scales="free_y")
palpha_line_cv


#pbias plots

df_alpha_est_redux<- df[df$parameter=="alpha"&df$variable=="mode"&
df$scenario%in%c("decLinearProd", "shiftProd", "sineProd")&
df$model%in%c("autocorr", "rwa", "hmma")&
df$method=="MLE",]



palpha_violin_abspbias<-ggplot(df_alpha_est_redux) + 
geom_violin(aes(x=model,y=(pbias), fill=model), scale="width", trim=TRUE, alpha=.7,adjust = 1.8)+
geom_boxplot(aes(x=model,y=(pbias), fill=model),outlier.shape = NA, width=0.1)+
geom_hline(yintercept=0,color="black", alpha=.6,linewidth=1.2)+
 scale_fill_viridis_d(begin=.1, end=.8,,option = "E") +
 facet_grid(scenario2~., scales="free_y")+
 ylab("absolute bias in log(alpha)") +
mytheme
palpha_violin_abspbias


multi.page.abs <- ggarrange(palpha_line_cv, palpha_violin_abs,
                        nrow = 1, ncol = 2,
                        common.legend = TRUE,
                        legend="bottom")
multi.page.abs
#==================================================================
# Bias in beta

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


unique(summarydf_smax_sim$scenario)
summarydf_smax_sim_redux<-summarydf_smax_sim[summarydf_smax_sim$scenario%in%c("decLinearCap", "regimeCap", "shiftCap")&
summarydf_smax$model%in%c("autocorr", "rwb", "hmmb")&summarydf_alpha$method=="MLE",]

summarydf_smax_sim_redux$scenario2<-case_match(summarydf_smax_sim_redux$scenario,
    "decLinearCap"~ "linear decline",
     "regimeCap" ~ "regime increase" , 
     "shiftCap" ~ "shift decline")

summarydf_smax_redux<-summarydf_smax[summarydf_smax$scenario%in%c( "decLinearCap", "regimeCap", "shiftCap" )&
  summarydf_smax$model%in%c("autocorr", "rwb", "hmmb") &summarydf_smax$method=="MLE",]


summarydf_smax_redux$scenario2<-case_match(summarydf_smax_redux$scenario,
    "decLinearCap"~ "linear decline",
     "regimeCap" ~ "regime increase" , 
     "shiftCap" ~ "shift decline")



psmax_line<-ggplot() + 
geom_pointrange(data=summarydf_smax_redux,aes(x=by,y= x.50.,ymin = x.2.5., ymax = x.97.5., color=model),alpha=.9)+
geom_line(data=summarydf_smax_sim_redux,aes(x=by,y= x),color="black", alpha=.8,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8,option = "E") +
scale_fill_viridis_d(begin=.1, end=.8,option = "E") +
coord_cartesian(ylim = c(50000,400000))+ 
mytheme+ 
ylab("Smax") +
xlab("year") +
theme( strip.text.y = element_blank(),strip.text.x = element_blank())+
facet_grid(scenario2~model, scales="free_y")
psmax_line


df_smax_est_redux<- df[df$parameter=="smax"&df$variable=="mode"&
df$scenario%in%c("decLinearCap", "regimeCap", "shiftCap")&
df$model%in%c("autocorr", "rwb", "hmmb")&
df$method=="MLE",]
head(df_smax_est_redux)

df_smax_est_redux$scenario2<-case_match(df_smax_est_redux$scenario,
    "decLinearCap"~ "linear decline",
     "regimeCap" ~ "regime increase" , 
     "shiftCap" ~ "shift decline")

head(df)
psmax_violin<-ggplot(df_smax_est_redux) + 
geom_violin(aes(x=model,y=bias, fill=model), scale="width", trim=TRUE, alpha=.7,adjust = 1.8)+
geom_boxplot(aes(x=model,y=bias, fill=model),outlier.shape = NA, width=0.1)+
geom_hline(yintercept=0,color="black", alpha=.6,linewidth=1.2)+
 scale_fill_viridis_d(begin=.1, end=.8,,option = "E") +
 facet_grid(scenario2~., scales="free_y")+
 ggtitle("estimation model")+
 coord_cartesian(ylim = c(-350000,350000))+ 
mytheme
psmax_violin


multi.page <- ggarrange(psmax_line, psmax_violin,
                        nrow = 1, ncol = 2,
                        common.legend = TRUE,
                        legend="bottom")
multi.page




psmax_violin_abs<-ggplot(df_smax_est_redux) + 
geom_violin(aes(x=model,y=abs(bias), fill=model), scale="width", trim=TRUE, alpha=.7,adjust = 1.8)+
geom_boxplot(aes(x=model,y=abs(bias), fill=model),outlier.shape = NA, width=0.1)+
geom_hline(yintercept=0,color="black", alpha=.6,linewidth=1.2)+
 scale_fill_viridis_d(begin=.1, end=.8,,option = "E") +
 facet_grid(scenario2~., scales="free_y")+
 ylab("absolute bias in smax") +
 coord_cartesian(ylim = c(0,450000))+ 
mytheme
psmax_violin_abs

multi.page.abs.smax <- ggarrange(psmax_line, psmax_violin_abs,
                        nrow = 1, ncol = 2,
                        common.legend = TRUE,
                        legend="bottom")
multi.page.abs.smax

ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/summary_figs/bias_trends_violins_smax_filterbounds.png",
    plot=multi.page.abs.smax)



head(df_smax_sim)
head(df_smax_est_redux)

psmax_trends_abs<-ggplot(df_smax_est_redux) + 
geom_line(aes(x=by,y=value, color=model, group=iteration), alpha=.2)+
geom_line(data=summarydf_smax_sim_redux,aes(x=by,y= x),color="black", alpha=.8,linewidth=1.2)+
#geom_boxplot(aes(x=model,y=abs(bias), fill=model),outlier.shape = NA, width=0.1)+
 scale_color_viridis_d(begin=.1, end=.8) +
 facet_grid(scenario2~model, scales="free_y")+
 ylab("absolute bias in smax") +
 ggtitle("estimation model")+
 coord_cartesian(ylim = c(0,750000))+ 
mytheme
psmax_trends_abs


#==================================================================
# Bias in alpha and beta




df_alpha.smax_sim<- df[df$parameter%in%c("alpha","smax")&df$variable=="sim",]
df_alpha.smax_est<- df[df$parameter%in%c("alpha","smax")&df$variable=="mode",]


summarydf_alpha.smax<-aggregate(df_alpha.smax_est$value,by=list(scenario=df_alpha.smax_est$scenario, 
    method=df_alpha.smax_est$method, 
    model=df_alpha.smax_est$model,
    by=df_alpha.smax_est$by, 
    parameter=df_alpha.smax_est$parameter),
    function(x) {quantile(x,probs = c(0.025, .5, 0.975))})
summarydf_alpha.smax<-do.call(data.frame, summarydf_alpha.smax)

summarydf_alpha.smax_sim<-aggregate(df_alpha.smax_sim$value,by=list(scenario=df_alpha.smax_sim$scenario, 
    method=df_alpha.smax_sim$method, 
    model=df_alpha.smax_sim$model,
    by=df_alpha.smax_sim$by,
     parameter=df_alpha.smax_sim$parameter ),
    function(x) {unique(x)})


unique(summarydf_alpha.smax_sim$scenario)
summarydf_alpha.smax_sim_redux<-summarydf_alpha.smax_sim[summarydf_alpha.smax_sim$scenario%in%c("regimeProdCap",
         "decLinearProdshiftCap")&
summarydf_alpha.smax_sim$model%in%c("autocorr", "rwab", "hmmab")&summarydf_alpha.smax_sim$method=="MLE",]

summarydf_alpha.smax_sim_redux$scenario2<-case_match(summarydf_alpha.smax_sim_redux$scenario,
    "regimeProdCap"~ "regime",
     "decLinearProdshiftCap" ~ "mixed trends" )

summarydf_alpha.smax_redux<-summarydf_alpha.smax[summarydf_alpha.smax$scenario%in%c( "regimeProdCap",
         "decLinearProdshiftCap" )&
  summarydf_alpha.smax$model%in%c("autocorr", "rwab", "hmmab") &summarydf_alpha.smax$method=="MLE",]


summarydf_alpha.smax_redux$scenario2<-case_match(summarydf_alpha.smax_redux$scenario,
      "regimeProdCap"~ "regime",
     "decLinearProdshiftCap" ~ "mixed trends")

head(summarydf_alpha.smax_redux)
head(summarydf_alpha.smax_sim_redux)

palpha.smax_line<-ggplot() + 
geom_pointrange(data=summarydf_alpha.smax_redux,aes(x=by,y= x.50.,ymin = x.2.5., ymax = x.97.5., color=model),alpha=.9)+
geom_line(data=summarydf_alpha.smax_sim_redux,aes(x=by,y= x),color="black", alpha=.8,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8,option = "E") +
scale_fill_viridis_d(begin=.1, end=.8,option = "E") +
#coord_cartesian(ylim = c(50000,400000))+ 
mytheme+ 
ylab("parameter value") +
xlab("year") +
facet_grid(scenario2+parameter~model, scales="free_y")+
theme( strip.text.y = element_blank(),strip.text.x = element_blank())
palpha.smax_line


df_alpha.smax_est_redux<- df[df$parameter%in%c("alpha","smax")&df$variable=="mode"&
df$scenario%in%c("regimeProdCap",
         "decLinearProdshiftCap")&
df$model%in%c("autocorr", "rwab", "hmmab")&
df$method=="MLE",]
head(df_alpha.smax_est_redux)

df_alpha.smax_est_redux$scenario2<-case_match(df_alpha.smax_est_redux$scenario,
    "regimeProdCap"~ "regime",
     "decLinearProdshiftCap" ~ "mixed trends")

head(df_alpha.smax_est_redux)
head(df)
palpha.smax_violin<-ggplot(df_alpha.smax_est_redux) + 
geom_violin(aes(x=model,y=bias, fill=model), scale="width", trim=TRUE, alpha=.7,adjust = 1.8)+
geom_boxplot(aes(x=model,y=bias, fill=model),outlier.shape = NA, width=0.1)+
geom_hline(yintercept=0,color="black", alpha=.6,linewidth=1.2)+
 scale_fill_viridis_d(begin=.1, end=.8,,option = "E") +
 facet_grid(scenario2+parameter~., scales="free_y")+
 #ggtitle("estimation model")+
 #coord_cartesian(ylim = c(-350000,350000))+ 
mytheme
palpha.smax_violin


multi.page.alpha.smax <- ggarrange(palpha.smax_line, palpha.smax_violin,
                        nrow = 1, ncol = 2,
                        common.legend = TRUE,
                        legend="bottom")
multi.page.alpha.smax

summary(log10(abs(df_smax_est_redux$bias)))




palpha.smax_violin_abs<-ggplot(df_alpha.smax_est_redux) + 
geom_violin(aes(x=model,y=abs(bias), fill=model), scale="width", trim=TRUE, alpha=.7,adjust = 1.8)+
geom_boxplot(aes(x=model,y=abs(bias), fill=model),outlier.shape = NA, width=0.1)+
geom_hline(yintercept=0,color="black", alpha=.6,linewidth=1.2)+
 scale_fill_viridis_d(begin=.1, end=.8,,option = "E") +
 facet_grid(scenario2+parameter~., scales="free_y")+
 ylab("absolute bias in parameter") +
 #ggtitle("estimation model")+
 #coord_cartesian(ylim = c(0,350000))+ 
 
mytheme
#theme( strip.text.y = element_blank())
palpha.smax_violin_abs

multi.page.abs.alpha.smax <- ggarrange(palpha.smax_line, palpha.smax_violin_abs,
                        nrow = 1, ncol = 2,
                        common.legend = TRUE,
                        legend="bottom")
multi.page.abs.alpha.smax

ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/summary_figs/bias_trends_violins_alphasmax.png",
    plot=multi.page.abs.alpha.smax)




#=========================================================================================================\
#Same as baove but smsy-focused

# alpha varies


df_smsy_sim<- df[df$parameter%in%c("smsy")&df$variable=="sim",]
df_smsy_est<- df[df$parameter%in%c("smsy")&df$variable=="mode",]


summarydf_smsy<-aggregate(df_smsy_est$value,by=list(scenario=df_smsy_est$scenario, 
    method=df_smsy_est$method, 
    model=df_smsy_est$model,
    by=df_smsy_est$by ),
    function(x) {quantile(x,probs = c(0.025, .5, 0.975))})
summarydf_smsy<-do.call(data.frame, summarydf_smsy)

summarydf_smsy_sim<-aggregate(df_smsy_sim$value,by=list(scenario=df_smsy_sim$scenario, 
    method=df_smsy_sim$method, 
    model=df_smsy_sim$model,
    by=df_smsy_sim$by ),
    function(x) {unique(x)})

head(summarydf_smsy_sim_redux_alpha)

summarydf_smsy_sim_redux_alpha<-summarydf_smsy_sim[summarydf_smsy_sim$scenario%in%c("decLinearProd", "shiftProd", "sineProd")&
summarydf_smsy_sim$model%in%c("autocorr", "rwa", "hmma")&summarydf_smsy_sim$method=="MLE",]

summarydf_smsy_sim_redux_alpha$scenario2<-case_match(summarydf_smsy_sim_redux_alpha$scenario,
    "decLinearProd"~ "linear decline",
     "shiftProd" ~ "shift decline" , 
     "sineProd" ~ "sine fluctuation")

summarydf_smsy_redux_alpha<-summarydf_smsy[summarydf_smsy$scenario%in%c( "decLinearProd", "shiftProd", "sineProd" )&
  summarydf_smsy$model%in%c("autocorr", "rwa", "hmma") &summarydf_smsy$method=="MLE",]


summarydf_smsy_redux_alpha$scenario2<-case_match(summarydf_smsy_redux_alpha$scenario,
    "decLinearProd"~ "linear decline",
     "shiftProd" ~ "shift decline" , 
     "sineProd" ~ "sine fluctuation")
head(summarydf_smsy)
head(summarydf_smsy_redux_alpha)

psmsy_alphascn_line<-ggplot() + 
geom_pointrange(data=summarydf_smsy_redux_alpha,aes(x=by,y= x.50.,ymin = x.2.5., ymax = x.97.5., color=model),alpha=.9)+
geom_line(data=summarydf_smsy_sim_redux_alpha,aes(x=by,y= x),color="black", alpha=.8,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8,option = "E") +
scale_fill_viridis_d(begin=.1, end=.8,option = "E") +
coord_cartesian(ylim = c(20000,170000))+ 
mytheme+ 
ylab("Smsy") +
xlab("year") +
theme( strip.text.y = element_blank(),strip.text.x = element_blank())+
facet_grid(scenario2~model, scales="free_y")
psmsy_alphascn_line


df_smsy_est_redux_alpha<- df[df$parameter=="smsy"&df$variable=="mode"&
df$scenario%in%c("decLinearProd", "shiftProd", "sineProd")&
df$model%in%c("autocorr", "rwa", "hmma")&
df$method=="MLE",]


df_smsy_est_redux_alpha$scenario2<-case_match(df_smsy_est_redux_alpha$scenario,
    "decLinearProd"~ "linear decline",
     "shiftProd" ~ "shift decline" , 
 



psmsy_alphascn_violin_abs<-ggplot(df_smsy_est_redux_alpha) + 
geom_violin(aes(x=model,y=abs(bias), fill=model), scale="width", trim=TRUE, alpha=.7,adjust = 1.8)+
geom_boxplot(aes(x=model,y=abs(bias), fill=model),outlier.shape = NA, width=0.1)+
geom_hline(yintercept=0,color="black", alpha=.6,linewidth=1.2)+
 scale_fill_viridis_d(begin=.1, end=.8,,option = "E") +
 facet_grid(scenario2~., scales="free_y")+
 ylab("absolute bias in log(alpha)") +
 coord_cartesian(ylim = c(0,250000))+ 
mytheme
psmsy_alphascn_violin_abs

multi.page.abs.smsy_alphascn <- ggarrange(psmsy_alphascn_line, psmsy_alphascn_violin_abs,
                        nrow = 1, ncol = 2,
                        common.legend = TRUE,
                        legend="bottom")
multi.page.abs.smsy_alphascn

ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/summary_figs/bias_trends_violins_smsy_alphascn.png",
    plot=multi.page.abs.smsy_alphascn)

#==================================================================
# Bias in beta


unique(summarydf_smax_sim$scenario)
summarydf_smsy_sim_redux_smax<-summarydf_smsy_sim[summarydf_smsy_sim$scenario%in%c("decLinearCap", "regimeCap", "shiftCap")&
summarydf_smsy$model%in%c("autocorr", "rwb", "hmmb")&summarydf_smsy$method=="MLE",]

summarydf_smsy_sim_redux_smax$scenario2<-case_match(summarydf_smsy_sim_redux_smax$scenario,
    "decLinearCap"~ "linear decline",
     "regimeCap" ~ "regime increase" , 
     "shiftCap" ~ "shift decline")

summarydf_smsy_redux_smax<-summarydf_smsy[summarydf_smsy$scenario%in%c( "decLinearCap", "regimeCap", "shiftCap" )&
  summarydf_smsy$model%in%c("autocorr", "rwb", "hmmb") &summarydf_smsy$method=="MLE",]


summarydf_smsy_redux_smax$scenario2<-case_match(summarydf_smsy_redux_smax$scenario,
    "decLinearCap"~ "linear decline",
     "regimeCap" ~ "regime increase" , 
     "shiftCap" ~ "shift decline")



psmsy_smaxscn_line<-ggplot() + 
geom_pointrange(data=summarydf_smsy_redux_smax,aes(x=by,y= x.50.,ymin = x.2.5., ymax = x.97.5., color=model),alpha=.9)+
geom_line(data=summarydf_smsy_sim_redux_smax,aes(x=by,y= x),color="black", alpha=.8,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8,option = "E") +
scale_fill_viridis_d(begin=.1, end=.8,option = "E") +
coord_cartesian(ylim = c(50000,170000))+ 
mytheme+ 
ylab("Smax") +
xlab("year") +
theme( strip.text.y = element_blank(),strip.text.x = element_blank())+
facet_grid(scenario2~model, scales="free_y")
psmsy_smaxscn_line


df_smsy_est_redux_smax<- df[df$parameter=="smsy"&df$variable=="mode"&
df$scenario%in%c("decLinearCap", "regimeCap", "shiftCap")&
df$model%in%c("autocorr", "rwb", "hmmb")&
df$method=="MLE",]
head(df_smax_est_redux)

df_smsy_est_redux_smax$scenario2<-case_match(df_smsy_est_redux_smax$scenario,
    "decLinearCap"~ "linear decline",
     "regimeCap" ~ "regime increase" , 
     "shiftCap" ~ "shift decline")



psmsy_smaxscn_violin_abs<-ggplot(df_smsy_est_redux_smax) + 
geom_violin(aes(x=model,y=abs(bias), fill=model), scale="width", trim=TRUE, alpha=.7,adjust = 1.8)+
geom_boxplot(aes(x=model,y=abs(bias), fill=model),outlier.shape = NA, width=0.1)+
geom_hline(yintercept=0,color="black", alpha=.6,linewidth=1.2)+
 scale_fill_viridis_d(begin=.1, end=.8,,option = "E") +
 facet_grid(scenario2~., scales="free_y")+
 ylab("absolute bias in smax") +
 coord_cartesian(ylim = c(0,200000))+ 
mytheme
psmsy_smaxscn_violin_abs

multi.page.abs.smsy.smaxscn <- ggarrange(psmsy_smaxscn_line, psmsy_smaxscn_violin_abs,
                        nrow = 1, ncol = 2,
                        common.legend = TRUE,
                        legend="bottom")
multi.page.abs.smsy.smaxscn

ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/summary_figs/bias_trends_violins_smsy_smaxscn.png",
    plot=multi.page.abs.smsy.smaxscn)


#==================================================================
# Bias in alpha and beta

summarydf_smsy_sim_redux_bothscn<-summarydf_smsy_sim[summarydf_smsy_sim$scenario%in%c("regimeProdCap",
         "decLinearProdshiftCap")&
summarydf_smsy_sim$model%in%c("autocorr", "rwab", "hmmab")&summarydf_smsy_sim$method=="MLE",]

summarydf_smsy_sim_redux_bothscn$scenario2<-case_match(summarydf_smsy_sim_redux_bothscn$scenario,
    "regimeProdCap"~ "regime",
     "decLinearProdshiftCap" ~ "mixed trends" )

summarydf_smsy_redux_bothscn<-summarydf_smsy[summarydf_smsy$scenario%in%c( "regimeProdCap",
         "decLinearProdshiftCap" )&
  summarydf_smsy$model%in%c("autocorr", "rwab", "hmmab") &summarydf_smsy$method=="MLE",]


summarydf_smsy_redux_bothscn$scenario2<-case_match(summarydf_smsy_redux_bothscn$scenario,
      "regimeProdCap"~ "regime",
     "decLinearProdshiftCap" ~ "mixed trends")
head(summarydf_smsy_sim_redux_bothscn)
psmsy_bothscn_line<-ggplot() + 
geom_pointrange(data=summarydf_smsy_redux_bothscn,aes(x=by,y= x.50.,ymin = x.2.5., ymax = x.97.5., color=model),alpha=.9)+
geom_line(data=summarydf_smsy_sim_redux_bothscn,aes(x=by,y= x),color="black", alpha=.8,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8,option = "E") +
scale_fill_viridis_d(begin=.1, end=.8,option = "E") +
#coord_cartesian(ylim = c(50000,400000))+ 
mytheme+ 
ylab("parameter value") +
xlab("year") +
theme( strip.text.y = element_blank(),strip.text.x = element_blank())+
facet_grid(scenario2~model, scales="free_y")
psmsy_bothscn_line


df_smsy_est_redux_both<- df[df$parameter%in%c("smsy")&df$variable=="mode"&
df$scenario%in%c("regimeProdCap",
         "decLinearProdshiftCap")&
df$model%in%c("autocorr", "rwab", "hmmab")&
df$method=="MLE",]

df_smsy_est_redux_both$scenario2<-case_match(df_smsy_est_redux_both$scenario,
    "regimeProdCap"~ "regime",
     "decLinearProdshiftCap" ~ "mixed trends")



psmsy_bothscn_violin_abs<-ggplot(df_smsy_est_redux_both) + 
geom_violin(aes(x=model,y=abs(bias), fill=model), scale="width", trim=TRUE, alpha=.7,adjust = 1.8)+
geom_boxplot(aes(x=model,y=abs(bias), fill=model),outlier.shape = NA, width=0.1)+
geom_hline(yintercept=0,color="black", alpha=.6,linewidth=1.2)+
 scale_fill_viridis_d(begin=.1, end=.8,,option = "E") +
 facet_grid(scenario2~., scales="free_y")+
 ylab("absolute bias in parameter") +
 #ggtitle("estimation model")+
 #coord_cartesian(ylim = c(0,350000))+ 
 
mytheme
#theme( strip.text.y = element_blank())
psmsy_bothscn_violin_abs

multi.page.abs.smsy.bothscn <- ggarrange(psmsy_bothscn_line,psmsy_bothscn_violin_abs,
                        nrow = 1, ncol = 2,
                        common.legend = TRUE,
                        legend="bottom")
multi.page.abs.smsy.bothscn 

ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/summary_figs/bias_trends_violins_smsy_bothscn.png",
    plot=multi.page.abs.smsy.bothscn )






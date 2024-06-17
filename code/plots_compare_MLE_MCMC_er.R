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
        theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
              legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="bold", size=12),
              axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
)


source("code/read_er_data.R")

#========================================================================================================
unique(resparam$scenario)

df<-reshape2::melt(resparam, id.vars=c("parameter","iteration","scenario","method","model","by", 
                                      "convergence","conv_warning","pbias","bias"))

#df_alpha<-df[df$parameter%in%c("alpha"),]
df$variable<-factor(df$variable,levels=c("median","mode", "sim"))

df$method[df$method=="MLE"]<-"MAP"
unique(df$method)
df$scenario<-factor(df$scenario,levels=c("decLinearProd_highERLowError",
                                         "decLinearProd_ShiftERLowError",
                                         "decLinearProd_lowERLowError", 
                                         "incLinearProd_highERLowError",  
                                         "incLinearProd_ShiftERLowError",
                                         "incLinearProd_lowERLowError", 
                                         "highERHighError",               
                                         "highERLowError",               
                                         "lowERHighError",                
                                         "lowERLowError",                
                                          "ShiftERHighError",              
                                          "ShiftERLowError",               
                                          "trendERHighError",              
                                          "trendERLowError" )  ) 


df$ERtrend<-case_match(df$scenario,"decLinearProd_highERLowError"~"highER",
                                   "decLinearProd_ShiftERLowError"~"ShiftER", 
                                   "decLinearProd_lowERLowError"~"lowER",
                                   "incLinearProd_highERLowError"~"highER",  
                                   "incLinearProd_ShiftERLowError"~"ShiftER",
                                   "incLinearProd_lowERLowError"~"lowER",  
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
                                   "decLinearProd_lowERLowError"~"LowError",
                                   "incLinearProd_lowERLowError"~"LowError",
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
                                   "decLinearProd_lowERLowError"~"decLinear",
                                   "incLinearProd_highERLowError"~"incLinear",  
                                   "incLinearProd_ShiftERLowError"~"incLinear", 
                                   "incLinearProd_lowERLowError"~"incLinear", 
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


summarydf  <- df %>%
   group_by(scenario,parameter,
    method,model,by,variable,dynamics,
    ERerror,ERtrend) %>%
   reframe(qs = quantile(value, c(0.025, .5, 0.975),na.rm=T), prob = c("lower","median", "upper"))



summarydf <- reshape2::dcast(data=summarydf,  
    scenario + parameter + method + model + by + variable + dynamics + ERerror + ERtrend  ~prob, 
    value.var= "qs",fun.aggregate=mean)


summarydf$model2<-case_match(summarydf$model,
     "simple"~"stationary",
    "autocorr"~"autocorr", 
     "rwa"~ "rw.a",
     "hmma" ~"hmm.a",
      "rwb" ~"rw.b",
      "hmmb" ~"hmm.b", 
       "rwab"~"rw.ab",
        "hmmab"~ "hmm.ab" )
summarydf$model2<-factor(summarydf$model2, levels=c("stationary","autocorr","rw.a","hmm.a","rw.b","hmm.b","rw.ab", "hmm.ab"))

summarydf_alpha_sim1<-summarydf[summarydf$parameter=="alpha"&
                                summarydf$variable=="sim"&
                                summarydf$scenario%in%c("decLinearProd_highERLowError",
                                  "decLinearProd_ShiftERLowError", 
                                  "incLinearProd_highERLowError",  
                                  "incLinearProd_ShiftERLowError",
                                  "decLinearProd_lowERLowError",
                                         "incLinearProd_lowERLowError"),
                                ]

summarydf_alpha1<-summarydf[summarydf$parameter=="alpha"&
                            summarydf$variable=="mode"&
                            summarydf$scenario%in%c("decLinearProd_highERLowError",
                                         "decLinearProd_ShiftERLowError", 
                                         "incLinearProd_highERLowError",  
                                         "incLinearProd_ShiftERLowError",
                                         "decLinearProd_lowERLowError",
                                         "incLinearProd_lowERLowError"),
                            ]

scenlab1<-c("highER","ShiftER","lowER", "highER", "ShiftER", "lowER")
names(scenlab1) <- c("decLinearProd_highERLowError",
                                         "decLinearProd_ShiftERLowError",
                                         "decLinearProd_lowERLowError", 
                                         "incLinearProd_highERLowError",  
                                         "incLinearProd_ShiftERLowError",
                                         "incLinearProd_lowERLowError")


er_alpha1<-ggplot() + 
geom_pointrange(data=summarydf_alpha1,aes(x=by-54,y= median,ymin = lower, ymax = upper, col=method),alpha=.6)+
geom_line(data=summarydf_alpha_sim1,aes(x=by-54,y= median),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
coord_cartesian(ylim = c(0.2,2.7))+ 
mytheme+ 
ylab(expression(log(alpha))) +
xlab("year") +
facet_grid(dynamics+scenario~model2, scales="free_y",labeller = labeller(scenario= scenlab1))
er_alpha1
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/ER/compareMCMC_MLE_alpha_trendascn.png",
    plot=er_alpha1, width = 15,height = 8)


summarydf_alpha_sim2<-summarydf[summarydf$parameter=="alpha"&
                                summarydf$variable=="sim"&
                                summarydf$scenario%in%c("highERHighError",                              
                                         "lowERHighError",                           
                                          "ShiftERHighError",                                                                   
                                          "trendERHighError"),
                                ]

summarydf_alpha2<-summarydf[summarydf$parameter=="alpha"&
                            summarydf$variable=="mode"&
                            summarydf$scenario%in%c("highERHighError",                              
                                         "lowERHighError",                           
                                          "ShiftERHighError",                                                                   
                                          "trendERHighError"),
                            ]

                                      
scenlab2<-c("highER","lowER", "ShiftER","trendER")
names(scenlab2) <- c( "highERHighError",                              
                                         "lowERHighError",                           
                                          "ShiftERHighError",                                                                   
                                          "trendERHighError")

er_alpha2<-ggplot() + 
geom_pointrange(data=summarydf_alpha2,aes(x=by-54,y= median,ymin = lower, ymax = upper, col=method),alpha=.6)+
geom_line(data=summarydf_alpha_sim2,aes(x=by-54,y= median),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
coord_cartesian(ylim = c(0.2,2.7))+ 
mytheme + 
ylab(expression(log(alpha))) +
xlab("year") +
facet_grid(dynamics+scenario~model2, scales="free_y",labeller = labeller(scenario= scenlab2))
er_alpha2
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/ER/compareMCMC_MLE_alpha_higherscn.png",
    plot=er_alpha2, width = 15,height = 6)
#MLE estimates are less biased and higher than MCMC



summarydf_alpha_sim3<-summarydf[summarydf$parameter=="alpha"&
                                summarydf$variable=="sim"&
                                summarydf$scenario%in%c("highERLowError", 
                                         "lowERLowError",  
                                          "ShiftERLowError",      
                                          "trendERLowError"),
                                ]

summarydf_alpha3<-summarydf[summarydf$parameter=="alpha"&
                            summarydf$variable=="mode"&
                            summarydf$scenario%in%c("highERLowError", 
                                         "lowERLowError",  
                                          "ShiftERLowError",      
                                          "trendERLowError"),
                            ]



                                     
scenlab3<-c("highER","lowER", "ShiftER","trendER")
names(scenlab3) <- c( "highERLowError", 
                                         "lowERLowError",  
                                          "ShiftERLowError",      
                                          "trendERLowError" )



er_alpha3<-ggplot() + 
geom_pointrange(data=summarydf_alpha3,aes(x=by-54,y= median,ymin = lower, ymax = upper, col=method),alpha=.6)+
geom_line(data=summarydf_alpha_sim3,aes(x=by-54,y= median),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
coord_cartesian(ylim = c(0.2,2.7))+ 
mytheme + 
ylab(expression(log(alpha))) +
xlab("year") +
facet_grid(dynamics+scenario~model2, scales="free_y",labeller = labeller(scenario= scenlab3))
er_alpha3
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/ER/compareMCMC_MLE_alpha_lowerscn.png",
    plot=er_alpha3, width = 15,height = 6)


#=======================================================
#smax estimates


summarydf_smax_sim1<-summarydf[summarydf$parameter=="smax"&
                                summarydf$variable=="sim"&
                                summarydf$scenario%in%c("decLinearProd_highERLowError",
                                  "decLinearProd_ShiftERLowError", 
                                  "incLinearProd_highERLowError",  
                                  "incLinearProd_ShiftERLowError",
                                  "decLinearProd_lowERLowError",
                                         "incLinearProd_lowERLowError"),
                                ]

summarydf_smax1<-summarydf[summarydf$parameter=="smax"&
                            summarydf$variable=="mode"&
                            summarydf$scenario%in%c("decLinearProd_highERLowError",
                                         "decLinearProd_ShiftERLowError", 
                                         "incLinearProd_highERLowError",  
                                         "incLinearProd_ShiftERLowError",
                                         "decLinearProd_lowERLowError",
                                         "incLinearProd_lowERLowError"),
                            ]

scenlab1<-c("highER","ShiftER","lowER", "highER", "ShiftER", "lowER")
names(scenlab1) <- c("decLinearProd_highERLowError",
                                         "decLinearProd_ShiftERLowError",
                                         "decLinearProd_lowERLowError", 
                                         "incLinearProd_highERLowError",  
                                         "incLinearProd_ShiftERLowError",
                                         "incLinearProd_lowERLowError")


er_smax1<-ggplot() + 
geom_pointrange(data=summarydf_smax1,aes(x=by-54,y=median,ymin = lower, ymax = upper,  col=method),alpha=.6)+
geom_line(data=summarydf_smax_sim1,aes(x=by-54,y= median),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
mytheme + 
ylab(expression(S[max])) +
xlab("year") +
coord_cartesian(ylim = c(60000,400000))+ 
facet_grid(dynamics+scenario~model2, scales="free_y",labeller = labeller(scenario= scenlab1))
er_smax1
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/ER/compareMCMC_MLE_smax_trendascn.png",
    plot=er_smax1, width = 15,height = 8)



summarydf_smax_sim2<-summarydf[summarydf$parameter=="smax"&
                                summarydf$variable=="sim"&
                                summarydf$scenario%in%c("highERHighError",                              
                                         "lowERHighError",                           
                                          "ShiftERHighError",                                                                   
                                          "trendERHighError"),
                                ]

summarydf_smax2<-summarydf[summarydf$parameter=="smax"&
                            summarydf$variable=="mode"&
                            summarydf$scenario%in%c("highERHighError",                              
                                         "lowERHighError",                           
                                          "ShiftERHighError",                                                                   
                                          "trendERHighError"),
                            ]


er_smax2<-ggplot() + 
geom_pointrange(data=summarydf_smax2,aes(x=by-54,y=median,ymin = lower, ymax = upper, col=method),alpha=.6)+
geom_line(data=summarydf_smax_sim2,aes(x=by-54,y=median),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
mytheme + 
ylab(expression(S[max])) +
xlab("year") +
coord_cartesian(ylim = c(60000,400000))+ 
facet_grid(dynamics+scenario~model2, scales="free_y",labeller = labeller(scenario= scenlab2))
er_smax2
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/ER/compareMCMC_MLE_smax_higherscn.png",
    plot=er_smax2, width = 15,height = 6)




summarydf_smax_sim3<-summarydf[summarydf$parameter=="smax"&
                                summarydf$variable=="sim"&
                                summarydf$scenario%in%c("highERLowError",                              
                                         "lowERLowError",                           
                                          "ShiftERLowError",                                                                   
                                          "trendERLowError"),
                                ]

summarydf_smax3<-summarydf[summarydf$parameter=="smax"&
                            summarydf$variable=="mode"&
                            summarydf$scenario%in%c("highERLowError",                              
                                         "lowERLowError",                           
                                          "ShiftERLowError",                                                                   
                                          "trendERLowError"),
                            ]



er_smax3<-ggplot() + 
geom_pointrange(data=summarydf_smax3,aes(x=by-54,y=median,ymin =lower, ymax =upper, col=method),alpha=.6)+
geom_line(data=summarydf_smax_sim3,aes(x=by-54,y=median),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
mytheme + 
ylab(expression(S[max])) +
xlab("year") +
coord_cartesian(ylim = c(60000,400000))+ 
facet_grid(dynamics+scenario~model2, scales="free_y",labeller = labeller(scenario= scenlab3))
er_smax3
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/ER/compareMCMC_MLE_smax_lowerscn.png",
    plot=er_smax3, width = 15,height = 6)


#=======================================================
#smsy estimates




summarydf_smsy_sim1<-summarydf[summarydf$parameter=="smsy"&
                                summarydf$variable=="sim"&
                                summarydf$scenario%in%c("decLinearProd_highERLowError",
                                  "decLinearProd_ShiftERLowError", 
                                  "incLinearProd_highERLowError",  
                                  "incLinearProd_ShiftERLowError",
                                  "decLinearProd_lowERLowError",
                                         "incLinearProd_lowERLowError"),
                                ]

summarydf_smsy1<-summarydf[summarydf$parameter=="smsy"&
                            summarydf$variable=="mode"&
                            summarydf$scenario%in%c("decLinearProd_highERLowError",
                                         "decLinearProd_ShiftERLowError", 
                                         "incLinearProd_highERLowError",  
                                         "incLinearProd_ShiftERLowError",
                                         "decLinearProd_lowERLowError",
                                         "incLinearProd_lowERLowError"),
                            ]



er_smsy1<-ggplot() + 
geom_pointrange(data=summarydf_smsy1,aes(x=by-54,y=median,ymin = lower, ymax = upper,  col=method),alpha=.6)+
geom_line(data=summarydf_smsy_sim1,aes(x=by-54,y= median),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
mytheme + 
ylab(expression(S[MSY])) +
xlab("year") +
coord_cartesian(ylim = c(20000,150000))+ 
facet_grid(dynamics+scenario~model2, scales="free_y",labeller = labeller(scenario= scenlab1))
er_smsy1
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/ER/compareMCMC_MLE_smsy_trendascn.png",
    plot=er_smsy1, width = 15,height = 8)




summarydf_smsy_sim2<-summarydf[summarydf$parameter=="smsy"&
                                summarydf$variable=="sim"&
                                summarydf$scenario%in%c( "highERHighError",                              
                                         "lowERHighError",                           
                                          "ShiftERHighError",                                                                   
                                          "trendERHighError"),
                                ]

summarydf_smsy2<-summarydf[summarydf$parameter=="smsy"&
                            summarydf$variable=="mode"&
                            summarydf$scenario%in%c( "highERHighError",                              
                                         "lowERHighError",                           
                                          "ShiftERHighError",                                                                   
                                          "trendERHighError"),
                            ]

er_smsy2<-ggplot() + 
geom_pointrange(data=summarydf_smsy2,aes(x=by-54,y=median,ymin = lower, ymax = upper,  col=method),alpha=.6)+
geom_line(data=summarydf_smsy_sim2,aes(x=by-54,y= median),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
mytheme + 
ylab(expression(S[MSY])) +
xlab("year") +
coord_cartesian(ylim = c(20000,150000))+ 
facet_grid(dynamics+scenario~model2, scales="free_y",labeller = labeller(scenario= scenlab2))
er_smsy2
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/ER/compareMCMC_MLE_smsy_higherscn.png",
    plot=er_smsy2, width = 15,height = 6)



summarydf_smsy_sim3<-summarydf[summarydf$parameter=="smsy"&
                                summarydf$variable=="sim"&
                                summarydf$scenario%in%c( "highERLowError",                              
                                         "lowERLowError",                           
                                          "ShiftERLowError",                                                                   
                                          "trendERLowError"),
                                ]

summarydf_smsy3<-summarydf[summarydf$parameter=="smsy"&
                            summarydf$variable=="mode"&
                            summarydf$scenario%in%c( "highERLowError",                              
                                         "lowERLowError",                           
                                          "ShiftERLowError",                                                                   
                                          "trendERLowError"),
                            ]

er_smsy3<-ggplot() + 
geom_pointrange(data=summarydf_smsy3,aes(x=by-54,y=median,ymin = lower, ymax = upper,  col=method),alpha=.6)+
geom_line(data=summarydf_smsy_sim3,aes(x=by-54,y= median),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
mytheme + 
ylab(expression(S[MSY])) +
xlab("year") +
coord_cartesian(ylim = c(20000,150000))+ 
facet_grid(dynamics+scenario~model2, scales="free_y",labeller = labeller(scenario= scenlab3))
er_smsy3
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/ER/compareMCMC_MLE_smsy_lowerscn.png",
    plot=er_smsy3, width = 15,height = 6)





#-----------------------------------------------------------------------------------
#summary plots 


head(summarydf_smsy_sim_redux)
unique(summarydf$ERtrend)



summarydf_smsy_sim_redux<-summarydf[summarydf$parameter=="smsy"&
                                summarydf$variable=="sim"&
                                summarydf$scenario%in%c("highERLowError", 
                                                        "ShiftERLowError",
                                                        "lowERLowError", 
                                                        "decLinearProd_highERLowError",
                                                        "incLinearProd_highERLowError", 
                                                        "decLinearProd_ShiftERLowError", 
                                                        "incLinearProd_ShiftERLowError",
                                                        "decLinearProd_lowERLowError", 
                                                        "incLinearProd_lowERLowError")&
                                summarydf$model%in%c("autocorr", "rwa")&
                                summarydf$method=="MAP",]

summarydf_smsy_sim_redux$ERtrend<-factor(summarydf_smsy_sim_redux$ERtrend,levels=c("highER", "ShiftER", "lowER"))

summarydf_smsy_sim_redux$dynamics2<-case_match(summarydf_smsy_sim_redux$dynamics,
    "decLinear"~"decline",
    "incLinear"~"increase",
    "stationary"~"stationary")



summarydf_smsy_redux<-summarydf[summarydf$variable=="mode"&
                                summarydf$parameter=="smsy"&
                                summarydf$scenario%in%c("highERLowError", 
                                                        "ShiftERLowError", 
                                                        "lowERLowError",
                                                        "decLinearProd_highERLowError",
                                                        "incLinearProd_highERLowError", 
                                                        "decLinearProd_ShiftERLowError", 
                                                        "incLinearProd_ShiftERLowError",
                                                        "decLinearProd_lowERLowError", 
                                                        "incLinearProd_lowERLowError")&
                                    summarydf$model%in%c("autocorr", "rwa")&
                                    summarydf$method=="MAP",]


summarydf_smsy_redux$ERtrend<-factor(summarydf_smsy_redux$ERtrend,levels=c("highER", "ShiftER", "lowER"))


head(summarydf_smsy_redux)
unique(summarydf_smsy_redux$dynamics)

summarydf_smsy_redux$dynamics2<-case_match(summarydf_smsy_redux$dynamics,
    "decLinear"~"decline",
    "incLinear"~"increase",
    "stationary"~"stationary")

summarydf_smsy_redux$dynamics


psmsy_highERscn_line<-ggplot() + 
geom_pointrange(data=summarydf_smsy_redux,aes(x=by-54,y= median,ymin = lower, ymax = upper, color=model2),alpha=.9)+
geom_line(data=summarydf_smsy_sim_redux,aes(x=by-54,y=median),color="black", alpha=.8,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8,option = "E") +
scale_fill_viridis_d(begin=.1, end=.8,option = "E") +
coord_cartesian(ylim = c(20000,180000))+ 
mytheme+ 
ylab(expression(S[MSY])) +
xlab("year") +
facet_grid(dynamics2+ERtrend ~model2, scales="free_y")+
theme(strip.text = element_text(
    size =11))
psmsy_highERscn_line



df$model2<-case_match(df$model,
     "simple"~"stationary",
    "autocorr"~"autocorr", 
     "rwa"~ "rw.a",
     "hmma" ~"hmm.a",
      "rwb" ~"rw.b",
      "hmmb" ~"hmm.b", 
       "rwab"~"rw.ab",
        "hmmab"~ "hmm.ab" )
df$model2<-factor(df$model2, levels=c("stationary","autocorr","rw.a","hmm.a","rw.b","hmm.b","rw.ab", "hmm.ab"))

df_smsy_est_redux<- df[df$parameter=="smsy"&
                    df$variable=="mode"&
                    df$scenario%in%c("highERLowError",  
                                     "ShiftERLowError", 
                                     "lowERLowError",
                                     "decLinearProd_highERLowError",
                                     "incLinearProd_highERLowError",
                                     "decLinearProd_ShiftERLowError", 
                                     "incLinearProd_ShiftERLowError",
                                     "decLinearProd_lowERLowError", 
                                     "incLinearProd_lowERLowError")&
                    df$model%in%c("autocorr", "rwa")&
                    df$method=="MAP",]
df_smsy_est_redux$ERtrend<-factor(df_smsy_est_redux$ERtrend,levels=c("highER", "ShiftER", "lowER"))


df_smsy_est_redux$dynamics2<-case_match(df_smsy_est_redux$dynamics,
    "decLinear"~"decline",
    "incLinear"~"increase",
    "stationary"~"stationary")


psmsy_erscn_violin_abs<-ggplot(df_smsy_est_redux) + 
geom_violin(aes(x=model2,y=abs(bias), fill=model2), scale="width", trim=TRUE, alpha=.7,adjust = 1.8)+
geom_boxplot(aes(x=model2,y=abs(bias), fill=model2),outlier.shape = NA, width=0.1)+
geom_hline(yintercept=0,color="black", alpha=.6,linewidth=1.2)+
 scale_fill_viridis_d(begin=.1, end=.8,,option = "E") +
 facet_grid(dynamics2 ~ERtrend, scales="free_y")+
 ylab(expression("absolute bias in"~S[MSY])) +
 coord_cartesian(ylim = c(0,75000))+ 
mytheme
psmsy_erscn_violin_abs

multi.page.abs.smsy_erscn <- ggarrange(psmsy_highERscn_line, psmsy_erscn_violin_abs,
                        nrow = 1, ncol = 2,
                        common.legend = TRUE,
                        legend="bottom")
multi.page.abs.smsy_erscn


multi.page.abs.smsy_erscn_title<-annotate_figure(multi.page.abs.smsy_erscn, 
    top = text_grob(expression("Sensitivity to exploitation history and time-varying"~log(alpha)), 
               face = "bold" , size = 14))




ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/summary_figs/bias_trends_violins_smsy_ERscn.png",
    plot=multi.page.abs.smsy_erscn_title, width = 12,height = 9, bg = "white")


#=======================================================================================================
#summary plots 
#with ow ER scenarios
unique(summarydf_smsy_sim$scenario)



summarydf_smsy_sim_redux2<-summarydf[summarydf$parameter=="smsy"&
                                     summarydf$variable=="sim"&
                                    summarydf$scenario%in%c("highERLowError",  
                                                             "ShiftERLowError", 
                                                             "lowERLowError",
                                                             "decLinearProd_highERLowError",
                                                              "incLinearProd_highERLowError", 
                                                              "decLinearProd_ShiftERLowError", 
                                                              "incLinearProd_ShiftERLowError",
                                                              "decLinearProd_lowERLowError", 
                                                              "incLinearProd_lowERLowError")&
                                    summarydf$model%in%c("autocorr", "rwa")&
                                    summarydf$method=="MLE",]


summarydf_smsy_sim_redux2$ERtrend<-factor(summarydf_smsy_sim_redux2$ERtrend,levels=c("highER", "ShiftER", "lowER"))


summarydf_smsy_redux2<-summarydf[summarydf$parameter=="smsy"&
                                summarydf$variable=="mode"&
                                summarydf$scenario%in%c("highERLowError", 
                                                        "ShiftERLowError", 
                                                        "lowERLowError",
                                                             "decLinearProd_highERLowError",
                                                             "incLinearProd_highERLowError",
                                                             "decLinearProd_ShiftERLowError", 
                                                             "incLinearProd_ShiftERLowError",
                                                             "decLinearProd_lowERLowError", 
                                                             "incLinearProd_lowERLowError" )&
                                summarydf$model%in%c("autocorr", "rwa")&
                                summarydf$method=="MLE",]


summarydf_smsy_redux2$ERtrend<-factor(summarydf_smsy_redux2$ERtrend,levels=c("highER", "ShiftER", "lowER"))


psmsy_highERscn_line2<-ggplot() + 
geom_pointrange(data=summarydf_smsy_redux2,aes(x=by,y=median,ymin =lower, ymax =upper, color=model),alpha=.9)+
geom_line(data=summarydf_smsy_sim_redux2,aes(x=by,y=median),color="black", alpha=.8,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8,option = "E") +
scale_fill_viridis_d(begin=.1, end=.8,option = "E") +
coord_cartesian(ylim = c(20000,180000))+ 
mytheme+ 
ylab("Smsy") +
xlab("year") +
facet_grid(dynamics+ERtrend ~model, scales="free_y")
psmsy_highERscn_line2


head(df$scenario)
df_smsy_est_redux2<- df[df$parameter=="smsy"&
                        df$variable=="mode"&
                        df$scenario%in%c("highERLowError", 
                                         "ShiftERLowError", 
                                         "lowERLowError",
                                         "decLinearProd_highERLowError",
                                         "incLinearProd_highERLowError",
                                         "decLinearProd_ShiftERLowError", 
                                         "incLinearProd_ShiftERLowError",
                                         "decLinearProd_lowERLowError", 
                                         "incLinearProd_lowERLowError" )&
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


summarydf_smax_sim_redux<-summarydf[summarydf$parameter=="smax"&
                                    summarydf$variable=="mode"&
                                    summarydf$scenario%in%c("highERLowError", 
                                                            "ShiftERLowError", 
                                                            "lowERLowError",
                                                            "decLinearProd_highERLowError",
                                                            "incLinearProd_highERLowError",
                                                            "decLinearProd_ShiftERLowError", 
                                                            "incLinearProd_ShiftERLowError",
                                                            "decLinearProd_lowERLowError", 
                                                            "incLinearProd_lowERLowError" )&
                                    summarydf$model%in%c("autocorr", "rwa")&
                                    summarydf$method=="MLE",]


summarydf_smax_sim_redux$ERtrend<-factor(summarydf_smax_sim_redux$ERtrend,levels=c("highER", "ShiftER", "lowER"))

summarydf_smax_redux<-summarydf[summarydf$parameter=="smax"&
                                summarydf$variable=="mode"&
                                summarydf$scenario%in%c("highERLowError", 
                                                        "ShiftERLowError", 
                                                        "lowERLowError",
                                                        "decLinearProd_highERLowError",
                                                        "incLinearProd_highERLowError",
                                                        "decLinearProd_ShiftERLowError", 
                                                        "incLinearProd_ShiftERLowError",
                                                        "decLinearProd_lowERLowError", 
                                                        "incLinearProd_lowERLowError" )&
                                summarydf$model%in%c("autocorr", "rwa") &
                                summarydf$method=="MLE",]




summarydf_smax_redux$ERtrend<-factor(summarydf_smax_redux$ERtrend,levels=c("highER", "ShiftER", "lowER"))


psmax_highERscn_line<-ggplot() + 
geom_pointrange(data=summarydf_smax_redux,aes(x=by,y=median,ymin =lower, ymax = upper, color=model),alpha=.9)+
geom_line(data=summarydf_smax_sim_redux,aes(x=by,y=median),color="black", alpha=.8,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8,option = "E") +
scale_fill_viridis_d(begin=.1, end=.8,option = "E") +
coord_cartesian(ylim = c(50000,400000))+ 
mytheme+ 
ylab("Smax") +
xlab("year") +
facet_grid(dynamics+ERtrend ~model, scales="free_y")
psmax_highERscn_line


head(df$scenario)
df_smax_est_redux<- df[df$parameter=="smax"&
                       df$variable=="mode"&
                       df$scenario%in%c("highERLowError", 
                                        "ShiftERLowError", 
                                        "lowERLowError",
                                        "decLinearProd_highERLowError",
                                        "incLinearProd_highERLowError",
                                        "decLinearProd_ShiftERLowError", 
                                        "incLinearProd_ShiftERLowError",
                                        "decLinearProd_lowERLowError", 
                                        "incLinearProd_lowERLowError")&
                        df$model%in%c("autocorr", "rwa")&
                        df$method=="MLE",]


df_smax_est_redux$ERtrend<-factor(df_smax_est_redux$ERtrend,levels=c("highER", "ShiftER", "lowER"))


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


summarydf_alpha_sim_redux<-summarydf[summarydf$parameter=="alpha"&
                                    summarydf$variable=="sim"&
                                    summarydf$scenario%in%c("highERLowError", 
                                                            "ShiftERLowError", 
                                                            "lowERLowError",
                                                            "decLinearProd_highERLowError",
                                                            "incLinearProd_highERLowError",
                                                            "decLinearProd_ShiftERLowError", 
                                                            "incLinearProd_ShiftERLowError",
                                                            "decLinearProd_lowERLowError", 
                                                            "incLinearProd_lowERLowError")&
                                    summarydf$model%in%c("autocorr", "rwa")&
                                    summarydf$method=="MLE",]

summarydf_alpha_sim_redux$ERtrend<-factor(summarydf_alpha_sim_redux$ERtrend,levels=c("highER", "ShiftER", "lowER"))

summarydf_alpha_redux<-summarydf[summarydf$parameter=="alpha"&
                                       summarydf$variable=="mode"&
                                       summarydf$scenario%in%c("highERLowError", 
                                                            "ShiftERLowError", 
                                                            "lowERLowError",
                                                            "decLinearProd_highERLowError",
                                                            "incLinearProd_highERLowError",
                                                            "decLinearProd_ShiftERLowError", 
                                                            "incLinearProd_ShiftERLowError",
                                                            "decLinearProd_lowERLowError", 
                                                            "incLinearProd_lowERLowError")&
                                       summarydf$model%in%c("autocorr", "rwa") &
                                       summarydf$method=="MLE",]


summarydf_alpha_redux$ERtrend<-factor(summarydf_alpha_redux$ERtrend,levels=c("highER", "ShiftER", "lowER"))



palpha_highERscn_line<-ggplot() + 
geom_pointrange(data=summarydf_alpha_redux,aes(x=by,y=median,ymin = lower, ymax = upper, color=model),alpha=.9)+
geom_line(data=summarydf_alpha_sim_redux,aes(x=by,y=median),color="black", alpha=.8,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8,option = "E") +
scale_fill_viridis_d(begin=.1, end=.8,option = "E") +
coord_cartesian(ylim = c(0,2.3))+ 
mytheme+ 
ylab("log(alpha)") +
xlab("year") +
facet_grid(dynamics+ERtrend ~model, scales="free_y")
palpha_highERscn_line


df_alpha_est_redux<- df[df$parameter=="alpha"&df$variable=="mode"&
                        df$scenario%in%c("highERLowError", 
                                        "ShiftERLowError", 
                                        "lowERLowError",
                                        "decLinearProd_highERLowError",
                                        "incLinearProd_highERLowError",
                                        "decLinearProd_ShiftERLowError", 
                                        "incLinearProd_ShiftERLowError",
                                        "decLinearProd_lowERLowError", 
                                        "incLinearProd_lowERLowError")&
                        df$model%in%c("autocorr", "rwa")&
                        df$method=="MLE",]



df_alpha_est_redux$ERtrend<-factor(df_alpha_est_redux$ERtrend,levels=c("highER", "ShiftER", "lowER"))


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

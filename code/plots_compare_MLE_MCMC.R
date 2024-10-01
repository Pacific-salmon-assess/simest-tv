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
    theme_classic(13)+
        theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),
            panel.border = element_rect(fill = NA, color = "black"), legend.title = element_blank(),
            legend.position="bottom", strip.text = element_text( size=13),
            axis.text=element_text(face="bold",size=13),axis.title = element_text(face="bold",size=13),
            legend.text=element_text(size=13),
            plot.title = element_text(face = "bold", hjust = 0.5,size=13))
)

source("code/read_base_data.R")

#========================================================================================================
#base case
#read in data

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
      "autocorr"~"autocorrelation",
      "sigmaShift"~"sigma", 
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


df$scentype<-factor(df$scentype, levels=c("stationary","autocorrelation","sigma","log(alpha)",
    "Smax","both" )  ) 


df$scendesc<-dplyr::case_match(df$scenario, 
      "stationary"~".",
      "autocorr"~".",
      "sigmaShift"~"shift up", 
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


df$scendesc<-factor(df$scendesc, levels=c(".","shift up","linear decline",
    "sine trend","shift up + down","shift down + up","shift down",
    "regime both","trend & shift") ) 
     

df$scentrend<-dplyr::case_match(df$scenario, 
      "stationary"~"stationary",
      "autocorr"~"autocorr",
      "sigmaShift"~"sigma shift", 
      "decLinearProd"~"trend",
      "sineProd"~"trend",
      "regimeProd"~"shift",
      "shiftProd"~"shift",
      "decLinearCap"~"trend",
      "regimeCap"~"shift",
      "shiftCap"~"shift", 
      "regimeProdCap"~"shift",
      "decLinearProdshiftCap"~"trend & shift"
      )   


df$model<-factor(df$model,levels=c("simple",
                                   "autocorr", 
                                   "rwa",
                                   "hmma",
                                   "rwb",
                                   "hmmb",
                                   "rwab",
                                   "hmmab"  ))


df$model2<-dplyr::case_match(df$model, 
     "simple"~"stationary",
     "autocorr"~"autocorr", 
     "rwa"~"rw.a",
     "hmma"~"hmm.a",
     "rwb"~"rw.b",
     "hmmb"~"hmm.b",
     "rwab"~"rw.ab",
     "hmmab"~"hmm.ab" 
      )   


df$model2<-factor(df$model2,levels=c("stationary",
                                   "autocorr", 
                                   "rw.a",
                                   "hmm.a",
                                   "rw.b",
                                   "hmm.b",
                                   "rw.ab",
                                   "hmm.ab"  ))



summarydf  <- df %>%
   group_by(scenario,parameter,
    method,model,model2,by,variable,scencode,scentype,scentrend,scendesc) %>%
   reframe(qs = quantile(value, c(0.025, .5, 0.975),na.rm=T), prob = c("lower","median", "upper"))



summarydf <- reshape2::dcast(data=summarydf,  
    scenario + parameter + method + model +model2 + by + variable + scencode + scentype +scentrend +scendesc~ prob, 
    value.var= "qs",fun.aggregate=mean)

#==================================================================================
#mean CV estimates
head(resparam)

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

meancvdf$model2<-case_match(meancvdf$model,
     "simple"~"stationary",
    "autocorr"~"autocorr", 
     "rwa"~ "rw.a",
     "hmma" ~"hmm.a",
      "rwb" ~"rw.b",
      "hmmb" ~"hmm.b", 
       "rwab"~"rw.ab",
        "hmmab"~ "hmm.ab" )


meancvdf$scendesc<-dplyr::case_match(meancvdf$scenario, 
      "stationary"~".",
      "autocorr"~".",
      "sigmaShift"~"shift up", 
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

meancvdf$scendesc <- factor(meancvdf$scendesc, levels=c(".","shift up","linear decline",
    "sine trend","shift up + down","shift down + up","shift down",
    "regime both","trend & shift") ) 
     

meancvdf$scentype<-dplyr::case_match(meancvdf$scenario, 
       "stationary"~"stationary",
      "autocorr"~"autocorrelation",
      "sigmaShift"~"sigma", 
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

meancvdf$scentype<-factor(meancvdf$scentype, levels=c("stationary","autocorrelation",
    "sigma","log(alpha)","Smax","both" )  ) 

meancvdf$labels<-round(meancvdf$x,2)

#mean absolute percent bias
abspbiasdf <- aggregate(resparam$pbias, list(parameter=resparam$parameter,
                               scenario=resparam$scenario,
                               method=resparam$method,
                               model=resparam$model,
                               by=resparam$by ), function(x){median(abs(x),na.rm=T)})

meanabspbiasdf <- aggregate(abspbiasdf$x, list(parameter=abspbiasdf$parameter,
                               scenario=abspbiasdf$scenario,
                               method=abspbiasdf$method,
                               model=abspbiasdf$model
                              ), function(x){mean(x,na.rm=T)})
         
meanabspbiasdf$model2<-case_match(meanabspbiasdf$model,
     "simple"~"stationary",
    "autocorr"~"autocorr", 
     "rwa"~ "rw.a",
     "hmma" ~"hmm.a",
      "rwb" ~"rw.b",
      "hmmb" ~"hmm.b", 
       "rwab"~"rw.ab",
        "hmmab"~ "hmm.ab" )


#===================================================================================

summarydf_alpha<-summarydf[summarydf$parameter=="logalpha"&
                            summarydf$variable=="mode",
                            ]

summarydf_alpha_sim<-summarydf[summarydf$parameter=="logalpha"&
                                summarydf$variable=="sim",
                                ]

meancvdf_alpha<-meancvdf[meancvdf$parameter=="logalpha", ]

alphabase<-ggplot() + 
geom_pointrange(data=summarydf_alpha,aes(x=by-54,y= median,ymin = lower, ymax = upper, col=method),alpha=.6)+
#geom_line(data=summarydf_alpha_sim,aes(x=by-54,y= median),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
coord_cartesian(ylim = c(0.2,2.3)) + 
mytheme+ 
ylab(expression(log(alpha))) +
xlab("year") +
facet_grid(scentype+scendesc~model2, scales="free_y")
alphabase
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/base/compareMCMC_MLE_alpha_base.png",
    plot=alphabase, width = 15,height = 18 )



#=======================================================
#b estimates

summarydf_smax<-summarydf[summarydf$parameter=="smax"&
                            summarydf$variable=="mode",
                            ]

summarydf_smax_sim<-summarydf[summarydf$parameter=="smax"&
                                summarydf$variable=="sim",
                                ]

smaxbase<-ggplot() + 
geom_pointrange(data=summarydf_smax,aes(x=by-54,y=median/1000,ymin =lower/1000,
   ymax = upper/1000, col=method),alpha=.6)+
geom_line(data=summarydf_smax_sim,aes(x=by-54,y= median/1000),color="black", 
    alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
mytheme + 
ylab(expression(S[max]~"(1000s)")) +
xlab("year") +
coord_cartesian(ylim = c(60000,400000)/1000)+ 
facet_grid(scentype+scendesc~model2, scales="free_y")
smaxbase
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/base/compareMCMC_MLE_smax_base.png",
    plot=smaxbase, width = 15,height = 18)




#=======================================================
#smsy estimates


summarydf_smsy<-summarydf[summarydf$parameter=="smsy"&
                            summarydf$variable=="mode",
                            ]

summarydf_smsy_sim<-summarydf[summarydf$parameter=="smsy"&
                                summarydf$variable=="sim",
                                ]


smsybase<-ggplot() + 
geom_pointrange(data=summarydf_smsy,aes(x=by-54,y=median/1000,ymin =lower/1000,
      ymax = upper/1000, col=method),alpha=.6)+
geom_line(data=summarydf_smsy_sim,aes(x=by-54,y= median/1000), color="black", 
    alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
mytheme + 
ylab(expression(S[MSY]~"(1000s)")) +
xlab("year") +
coord_cartesian(ylim = c(20000,150000)/1000)+ 
facet_grid(scentype+scendesc~model2, scales="free_y")
smsybase
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/base/compareMCMC_MLE_smsy_base.png",
    plot=smsybase, width = 15,height = 18)


#=========================================================================================================\
# smsy-focusedredux plots

# alpha varies


meancvdf_redux_alpha<-meancvdf[meancvdf$parameter=="smsy"&
meancvdf$scenario%in%c("autocorr","decLinearProd", "shiftProd", "sineProd")&
meancvdf$model%in%c("autocorr", "rwa", "hmma")&
meancvdf$method=="MLE",]


meancvdf_redux_alpha$meancv<-round(meancvdf_redux_alpha$x,2)
meancvdf_redux_alpha$scenario2<-case_match(meancvdf_redux_alpha$scenario,
     "autocorr"~ "stationary",
    "decLinearProd"~ "linear decline",
     "shiftProd" ~ "shift decline" , 
     "sineProd" ~ "sine fluctuation")

meancvdf_redux_alpha$model2<-factor(meancvdf_redux_alpha$model2, levels=c("autocorr", "rw.a", "hmm.a"))

meancvdf_redux_alpha$scenario2<-factor(meancvdf_redux_alpha$scenario2,levels=c("stationary",
    "linear decline" ,"sine fluctuation","shift decline" ))

#absolute median percent bias
abspbiasdf<-aggregate(resparam$pbias, list(parameter=resparam$parameter,
                               scenario=resparam$scenario,
                               method=resparam$method,
                               model=resparam$model,
                               by=resparam$by ), function(x){median(abs(x),na.rm=T)})

meanabspbiasdf<-aggregate(abspbiasdf$x, list(parameter=abspbiasdf$parameter,
                               scenario=abspbiasdf$scenario,
                               method=abspbiasdf$method,
                               model=abspbiasdf$model
                              ), function(x){mean(x,na.rm=T)})

meanabspbiasdf$model2<-case_match(meanabspbiasdf$model,
     "simple"~"stationary",
    "autocorr"~"autocorr", 
     "rwa"~ "rw.a",
     "hmma" ~"hmm.a",
      "rwb" ~"rw.b",
      "hmmb" ~"hmm.b", 
       "rwab"~"rw.ab",
        "hmmab"~ "hmm.ab" )

meanabspbiasdf_redux_alpha<-meanabspbiasdf[meanabspbiasdf$parameter=="smsy"&
meanabspbiasdf$scenario%in%c("autocorr","decLinearProd", "shiftProd", "sineProd")&
meanabspbiasdf$model%in%c("autocorr", "rwa", "hmma")&
meanabspbiasdf$method=="MLE",]

meanabspbiasdf_redux_alpha$meanabspbias<-round(meanabspbiasdf_redux_alpha$x,2)
meanabspbiasdf_redux_alpha$scenario2<-case_match(meanabspbiasdf_redux_alpha$scenario,
     "autocorr"~ "stationary",
    "decLinearProd"~ "linear decline",
     "shiftProd" ~ "shift decline" , 
     "sineProd" ~ "sine fluctuation")
meanabspbiasdf_redux_alpha$model2<-factor(meanabspbiasdf_redux_alpha$model2, levels=c("autocorr", "rw.a", "hmm.a"))

meanabspbiasdf_redux_alpha$scenario2<-factor(meanabspbiasdf_redux_alpha$scenario2,levels=c("stationary",
    "linear decline" ,"sine fluctuation","shift decline" ))

#estimates
df_smsy_sim<- df[df$parameter%in%c("smsy")&df$variable=="sim",]
df_smsy_est<- df[df$parameter%in%c("smsy")&df$variable=="mode",]


summarydf_smsy<-aggregate(df_smsy_est$value,by=list(scenario=df_smsy_est$scenario, 
    method=df_smsy_est$method, 
    model2=df_smsy_est$model2,
    by=df_smsy_est$by ),
    function(x) {quantile(x,probs = c(0.025, .5, 0.975))})
summarydf_smsy<-do.call(data.frame, summarydf_smsy)

summarydf_smsy_sim<-aggregate(df_smsy_sim$value,by=list(scenario=df_smsy_sim$scenario, 
    method=df_smsy_sim$method, 
    model2=df_smsy_sim$model2,
    by=df_smsy_sim$by ),
    function(x) {unique(x)})

summarydf_smsy_sim_redux_alpha<-summarydf_smsy_sim[summarydf_smsy_sim$scenario%in%c("autocorr","decLinearProd", "shiftProd", "sineProd")&
summarydf_smsy_sim$model2%in%c("autocorr", "rw.a", "hmm.a")&summarydf_smsy_sim$method=="MLE",]

summarydf_smsy_sim_redux_alpha$scenario2<-case_match(summarydf_smsy_sim_redux_alpha$scenario,
     "autocorr"~ "stationary",
    "decLinearProd"~ "linear decline",
     "shiftProd" ~ "shift decline" , 
     "sineProd" ~ "sine fluctuation")

summarydf_smsy_sim_redux_alpha$scenario2<-factor(summarydf_smsy_sim_redux_alpha$scenario2, 
    levels=c("stationary",
     "linear decline" ,
     "sine fluctuation",
     "shift decline" ))

summarydf_smsy_redux_alpha<-summarydf_smsy[summarydf_smsy$scenario%in%c("autocorr", "decLinearProd", "shiftProd", "sineProd" )&
  summarydf_smsy$model2%in%c("autocorr", "rw.a", "hmm.a") &summarydf_smsy$method=="MLE",]


summarydf_smsy_redux_alpha$scenario2<-case_match(summarydf_smsy_redux_alpha$scenario,
     "autocorr"~ "stationary",
    "decLinearProd"~ "linear decline",
     "shiftProd" ~ "shift decline" , 
     "sineProd" ~ "sine fluctuation")

summarydf_smsy_redux_alpha$scenario2<-factor(summarydf_smsy_redux_alpha$scenario2, 
    levels=c("stationary",
     "linear decline",
     "sine fluctuation",
     "shift decline" ))

meancvdf_redux_alpha$labels<-paste("CV = ",meancvdf_redux_alpha$meancv)


psmsy_alphascn_line<-ggplot() + 
geom_pointrange(data=summarydf_smsy_redux_alpha,aes(x=by-54,y= x.50./10000,ymin = x.2.5./10000, ymax = x.97.5./10000, color=model2),alpha=.9)+
geom_line(data=summarydf_smsy_sim_redux_alpha,aes(x=by-54,y= x/10000),color="black", alpha=.8,linewidth=1.2)+
scale_color_viridis_d(begin=.2, end=.8,option = "E") +
scale_fill_viridis_d(begin=.2, end=.8,option = "E") +
coord_cartesian(ylim = c(2,17))+ 
mytheme+ 
ylab(expression(paste(S[MSY]~("10,000s")))) +
xlab("year") +
theme( strip.text.y = element_blank(),strip.text.x = element_blank())+
geom_text(data = meancvdf_redux_alpha, aes(x=-Inf,y=Inf,hjust=-0.05,
                vjust=1.7,label=labels), size=4.5)+
facet_grid(scenario2~model2, scales="free_y")
psmsy_alphascn_line


#line plot
df_smsy_est_redux_alpha<-df_smsy_est[
           df_smsy_est$scenario%in%c("autocorr","decLinearProd", "shiftProd", "sineProd")&
           df_smsy_est$model2%in%c("autocorr", "rw.a", "hmm.a")&
           df_smsy_est$method=="MLE",]


df_smsy_est_redux_alpha$scenario2<-case_match(df_smsy_est_redux_alpha$scenario,
     "autocorr"~ "stationary",
    "decLinearProd"~ "linear decline",
     "shiftProd" ~ "shift decline" , 
     "sineProd" ~ "sine fluctuation")


df_smsy_est_redux_alpha$scenario2<-factor(df_smsy_est_redux_alpha$scenario2, 
    levels=c("stationary",
     "linear decline",
     "sine fluctuation",
     "shift decline" ))


df_smsy_sim_redux_alpha<-df_smsy_sim[
           df_smsy_sim$scenario%in%c("autocorr","decLinearProd", "shiftProd", "sineProd")&
           df_smsy_sim$model2%in%c("autocorr", "rw.a", "hmm.a")&
           df_smsy_sim$method=="MLE",]


df_smsy_est_redux_alpha[df_smsy_est_redux_alpha$iteration==1&
           df_smsy_est_redux_alpha$model=="autocorr",]


psmsy_alphascn_hair<-ggplot(df_smsy_est_redux_alpha) + 
geom_line(aes(x=by-54,y=value/10000, color=model2, group=iteration), alpha=.1)+
guides(color = guide_legend(override.aes = list(linewidth = 2, alpha=1)))+
geom_line(data=summarydf_smsy_sim_redux_alpha,aes(x=by-54,y= x/10000),color="black", alpha=.8,linewidth=1.2)+
geom_line(data=summarydf_smsy_redux_alpha,aes(x=by-54,y= x.50./10000),color="darkred", alpha=.8,linewidth=1.2)+
scale_color_viridis_d(begin=.2, end=.8,option = "E") +
scale_fill_viridis_d(begin=.2, end=.8,option = "E") +
 facet_grid(scenario2~model2, scales="free_y")+
coord_cartesian(ylim = c(2,19))+ 
mytheme+ 
ylab(expression(paste(S[MSY]~("10,000s")))) +
xlab("year") +
mytheme+
theme( strip.text.y = element_blank(),strip.text.x = element_blank(),
    legend.text = element_text(size=18),legend.key.width = unit(2, 'cm'))+
geom_text(data = meancvdf_redux_alpha, aes(x=-Inf,y=Inf,hjust=-0.05,
                vjust=1.7,label=labels), size=4.5)
psmsy_alphascn_hair



df_smsy_est_redux_alpha<- df[df$parameter=="smsy"&df$variable=="mode"&
df$scenario%in%c("autocorr","decLinearProd", "shiftProd", "sineProd")&
df$model%in%c("autocorr", "rwa", "hmma")&
df$method=="MLE",]


df_smsy_est_redux_alpha$scenario2<-case_match(df_smsy_est_redux_alpha$scenario,
    "autocorr"~ "stationary",
    "decLinearProd"~ "linear decline",
     "shiftProd" ~ "shift decline" , 
     "sineProd" ~ "sine fluctuation")


df_smsy_est_redux_alpha$scenario2<-factor(df_smsy_est_redux_alpha$scenario2, levels=c("stationary",
     "linear decline",
     "sine fluctuation",
     "shift decline" ))
head(meanabspbiasdf_redux_alpha)
psmsy_alphascn_violin_abs<-ggplot(df_smsy_est_redux_alpha) + 
geom_violin(aes(x=model2,y=abs(bias), fill=model2), scale="width", trim=TRUE, alpha=.7,adjust = 1.8)+
geom_boxplot(aes(x=model2,y=abs(bias), fill=model2), outlier.shape = NA, width=0.1)+
geom_hline(yintercept=0,color="black", alpha=.6,linewidth=1.2)+
 scale_fill_viridis_d(begin=.2, end=.8,,option = "E") +
 facet_grid(scenario2~., scales="free_y")+
 ylab(expression(absolute~bias~"in"~ S[MSY]~("10,000s"))) +
  coord_cartesian(ylim = c(0,150000))+ 
 geom_text(data = meanabspbiasdf_redux_alpha, aes(x=model2,y=Inf,hjust=0,
                vjust=1.0,label=meanabspbias), size=6)+

mytheme
psmsy_alphascn_violin_abs

multi.page.abs.smsy_alphascn <- ggarrange(psmsy_alphascn_line, psmsy_alphascn_violin_abs,
                        nrow = 1, ncol = 2,
                        common.legend = TRUE,
                        legend="bottom")
multi.page.abs.smsy_alphascn

head(meanabspbiasdf_redux_alpha)

meanabspbiasdf_redux_alpha$labels<-paste("|%bias| = ",meanabspbiasdf_redux_alpha$meanabspbias)

psmsy_alphascn_boxplot_abs<-ggplot(df_smsy_est_redux_alpha) + 
#geom_violin(aes(x=model2,y=abs(bias), fill=model2), scale="width", trim=TRUE, alpha=.7,adjust = 1.8)+
geom_boxplot(aes(x=model2,y=abs(bias)/10000, fill=model2), outlier.shape = NA, width=0.8)+
geom_hline(yintercept=0,color="black", alpha=.6,linewidth=1.2)+
 scale_fill_viridis_d(begin=.2, end=.8,,option = "E") +
 facet_grid(scenario2~., scales="free_y")+
 ylab(expression(absolute~bias~"in"~ S[MSY]~("10,000s"))) +
  coord_cartesian(ylim = c(0,5))+ 
   xlab("estimation model") +
 geom_text(data = meanabspbiasdf_redux_alpha, aes(x=model2,y=Inf,hjust=0.5,
                vjust=1.7,label=labels), size=4.5)+
mytheme
psmsy_alphascn_boxplot_abs

multi.page.abs.smsy_alphascn.bx <- ggarrange(psmsy_alphascn_line, psmsy_alphascn_boxplot_abs,
                        nrow = 1, ncol = 2,
                        common.legend = TRUE,
                        legend="bottom")
multi.page.abs.smsy_alphascn.bx

multi.page.abs.smsy_alphascn_title<-annotate_figure(multi.page.abs.smsy_alphascn.bx, 
    top = text_grob(expression("Stationary and time-varying"~log(alpha)), 
               face = "bold" , size = 14))

ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/summary_figs/bias_trends_violins_smsy_basealpha_scn_cv.png",
    plot=multi.page.abs.smsy_alphascn_title,width = 14,height = 8, bg = "white")


#hair plot version


multi.page.abs.smsy_alphascn.hairbx <- ggarrange(psmsy_alphascn_hair, psmsy_alphascn_boxplot_abs,
                        nrow = 1, ncol = 2,
                        common.legend = TRUE,
                        legend="bottom")
multi.page.abs.smsy_alphascn.hairbx

multi.page.abs.smsy_alphascn_title.hairbx<-annotate_figure(multi.page.abs.smsy_alphascn.hairbx, 
    top = text_grob(expression("Stationary and time-varying"~log(alpha)), 
               face = "bold" , size = 14))

ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/summary_figs/bias_trends_smsy_basealpha_scn_hairbx.png",
    plot=multi.page.abs.smsy_alphascn_title.hairbx,width = 14,height = 8, bg = "white")



#==================================================================
# Bias when beta or both parameters change


meancvdf_redux_smaxboth<-meancvdf[meancvdf$parameter=="smsy"&
meancvdf$scenario%in%c("decLinearCap",  "shiftCap","regimeProdCap",
         "decLinearProdshiftCap")&
meancvdf$model%in%c("autocorr","rwb", "rwab", "hmmb", "hmmab")&
meancvdf$method=="MLE",]

meancvdf_redux_smaxboth$meancv<-round(meancvdf_redux_smaxboth$x,2)
meancvdf_redux_smaxboth$scenario2<-case_match(meancvdf_redux_smaxboth$scenario,
     "decLinearCap"~ "linear decline - cap", 
     "shiftCap" ~ "shift decline - cap",
     "regimeProdCap" ~ "shift - both",
    "decLinearProdshiftCap"~ "mixed trend - both")
meancvdf_redux_smaxboth$model2<-factor(meancvdf_redux_smaxboth$model2, levels=c("autocorr","rw.b", "rw.ab","hmm.b", "hmm.ab"))

meancvdf_redux_smaxboth$scenario2<-factor(meancvdf_redux_smaxboth$scenario2,levels=c("linear decline - cap", 
    "shift decline - cap",
     "shift - both",
    "mixed trend - both" ))


#estimates

summarydf_smsy_sim_redux_smaxboth<-summarydf_smsy_sim[summarydf_smsy_sim$scenario%in%c("decLinearCap", 
      "shiftCap","regimeProdCap","decLinearProdshiftCap")&
         summarydf_smsy_sim$model2%in%c("autocorr","rw.b", "rw.ab", "hmm.b", "hmm.ab")&
         summarydf_smsy_sim$method=="MLE",]

summarydf_smsy_sim_redux_smaxboth$scenario2<-case_match(summarydf_smsy_sim_redux_smaxboth$scenario,
    "decLinearCap"~ "linear decline - cap", 
     "shiftCap" ~ "shift decline - cap",
     "regimeProdCap" ~ "shift - both",
    "decLinearProdshiftCap"~ "mixed trend - both")


summarydf_smsy_sim_redux_smaxboth$scenario2<-factor(summarydf_smsy_sim_redux_smaxboth$scenario2, levels=c(
    "linear decline - cap", 
    "shift decline - cap",
     "shift - both",
    "mixed trend - both"))


summarydf_smsy_redux_smaxboth<-summarydf_smsy[summarydf_smsy$scenario%in%c( "decLinearCap",  "shiftCap","regimeProdCap",
         "decLinearProdshiftCap" )&
  summarydf_smsy$model%in%c("autocorr","rw.b", "rw.ab", "hmm.b", "hmm.ab") &summarydf_smsy$method=="MLE",]
#"autocorr"

summarydf_smsy_redux_smaxboth$scenario2<-case_match(summarydf_smsy_redux_smaxboth$scenario,
    "decLinearCap"~ "linear decline - cap", 
     "shiftCap" ~ "shift decline - cap",
     "regimeProdCap" ~ "shift - both",
    "decLinearProdshiftCap"~ "mixed trend - both")

summarydf_smsy_redux_smaxboth$scenario2<-factor(summarydf_smsy_redux_smaxboth$scenario2, levels=c(
    "linear decline - cap", 
    "shift decline - cap",
     "shift - both",
    "mixed trend - both"))


meancvdf_redux_smaxboth$labels<-paste("CV = ",meancvdf_redux_smaxboth$meancv)


psmsy_smaxscn_line<-ggplot() + 
geom_pointrange(data=summarydf_smsy_redux_smaxboth,aes(x=by-54,y= x.50./10000,ymin = x.2.5./10000,
 ymax = x.97.5./10000, color=model2),alpha=.9)+
geom_line(data=summarydf_smsy_sim_redux_smaxboth,aes(x=by-54,y= x/10000),color="black", alpha=.8,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8,option = "E") +
scale_fill_viridis_d(begin=.1, end=.8,option = "E") +
coord_cartesian(ylim = c(1,14))+ 
mytheme+ 
ylab(expression(paste(S[MSY]~("10,000s"))))  +
xlab("year") +
theme( strip.text.y = element_blank(),strip.text.x = element_blank())+
geom_text(data = meancvdf_redux_smaxboth, aes(x=-Inf,y=Inf,hjust=-0.05,
                vjust=1.7,label=labels), size=4.5)+
facet_grid(scenario2~model2, scales="free_y")
psmsy_smaxscn_line
head(meancvdf_redux_smaxboth)




#hair plot
df_smsy_est_redux_smaxboth<-df_smsy_est[df_smsy_est$scenario%in%c("decLinearCap", 
      "shiftCap","regimeProdCap","decLinearProdshiftCap")&
df_smsy_est$model2%in%c("autocorr","rw.b", "rw.ab", "hmm.b", "hmm.ab")&
df_smsy_est$method=="MLE",]


df_smsy_est_redux_smaxboth$scenario2<-case_match(df_smsy_est_redux_smaxboth$scenario,
     "decLinearCap"~ "linear decline - cap", 
     "shiftCap" ~ "shift decline - cap",
     "regimeProdCap" ~ "shift - both",
    "decLinearProdshiftCap"~ "mixed trend - both")


df_smsy_est_redux_smaxboth$scenario2<-factor(df_smsy_est_redux_smaxboth$scenario2, 
    levels=c("linear decline - cap", 
    "shift decline - cap",
     "shift - both",
    "mixed trend - both"))


df_smsy_sim_redux_smaxboth<-df_smsy_sim[df_smsy_sim$scenario%in%c("decLinearCap", 
      "shiftCap","regimeProdCap","decLinearProdshiftCap")&
df_smsy_sim$model2%in%c("autocorr","rw.b", "rw.ab", "hmm.b", "hmm.ab"))&
df_smsy_sim$method=="MLE",]

head(df_smsy_est_redux_alpha)
unique(df_smsy_est_redux_alpha$parameter)


psmsy_smaxbothscn_hair<-ggplot(df_smsy_est_redux_smaxboth) + 
geom_line(aes(x=by-54,y=value/10000, color=model2, group=iteration), alpha=.1)+
guides(color = guide_legend(override.aes = list(linewidth = 2, alpha=1)))+
geom_line(data=summarydf_smsy_sim_redux_smaxboth,aes(x=by-54,y= x/10000),color="black", alpha=.8,linewidth=1.2)+
geom_line(data=summarydf_smsy_redux_smaxboth,aes(x=by-54,y= x.50./10000),color="darkred", alpha=.8,linewidth=1.2)+
scale_color_viridis_d(begin=.2, end=.8,option = "E") +
scale_fill_viridis_d(begin=.2, end=.8,option = "E") +
 facet_grid(scenario2~model2, scales="free_y")+
coord_cartesian(ylim = c(2,15))+ 
mytheme+ 
ylab(expression(paste(S[MSY]~("10,000s")))) +
xlab("year") +
theme(strip.text.y = element_blank(),strip.text.x = element_blank(),
    legend.text = element_text(size=18),legend.key.width = unit(2, 'cm'))+
geom_text(data = meancvdf_redux_smaxboth, aes(x=-Inf,y=Inf,hjust=-0.05,
                vjust=1.7,label=labels), size=4.5)
psmsy_smaxbothscn_hair



#calculate pbias

meanabspbiasdf_redux_smaxboth<-meanabspbiasdf[meanabspbiasdf$parameter=="smsy"&
meanabspbiasdf$scenario%in%c("decLinearCap",  "shiftCap","regimeProdCap",
         "decLinearProdshiftCap")&
meanabspbiasdf$model%in%c("autocorr","rwb", "rwab","hmmb", "hmmab")&
meanabspbiasdf$method=="MLE",]


meanabspbiasdf_redux_smaxboth$meanabspbias<-round(meanabspbiasdf_redux_smaxboth$x,2)
meanabspbiasdf_redux_smaxboth$scenario2<-case_match(meanabspbiasdf_redux_smaxboth$scenario,
     "decLinearCap"~ "linear decline - cap", 
     "shiftCap" ~ "shift decline - cap",
     "regimeProdCap" ~ "shift - both",
    "decLinearProdshiftCap"~ "mixed trend - both")
unique(meanabspbiasdf_redux_smaxboth$model2)
meanabspbiasdf_redux_smaxboth$model2<-factor(meanabspbiasdf_redux_smaxboth$model2, 
    levels=c("autocorr","rw.b", "hmm.b", "rw.ab", "hmm.ab"))

meanabspbiasdf_redux_smaxboth$scenario2<-factor(meanabspbiasdf_redux_smaxboth$scenario2,
    levels=c("linear decline - cap", 
    "shift decline - cap",
     "shift - both",
    "mixed trend - both"))




df_smsy_est_redux_smaxboth<- df[df$parameter=="smsy"&df$variable=="mode"&
df$scenario%in%c("decLinearCap",  "shiftCap","regimeProdCap",
         "decLinearProdshiftCap")&
df$model%in%c("autocorr","rwb", "rwab","hmmb", "hmmab")&
df$method=="MLE",]
#"autocorr"

df_smsy_est_redux_smaxboth$scenario2<-case_match(df_smsy_est_redux_smaxboth$scenario,
     "decLinearCap"~ "linear decline - cap", 
     "shiftCap" ~ "shift decline - cap",
     "regimeProdCap" ~ "shift - both",
    "decLinearProdshiftCap"~ "mixed trend - both")

df_smsy_est_redux_smaxboth$scenario2<-factor(df_smsy_est_redux_smaxboth$scenario2, levels=c(
    "linear decline - cap", 
    "shift decline - cap",
     "shift - both",
    "mixed trend - both"))

head(df_smsy_est_redux_smaxboth)
psmsy_smaxscn_violin_abs<-ggplot(df_smsy_est_redux_smaxboth) + 
geom_violin(aes(x=model2,y=abs(bias), fill=model2), scale="width", trim=TRUE, alpha=.7,adjust = 1.8)+
geom_boxplot(aes(x=model2,y=abs(bias), fill=model2),outlier.shape = NA, width=0.1)+
geom_hline(yintercept=0,color="black", alpha=.6,linewidth=1.2)+
 scale_fill_viridis_d(begin=.1, end=.8,,option = "E") +
 facet_grid(scenario2~., scales="free_y")+
 ylab(expression(absolute~bias~"in" ~ S[MSY])) +
 coord_cartesian(ylim = c(0,100000))+ 
 xlab("estimation model")+
 mytheme
psmsy_smaxscn_violin_abs

multi.page.abs.smsy.smaxscn <- ggarrange(psmsy_smaxscn_line, psmsy_smaxscn_violin_abs,
                        nrow = 1, ncol = 2,
                        common.legend = TRUE,
                        legend="bottom")
multi.page.abs.smsy.smaxscn


multi.page.abs.smsy.smaxscn_title<-annotate_figure(multi.page.abs.smsy.smaxscn, 
    top = text_grob(expression("Time-varying"~S[max]~"or both parameters"), 
               face = "bold" , size = 14))



ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/summary_figs/bias_trends_violins_smsy_smaxbothscn.png",
    plot=multi.page.abs.smsy.smaxscn_title,width = 12,height = 8)





meanabspbiasdf_redux_smaxboth$labels<-paste("|%bias| =",meanabspbiasdf_redux_smaxboth$meanabspbias)

psmsy_smaxscn_boxplot_abs<-ggplot(df_smsy_est_redux_smaxboth) + 
geom_boxplot(aes(x=model2,y=abs(bias)/10000, fill=model2), outlier.shape = NA, width=0.8)+
geom_hline(yintercept=0,color="black", alpha=.6,linewidth=1.2)+
 scale_fill_viridis_d(begin=.2, end=.8,,option = "E") +
 facet_grid(scenario2~., scales="free_y")+
 ylab(expression(absolute~bias~"in"~ S[MSY]~("10,000s"))) +
  xlab("estimation model") +
  coord_cartesian(ylim = c(0,5))+ 
 geom_text(data = meanabspbiasdf_redux_smaxboth, aes(x=model2,y=Inf,hjust=.6,
                vjust=1.7,label=labels), size=4.5)+
mytheme
psmsy_smaxscn_boxplot_abs

multi.page.abs.smsy_smaxscn.bx <- ggarrange(psmsy_smaxscn_line, psmsy_smaxscn_boxplot_abs,
                        nrow = 1, ncol = 2,
                        common.legend = TRUE,
                        legend="bottom")
multi.page.abs.smsy_smaxscn.bx

multi.page.abs.smsy_smaxscn.bx_title<-annotate_figure(multi.page.abs.smsy_smaxscn.bx, 
    top = text_grob(expression("Time-varying"~S[max]~"or both parameters"), 
               face = "bold" , size = 14))
multi.page.abs.smsy_smaxscn.bx_title


ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/summary_figs/bias_trends_boxplot_smsy_smaxbothscn_scn_cv.png",
    plot=multi.page.abs.smsy_smaxscn.bx_title,width = 21,height = 8, bg = "white")

#hair plot version



multi.page.abs.smsy_smaxbothscn.hairbx <- ggarrange(psmsy_smaxbothscn_hair, psmsy_smaxscn_boxplot_abs,
                        nrow = 1, ncol = 2,
                        common.legend = TRUE,
                        legend="bottom")
multi.page.abs.smsy_smaxbothscn.hairbx

multi.page.abs.smsy_smaxbothscn_title.hairbx<-annotate_figure(multi.page.abs.smsy_smaxbothscn.hairbx, 
    top = text_grob(expression("Time-varying"~S[max]~"or both parameters"), 
               face = "bold" , size = 14))
multi.page.abs.smsy_smaxbothscn_title.hairbx

ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/summary_figs/bias_trends_boxplot_smsy_smaxbothscn_scn_hairbx.png",
    plot=multi.page.abs.smsy_smaxbothscn_title.hairbx,width = 22,height = 8, bg = "white")



##Old code??###
#==================================================================
# Bias in alpha and beta for hershey and quantile plots, 
#need to update it to hairplot and boxplots
#redux other parameter plots

# Bias in alpha
summarydf_alpha_sim_redux<-summarydf_alpha_sim[summarydf_alpha_sim$scenario%in%c("autocorr","decLinearProd", "regimeProd", "sineProd")&
summarydf_alpha$model%in%c("autocorr", "rwa", "hmma")&summarydf_alpha$method=="MLE",]
head(summarydf_alpha_sim_redux)

summarydf_alpha_sim_redux$scenario2<-case_match(summarydf_alpha_sim_redux$scenario,
    "autocorr"~ "stationary",
    "decLinearProd"~ "linear decline",
     "regimeProd" ~ "shift increase" , 
     "sineProd" ~ "sine fluctuation")

summarydf_alpha_sim_redux$scenario2<-factor(summarydf_alpha_sim_redux$scenario2,levels=c("stationary",
    "linear decline","sine fluctuation","shift increase" ))


summarydf_alpha_redux<-summarydf_alpha[summarydf_alpha$scenario%in%c("autocorr", "decLinearProd", 
    "regimeProd", "sineProd" )&
  summarydf_alpha$model%in%c("autocorr", "rwa", "hmma") &summarydf_alpha$method=="MLE",]


summarydf_alpha_redux$scenario2<-case_match(summarydf_alpha_redux$scenario,
     "autocorr"~ "stationary",
    "decLinearProd"~ "linear decline",
     "regimeProd" ~ "shift increase" , 
     "sineProd" ~ "sine fluctuation")

summarydf_alpha_redux$scenario2<-factor(summarydf_alpha_redux$scenario2,levels=c("stationary",
    "linear decline","sine fluctuation","shift increase" ))


palpha_line<-ggplot() + 
geom_pointrange(data=summarydf_alpha_redux,aes(x=by,y= median,ymin =lower, ymax = upper, color=model),alpha=.9)+
geom_line(data=summarydf_alpha_sim_redux,aes(x=by,y=median),color="black", alpha=.8,linewidth=1.2)+
scale_color_viridis_d(begin=.2, end=.8,option = "E") +
scale_fill_viridis_d(begin=.2, end=.8,option = "E") +
mytheme+ 
ylab(expression(log(alpha))) +
xlab("year") +
theme( strip.text.y = element_blank(),strip.text.x = element_blank())+
facet_grid(scenario2~model, scales="free_y")
palpha_line


df_alpha_est_redux<- df[df$parameter=="logalpha"&df$variable=="mode"&
df$scenario%in%c("autocorr","decLinearProd", "shiftProd", "sineProd")&
df$model%in%c("autocorr", "rwa", "hmma")&
df$method=="MLE",]


df_alpha_est_redux$scenario2<-case_match(df_alpha_est_redux$scenario,
    "autocorr"~ "stationary",
    "decLinearProd"~ "linear decline",
     "shiftProd" ~ "shift increase" , 
     "sineProd" ~ "sine fluctuation")

df_alpha_est_redux$scenario2<-factor(df_alpha_est_redux$scenario2,levels=c("stationary",
    "linear decline","sine fluctuation","shift increase" ))


head(df_alpha_est_redux)


palpha_violin_abs<-ggplot(df_alpha_est_redux) + 
geom_violin(aes(x=model,y=abs(bias), fill=model), scale="width", trim=TRUE, alpha=.7,adjust = 1.8)+
geom_boxplot(aes(x=model,y=abs(bias), fill=model),outlier.shape = NA, width=0.1, alpha=.85)+
geom_hline(yintercept=0,color="black", alpha=.6,linewidth=1.2)+
 scale_fill_viridis_d(begin=.2, end=.8,,option = "E") +
 facet_grid(scenario2~., scales="free_y")+
 ylab(expression(paste("absolute bias in ",~log(alpha)))) +
 xlab("estimation model") +
mytheme
palpha_violin_abs

multi.page.abs <- ggarrange(palpha_line, palpha_violin_abs,
                        nrow = 1, ncol = 2,
                        common.legend = TRUE,
                        legend="bottom")
multi.page.abs

ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/summary_figs/bias_trends_violins_log_alpha.png",
    plot=multi.page.abs,width = 12,height = 8)


#GUIDELINES
head(resparam)
#cv for line plot
#calculate mean absolute bias across time series 
absbiasdf<-aggregate(resparam$bias, list(parameter=resparam$parameter,
                               scenario=resparam$scenario,
                               method=resparam$method,
                               model=resparam$model,
                               by=resparam$by ), function(x){mean(abs(x),na.rm=T)})

meanabsbiasdf<-aggregate(absbiasdf$x, list(parameter=absbiasdf$parameter,
                               scenario=absbiasdf$scenario,
                               method=absbiasdf$method,
                               model=absbiasdf$model
                              ), function(x){mean(x,na.rm=T)})

head(meanabsbiasdf)
meanabsbiasdf[meanabsbiasdf$method=="MLE"&meanabsbiasdf$parameter=="smsy"&meanabsbiasdf$scenario=="decLinearProd",]

meanabsbiasdf[meanabsbiasdf$method=="MLE"&meanabsbiasdf$parameter=="smsy"&meanabsbiasdf$scenario=="stationary",]

meanabsbiasdf_redux<-meanabsbiasdf[meanabsbiasdf$parameter=="logalpha"&
meanabsbiasdf$scenario%in%c("autocorr","decLinearProd", "regimeProd", "sineProd")&
meanabsbiasdf$model%in%c("autocorr", "rwa", "hmma")&
meanabsbiasdf$method=="MLE",]

head(meanabsbiasdf_redux)
meanabsbiasdf_redux$meanabsbias<-round(meanabsbiasdf_redux$x,2)
meanabsbiasdf_redux$scenario2<-case_match(meanabsbiasdf_redux$scenario,
     "autocorr"~ "stationary",
    "decLinearProd"~ "linear decline",
     "regimeProd" ~ "shift increase" , 
     "sineProd" ~ "sine fluctuation")
meanabsbiasdf_redux$model<-factor(meanabsbiasdf_redux$model, levels=c("autocorr", "rwa", "hmma"))

meanabsbiasdf_redux$scenario2<-factor(meanabsbiasdf_redux$scenario2,levels=c("stationary",
    "linear decline","sine fluctuation","shift increase" ))




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
meancvdf$scenario%in%c("autocorr","decLinearProd", "regimeProd", "sineProd")&
meancvdf$model%in%c("autocorr", "rwa", "hmma")&
meancvdf$method=="MLE",]

meancvdf_redux$meancv<-round(meancvdf_redux$x,2)
meancvdf_redux$scenario2<-case_match(meancvdf_redux$scenario,
     "autocorr"~ "stationary",
    "decLinearProd"~ "linear decline",
     "regimeProd" ~ "shift increase" , 
     "sineProd" ~ "sine fluctuation")
meancvdf_redux$model<-factor(meancvdf_redux$model, levels=c("autocorr", "rwa", "hmma"))

meancvdf_redux$scenario2<-factor(meancvdf_redux$scenario2,levels=c("stationary",
    "linear decline","sine fluctuation","shift increase" ))


palpha_line_cv<-ggplot() + 
geom_pointrange(data=summarydf_alpha_redux,aes(x=by,y= median,ymin = lower, ymax = upper, color=model),alpha=.9)+
geom_line(data=summarydf_alpha_sim_redux,aes(x=by,y= median),color="black", alpha=.8,linewidth=1.2)+
scale_color_viridis_d(begin=.2, end=.8,option = "E") +
scale_fill_viridis_d(begin=.2, end=.8,option = "E") +
#coord_cartesian(ylim = c(0.2,3.0))+ 
mytheme+ 
ylab(expression(log(alpha))) +
xlab("year") +
theme( strip.text.y = element_blank(),strip.text.x = element_blank())+
geom_text(data = meancvdf_redux, aes(x=-Inf,y=Inf,hjust=0,
                vjust=1.0,label=meancv), size=6)+
facet_grid(scenario2~model, scales="free_y")
palpha_line_cv


#pbias plots

head(df_alpha_est_redux)

df_alpha_est_redux<- df[df$parameter=="logalpha"&df$variable=="mode"&
df$scenario%in%c("autocorr","decLinearProd", "regimeProd", "sineProd")&
df$model%in%c("autocorr", "rwa", "hmma")&
df$method=="MLE",]

df_alpha_est_redux$scenario2<-case_match(df_alpha_est_redux$scenario,
    "autocorr"~ "stationary",
    "decLinearProd"~ "linear decline",
     "regimeProd" ~ "shift increase" , 
     "sineProd" ~ "sine fluctuation")

df_alpha_est_redux$scenario2<-factor(df_alpha_est_redux$scenario2,levels=c("stationary",
    "linear decline","sine fluctuation","shift increase" ))

unique(df_alpha_est_redux$scenario2)
head(df_alpha_est_redux)

palpha_violin_abspbias<-ggplot(df_alpha_est_redux) + 
geom_violin(aes(x=model,y=abs(pbias), fill=model), scale="width", trim=TRUE, alpha=.7,adjust = 1.8)+
geom_boxplot(aes(x=model,y=abs(pbias), fill=model),outlier.shape = NA, width=0.1)+
geom_hline(yintercept=0,color="black", alpha=.6,linewidth=1.2)+
 scale_fill_viridis_d(begin=.2, end=.8,,option = "E") +
 facet_grid(scenario2~., scales="free_y")+
 ylab(expression(paste("absolute bias in ",~log(alpha)))) +
 xlab("estimation model")+
mytheme
palpha_violin_abspbias


multi.page.abs.cv <- ggarrange(palpha_line_cv, palpha_violin_abs,
                        nrow = 1, ncol = 2,
                        common.legend = TRUE,
                        legend="bottom")
multi.page.abs.cv


ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/summary_figs/bias_trends_violins_cv_log_alpha.png",
    plot=multi.page.abs.cv,width = 12,height = 8)


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

summarydf_smax_sim_redux<-summarydf_smax_sim[summarydf_smax_sim$scenario%in%c("decLinearCap", "regimeCap", "shiftCap")&
summarydf_smax$model%in%c("autocorr", "rwb", "hmmb")&summarydf_alpha$method=="MLE",]


summarydf_smax_sim_redux$scenario2<-case_match(summarydf_smax_sim_redux$scenario,
    "decLinearCap"~ "linear decline",
     "regimeCap" ~ "regime increase" , 
     "shiftCap" ~ "shift decline")

summarydf_smax_sim_redux$scenario2<-factor(summarydf_smax_sim_redux$scenario2,levels=c(
    "linear decline","regime increase","shift decline" ))


summarydf_smax_redux<-summarydf_smax[summarydf_smax$scenario%in%c( "decLinearCap", "regimeCap", "shiftCap" )&
  summarydf_smax$model%in%c("autocorr", "rwb", "hmmb") &summarydf_smax$method=="MAP",]


summarydf_smax_redux$scenario2<-case_match(summarydf_smax_redux$scenario,
    "decLinearCap"~ "linear decline",
     "regimeCap" ~ "regime increase" , 
     "shiftCap" ~ "shift decline")


summarydf_smax_redux$scenario2<-factor(summarydf_smax_redux$scenario2,levels=c("stationary",
    "linear decline","regime increase","shift decline"  ))



meancvdf_redux_smax<-meancvdf[meancvdf$parameter=="smax"&
meancvdf$scenario%in%c("decLinearCap", "regimeCap", "shiftCap")&
meancvdf$model%in%c("autocorr", "rwb", "hmmb")&
meancvdf$method=="MLE",]

meancvdf_redux_smax$meancv<-round(meancvdf_redux_smax$x,2)
meancvdf_redux_smax$scenario2<-case_match(meancvdf_redux_smax$scenario,
     "decLinearCap"~ "linear decline",
     "regimeCap" ~ "regime increase" , 
     "shiftCap" ~ "shift decline")
meancvdf_redux_smax$model<-factor(meancvdf_redux_smax$model, levels=c("autocorr", "rwb", "hmmb"))

meancvdf_redux_smax$scenario2<-factor(meancvdf_redux_smax$scenario2,levels=c("stationary",
    "linear decline","regime increase","shift decline"))



psmax_line<-ggplot() + 
geom_pointrange(data=summarydf_smax_redux,aes(x=by,y= x.50.,ymin = x.2.5., ymax = x.97.5., color=model),alpha=.9)+
geom_line(data=summarydf_smax_sim_redux,aes(x=by,y= x),color="black", alpha=.8,linewidth=1.2)+
scale_color_viridis_d(begin=.2, end=.8,option = "E") +
scale_fill_viridis_d(begin=.2, end=.8,option = "E") +
coord_cartesian(ylim = c(50000,400000))+ 
mytheme+ 
ylab(expression(S[max])) +
xlab("year") +
theme( strip.text.y = element_blank(),strip.text.x = element_blank())+
geom_text(data = meancvdf_redux_smax, aes(x=-Inf,y=Inf,hjust=0,
                vjust=1.0,label=meancv), size=6)+
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



psmax_violin_abs<-ggplot(df_smax_est_redux) + 
geom_violin(aes(x=model,y=abs(bias), fill=model), scale="width", trim=TRUE, alpha=.7,adjust = 1.8)+
geom_boxplot(aes(x=model,y=abs(bias), fill=model),outlier.shape = NA, width=0.1)+
geom_hline(yintercept=0,color="black", alpha=.6,linewidth=1.2)+
 scale_fill_viridis_d(begin=.2, end=.8,,option = "E") +
 facet_grid(scenario2~., scales="free_y")+
 ylab("absolute bias in smax") +
 coord_cartesian(ylim = c(0,450000))+
 xlab("estimation model") +
mytheme
psmax_violin_abs

multi.page.abs.smax <- ggarrange(psmax_line, psmax_violin_abs,
                        nrow = 1, ncol = 2,
                        common.legend = TRUE,
                        legend="bottom")
multi.page.abs.smax

ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/summary_figs/bias_trends_violins_cv_smax.png",
    plot=multi.page.abs.smax,width = 12,height = 8)



head(df_smax_sim)
head(df_smax_est_redux)

psmax_trends_abs<-ggplot(df_smax_est_redux) + 
geom_line(aes(x=by,y=value, color=model, group=iteration), alpha=.1)+
geom_line(data=summarydf_smax_sim_redux,aes(x=by,y= x),color="black", alpha=.8,linewidth=1.2)+
 scale_color_viridis_d(begin=.1, end=.8) +
 facet_grid(scenario2~model, scales="free_y")+
 ylab("absolute bias in smax") +
 ggtitle("estimation model")+
 coord_cartesian(ylim = c(0,750000))+ 
mytheme
psmax_trends_abs


#==================================================================
# Bias in alpha and beta



meancvdf_redux_both<-meancvdf[meancvdf$parameter%in%c("alpha","smax")&
meancvdf$scenario%in%c("regimeProdCap",
         "decLinearProdshiftCap")&
meancvdf$model%in%c("autocorr", "rwab", "hmmab")&
meancvdf$method=="MLE",]

meancvdf_redux_both$meancv<-round(meancvdf_redux_both$x,2)
meancvdf_redux_both$scenario2<-case_match(meancvdf_redux_both$scenario,
      "regimeProdCap"~ "regime",
     "decLinearProdshiftCap" ~ "mixed trends")
meancvdf_redux_both$model<-factor(meancvdf_redux_both$model, levels=c("autocorr", "rwab", "hmmab"))

#meancvdf_redux_both$scenario2<-factor(meancvdf_redux_both$scenario2,levels=c("regimeProdCap",
#         "decLinearProdshiftCap"))




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



summarydf_alpha.smax_sim_redux<-summarydf_alpha.smax_sim[summarydf_alpha.smax_sim$scenario%in%c("regimeProdCap",
         "decLinearProdshiftCap")&
summarydf_alpha.smax_sim$model%in%c("autocorr", "rwab", "hmmab")&summarydf_alpha.smax_sim$method=="MAP",]

summarydf_alpha.smax_sim_redux$scenario2<-case_match(summarydf_alpha.smax_sim_redux$scenario,
    "regimeProdCap"~ "regime",
     "decLinearProdshiftCap" ~ "mixed trends" )

summarydf_alpha.smax_redux<-summarydf_alpha.smax[summarydf_alpha.smax$scenario%in%c( "regimeProdCap",
         "decLinearProdshiftCap" )&
  summarydf_alpha.smax$model%in%c("autocorr", "rwab", "hmmab") &summarydf_alpha.smax$method=="MAP",]


summarydf_alpha.smax_redux$scenario2<-case_match(summarydf_alpha.smax_redux$scenario,
      "regimeProdCap"~ "regime",
     "decLinearProdshiftCap" ~ "mixed trends")

head(summarydf_alpha.smax_redux)

palpha.smax_line<-ggplot() + 
geom_pointrange(data=summarydf_alpha.smax_redux,aes(x=by-54,y= x.50.,ymin = x.2.5., ymax = x.97.5., color=model),alpha=.9)+
geom_line(data=summarydf_alpha.smax_sim_redux,aes(x=by-54,y= x),color="black", alpha=.8,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8,option = "E") +
scale_fill_viridis_d(begin=.1, end=.8,option = "E") +
mytheme+ 
ylab("parameter value") +
xlab("year") +
facet_grid(scenario2+parameter~model, scales="free_y")+
geom_text(data = meancvdf_redux_both, aes(x=-Inf,y=Inf,hjust=0,
                vjust=1.0,label=meancv), size=6)+
theme( strip.text.y = element_blank(),strip.text.x = element_blank())
palpha.smax_line






df_alpha.smax_est_redux<- df[df$parameter%in%c("alpha","smax")&df$variable=="mode"&
df$scenario%in%c("regimeProdCap",
         "decLinearProdshiftCap")&
df$model%in%c("autocorr", "rwab", "hmmab")&
df$method=="MAP",]
head(df_alpha.smax_est_redux)

df_alpha.smax_est_redux$scenario2<-case_match(df_alpha.smax_est_redux$scenario,
    "regimeProdCap"~ "regime",
     "decLinearProdshiftCap" ~ "mixed trends")



palpha.smax_violin_abs<-ggplot(df_alpha.smax_est_redux) + 
geom_violin(aes(x=model,y=abs(bias), fill=model), scale="width", trim=TRUE, alpha=.7,adjust = 1.8)+
geom_boxplot(aes(x=model,y=abs(bias), fill=model),outlier.shape = NA, width=0.1)+
geom_hline(yintercept=0,color="black", alpha=.6,linewidth=1.2)+
 scale_fill_viridis_d(begin=.1, end=.8,,option = "E") +
 facet_grid(scenario2+parameter~., scales="free_y")+
 ylab("absolute bias in parameter") + 
 xlab("estimation model")+
mytheme

palpha.smax_violin_abs

multi.page.abs.alpha.smax <- ggarrange(palpha.smax_line, palpha.smax_violin_abs,
                        nrow = 1, ncol = 2,
                        common.legend = TRUE,
                        legend="bottom")
multi.page.abs.alpha.smax

ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/summary_figs/bias_trends_violins_cv_both.png",
    plot=multi.page.abs.alpha.smax,width = 12,height = 8)



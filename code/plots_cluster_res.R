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

#res<-readRDS(file = "outs/simest/generic/res.rds")
res1<-readRDS(file = "outs/simest/generic/res1.rds")
res2<-readRDS(file = "outs/simest/generic/res2.rds")
restmb<-rbind(res1,res2)
#
resstan1<-readRDS(file = "outs/simest/generic/resstan1.rds")
resstan2<-readRDS(file = "outs/simest/generic/resstan2.rds")
resstan<-rbind(resstan1,resstan2)
#resstan<-readRDS(file = "outs/simest/generic/resstan.rds")

res<-rbind(restmb,resstan)

res$parameter[res$parameter=="smax"]<-"Smax"

resparam<-res[res$parameter%in%c("alpha","Smax","sigma","smsy","sgen","umsy"),]

names(resparam)

#mean by iteration
dfmpbias<-aggregate(resparam$pbias, 
    list(model=resparam$model,
         method=resparam$method,
         iteration=resparam$iteration,
         parameter=resparam$parameter,
         scenario=resparam$scenario),
    mean)

unique(dfmpbias$method)
dfmpbias$model <- factor(dfmpbias$model, 
    levels=c("simple",
             "autocorr",
             "rwa",
             "rwb",
             "rwab",
             "hmma",
             "hmmb",
             "hmmab"))
#=================================
#plots
#dfmpbias$method <- factor(dfmpbias$method, levels=c("MLE","MCMC"))
dfmpbias <- dfmpbias[!is.na(dfmpbias$x),] 

dfmpbias$simulated <- dplyr::recode(dfmpbias$scenario, 
      "stationary"="simple",
      "autocorr"="autocorr",
      "sigmaShift"="simple", 
      "decLinearProd"="rwa",
      "sineProd"="rwa",
      "regimeProd"="hmma",
      "decLinearCap"="rwb",
      "regimeCap"="hmmb",
      "shiftCap"="hmmb", 
      "shiftProd"="hmma",
      "regimeProdCap"="hmmab",
      "decLinearProdshiftCap"="rwab"
      )   

dfmpbias$simulated_f<-factor(dfmpbias$simulated, levels=c("simple",
                                            "autocorr",
                                            "rwa",
                                            "rwb",
                                            "rwab", 
                                            "hmma",
                                            "hmmb",
                                            "hmmab"))



dfmpbias$model <- factor(dfmpbias$model, levels=c(
       "simple", 
       "autocorr",
       "rwa",
       "rwb",
       "rwab",
       "hmma",  
        "hmmb",    
        "hmmab") )

     
dfmpbias$scenario_f<- factor(dfmpbias$scenario, levels=c(
      "stationary",
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
      "decLinearProdshiftCap"
      )  )


dfmpbias$model_agg <- dplyr::recode(dfmpbias$model, 
       "simple"="simple", 
       "autocorr"="simple",
       "rwa"="rw",
       "rwb"="rw",
       "rwab"="rw",
       "hmma"="hmm",  
       "hmmb"="hmm",    
       "hmmab"="hmm")

dfmpbias$simulated_agg <- dplyr::recode(dfmpbias$simulated, 
       "simple"="simple", 
       "autocorr"="simple",
       "rwa"="rw",
       "rwb"="rw",
       "rwab"="rw",
       "hmma"="hmm",  
       "hmmb"="hmm",    
       "hmmab"="hmm")

dfmpbias$model_paragg <- dplyr::recode(dfmpbias$model, 
       "simple"="simple", 
       "autocorr"="simple",
       "rwa"="tva",
       "rwb"="tvb",
       "rwab"="tvab",
       "hmma"="tva",  
       "hmmb"="tvb",    
       "hmmab"="tvab")

dfmpbias$simulated_paragg <- dplyr::recode(dfmpbias$simulated, 
       "simple"="simple", 
       "autocorr"="simple",
       "rwa"="tva",
       "rwb"="tvb",
       "rwab"="tvab",
       "hmma"="tva",  
       "hmmb"="tvb",    
       "hmmab"="tvab")

dfscen<-aggregate(dfmpbias$simulated,list(parameter=dfmpbias$parameter,
  scenario_f=dfmpbias$scenario_f, model=dfmpbias$model),unique)
dfscen$x<-factor(dfscen$x, levels=c("simple", 
       "autocorr",
       "rwa",
       "rwb",
       "rwab",
       "hmma",  
        "hmmb",    
        "hmmab"))
dfscen$simulated_f<-as.numeric(dfscen$x)

dfmpbias_main<-dfmpbias[dfmpbias$parameter%in%c("alpha",
      "Smax",
      "sigma"),]

dfscen_main<-dfscen[dfscen$parameter%in%c("alpha",
      "Smax",
      "sigma"),]

boxplot_allscn(dfmpbias_main,dfscen_main)


#========================
#last 5 years

unique(resparam$by)
resparaml5<-resparam[resparam$by%in%c(90,91,92,93,94),]
dfmpbiasl5<-aggregate(resparaml5$pbias, 
    list(model=resparaml5$model,
         method=resparaml5$method,
         iteration=resparaml5$iteration,
         parameter=resparaml5$parameter,
         scenario=resparaml5$scenario),
    mean)


dfmpbiasl5$model <- factor(dfmpbiasl5$model, 
    levels=c("simple",
             "autocorr",
             "rwa",
             "rwb",
             "rwab",
             "hmma",
             "hmmb",
             "hmmab"))

dfmpbiasl5 <- dfmpbiasl5[!is.na(dfmpbiasl5$x),] 

dfmpbiasl5$simulated <- dplyr::recode(dfmpbiasl5$scenario, 
      "stationary"="simple",
      "autocorr"="autocorr",
      "sigmaShift"="simple", 
      "decLinearProd"="rwa",
      "sineProd"="rwa",
      "regimeProd"="hmma",
      "decLinearCap"="rwb",
      "regimeCap"="hmmb",
      "shiftCap"="hmmb", 
      "shiftProd"="hmma",
      "regimeProdCap"="hmmab",
      "decLinearProdshiftCap"="rwab"
      )   

dfmpbiasl5$simulated_f<-factor(dfmpbiasl5$simulated, levels=c("simple",
                                            "autocorr",
                                            "rwa",
                                            "rwb",
                                            "rwab", 
                                            "hmma",
                                            "hmmb",
                                            "hmmab"))



dfmpbiasl5$model <- factor(dfmpbiasl5$model, levels=c(
       "simple", 
       "autocorr",
       "rwa",
       "rwb",
       "rwab",
       "hmma",  
        "hmmb",    
        "hmmab") )

     
dfmpbiasl5$scenario_f<- factor(dfmpbiasl5$scenario, levels=c(
      "stationary",
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
      "decLinearProdshiftCap"
      )  )


dfmpbiasl5$model_agg <- dplyr::recode(dfmpbiasl5$model, 
       "simple"="simple", 
       "autocorr"="simple",
       "rwa"="rw",
       "rwb"="rw",
       "rwab"="rw",
       "hmma"="hmm",  
       "hmmb"="hmm",    
       "hmmab"="hmm")

dfmpbiasl5$simulated_agg <- dplyr::recode(dfmpbiasl5$simulated, 
       "simple"="simple", 
       "autocorr"="simple",
       "rwa"="rw",
       "rwb"="rw",
       "rwab"="rw",
       "hmma"="hmm",  
       "hmmb"="hmm",    
       "hmmab"="hmm")

dfmpbiasl5$model_paragg <- dplyr::recode(dfmpbiasl5$model, 
       "simple"="simple", 
       "autocorr"="simple",
       "rwa"="tva",
       "rwb"="tvb",
       "rwab"="tvab",
       "hmma"="tva",  
       "hmmb"="tvb",    
       "hmmab"="tvab")

dfmpbiasl5$simulated_paragg <- dplyr::recode(dfmpbiasl5$simulated, 
       "simple"="simple", 
       "autocorr"="simple",
       "rwa"="tva",
       "rwb"="tvb",
       "rwab"="tvab",
       "hmma"="tva",  
       "hmmb"="tvb",    
       "hmmab"="tvab")

dfscenl5<-aggregate(dfmpbiasl5$simulated,list(parameter=dfmpbiasl5$parameter,
  scenario_f=dfmpbiasl5$scenario_f, model=dfmpbiasl5$model),unique)
dfscenl5$x<-factor(dfscenl5$x, levels=c("simple", 
       "autocorr",
       "rwa",
       "rwb",
       "rwab",
       "hmma",  
        "hmmb",    
        "hmmab"))
dfscenl5$simulated_f<-as.numeric(dfscenl5$x)

dfmpbiasl5_main<-dfmpbiasl5[dfmpbiasl5$parameter%in%c("alpha",
      "Smax",
      "sigma"),]

dfscen_mainl5<-dfscen[dfscenl5$parameter%in%c("alpha",
      "Smax",
      "sigma"),]

boxplot_allscn(dfmpbiasl5_main,dfscen_main)



#=================================
#plots
#dfmpbias$method <- factor(dfmpbias$method, levels=c("MLE","MCMC"))
dfmpbias <- dfmpbias[!is.na(dfmpbias$x),] 

dfmpbias$simulated <- dplyr::recode(dfmpbias$scenario, 
      "stationary"="simple",
      "autocorr"="autocorr",
      "sigmaShift"="simple", 
      "decLinearProd"="rwa",
      "sineProd"="rwa",
      "regimeProd"="hmma",
      "decLinearCap"="rwb",
      "regimeCap"="hmmb",
      "shiftCap"="hmmb", 
      "shiftProd"="hmma",
      "regimeProdCap"="hmmab",
      "decLinearProdshiftCap"="rwab"
      )   

dfmpbias$simulated_f<-factor(dfmpbias$simulated, levels=c("simple",
                                            "autocorr",
                                            "rwa",
                                            "rwb",
                                            "rwab", 
                                            "hmma",
                                            "hmmb",
                                            "hmmab"))



dfmpbias$model <- factor(dfmpbias$model, levels=c(
       "simple", 
       "autocorr",
       "rwa",
       "rwb",
       "rwab",
       "hmma",  
        "hmmb",    
        "hmmab") )

     
dfmpbias$scenario_f<- factor(dfmpbias$scenario, levels=c(
      "stationary",
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
      "decLinearProdshiftCap"
      )  )


dfmpbias$model_agg <- dplyr::recode(dfmpbias$model, 
       "simple"="simple", 
       "autocorr"="simple",
       "rwa"="rw",
       "rwb"="rw",
       "rwab"="rw",
       "hmma"="hmm",  
       "hmmb"="hmm",    
       "hmmab"="hmm")

dfmpbias$simulated_agg <- dplyr::recode(dfmpbias$simulated, 
       "simple"="simple", 
       "autocorr"="simple",
       "rwa"="rw",
       "rwb"="rw",
       "rwab"="rw",
       "hmma"="hmm",  
       "hmmb"="hmm",    
       "hmmab"="hmm")

dfmpbias$model_paragg <- dplyr::recode(dfmpbias$model, 
       "simple"="simple", 
       "autocorr"="simple",
       "rwa"="tva",
       "rwb"="tvb",
       "rwab"="tvab",
       "hmma"="tva",  
       "hmmb"="tvb",    
       "hmmab"="tvab")

dfmpbias$simulated_paragg <- dplyr::recode(dfmpbias$simulated, 
       "simple"="simple", 
       "autocorr"="simple",
       "rwa"="tva",
       "rwb"="tvb",
       "rwab"="tvab",
       "hmma"="tva",  
       "hmmb"="tvb",    
       "hmmab"="tvab")

dfscen<-aggregate(dfmpbias$simulated,list(parameter=dfmpbias$parameter,
  scenario_f=dfmpbias$scenario_f, model=dfmpbias$model),unique)
dfscen$x<-factor(dfscen$x, levels=c("simple", 
       "autocorr",
       "rwa",
       "rwb",
       "rwab",
       "hmma",  
        "hmmb",    
        "hmmab"))
dfscen$simulated_f<-as.numeric(dfscen$x)

dfmpbias_main<-dfmpbias[dfmpbias$parameter%in%c("alpha",
      "Smax",
      "sigma"),]

dfscen_main<-dfscen[dfscen$parameter%in%c("alpha",
      "Smax",
      "sigma"),]

boxplot_allscn(dfmpbias_main,dfscen_main)

head(dfmpbias_main)
dfmpbias_alpha<-dfmpbias_main[dfmpbias_main$parameter=="alpha"&dfmpbias_main$method=="MCMC",]
 

#========================




dfmpbias_deriv<-dfmpbias[dfmpbias$parameter%in%c("sgen",
      "smsy",
      "umsy"),]


dfscen_deriv<-dfscen[dfscen$parameter%in%c("sgen",
      "smsy",
      "umsy"),]

boxplot_allscn(dfmpbias_deriv,dfscen_deriv)

#aggregate plots
head(dfmpbias_deriv)

dfscen_agg<-aggregate(dfmpbias$simulated_agg,list(parameter=dfmpbias$parameter,
  scenario_f=dfmpbias$scenario_f, model=dfmpbias$model_agg),unique)
dfscen_agg$x<-factor(dfscen_agg$x, levels=c("simple", 
       "rw",
       "hmm"))
dfscen_agg$simulated_f<-as.numeric(dfscen_agg$x)

dfscen_main_agg<-dfscen_agg[dfscen_agg$parameter%in%c("alpha",
      "Smax",
      "sigma"),]

boxplot_agg(dfmpbias_main,dfscen_main_agg)



dfscen_paragg<-aggregate(dfmpbias$simulated_paragg,list(parameter=dfmpbias$parameter,
  scenario_f=dfmpbias$scenario_f, model=dfmpbias$model_agg),unique)
dfscen_paragg$x<-factor(dfscen_paragg$x, levels=c("simple", 
       "tva",
       "tvb",
       "tvab"))
dfscen_paragg$simulated_f<-as.numeric(dfscen_paragg$x)

dfscen_main_paragg<-dfscen_paragg[dfscen_paragg$parameter%in%c("alpha",
      "Smax",
      "sigma"),]

boxplot_paragg(dfmpbias_main,dfscen_main_paragg)

ggplot(dfmpbias_main) +   
  geom_boxplot(data=dfmpbias_main,aes(fill=method,x=model_paragg,y=x), outlier.shape = NA,alpha=.8) +         
  coord_cartesian(ylim = quantile(dfmpbias_main$x, c(0.025, 0.975)))+
  geom_hline(yintercept=0) +
  mytheme+ 
  ylab("mean % bias") +
  xlab("estimation model") +
  facet_grid(parameter~scenario_f, scales="free_y")+
  scale_fill_viridis_d(begin=.3, end=.8,option = "D") +
  scale_color_viridis_d(begin=.4, end=.8) +
  geom_rect(data=dfscen_main_paragg,aes(xmin=as.numeric(simulated_f)-.5,
                         xmax=as.numeric(simulated_f)+.5,
                         ymin=-Inf,ymax=Inf),
                    color="gray90",alpha=0.05)+
  #stat_summary(aes(color=method,x=(model),y=x),fun.data = give.n, geom = "text", 
  #        vjust = -2,position = position_dodge(1))+ 
  theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1))


#close look plot


dfmpbias_ss<-dfmpbias[dfmpbias$scenario_f%in%c("stationary",
      "decLinearProd",
      "regimeProd",
       "decLinearCap",
      "shiftCap"),]


dfmpbias_ss<-dfmpbias_ss[dfmpbias_ss$parameter%in%c("alpha",
      "Smax",
      "smsy"),]


dfscen_ss<-aggregate(dfmpbias_ss$simulated_agg,list(parameter=dfmpbias_ss$parameter,
  scenario_f=dfmpbias_ss$scenario_f, model=dfmpbias_ss$model_agg),unique)
dfscen_ss$x<-factor(dfscen_ss$x, levels=c("simple","rw","hmm"))
dfscen_ss$simulated_f<-as.numeric(dfscen_ss$x)






ggplot(dfmpbias_ss) +   
  geom_boxplot(data=dfmpbias_ss,aes(fill=method,x=model_agg,y=x), outlier.shape = NA) +         
  coord_cartesian(ylim = quantile(dfmpbias_ss$x, c(0.025, 0.975)))+
  geom_hline(yintercept=0) +
  mytheme+ 
  ylab("mean % bias") +
  xlab("estimation model") +
  facet_grid(parameter~scenario_f, scales="free_y")+
  scale_fill_viridis_d(begin=.3, end=.8,option = "D") +
  scale_color_viridis_d(begin=.3, end=.8) +
  geom_rect(data=dfscen_ss,aes(xmin=as.numeric(simulated_f)-.5,
                         xmax=as.numeric(simulated_f)+.5,
                         ymin=-Inf,ymax=Inf),
                    color="gray90",alpha=0.05)+
  #stat_summary(aes(color=method,x=(model),y=x),fun.data = give.n, geom = "text", 
  #        vjust = -2,position = position_dodge(1))+ 
  theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1))




#========================================================================================================
#sensitivity a scenario
#read in data
simPar <- read.csv("data/sensitivity/SimPars.csv")

## Store relevant object names to help run simulation 
scenNames <- unique(simPar$scenario)

res_a<-readRDS(file = "outs/simest/sensitivity/res_a.rds")
res_a_stan<-readRDS(file = "outs/simest/sensitivity/resstan_a.rds")


res_a<-rbind(res_a,res_a_stan)
res_a$parameter[res_a$parameter=="smax"]<-"Smax"

resparam<-res_a[res_a$parameter%in%c("alpha","Smax","sigma","smsy","sgen","umsy"),]
head(resparam)
resparam$bias<-(resparam$pbias/100)*resparam$sim
df<-reshape2::melt(resparam, id.vars=c("parameter","iteration","scenario","method","model","by", "convergence","pbias","bias"))
head(df)
df_alpha<-df[df$parameter%in%c("alpha"),]
df_alpha$col<-factor(df_alpha$variable,levels=c("est","sim"))

df_alpha_mle<-df_alpha[df_alpha$method=="MLE",]

param_traj(df_alpha_mle)
 
df_alpha_mcmc<-df_alpha[df_alpha$method=="MCMC",]

param_traj(df_alpha_mcmc)

head(df_alpha_mcmc)

df_alpha_pbias<-df_alpha_mcmc[df_alpha_mcmc$variable%in%c("est"),]

df_alpha_pbias$model<-factor(df_alpha_pbias$model, levels=c("simple", 
       "autocorr",
       "rwa",
       "rwb",
       "rwab",
       "hmma",  
        "hmmb",    
        "hmmab"))

unique(df_alpha_pbias$scenario)
head(df_alpha_pbias)

df_alpha_pbias1<-df_alpha_pbias[df_alpha_pbias$scenario%in%c("trendLinearProd1",  "trendLinearProd2", 
 "trendLinearProd5", "trendLinearProd7",  "trendLinearProd10"),]

ggplot(df_alpha_pbias1) +  
coord_cartesian(ylim = c(-100,100))+ 
       geom_boxplot(aes(x=as.factor(by),y=pbias ),
        outlier.shape = NA, alpha=.9) +  
       geom_hline(yintercept=0, color="darkred",alpha=.7,linewidth=1.2) +
       mytheme+ 
       ylab("parameters") +
       xlab("year") +
       facet_grid(scenario~model, scales="free_y")

ggplot(df_alpha_pbias1) +  
coord_cartesian(ylim = c(-1,1))+ 
       geom_boxplot(aes(x=as.factor(by),y=bias ),
        outlier.shape = NA, alpha=.9) +  
       geom_hline(yintercept=0, color="darkred",alpha=.7,linewidth=1.2) +
       mytheme+ 
       ylab("parameters") +
       xlab("year") +
       facet_grid(scenario~model, scales="free_y")



df_alpha_pbias2<-df_alpha_pbias[df_alpha_pbias$scenario%in%c("regimeProd1",      
 "regimeProd2",       "regimeProd5",       "regimeProd7",       "regimeProd10" ),]

ggplot(df_alpha_pbias2) +  
coord_cartesian(ylim = c(-100,100))+ 
       geom_boxplot(aes(x=as.factor(by),y=pbias ),
        outlier.shape = NA, alpha=.9) +  
       geom_hline(yintercept=0, color="darkred",alpha=.7,linewidth=1.2) +
       mytheme+ 
       ylab("parameters") +
       xlab("year") +
       facet_grid(scenario~model, scales="free_y")

#=========================================================
aggregate((resparam$pbias/100)*resparam$sim, 
    list(model=resparam$model,
         method=resparam$method,
         iteration=resparam$iteration,
         parameter=resparam$parameter,



#pbias boxplots

resparam<-resparam[resparam$convergence<1,]

dfmpbias<-aggregate(resparam$pbias, 
    list(model=resparam$model,
         method=resparam$method,
         iteration=resparam$iteration,
         parameter=resparam$parameter,
         scenario=resparam$scenario),
    mean)

dfmpbias$model <- factor(dfmpbias$model, 
    levels=c("simple",
             "autocorr",
             "rwa",
             "rwb",
             "rwab",
             "hmma",
             "hmmb",
             "hmmab"))

dfmpbias <- dfmpbias[!is.na(dfmpbias$x),] 



dfmpbias$simulated <- dplyr::recode(dfmpbias$scenario, 
      "trendLinearProd1"="rwa",
      "trendLinearProd2"="rwa",
      "trendLinearProd5"="rwa",
      "trendLinearProd7"="rwa",
      "trendLinearProd10"="rwa", 
      "regimeProd1"="hmma",
      "regimeProd2"="hmma",
      "regimeProd5"="hmma", 
      "regimeProd7"="hmma",
      "regimeProd10"="hmma"
      )   

dfmpbias$simulated_f<-factor(dfmpbias$simulated, levels=c("simple",
                                            "autocorr",
                                            "rwa",
                                            "rwb",
                                            "rwab", 
                                            "hmma",
                                            "hmmb",
                                            "hmmab"))


     
dfmpbias$scenario_f<- factor(dfmpbias$scenario, levels=c(
      "trendLinearProd1",
      "trendLinearProd2",
      "trendLinearProd5",
      "trendLinearProd7", 
      "trendLinearProd10",
      "regimeProd1",
      "regimeProd2",
      "regimeProd5", 
      "regimeProd7",
      "regimeProd10"
      )  )


dfmpbias$model_agg <- dplyr::recode(dfmpbias$model, 
       "simple"="simple", 
       "autocorr"="simple",
       "rwa"="rw",
       "rwb"="rw",
       "rwab"="rw",
       "hmma"="hmm",  
       "hmmb"="hmm",    
       "hmmab"="hmm")



dfmpbias$simulated_agg <- dplyr::recode(dfmpbias$simulated, 
       "simple"="simple", 
       "autocorr"="simple",
       "rwa"="rw",
       "rwb"="rw",
       "rwab"="rw",
       "hmma"="hmm",  
       "hmmb"="hmm",    
       "hmmab"="hmm")


dfmpbias$model_paragg <- dplyr::recode(dfmpbias$model, 
       "simple"="simple", 
       "autocorr"="simple",
       "rwa"="tva",
       "rwb"="tvb",
       "rwab"="tvab",
       "hmma"="tva",  
       "hmmb"="tvb",    
       "hmmab"="tvab")



dfmpbias$simulated_paragg <- dplyr::recode(dfmpbias$simulated, 
       "simple"="simple", 
       "autocorr"="simple",
       "rwa"="tva",
       "rwb"="tvb",
       "rwab"="tvab",
       "hmma"="tva",  
       "hmmb"="tvb",    
       "hmmab"="tvab")



#all scn pbias plot

dfscen<-aggregate(dfmpbias$simulated,list(parameter=dfmpbias$parameter,
  scenario_f=dfmpbias$scenario_f, model=dfmpbias$model),unique)
dfscen$x<-factor(dfscen$x, levels=c("simple", 
       "autocorr",
       "rwa",
       "rwb",
       "rwab",
       "hmma",  
        "hmmb",    
        "hmmab"))
dfscen$simulated_f<-as.numeric(dfscen$x)
dfmpbias <- dfmpbias[!is.na(dfmpbias$x),] 
head(dfmpbias)

boxplot_allscn(dfmpbias,dfscen)


#split parameter plot

dfmpbias_main<-dfmpbias[dfmpbias$parameter%in%c("alpha",
      "Smax",
      "sigma"),]

head(dfmpbias_main)
dfscen_main<-dfscen[dfscen$parameter%in%c("alpha",
      "Smax",
      "sigma"),]

boxplot_allscn(dfmpbias_main,dfscen_main)

dfmpbias_deriv<-dfmpbias[dfmpbias$parameter%in%c("sgen",
      "smsy",
      "umsy"),]


dfscen_deriv<-dfscen[dfscen$parameter%in%c("sgen",
      "smsy",
      "umsy"),]

boxplot_allscn(dfmpbias_deriv,dfscen_deriv)



#aggegated plot

dfscen_agg<-aggregate(dfmpbias$simulated_agg,list(parameter=dfmpbias$parameter,
  scenario_f=dfmpbias$scenario_f, model=dfmpbias$model_agg),unique)
dfscen_agg$x<-factor(dfscen_agg$x, levels=c("simple", 
       "rw",
       "hmm"))
dfscen_agg$simulated_f<-as.numeric(dfscen_agg$x)

dfscen_main_agg<-dfscen_agg[dfscen_agg$parameter%in%c("alpha",
      "Smax",
      "sigma"),]


boxplot_agg(dfmpbias=dfmpbias_main,dfscen_agg=dfscen_main_agg)





#absolute bias
dfmbias<-aggregate((resparam$pbias/100)*resparam$sim, 
    list(model=resparam$model,
         method=resparam$method,
         iteration=resparam$iteration,
         parameter=resparam$parameter,
         scenario=resparam$scenario),
    mean)

unique(dfmpbias$model)
dfmbias$model <- factor(dfmbias$model, 
    levels=c("simple",
             "autocorr",
             "rwa",
             "rwb",
             "rwab",
             "hmma",
             "hmmb",
             "hmmab"))

dfmbias <- dfmbias[!is.na(dfmbias$x),] 


dfmbias$simulated <- dplyr::recode(dfmbias$scenario, 
      "trendLinearProd1"="rwa",
      "trendLinearProd2"="rwa",
      "trendLinearProd5"="rwa",
      "trendLinearProd7"="rwa",
      "trendLinearProd10"="rwa", 
      "regimeProd1"="hmma",
      "regimeProd2"="hmma",
      "regimeProd5"="hmma", 
      "regimeProd7"="hmma",
       "regimeProd10"="hmma"
      )   

dfmbias$simulated_f<-factor(dfmbias$simulated, levels=c("simple",
                                            "autocorr",
                                            "rwa",
                                            "rwb",
                                            "rwab", 
                                            "hmma",
                                            "hmmb",
                                            "hmmab"))


     
dfmbias$scenario_f<- factor(dfmbias$scenario, levels=c(
      "trendLinearProd1",
      "trendLinearProd2",
      "trendLinearProd5",
      "trendLinearProd7",
      "trendLinearProd10",  
      "regimeProd1",
      "regimeProd2",
      "regimeProd5", 
      "regimeProd7",
      "regimeProd10"
      )  )


dfmbias$model_agg <- dplyr::recode(dfmbias$model, 
       "simple"="simple", 
       "autocorr"="simple",
       "rwa"="rw",
       "rwb"="rw",
       "rwab"="rw",
       "hmma"="hmm",  
       "hmmb"="hmm",    
       "hmmab"="hmm")



dfmbias$simulated_agg <- dplyr::recode(dfmbias$simulated, 
       "simple"="simple", 
       "autocorr"="simple",
       "rwa"="rw",
       "rwb"="rw",
       "rwab"="rw",
       "hmma"="hmm",  
       "hmmb"="hmm",    
       "hmmab"="hmm")


dfmbias$model_paragg <- dplyr::recode(dfmbias$model, 
       "simple"="simple", 
       "autocorr"="simple",
       "rwa"="tva",
       "rwb"="tvb",
       "rwab"="tvab",
       "hmma"="tva",  
       "hmmb"="tvb",    
       "hmmab"="tvab")



dfmbias$simulated_paragg <- dplyr::recode(dfmbias$simulated, 
       "simple"="simple", 
       "autocorr"="simple",
       "rwa"="tva",
       "rwb"="tvb",
       "rwab"="tvab",
       "hmma"="tva",  
       "hmmb"="tvb",    
       "hmmab"="tvab")



#all scn pbias plot

dfscen<-aggregate(dfmbias$simulated,list(parameter=dfmbias$parameter,
  scenario_f=dfmbias$scenario_f, model=dfmbias$model),unique)
dfscen$x<-factor(dfscen$x, levels=c("simple", 
       "autocorr",
       "rwa",
       "rwb",
       "rwab",
       "hmma",  
        "hmmb",    
        "hmmab"))
dfscen$simulated_f<-as.numeric(dfscen$x)


dfmbias_alpha<-dfmbias[dfmbias$parameter%in%c("alpha"),]

dfscen_alpha<-dfscen[dfscen$parameter%in%c("alpha"),]

head(dfmbias_alpha)

ggplot(dfmbias_alpha) +   
  geom_boxplot(aes(x=model,y=x), outlier.shape = NA, alpha=.7) + 
  coord_cartesian(ylim = quantile(dfmbias_alpha$x, c(0.025, 0.975)))+        
  geom_hline(yintercept=0) +
  mytheme+ 
  ylab("mean bias") +
  xlab("estimation model") +
  facet_grid(parameter~scenario_f, scales="free_y")+
  scale_fill_viridis_d(begin=.1, end=.8) +
  scale_color_viridis_d(begin=.1, end=.8) +
  geom_rect(data=dfscen_alpha,aes(xmin=as.numeric(simulated_f)-.5,
                         xmax=as.numeric(simulated_f)+.5,
                         ymin=-200,ymax=Inf),
                    color="gray90",alpha=0.05)+
  theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1))




dfmpbias_alpha<-dfmpbias[dfmpbias$parameter%in%c("sgen",
      "smsy",    
      "Smax"),]

dfscen_main<-dfscen[dfscen$parameter%in%c("alpha",
      "Smax",
      "sigma"),]


ggplot(dfmbias) +   
  geom_boxplot(aes(fill=method,x=model,y=x), outlier.shape = NA, alpha=.7) +         
  geom_hline(yintercept=0) +
  mytheme+ 
  ylab("mean bias") +
  xlab("estimation model") +
  facet_grid(parameter~scenario_f, scales="free_y")+
  scale_fill_viridis_d(begin=.1, end=.8) +
  scale_color_viridis_d(begin=.1, end=.8) +
  geom_rect(data=dfscen,aes(xmin=as.numeric(simulated_f)-.5,
                         xmax=as.numeric(simulated_f)+.5,
                         ymin=-Inf,ymax=Inf),
                    color="gray90",alpha=0.05)+
  #stat_summary(aes(color=method,x=(model),y=x),fun.data = give.n, geom = "text", 
  #        vjust = -2,position = position_dodge(1))+ 
  theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1))






#========================================================================================================
#sensitivity a scenario
#read in data
simPar <- read.csv("data/Smax_sensitivity/SimPars.csv")

## Store relevant object names to help run simulation 
scenNames <- unique(simPar$scenario)

res_smax<-readRDS(file = "outs/simest/Smax_sensitivity/res_smax.rds")

res_smax_stan<-readRDS(file = "outs/simest/resstan_smax.rds")


res_smax<-rbind(res_smax,res_smax_stan)

res_smax$parameter[res_smax$parameter=="smax"]<-"Smax"

resparam<-res_smax[res_smax$parameter%in%c("alpha","Smax","sigma","smsy","sgen","umsy"),]

df<-reshape2::melt(resparam, id.vars=c("parameter","iteration","scenario","method","model","by", "convergence","pbias"))

df_Smax<-df[df$parameter%in%c("Smax"),]

df_Smax$col<-factor(df_Smax$variable,levels=c("est","sim"))

df_Smax_mle<-df_Smax[df_Smax$method=="MLE",]


unique(df_Smax_mle$scenario)
df_Smax_mle1<-df_Smax_mle[df_Smax_mle$model!="hmmb",]

param_traj(df_Smax_mle1)

df_Smax_mcmc<-df_Smax[df_Smax$method=="MCMC",]

param_traj(df_Smax_mcmc)



#pbias boxplots -- need to make these plots

resparam<-resparam[resparam$convergence<1,]

dfmpbias<-aggregate(resparam$pbias, 
    list(model=resparam$model,
         method=resparam$method,
         iteration=resparam$iteration,
         parameter=resparam$parameter,
         scenario=resparam$scenario),
    mean)

dfmpbias$model <- factor(dfmpbias$model, 
    levels=c("simple",
             "autocorr",
             "rwa",
             "rwb",
             "rwab",
             "hmma",
             "hmmb",
             "hmmab"))

dfmpbias <- dfmpbias[!is.na(dfmpbias$x),] 

unique(dfmpbias$scenario)

dfmpbias$simulated <- dplyr::recode(dfmpbias$scenario, 
    "regimeSmax025"="hmmb", 
    "regimeSmax050"="hmmb",      
    "regimeSmax150"="hmmb",      
    "regimeSmax200"="hmmb",      
    "regimeSmax300"="hmmb",      
    "trendLinearSmax025"="rwb",
    "trendLinearSmax050"="rwb", 
    "trendLinearSmax150"="rwb", 
    "trendLinearSmax200"="rwb",
    "trendLinearSmax300"="rwb"
      )   

dfmpbias$simulated_f<-factor(dfmpbias$simulated, levels=c("simple",
                                            "autocorr",
                                            "rwa",
                                            "rwb",
                                            "rwab", 
                                            "hmma",
                                            "hmmb",
                                            "hmmab"))


     
dfmpbias$scenario_f<- factor(dfmpbias$scenario, levels=c(
      "regimeSmax025", 
    "regimeSmax050",      
    "regimeSmax150",      
    "regimeSmax200",      
    "regimeSmax300",      
    "trendLinearSmax025",
    "trendLinearSmax050", 
    "trendLinearSmax150", 
    "trendLinearSmax200",
    "trendLinearSmax300"
      )  )


dfmpbias$model_agg <- dplyr::recode(dfmpbias$model, 
       "simple"="simple", 
       "autocorr"="simple",
       "rwa"="rw",
       "rwb"="rw",
       "rwab"="rw",
       "hmma"="hmm",  
       "hmmb"="hmm",    
       "hmmab"="hmm")



dfmpbias$simulated_agg <- dplyr::recode(dfmpbias$simulated, 
       "simple"="simple", 
       "autocorr"="simple",
       "rwa"="rw",
       "rwb"="rw",
       "rwab"="rw",
       "hmma"="hmm",  
       "hmmb"="hmm",    
       "hmmab"="hmm")


dfmpbias$model_paragg <- dplyr::recode(dfmpbias$model, 
       "simple"="simple", 
       "autocorr"="simple",
       "rwa"="tva",
       "rwb"="tvb",
       "rwab"="tvab",
       "hmma"="tva",  
       "hmmb"="tvb",    
       "hmmab"="tvab")



dfmpbias$simulated_paragg <- dplyr::recode(dfmpbias$simulated, 
       "simple"="simple", 
       "autocorr"="simple",
       "rwa"="tva",
       "rwb"="tvb",
       "rwab"="tvab",
       "hmma"="tva",  
       "hmmb"="tvb",    
       "hmmab"="tvab")



#all scn pbias plot

dfscen<-aggregate(dfmpbias$simulated,list(parameter=dfmpbias$parameter,
  scenario_f=dfmpbias$scenario_f, model=dfmpbias$model),unique)
dfscen$x<-factor(dfscen$x, levels=c("simple", 
       "autocorr",
       "rwa",
       "rwb",
       "rwab",
       "hmma",  
        "hmmb",    
        "hmmab"))
dfscen$simulated_f<-as.numeric(dfscen$x)
dfmpbias <- dfmpbias[!is.na(dfmpbias$x),] 
head(dfmpbias)

boxplot_allscn(dfmpbias,dfscen)


#split parameter plot

dfmpbias_main<-dfmpbias[dfmpbias$parameter%in%c("alpha",
      "Smax",
      "sigma"),]

head(dfmpbias_main)
dfscen_main<-dfscen[dfscen$parameter%in%c("alpha",
      "Smax",
      "sigma"),]

boxplot_allscn(dfmpbias_main,dfscen_main)

dfmpbias_deriv<-dfmpbias[dfmpbias$parameter%in%c("sgen",
      "smsy",
      "umsy"),]


dfscen_deriv<-dfscen[dfscen$parameter%in%c("sgen",
      "smsy",
      "umsy"),]

boxplot_allscn(dfmpbias_deriv,dfscen_deriv)



#aggegated plot

dfscen_agg<-aggregate(dfmpbias$simulated_agg,list(parameter=dfmpbias$parameter,
  scenario_f=dfmpbias$scenario_f, model=dfmpbias$model_agg),unique)
dfscen_agg$x<-factor(dfscen_agg$x, levels=c("simple", 
       "rw",
       "hmm"))
dfscen_agg$simulated_f<-as.numeric(dfscen_agg$x)

dfscen_main_agg<-dfscen_agg[dfscen_agg$parameter%in%c("alpha",
      "Smax",
      "sigma"),]


boxplot_agg(dfmpbias=dfmpbias_main,dfscen_agg=dfscen_main_agg)





#absolute bias
dfmbias<-aggregate((resparam$pbias/100)*resparam$sim, 
    list(model=resparam$model,
         method=resparam$method,
         iteration=resparam$iteration,
         parameter=resparam$parameter,
         scenario=resparam$scenario),
    mean)
summary(dfmpbias$x)
unique(dfmpbias$model)
dfmbias$model <- factor(dfmbias$model, 
    levels=c("simple",
             "autocorr",
             "rwa",
             "rwb",
             "rwab",
             "hmma",
             "hmmb",
             "hmmab"))

dfmbias <- dfmbias[!is.na(dfmbias$x),] 


dfmbias$simulated <- dplyr::recode(dfmbias$scenario, 
      "trendLinearProd1"="rwa",
      "trendLinearProd2"="rwa",
      "trendLinearProd5"="rwa",
      "trendLinearProd7"="rwa",
      "trendLinearProd10"="rwa", 
      "regimeProd1"="hmma",
      "regimeProd2"="hmma",
      "regimeProd5"="hmma", 
      "regimeProd7"="hmma",
       "regimeProd10"="hmma"
      )   

dfmbias$simulated_f<-factor(dfmbias$simulated, levels=c("simple",
                                            "autocorr",
                                            "rwa",
                                            "rwb",
                                            "rwab", 
                                            "hmma",
                                            "hmmb",
                                            "hmmab"))


     
dfmbias$scenario_f<- factor(dfmbias$scenario, levels=c(
      "trendLinearProd1",
      "trendLinearProd2",
      "trendLinearProd5",
      "trendLinearProd7",
      "trendLinearProd10",  
      "regimeProd1",
      "regimeProd2",
      "regimeProd5", 
      "regimeProd7",
      "regimeProd10"
      )  )


dfmbias$model_agg <- dplyr::recode(dfmbias$model, 
       "simple"="simple", 
       "autocorr"="simple",
       "rwa"="rw",
       "rwb"="rw",
       "rwab"="rw",
       "hmma"="hmm",  
       "hmmb"="hmm",    
       "hmmab"="hmm")



dfmbias$simulated_agg <- dplyr::recode(dfmbias$simulated, 
       "simple"="simple", 
       "autocorr"="simple",
       "rwa"="rw",
       "rwb"="rw",
       "rwab"="rw",
       "hmma"="hmm",  
       "hmmb"="hmm",    
       "hmmab"="hmm")


dfmbias$model_paragg <- dplyr::recode(dfmbias$model, 
       "simple"="simple", 
       "autocorr"="simple",
       "rwa"="tva",
       "rwb"="tvb",
       "rwab"="tvab",
       "hmma"="tva",  
       "hmmb"="tvb",    
       "hmmab"="tvab")



dfmbias$simulated_paragg <- dplyr::recode(dfmbias$simulated, 
       "simple"="simple", 
       "autocorr"="simple",
       "rwa"="tva",
       "rwb"="tvb",
       "rwab"="tvab",
       "hmma"="tva",  
       "hmmb"="tvb",    
       "hmmab"="tvab")



#all scn pbias plot

dfscen<-aggregate(dfmbias$simulated,list(parameter=dfmbias$parameter,
  scenario_f=dfmbias$scenario_f, model=dfmbias$model),unique)
dfscen$x<-factor(dfscen$x, levels=c("simple", 
       "autocorr",
       "rwa",
       "rwb",
       "rwab",
       "hmma",  
        "hmmb",    
        "hmmab"))
dfscen$simulated_f<-as.numeric(dfscen$x)


dfmbias_alpha<-dfmbias[dfmbias$parameter%in%c("alpha"),]

dfscen_alpha<-dfscen[dfscen$parameter%in%c("alpha"),]


ggplot(dfmbias_alpha) +   
  geom_boxplot(aes(x=model,y=x), outlier.shape = NA, alpha=.7) + 
  coord_cartesian(ylim = quantile(dfmbias_alpha$x, c(0.025, 0.975)))+        
  geom_hline(yintercept=0) +
  mytheme+ 
  ylab("mean bias") +
  xlab("estimation model") +
  facet_grid(parameter~scenario_f, scales="free_y")+
  scale_fill_viridis_d(begin=.1, end=.8) +
  scale_color_viridis_d(begin=.1, end=.8) +
  geom_rect(data=dfscen_alpha,aes(xmin=as.numeric(simulated_f)-.5,
                         xmax=as.numeric(simulated_f)+.5,
                         ymin=-200,ymax=Inf),
                    color="gray90",alpha=0.05)+
  theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1))




dfmpbias_alpha<-dfmpbias[dfmpbias$parameter%in%c("sgen",
      "smsy",    
      "Smax"),]

dfscen_main<-dfscen[dfscen$parameter%in%c("alpha",
      "Smax",
      "sigma"),]


ggplot(dfmbias) +   
  geom_boxplot(aes(fill=method,x=model,y=x), outlier.shape = NA, alpha=.7) +         
  geom_hline(yintercept=0) +
  mytheme+ 
  ylab("mean bias") +
  xlab("estimation model") +
  facet_grid(parameter~scenario_f, scales="free_y")+
  scale_fill_viridis_d(begin=.1, end=.8) +
  scale_color_viridis_d(begin=.1, end=.8) +
  geom_rect(data=dfscen,aes(xmin=as.numeric(simulated_f)-.5,
                         xmax=as.numeric(simulated_f)+.5,
                         ymin=-Inf,ymax=Inf),
                    color="gray90",alpha=0.05)+
  #stat_summary(aes(color=method,x=(model),y=x),fun.data = give.n, geom = "text", 
  #        vjust = -2,position = position_dodge(1))+ 
  theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1))






#========================================================================================================
#base case
#read in data
simPar <- read.csv("data/genericER/SimPars_ER.csv")

## Store relevant object names to help run simulation 
scenNames <- unique(simPar$scenario)

res<-readRDS(file = "outs/simest/generic/res.rds")

#
#resstan1<-readRDS(file = "outs/simest/generic/resstan1.rds")
#resstan2<-readRDS(file = "outs/simest/generic/resstan2.rds")
#resstan<-rbind(resstan1,resstan2)
resstan<-readRDS(file = "outs/simest/generic/resstan.rds")

res<-rbind(res,resstan)

res$parameter[res$parameter=="smax"]<-"Smax"

resparam<-res[res$parameter%in%c("alpha","Smax","sigma","smsy","sgen","umsy"),]

#mean by iteration
dfmpbias<-aggregate(resparam$pbias, 
    list(model=resparam$model,
         method=resparam$method,
         iteration=resparam$iteration,
         parameter=resparam$parameter,
         scenario=resparam$scenario),
    mean)

unique(dfmpbias$method)
dfmpbias$model <- factor(dfmpbias$model, 
    levels=c("simple",
             "autocorr",
             "rwa",
             "rwb",
             "rwab",
             "hmma",
             "hmmb",
             "hmmab"))
#=================================
#plots
#dfmpbias$method <- factor(dfmpbias$method, levels=c("MLE","MCMC"))
dfmpbias <- dfmpbias[!is.na(dfmpbias$x),] 

dfmpbias$simulated <- dplyr::recode(dfmpbias$scenario, 
      "stationary"="simple",
      "autocorr"="autocorr",
      "sigmaShift"="simple", 
      "decLinearProd"="rwa",
      "sineProd"="rwa",
      "regimeProd"="hmma",
      "decLinearCap"="rwb",
      "regimeCap"="hmmb",
      "shiftCap"="hmmb", 
      "shiftProd"="hmma",
      "regimeProdCap"="hmmab",
      "decLinearProdshiftCap"="rwab"
      )   

dfmpbias$simulated_f<-factor(dfmpbias$simulated, levels=c("simple",
                                            "autocorr",
                                            "rwa",
                                            "rwb",
                                            "rwab", 
                                            "hmma",
                                            "hmmb",
                                            "hmmab"))



dfmpbias$model <- factor(dfmpbias$model, levels=c(
       "simple", 
       "autocorr",
       "rwa",
       "rwb",
       "rwab",
       "hmma",  
        "hmmb",    
        "hmmab") )

     
dfmpbias$scenario_f<- factor(dfmpbias$scenario, levels=c(
      "stationary",
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
      "decLinearProdshiftCap"
      )  )


dfmpbias$model_agg <- dplyr::recode(dfmpbias$model, 
       "simple"="simple", 
       "autocorr"="simple",
       "rwa"="rw",
       "rwb"="rw",
       "rwab"="rw",
       "hmma"="hmm",  
       "hmmb"="hmm",    
       "hmmab"="hmm")

dfmpbias$simulated_agg <- dplyr::recode(dfmpbias$simulated, 
       "simple"="simple", 
       "autocorr"="simple",
       "rwa"="rw",
       "rwb"="rw",
       "rwab"="rw",
       "hmma"="hmm",  
       "hmmb"="hmm",    
       "hmmab"="hmm")

dfmpbias$model_paragg <- dplyr::recode(dfmpbias$model, 
       "simple"="simple", 
       "autocorr"="simple",
       "rwa"="tva",
       "rwb"="tvb",
       "rwab"="tvab",
       "hmma"="tva",  
       "hmmb"="tvb",    
       "hmmab"="tvab")

dfmpbias$simulated_paragg <- dplyr::recode(dfmpbias$simulated, 
       "simple"="simple", 
       "autocorr"="simple",
       "rwa"="tva",
       "rwb"="tvb",
       "rwab"="tvab",
       "hmma"="tva",  
       "hmmb"="tvb",    
       "hmmab"="tvab")

dfscen<-aggregate(dfmpbias$simulated,list(parameter=dfmpbias$parameter,
  scenario_f=dfmpbias$scenario_f, model=dfmpbias$model),unique)
dfscen$x<-factor(dfscen$x, levels=c("simple", 
       "autocorr",
       "rwa",
       "rwb",
       "rwab",
       "hmma",  
        "hmmb",    
        "hmmab"))
dfscen$simulated_f<-as.numeric(dfscen$x)

dfmpbias_main<-dfmpbias[dfmpbias$parameter%in%c("alpha",
      "Smax",
      "sigma"),]

dfscen_main<-dfscen[dfscen$parameter%in%c("alpha",
      "Smax",
      "sigma"),]

boxplot_allscn(dfmpbias_main,dfscen_main)



dfmpbias_deriv<-dfmpbias[dfmpbias$parameter%in%c("sgen",
      "smsy",
      "umsy"),]


dfscen_deriv<-dfscen[dfscen$parameter%in%c("sgen",
      "smsy",
      "umsy"),]

boxplot_allscn(dfmpbias_deriv,dfscen_deriv)

#aggregate plots




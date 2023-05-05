#routines to visualize smulation estimation Results
#Catarina Wor
#November 2022
#============================================



library(ggplot2)
library(gridExtra)
source("code/utils.R")


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

res<-readRDS(file = "outs/simest/res.rds")
res14<-readRDS(file = "outs/simest/res14.rds")



#
resstan<-readRDS(file = "outs/simest/resstan.rds")


res<-rbind(res,res14,resstan)#,resstan16,resstan712)


res$parameter[res$parameter=="smax"]<-"Smax"
#mean by iteration
resparam<-res[res$parameter%in%c("alpha","Smax","sigma","smsy","sgen","umsy"),]

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




ggplot(dfmpbias) +   
  geom_boxplot(data=dfmpbias,aes(fill=method,x=model,y=x), outlier.shape = NA, alpha=.7) +         
  coord_cartesian(ylim = quantile(dfmpbias$x, c(0.025, 0.975)))+
  geom_hline(yintercept=0) +
  mytheme+ 
  ylab("mean % bias") +
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


#split parameter plot

dfmpbias_main<-dfmpbias[dfmpbias$parameter%in%c("alpha",
      "Smax",
      "sigma"),]


dfscen_main<-dfscen[dfscen$parameter%in%c("alpha",
      "Smax",
      "sigma"),]

ggplot(dfmpbias_main) +   
  geom_boxplot(data=dfmpbias_main,aes(fill=method,x=model,y=x), outlier.shape = NA,alpha=.9) +         
  coord_cartesian(ylim = quantile(dfmpbias_main$x, c(0.025, 0.975)))+
  geom_hline(yintercept=0) +
  mytheme+ 
  ylab("mean % bias") +
  xlab("estimation model") +
  facet_grid(parameter~scenario_f, scales="free_y")+
  scale_fill_viridis_d(begin=.1, end=.8,option = "D",) +
  scale_color_viridis_d(begin=.4, end=.8) +
  geom_rect(data=dfscen_main,aes(xmin=as.numeric(simulated_f)-.5,
                         xmax=as.numeric(simulated_f)+.5,
                         ymin=-Inf,ymax=Inf),
                    color="gray90",alpha=0.05)+
  #stat_summary(aes(color=method,x=(model),y=x),fun.data = give.n, geom = "text", 
  #        vjust = -2,position = position_dodge(1))+ 
  theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1))


dfmpbias_deriv<-dfmpbias[dfmpbias$parameter%in%c("sgen",
      "smsy",
      "umsy"),]


dfscen_deriv<-dfscen[dfscen$parameter%in%c("sgen",
      "smsy",
      "umsy"),]


ggplot(dfmpbias_deriv) +   
  geom_boxplot(data=dfmpbias_deriv,aes(fill=method,x=model,y=x), outlier.shape = NA,alpha=.3) +         
  coord_cartesian(ylim = quantile(dfmpbias_main$x, c(0.025, 0.975)))+
  geom_hline(yintercept=0) +
  mytheme+ 
  ylab("mean % bias") +
  xlab("estimation model") +
  facet_grid(parameter~scenario_f, scales="free_y")+
  scale_fill_viridis_d(begin=.1, end=.8,option = "D",) +
  scale_color_viridis_d(begin=.4, end=.8) +
  geom_rect(data=dfscen_deriv,aes(xmin=as.numeric(simulated_f)-.5,
                         xmax=as.numeric(simulated_f)+.5,
                         ymin=-Inf,ymax=Inf),
                    color="gray90",alpha=0.05)+
  #stat_summary(aes(color=method,x=(model),y=x),fun.data = give.n, geom = "text", 
  #        vjust = -2,position = position_dodge(1))+ 
  theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1))

#close look plot


dfmpbias_ss<-dfmpbias[dfmpbias$scenario_f%in%c("stationary",
      "decLinearProd",
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
  scale_fill_viridis_d(begin=.3, end=.8) +
  scale_color_viridis_d(begin=.3, end=.8) +
  geom_rect(data=dfscen_ss,aes(xmin=as.numeric(simulated_f)-.5,
                         xmax=as.numeric(simulated_f)+.5,
                         ymin=-Inf,ymax=Inf),
                    color="gray90",alpha=0.05)+
  #stat_summary(aes(color=method,x=(model),y=x),fun.data = give.n, geom = "text", 
  #        vjust = -2,position = position_dodge(1))+ 
  theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1))






ggsave(
      filename = "outs/SamSimOutputs/plotcheck/pbiasbymodel.pdf", 
      plot = marrangeGrob(pbiasmodel, nrow=1, ncol=1), 
      width = 12, height = 8
    )


ggsave(
      filename = "outs/SamSimOutputs/plotcheck/rmsebymodel.pdf", 
      plot = marrangeGrob(rmsemodel, nrow=1, ncol=1), 
      width = 12, height = 8
    )


#===========================================
#confusion matrices 




#summarize model selection
aic_set=list()
bic_set=list()
lfo_set=list()

aic=subset(res,parameter=='AIC')
bic=subset(res,parameter=='BIC')
lfo=subset(res,parameter=='LFO')
dim(lfo)
summary(aic)
aggregate(lfo$iteration,list(lfo$iteration),length)

lfo<-lfo[lfo$model %in% c("simple","autocorr","rwa","rwb","rwab","hmma","hmmb","hmmab"),]
unique(lfo$model)
lfo[is.na(lfo$est),]<--Inf

scn<-factor(unique(aic$scenario), levels=c(
  "stationary", 
  "sigmaShift",
  "autocorr",
  "decLinearProd",  
  "sineProd",
  "decLinearCap",
  "decLinearProdshiftCap",      
 "regimeProd",
 "shiftProd", 
 "regimeCap", 
 "shiftCap",
 "regimeProdCap" ) )


EM=c("stationary",
     "autocorr",
     "dynamic.a","dynamic.b","dynamic.ab",
     "regime.a","regime.b","regime.ab")
##Confusion matrices
conf_matrix<-expand.grid(EM=EM,OM=scn)
conf_matrix$w_AIC=NA
conf_matrix$BIC=NA
conf_matrix$LFO=NA

cn1<-list()
cn2<-list()
cn3<-list()

o=0
for(a in seq_along(scn)){

  #AIC
  aica<-subset(aic,scenario==scn[a])
  aic_set[[a]]=tidyr::spread(aica[,-9],key=model,value=est)
  aic_set[[a]]=aic_set[[a]][c(15,8,12,14,13,9,11,10)] #reorder estimation models

  sc1=apply(aic_set[[a]],1,which.min)
  cn1[[a]]=summary(factor(sc1,levels=seq(1:ncol(aic_set[[a]]))))/1000

  bica=subset(bic,scenario==scn[a])
  bic_set[[a]]=tidyr::spread(bica[,-9],key=model,value=est)
  bic_set[[a]]=bic_set[[a]][c(15,8,12,14,13,9,11,10)] #reorder estimation models

  sc2=apply(bic_set[[a]],1,which.min)
  cn2[[a]]=summary(factor(sc2,levels=seq(1:ncol(bic_set[[a]]))))/1000

  lfoa=subset(lfo,scenario==scn[a])
  lfo_set[[a]]=tidyr::spread(lfoa[,-9],key=model,value=est)
  lfo_set[[a]]=lfo_set[[a]][c(15,8,12,14,13,9,11,10)] #reorder estimation models

  sc3=apply(lfo_set[[a]],1,which.max)
  cn3[[a]]=summary(factor(sc3,levels=seq(1:ncol(lfo_set[[a]]))))/1000


  myseq<-seq(from=o+1, length.out=length(EM))
  conf_matrix$w_AIC[myseq]<-cn1[[a]]
  conf_matrix$BIC[myseq]<-cn2[[a]]
  conf_matrix$LFO[myseq]<-cn3[[a]]
  o=max(myseq)

}


conf_matrix$eqem_om <- dplyr::recode(conf_matrix$OM, 
      "stationary"="stationary",
      "autocorr"="autocorr",
      "sigmaShift"="stationary", 
      "decLinearProd"="dynamic.a",
      "sineProd"="dynamic.a",
      "decLinearCap"="dynamic.b",
      "decLinearProdshiftCap"="dynamic.ab",
      "regimeProd"="regime.a",
      "shiftProd"="regime.a",
      "regimeCap"="regime.b",
      "shiftCap"="regime.b", 
      "regimeProdCap"="regime.ab",
      )   
conf_matrix$diag<-conf_matrix$eqem_om==conf_matrix$EM

conf_matrix$EM


mytheme = list(
  theme_classic(14)+
    theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
          legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="bold", size=13),
          axis.text=element_text(face="bold"),axis.title = element_text(face="bold",size=16),plot.title = element_text(face = "bold", hjust = 0.5,size=16))
)





p=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = w_AIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_AIC,2)), vjust = 1,size=6) +
  ggtitle("AIC")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray90", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1)) +
  xlab("Simulation Scenario")+ylab("Estimation Model")
p



p=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = BIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(BIC,2)), vjust = 1,size=6) +
  ggtitle("BIC")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray90", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
p



p=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = LFO), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(LFO,2)), vjust = 1, size=6) +
  ggtitle("LFO")+
  scale_fill_gradient(low="white", high="#009194")  +
  geom_segment(data=transform(subset(conf_matrix, !!diag), 
                    simulated=as.numeric(OM), 
                    estimated=as.numeric(EM)), 
               aes(x=simulated-.49, xend=simulated+.49, y=estimated-.49, yend=estimated+.49), 
               color="gray90", linewidth=2)+
  mytheme + theme(legend.position="none", axis.text.x = element_text(angle = 45,  hjust=1))+
  xlab("Simulation Scenario")+ylab("Estimation Model")
p




#========================================================================================================
#sensitivity a scenario
#read in data
simPar <- read.csv("data/sensitivity/SimPars.csv")

## Store relevant object names to help run simulation 
scenNames <- unique(simPar$scenario)

res_a<-readRDS(file = "outs/simest/res_a.rds")

unique(res_a$parameter)
head(res_a[res_a$sim==0,])




#mean by iteration
resparam<-res_a[res_a$parameter%in%c("alpha","Smax","sigma","smsy","sgen","umsy"),]





resparam<-resparam[resparam$scenario!="regimeProd1",]
resparam<-resparam[resparam$convergence<1,]
summary(resparam)

 df<-reshape2::melt(resparam, id.vars=c("parameter","iteration","scenario","method","model","by", "convergence","pbias"))
 df_alpha<-df[df$parameter%in%c("alpha"),]

df_alpha$col<-factor(df_alpha$variable,levels=c("est","sim"))

ggplot(df_alpha) +   
  geom_line(aes(x=by,y=value,color=col,
    group= interaction(iteration,col) ),
    size=1.4,alpha=.5) +  
  mytheme+ 
  ylab("parameters") +
  xlab("year") +
  facet_grid(scenario~model, scales="free_y")+
   scale_color_manual(values = c(sim = "darkred", est = "gray50"))
  

df_alpha<-df_alpha[df_alpha$scenario%in%c("trendLinearProd1"),]
head(df_alpha)

df_alpha$col<-factor(df_alpha$variable,levels=c("est","sim"))

ggplot(df_alpha) +   
  geom_line(aes(x=by,y=value,color=col,
    group= interaction(iteration,col) ),
    size=1.4,alpha=.5) +  
  mytheme+ 
  ylab("parameters") +
  xlab("year") +
  facet_grid(~model, scales="free_y")+
   scale_color_manual(values = c(sim = "darkred", est = "gray50"))
  
#scale_color_viridis_d(begin=.1, end=.8)


dfmpbias<-aggregate(resparam$pbias, 
    list(model=resparam$model,
         method=resparam$method,
         iteration=resparam$iteration,
         parameter=resparam$parameter,
         scenario=resparam$scenario),
    mean)
summary(dfmpbias$x)
unique(dfmpbias$model)
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
summary(dfmpbias)


dfmpbias$simulated <- dplyr::recode(dfmpbias$scenario, 
      "trendLinearProd1"="rwa",
      "trendLinearProd2"="rwa",
      "trendLinearProd5"="rwa",
      "trendLinearProd7"="rwa", 
      #"regimeProd1"="hmma",
      "regimeProd2"="hmma",
      "regimeProd5"="hmma", 
      "regimeProd7"="hmma"
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
      #"regimeProd1",
      "regimeProd2",
      "regimeProd5", 
      "regimeProd7"
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




ggplot(dfmpbias) +   
  geom_boxplot(data=dfmpbias,aes(fill=method,x=model,y=x), outlier.shape = NA, alpha=.7) +         
  coord_cartesian(ylim = quantile(dfmpbias$x, c(0.025, 0.975)))+
  geom_hline(yintercept=0) +
  mytheme+ 
  ylab("mean % bias") +
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


#split parameter plot

dfmpbias_main<-dfmpbias[dfmpbias$parameter%in%c("alpha",
      "Smax",
      "sigma"),]

head(dfmpbias_main)
dfscen_main<-dfscen[dfscen$parameter%in%c("alpha",
      "Smax",
      "sigma"),]

ggplot(dfmpbias_main) +   
  geom_boxplot(data=dfmpbias_main,aes(fill=method,x=model,y=x), outlier.shape = NA,alpha=.3) +         
  coord_cartesian(ylim = quantile(dfmpbias_main$x, c(0.025, 0.975)))+
  geom_hline(yintercept=0) +
  mytheme+ 
  ylab("mean % bias") +
  xlab("estimation model") +
  facet_grid(parameter~scenario_f, scales="free_y")+
  scale_fill_viridis_d(begin=.1, end=.8,option = "D",) +
  scale_color_viridis_d(begin=.4, end=.8) +
  geom_rect(data=dfscen_main,aes(xmin=as.numeric(simulated_f)-.5,
                         xmax=as.numeric(simulated_f)+.5,
                         ymin=-Inf,ymax=Inf),
                    color="gray90",alpha=0.05)+
  #stat_summary(aes(color=method,x=(model),y=x),fun.data = give.n, geom = "text", 
  #        vjust = -2,position = position_dodge(1))+ 
  theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1))


dfmpbias_deriv<-dfmpbias[dfmpbias$parameter%in%c("sgen",
      "smsy",
      "umsy"),]


dfscen_deriv<-dfscen[dfscen$parameter%in%c("sgen",
      "smsy",
      "umsy"),]


ggplot(dfmpbias_deriv) +   
  geom_boxplot(data=dfmpbias_deriv,aes(fill=method,x=model,y=x), outlier.shape = NA,alpha=.3) +         
  coord_cartesian(ylim = quantile(dfmpbias_main$x, c(0.025, 0.975)))+
  geom_hline(yintercept=0) +
  mytheme+ 
  ylab("mean % bias") +
  xlab("estimation model") +
  facet_grid(parameter~scenario_f, scales="free_y")+
  scale_fill_viridis_d(begin=.1, end=.8,option = "D",) +
  scale_color_viridis_d(begin=.4, end=.8) +
  geom_rect(data=dfscen_deriv,aes(xmin=as.numeric(simulated_f)-.5,
                         xmax=as.numeric(simulated_f)+.5,
                         ymin=-Inf,ymax=Inf),
                    color="gray90",alpha=0.05)+
  #stat_summary(aes(color=method,x=(model),y=x),fun.data = give.n, geom = "text", 
  #        vjust = -2,position = position_dodge(1))+ 
  theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1))

#close look plot


dfmpbias_ss<-dfmpbias[dfmpbias$scenario_f%in%c("stationary",
      "decLinearProd",
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
  scale_fill_viridis_d(begin=.3, end=.8) +
  scale_color_viridis_d(begin=.3, end=.8) +
  geom_rect(data=dfscen_ss,aes(xmin=as.numeric(simulated_f)-.5,
                         xmax=as.numeric(simulated_f)+.5,
                         ymin=-Inf,ymax=Inf),
                    color="gray90",alpha=0.05)+
  #stat_summary(aes(color=method,x=(model),y=x),fun.data = give.n, geom = "text", 
  #        vjust = -2,position = position_dodge(1))+ 
  theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1))

head(resparam)

#absolute bias
dfmbias<-aggregate(resparam$pbias*resparam$sim, 
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
      #"regimeProd1"="hmma",
      "regimeProd2"="hmma",
      "regimeProd5"="hmma", 
      "regimeProd7"="hmma"
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
      #"regimeProd1",
      "regimeProd2",
      "regimeProd5", 
      "regimeProd7"
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
  #stat_summary(aes(color=method,x=(model),y=x),fun.data = give.n, geom = "text", 
  #        vjust = -2,position = position_dodge(1))+ 
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




#===========================================
#confusion matrices 




#summarize model selection
aic_set=list()
bic_set=list()
lfo_set=list()

#res_a <- res_a[res_a$convergence==0,]

aic=subset(res_a,parameter=='AIC')
bic=subset(res_a,parameter=='BIC')
lfo=subset(res_a,parameter=='LFO')


aggregate(aic$scenario,list(aic$scenario),length)

lfo<-lfo[lfo$model %in% c("simple","autocorr","rwa","rwb","rwab","hmma","hmmb","hmmab"),]
lfo[is.na(lfo$est),]<--Inf


aic$est[aic$convergence>0]<-Inf
bic$est[aic$convergence>0]<-Inf

scn<-factor(unique(aic$scenario), levels=c(
  "trendLinearProd1", 
  "trendLinearProd2", 
  "trendLinearProd5",
  "trendLinearProd7",
  "regimeProd1",
  "regimeProd2",      
  "regimeProd5",      
  "regimeProd7" ) )


EM=c("stationary",
     "autocorr",
     "dynamic.a","dynamic.b","dynamic.ab",
     "regime.a","regime.b","regime.ab")
##Confusion matrices
conf_matrix<-expand.grid(EM=EM,OM=scn)
conf_matrix$w_AIC=NA
conf_matrix$BIC=NA
conf_matrix$LFO=NA

cn1<-list()
cn2<-list()
cn3<-list()

  aic3<-subset(aic,scenario==scn[2])
  summary(aic3)
  head(tidyr::spread(aic3[,-c(9)],key=model,value=est))

  aic4<-subset(aic,scenario==scn[4])
  summary(aic4)
    head(tidyr::spread(aic4[,-c(9)],key=model,value=est))

o=0
for(a in seq_along(scn)){

  #AIC
  aica<-subset(aic,scenario==scn[a])
  dim(aica)
  unique(aica$scenario)

  #aica[aica$model=="simple",][990:1100,]
  #aicset<-data.frame(simple=aica$est[aica$model=="simple"],
  #                   autocorr=aica$est[aica$model=="autocorr"],
  #                   rwa=aica$est[aica$model=="rwa"],
  #                   rwb=aica$est[aica$model=="rwb"],
  #                   rwab=aica$est[aica$model=="rwab"],
  #                   hmma=aica$est[aica$model=="hmma"],
  #                   hmmb=aica$est[aica$model=="hmmb"],
  #                   hmmab=aica$est[aica$model=="hmmab"])
     



  aic_set[[a]]=tidyr::spread(aica[,-c(9)],key=model,value=est)
  aic_set[[a]]=aic_set[[a]][c(15,8,12,14,13,9,11,10)] #reorder estimation models

  sc1=apply(aic_set[[a]],1,which.min)
  cn1[[a]]=summary(factor(sc1,levels=seq(1:ncol(aic_set[[a]]))))/1000

  bica=subset(bic,scenario==scn[a])
  bic_set[[a]]=tidyr::spread(bica[,-9],key=model,value=est)
  bic_set[[a]]=bic_set[[a]][c(15,8,12,14,13,9,11,10)] #reorder estimation models

  sc2=apply(bic_set[[a]],1,which.min)
  cn2[[a]]=summary(factor(sc2,levels=seq(1:ncol(bic_set[[a]]))))/1000

  lfoa=subset(lfo,scenario==scn[a])
  lfo_set[[a]]=tidyr::spread(lfoa[,-9],key=model,value=est)
  lfo_set[[a]]=lfo_set[[a]][c(15,8,12,14,13,9,11,10)] #reorder estimation models

  sc3=apply(lfo_set[[a]],1,which.max)
  cn3[[a]]=summary(factor(sc3,levels=seq(1:ncol(lfo_set[[a]]))))/1000


  myseq<-seq(from=o+1, length.out=length(EM))
  conf_matrix$w_AIC[myseq]<-cn1[[a]]
  conf_matrix$BIC[myseq]<-cn2[[a]]
  conf_matrix$LFO[myseq]<-cn3[[a]]
  o=max(myseq)

}


conf_matrix$eqem_om <- dplyr::recode(conf_matrix$OM, 
      "stationary"="stationary",
      "autocorr"="autocorr",
      "sigmaShift"="stationary", 
      "decLinearProd"="dynamic.a",
      "sineProd"="dynamic.a",
      "decLinearCap"="dynamic.b",
      "decLinearProdshiftCap"="dynamic.ab",
      "regimeProd"="regime.a",
      "shiftProd"="regime.a",
      "regimeCap"="regime.b",
      "shiftCap"="regime.b", 
      "regimeProdCap"="regime.ab",
      )   
conf_matrix$diag<-conf_matrix$eqem_om==conf_matrix$EM

conf_matrix$EM


mytheme = list(
  theme_classic(14)+
    theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
          legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="bold", size=13),
          axis.text=element_text(face="bold"),axis.title = element_text(face="bold",size=16),plot.title = element_text(face = "bold", hjust = 0.5,size=16))
)








#routines to visualize smulation estimation Results
#Catarina Wor
#November 2022
#============================================
library(ggplot2)
library(gridExtra)
source("code/utils.R")
#read in data
simPar <- read.csv("data/HarCk/harcnkSimPars.csv")

## Store relevant object names to help run simulation 
scenNames <- unique(simPar$scenario)

load("outs/simest/simest_prodcapscenarios1_4.Rdata") 
allsimest14<-allsimest
allrmse14<-allrmse
load("outs/simest/simest_prodcapscenarios5_11.Rdata") 
allsimest511<-allsimest
allrmse511<-allrmse
#combine lists
allsimest<-list()
allrmse<-list()
for(i in seq_along(scenNames)){
    if(i<5){
      allsimest[[i]]<-allsimest14[[i]]
      allrmse[[i]]<-allrmse14[[i]]
    }else{
      allsimest[[i]]<-allsimest511[[i]]
      allrmse[[i]]<-allrmse511[[i]]
    }
}

#=================================
#plots


pbiasplot<-list()
rmseplot<-list()
for(a in seq_along(scenNames)){
  
#pbias plot
#todo - summarize pbias for each iteration 
  scna<-allsimest[[a]]
  dfpbias<-do.call("rbind",scna)

  dfpbias<- dfpbias[dfpbias$convergence==0,]
  #mean by iteration
  dfmpbias<-aggregate(dfpbias$pbias, 
    list(model=dfpbias$model,
         method=dfpbias$method,
         iteration=dfpbias$iteration,
         parameter=dfpbias$parameter),
    mean)

  dfmpbias$convergence <- aggregate(dfpbias$convergence, 
    list(model=dfpbias$model,
         method=dfpbias$method,
         iteration=dfpbias$iteration,
         parameter=dfpbias$parameter),
    function(x)sum(x,na.rm=T))$x
  

  dfmpbias$model <- factor(dfmpbias$model, 
    levels=c("simple",
             "autocorr",
             "rwa",
             "rwb",
             "rwab",
             "hmma_regime",
             "hmmb_regime",
             "hmmab_regime"))
  
  dfmpbias$method <- factor(dfmpbias$method, levels=c("MLE","MCMC"))
  dfmpbias <- dfmpbias[!is.na(dfmpbias$x),] 
  dfmpbias$vjust <- -2
  dfmpbias$vjust <- -2
   
  pbiasplot[[a]] <- ggplot(dfmpbias,aes(x=model,y=x)) +
      geom_boxplot(aes(fill=method),outlier.shape = NA) +
      coord_cartesian(ylim = c(-100,100)) +
      geom_hline(yintercept=0) +
      theme_bw(14)+ 
      facet_wrap(~parameter, scales="free_y") +
      scale_fill_viridis_d(begin=.3, end=.9) +
      stat_summary(aes(hjust = method),fun.data = give.n, 
                   geom = "text", 
                   vjust = 4,
                   fontface="bold") + 
      labs(title = simPar$nameOM[a]) +
      theme(axis.text.x = element_text(angle = 45,
                                       vjust = 0.5, 
                                       hjust=1))

      

  #rmse
  rmsea<-allrmse[[a]]
  dfrmse<-do.call("rbind",rmsea)

  dfrmse<- dfrmse[dfrmse$convergence==0,]
  

  dfrmse$model <- factor(dfrmse$model, 
                    levels=c("simple",
                             "autocorr",
                             "rwa",
                             "rwb",
                             "rwab",
                             "hmma_regime",
                             "hmmb_regime",
                             "hmmab_regime"))

  #c("simple","autocorr","rwa",
  #"rwb","rwab","hmma_regime","hmma_average","hmmb_regime", "hmmb_average",
  #"hmmab_regime", "hmmab_average",  "hmmabhc_regime","hmmabhc_average" )

   rmseplot[[a]] <- ggplot(dfrmse, aes(x=model,y=x)) +
                    geom_boxplot(aes(fill=method),outlier.shape = NA ) +
                    theme_bw(14)+ 
                    facet_wrap(~parameter, scales="free_y")+
                    scale_fill_viridis_d(begin=.3, end=.9)+
                    ylab("RMSE")+ labs(title = simPar$nameOM[a])+
                    stat_summary( aes(vjust=method,hjust=method),
                      fun.data = give.n, 
                                 geom = "text"
                                )+
                    theme(axis.text.x = element_text(angle = 45, 
                                                     vjust = 0.5, 
                                                     hjust=1))

}

ggsave(
      filename = "outs/SamSimOutputs/plotcheck/pbias.pdf", 
      plot = marrangeGrob(pbiasplot, nrow=1, ncol=1), 
      width = 12, height = 8
    )


ggsave(
      filename = "outs/SamSimOutputs/plotcheck/rmse.pdf", 
      plot = marrangeGrob(rmseplot, nrow=1, ncol=1), 
      width = 12, height = 8
    )



allscnest<-list()
allscnrmse<-list()

for(a in seq_along(scenNames)){
  
#pbias plot
#todo - summarize pbias for each iteration 
  scna<-allsimest[[a]]
  dfpbias<-do.call("rbind",scna)
  dfpbias$scn<-scenNames[a]
  dfpbias<- dfpbias[dfpbias$convergence==0,]
  

  allscnest[[a]]<-dfpbias


  rmsea<-allrmse[[a]]
  dfrmse<-do.call("rbind",rmsea)
  dfrmse$scn<-scenNames[a]
  dfrmse<- dfrmse[dfrmse$convergence==0,]
  
  

  allscnrmse[[a]]<-dfrmse

}



scnest <- do.call("rbind",allscnest)
scnrmse <- do.call("rbind",allscnrmse)
head(scnest)
unique(scnest$model)

pbiasmodel<-list()
dfpbiasmodel<-list()
rmsemodel<-list()
scnest<-scnest[!is.na(scnest$model),]
scnrmse<-scnrmse[!is.na(scnrmse$model),]


for(n in seq_along(unique(scnest$model))){

  df<-scnest[scnest$model==unique(scnest$model)[n],]

  #mean by iteration
  dfmpbias<-aggregate(df$pbias, 
    list(model=df$model,
         method=df$method,
         scenario=df$scn,
         iteration=df$iteration,
         parameter=df$parameter),
    mean)

  dfmpbias$convergence <- aggregate(df$convergence, 
    list(model=df$model,
         method=df$method,
         scenario=df$scn,
         iteration=df$iteration,
         parameter=df$parameter),
    function(x)sum(x,na.rm=T))$x
  
  

  
  dfmpbias$method <- factor(dfmpbias$method, levels=c("MLE","MCMC"))
  dfmpbias <- dfmpbias[!is.na(dfmpbias$x),] 

  pbiasmodel[[n]]<-dfmpbias 
 
  pbiasmodel[[n]] <- ggplot(dfmpbias,aes(x=scenario,y=x)) +
      geom_boxplot(aes(fill=method)) +
      coord_cartesian(ylim = c(-100,100)) +
      geom_hline(yintercept=0) +
      theme_bw(14)+ 
      facet_wrap(~parameter, scales="free_y") +
      scale_fill_viridis_d(begin=.3, end=.9) +
      stat_summary(fun.data = give.n, 
                   geom = "text", 
                   hjust = 0.5,
                   vjust = -2) + 
      labs(title = unique(scnest$model)[n]) +
      theme(axis.text.x = element_text(angle = 45,
                                       vjust = 0.5, 
                                       hjust=1))

  dfrmse<-scnrmse[scnrmse$model==unique(scnrmse$model)[n],]


  rmsemodel[[n]] <- ggplot(dfrmse, aes(x=scn,y=x)) +
                    geom_boxplot(aes(fill=method)) +
                    theme_bw(14)+ 
                    facet_wrap(~parameter, scales="free_y")+
                    scale_fill_viridis_d(begin=.3, end=.9)+
                    ylab("RMSE")+ 
                    labs(title = unique(scnest$model)[n])+
                    stat_summary(fun.data = give.n, 
                                 geom = "text", 
                                 hjust = 0.5,
                                 vjust = -2)+
                    theme(axis.text.x = element_text(angle = 45, 
                                                     vjust = 0.5, 
                                                     hjust=1))
}


df <- do.call("rbind",allscnest)

dfmpbias<-aggregate(df$pbias, 
    list(model=df$model,
         method=df$method,
         scenario=df$scn,
         iteration=df$iteration,
         parameter=df$parameter),
    mean)

  dfmpbias$convergence <- aggregate(df$convergence, 
    list(model=df$model,
         method=df$method,
         scenario=df$scn,
         iteration=df$iteration,
         parameter=df$parameter),
    function(x)sum(x,na.rm=T))$x
  
  

  
  dfmpbias$method <- factor(dfmpbias$method, levels=c("MLE","MCMC"))
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
      "regimeProdCap"="hmm",
      "decLinearProdshiftCap"="rwab")   

dfmpbias$simulated_f<-factor(dfmpbias$simulated, levels=c("simple",
                                            "autocorr",
                                            "rwa",
                                            "rwb",
                                            "rwab", 
                                            "hmma",
                                            "hmmb",
                                            "hmm"))


dfmpbias$model <- dplyr::recode(dfmpbias$model, 
       "simple"="simple", 
       "autocorr"="autocorr",
       "rwa"="rwa",
       "rwb"="rwb",
       "rwab"="rwab",
       "hmma"="hmma_regime",  
       "hmmb"="hmmb_regime",    
       "hmmab"="hmmab_regime")

dfmpbias$model <- factor(dfmpbias$model, levels=c(
       "simple", 
       "autocorr",
       "rwa",
       "rwb",
       "rwab",
       "hmma_regime",  
        "hmmb_regime",    
        "hmmab_regime") )

     
dfmpbias$scenario_f<- factor(dfmpbias$scenario, levels=c(
      "stationary",
      "autocorr",
      "sigmaShift", 
      "decLinearProd",
      "sineProd",
      "regimeProd",
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
       "hmma_regime"="hmm",  
       "hmmb_regime"="hmm",    
       "hmmab_regime"="hmm")



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
       "hmma_regime"="tva",  
       "hmmb_regime"="tvb",    
       "hmmab_regime"="tvab")



dfmpbias$simulated_paragg <- dplyr::recode(dfmpbias$simulated, 
       "simple"="simple", 
       "autocorr"="simple",
       "rwa"="tva",
       "rwb"="tvb",
       "rwab"="tvab",
       "hmma"="tva",  
       "hmmb"="tvb",    
       "hmmab"="tvab")



head(dfmpbias)

dfmpbias_ss<-dfmpbias[dfmpbias$scenario_f%in%c("stationary",
      "decLinearProd",
      "shiftCap"),]




mytheme = list(
    theme_classic(16)+
        theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
              legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="bold", size=12),
              axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
)

dfmpbias15<-dfmpbias[dfmpbias$scenario_f%in%c("stationary",
      "autocorr",
      "sigmaShift", 
      "decLinearProd",
      "sineProd"),]

dfmpbias15<-dfmpbias15[dfmpbias15$parameter%in%c("alpha",
      "Smax",
      "smsy"),]

dfscen<-aggregate(dfmpbias15$simulated_f,list(parameter=dfmpbias15$parameter,
  scenario_f=dfmpbias15$scenario_f, model=dfmpbias15$model),unique)
dfscen$simulated_f<-as.numeric(dfscen$x)


#kitchen sink plot
ggplot(dfmpbias15) +   
  geom_boxplot(data=dfmpbias15,aes(fill=method,x=model,y=x), outlier.shape = NA) +         
  coord_cartesian(ylim = quantile(dfmpbias15$x, c(0.1, 0.9)))+
  geom_hline(yintercept=0) +
  mytheme+ 
  ylab("mean % bias") +
  facet_grid(parameter~scenario_f, scales="free_y")+
  scale_fill_viridis_d(begin=.3, end=.8) +
  scale_color_viridis_d(begin=.3, end=.8) +
  geom_rect(data=dfscen,aes(xmin=as.numeric(simulated_f)-.5,
                         xmax=as.numeric(simulated_f)+.5,
                         ymin=-Inf,ymax=Inf),
                    color="gray90",alpha=0.02)+
  #stat_summary(aes(color=method,x=(model),y=x),fun.data = give.n, geom = "text", 
  #        vjust = -2,position = position_dodge(1))+ 
  theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1))


#aggregate plot
dfscen_agg<-aggregate(dfmpbias15$simulated_agg,list(parameter=dfmpbias15$parameter,
  scenario_f=dfmpbias15$scenario_f, model=dfmpbias15$model_par),unique)
dfscen_agg$x<-factor(dfscen_agg$x, levels=c("simple","rw","hmm"))
dfscen_agg$simulated_f<-as.numeric(dfscen_agg$x)


ggplot(dfmpbias15) +   
  geom_boxplot(data=dfmpbias15,aes(fill=method,x=model_agg,y=x), outlier.shape = NA) +         
  coord_cartesian(ylim = c(-90,90))+
  #coord_cartesian(ylim = quantile(dfmpbias15$x, c(0.1, 0.9)))+
  geom_hline(yintercept=0) +
  mytheme+ 
  ylab("mean % bias") +
  facet_grid(parameter~scenario_f, scales="free_y")+
  scale_fill_viridis_d(begin=.3, end=.8) +
  scale_color_viridis_d(begin=.3, end=.8) +
  geom_rect(data=dfscen_agg,aes(xmin=as.numeric(simulated_f)-.5,
                         xmax=as.numeric(simulated_f)+.5,
                         ymin=-Inf,ymax=Inf),
                    color="gray90",alpha=0.02)+
  #stat_summary(aes(color=method,x=(model),y=x),fun.data = give.n, geom = "text", 
  #        vjust = -2,position = position_dodge(1))+ 
  theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1))



summary(dfmpbias15)
dfmpbias15MLE<-dfmpbias15[dfmpbias15$method=="MLE",]


ggplot(dfmpbias15MLE) +   
  geom_violin(trim=FALSE, aes(x=model_agg,y=x),fill="#009194",alpha=.3)+ 
  geom_boxplot(width=0.2,position = position_dodge(1),aes(x=model_agg,y=x),fill="#009194", outlier.shape = NA) +        
  coord_cartesian(ylim = quantile(dfmpbias15$x, c(0.025, 0.975)))+
  geom_hline(yintercept=0) +
  mytheme+ 
  ylab("mean % bias") +
  facet_grid(parameter~scenario_f, scales="free_y")+
  scale_fill_viridis_d(begin=.3, end=.8) +
  scale_color_viridis_d(begin=.3, end=.8) +
  theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1))+
  geom_rect(data=dfscen_agg,aes(xmin=as.numeric(simulated_f)-.5,
                         xmax=as.numeric(simulated_f)+.5,
                         ymin=-Inf,ymax=Inf),
                    color="gray90",alpha=0.02)


dfmpbias15MCMC<-dfmpbias15[dfmpbias15$method=="MCMC",]


ggplot(dfmpbias15MCMC) +   
  geom_violin(trim=FALSE, aes(x=model_agg,y=x),fill="#009194",alpha=.3)+ 
  geom_boxplot(width=0.2,position = position_dodge(1),aes(x=model_agg,y=x),fill="#009194", outlier.shape = NA) +        
  coord_cartesian(ylim = quantile(dfmpbias15$x, c(0.025, 0.975)))+
  geom_hline(yintercept=0) +
  mytheme+ 
  ylab("mean % bias") +
   xlab("estimation model") +
  facet_grid(parameter~scenario_f, scales="free_y")+
  scale_fill_viridis_d(begin=.3, end=.8) +
  scale_color_viridis_d(begin=.3, end=.8) +
  theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1))+
  geom_rect(data=dfscen_agg,aes(xmin=as.numeric(simulated_f)-.5,
                         xmax=as.numeric(simulated_f)+.5,
                         ymin=-Inf,ymax=Inf),
                    color="gray90",alpha=0.05)
  #stat_summary(aes(color=method,x=(model),y=x),fun.data = give.n, geom = "text", 
  #        vjust = -2,position = position_dodge(1))+ 

#aggregate estimates




dfmpbias611<-dfmpbias[dfmpbias$scenario_f%in%c("regimeProd",
      "decLinearCap",
      "regimeCap",
      "shiftCap", 
      "regimeProdCap",
      "decLinearProdshiftCap"),]


dfmpbias611<-dfmpbias611[dfmpbias611$parameter%in%c("alpha",
      "Smax",
      "smsy"),]

dfscen611<-aggregate(dfmpbias611$simulated_f,list(parameter=dfmpbias611$parameter,
  scenario_f=dfmpbias611$scenario_f, model=dfmpbias611$model),unique)
dfscen611$simulated_f<-as.numeric(dfscen611$x)

ggplot(dfmpbias611) +   
  geom_boxplot(data=dfmpbias611,aes(fill=method,x=model,y=x)) +         
  coord_cartesian(ylim = c(-100,100))+
  geom_hline(yintercept=0) +
  mytheme+ 
  ylab("mean % bias") +
  facet_grid(parameter~scenario_f, scales="free_y")+
  scale_fill_viridis_d(begin=.3, end=.8) +
  scale_color_viridis_d(begin=.3, end=.8) +
  geom_rect(data=dfscen611,aes(xmin=as.numeric(simulated_f)-.5,
                         xmax=as.numeric(simulated_f)+.5,
                         ymin=-Inf,ymax=Inf),
                    color="gray90",alpha=0.02)+
  stat_summary(aes(color=method,x=(model),y=x),fun.data = give.n, geom = "text", 
          vjust = -2,position = position_dodge(1))+ 
  theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1))


ggplot(dfmpbias611,aes(x=model,y=x)) +
      geom_boxplot(aes(fill=method)) +
      coord_cartesian(ylim = c(-100,100))+
      geom_hline(yintercept=0) +
      geom_rect(aes(xmin=as.numeric(simulated_f)-.5,
                         xmax=as.numeric(simulated_f)+.5,
                         ymin=-Inf,ymax=Inf),
                    color="gray90",alpha=0.002)+
      mytheme+ 
      ylab("mean % bias") +
      facet_grid(parameter~scenario_f, scales="free_y")+
      scale_fill_viridis_d(begin=.3, end=.8) +
      scale_color_viridis_d(begin=.3, end=.8) +
      stat_summary(aes(color=method),fun.data = give.n, geom = "text", 
          vjust = -2, position = position_dodge(1))+ 
      theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1))




#aggregate plot
dfscen_agg<-aggregate(dfmpbias611$simulated_agg,list(parameter=dfmpbias611$parameter,
  scenario_f=dfmpbias611$scenario_f, model=dfmpbias611$model_par),unique)
dfscen_agg$x<-factor(dfscen_agg$x, levels=c("simple","rw","hmm"))
dfscen_agg$simulated_f<-as.numeric(dfscen_agg$x)


ggplot(dfmpbias611) +   
  geom_boxplot(data=dfmpbias611,aes(fill=method,x=model_agg,y=x), outlier.shape = NA) +         
  coord_cartesian(ylim = c(-90,90))+
  #coord_cartesian(ylim = quantile(dfmpbias15$x, c(0.1, 0.9)))+
  geom_hline(yintercept=0) +
  mytheme+ 
  ylab("mean % bias") +
  facet_grid(parameter~scenario_f, scales="free_y")+
  scale_fill_viridis_d(begin=.3, end=.8) +
  scale_color_viridis_d(begin=.3, end=.8) +
  geom_rect(data=dfscen_agg,aes(xmin=as.numeric(simulated_f)-.5,
                         xmax=as.numeric(simulated_f)+.5,
                         ymin=-Inf,ymax=Inf),
                    color="gray90",alpha=0.02)+
  #stat_summary(aes(color=method,x=(model),y=x),fun.data = give.n, geom = "text", 
  #        vjust = -2,position = position_dodge(1))+ 
  theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1))



#try violin plot
ggplot(dfmpbias611,aes(x=model,y=x)) +
       geom_violin(aes(color=method),trim=FALSE)+
      #stat_summary(fun.data="mean_sdl", mult=1, 
      #           geom="crossbar", width=0.2 )+
      geom_boxplot(aes(fill=method),width=0.2,position = position_dodge(1)) +
      coord_cartesian(ylim = c(-100,100))+
      geom_hline(yintercept=0) +
      geom_rect(aes(xmin=as.numeric(simulated_f)-.5,
                         xmax=as.numeric(simulated_f)+.5,
                         ymin=-Inf,ymax=Inf),
                    color="gray90",alpha=0.002)+
      mytheme+ 
      ylab("mean % bias") +
      facet_grid(parameter~scenario_f, scales="free_y")+
      scale_fill_viridis_d(begin=.3, end=.9) +
      scale_color_viridis_d(begin=.3, end=.9) +
      stat_summary(aes(color=method),fun.data = give.n, geom = "text", 
          vjust = -2,position = position_dodge(1))+ 
      theme(axis.text.x = element_text(angle = 45, vjust=0.5))




dfmpbias611MCMC<-dfmpbias611[dfmpbias611$method=="MCMC",]


ggplot(dfmpbias611MCMC) +   
  geom_violin(trim=FALSE, aes(x=model_agg,y=x),fill="#009194",alpha=.3)+ 
  geom_boxplot(width=0.2,position = position_dodge(1),aes(x=model_agg,y=x),fill="#009194", outlier.shape = NA) +        
  coord_cartesian(ylim = quantile(dfmpbias15$x, c(0.025, 0.975)))+
  geom_hline(yintercept=0) +
  mytheme+ 
  ylab("mean % bias") +
  xlab("estimation model") +
  facet_grid(parameter~scenario_f, scales="free_y")+
  scale_fill_viridis_d(begin=.3, end=.8) +
  scale_color_viridis_d(begin=.3, end=.8) +
  theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1))+
  geom_rect(data=dfscen_agg,aes(xmin=as.numeric(simulated_f)-.5,
                         xmax=as.numeric(simulated_f)+.5,
                         ymin=-Inf,ymax=Inf),
                    color="gray90",alpha=0.05)



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





dfmpbias<-aggregate(df$pbias, 
    list(model=df$model,
         method=df$method,
         scenario=df$scn,
         iteration=df$iteration,
         parameter=df$parameter),
    function(x){mean(abs(x))})

  dfmpbias$convergence <- aggregate(df$convergence, 
    list(model=df$model,
         method=df$method,
         scenario=df$scn,
         iteration=df$iteration,
         parameter=df$parameter),
    function(x)sum(x,na.rm=T))$x
  
  

  
  dfmpbias$method <- factor(dfmpbias$method, levels=c("MLE","MCMC"))
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
      "regimeProdCap"="hmm",
      "decLinearProdshiftCap"="rwab")   

dfmpbias$simulated_f<-factor(dfmpbias$simulated, levels=c("simple",
                                            "autocorr",
                                            "rwa",
                                            "rwb",
                                            "rwab", 
                                            "hmma",
                                            "hmmb",
                                            "hmm"))


dfmpbias$model <- dplyr::recode(dfmpbias$model, 
       "simple"="simple", 
       "autocorr"="autocorr",
       "rwa"="rwa",
       "rwb"="rwb",
       "rwab"="rwab",
       "hmma"="hmma_regime",  
       "hmmb"="hmmb_regime",    
       "hmmab"="hmmab_regime")

dfmpbias$model <- factor(dfmpbias$model, levels=c(
       "simple", 
       "autocorr",
       "rwa",
       "rwb",
       "rwab",
       "hmma_regime",  
        "hmmb_regime",    
        "hmmab_regime") )

     
dfmpbias$scenario_f<- factor(dfmpbias$scenario, levels=c(
      "stationary",
      "autocorr",
      "sigmaShift", 
      "decLinearProd",
      "sineProd",
      "regimeProd",
      "decLinearCap",
      "regimeCap",
      "shiftCap", 
      "regimeProdCap",
      "decLinearProdshiftCap"
      )  )



mytheme = list(
    theme_classic(16)+
        theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
              legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="bold", size=12),
              axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
)

dfmpbias15<-dfmpbias[dfmpbias$scenario_f%in%c("stationary",
      "autocorr",
      "sigmaShift", 
      "decLinearProd",
      "sineProd"),]

dfmpbias15<-dfmpbias15[dfmpbias15$parameter%in%c("alpha",
      "Smax",
      "smsy"),]

dfscen<-aggregate(dfmpbias15$simulated_f,list(parameter=dfmpbias15$parameter,
  scenario_f=dfmpbias15$scenario_f, model=dfmpbias15$model),unique)
dfscen$simulated_f<-as.numeric(dfscen$x)

ggplot(dfmpbias15) +   
  geom_boxplot(data=dfmpbias15,aes(fill=method,x=model,y=x)) +         
  coord_cartesian(ylim = c(0,100))+
  geom_hline(yintercept=0) +
  mytheme+ 
  ylab("mean % bias") +
  facet_grid(parameter~scenario_f, scales="free_y")+
  scale_fill_viridis_d(begin=.3, end=.8) +
  scale_color_viridis_d(begin=.3, end=.8) +
  geom_rect(data=dfscen,aes(xmin=as.numeric(simulated_f)-.5,
                         xmax=as.numeric(simulated_f)+.5,
                         ymin=-Inf,ymax=Inf),
                    color="gray90",alpha=0.02)+
  stat_summary(aes(color=method,x=(model),y=x),fun.data = give.n, geom = "text", 
          vjust = -2,position = position_dodge(1))+ 
  theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1))

dfmpbias611<-dfmpbias[dfmpbias$scenario_f%in%c("regimeProd",
      "decLinearCap",
      "regimeCap",
      "shiftCap", 
      "regimeProdCap",
      "decLinearProdshiftCap"),]


dfmpbias611<-dfmpbias611[dfmpbias611$parameter%in%c("alpha",
      "Smax",
      "smsy"),]

dfscen611<-aggregate(dfmpbias611$simulated_f,list(parameter=dfmpbias611$parameter,
  scenario_f=dfmpbias611$scenario_f, model=dfmpbias611$model),unique)
dfscen611$simulated_f<-as.numeric(dfscen611$x)

ggplot(dfmpbias611) +   
  geom_boxplot(data=dfmpbias611,aes(fill=method,x=model,y=x)) +         
  coord_cartesian(ylim = c(-100,100))+
  geom_hline(yintercept=0) +
  mytheme+ 
  ylab("mean % bias") +
  facet_grid(parameter~scenario_f, scales="free_y")+
  scale_fill_viridis_d(begin=.3, end=.8) +
  scale_color_viridis_d(begin=.3, end=.8) +
  geom_rect(data=dfscen611,aes(xmin=as.numeric(simulated_f)-.5,
                         xmax=as.numeric(simulated_f)+.5,
                         ymin=-Inf,ymax=Inf),
                    color="gray90",alpha=0.02)+
  stat_summary(aes(color=method,x=(model),y=x),fun.data = give.n, geom = "text", 
          vjust = -2,position = position_dodge(1))+ 
  theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1))


ggplot(dfmpbias611,aes(x=model,y=x)) +
      geom_boxplot(aes(fill=method)) +
      coord_cartesian(ylim = c(-100,100))+
      geom_hline(yintercept=0) +
      geom_rect(aes(xmin=as.numeric(simulated_f)-.5,
                         xmax=as.numeric(simulated_f)+.5,
                         ymin=-Inf,ymax=Inf),
                    color="gray90",alpha=0.002)+
      mytheme+ 
      ylab("mean % bias") +
      facet_grid(parameter~scenario_f, scales="free_y")+
      scale_fill_viridis_d(begin=.3, end=.8) +
      scale_color_viridis_d(begin=.3, end=.8) +
      stat_summary(aes(color=method),fun.data = give.n, geom = "text", 
          vjust = -2, position = position_dodge(1))+ 
      theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1))

#try violin plot
ggplot(dfmpbias611,aes(x=model,y=x)) +
       geom_violin(aes(color=method),trim=FALSE)+
      #stat_summary(fun.data="mean_sdl", mult=1, 
      #           geom="crossbar", width=0.2 )+
      geom_boxplot(aes(fill=method),width=0.2,position = position_dodge(1)) +
      coord_cartesian(ylim = c(-100,100))+
      geom_hline(yintercept=0) +
      geom_rect(aes(xmin=as.numeric(simulated_f)-.5,
                         xmax=as.numeric(simulated_f)+.5,
                         ymin=-Inf,ymax=Inf),
                    color="gray90",alpha=0.002)+
      mytheme+ 
      ylab("mean % bias") +
      facet_grid(parameter~scenario_f, scales="free_y")+
      scale_fill_viridis_d(begin=.3, end=.9) +
      scale_color_viridis_d(begin=.3, end=.9) +
      stat_summary(aes(color=method),fun.data = give.n, geom = "text", 
          vjust = -2,position = position_dodge(1))+ 
      theme(axis.text.x = element_text(angle = 45, vjust=0.5))





scnrmse <- do.call("rbind",allscnrmse)


ggplot(scnrmse) +   
  geom_boxplot(data=scnrmse,aes(fill=method,x=model,y=x),outlier.shape = NA) +         
  geom_hline(yintercept=0) +
  mytheme+ 
  ylab("mean % bias") +
  facet_grid(parameter~scn, scales="free_y")+
  scale_fill_viridis_d(begin=.3, end=.8) +
  scale_color_viridis_d(begin=.3, end=.8) +
  geom_rect(data=dfscen611,aes(xmin=as.numeric(simulated_f)-.5,
                         xmax=as.numeric(simulated_f)+.5,
                         ymin=-Inf,ymax=Inf),
                    color="gray90",alpha=0.02)+
  stat_summary(aes(color=method,x=(model),y=x),fun.data = give.n, geom = "text", 
          vjust = -2,position = position_dodge(1))+ 
  theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1))



#============================================================
#plot of pbiasER
#read in data
simParER <- read.csv("data/HarCkER/harcnkSimPars_ER.csv")
## Store relevant object names to help run simulation 
scenNamesER <- unique(simParER$scenario)

load("outs/simest/simest_ERcenarios1_3.Rdata") 
allsimestER13<-allsimest
allrmseER13<-allrmse
load("outs/simest/simest_ERcenarios4_11.Rdata") 
allsimestER58<-allsimest
allrmseER58<-allrmse
#combine lists
allsimestER<-list()
allrmseER<-list()
for(i in seq_along(scenNamesER)){
    if(i<4){
      allsimestER[[i]]<-allsimestER13[[i]]
      allrmseER[[i]]<-allrmseER13[[i]]
    }else{
      allsimestER[[i]]<-allsimestER58[[i]]
      allrmseER[[i]]<-allrmseER58[[i]]
    }
}






dfpbiasallER<-list()


for(a in seq_along(scenNamesER)){
  #a<-1
#pbias plot
#todo - summarize pbias for each iteration 
  scna<-allsimestER[[a]]
  dfpbiasER<-do.call("rbind",scna)

  dfpbiasER<- dfpbiasER[dfpbiasER$convergence==0,]
  #mean by iteration
  dfmpbiasER<-aggregate(dfpbiasER$pbias, 
    list(model=dfpbiasER$model,
         method=dfpbiasER$method,
         iteration=dfpbiasER$iteration,
         parameter=dfpbiasER$parameter),
    mean)

  dfmpbiasER$convergence <- aggregate(dfpbiasER$convergence, 
    list(model=dfpbiasER$model,
         method=dfpbiasER$method,
         iteration=dfpbiasER$iteration,
         parameter=dfpbiasER$parameter),
    function(x)sum(x,na.rm=T))$x
  

  dfmpbiasER$model <- factor(dfmpbiasER$model, 
    levels=c("simple",
             "autocorr",
             "rwa",
             "rwb",
             "rwab",
             "hmma_regime",
             "hmmb_regime",
             "hmmab_regime"))
  
  dfmpbiasER$method <- factor(dfmpbiasER$method, levels=c("MLE","MCMC"))
  dfmpbiasER <- dfmpbiasER[!is.na(dfmpbiasER$x),] 
  dfmpbiasER$vjust <- -2
  dfmpbiasER$vjust <- -2
  dfmpbiasER$scenario <-scenNamesER[a]
   
   dfpbiasallER[[a]]<-dfmpbiasER

  #rmse
  rmsea<-allrmseER[[a]]
  dfrmseER<-do.call("rbind",rmsea)

  dfrmseER<- dfrmseER[dfrmseER$convergence==0,]
   dfrmseER$scenario <-scenNamesER[a]
  

  dfrmseER$model <- factor(dfrmseER$model, 
                    levels=c("simple",
                             "autocorr",
                             "rwa",
                             "rwb",
                             "rwab",
                             "hmma_regime",
                             "hmmb_regime",
                             "hmmab_regime"))

 
   
}




dfer<-do.call(rbind,dfpbiasallER)
dfer$simulated <- "simple"  

dfer$simulated_f<-factor(dfer$simulated, levels=c("simple",
                                            "autocorr",
                                            "rwa",
                                            "rwb",
                                            "rwab", 
                                            "hmma",
                                            "hmmb",
                                            "hmm"))



dfer$trend<- "low ER"
dfer$trend[dfer$scenario == "highERLowError"|
dfer$scenario == "highERHighError"] <- "high ER"

dfer$trend[dfer$scenario == "ShiftERLowError"|
dfer$scenario == "ShiftERHighError"] <- "shift ER"

dfer$trend[dfer$scenario == "trendERLowError"|
dfer$scenario == "trendERHighError"] <- "trend ER"

dfer$ERvar<- "High CV"

dfer$ERvar[dfer$scenario == "highERLowError"|
dfer$scenario == "ShiftERLowError"|
dfer$scenario == "trendERLowError"|
dfer$scenario == "lowERLowError"] <- "Low CV"



summary(dfer)

dfera<-dfer[dfer$parameter=="alpha"|
dfer$parameter=="Smax"|
dfer$parameter=="smsy",]



mytheme = list(
    theme_classic(16)+
        theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
              legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="bold", size=12),
              axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
)

pbiasparamER<-ggplot(dfera,aes(x=model,y=x)) +
      geom_boxplot(aes(fill=method)) +
      coord_cartesian(ylim = c(-100,100))+
      geom_hline(yintercept=0) +
      geom_rect(aes(xmin=as.numeric(simulated_f)-.5,
                         xmax=as.numeric(simulated_f)+.5,
                         ymin=-Inf,ymax=Inf),
                    color="gray90",alpha=0.002)+
      mytheme+ 
      facet_grid(parameter~scenario, scales="free_y")+
      scale_fill_viridis_d(begin=.3, end=.9) +
      stat_summary(fun.data = give.n, geom = "text", 
          vjust = -2)+ 
      theme(axis.text.x = element_text(angle = 45, vjust=0.5))

pbiasparamER

ggsave(
      filename = "outs/SamSimOutputs/plotcheck/pbias_all_erscn.png", 
      plot = pbiasparamER, 
      width = 12, height = 7
    )



#============================================================
#plot of prior distributions


#dgamma(3,1)

al<-seq(0,10,by=.1)
be<-seq(-25,-1,by=.1)
si<-seq(0,3,by=.01)

par(mfrow=c(1,3),
  font.lab=2,
  cex.axis=1.5,
  las=1,
  cex.main=2)
plot(al,dgamma(al,3,1),type="l", lwd=2,ylab="",
  main=expression(alpha), 
  xlab="alpha",
  font=2
  )
plot(exp(be),dnorm(be,-12,3),type="l", lwd=2,ylab="",
  main=expression(beta), 
  xlab="beta", 
  xlim=c(1e-10,1e-1),
  log="x",
  font=2)
plot(si,dgamma(si,2,3),type="l", lwd=2,ylab="",
  main=expression(sigma),
  xlab="sigma",
  font=2)







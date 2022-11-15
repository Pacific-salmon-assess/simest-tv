#routines to visualize smulation estimation Results
#Catarina Wor
#November 2022
#============================================

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
  
  #Remove hmm average
  dfmpbias <- dfmpbias[!(dfmpbias$model %in% 
    c("hmma_average",
      "hmmb_average",
      "hmmab_average",
      "hmmabhc_average" )),]

  dfmpbias$model <- factor(dfmpbias$model, 
    levels=c("simple",
             "autocorr",
             "rwa",
             "rwb",
             "rwab",
             "hmma_regime",
             "hmmb_regime",
             "hmmab_regime",
             "hmmabhc_regime"))
  
  dfmpbias$method <- factor(dfmpbias$method, levels=c("MLE","MCMC"))
  dfmpbias <- dfmpbias[!is.na(dfmpbias$x),] 
 
  pbiasplot[[a]] <- ggplot(dfmpbias,aes(x=model,y=x)) +
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
      labs(title = simPar$nameOM[a]) +
      theme(axis.text.x = element_text(angle = 45,
                                       vjust = 0.5, 
                                       hjust=1))

  #rmse
  rmsea<-allrmse[[a]]
  dfrmse<-do.call("rbind",rmsea)

  dfrmse<- dfrmse[dfrmse$convergence==0,]
  
  dfrmse <- dfrmse[!(dfrmse$model %in% 
    c("hmma_average",
      "hmmb_average",
      "hmmab_average",
      "hmmabhc_average")),]

  dfrmse$model <- factor(dfrmse$model, 
                    levels=c("simple",
                             "autocorr",
                             "rwa",
                             "rwb",
                             "rwab",
                             "hmma_regime",
                             "hmmb_regime",
                             "hmmab_regime",
                             "hmmabhc_regime"))

  #c("simple","autocorr","rwa",
  #"rwb","rwab","hmma_regime","hmma_average","hmmb_regime", "hmmb_average",
  #"hmmab_regime", "hmmab_average",  "hmmabhc_regime","hmmabhc_average" )

   rmseplot[[a]] <- ggplot(dfrmse, aes(x=model,y=x)) +
                    geom_boxplot(aes(fill=method)) +
                    theme_bw(14)+ 
                    facet_wrap(~parameter, scales="free_y")+
                    scale_fill_viridis_d(begin=.3, end=.9)+
                    ylab("RMSE")+ labs(title = simPar$nameOM[a])+
                    stat_summary(fun.data = give.n, 
                                 geom = "text", 
                                 hjust = 0.5,
                                 vjust = -2)+
                    theme(axis.text.x = element_text(angle = 45, 
                                                     vjust = 0.5, 
                                                     hjust=1))

}

ggsave(
      filename = "outs/SamSimOutputs/plotcheck/pbias.pdf", 
      plot = marrangeGrob(pbiasplot, nrow=1, ncol=1), 
      width = 12, height = 5
    )


ggsave(
      filename = "outs/SamSimOutputs/plotcheck/rmse.pdf", 
      plot = marrangeGrob(rmseplot, nrow=1, ncol=1), 
      width = 12, height = 8
    )




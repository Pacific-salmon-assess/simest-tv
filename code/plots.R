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
  dfmpbias$vjust <- -2
  dfmpbias$vjust[] <- -2
   
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
                    geom_boxplot(aes(fill=method),outlier.shape = NA ) +
                    theme_bw(14)+ 
                    facet_wrap(~parameter, scales="free_y")+
                    scale_fill_viridis_d(begin=.3, end=.9)+
                    ylab("RMSE")+ labs(title = simPar$nameOM[a])+
                    stat_summary( aes(fill=method)
                      fun.data = give.n, 
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
  #Remove hmm average
  dfpbias <- dfpbias[!(dfpbias$model %in% 
    c("hmma_average",
      "hmmb_average",
      "hmmab_average",
      "hmmabhc_average" )),]

  allscnest[[a]]<-dfpbias


  rmsea<-allrmse[[a]]
  dfrmse<-do.call("rbind",rmsea)
  dfrmse$scn<-scenNames[a]
  dfrmse<- dfrmse[dfrmse$convergence==0,]
  
  dfrmse <- dfrmse[!(dfrmse$model %in% 
    c("hmma_average",
      "hmmb_average",
      "hmmab_average",
      "hmmabhc_average")),]

  allscnrmse[[a]]<-dfrmse

}



scnest <- do.call("rbind",allscnest)
scnrmse <- do.call("rbind",allscnrmse)
head(scnest)
unique(scnest$model)

pbiasmodel<-list()
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







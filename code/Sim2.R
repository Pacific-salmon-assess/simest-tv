#=============================================================
#Run samSim simulations 
#using the harrison chinook data as example
#samSim updates
#Catarina Wor
#October 2022 
#=============================================================

#install samsim 
#remotes::install_github("Pacific-salmon-assess/samSim", ref="timevar", force=TRUE)

#install samest
#remotes::install_git('https://github.com/Pacific-salmon-assess/samEst')

library(samSim)
library(ggplot2)
library(devtools)
library(gridExtra)
library(dplyr)
library(here)
source("code/utils.R")

## Load relevant input data
# Simulation run parameters describing different scenarios
simPar <- read.csv("data/generic/SimPars.csv")
# CU-specific parameters
cuPar <- read.csv("data/generic/CUPars.csv")

## Store relevant object names to help run simulation 
scenNames <- unique(simPar$scenario)

## First check to ensure that a single scenario can be run (only a small number
# of trials necessary)

#Run and save simulated data

for(a in seq_len(nrow(simPar))){
for(a in 5:nrow(simPar)){
   print(a)
   genericRecoverySim(simPar=simPar[a,], 
                      cuPar=cuPar, 
                      catchDat=NULL, 
                      srDat=NULL,
                      variableCU=FALSE, 
                      ricPars=NULL, 
                      larkPars=NULL, 
                      cuCustomCorrMat= NULL,
                      outDir="outs", 
                      nTrials=10000, 
                      makeSubDirs=TRUE, 
                      random=FALSE, 
                      uniqueProd=TRUE,
                      uniqueSurv=FALSE)

}




simData<-list()
srplot<-list()
recplot<-list()
spnplot<-list()
erplot<-list()
paramplot<-list()
stackcu<-list()



for(a in seq_len(nrow(simPar))){

  simData[[a]] <- readRDS(paste0("outs/SamSimOutputs/", 
                          simPar$nameOM[a],"/",
                          simPar$scenario[a],"/",
                          paste(simPar$nameOM[a],"_", 
                          simPar$nameMP[a], "_", 
                          "CUsrDat.RData",sep="")))$srDatout
  


  dat<-simData[[a]] 
  dat<-dat[dat$year>(max(dat$year)-46),]
  dat <- dat[!is.na(dat$obsRecruits),]
  


  stackcu[[a]]<-cbind(dat[,-c(9,10,11)],stack(dat[,9:11]))
  stackcu[[a]]$scenario<-simPar$nameOM[a]
  stackcu1<-stackcu[[a]]


  paramplot[[a]] <- ggplot(stackcu1) +
    geom_line(aes(x=year,y=values, col=as.factor(CU))) +
    geom_point(aes(x=year,y=values, col=as.factor(CU))) +
    theme_classic(14) + theme(legend.position="none") +
    facet_wrap(~ind, scales="free_y") +
    scale_colour_viridis_d() +
    labs(title = simPar$nameOM[a])

    
  S <- seq(0,max(dat$spawners),by=1000)
  R <- matrix(NA, ncol=length(unique(dat$year)),nrow=length(S))
  
  for(i in seq_along(unique(dat$year))){
    alpha<- mean(dat$alpha[dat$year==unique(dat$year)[i]])
    beta<- unique(dat$beta[dat$year==unique(dat$year)[i]])
    R[,i]<-S*exp(alpha-beta*S)
  }
    
  actualSR<-data.frame(year=rep(unique(dat$year),
      each=length(S)),
      spawners=S,
      recruits=c(R))
  
  srplot[[a]]<-ggplot( actualSR) +
    geom_line(aes(x=spawners,y=recruits, col=as.factor(year)),linewidth=2) +
    theme_bw(14) + 
    scale_colour_viridis_d(end=.7) +
    geom_point(data=dat,aes(x=spawners,y=recruits),alpha=.5) +
    labs(title = simPar$nameOM[a])

  recplot[[a]]<-ggplot(dat) +
    geom_boxplot(aes(x=as.factor(year),y=recruits),alpha=.5) +
    theme_bw(14)+ xlab("year") +
    labs(title = simPar$nameOM[a])
   
  spnplot[[a]]<-ggplot(dat) +
      geom_boxplot(aes(x=as.factor(year),y=spawners),alpha=.5) +
      geom_hline(data=dat,aes(yintercept=sMSY, col=as.factor(year))) +
      theme_bw(14) + 
      xlab("year") + 
      scale_color_discrete(name = "Smsy") +
      labs(title = simPar$nameOM[a])

  erplot[[a]]<-ggplot(dat) +
      geom_boxplot(aes(x=as.factor(year),y=ER),alpha=.5) +
      geom_hline(data=dat,aes(yintercept=uMSY, col=as.factor(year))) +
      theme_bw(14) + 
      xlab("year") + 
      scale_color_discrete(name = "Umsy") +
      labs(title = simPar$nameOM[a])
    
  }





ggsave(
      filename = "outs/SamSimOutputs/plotcheck/paramplots.pdf", 
      plot = marrangeGrob(paramplot, nrow=1, ncol=1), 
      width = 12, height = 5
    )

ggsave(
      filename = "outs/SamSimOutputs/plotcheck/srplots.pdf", 
      plot = marrangeGrob(srplot, nrow=1, ncol=1), 
      width = 12, height = 5
    )
  
ggsave(
      filename = "outs/SamSimOutputs/plotcheck/recplots.pdf", 
      plot = marrangeGrob(recplot, nrow=1, ncol=1), 
      width = 12, height = 5
    )

ggsave(
      filename = "outs/SamSimOutputs/plotcheck/spnplots.pdf", 
      plot = marrangeGrob(spnplot, nrow=1, ncol=1), 
      width = 12, height = 5
    )

ggsave(
      filename = "outs/SamSimOutputs/plotcheck/erplots.pdf", 
      plot = marrangeGrob(erplot, nrow=1, ncol=1), 
      width = 12, height = 5
    )

#==========================




#Scenario plots:
allscnsim <- do.call(rbind,stackcu)
#one example per scearnio
simpr=subset(allscnsim, iteration==1)
simpr$year=simpr$year-54 #just rescaling to years observed in the time-series, easier to interpret

summary(simpr)

unique(simpr$scenario)


simpr$scenario[simpr$scenario== "decLinearProdshiftCap"] <- "decLinProd shiftCap"
simpr$scenario[simpr$scenario== "decLinearProd"] <- "decLinProd"
simpr$scenario[simpr$scenario== "decLinearCap"] <- "decLinCap"

simpr <- simpr[simpr$scenario!="sigmaShift",]
simpr <- simpr[simpr$scenario!="sineProd",] 
simpr <- simpr[simpr$scenario!="shiftProd",] 
simpr <- simpr[simpr$scenario!="shiftCap",] 
#simpr <- simpr[simpr$scenario!="autocorr",] 
simpr <- simpr[simpr$ind!="sigma",] 

simpr$scenario_f = factor(simpr$scenario, 
  levels=c("stationary",
           "autocorr",
           "decLinProd", 
           "decLinCap",
           "decLinProd shiftCap",
           "regimeProd",
           "regimeCap",
           "regimeProdCap"
           ))

simpr$tipo <- "dynamic"
simpr$tipo[simpr$scenario == "regimeProd"|
simpr$scenario == "regimeCap"|
simpr$scenario == "regimeProdCap"] <- "regime"


simpr$tipo[simpr$scenario == "stationary" |
             simpr$scenario == "autocorr"] <- "astat"


 
simpr

mytheme = list(
    theme_classic(14)+
        theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
              legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="bold", size=12),
              axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
)



p <- ggplot(simpr) +
     geom_line(aes(x=year,y=values,col=tipo)) +
     geom_point(aes(x=year,y=values,col=tipo)) +
     mytheme + theme(legend.position="none") +
     facet_grid(ind~scenario_f, scales="free_y",
               labeller = labeller(scenario_f = label_wrap_gen(10))) +
     scale_colour_brewer(palette='Blues') 
p
ggsave(filename = "outs/par_sims.pdf",
        plot=p,
        width=12,height=6)

ggsave(
      filename = "outs/SamSimOutputs/plotcheck/erplots.pdf", 
      plot = marrangeGrob(erplot, nrow=1, ncol=1), 
      width = 12, height = 5
    )     

# add sceatios for changes of age at spawners - later


#Plot one iteration and recruitment curves per scenario



actualSR<-list()
alldat<-list()

for(a in seq_len(nrow(simPar))){

  simData[[a]] <- readRDS(paste0("outs/SamSimOutputs/", 
                          simPar$nameOM[a],"/",
                          simPar$scenario[a],"/",
                          paste(simPar$nameOM[a],"_", 
                          simPar$nameMP[a], "_", 
                          "CUsrDat.RData",sep="")))$srDatout

  dat<-simData[[a]] 
  dat<-dat[dat$year>(max(dat$year)-46),]
  dat <- dat[!is.na(dat$obsRecruits),]
  
  dat <- dat[dat$iteration==sample(unique(dat$iteration),1),]
  dat$scenario <- simPar$scenario[a]
  alldat[[a]]<-dat

  S <- seq(0,600000,by=1000)
  R <- matrix(NA, ncol=length(unique(dat$year)),nrow=length(S))
  
  for(i in unique(dat$year)){

    alpha<- dat$alpha[dat$year==i]
    beta<- dat$beta[dat$year==i]
    R[,which(unique(dat$year)==i)]<-S*exp(alpha-beta*S)
  }
    
  actualSR[[a]]<-data.frame(year=rep(unique(dat$year),
      each=length(S)),
      spawners=S,
      recruits=c(R),
      scenario=simPar$scenario[a])

}

SRdf<-do.call(rbind,actualSR)
datdf<-do.call(rbind,alldat)


scn=c("stationary","autocorr",
             "decLinearProd","decLinearCap",  "decLinearProdshiftCap", "regimeProd",
             "regimeCap",
             "regimeProdCap")
SRdf=subset(SRdf,scenario %in% scn)

SRdf$scenario_f=SRdf$scenario
SRdf$scenario_f=ifelse(SRdf$scenario_f=="decLinearProd",'dec.prod',SRdf$scenario_f)
SRdf$scenario_f=ifelse(SRdf$scenario_f=="decLinearCap",'dec.cap',SRdf$scenario_f)
SRdf$scenario_f=ifelse(SRdf$scenario_f=="decLinearProdshiftCap",'dec.prodcap',SRdf$scenario_f)
SRdf$scenario_f=ifelse(SRdf$scenario_f=="regimeProd",'regime.prod',SRdf$scenario_f)
SRdf$scenario_f=ifelse(SRdf$scenario_f=="regimeCap",'regime.cap',SRdf$scenario_f)
SRdf$scenario_f=ifelse(SRdf$scenario_f=="regimeProdCap",'regime.prodcap',SRdf$scenario_f)

SRdf$scenario_f <-factor(SRdf$scenario_f, levels=c("stationary","autocorr",
                                                        "dec.prod","dec.cap","dec.prodcap","regime.prod",
                                                        "regime.cap",
                                                        "regime.prodcap"))


datdf$scenario_f <-factor(datdf$scenario_f, levels=c("stationary","autocorr",
                                                   "dec.prod","dec.cap","dec.prodcap","regime.prod",
                                                   "regime.cap",
                                                   "regime.prodcap"))




mytheme = list(
    theme_classic(14)+
        theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
              legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="bold", size=12),
              axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
)

SRexample<-  ggplot(SRdf) +
    geom_line(aes(x=spawners,y=recruits, col=as.factor(year)),linewidth=2) +
    mytheme + 
    theme(legend.position="right") +
    scale_colour_viridis_d(end=.85) +
    labs(col = "year") +
    geom_point(data=datdf,aes(x=spawners,y=recruits,col=as.factor(year)),alpha=.5) +
    facet_wrap(~scenario_f)
pdf('sr_examples.pdf',width=12,height=10)
SRexample    
dev.off()

ggsave(
      filename = "outs/SamSimOutputs/plotcheck/srexample.png", 
      plot = SRexample, 
      width = 12, height = 6
    )








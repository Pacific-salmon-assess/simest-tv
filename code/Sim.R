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
simPar <- read.csv("data/harck/harcnkSimPars.csv")
# CU-specific parameters
cuPar <- read.csv("data/harck/harcnkCUPars.csv")

## Store relevant object names to help run simulation 
scenNames <- unique(simPar$scenario)
dirNames <- sapply(scenNames, function(x) paste(x, unique(simPar$species),sep = "_"))

## First check to ensure that a single scenario can be run (only a small number
# of trials necessary)

#Run and save simulated data

for(a in seq_len(nrow(simPar))){
  
   genericRecoverySim(simPar=simPar[a,], cuPar=cuPar, catchDat=NULL, srDat=NULL,
            variableCU=FALSE, ricPars=NULL, larkPars=NULL, cuCustomCorrMat= NULL,
            outDir="outs", nTrials=100, makeSubDirs=TRUE, random=FALSE, uniqueProd=TRUE,
                               uniqueSurv=FALSE)

}




simData<-list()
srplot<-list()
recplot<-list()
spnplot<-list()
erplot<-list()
paramplot<-list()


for(a in seq_len(nrow(simPar))){

  simData[[a]] <- readRDS(paste0("outs/SamSimOutputs/simData/", simPar$nameOM[a],"/",simPar$scenario[a],"/",
                         paste(simPar$nameOM[a],"_", simPar$nameMP[a], "_", "CUsrDat.RData",sep="")))$srDatout
  


  dat<-simData[[a]] 
  dat<-dat[dat$year>(max(dat$year)-46),]
  dat <- dat[!is.na(dat$obsRecruits),]
  


  stackcu1<-cbind(dat[,-c(9,10,11)],stack(dat[,9:11]))

  paramplot[[i]] <- ggplot(stackcu1) +
    geom_line(aes(x=year,y=values, col=as.factor(CU)))+
    geom_point(aes(x=year,y=values, col=as.factor(CU)))+
    theme_bw(14)+ theme(legend.position="none")+
    facet_wrap(~ind, scales="free_y")+
    scale_colour_viridis_d() +
    labs(title = simPar$nameOM[i])

    
  S<-seq(0,max(dat$spawners),by=1000)
  R<-matrix(NA, ncol=length(unique(dat$year)),nrow=length(S))
  for(i in seq_along(unique(dat$year))){
    alpha<- mean(dat$alpha[dat$year==unique(dat$year)[i]])

    beta<- unique(dat$beta[dat$year==unique(dat$year)[i]])
    R[,i]<-S*exp(alpha-beta*S)
  }
    
  actualSR<-data.frame(year=rep(unique(dat$year),each=length(S)),
      spawners=S,
      recruits=c(R))
  
 

  srplot[[a]]<-ggplot( actualSR) +
    geom_line(aes(x=spawners,y=recruits, col=as.factor(year)),size=2)+
    theme_bw(14) + scale_colour_viridis_d(end=.7)+
    geom_point(data=dat,aes(x=spawners,y=recruits),alpha=.5)+
    labs(title = simPar$nameOM[a])

  
  
   

  recplot[[a]]<-ggplot(dat) +
    geom_boxplot(aes(x=as.factor(year),y=recruits),alpha=.5)+
    theme_bw(14)+ xlab("year")+
      labs(title = simPar$nameOM[a])

   
  spnplot[[a]]<-ggplot(dat) +
      geom_boxplot(aes(x=as.factor(year),y=spawners),alpha=.5)+
      geom_hline(data=dat,aes(yintercept=sMSY, col=as.factor(year)))+
      theme_bw(14)+ xlab("year")+ scale_color_discrete(name = "Smsy")+
      labs(title = simPar$nameOM[a])

  erplot[[a]]<-ggplot(dat) +
      geom_boxplot(aes(x=as.factor(year),y=ER),alpha=.5)+
      geom_hline(data=dat,aes(yintercept=uMSY, col=as.factor(year)))+
      theme_bw(14)+ xlab("year")+ scale_color_discrete(name = "Umsy")+
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
#todo





# add sceatios for changes of age at spawners - later



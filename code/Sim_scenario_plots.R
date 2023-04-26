#===============================================
#Ricker curve scenario plots
# april 2023
#===============================================
library(ggplot2)
library(dplyr)




mytheme = list(
    theme_classic(14)+
        theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
              legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="bold", size=12),
              axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
)


#----------------------------------------
#base case                              |
#----------------------------------------
simPar <- read.csv("data/generic/SimPars.csv")

simData<-list()
actualSR<-list()
alldat<-list()

for(a in seq_len(nrow(simPar))){

  simData[[a]] <- readRDS(paste0("outs/SamSimOutputs/simData/", 
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

SRdf$scenario_f <-factor(SRdf$scenario, levels=c("stationary", "autocorr","sigmaShift",
                                                "decLinearProd", "regimeProd", "sineProd", "shiftProd",
                                                "regimeCap", "decLinearCap", "shiftCap", 
                                                "regimeProdCap", "decLinearProdshiftCap"))   



datdf$scenario_f <-factor(datdf$scenario, levels=c("stationary","autocorr","sigmaShift",
              "decLinearProd", "regimeProd", "sineProd","shiftProd",
               "regimeCap", "decLinearCap", "shiftCap",
               "regimeProdCap",  "decLinearProdshiftCap"))




SRexample<-  ggplot(SRdf) +
    geom_line(aes(x=spawners,y=recruits, col=as.factor(year)),linewidth=2) +
    mytheme + 
    theme(legend.position="right") +
    scale_colour_viridis_d(end=.85) +
    labs(col = "year") +
    geom_point(data=datdf,aes(x=spawners,y=recruits,col=as.factor(year)),alpha=.5) +
    facet_wrap(~scenario_f)
SRexample    
 
ggsave(
      filename = "outs/SamSimOutputs/plotcheck/srcurve_basecase.png", 
      plot = SRexample, 
      width = 12, height = 6
    )


#----------------------------------------
#sensitivity alpha case                              |
#----------------------------------------

simPar_sena <- read.csv("data/sensitivity/SimPars.csv")


simData_sena<-list()
actualSR_sena<-list()
alldat_sena<-list()

for(a in seq_len(nrow(simPar_sena))){

  simData_sena[[a]] <- readRDS(paste0("outs/SamSimOutputs/simData/", 
                          simPar_sena$nameOM[a],"/",
                          simPar_sena$scenario[a],"/",
                          paste(simPar_sena$nameOM[a],"_", 
                          simPar_sena$nameMP[a], "_", 
                          "CUsrDat.RData",sep="")))$srDatout

  dat<-simData_sena[[a]] 
  dat<-dat[dat$year>(max(dat$year)-46),]
  dat <- dat[!is.na(dat$obsRecruits),]
  
  dat <- dat[dat$iteration==sample(unique(dat$iteration),1),]
  dat$scenario <- simPar_sena$scenario[a]
  alldat_sena[[a]]<-dat

  S <- seq(0,600000,by=1000)
  R <- matrix(NA, ncol=length(unique(dat$year)),nrow=length(S))
  
  for(i in unique(dat$year)){

    alpha<- dat$alpha[dat$year==i]
    beta<- dat$beta[dat$year==i]
    R[,which(unique(dat$year)==i)]<-S*exp(alpha-beta*S)
  }
    
  actualSR_sena[[a]]<-data.frame(year=rep(unique(dat$year),
      each=length(S)),
      spawners=S,
      recruits=c(R),
      scenario=simPar_sena$scenario[a])

}

SRdf_sena<-do.call(rbind,actualSR_sena)
datdf_sena<-do.call(rbind,alldat_sena)
unique(SRdf_sena$scenario)
SRdf_sena$scenario_f <-factor(SRdf_sena$scenario, levels=c("trendLinearProd1", "trendLinearProd2",
                                                           "trendLinearProd5", "trendLinearProd7",
                                                            "regimeProd1",      "regimeProd2",
                                                            "regimeProd5",      "regimeProd7"))   



datdf_sena$scenario_f <-factor(datdf_sena$scenario, levels=c("trendLinearProd1", "trendLinearProd2",
                                                           "trendLinearProd5", "trendLinearProd7",
                                                            "regimeProd1",      "regimeProd2",
                                                            "regimeProd5",      "regimeProd7"))




SRexample_sena<-  ggplot(SRdf_sena) +
    geom_line(aes(x=spawners,y=recruits, col=as.factor(year)),linewidth=2) +
    mytheme + 
    theme(legend.position="right") +
    scale_colour_viridis_d(end=.85) +
    labs(col = "year") +
    geom_point(data=datdf_sena,aes(x=spawners,y=recruits,col=as.factor(year)),alpha=.5) +
    facet_wrap(~scenario_f)
SRexample_sena    
 
ggsave(
      filename = "outs/SamSimOutputs/plotcheck/srcurve_sensitivitya.png", 
      plot = SRexample_sena, 
      width = 12, height = 6
    )

#----------------------------------------
#sensitivity alpha - half Smax case                              |
#----------------------------------------


simPar_sena_hsmax <- read.csv("data/sensitivity_halfSmax/SimPars.csv")


simData_sena_hsmax<-list()
actualSR_sena_hsmax<-list()
alldat_sena_hsmax<-list()

for(a in seq_len(nrow(simPar_sena_hsmax))){

  simData_sena_hsmax[[a]] <- readRDS(paste0("outs/SamSimOutputs/simData/", 
                          simPar_sena_hsmax$nameOM[a],"/",
                          simPar_sena_hsmax$scenario[a],"/",
                          paste(simPar_sena_hsmax$nameOM[a],"_", 
                          simPar_sena_hsmax$nameMP[a], "_", 
                          "CUsrDat.RData",sep="")))$srDatout

  dat<-simData_sena_hsmax[[a]] 
  dat<-dat[dat$year>(max(dat$year)-46),]
  dat <- dat[!is.na(dat$obsRecruits),]
  
  dat <- dat[dat$iteration==sample(unique(dat$iteration),1),]
  dat$scenario <- simPar_sena_hsmax$scenario[a]
  alldat_sena_hsmax[[a]]<-dat

  S <- seq(0,600000,by=1000)
  R <- matrix(NA, ncol=length(unique(dat$year)),nrow=length(S))
  
  for(i in unique(dat$year)){

    alpha<- dat$alpha[dat$year==i]
    beta<- dat$beta[dat$year==i]
    R[,which(unique(dat$year)==i)]<-S*exp(alpha-beta*S)
  }
    
  actualSR_sena_hsmax[[a]]<-data.frame(year=rep(unique(dat$year),
      each=length(S)),
      spawners=S,
      recruits=c(R),
      scenario=simPar_sena_hsmax$scenario[a])

}

SRdf_sena_hsmax<-do.call(rbind,actualSR_sena_hsmax)
datdf_sena_hsmax<-do.call(rbind,alldat_sena_hsmax)
unique(SRdf_sena_hsmax$scenario)
SRdf_sena_hsmax$scenario_f <-factor(SRdf_sena_hsmax$scenario, levels=c("trendLinearProd1_halfbeta",
                                                                       "trendLinearProd2_halfbeta", 
                                                                       "trendLinearProd5_halfbeta",
                                                                    "trendLinearProd7_halfbeta",
                                                                    "regimeProd1_halfbeta",
                                                                    "regimeProd2_halfbeta",     
                                                                    "regimeProd5_halfbeta",
                                                                    "regimeProd7_halfbeta"))   



datdf_sena_hsmax$scenario_f <-factor(datdf_sena_hsmax$scenario, levels=c("trendLinearProd1_halfbeta",
                                                            "trendLinearProd2_halfbeta", 
                                                            "trendLinearProd5_halfbeta",
                                                            "trendLinearProd7_halfbeta",
                                                            "regimeProd1_halfbeta",
                                                            "regimeProd2_halfbeta",     
                                                            "regimeProd5_halfbeta",
                                                            "regimeProd7_halfbeta"))




SRexample_sena_hsmax<-  ggplot(SRdf_sena_hsmax) +
    geom_line(aes(x=spawners,y=recruits, col=as.factor(year)),linewidth=2) +
    mytheme + 
    theme(legend.position="right") +
    scale_colour_viridis_d(end=.85) +
    labs(col = "year") +
    geom_point(data=datdf_sena_hsmax,aes(x=spawners,y=recruits,col=as.factor(year)),alpha=.5) +
    facet_wrap(~scenario_f)
SRexample_sena_hsmax    
 
ggsave(
      filename = "outs/SamSimOutputs/plotcheck/srcurve_sensitivitya_halfSmax.png", 
      plot = SRexample_sena, 
      width = 12, height = 6
    )



#----------------------------------------
#sensitivity Smax case                              |
#----------------------------------------


simPar_senSmax <- read.csv("data/Smax_sensitivity/SimPars.csv")


simData_senSmax<-list()
actualSR_senSmax<-list()
alldat_senSmax<-list()

for(a in seq_len(nrow(simPar_senSmax))){

  simData_senSmax[[a]] <- readRDS(paste0("outs/SamSimOutputs/simData/", 
                          simPar_senSmax$nameOM[a],"/",
                          simPar_senSmax$scenario[a],"/",
                          paste(simPar_senSmax$nameOM[a],"_", 
                          simPar_senSmax$nameMP[a], "_", 
                          "CUsrDat.RData",sep="")))$srDatout

  dat<-simData_senSmax[[a]] 
  dat<-dat[dat$year>(max(dat$year)-46),]
  dat <- dat[!is.na(dat$obsRecruits),]
  
  dat <- dat[dat$iteration==sample(unique(dat$iteration),1),]
  dat$scenario <- simPar_senSmax$scenario[a]
  alldat_senSmax[[a]]<-dat

  S <- seq(0,750000,by=1000)
  R <- matrix(NA, ncol=length(unique(dat$year)),nrow=length(S))
  
  for(i in unique(dat$year)){

    alpha<- dat$alpha[dat$year==i]
    beta<- dat$beta[dat$year==i]
    R[,which(unique(dat$year)==i)]<-S*exp(alpha-beta*S)
  }
    
  actualSR_senSmax[[a]]<-data.frame(year=rep(unique(dat$year),
      each=length(S)),
      spawners=S,
      recruits=c(R),
      scenario=simPar_senSmax$scenario[a])

}

SRdf_senSmax<-do.call(rbind,actualSR_senSmax)
datdf_senSmax<-do.call(rbind,alldat_senSmax)
unique(SRdf_senSmax$scenario)
SRdf_senSmax$scenario_f <-factor(SRdf_senSmax$scenario, levels=c("trendLinearSmax025",
                                                                "trendLinearSmax050",
                                                                 "trendLinearSmax150",
                                                                "trendLinearSmax200",
                                                                "trendLinearSmax300",
                                                                 "regimeSmax025",
                                                                 "regimeSmax050",
                                                                "regimeSmax150",
                                                                "regimeSmax200",      
                                                                "regimeSmax300"))   



datdf_senSmax$scenario_f <-factor(datdf_senSmax$scenario, levels=c("trendLinearSmax025",
                                                                "trendLinearSmax050",
                                                                 "trendLinearSmax150",
                                                                "trendLinearSmax200",
                                                                "trendLinearSmax300",
                                                                 "regimeSmax025",
                                                                 "regimeSmax050",
                                                                "regimeSmax150",
                                                                "regimeSmax200",      
                                                                "regimeSmax300"))




SRexample_senSmax<-  ggplot(SRdf_senSmax) +
    geom_line(aes(x=spawners,y=recruits, col=as.factor(year)),linewidth=2) +
    mytheme + 
    theme(legend.position="right") +
    scale_colour_viridis_d(end=.85) +
    labs(col = "year") +
    geom_point(data=datdf_senSmax,aes(x=spawners,y=recruits,col=as.factor(year)),alpha=.5) +
    facet_wrap(~scenario_f)
SRexample_senSmax    
 
ggsave(
      filename = "outs/SamSimOutputs/plotcheck/srcurve_sensitivity_Smax.png", 
      plot = SRexample_sena, 
      width = 12, height = 6
    )



#----------------------------------------
#sensitivity Smax case                              |
#----------------------------------------


simPar_senSmax_da <- read.csv("data/Smax_sensitivity_doublealpha/SimPars.csv")


simData_senSmax_da<-list()
actualSR_senSmax_da<-list()
alldat_senSmax_da<-list()

for(a in seq_len(nrow(simPar_senSmax_da))){

  simData_senSmax_da[[a]] <- readRDS(paste0("outs/SamSimOutputs/simData/", 
                          simPar_senSmax_da$nameOM[a],"/",
                          simPar_senSmax_da$scenario[a],"/",
                          paste(simPar_senSmax_da$nameOM[a],"_", 
                          simPar_senSmax_da$nameMP[a], "_", 
                          "CUsrDat.RData",sep="")))$srDatout

  dat<-simData_senSmax_da[[a]] 
  dat<-dat[dat$year>(max(dat$year)-46),]
  dat <- dat[!is.na(dat$obsRecruits),]
  
  dat <- dat[dat$iteration==sample(unique(dat$iteration),1),]
  dat$scenario <- simPar_senSmax_da$scenario[a]
  alldat_senSmax_da[[a]]<-dat

  S <- seq(0,750000,by=1000)
  R <- matrix(NA, ncol=length(unique(dat$year)),nrow=length(S))
  
  for(i in unique(dat$year)){

    alpha<- dat$alpha[dat$year==i]
    beta<- dat$beta[dat$year==i]
    R[,which(unique(dat$year)==i)]<-S*exp(alpha-beta*S)
  }
    
  actualSR_senSmax_da[[a]]<-data.frame(year=rep(unique(dat$year),
      each=length(S)),
      spawners=S,
      recruits=c(R),
      scenario=simPar_senSmax_da$scenario[a])

}

SRdf_senSmax_da<-do.call(rbind,actualSR_senSmax_da)
datdf_senSmax_da<-do.call(rbind,alldat_senSmax_da)
unique(SRdf_senSmax_da$scenario)
SRdf_senSmax_da$scenario_f <-factor(SRdf_senSmax_da$scenario, levels=c("trendLinearSmax025",
                                                                "trendLinearSmax050",
                                                                 "trendLinearSmax150",
                                                                "trendLinearSmax200",
                                                                "trendLinearSmax300",
                                                                 "regimeSmax025",
                                                                 "regimeSmax050",
                                                                "regimeSmax150",
                                                                "regimeSmax200",      
                                                                "regimeSmax300"))   



datdf_senSmax_da$scenario_f <-factor(datdf_senSmax_da$scenario, levels=c("trendLinearSmax025",
                                                                "trendLinearSmax050",
                                                                 "trendLinearSmax150",
                                                                "trendLinearSmax200",
                                                                "trendLinearSmax300",
                                                                 "regimeSmax025",
                                                                 "regimeSmax050",
                                                                "regimeSmax150",
                                                                "regimeSmax200",      
                                                                "regimeSmax300"))




SRexample_senSmax<-  ggplot(SRdf_senSmax_da) +
    geom_line(aes(x=spawners,y=recruits, col=as.factor(year)),linewidth=2) +
    mytheme + 
    theme(legend.position="right") +
    scale_colour_viridis_d(end=.85) +
    labs(col = "year") +
    geom_point(data=datdf_senSmax_da,aes(x=spawners,y=recruits,col=as.factor(year)),alpha=.5) +
    facet_wrap(~scenario_f)
SRexample_senSmax    
 
ggsave(
      filename = "outs/SamSimOutputs/plotcheck/srcurve_sensitivity_Smax.png", 
      plot = SRexample_sena, 
      width = 12, height = 6
    )




#----------------------------------------
#sigmalow sensitivity                   |
#----------------------------------------



#----------------------------------------
#sigmamed sensitivity                   |
#----------------------------------------



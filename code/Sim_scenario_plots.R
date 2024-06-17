#===============================================
#Ricker curve scenario plots
# april 2023
#===============================================
library(ggplot2)
library(dplyr)
library(samEst)
library(cowplot)



mytheme = list(
    theme_classic(14)+
        theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
              legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="bold", size=13),
              axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=13))
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





SRdf$scenario_f <-factor(SRdf$scenario, levels=c("stationary",  
                                                  "autocorr",
                                                  "sigmaShift",
                                                  "decLinearProd", 
                                                  "regimeProd", 
                                                  "sineProd",  
                                                  "shiftProd",
                                                  "regimeCap", 
                                                  "decLinearCap", 
                                                  "shiftCap", 
                                                  "regimeProdCap", 
                                                  "decLinearProdshiftCap"))   
                                                                  

SRdf$scencode<-dplyr::case_match(SRdf$scenario, 
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



SRdf$scencode <-factor(SRdf$scencode, levels=c("Base1","Base2","Base3",
             "Base4","Base5","Base6",
              "Base7","Base8","Base9",
               "Base10","Base11","Base12"))



datdf$scenario_f <-factor(datdf$scenario, levels=c("stationary","autocorr","sigmaShift",
              "decLinearProd", "regimeProd", "sineProd","shiftProd",
               "regimeCap", "decLinearCap", "shiftCap",
               "regimeProdCap",  "decLinearProdshiftCap"))


datdf$scencode<-dplyr::case_match(datdf$scenario, 
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


datdf$scencode <-factor(datdf$scencode, levels=c("Base1","Base2","Base3",
             "Base4","Base5","Base6",
              "Base7","Base8","Base9",
               "Base10","Base11","Base12"))




SRexample<-  ggplot(SRdf) +
    geom_line(aes(x=spawners,y=recruits, col=as.factor(year)),linewidth=2) +
    mytheme + 
    theme(legend.position="right") +
    scale_colour_viridis_d(end=.85) +
    labs(col = "year") +
    geom_point(data=datdf,aes(x=spawners,y=recruits,col=as.factor(year)),alpha=.5) +
    facet_wrap(~scencode)
SRexample    
 
ggsave(
      filename = "C:/Users/worc/Documents/timevar/Best-Practices-time-varying-salmon-SR-models/figures/scenarios/srcurve_basecase.png",
      plot = SRexample, 
      width = 12, height = 6
    )


#select scenarios for main paper
datdfselec<-datdf[datdf$scenario%in%c("autocorr","decLinearProd","sineProd","shiftProd",
    "decLinearCap","shiftCap","regimeProdCap","decLinearProdshiftCap"),]



datdfselec$genericscenario<-dplyr::case_match(datdfselec$scenario, 
     "autocorr"~"stationary",
      "decLinearProd"~"linear~decline",
      "sineProd"~"sine~fluctuation",
      "shiftProd"~"shift~decline",
    "decLinearCap"~"linear~decline",
    "shiftCap"~"shift~decline",
    "regimeProdCap"~ "shift~-~both",
    "decLinearProdshiftCap"~"mixed~trend"
      )   

datdfselec$paramvary<-dplyr::case_match(datdfselec$scenario, 
      "autocorr"~"stationary",
      "decLinearProd"~"log(alpha)",
      "sineProd"~"log(alpha)",
      "shiftProd"~"log(alpha)",
    "decLinearCap"~"S[max]",
    "shiftCap"~"S[max]",
    "regimeProdCap"~ "both",
    "decLinearProdshiftCap"~"both"
      ) 



datdfselec$paramvary<-factor(datdfselec$paramvary, 
    levels=c("stationary","log(alpha)","S[max]","both"))


SRdfselec<-SRdf[SRdf$scenario%in%c("autocorr","decLinearProd","sineProd","shiftProd",
    "decLinearCap","shiftCap","regimeProdCap","decLinearProdshiftCap"),]



SRdfselec$genericscenario<-dplyr::case_match(SRdfselec$scenario, 
      "autocorr"~"stationary",
      "decLinearProd"~"linear~decline",
      "sineProd"~"sine~fluctuation",
      "shiftProd"~"shift~decline",
    "decLinearCap"~"linear~decline",
    "shiftCap"~"shift~decline",
    "regimeProdCap"~ "shift~-~both",
    "decLinearProdshiftCap"~"mixed~trend"
      )   

SRdfselec$paramvary<-dplyr::case_match(SRdfselec$scenario, 
      "autocorr"~"stationary",
      "decLinearProd"~"log(alpha)",
      "sineProd"~"log(alpha)",
      "shiftProd"~"log(alpha)",
    "decLinearCap"~"S[max]",
    "shiftCap"~"S[max]",
    "regimeProdCap"~ "both",
    "decLinearProdshiftCap"~"both"
      ) 


SRdfselec$paramvary<-factor(SRdfselec$paramvary, 
    levels=c("stationary","log(alpha)","S[max]","both"))


SRexampleselec<-  ggplot(SRdfselec) +
    geom_line(aes(x=spawners,y=recruits, col=as.factor(year)),linewidth=2) +
    mytheme + 
    coord_cartesian(ylim = c(-0, 600000))+
    scale_colour_viridis_d(end=.85) +
    labs(col = "year") +
    geom_point(data=datdfselec,aes(x=spawners,y=recruits,col=as.factor(year)),alpha=.5) +
    facet_wrap(paramvary~genericscenario, ncol=8,  labeller=label_parsed)+ 
    theme(legend.position = "none",axis.text.x = element_text(angle = 45, vjust = 1.0, hjust=1))
SRexampleselec

unique(SRdfselec$genericscenario)

ggsave(
      filename = "C:/Users/worc/Documents/timevar/Best-Practices-time-varying-salmon-SR-models/figures/scenarios/srcurve_selecmain.png",
      plot = SRexampleselec, 
      width = 8, height = 4
    )


paramdf<-reshape2::melt(datdf,id.vars=c("iteration","year","CU","spawners","recruits","obsSpawners","obsRecruits",
                                "ER", "obsER", "targetER",   
                                "sMSY", "sGen", "uMSY", "scenario", "scenario_f","scencode"))

summary(paramdf)

paramdf$simulated<-dplyr::recode(paramdf$scenario, 
      "stationary"="simple",
      "autocorr"="simple",
      "sigmaShift"="simple", 
      "decLinearProd"="rw",
      "sineProd"="rw",
      "regimeProd"="hmm",
      "decLinearCap"="rw",
      "regimeCap"="hmm",
      "shiftCap"="hmm", 
      "shiftProd"="hmm",
      "regimeProdCap"="hmm",
      "decLinearProdshiftCap"="rw"
      )   



paramdf$paramch<-dplyr::recode(paramdf$scenario, 
      "stationary"="simple",
      "autocorr"="simple",
      "sigmaShift"="simple", 
      "decLinearProd"="tva",
      "sineProd"="tva",
      "regimeProd"="tva",
      "shiftProd"="tva",
      "decLinearCap"="tvb",
      "regimeCap"="tvb",
      "shiftCap"="tvb", 
      "regimeProdCap"="both",
      "decLinearProdshiftCap"="both"
      )   


paramdf<-paramdf[paramdf$variable!="beta",]

head(paramdf)
 
paramdf$scencode<-dplyr::case_match(paramdf$scenario, 
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



paramdf$scenario <-factor(paramdf$scenario, levels=c("stationary",
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
      "decLinearProdshiftCap"))




paramdf$scencode <-factor(paramdf$scencode, levels=c("Base1","Base2","Base3",
             "Base4","Base5","Base6",
              "Base7","Base8","Base9",
               "Base10","Base11","Base12"))

head(paramdf)


paramtraj<-ggplot(paramdf) +
    geom_line(aes(x=year,y=value,col=paramch), linewidth=2)+
    mytheme + 
    theme(legend.position="bottom") +
    scale_colour_viridis_d(end=.85) +
    labs(col = "year") +
    facet_grid(variable~scencode, scales="free")
paramtraj

ggsave(
      filename = "C:/Users/worc/Documents/timevar/Best-Practices-time-varying-salmon-SR-models/figures/scenarios/params_basecase.png",
      plot = paramtraj, 
      width = 12, height = 6
    )

paramtraj2<-ggplot(paramdf) +
    geom_line(aes(x=year,y=value,col=paramch), linewidth=2)+
    mytheme + 
    theme(legend.position="bottom") +
    scale_colour_viridis_d(end=.85) +
    labs(col = "year") +
    facet_grid(variable~scenario_f, scales="free")
paramtraj2

ggsave(
      filename = "C:/Users/worc/Documents/timevar/Best-Practices-time-varying-salmon-SR-models/figures/scenarios/params_basecase_scnname.png",
      plot = paramtraj2, 
      width = 19, height = 6
    )

paramdf[paramdf$scencode=="Base5"&paramdf$variable=="alpha",]


#selected scenarios


paramdfselec <- paramdf[paramdf$scenario%in%c("autocorr","decLinearProd","sineProd","shiftProd",
    "decLinearCap","shiftCap","regimeProdCap","decLinearProdshiftCap")&
     paramdf$variable!="sigma",]

paramdfselec$variable<-as.character(paramdfselec$variable)
paramdfselec$variable[paramdfselec$variable=="capacity"]<-paste("S[max]")
paramdfselec$variable[paramdfselec$variable=="alpha"]<-paste("log(alpha)")

unique(paramdfselec$variable)

paramdfselec$genericscenario<-dplyr::case_match(paramdfselec$scenario, 
      "autocorr"~"stationary",
      "decLinearProd"~"linear~decline",
      "sineProd"~"sine~fluctuation",
      "shiftProd"~"shift~decline",
    "decLinearCap"~"linear~decline",
    "shiftCap"~"shift~decline",
    "regimeProdCap"~ "shift~-~both",
    "decLinearProdshiftCap"~"mixed~trend"
      )   

paramdfselec$paramvary<-dplyr::case_match(paramdfselec$scenario, 
      "autocorr"~".",
      "decLinearProd"~"log(alpha)",
      "sineProd"~"log(alpha)",
      "shiftProd"~"log(alpha)",
    "decLinearCap"~"S[max]",
    "shiftCap"~"S[max]",
    "regimeProdCap"~ "both",
    "decLinearProdshiftCap"~"both"
      ) 



paramdfselec$paramvary<-factor(paramdfselec$paramvary, 
    levels=c(".","log(alpha)","S[max]","both"))

unique(paramdfselec$paramvary)


paramtrajselec<-ggplot(paramdfselec) +
    geom_line(aes(x=year-54,y=value,col=paramch), linewidth=2)+
    mytheme + 
    theme(legend.position="bottom") +
    scale_colour_viridis_d(end=.85,option="A") +
    labs(col = "year") +
    ylab("")+
    theme(legend.position="none")+
    facet_grid(variable~paramvary+genericscenario, scales="free", labeller=label_parsed)
paramtrajselec









ggsave(
      filename = "C:/Users/worc/Documents/timevar/Best-Practices-time-varying-salmon-SR-models/figures/scenarios/params_selecmain.png",
      plot = SRexampleselec, 
      width = 8, height = 4
    )



multi.page.scenario <- ggarrange( paramtrajselec,SRexampleselec,
                        nrow = 2, ncol = 1,
                        legend="none",
                        heights=c(1,1),
                        align="v")
multi.page.scenario




#use cowplot



logatrajselec<-ggplot(paramdfselec[paramdfselec$variable=="log(alpha)",]) +
    geom_line(aes(x=year-54,y=value,col=paramch), linewidth=2)+
    mytheme + 
    theme(legend.position="bottom") +
    scale_colour_viridis_d(end=.85,option="A") +
    xlab( "year") +
    ylab(expression(log(alpha)))+
    theme(legend.position="none",axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
    facet_grid(~paramvary+genericscenario, scales="free", labeller=label_parsed)
logatrajselec



yrcol<-viridis::viridis(n = 5)

length(unique(paramdfselec$year)))

smaxtrajselec<-ggplot(paramdfselec[paramdfselec$variable=="S[max]",]) +
    geom_line(aes(x=year-54,y=value/1000,col=paramch), linewidth=2)+
    mytheme + 
    theme(legend.position="bottom") +
    scale_colour_viridis_d(end=.85,option="A") +
    xlab( "year") +
    ylab(expression(paste(S[max]~"(1000s)")))+
    theme(legend.position="none",  strip.text.x = element_blank())+
    facet_grid(~paramvary+genericscenario, scales="free")
smaxtrajselec



#xLabVals <- as.numeric(ggplot_build(smaxtrajselec)$layout$panel_params[[1]]$x$get_labels()) 
#yrcol<-viridis::viridis(length(unique(xLabVals )),end=.85)


#smaxtrajselec<-smaxtrajselec+theme(axis.text.x = element_text(angle = 0, hjust = 1, colour = yrcol))




SRexampleselec<-  ggplot(SRdfselec) +
    geom_line(aes(x=spawners/1000,y=recruits/1000, col=as.factor(year)),linewidth=2) +
    mytheme + 
    coord_cartesian(ylim = c(0, 570))+
    scale_colour_viridis_d(end=.85) +
    xlab("spawners") +
    ylab("recruits") +
    geom_point(data=datdfselec,aes(x=spawners/1000,y=recruits/1000,col=as.factor(year)),alpha=.5) +
    facet_wrap(paramvary~genericscenario, ncol=8)+ 
    theme(legend.position = "none",axis.text.x = element_text(angle = 45, vjust = 1.0, hjust=1), strip.text.x = element_blank())
SRexampleselec


plot2<-logatrajselec/ smaxtrajselec/SRexampleselec
plot2

ggsave(
      filename = "C:/Users/worc/Documents/timevar/Best-Practices-time-varying-salmon-SR-models/figures/scenarios/multi.page.scenario_basecase.png",
      plot = plot2, 
      width = 12, height = 6
    )

#----------------------------------------
#equilibrium parameter plots                            |
#----------------------------------------


names(datdf)
eqdf<-datdf

eqdf$sMSY<-smsyCalc(eqdf$alpha,eqdf$beta)
eqdf$uMSY<-umsyCalc(eqdf$alpha)
eqdf$sGen<-unlist(mapply(sGenCalc,a=eqdf$alpha,
          Smsy=eqdf$sMSY, 
          b=eqdf$beta))


pareqdf<-reshape2::melt(eqdf,id.vars=c("iteration","year","CU","spawners","recruits","obsSpawners","obsRecruits",
                                "ER", "obsER", "targetER", "scenario", "scenario_f"))



pareqdf<-pareqdf[pareqdf$variable!="sigma"&,]


 ggplot(pareqdf) +
    geom_line(aes(x=year,y=value), linewidth=2)+
    mytheme + 
    theme(legend.position="right") +
    scale_colour_viridis_d(end=.85) +
    labs(col = "year") +
    facet_grid(variable~scenario_f, scales="free")

#----------------------------------------
#sensitivity alpha case                              |
#----------------------------------------

simPar_sena <- read.csv("data/sensitivity/SimPars.csv")

simPar_sena$scenario
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
                                                           "trendLinearProd10",
                                                            "regimeProd1",      "regimeProd2",
                                                            "regimeProd5",      "regimeProd7",
                                                            "regimeProd10"))   



datdf_sena$scenario_f <-factor(datdf_sena$scenario, levels=c("trendLinearProd1", "trendLinearProd2",
                                                           "trendLinearProd5", "trendLinearProd7",
                                                           "trendLinearProd10",
                                                            "regimeProd1",      "regimeProd2",
                                                            "regimeProd5",      "regimeProd7",
                                                            "regimeProd10"))




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
SRdf_sena_hsmax$scenario_f <-factor(SRdf_sena_hsmax$scenario, levels=c("trendLinearProd1_halfSmax",
                                                                       "trendLinearProd2_halfSmax", 
                                                                       "trendLinearProd5_halfSmax",
                                                                    "trendLinearProd7_halfSmax",
                                                                    "trendLinearProd10_halfSmax",
                                                                    "regimeProd1_halfSmax",
                                                                    "regimeProd2_halfSmax",     
                                                                    "regimeProd5_halfSmax",
                                                                    "regimeProd7_halfSmax",
                                                                    "regimeProd10_halfSmax"))   



datdf_sena_hsmax$scenario_f <-factor(datdf_sena_hsmax$scenario, levels=c("trendLinearProd1_halfSmax",
                                                                       "trendLinearProd2_halfSmax", 
                                                                       "trendLinearProd5_halfSmax",
                                                                    "trendLinearProd7_halfSmax",
                                                                    "trendLinearProd10_halfSmax",
                                                                    "regimeProd1_halfSmax",
                                                                    "regimeProd2_halfSmax",     
                                                                    "regimeProd5_halfSmax",
                                                                    "regimeProd7_halfSmax",
                                                                    "regimeProd10_halfSmax"))




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

  S <- seq(0,1300000,by=1000)
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


names(datdf)

datdf[datdf$scenario_f=="stationary",]


paramdf<-reshape2::melt(datdf,id.vars=c("iteration","year","CU","spawners","recruits","obsSpawners","obsRecruits",
                                "ER", "obsER", "targetER",   
                                "sMSY", "sGen", "uMSY", "scenario", "scenario_f"))

summary(paramdf)

paramdf$simulated<-dplyr::recode(paramdf$scenario, 
      "stationary"="simple",
      "autocorr"="simple",
      "sigmaShift"="simple", 
      "decLinearProd"="rw",
      "sineProd"="rw",
      "regimeProd"="hmm",
      "decLinearCap"="rw",
      "regimeCap"="hmm",
      "shiftCap"="hmm", 
      "shiftProd"="hmm",
      "regimeProdCap"="hmm",
      "decLinearProdshiftCap"="rw"
      )   

paramdf<-paramdf[paramdf$variable!="beta",]

 ggplot(paramdf) +
    geom_line(aes(x=year,y=value,col=simulated), linewidth=2)+
    mytheme + 
    theme(legend.position="right") +
    scale_colour_viridis_d(end=.85) +
    labs(col = "year") +
    facet_grid(variable~scenario_f, scales="free")


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

  S <- seq(0,1300000,by=1000)
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
SRdf_senSmax_da$scenario_f <-factor(SRdf_senSmax_da$scenario, levels=c("trendLinearSmax025_da",
                                                                "trendLinearSmax050_da",
                                                                 "trendLinearSmax150_da",
                                                                "trendLinearSmax200_da",
                                                                "trendLinearSmax300_da",
                                                                 "regimeSmax025_da",
                                                                 "regimeSmax050_da",
                                                                "regimeSmax150_da",
                                                                "regimeSmax200_da",      
                                                                "regimeSmax300_da"))   



datdf_senSmax_da$scenario_f <-factor(datdf_senSmax_da$scenario, levels=c("trendLinearSmax025_da",
                                                                "trendLinearSmax050_da",
                                                                 "trendLinearSmax150_da",
                                                                "trendLinearSmax200_da",
                                                                "trendLinearSmax300_da",
                                                                 "regimeSmax025_da",
                                                                 "regimeSmax050_da",
                                                                "regimeSmax150_da",
                                                                "regimeSmax200_da",      
                                                                "regimeSmax300_da"))




SRexample_senSmax_da<-  ggplot(SRdf_senSmax_da) +
    geom_line(aes(x=spawners,y=recruits, col=as.factor(year)),linewidth=2) +
    mytheme + 
    theme(legend.position="right") +
    scale_colour_viridis_d(end=.85) +
    labs(col = "year") +
    geom_point(data=datdf_senSmax_da,aes(x=spawners,y=recruits,col=as.factor(year)),alpha=.5) +
    facet_wrap(~scenario_f,scales = "free")
SRexample_senSmax_da    
 
ggsave(
      filename = "outs/SamSimOutputs/plotcheck/srcurve_sensitivity_Smax.png", 
      plot = SRexample_sena, 
      width = 12, height = 6
    )




#----------------------------------------
#sigmalow sensitivity                   |
#----------------------------------------

simPar_siglow <- read.csv("data/sigmalow_sensitivity/SimPars.csv")



simData_siglow<-list()
actualSR_siglow<-list()
alldat_siglow<-list()

for(a in seq_len(nrow(simPar_siglow))){

  simData_siglow[[a]] <- readRDS(paste0("outs/SamSimOutputs/simData/", 
                          simPar_siglow$nameOM[a],"/",
                          simPar_siglow$scenario[a],"/",
                          paste(simPar_siglow$nameOM[a],"_", 
                          simPar_siglow$nameMP[a], "_", 
                          "CUsrDat.RData",sep="")))$srDatout

  dat<-simData_siglow[[a]] 
  dat<-dat[dat$year>(max(dat$year)-46),]
  dat <- dat[!is.na(dat$obsRecruits),]
  
  dat <- dat[dat$iteration==sample(unique(dat$iteration),1),]
  dat$scenario <- simPar_siglow$scenario[a]
  alldat_siglow[[a]]<-dat

  S <- seq(0,750000,by=1000)
  R <- matrix(NA, ncol=length(unique(dat$year)),nrow=length(S))
  
  for(i in unique(dat$year)){

    alpha<- dat$alpha[dat$year==i]
    beta<- dat$beta[dat$year==i]
    R[,which(unique(dat$year)==i)]<-S*exp(alpha-beta*S)
  }
    
  actualSR_siglow[[a]]<-data.frame(year=rep(unique(dat$year),
      each=length(S)),
      spawners=S,
      recruits=c(R),
      scenario=simPar_siglow$scenario[a])

}

SRdf_siglow<-do.call(rbind,actualSR_siglow)
datdf_siglow<-do.call(rbind,alldat_siglow)
unique(SRdf_siglow$scenario)
SRdf_siglow$scenario_f <-factor(SRdf_siglow$scenario, levels=c("sigmalow_stationary",
                                                               "sigmalow_decLinearProd",
                                                               "sigmalow_regimeProd",           
                                                               "sigmalow_sineProd",
                                                               "sigmalow_regimeCap",
                                                               "sigmalow_decLinearCap",         
                                                               "sigmalow_regimeProdCap",
                                                               "sigmalow_shiftCap",
                                                               "sigmalow_decLinearProdshiftCap"))   



datdf_siglow$scenario_f <-factor(datdf_siglow$scenario, levels=c("sigmalow_stationary",
                                                               "sigmalow_decLinearProd",
                                                               "sigmalow_regimeProd",           
                                                               "sigmalow_sineProd",
                                                               "sigmalow_regimeCap",
                                                               "sigmalow_decLinearCap",         
                                                               "sigmalow_regimeProdCap",
                                                               "sigmalow_shiftCap",
                                                               "sigmalow_decLinearProdshiftCap"))




SRexample_siglow <- ggplot(SRdf_siglow) +
    geom_line(aes(x=spawners,y=recruits, col=as.factor(year)),linewidth=2) +
    mytheme + 
    theme(legend.position="right") +
    scale_colour_viridis_d(end=.85) +
    labs(col = "year") +
    geom_point(data=datdf_siglow,aes(x=spawners,y=recruits,col=as.factor(year)),alpha=.5) +
    facet_wrap(~scenario_f)
SRexample_siglow    
 
ggsave(
      filename = "outs/SamSimOutputs/plotcheck/srcurve_sensitivity_siglow.png", 
      plot = SRexample_siglow, 
      width = 12, height = 6
    )




#----------------------------------------
#sigmamed sensitivity                   |
#----------------------------------------


simPar_sigmed <- read.csv("data/sigmamed_sensitivity/SimPars.csv")



simData_sigmed <- list()
actualSR_sigmed <- list()
alldat_sigmed <- list()

for(a in seq_len(nrow(simPar_sigmed))){
  
  simData_sigmed[[a]] <- readRDS(paste0("outs/SamSimOutputs/simData/", 
                                        simPar_sigmed$nameOM[a],"/",
                                        simPar_sigmed$scenario[a],"/",
                                        paste(simPar_sigmed$nameOM[a],"_", 
                                              simPar_sigmed$nameMP[a], "_", 
                                              "CUsrDat.RData",sep="")))$srDatout
  
  dat <- simData_sigmed[[a]] 
  dat <- dat[dat$year>(max(dat$year)-46),]
  dat <- dat[!is.na(dat$obsRecruits),]
  
  dat <- dat[dat$iteration==sample(unique(dat$iteration),1),]
  dat$scenario <- simPar_sigmed$scenario[a]
  alldat_sigmed[[a]]<-dat
  
  S <- seq(0,750000,by=1000)
  R <- matrix(NA, ncol=length(unique(dat$year)),nrow=length(S))
  
  for(i in unique(dat$year)){
    
    alpha<- dat$alpha[dat$year==i]
    beta<- dat$beta[dat$year==i]
    R[,which(unique(dat$year)==i)]<-S*exp(alpha-beta*S)
  }
  
  actualSR_sigmed[[a]]<-data.frame(year=rep(unique(dat$year),
                                            each=length(S)),
                                   spawners=S,
                                   recruits=c(R),
                                   scenario=simPar_sigmed$scenario[a])
  
}

SRdf_sigmed <- do.call(rbind,actualSR_sigmed)
datdf_sigmed <- do.call(rbind,alldat_sigmed)
unique(SRdf_sigmed$scenario)
SRdf_sigmed$scenario_f <-factor(SRdf_sigmed$scenario, levels=c("sigmamed_stationary",
                                                               "sigmamed_decLinearProd",
                                                               "sigmamed_regimeProd",           
                                                               "sigmamed_sineProd",
                                                               "sigmamed_regimeCap",
                                                               "sigmamed_decLinearCap",         
                                                               "sigmamed_regimeProdCap",
                                                               "sigmamed_shiftCap",
                                                               "sigmamed_decLinearProdshiftCap"))   



datdf_sigmed$scenario_f <-factor(datdf_sigmed$scenario, levels=c("sigmamed_stationary",
                                                                 "sigmamed_decLinearProd",
                                                                 "sigmamed_regimeProd",           
                                                                 "sigmamed_sineProd",
                                                                 "sigmamed_regimeCap",
                                                                 "sigmamed_decLinearCap",         
                                                                 "sigmamed_regimeProdCap",
                                                                 "sigmamed_shiftCap",
                                                                 "sigmamed_decLinearProdshiftCap"))




SRexample_sigmed <- ggplot(SRdf_sigmed) +
  geom_line(aes(x=spawners,y=recruits, col=as.factor(year)),linewidth=2) +
  mytheme + 
  theme(legend.position="right") +
  scale_colour_viridis_d(end=.85) +
  labs(col = "year") +
  geom_point(data=datdf_sigmed,aes(x=spawners,y=recruits,col=as.factor(year)),alpha=.5) +
  facet_wrap(~scenario_f)
SRexample_sigmed    

ggsave(
  filename = "outs/SamSimOutputs/plotcheck/srcurve_sensitivity_siglow.png", 
  plot = SRexample_siglow, 
  width = 12, height = 6
)





#-----------------------------------------------------------------
#ER trend plots



simPar_ER <- read.csv("data/genericER/SimPars_ER.csv")



simData_ER <- list()
actualSR_ER <- list()
alldat_ER <- list()

for(a in seq_len(nrow(simPar_ER))){


  simData_ER[[a]] <- readRDS(paste0("outs/SamSimOutputs/", 
                          simPar_ER$nameOM[a],"/",
                          simPar_ER$scenario[a],"/",
                          paste(simPar_ER$nameOM[a],"_", 
                          simPar_ER$nameMP[a], "_", 
                          "CUsrDat.RData",sep="")))$srDatout

  dat <- simData_ER[[a]] 
  dat <- dat[dat$year>(max(dat$year)-46),]
  dat <- dat[!is.na(dat$obsRecruits),]
  
  dat <- dat[dat$iteration==sample(unique(dat$iteration),1),]
  dat$scenario <- simPar_ER$scenario[a]
  alldat_ER[[a]]<-dat
  
  S <- seq(0,750000,by=1000)
  R <- matrix(NA, ncol=length(unique(dat$year)),nrow=length(S))
  
  for(i in unique(dat$year)){
    
    alpha<- dat$alpha[dat$year==i]
    beta<- dat$beta[dat$year==i]
    R[,which(unique(dat$year)==i)]<-S*exp(alpha-beta*S)
  }
  
  actualSR_ER[[a]]<-data.frame(year=rep(unique(dat$year),
                                            each=length(S)),
                                   spawners=S,
                                   recruits=c(R),
                                   scenario=simPar_ER$scenario[a])
  


  

  
  
}


ggplot(dat) +
      geom_boxplot(aes(x=as.factor(year),y=ER),alpha=.5) +
      geom_hline(data=dat,aes(yintercept=uMSY, col=as.factor(year))) +
      theme_bw(14) + 
      xlab("year") + 
      scale_color_discrete(name = "Umsy") +
      labs(title = simPar$nameOM[a])
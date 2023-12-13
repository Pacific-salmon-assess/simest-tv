#routines to visualize smulation estimation Results
#Catarina Wor
#September 2023
#============================================


#remotes::install_git("https://gitlab.com/michaelfolkes/performancer") 

library(ggplot2)
library(gridExtra)
library(dplyr)
library(cowplot)
library(performancer)
library(ggpattern)

source("code/utils.R")
source("code/cluster_func_plots.R")

mytheme = list(
    theme_classic(16)+
        theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
              legend.position="bottom", strip.text = element_text(face="bold", size=12),
              axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
)
#legend.title = element_blank(),


#========================================================================================================
#sensitivity Alpha cases
#read in data
simPar <- read.csv("data/Smax_sensitivity/SimPars.csv")

## Store relevant object names to help run simulation 
scenNames <- unique(simPar$scenario)
simPar <- read.csv("data/Smax_sensitivity/SimPars.csv")

## Store relevant object names to help run simulation 
scenNames <- unique(simPar$scenario)

ressmax1<-readRDS(file = "outs/simest/Smax_sensitivity/res_smax1.rds")
ressmax2<-readRDS(file = "outs/simest/Smax_sensitivity/res_smax2.rds")
ressmax240<-readRDS(file = "outs/simest/Smax_sensitivity/res_smax_240.rds")
restmb<-rbind(ressmax1,ressmax2,ressmax240)

head(restmb)

resstansmax1<-readRDS(file = "outs/simest/Smax_sensitivity/resstan_smax1.rds")
resstansmax2<-readRDS(file = "outs/simest/Smax_sensitivity/resstan_smax2.rds")

resstan<-rbind(resstansmax1,resstansmax2)

#resstan<-readRDS(file = "outs/simest/generic/resstan.rds")



res<-rbind(restmb,resstan)


res$parameter[res$parameter=="Smax"]<-"smax"

#===========================================================
#MCMC results  

resparam<-res[res$parameter%in%c("alpha","smax","smsy"),]
resparam<-resparam[resparam$convergence==0&resparam$method=="MCMC",]


scns<-unique(resparam$scenario)
mods<-unique(resparam$model)
paras<-unique(resparam$parameter)

pml<-list()
a=1

for(i in seq_along(scns)){
    for(j in seq_along(mods)){
        for(u in seq_along(paras)){
          dd <- filter(resparam,  scenario == scns[i] & model == mods[j],parameter==paras[u])
          aa<-calcMetrics(
            expect=dd$mode,
            obs=dd$sim,
            metrics=c("mae","rmse","mape","sd2"))
            #"mre",excluded bias only
            #"mae", mean bsolute error - bias and precision, no direction
            #"mpe",excluded, bias
            #"mape",mean absolute percent error. bias and precision, no direction - will change as variable value changes maybe not ideal
            #"mse" - over values outliers -excluded
            #?rmse -  precision
            #"sd2"
            #"all", "mre", "mae", "mpe", "mape", "mrpd", "mse", "rmse", "sd2", "ustat2",
            #"ustat2mean", "mase", "maz"))
          pml[[a]]<-cbind(scns[i],mods[j],paras[u],aa$m)
          a<-a+1
        }
    }
}

performance.results.df <- do.call(rbind,pml)
head(performance.results.df )
ranksbase<-list()
n=1
for(i in seq_along(scns)){
    for(u in seq_along(paras)){
    tmp.results.df<-performance.results.df[performance.results.df$"scns[i]"==scns[i]&performance.results.df$"paras[u]"==paras[u],]
    ranksbase[[n]]<-calcRanks(tmp.results.df, columnToRank = 4:7, ranking ="relative1")$ranks
    n<-n+1
}}
unique(ranksbaseall$"scns[i]")

ranksbaseall<-do.call(rbind,ranksbase)
ranksbaseall$model=factor(ranksbaseall$"mods[j]", levels=c("simple",
   "autocorr", "rwa", "hmma","rwb",  "hmmb","rwab", "hmmab"))
ranksbaseall$parameter=factor(ranksbaseall$"paras[u]", levels=c("alpha","smax","smsy"))
ranksbaseall$scenario=factor(ranksbaseall$"scns[i]", levels=c("trendLinearSmax025",
                            "trendLinearSmax050",
                            "trendLinearSmax150",
                            "trendLinearSmax200",
                            "trendLinearSmax300",
                            "regimeSmax025",
                            "regimeSmax050",
                            "regimeSmax150",      
                            "regimeSmax200",      
                            "regimeSmax300" ))

ranksbaseall$parameter=factor(ranksbaseall$"paras[u]", levels=c("alpha","smax","smsy"))

ranksbaseall$paramch<-"a"


ranksbaseall$modtype<-recode(ranksbaseall$model, 
    "simple"="stationary",
   "autocorr"="stationary", 
   "rwa"="tva",
    "rwb"="tvb",
   "rwab"="tvab", 
   "hmma"="tva",
     "hmmb"="tvb",
     "hmmab"="tvab")

ranksbaseall$scnmodmatch<-0
unique(ranksbaseall$scenario)


ranksbaseall$scnmodmatch[ranksbaseall$scenario %in% c("trendLinearSmax025",
                            "trendLinearSmax050",
                            "trendLinearSmax150",
                            "trendLinearSmax200",
                            "trendLinearSmax300" )&
                        ranksbaseall$model=="rwb"]<-1

ranksbaseall$scnmodmatch[ranksbaseall$scenario %in% c( "regimeSmax025",
                            "regimeSmax050",
                            "regimeSmax150",      
                            "regimeSmax200",      
                            "regimeSmax300" )&
                        ranksbaseall$model=="hmmb"]<-1

ranksbaseall$scnmodmatch<-as.factor(ranksbaseall$scnmodmatch)


ranksbaseall$scnmodmatch2<-"0"

for(i in seq_along(scns)){
    for(u in seq_along(paras)){
    
    tmp<-min(ranksbaseall[ranksbaseall$scenario==scns[i]&
                        ranksbaseall$parameter==paras[u],"average.rank"])

    ranksbaseall$scnmodmatch2[ranksbaseall$scenario==scns[i]&
                        ranksbaseall$parameter==paras[u]&
                        ranksbaseall$average.rank==tmp]<-
                        "1"
  }
}



ranksbaseall$scnmodmatch2<-as.factor(ranksbaseall$scnmodmatch2)

ranksbaseall$scenario2 <- 
case_match(
  ranksbaseall$scenario,
   "trendLinearSmax025"~"trend smax*0.25"  ,
    "trendLinearSmax050"~"trend smax*0.50",
    "trendLinearSmax150"~"trend smax*1.5",
    "trendLinearSmax200"~"trend smax*2.0",
    "trendLinearSmax300"~"trend smax*3.0",
    "regimeSmax025"~"regime smax*0.25",
    "regimeSmax050"~"regime smax*0.50",
    "regimeSmax150"~"regime smax*1.5",      
    "regimeSmax200"~"regime smax*2.0",      
    "regimeSmax300"~"regime smax*3.0"
)

ranksbaseall$scenario2<-factor(ranksbaseall$scenario2, levels=c( "trend smax*0.25",
"trend smax*0.50",
  "trend smax*1.5", 
 "trend smax*2.0", 
  "trend smax*3.0",
  "regime smax*0.25",      
  "regime smax*0.50",      
  "regime smax*1.5",      
  "regime smax*2.0",       
  "regime smax*3.0" ))

unique(ranksbaseall$scenario2)


pranks_senssmax<-ggplot(data=ranksbaseall, aes(x=model, y=average.rank,fill=modtype,color=scnmodmatch ))+
     geom_bar_pattern(stat="identity", linewidth=1.3, aes(pattern_shape=scnmodmatch2 ),pattern = 'pch')+
    scale_pattern_shape_manual(values = c(NA, 16))+
   scale_fill_grey(start=0.4)+
   scale_color_manual(values=c("transparent",                             
                             "firebrick2"))+
    mytheme+
    theme(axis.text.x = element_text(angle = 90,vjust=-0), legend.position="none") +
   facet_grid(parameter~scenario2)
pranks_senssmax


ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/ranks/sens_smax/senssmax_average_rank.png", plot=pranks_sensa)




#=========================================================================
#Rmse
#c("rmse","mae","mape","mpe","mre","sd2")

ranksbaseall$scnmodmatch_rmse<-"0"

for(i in seq_along(scns)){
    for(u in seq_along(paras)){
    
    tmp<-min(ranksbaseall[ranksbaseall$scenario==scns[i]&
                        ranksbaseall$parameter==paras[u],"rmse.rank"])

    ranksbaseall$scnmodmatch_rmse[ranksbaseall$scenario==scns[i]&
                        ranksbaseall$parameter==paras[u]&
                        ranksbaseall$rmse.rank==tmp]<-
                        "1"
  }
}


pranks_rmse<-ggplot(data=ranksbaseall, aes(x=model, y=rmse.rank,fill=modtype,color=scnmodmatch ))+
     geom_bar_pattern(stat="identity", linewidth=1.3, aes(pattern_shape=scnmodmatch_rmse),pattern = 'pch', 
        pattern_aspect_ratio = 1, 
    pattern_density      = 0.5)+
    scale_pattern_shape_manual(values = c(NA, 16
                                 ))+
    scale_pattern_size_manual(values = c(NA, 2))+
   scale_fill_grey(start=0.4)+
   scale_color_manual(values=c("transparent",                           
                             "firebrick2"#,
                             ))+
    mytheme+
    theme(axis.text.x = element_text(angle = 90,vjust=-0), legend.position="none") +
   facet_grid(parameter~scenario2)
pranks_rmse
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/ranks/sens_smax/senssmax_rmse_rank.png", plot=pranks_rmse)

#=================================================================================
#"mae"
#c("rmse","mae","mape","mpe","mre","sd2")

ranksbaseall$scnmodmatch_mae<-"0"

for(i in seq_along(scns)){
    for(u in seq_along(paras)){
    
    tmp<-min(ranksbaseall[ranksbaseall$scenario==scns[i]&
                        ranksbaseall$parameter==paras[u],"mae.rank"])

    ranksbaseall$scnmodmatch_mae[ranksbaseall$scenario==scns[i]&
                        ranksbaseall$parameter==paras[u]&
                        ranksbaseall$mae.rank==tmp]<-
                        "1"
  }
}


pranks_mae<-ggplot(data=ranksbaseall, aes(x=model, y=mae.rank,fill=modtype,color=scnmodmatch ))+
     geom_bar_pattern(stat="identity", linewidth=1.3, aes(pattern_shape=scnmodmatch_mae),pattern = 'pch',
         pattern_aspect_ratio = 1, 
    pattern_density      = 0.5)+
    scale_pattern_shape_manual(values = c(NA, 16
                                 ))+
   scale_fill_grey(start=0.4)+
   scale_color_manual(values=c("transparent",                           
                             "firebrick2"#,
                             ))+
    mytheme+
    theme(axis.text.x = element_text(angle = 90,vjust=-0)) +
   facet_grid(parameter~scenario2)
pranks_mae

ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/ranks/sens_smax/senssmax_mae_rank.png", plot=pranks_mae)



#=================================================================================
#"mape
#c("rmse","mae","mape","mpe","mre","sd2")

ranksbaseall$scnmodmatch_mape<-"0"

for(i in seq_along(scns)){
    for(u in seq_along(paras)){
    
    tmp<-min(ranksbaseall[ranksbaseall$scenario==scns[i]&
                        ranksbaseall$parameter==paras[u],"mape.rank"])

    ranksbaseall$scnmodmatch_mape[ranksbaseall$scenario==scns[i]&
                        ranksbaseall$parameter==paras[u]&
                        ranksbaseall$mape.rank==tmp]<-
                        "1"
  }
}


pranks_mape<-ggplot(data=ranksbaseall, aes(x=model, y=mape.rank,fill=modtype,color=scnmodmatch ))+
     geom_bar_pattern(stat="identity", linewidth=1.3, aes(pattern_shape=scnmodmatch_mape),pattern = 'pch',
        pattern_aspect_ratio = 1, 
    pattern_density      = 0.5)+
    scale_pattern_shape_manual(values = c(NA, 16 ))+
   scale_fill_grey(start=0.4)+
   scale_color_manual(values=c("transparent",                           
                             "firebrick2"#,
                             ))+
    mytheme+
    theme(axis.text.x = element_text(angle = 90,vjust=-0)) +
   facet_grid(parameter~scenario2)
pranks_mape
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/ranks/sens_smax/senssmax_mape_rank.png", plot=pranks_mape)


#=================================================================================
#sd2
#c("rmse","mae","mape","mpe","mre","sd2")

ranksbaseall$scnmodmatch_sd2<-"0"

for(i in seq_along(scns)){
    for(u in seq_along(paras)){
    
    tmp<-min(ranksbaseall[ranksbaseall$scenario==scns[i]&
                        ranksbaseall$parameter==paras[u],"sd2.rank"])

    ranksbaseall$scnmodmatch_sd2[ranksbaseall$scenario==scns[i]&
                        ranksbaseall$parameter==paras[u]&
                        ranksbaseall$sd2.rank==tmp]<-
                        "1"
  }
}


pranks_sd2<-ggplot(data=ranksbaseall, aes(x=model, y=sd2.rank,fill=modtype,color=scnmodmatch ))+
     geom_bar_pattern(stat="identity", linewidth=1.3, aes(pattern_shape=scnmodmatch_sd2),pattern = 'pch',
        pattern_aspect_ratio = 1, 
    pattern_density      = 0.5)+
    scale_pattern_shape_manual(values = c(NA, 16))+
   scale_fill_grey(start=0.4)+
   scale_color_manual(values=c("transparent",                           
                             "firebrick2"))+
    mytheme+
    theme(axis.text.x = element_text(angle = 90,vjust=-0)) +
   facet_grid(parameter~scenario2)
pranks_sd2
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/ranks/sens_smax/senssmax_sd2_rank.png", plot=pranks_sd2)



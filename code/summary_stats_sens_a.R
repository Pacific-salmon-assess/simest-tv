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
simPar <- read.csv("data/sensitivity/SimPars.csv")

## Store relevant object names to help run simulation 
scenNames <- unique(simPar$scenario)

resa1<-readRDS(file = "outs/simest/sensitivity/res_a1.rds")
resa2<-readRDS(file = "outs/simest/sensitivity/res_a2.rds")


restmb<-rbind(resa1,resa2)

resstana1<-readRDS(file = "outs/simest/sensitivity/resstan_a1.rds")
resstana2<-readRDS(file = "outs/simest/sensitivity/resstan_a2.rds")
resstan<-rbind(resstana1,resstana2)

 names(resstan)

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


ranksbaseall<-do.call(rbind,ranksbase)
ranksbaseall$model=factor(ranksbaseall$"mods[j]", levels=c("simple",
   "autocorr", "rwa", "hmma","rwb",  "hmmb","rwab", "hmmab"))
ranksbaseall$parameter=factor(ranksbaseall$"paras[u]", levels=c("alpha","smax","smsy"))
ranksbaseall$scenario=factor(ranksbaseall$"scns[i]", levels=c("trendLinearProd1" ,
                                                        "trendLinearProd2", 
                                                        "trendLinearProd5", 
                                                        "trendLinearProd7",
                                                        "trendLinearProd10" ,    
                                                        "regimeProd1",  
                                                        "regimeProd2" ,     
                                                        "regimeProd5",     
                                                        "regimeProd7",
                                                        "regimeProd10" ))

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

ranksbaseall$scnmodmatch[ranksbaseall$scenario %in% c("trendLinearProd1" ,
                                                        "trendLinearProd2", 
                                                        "trendLinearProd5", 
                                                        "trendLinearProd7",
                                                        "trendLinearProd10" )&
                        ranksbaseall$model=="rwa"]<-1

ranksbaseall$scnmodmatch[ranksbaseall$scenario %in% c( "regimeProd1",  
                                                        "regimeProd2" ,     
                                                        "regimeProd5",     
                                                        "regimeProd7",
                                                        "regimeProd10" )&
                        ranksbaseall$model=="hmma"]<-1

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




pranks_sensa<-ggplot(data=ranksbaseall, aes(x=model, y=average.rank,fill=modtype,color=scnmodmatch ))+
     geom_bar_pattern(stat="identity", size=1.3, aes(pattern_shape=scnmodmatch2 ),pattern = 'pch')+
    scale_pattern_shape_manual(values = c(NA, 16))+
   scale_fill_grey(start=0.3)+
   scale_color_manual(values=c("transparent",                             
                             "firebrick2"))+
    mytheme+
    theme(axis.text.x = element_text(angle = 90,vjust=-0), legend.position="none") +
   facet_grid(parameter~scenario)
pranks_sensa
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/ranks/sens_a/sensa_average_rank.png", plot=pranks_sensa)




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
     geom_bar_pattern(stat="identity", size=1.3, aes(pattern_shape=scnmodmatch_rmse),pattern = 'pch')+
    scale_pattern_shape_manual(values = c(NA, 16
                                 ))+
   scale_fill_grey(start=0.3)+
   scale_color_manual(values=c("transparent",                           
                             "firebrick2"#,
                             ))+
    mytheme+
    theme(axis.text.x = element_text(angle = 90,vjust=-0), legend.position="none") +
   facet_grid(parameter~scenario)
pranks_rmse
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/ranks/base_rmse_rank.png", plot=pranks_rmse)

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
     geom_bar_pattern(stat="identity", size=1.3, aes(pattern_shape=scnmodmatch_mae),pattern = 'pch')+
    scale_pattern_shape_manual(values = c(NA, 16
                                 ))+
   scale_fill_grey(start=0.4)+
   scale_color_manual(values=c("transparent",                           
                             "firebrick2"#,
                             ))+
    mytheme+
    theme(axis.text.x = element_text(angle = 90,vjust=-0)) +
   facet_grid(parameter~scenario)
pranks_mae

ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/ranks/base_mae_rank.png", plot=pranks_mae)



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
     geom_bar_pattern(stat="identity", size=1.3, aes(pattern_shape=scnmodmatch_mape),pattern = 'pch')+
    scale_pattern_shape_manual(values = c(NA, 16 ))+
   scale_fill_grey(start=0.4)+
   scale_color_manual(values=c("transparent",                           
                             "firebrick2"#,
                             ))+
    mytheme+
    theme(axis.text.x = element_text(angle = 90,vjust=-0)) +
   facet_grid(parameter~scenario)
pranks_mape
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/ranks/base_mape_rank.png", plot=pranks_mape)


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
     geom_bar_pattern(stat="identity", size=1.3, aes(pattern_shape=scnmodmatch_sd2),pattern = 'pch')+
    scale_pattern_shape_manual(values = c(NA, 16))+
   scale_fill_grey(start=0.4)+
   scale_color_manual(values=c("transparent",                           
                             "firebrick2"))+
    mytheme+
    theme(axis.text.x = element_text(angle = 90,vjust=-0)) +
   facet_grid(parameter~scenario)
pranks_sd2
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/ranks/base_sd2_rank.png", plot=pranks_sd2)



#======================================================================
#aggregate averages by 
aggrankbyscn<-aggregate(ranksbaseall$average.rank, list(parameter=ranksbaseall$parameter,
                                                       paramch=ranksbaseall$paramch, 
                                                       model=ranksbaseall$model,
                                                       modtype=ranksbaseall$modtype), mean)

aggrankbyscn$scnmodmatch<-FALSE
aggrankbyscn$scnmodmatch[aggrankbyscn$paramch %in% c("stationary")&
                        aggrankbyscn$model%in%c("simple","autocorr")]<-TRUE

aggrankbyscn$scnmodmatch[aggrankbyscn$paramch %in% c("a")&
                        aggrankbyscn$model%in%c("rwa","hmma")]<-TRUE

aggrankbyscn$scnmodmatch[aggrankbyscn$paramch %in% c("b")&
                        aggrankbyscn$model%in%c("rwb","hmmb")]<-TRUE

aggrankbyscn$scnmodmatch[aggrankbyscn$paramch %in% c("both")&
                        aggrankbyscn$model%in%c("rwab","hmmab")]<-TRUE

head(aggrankbyscn)



paramchs<-unique(aggrankbyscn$paramch)

aggrankbyscn$scnmodmatch2<-"0"

for(i in seq_along(paramchs)){
    for(u in seq_along(paras)){

    tmp<-min(aggrankbyscn[aggrankbyscn$paramch==paramchs[i]&
                        aggrankbyscn$parameter==paras[u],"x"])

     aggrankbyscn$scnmodmatch2[ aggrankbyscn$paramch==paramchs[i]&
                        aggrankbyscn$parameter==paras[u]&
                        aggrankbyscn$x==tmp]<-
                        "1"
  }
}



aggscn_base_rank<-ggplot(data=aggrankbyscn, aes(x=model, y=x,fill=modtype,color=scnmodmatch ))+
    geom_bar(stat="identity", size=1.3)+
    geom_bar_pattern(stat="identity", size=1.3, aes(pattern_shape=scnmodmatch2),pattern = 'pch')+
    scale_pattern_shape_manual(values = c(NA, 16))+
   scale_fill_grey(start=0.4)+
   scale_color_manual(values=c("transparent",
                             "darkred"))+
    mytheme+
    theme(axis.text.x = element_text(angle = 90,vjust=-0)) +
   facet_grid(parameter~paramch)
#Does not look representative of split by scenario
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/ranks/base_rank_aggscn.png", plot=aggscn_base_rank)




#routines to visualize smulation estimation Results
#Catarina Wor
#November 2022
#============================================




library(ggplot2)
library(gridExtra)
library(dplyr)
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
#base case
#read in data
simPar <- read.csv("data/generic/SimPars.csv")

## Store relevant object names to help run simulation 
scenNames <- unique(simPar$scenario)

res1<-readRDS(file = "outs/simest/generic/resbase1.rds")
res2<-readRDS(file = "outs/simest/generic/resbase2.rds")


restmb<-rbind(res1,res2)


#

resstan1<-readRDS(file = "outs/simest/generic/resstan1.rds")
resstan2<-readRDS(file = "outs/simest/generic/resstan2.rds")
resstan<-rbind(resstan1,resstan2)

head(resstan)
head(resstan[resstan$scenario=="decLinearProd"&resstan$model=="simple"&resstan$parameter=="smsy",],10)

#resstan<-readRDS(file = "outs/simest/generic/resstan.rds")

res<-rbind(restmb,resstan)

res$parameter[res$parameter=="Smax"]<-"smax"

resparam<-res[res$parameter%in%c("alpha","smax","sigma","smsy","sgen","umsy"),]



convstat<-aggregate(resparam$convergence,
    list(scenario=resparam$scenario,
        model=resparam$model,
        method=resparam$method,
        iteration=resparam$iteration),
    function(x){sum(x)})
convstatMLE<-convstat[convstat$x==0&convstat$method=="MLE",]
convstatMCMC<-convstat[convstat$x==0&convstat$method=="MCMC",]



allconv<-inner_join(convstatMLE[,-3], convstatMCMC[,-3])


convsum<-aggregate(allconv$iteration,
    list(model=allconv$model,scenario=allconv$scenario),
    function(x){length(unique(x))})

convsumMLE<-aggregate(convstatMLE$iteration,
    list(model=convstatMLE$model,scenario=convstatMLE$scenario),
    function(x){length(unique(x))})


concsumMCMC<-aggregate(convstatMCMC$iteration,
    list(model=convstatMCMC$model,scenario=convstatMCMC$scenario),
    function(x){length(unique(x))})


convsnc<-as.numeric(rownames(convsum))

conv_iter<-aggregate(allconv$iteration,
    list(model=allconv$model,scenario=allconv$scenario),
    function(x){(unique(x))})

resl<-list()
for(i in seq_along(convsnc)){

    sel<-conv_iter[convsnc[i],]
    resl[[i]]<-resparam %>% filter(model==sel$model&
                            scenario==sel$scenario&
                            iteration%in%sel$x[[1]])
    
}

resparam<-as.data.frame(data.table::rbindlist(resl))

#resparam<-resparam[resparam$convergence==0,]






df<-reshape2::melt(resparam, id.vars=c("parameter","iteration","scenario","method","model","by", "convergence","pbias","bias"))

#df_alpha<-df[df$parameter%in%c("alpha"),]
df$col<-factor(df$variable,levels=c("est","sim"))

#df_pbias<-df_mcmc[df_mcmc$variable%in%c("est"),]
df_pbias<-df[df$variable%in%c("mode"),]


df_pbias<-df_pbias[df_pbias$parameter!="sigma"&!is.na(df_pbias$bias),]



df_biasq<-aggregate(df_pbias$bias,list(parameter=df_pbias$parameter,
    scenario=df_pbias$scenario,
            method =df_pbias$method,
            model =df_pbias$model,
            by =df_pbias$by),function (x) quantile(x, c(0.1, 0.5, 0.9)))
head(df_biasq)
df_biasq$model<-factor(df_biasq$model, levels=c("simple", 
       "autocorr",
       "rwa",
       "hmma", 
       "rwb",
       "hmmb",
       "rwab",    
        "hmmab"))
df_biasq<-do.call(data.frame, df_biasq)

head(df_biasq)
df_biasq$paramch<-"stationary"
df_biasq$paramch[df_biasq$scenario%in%c("decLinearCap" , "regimeCap", 
 "shiftCap")]<-"b"
df_biasq$paramch[df_biasq$scenario%in%c( "decLinearProd"  , "regimeProd" , "shiftProd",
  "sineProd")]<-"a"
df_biasq$paramch[df_biasq$scenario%in%c("decLinearProdshiftCap" ,    "regimeProdCap"  )]<-"both"
   
df_biasq$type<-"none"
df_biasq$type[df_biasq$scenario%in%c("decLinearProd" ,"decLinearCap" , "sineProd")]<-"trend"
df_biasq$type[df_biasq$scenario%in%c( "regimeCap"  , "regimeProd" , "regimeProdCap", "shiftCap",
  "shiftProd") ]<-"regime"
df_biasq$type[df_biasq$scenario%in%c("decLinearProdshiftCap" )]<-"combo"
       

df_biasq<-df_biasq[df_biasq$method=="MLE",]


df_alpha_biasq<-df_biasq[df_biasq$parameter=="alpha",]
basea<-ggplot(df_alpha_biasq) + 
coord_cartesian(ylim = c(-1,1))+ 
geom_ribbon(aes(x=by,ymin =  x.10., ymax=x.90.,fill=type,group=scenario),alpha=0.2)+
geom_line(aes(x=by,y= x.50.,color=type,group=scenario),linewidth=1.2)+
       scale_color_viridis_d("type of trend:",begin=.1, end=.8) +
       scale_fill_viridis_d("type of trend:",begin=.1, end=.8) +
       geom_hline(yintercept=0, color="black",alpha=.7,linewidth=1.2) +
       mytheme+ 
       scale_y_continuous(name = "MLE bias in log(alpha)", sec.axis = sec_axis(~., name = "varying parameters")) +
       theme(axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank()
       )+
       #ylab("MLE bias in log(alpha)") +
       xlab("year") +

       facet_grid(paramch~model, scales="free_y")
basea  
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/trends/base/basea_MLE.png")



df_smax_biasq<-df_biasq[df_biasq$parameter=="smax",]
baseb<-ggplot(df_smax_biasq) + 
coord_cartesian(ylim = c(-100000,150000))+ 
geom_ribbon(aes(x=by,ymin =  x.10., ymax=x.90.,fill=type,group=scenario),alpha=0.2)+
geom_line(aes(x=by,y= x.50.,color=type,group=scenario),linewidth=1.2)+
       scale_color_viridis_d(begin=.1, end=.8) +
       scale_fill_viridis_d(begin=.1, end=.8) +
       geom_hline(yintercept=0, color="black",alpha=.7,linewidth=1.2) +
       mytheme+ 
       scale_y_continuous(name = "MLE bias in Smax", sec.axis = sec_axis(~., name = "varying parameters")) +
       theme(axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank()
       )+
       xlab("year") +
       facet_grid(paramch~model, scales="free_y")
baseb
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/trends/base/baseSmax_MLE.png")



df_smsy_biasq<-df_biasq[df_biasq$parameter=="smsy",]
head(df_smsy_biasq)
basesmsy<-ggplot(df_smsy_biasq) + 
coord_cartesian(ylim = c(-80000,80000))+ 
geom_ribbon(aes(x=by,ymin =  x.10., ymax=x.90.,fill=type,group=scenario),alpha=0.2)+
geom_line(aes(x=by,y= x.50.,color=type,group=scenario),linewidth=1.2)+
       scale_color_viridis_d(begin=.1, end=.8) +
       scale_fill_viridis_d(begin=.1, end=.8) +
       geom_hline(yintercept=0, color="black",alpha=.7,linewidth=1.2) +
       mytheme+ 
      scale_y_continuous(name = "MLE bias in Smsy", sec.axis = sec_axis(~., name = "varying parameters")) +
       theme(axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank()
       )+
       xlab("year") +
       facet_grid(paramch~model, scales="free_y")
basesmsy
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/trends/base/basesmsy_MLE.png",plot=basesmsy)



df_sgen_biasq<-df_biasq[df_biasq$parameter=="sgen",]
basesgen<-ggplot(df_sgen_biasq) + 
coord_cartesian(ylim = c(-50000,50000))+ 
geom_ribbon(aes(x=by,ymin =  x.10., ymax=x.90.,fill=type,group=scenario),alpha=0.2)+
geom_line(aes(x=by,y= x.50.,color=type,group=scenario),linewidth=1.2)+
       scale_color_viridis_d(begin=.1, end=.8) +
       scale_fill_viridis_d(begin=.1, end=.8) +
       geom_hline(yintercept=0, color="black",alpha=.7,linewidth=1.2) +
       mytheme+ 
      scale_y_continuous(name = "MLE bias in Sgen", sec.axis = sec_axis(~., name = "varying parameters")) +
       theme(axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank()
       )+
       xlab("year") +
       facet_grid(paramch~model, scales="free_y")
basesgen
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/trends/base/basesgen_MLE.png", plot=basesgen)





df_umsy_biasq<-df_biasq[df_biasq$parameter=="umsy",]
head(df_umsy_biasq)
baseumsy<-ggplot(df_umsy_biasq) + 
coord_cartesian(ylim = c(-.35,.35))+ 
geom_ribbon(aes(x=by,ymin =  x.10., ymax=x.90.,fill=type,group=scenario),alpha=0.2)+
geom_line(aes(x=by,y= x.50.,color=type,group=scenario),linewidth=1.2)+
       scale_color_viridis_d(begin=.1, end=.8) +
       scale_fill_viridis_d(begin=.1, end=.8) +
       geom_hline(yintercept=0, color="black",alpha=.7,linewidth=1.2) +
       mytheme+ scale_y_continuous(name = "MLE bias in Sgen", sec.axis = sec_axis(~., name = "varying parameters")) +
       theme(axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank()
       )+
       xlab("year") +
       facet_grid(paramch~model, scales="free_y")
baseumsy
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/trends/base/baseumsy_MLE.png", plot=baseumsy)


#Same with MCMC
#need to re-work!
df_mcmc<-df[df$method=="MCMC",]
unique(df_mcmc$variable)
df_pbias_mcmc<-df_mcmc[df_mcmc$variable%in%c("mode"),]
head(df_pbias_mcmc)

df_pbias_mcmc<-df_pbias_mcmc[df_pbias_mcmc$parameter!="sigma"&!is.na(df_pbias_mcmc$bias),]




df_biasq_mcmc<-aggregate(df_pbias_mcmc$bias,list(parameter=df_pbias_mcmc$parameter,
    scenario=df_pbias_mcmc$scenario,
            method =df_pbias_mcmc$method,
            model =df_pbias_mcmc$model,
            by =df_pbias_mcmc$by),function (x) quantile(x, c(0.1, 0.5, 0.9)))
head(df_biasq_mcmc)
df_biasq_mcmc$model<-factor(df_biasq$model, levels=c("simple", 
       "autocorr",
       "rwa",
       "hmma", 
       "rwb",
        "hmmb",  
       "rwab", 
        "hmmab"))

df_biasq_mcmc<-do.call(data.frame, df_biasq_mcmc)
head(df_biasq_mcmc)



df_biasq_mcmc$paramch<-"stationary"
df_biasq_mcmc$paramch[df_biasq_mcmc$scenario%in%c("decLinearCap" , "regimeCap", 
 "shiftCap")]<-"b"
df_biasq_mcmc$paramch[df_biasq_mcmc$scenario%in%c( "decLinearProd"  , "regimeProd" , "shiftProd",
  "sineProd")]<-"a"
df_biasq_mcmc$paramch[df_biasq_mcmc$scenario%in%c("decLinearProdshiftCap" ,    "regimeProdCap"  )]<-"both"
   
df_biasq_mcmc$type<-"none"
df_biasq_mcmc$type[df_biasq_mcmc$scenario%in%c("decLinearProd" ,"decLinearCap" , "sineProd")]<-"trend"
df_biasq_mcmc$type[df_biasq_mcmc$scenario%in%c( "regimeCap"  , "regimeProd" , "regimeProdCap", "shiftCap",
  "shiftProd") ]<-"regime"
df_biasq_mcmc$type[df_biasq_mcmc$scenario%in%c("decLinearProdshiftCap" )]<-"combo"

df_alpha_biasq_mcmc<-df_biasq_mcmc[df_biasq_mcmc$parameter=="alpha",]
head(df_alpha_biasq_mcmc)
basea_mcmc<-ggplot(df_alpha_biasq_mcmc) + 
coord_cartesian(ylim = c(-1,1))+ 
geom_ribbon(aes(x=by,ymin =  x.10., ymax=x.90.,fill=type,group=scenario),alpha=0.2)+
geom_line(aes(x=by,y= x.50.,color=type,group=scenario),linewidth=1.2)+
       scale_color_viridis_d("type of trend:",begin=.1, end=.8) +
       scale_fill_viridis_d("type of trend:",begin=.1, end=.8) +
       geom_hline(yintercept=0, color="black",alpha=.7,linewidth=1.2) +
       mytheme+ 
       scale_y_continuous(name = "MCMC bias in log(alpha)", sec.axis = sec_axis(~., name = "varying parameters")) +
       theme(axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank()
       )+
       ylab("MCMC bias in log(alpha)") +
       xlab("year") +

       facet_grid(paramch~model, scales="free_y")
basea_mcmc  
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/trends/base/basea_MCMC.png")






df_smax_biasq_mcmc<-df_biasq_mcmc[df_biasq_mcmc$parameter=="smax",]
head(df_smax_biasq_mcmc)
baseb_mcmc<-ggplot(df_smax_biasq_mcmc) + 
coord_cartesian(ylim = c(-100000,150000))+ 
geom_ribbon(aes(x=by,ymin =  x.10., ymax=x.90.,fill=type,group=scenario),alpha=0.2)+
geom_line(aes(x=by,y= x.50.,color=type,group=scenario),linewidth=1.2)+
       scale_color_viridis_d(begin=.1, end=.8) +
       scale_fill_viridis_d(begin=.1, end=.8) +
       geom_hline(yintercept=0, color="black",alpha=.7,linewidth=1.2) +
       mytheme+ 
       scale_y_continuous(name = "MCMC bias in Smax", sec.axis = sec_axis(~., name = "varying parameters")) +
       theme(axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank()
       )+
       xlab("year") +
       facet_grid(paramch~model, scales="free_y")
baseb_mcmc
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/trends/base/baseSmax_MCMC.png")


head(df_smsy_biasq_mcmc)
df_smsy_biasq_mcmc<-df_biasq_mcmc[df_biasq_mcmc$parameter=="smsy",]
basesmsy_mcmc<-ggplot(df_smsy_biasq_mcmc) + 
coord_cartesian(ylim = c(-60000,70000))+ 
geom_ribbon(aes(x=by,ymin =  x.10., ymax=x.90.,fill=type,group=scenario),alpha=0.2)+
geom_line(aes(x=by,y= x.50.,color=type,group=scenario),linewidth=1.2)+
       scale_color_viridis_d(begin=.1, end=.8) +
       scale_fill_viridis_d(begin=.1, end=.8) +
       geom_hline(yintercept=0, color="black",alpha=.7,linewidth=1.2) +
       mytheme+ 
      scale_y_continuous(name = "MLE bias in Smsy", sec.axis = sec_axis(~., name = "varying parameters")) +
       theme(axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank()
       )+
       xlab("year") +
       facet_grid(paramch~model, scales="free_y")
basesmsy_mcmc
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/trends/base/basesmsy_MCMC.png",plot=basesmsy_mcmc)



df_sgen_biasq_mcmc<-df_biasq_mcmc[df_biasq_mcmc$parameter=="sgen",]
basesgen_mcmc<-ggplot(df_sgen_biasq_mcmc) + 
coord_cartesian(ylim = c(-25000,40000))+ 
geom_ribbon(aes(x=by,ymin =  x.10., ymax=x.90.,fill=type,group=scenario),alpha=0.2)+
geom_line(aes(x=by,y= x.50.,color=type,group=scenario),linewidth=1.2)+
       scale_color_viridis_d(begin=.1, end=.8) +
       scale_fill_viridis_d(begin=.1, end=.8) +
       geom_hline(yintercept=0, color="black",alpha=.7,linewidth=1.2) +
       mytheme+ 
      scale_y_continuous(name = "MLE bias in Sgen", sec.axis = sec_axis(~., name = "varying parameters")) +
       theme(axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank()
       )+
       xlab("year") +
       facet_grid(paramch~model, scales="free_y")
basesgen_mcmc
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/trends/base/basesgen_MCMC.png", plot=basesgen)


#=========================================================================================
#Percent bias
df_pbiasq<-aggregate(df_pbias$pbias,list(parameter=df_pbias$parameter,
    scenario=df_pbias$scenario,
            method =df_pbias$method,
            model =df_pbias$model,
            by =df_pbias$by),function (x) quantile(x, c(0.1, 0.5, 0.9)))

df_pbiasq$model<-factor(df_pbiasq$model, levels=c("simple", 
       "autocorr",
       "rwa",
       "rwb",
       "rwab",
       "hmma",  
        "hmmb",    
        "hmmab"))
df_pbiasq<-do.call(data.frame, df_pbiasq)

unique(df_biasq$scenario)

df_pbiasq$paramch<-"stationary"
df_pbiasq$paramch[df_pbiasq$scenario%in%c("decLinearCap" , "regimeCap", "shiftCap")]<-"b"
df_pbiasq$paramch[df_pbiasq$scenario%in%c( "decLinearProd"  , "regimeProd" , "shiftProd", "sineProd")]<-"a"
df_pbiasq$paramch[df_pbiasq$scenario%in%c("decLinearProdshiftCap" ,    "regimeProdCap"  )]<-"both"
   
df_pbiasq$type<-"none"
df_pbiasq$type[df_pbiasq$scenario%in%c("decLinearProd" ,"decLinearCap" , "sineProd")]<-"trend"
df_pbiasq$type[df_pbiasq$scenario%in%c( "regimeCap"  , "regimeProd" , "regimeProdCap", "shiftCap","shiftProd") ]<-"regime"
df_pbiasq$type[df_pbiasq$scenario%in%c("decLinearProdshiftCap" )]<-"combo"
       



df_alpha_pbiasq<-df_pbiasq[df_pbiasq$parameter=="alpha",]
ggplot(df_alpha_pbiasq) + 
coord_cartesian(ylim = c(-100,100))+ 
geom_ribbon(aes(x=by,ymin =  x.10., ymax=x.90.,fill=type,group=scenario),alpha=0.2)+
geom_line(aes(x=by,y= x.50.,color=type,group=scenario),linewidth=1.2)+
       scale_color_viridis_d(begin=.1, end=.8) +
       scale_fill_viridis_d(begin=.1, end=.8) +
       geom_hline(yintercept=0, color="black",alpha=.7,linewidth=1.2) +
       mytheme+ 
       ylab("Bias in log(alpha)") +
       xlab("year") +
       facet_grid(paramch~model, scales="free_y")




df_umsy_pbiasq<-df_pbiasq[df_pbiasq$parameter=="umsy",]
ggplot(df_umsy_pbiasq) + 
coord_cartesian(ylim = c(-100,100))+ 
geom_ribbon(aes(x=by,ymin =  x.10., ymax=x.90.,fill=type,group=scenario),alpha=0.2)+
geom_line(aes(x=by,y= x.50.,color=type,group=scenario),linewidth=1.2)+
       scale_color_viridis_d(begin=.1, end=.8) +
       scale_fill_viridis_d(begin=.1, end=.8) +
       geom_hline(yintercept=0, color="black",alpha=.7,linewidth=1.2) +
       mytheme+ 
       ylab("Bias in log(alpha)") +
       xlab("year") +
       facet_grid(paramch~model, scales="free_y")


#========================================================================================================
#sensitivity a scenario
#read in data
simPar <- read.csv("data/sensitivity/SimPars.csv")

## Store relevant object names to help run simulation 
scenNames <- unique(simPar$scenario)

resa1<-readRDS(file = "outs/simest/sensitivity/res_a1.rds")
resa2<-readRDS(file = "outs/simest/sensitivity/res_a2.rds")
restmb_a<-rbind(resa1,resa2)


resstana1<-readRDS(file = "outs/simest/sensitivity/resstan_a1.rds")
resstana2<-readRDS(file = "outs/simest/sensitivity/resstan_a2.rds")
resstan_a<-rbind(resstana1,resstana2)

res_a<-rbind(restmb_a,resstan_a)

res_a$parameter[res_a$parameter=="Smax"]<-"smax"
unique(res_a$parameter)
resparam<-res_a[res_a$parameter%in%c("alpha","smax","sigma","smsy","sgen","umsy"),]


convstat<-aggregate(resparam$convergence,
    list(scenario=resparam$scenario,
        model=resparam$model,
        method=resparam$method,
        iteration=resparam$iteration),
    function(x){sum(x)})


allconv<-inner_join(convstatMLE[,-3], convstatMCMC[,-3])


convsum<-aggregate(allconv$iteration,
    list(model=allconv$model,scenario=allconv$scenario),
    function(x){length(unique(x))})

convsnc<-as.numeric(rownames(convsum)) 

conv_iter<-aggregate(allconv$iteration,
    list(model=allconv$model,scenario=allconv$scenario),
    function(x){(unique(x))})


head(resparam )
resl<-list()
for(i in seq_along(convsnc)){

    sel<-conv_iter[convsnc[i],]
    resl[[i]]<-resparam %>% filter(model==sel$model&
                            scenario==sel$scenario&
                            iteration%in%sel$x[[1]])
    
}

resparam<-as.data.frame(data.table::rbindlist(resl))

head(resparam)

exp(unique(resparam$sim[resparam$parameter=="alpha"&
    resparam$scenario=="regimeProd1" ]))

exp(unique(resparam$sim[resparam$parameter=="alpha"&
    resparam$scenario=="regimeProd2" ]))


df<-reshape2::melt(resparam, id.vars=c("parameter","iteration","scenario","method","model","by", "convergence","pbias","bias"))

#df_alpha<-df[df$parameter%in%c("alpha"),]
df$col<-factor(df$variable,levels=c("est","sim"))

#df_pbias<-df_mcmc[df_mcmc$variable%in%c("est"),]
df_pbias<-df[df$variable%in%c("mode"),]


df_pbias<-df_pbias[df_pbias$parameter!="sigma"&!is.na(df_pbias$bias),]



df_biasq<-aggregate(df_pbias$bias,list(parameter=df_pbias$parameter,
    scenario=df_pbias$scenario,
            method =df_pbias$method,
            model =df_pbias$model,
            by =df_pbias$by),function (x) quantile(x, c(0.1, 0.5, 0.9)))
head(df_biasq)
df_biasq$model<-factor(df_biasq$model, levels=c("simple", 
       "autocorr",
       "rwa",
       "hmma", 
       "rwb",
       "hmmb",
       "rwab",    
        "hmmab"))
df_biasq<-do.call(data.frame, df_biasq)





df_biasq$paramch<-"a"
 
df_biasq$type<-"regime"
df_biasq$type[df_biasq$scenario%in%c("trendLinearProd10", "trendLinearProd2", 
 "trendLinearProd5", "trendLinearProd7","trendLinearProd1")]<-"trend"

df_biasq$magnitude_ch<-"3.7 -> 1.03"  
df_biasq$magnitude_ch[df_biasq$scenario%in%c("trendLinearProd10", "regimeProd10")]<-"3.7 -> 10"  
df_biasq$magnitude_ch[df_biasq$scenario%in%c("trendLinearProd2", "regimeProd2")]<-"3.7 -> 2"   
df_biasq$magnitude_ch[df_biasq$scenario%in%c("trendLinearProd5", "regimeProd5")]<-"3.7 -> 5"  
df_biasq$magnitude_ch[df_biasq$scenario%in%c("trendLinearProd7", "regimeProd7")]<-"3.7 -> 7"               
df_biasq$magnitude_ch<-factor(df_biasq$magnitude_ch,levels=c("3.7 -> 1.03", 
       "3.7 -> 2",
       "3.7 -> 5",
       "3.7 -> 7", 
       "3.7 -> 10"))



df_biasq_MCMC<-df_biasq[df_biasq$method=="MCMC",]


df_alpha_biasq_MCMC<-df_biasq_MCMC[df_biasq_MCMC$parameter=="alpha",]
sensa_alpha<-ggplot(df_alpha_biasq_MCMC) + 
coord_cartesian(ylim = c(-2,1))+ 
geom_ribbon(aes(x=by,ymin =  x.10., ymax=x.90.,fill=magnitude_ch,group=scenario),alpha=0.2)+
geom_line(aes(x=by,y= x.50.,color=magnitude_ch,group=scenario),linewidth=1.2)+
       scale_color_viridis_d("magnitude of change:",begin=.1, end=.8) +
       scale_fill_viridis_d("magnitude of change:",begin=.1, end=.8) +
       geom_hline(yintercept=0, color="black",alpha=.7,linewidth=1.2) +
       mytheme+ 
       scale_y_continuous(name = "MLE bias in log(alpha)", sec.axis = sec_axis(~., name = "varying parameters")) +
       theme(axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank()
       )+
       
       xlab("year") +
       facet_grid(type~model, scales="free_y")
sensa_alpha  
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/trends/sens_a/bias_a_MCMC.png", plot=sensa_alpha )



df_smax_biasq_MCMC<-df_biasq_MCMC[df_biasq_MCMC$parameter=="smax",]
sensa_smax<-ggplot(df_smax_biasq_MCMC) + 
coord_cartesian(ylim = c(-100000,150000))+ 
geom_ribbon(aes(x=by,ymin =  x.10., ymax=x.90.,fill=magnitude_ch,group=scenario),alpha=0.2)+
geom_line(aes(x=by,y= x.50.,color=magnitude_ch,group=scenario),linewidth=1.2)+
       scale_color_viridis_d("magnitude of change:",begin=.1, end=.8) +
       scale_fill_viridis_d("magnitude of change:",begin=.1, end=.8) +
       geom_hline(yintercept=0, color="black",alpha=.7,linewidth=1.2) +
       mytheme+ 
       scale_y_continuous(name = "MLE bias in Smax", sec.axis = sec_axis(~., name = "varying parameters")) +
       theme(axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank()
       )+
       #ylab("MLE bias in log(alpha)") +
       xlab("year") +
       facet_grid(type~model, scales="free_y")
sensa_smax  
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/trends/sens_a/bias_smax_MCMC.png", plot=sensa_smax   )



df_smsy_biasq_MCMC<-df_biasq_MCMC[df_biasq_MCMC$parameter=="smsy",]
sensa_smsy<-ggplot(df_smsy_biasq_MCMC) + 
coord_cartesian(ylim = c(-70000,70000))+ 
geom_ribbon(aes(x=by,ymin =  x.10., ymax=x.90.,fill=magnitude_ch,group=scenario),alpha=0.2)+
geom_line(aes(x=by,y= x.50.,color=magnitude_ch,group=scenario),linewidth=1.2)+
       scale_color_viridis_d("magnitude of change:",begin=.1, end=.8) +
       scale_fill_viridis_d("magnitude of change:",begin=.1, end=.8) +
       geom_hline(yintercept=0, color="black",alpha=.7,linewidth=1.2) +
       mytheme+ 
       scale_y_continuous(name = "MLE bias in Smsy", sec.axis = sec_axis(~., name = "varying parameters")) +
       theme(axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank()
       )+
       #ylab("MLE bias in log(alpha)") +
       xlab("year") +
       facet_grid(type~model, scales="free_y")
sensa_smsy  
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/trends/sens_a/bias_smsy_MCMC.png", plot=sensa_smsy   )




#=================================================

df_biasq_MLE<-df_biasq[df_biasq$method=="MLE",]
head(df_biasq_MLE)
df_alpha_biasq_MLE<-df_biasq_MLE[df_biasq_MLE$parameter=="alpha",]
head(df_alpha_biasq_MLE)
sensa_alpha_MLE<-ggplot(df_alpha_biasq_MLE) + 
coord_cartesian(ylim = c(-2,1))+ 
geom_ribbon(aes(x=by,ymin =  x.10., ymax=x.90.,fill=magnitude_ch,group=scenario),alpha=0.2)+
geom_line(aes(x=by,y= x.50.,color=magnitude_ch,group=scenario),linewidth=1.2)+
       scale_color_viridis_d("magnitude of change:",begin=.1, end=.8) +
       scale_fill_viridis_d("magnitude of change:",begin=.1, end=.8) +
       geom_hline(yintercept=0, color="black",alpha=.7,linewidth=1.2) +
       mytheme+ 
       scale_y_continuous(name = "MLE bias in log(alpha)", sec.axis = sec_axis(~., name = "varying parameters")) +
       theme(axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank()
       )+
       
       xlab("year") +
       facet_grid(type~model, scales="free_y")
sensa_alpha_MLE  
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/trends/sens_a/bias_a_MLE.png", plot=sensa_alpha_MLE )



df_smax_biasq_MLE<-df_biasq_MLE[df_biasq_MLE$parameter=="smax",]
sensa_smax_MLE<- ggplot(df_smax_biasq_MLE) + 
coord_cartesian(ylim = c(-120000,300000))+ 
geom_ribbon(aes(x=by,ymin =  x.10., ymax=x.90.,fill=magnitude_ch,group=scenario),alpha=0.2)+
geom_line(aes(x=by,y= x.50.,color=magnitude_ch,group=scenario),linewidth=1.2)+
       scale_color_viridis_d("magnitude of change:",begin=.1, end=.8) +
       scale_fill_viridis_d("magnitude of change:",begin=.1, end=.8) +
       geom_hline(yintercept=0, color="black",alpha=.7,linewidth=1.2) +
       mytheme+ 
       scale_y_continuous(name = "MLE bias in Smax", sec.axis = sec_axis(~., name = "varying parameters")) +
       theme(axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank()
       )+
       #ylab("MLE bias in log(alpha)") +
       xlab("year") +
       facet_grid(type~model, scales="free_y")
sensa_smax_MLE  
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/trends/sens_a/bias_smax_MLE.png", plot=sensa_smax_MLE  )


df_smsy_biasq_MLE<-df_biasq_MLE[df_biasq_MLE$parameter=="smsy",]
sensa_smsy_MLE <-ggplot(df_smsy_biasq_MLE) +  
coord_cartesian(ylim = c(-50000,100000))+ 
geom_ribbon(aes(x=by,ymin =  x.10., ymax=x.90.,fill=magnitude_ch,group=scenario),alpha=0.2)+
geom_line(aes(x=by,y= x.50.,color=magnitude_ch,group=scenario),linewidth=1.2)+
       scale_color_viridis_d(begin=.1, end=.8) +
       scale_fill_viridis_d(begin=.1, end=.8) +
       geom_hline(yintercept=0, color="black",alpha=.7,linewidth=1.2) +
       mytheme+ 
       ylab("Bias in Smsy") +
       xlab("year") +
       facet_grid(type~model, scales="free_y")
sensa_smsy_MLE
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/trends/sens_a/bias_smsy_MLE.png", plot=sensa_smsy_MLE  )

#========================================================================================================
#sensitivity smax scenario






#========================================================================================================
#sensitivity sigma scenario



#========================================================================================================
#sensitivity ER scenario
#read in data
simPar <- read.csv("data/genericER/SimPars_ER.csv")

## Store relevant object names to help run simulation 
scenNames <- unique(simPar$scenario)

res_er<-readRDS(file = "outs/simest/genericER/resbase_ER.rds")
#res_er_stan<-readRDS(file = "outs/simest/sensitivity/resstan_a.rds")


#res_er<-rbind(res_er,res_er_stan)
#res_er$parameter[res_er$parameter=="smax"]<-"Smax"
#unique(res_er$parameter)
resparam<-res_er[res_er$parameter%in%c("alpha","Smax","sigma","smsy","sgen","umsy"),]

resparam$bias<-(resparam$pbias/100)*resparam$sim
summary(resparam[resparam$parameter%in%c("Smax"),]$sim)


df<-reshape2::melt(resparam, id.vars=c("parameter","iteration","scenario","method","model","by", "convergence","pbias","bias"))

#df_alpha<-df[df$parameter%in%c("alpha"),]
df$col<-factor(df$variable,levels=c("est","sim"))

df_mle<-df[df$method=="MLE",]
df_mcmc<-df[df$method=="MCMC",]


df_pbias<-df_mle[df_mle$variable%in%c("est"),]


df_pbias<-df_pbias[df_pbias$parameter!="sigma"&df_pbias$convergence<1&!is.na(df_pbias$bias),]


df_biasq<-aggregate(df_pbias$bias,list(parameter=df_pbias$parameter,
    scenario=df_pbias$scenario,
            method =df_pbias$method,
            model =df_pbias$model,
            by =df_pbias$by),function (x) quantile(x, c(0.1, 0.5, 0.9)))

df_biasq$model<-factor(df_biasq$model, levels=c("simple", 
       "autocorr",
       "rwa",
       "rwb",
       "rwab",
       "hmma",  
        "hmmb",    
        "hmmab"))
df_biasq<-do.call(data.frame, df_biasq)

unique(df_biasq$scenario)

[1] "decLinearProd_highERLowError"  "decLinearProd_ShiftERLowError"
 [3] "highERHighError"               "highERLowError"               
 [5] "incLinearProd_highERLowError"  "incLinearProd_ShiftERLowError"
 [7] "lowERHighError"                "lowERLowError"                
 [9] "ShiftERHighError"              "ShiftERLowError"              
[11] "trendERHighError"              "trendERLowError"              
> 

df_biasq$er<-"high"
df_biasq$er[df_biasq$scenario%in%c("lowERHighError" ,"lowERLowError" )]<-"low"
df_biasq$er[df_biasq$scenario%in%c("decLinearProd_ShiftERLowError","incLinearProd_ShiftERLowError","ShiftERHighError", 
    "ShiftERLowError")]<-"shift"
df_biasq$er[df_biasq$scenario%in%c("trendERHighError","trendERLowError")]<-"trend"

df_biasq$ervar<-"high"
df_biasq$ervar[df_biasq$scenario%in%c("lowERLowError", "decLinearProd_ShiftERLowError","highERLowError",
"incLinearProd_highERLowError",  "incLinearProd_ShiftERLowError","lowERLowError","ShiftERLowError","trendERLowError"  )]<-"low"


df_biasq$prod<-"stationary"
df_biasq$prod[df_biasq$scenario%in%c("decLinearProd_highERLowError",  "decLinearProd_ShiftERLowError",
 "incLinearProd_highERLowError",  "incLinearProd_ShiftERLowError")]<-"trend in a"

library(ggnewscale)
library(dplyr)
library(grid)
library(gridExtra)

df_alpha_biasq<-df_biasq[df_biasq$parameter=="alpha",]
head(df_alpha_biasq)

ggplot(df_alpha_biasq) + 
coord_cartesian(ylim = c(-1,1))+ 
geom_ribbon(aes(x=by,ymin =  x.10., ymax=x.90.,fill=prod,group=scenario),
     alpha=0.2)+
geom_line(aes(x=by,y= x.50.,color=prod,group=scenario),
     linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
    scale_fill_viridis_d(begin=.1, end=.8) +
#
       geom_hline(yintercept=0, color="black",alpha=.7,linewidth=1.2) +
       mytheme+ 
       ylab("Bias in log(alpha)") +
       xlab("year") +
       facet_grid(er~model, scales="free_y")


df_smsy_biasq<-df_biasq[df_biasq$parameter=="smsy",]
ggplot(df_smsy_biasq) +  
coord_cartesian(ylim = c(-50000,100000))+ 
geom_ribbon(aes(x=by,ymin =  x.10., ymax=x.90.,fill=prod,group=scenario),alpha=0.2)+
geom_line(aes(x=by,y= x.50.,color=prod,group=scenario),linewidth=1.2)+
       scale_color_viridis_d(begin=.1, end=.8) +
       scale_fill_viridis_d(begin=.1, end=.8) +
       geom_hline(yintercept=0, color="black",alpha=.7,linewidth=1.2) +
       mytheme+ 
       ylab("Bias in Smsy") +
       xlab("year") +
       facet_grid(er~model, scales="free_y")



unique(df_biasq$parameter)
df_smax_biasq<-df_biasq[df_biasq$parameter=="Smax",]
ggplot(df_smax_biasq) +  
coord_cartesian(ylim = c(-50000,100000))+ 
geom_ribbon(aes(x=by,ymin =  x.10., ymax=x.90.,fill=prod,group=scenario),alpha=0.2)+
geom_line(aes(x=by,y= x.50.,color=prod,group=scenario),linewidth=1.2)+
       scale_color_viridis_d(begin=.1, end=.8) +
       scale_fill_viridis_d(begin=.1, end=.8) +
       geom_hline(yintercept=0, color="black",alpha=.7,linewidth=1.2) +
       mytheme+ 
       ylab("Bias in Smsy") +
       xlab("year") +
       facet_grid(er~model, scales="free_y")


#routines to visualize smulation estimation Results
#Catarina Wor
#November 2022
#============================================




library(ggplot2)
library(gridExtra)
library(dplyr)
library(cowplot)
source("code/utils.R")
source("code/cluster_func_plots.R")
source("code/scale_inidividual_facet_y_axes.R")

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



#res<-readRDS(file = "outs/simest/generic/res.rds")
#restmb<-restmb[restmb$convergence==0,]

#

resstan1<-readRDS(file = "outs/simest/generic/resstan1.rds")
resstan2<-readRDS(file = "outs/simest/generic/resstan2.rds")
resstan<-rbind(resstan1,resstan2)

#resstan<-readRDS(file = "outs/simest/generic/resstan.rds")

res<-rbind(restmb,resstan)

res$parameter[res$parameter=="Smax"]<-"smax"

resparam<-res[res$parameter%in%c("alpha","smax","sigma","smsy","sgen","umsy"),]

#use only the variables that converged. Decided not to use it vecause of the various instances when 
#TMB rwb version did not converge for stationary scenarios. 
convstat<-aggregate(resparam$convergence,
    list(scenario=resparam$scenario,
        model=resparam$model,
        method=resparam$method,
        iteration=resparam$iteration),
    function(x){sum(x)})
convstatMLE<-convstat[convstat$x==0&convstat$method=="MLE",]
convstatMCMC<-convstat[convstat$x==0&convstat$method=="MCMC",]

head(convstatMCMC[,-3])

(convstat[convstat$scenario=="stationary",][1:20,])


(convstat[convstat$scenario=="stationary"&convstat$model=="rwb",])
resparam[resparam$scenario=="stationary"&resparam$model=="rwb"&resparam$iteration==1&resparam$method=="MLE",]

allconv<-inner_join(convstatMLE[,-3], convstatMCMC[,-3])

head(allconv)

convsum<-aggregate(allconv$iteration,
    list(model=allconv$model,scenario=allconv$scenario),
    function(x){length(unique(x))})

convsumMLE<-aggregate(convstatMLE$iteration,
    list(model=convstatMLE$model,scenario=convstatMLE$scenario),
    function(x){length(unique(x))})


convsumMCMC<-aggregate(convstatMCMC$iteration,
    list(model=convstatMCMC$model,scenario=convstatMCMC$scenario),
    function(x){length(unique(x))})


conv_iter<-aggregate(allconv$iteration,
    list(model=allconv$model,scenario=allconv$scenario),
    function(x){(unique(x))})


conv_iter_MLE<-aggregate(convstatMLE$iteration,
    list(model=convstatMLE$model,scenario=convstatMLE$scenario),
    function(x){(unique(x))})

conv_iter_MCMC<-aggregate(convstatMCMC$iteration,
    list(model=convstatMCMC$model,scenario=convstatMCMC$scenario),
    function(x){(unique(x))})

convsnc<-as.numeric(rownames(convsum))

resl<-list()
for(i in seq_along(convsnc)){

    sel<-conv_iter[convsnc[i],]
    resl[[i]]<-resparam %>% filter(model==sel$model&
                            scenario==sel$scenario&
                            iteration%in%sel$x[[1]])
    
}

resparam<-as.data.frame(data.table::rbindlist(resl))




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

df_biasq$model<-factor(df_biasq$model, levels=c("simple", 
       "autocorr",
       "rwa",
       "rwb",
       "rwab",
       "hmma",  
        "hmmb",    
        "hmmab"))
df_biasq<-do.call(data.frame, df_biasq)


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




#===========================================================
#MCMC results -all conv. 

resparam<-res[res$parameter%in%c("alpha","smax","sigma","smsy","sgen","umsy"),]
resparam<-resparam[resparam$convergence==0,]

df_amc<-reshape2::melt(resparam, id.vars=c("parameter","iteration","scenario","method","model","by", "convergence","pbias","bias"))

#df_alpha<-df[df$parameter%in%c("alpha"),]
df_amc$col<-factor(df_amc$variable,levels=c("est","sim"))

df_mcmc_all<-df_amc[df_amc$method=="MCMC",]

df_pbias_mcmc_all<-df_mcmc_all[df_mcmc_all$variable%in%c("mode"),]


df_pbias_mcmc_all<-df_pbias_mcmc_all[df_pbias_mcmc_all$parameter!="sigma"&!is.na(df_pbias_mcmc_all$bias),]


df_pbias_mcmc_all$paramch<-"stationary"
df_pbias_mcmc_all$paramch[df_pbias_mcmc_all$scenario%in%c("decLinearCap" , "regimeCap", 
 "shiftCap")]<-"b"
df_pbias_mcmc_all$paramch[df_pbias_mcmc_all$scenario%in%c( "decLinearProd"  , "regimeProd" , "shiftProd",
  "sineProd")]<-"a"
df_pbias_mcmc_all$paramch[df_pbias_mcmc_all$scenario%in%c("decLinearProdshiftCap" ,    "regimeProdCap"  )]<-"both"
   
df_pbias_mcmc_all$type<-"none"
df_pbias_mcmc_all$type[df_pbias_mcmc_all$scenario%in%c("decLinearProd" ,"decLinearCap" , "sineProd")]<-"trend"
df_pbias_mcmc_all$type[df_pbias_mcmc_all$scenario%in%c( "regimeCap"  , "regimeProd" , "regimeProdCap", "shiftCap",
  "shiftProd") ]<-"regime"
df_pbias_mcmc_all$type[df_pbias_mcmc_all$scenario%in%c("decLinearProdshiftCap" )]<-"combo"


df_pbias_mcmc_all$modtype<-"simple"
df_pbias_mcmc_all$modtype[df_pbias_mcmc_all$model%in%c("rwa" ,"hmma")]<-"tva"
df_pbias_mcmc_all$modtype[df_pbias_mcmc_all$model%in%c("rwb" ,"hmmb") ]<-"tvb"
df_pbias_mcmc_all$modtype[df_pbias_mcmc_all$model%in%c("rwab" ,"hmmab")]<-"tvab"



#time-varying alpha and Smax scenarios
tvbxp<-df_pbias_mcmc_all %>% filter( paramch%in%c("a","b")&
                            model %in%c("rwa","hmma","simple", "rwb","hmmb")&
                            parameter %in% c("alpha","smax","smsy")
                            )

tvbxp$model<-factor(tvbxp$model, levels=c("rwa","hmma",
        "simple", 
       "rwb",
        "hmmb"))

tvbxp$parameter2<-recode(tvbxp$parameter, "alpha"="parameter: alpha",
    "smax"="parameter: smax",
    "smsy"="parameter: smsy")

tvbxp$paramch2<-recode(tvbxp$paramch, "a"="scenario: time-varying alpha", "b"="scenario: time-varying Smax")
?geom_violin

pbias_tv_v<-ggplot(tvbxp) + 
geom_violin(aes(x=model,y=pbias, fill=modtype), scale="width", trim=TRUE, alpha=.7,adjust = 1.8)+
geom_boxplot(aes(x=model,y=pbias, fill=modtype),outlier.shape = NA,width=0.1)+
coord_cartesian(ylim =  c(-80, 80))+
 scale_fill_viridis_d("estimation:",option = "E") +
mytheme+
facet_grid(parameter2~paramch2)+
geom_hline(yintercept=0, linewidth=1.2) +
 ylab(" % bias")

pbias_tv_v

ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/violins/pbias_tva_tvb_allpars.png", plot=pbias_tv_v)


tvbxp$scenario<-factor(tvbxp$scenario, levels=c("decLinearProd",
                                                 "sineProd",
                                                "regimeProd",
                                                "shiftProd",
                                                "decLinearCap",
                                                "regimeCap",
                                                "shiftCap"      ))

bias_tva_v_scn<-ggplot(tvbxp) + 
geom_violin(aes(x=model,y=bias, fill=modtype),scale="width", trim=TRUE, alpha=.7,adjust = 1.8)+
geom_boxplot(aes(x=model,y=bias, fill=modtype),outlier.shape = NA,width=0.1)+
 scale_fill_viridis_d("estimation:",option = "E") +
mytheme+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
geom_hline(yintercept=0, linewidth=1.2) +
 ylab("parameter bias")+
facet_wrap(parameter~scenario, scales="free", ncol=7)
 #annotate("text", x = Inf, y = Inf, label = paste("time-varying alpha"), vjust = 1, hjust = 1.2, size=5)
#bias_tva_v_scn


init_scales_orig =bias_tva_v_scn$facet$init_scales

tvbxp_a <-tvbxp %>% filter( parameter %in% c("alpha"))
alpha_bias_rng<-as.numeric(quantile(tvbxp_a$bias, c(0.025, 0.975)))


tvbxp_smax <-tvbxp %>% filter( parameter %in% c("smax"))
smax_bias_rng<-as.numeric(quantile(tvbxp_smax$bias, c(0.025, 0.975)))


tvbxp_smsy <-tvbxp %>% filter( parameter %in% c("smsy"))
smsy_bias_rng<-as.numeric(quantile(tvbxp_smsy$bias, c(0.025, 0.975)))

bias_tva_v_ylims = list(alpha_bias_rng, alpha_bias_rng, alpha_bias_rng, alpha_bias_rng, alpha_bias_rng, alpha_bias_rng, alpha_bias_rng,
    smax_bias_rng, smax_bias_rng, smax_bias_rng, smax_bias_rng, smax_bias_rng, smax_bias_rng, smax_bias_rng,
    smsy_bias_rng, smsy_bias_rng, smsy_bias_rng, smsy_bias_rng, smsy_bias_rng, smsy_bias_rng, smsy_bias_rng)


bias_v_scn_ylim<-scale_inidividual_facet_y_axes(bias_tva_v_scn, ylims = bias_tva_v_ylims)
bias_v_scn_ylim
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/violins/bias_v_scn_ylim.png", plot=bias_v_scn_ylim)


#==================================================


bias_tva_tvb_v<-ggplot(tvbxp) + 
geom_violin(aes(x=model,y=bias, fill=modtype),scale="width", trim=TRUE,alpha=.7,adjust = 1.8)+
geom_boxplot(aes(x=model,y=bias, fill=modtype),outlier.shape = NA,width=0.1)+
 scale_fill_viridis_d("estimation:",option = "E") +
mytheme+theme(axis.text.x = element_text(angle = 45,  hjust=1))+
geom_hline(yintercept=0, linewidth=1.2) +
 ylab("parameter bias")+
#facet_grid(parameter2~paramch2)+
facet_wrap(parameter2~paramch2, scales="free", ncol=2)
 #annotate("text", x = Inf, y = Inf, label = paste("time-varying alpha"), vjust = 1, hjust = 1.2, size=5)
bias_tva_tvb_v



bias_tva_tvb_ylims = list(alpha_bias_rng, alpha_bias_rng, 
    smax_bias_rng, smax_bias_rng, 
    smsy_bias_rng, smsy_bias_rng)

bias__tvb_v_ylim<-scale_inidividual_facet_y_axes(bias_tva_tvb_v, ylims = bias_tva_tvb_ylims)
bias_tva_tvb_v_ylim


ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/violins/bias_v_ylim.png", plot=bias_tva_tvb_v_ylim)
 
#==================================================




#alpha
tvabxp_a<-df_pbias_mcmc_all %>% filter( paramch=="a"&
                            model %in%c("rwa","hmma","simple", "rwb","hmmb")&
                            parameter %in% c("alpha"))

tvabxp_a$model<-factor(tvabxp_a$model, levels=c("rwa","hmma",
        "simple", 
       "rwb",
        "hmmb"))



pl_tva_a_v<-ggplot(tvabxp_a) + 
geom_violin(aes(x=model,y=bias, fill=modtype),scale="width", trim=TRUE,alpha=.7,adjust = 1.8)+
geom_boxplot(aes(x=model,y=bias, fill=modtype),outlier.shape = NA,width=0.1)+
coord_cartesian(ylim = quantile(tvabxp_a$bias, c(0.05, 0.95)))+
 scale_fill_viridis_d("estimation:",option = "E") +
mytheme+
geom_hline(yintercept=0, linewidth=1.2) +
 ylab("log(alpha) bias")+
 ggtitle("time-varying alpha")
#annotate("text", x = Inf, y = Inf, label = paste("time-varying alpha"), vjust = 1, hjust = 1.2, size=5)
pl_tva_a_v



pl_tva_a_v_scn<-ggplot(tvabxp_a) + 
geom_violin(aes(x=model,y=bias, fill=modtype),scale="width", trim=TRUE,alpha=.7,adjust = 1.8)+
geom_boxplot(aes(x=model,y=bias, fill=modtype),outlier.shape = NA,width=0.1)+
coord_cartesian(ylim = quantile(tvabxp_a$bias, c(0.05, 0.95)))+
 scale_fill_viridis_d("estimation:",option = "E") +
mytheme+
geom_hline(yintercept=0, linewidth=1.2) +
 ylab("log(alpha) bias")+
facet_wrap(~scenario)
#annotate("text", x = Inf, y = Inf, label = paste("time-varying alpha"), vjust = 1, hjust = 1.2, size=5)
pl_tva_a_v_scn
#ggsave("figures/pl_tva_a_v_scn.png", plot=pl_tva_a_v_scn)







tvabxp_smax<-df_pbias_mcmc_all %>% filter( paramch=="a"&
                            model %in%c("rwa","hmma","simple", "rwb","hmmb")&
                            parameter %in% c("smax"))

tvabxp_smax$model<-factor(tvabxp_smax$model, levels=c("rwa","hmma",
        "simple", 
       "rwb",
        "hmmb"))




pl_tva_smax_v<-ggplot(tvabxp_smax) + 
geom_violin(aes(x=model,y=bias, fill=modtype),scale="width", trim=TRUE,alpha=.7,adjust = 1.8)+
geom_boxplot(aes(x=model,y=bias, fill=modtype),outlier.shape = NA, width=0.1)+
coord_cartesian(ylim = quantile(tvabxp_smax$bias, c(0.05, 0.95)))+
 scale_fill_viridis_d("estimation:",option = "E") +
mytheme+
geom_hline(yintercept=0, color="black",linewidth=1.2) +
 ylab("Smax bias")+
 ggtitle("time-varying alpha")
 #annotate("text", x = Inf, y = Inf, label = paste("time-varying alpha"), vjust = 1, hjust = 1.2, size=5)
pl_tva_smax_v



pl_tva_smax_v_scn<-ggplot(tvabxp_smax) + 
geom_violin(aes(x=model,y=bias, fill=modtype),scale="width", trim=TRUE,alpha=.7,adjust = 1.8)+
geom_boxplot(aes(x=model,y=bias, fill=modtype),outlier.shape = NA, width=0.1)+
coord_cartesian(ylim = quantile(tvabxp_smax$bias, c(0.05, 0.95)))+
 scale_fill_viridis_d("estimation:",option = "E") +
mytheme+
geom_hline(yintercept=0, color="black",linewidth=1.2) +
 ylab("Smax bias")+
 facet_wrap(~scenario)
 #annotate("text", x = Inf, y = Inf, label = paste("time-varying alpha"), vjust = 1, hjust = 1.2, size=5)
pl_tva_smax_v_scn




tvabxp_smsy<-df_pbias_mcmc_all %>% filter( paramch=="a"&
                            model %in%c("rwa","hmma","simple", "rwb","hmmb")&
                            parameter %in% c("smsy"))

tvabxp_smsy$model<-factor(tvabxp_smsy$model, levels=c("rwa","hmma",
        "simple", 
       "rwb",
        "hmmb"))




pl_tva_smsy_v<-ggplot(tvabxp_smsy) + 
geom_violin(aes(x=model,y=bias, fill=modtype),scale="width", trim=TRUE,alpha=.7,adjust = 1.8)+
geom_boxplot(aes(x=model,y=bias, fill=modtype),outlier.shape = NA,width=0.1)+
coord_cartesian(ylim = quantile(tvabxp_smsy$bias, c(0.05, 0.95)))+
 scale_fill_viridis_d("estimation:",option = "E") +
mytheme+
geom_hline(yintercept=0, color="black",linewidth=1.2) +
 ylab("Smsy bias")+
 ggtitle("time-varying alpha")
#annotate("text", x = Inf, y = Inf, label = paste("time-varying alpha"), vjust = 1, hjust = 1.2, size=5)
pl_tva_smsy_v





pl_tva_smsy_v_scn<-ggplot(tvabxp_smsy) + 
geom_violin(aes(x=model,y=bias, fill=modtype),scale="width", trim=TRUE,alpha=.7,adjust = 1.8)+
geom_boxplot(aes(x=model,y=bias, fill=modtype),outlier.shape = NA,width=0.1)+
coord_cartesian(ylim = quantile(tvabxp_smsy$bias, c(0.05, 0.95)))+
 scale_fill_viridis_d("estimation:",option = "E") +
mytheme+
geom_hline(yintercept=0, color="black",linewidth=1.2) +
 ylab("Smsy bias")+
  facet_wrap(~scenario)
 #annotate("text", x = Inf, y = Inf, label = paste("time-varying alpha"), vjust = 1, hjust = 1.2, size=5)
pl_tva_smsy_v_scn
#ggsave("figures/pl_tva_smsy_v_scn.png", pl_tva_smsy_v_scn)




tva_bias_allpars<-plot_grid(pl_tva_a_v,pl_tva_smax_v,pl_tva_smsy_v, 
    ncol=3,label_size = 10)


#Scenario specific para
tva_bias_scn_allpars<-plot_grid(pl_tva_a_v_scn,pl_tva_smax_v_scn,pl_tva_smsy_v_scn, 
    ncol=3,label_size = 10)


#ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/violins/bias_v_ylim.png",
$ plot=tva_bias_scn_allpars)
 

#"Smax" scenarios


tvbbxp_a<-df_pbias_mcmc_all %>% filter( paramch=="b"&
                            model %in%c("rwa","hmma","simple", "rwb","hmmb")&
                            parameter %in% c("alpha"))

tvbbxp_a$model<-factor(tvbbxp_a$model, levels=c("rwa","hmma",
        "simple", 
       "rwb",
        "hmmb"))



pl_tvb_a_v<-ggplot(tvbbxp_a) + 
geom_violin(aes(x=model,y=bias, fill=modtype),scale="width", trim=TRUE,alpha=.7,adjust = 1.8)+
geom_boxplot(aes(x=model,y=bias, fill=modtype),outlier.shape = NA,width=0.1)+
coord_cartesian(ylim =  quantile(tvbbxp_a$bias, c(0.05, 0.95)))+
 scale_fill_viridis_d("estimation:",option = "E") +
mytheme+
geom_hline(yintercept=0, color="black",alpha=.7,linewidth=1.2) +
 ylab("log(alpha) bias")+
  ggtitle("time-varying Smax")
#annotate("text", x = Inf, y = Inf, label = paste("time-varying Smax"), vjust = 1, hjust = 1.2, size=5)
pl_tvb_a_v



pl_tvb_a_v_scn<-ggplot(tvbbxp_a) + 
geom_violin(aes(x=model,y=bias, fill=modtype),scale="width", trim=TRUE,alpha=.7,adjust = 1.8)+
geom_boxplot(aes(x=model,y=bias, fill=modtype),outlier.shape = NA,width=0.1)+
coord_cartesian(ylim =  quantile(tvbbxp_a$bias, c(0.05, 0.95)))+
 scale_fill_viridis_d("estimation:",option = "E") +
mytheme+
geom_hline(yintercept=0, color="black",alpha=.7,linewidth=1.2) +
 ylab("log(alpha) bias")+
   facet_wrap(~scenario,ncol=2)
pl_tvb_a_v_scn




tvbbxp_smax<-df_pbias_mcmc_all %>% filter( paramch=="b"&
                            model %in%c("rwa","hmma","simple", "rwb","hmmb")&
                            parameter %in% c("smax"))

tvbbxp_smax$model<-factor(tvbbxp_smax$model, levels=c("rwa","hmma",
        "simple", 
       "rwb",
        "hmmb"))


pl_tvb_smax_v<-ggplot(tvbbxp_smax) + 
geom_violin(aes(x=model,y=bias, fill=modtype),scale="width", trim=TRUE,alpha=.7,adjust = 1.8)+
geom_boxplot(aes(x=model,y=bias, fill=modtype),outlier.shape = NA,width=0.1)+
coord_cartesian(ylim = quantile(tvbbxp_smax$bias, c(0.05, 0.95)))+
 scale_fill_viridis_d("estimation:",option = "E") +
mytheme+
geom_hline(yintercept=0, color="black",linewidth=1.2) +
 ylab("Smax bias")+
 ggtitle("time-varying Smax")
#annotate("text", x = Inf, y = Inf, label = paste("time-varying Smax"), vjust = 1, hjust = 1.2, size=5)
pl_tvb_smax_v


pl_tvb_smax_v_scn<-ggplot(tvbbxp_smax) + 
geom_violin(aes(x=model,y=bias, fill=modtype),scale="width", trim=TRUE,alpha=.7,adjust = 1.8)+
geom_boxplot(aes(x=model,y=bias, fill=modtype),outlier.shape = NA,width=0.1)+
coord_cartesian(ylim = quantile(tvbbxp_smax$bias, c(0.05, 0.95)))+
 scale_fill_viridis_d("estimation:",option = "E") +
mytheme+
geom_hline(yintercept=0, color="black",linewidth=1.2) +
   facet_wrap(~scenario,ncol=2)+
 ylab("Smax bias")
 pl_tvb_smax_v_scn


tvbbxp_smsy<-df_pbias_mcmc_all %>% filter( paramch=="b"&
                            model %in%c("rwa","hmma","simple", "rwb","hmmb")&
                            parameter %in% c("smsy"))

tvbbxp_smsy$model<-factor(tvbbxp_smsy$model, levels=c("rwa","hmma",
        "simple", 
       "rwb",
        "hmmb"))




pl_tvb_smsy_v<-ggplot(tvbbxp_smsy) + 
geom_violin(aes(x=model,y=bias, fill=modtype),scale="width", trim=TRUE,alpha=.7,adjust = 1.8)+
geom_boxplot(aes(x=model,y=bias, fill=modtype),outlier.shape = NA,width=0.1)+
coord_cartesian(ylim = quantile(tvbbxp_smsy$bias, c(0.05, 0.95)))+
 scale_fill_viridis_d("estimation:",option = "E") +
mytheme+
geom_hline(yintercept=0, color="black",alpha=.7,linewidth=1.2) +
 ylab("Smsy bias")+
 ggtitle("time-varying Smax")
#annotate("text", x = Inf, y = Inf, label = paste("time-varying Smax"), vjust = 1, hjust = 1.2, size=5)
pl_tvb_smsy_v





pl_tvb_smsy_v_scn<-ggplot(tvbbxp_smsy) + 
geom_violin(aes(x=model,y=bias, fill=modtype),scale="width", trim=TRUE,alpha=.7,adjust = 1.8)+
geom_boxplot(aes(x=model,y=bias, fill=modtype),outlier.shape = NA,width=0.1)+
coord_cartesian(ylim = quantile(tvbbxp_smsy$bias, c(0.05, 0.95)))+
 scale_fill_viridis_d("estimation:",option = "E") +
mytheme+
geom_hline(yintercept=0, color="black",alpha=.7,linewidth=1.2) +
 ylab("Smsy bias")+
   facet_wrap(~scenario, ncol=2)
   pl_tvb_smsy_v_scn



plot_grid(pl_tva_a_v,pl_tva_smax_v,pl_tva_smsy_v,pl_tvb_a_v,pl_tvb_smax_v,pl_tvb_smsy_v, 
    ncol=3,label_size = 8)




pl_tva_tvb_scn_allpars<-plot_grid(pl_tva_a_v_scn,pl_tva_smax_v_scn,pl_tva_smsy_v_scn,
          pl_tvb_a_v_scn,pl_tvb_smax_v_scn,pl_tvb_smsy_v_scn, 
    ncol=3,label_size = 8)

ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/violins/bias_tva_tvb_scn_allpars.png",
 plot=pl_tva_tvb_scn_allpars)
 
#?plot_grid
#====================================================================================

# scenarios with both parameters are time-varying 

tvabbxp<-df_pbias_mcmc_all %>% filter( paramch%in%c("both")&
                            model %in%c("rwab","hmmab", "simple","rwa","hmma", "rwb","hmmb")&
                            parameter %in% c("alpha","smax","smsy")
                            )

tvabbxp$model<-factor(tvabbxp$model, levels=c("rwab","hmmab",
        "simple","rwa","hmma", "rwb","hmmb"))

tvabbxp$parameter2<-recode(tvabbxp$parameter, "alpha"="parameter: alpha",
    "smax"="parameter: smax",
    "smsy"="parameter: smsy")

#tvabbxp$paramch2<-recode(tvabbxp$paramch, "both"="scenario: time-varying alpha", "b"="scenario: time-varying Smax")
unique(tvabbxp$paramch)

pbias_tvab_v<-ggplot(tvabbxp) + 
geom_violin(aes(x=model,y=pbias, fill=modtype),scale="width", trim=TRUE,alpha=.7,adjust = 1.8)+
geom_boxplot(aes(x=model,y=pbias, fill=modtype),outlier.shape = NA,width=0.1)+
coord_cartesian(ylim =  c(-80, 80))+
 scale_fill_viridis_d("estimation:",option = "E") +
mytheme+
facet_grid(parameter2~paramch)+
geom_hline(yintercept=0, linewidth=1.2) +
 ylab(" % bias")

pbias_tvab_v

ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/violins/pbias_tvab_allpars.png",
 plot=pbias_tvab_v)
 





tvabbxp_a<-df_pbias_mcmc_all %>% filter( paramch=="both"&
                            model %in%c("rwab","hmmab","simple",
                                "rwa","hmma", "rwb","hmmb")&
                            parameter %in% c("alpha"))

tvabbxp_a$model<-factor(tvabbxp_a$model, levels=c("rwab","hmmab",
        "simple","rwa","hmma", "rwb","hmmb"))


pl_tvab_a_v<-ggplot(tvabbxp_a) + 
geom_violin(aes(x=model,y=bias, fill=modtype),scale="width", trim=TRUE,alpha=.7,adjust = 1.8)+
geom_boxplot(aes(x=model,y=bias, fill=modtype),outlier.shape = NA,width=0.1)+
coord_cartesian(ylim = quantile(tvabbxp_a$bias, c(0.05, 0.95)))+
 scale_fill_viridis_d("estimation:",option = "E") +
mytheme+
geom_hline(yintercept=0, color="black",linewidth=1.2) +
 ylab("log(alpha) bias")+
 ggtitle("time-varying alpha and Smax")
 #annotate("text", x = Inf, y = Inf, label = paste("time-varying alpha and Smax "), vjust = 1, hjust = 1.2, size=6)
 pl_tvab_a_v


pl_tvab_a_v_scn<-ggplot(tvabbxp_a) + 
geom_violin(aes(x=model,y=bias, fill=modtype),scale="width", trim=TRUE,alpha=.7,adjust = 1.8)+
geom_boxplot(aes(x=model,y=bias, fill=modtype),outlier.shape = NA,width=0.1)+
coord_cartesian(ylim = quantile(tvabbxp_a$bias, c(0.05, 0.95)))+
 scale_fill_viridis_d("estimation:",option = "E") +
mytheme+
   facet_wrap(~scenario)+
geom_hline(yintercept=0, color="black",linewidth=1.2) +
 ylab("log(alpha) bias")
pl_tvab_a_v_scn



tvabbxp_smax<-df_pbias_mcmc_all %>% filter( paramch=="both"&
                            model %in%c("rwab","hmmab","simple",
                                "rwa","hmma", "rwb","hmmb")&
                            parameter %in% c("smax"))

tvabbxp_smax$model<-factor(tvabbxp_smax$model, levels=c("rwab","hmmab",
        "simple","rwa","hmma", "rwb","hmmb"))


pl_tvab_smax_v<-ggplot(tvabbxp_smax) + 
geom_violin(aes(x=model,y=bias, fill=modtype),scale="width", trim=TRUE,alpha=.7,adjust = 1.8)+
geom_boxplot(aes(x=model,y=bias, fill=modtype),outlier.shape = NA,width=0.1)+
coord_cartesian(ylim = quantile(tvabbxp_smax$bias, c(0.05, 0.95)))+
 scale_fill_viridis_d("estimation:",option = "E") +
mytheme+
geom_hline(yintercept=0, color="black",linewidth=1.2) +
 ylab("Smax bias")+
  ggtitle("time-varying alpha and Smax")
 #annotate("text", x = Inf, y = Inf, label = paste("time-varying alpha and Smax"), vjust = 1, hjust = 1.2, size=6)
 pl_tvab_smax_v


pl_tvab_smax_v_scn<-ggplot(tvabbxp_smax) + 
geom_violin(aes(x=model,y=bias, fill=modtype),scale="width", trim=TRUE,alpha=.7,adjust = 1.8)+
geom_boxplot(aes(x=model,y=bias, fill=modtype),outlier.shape = NA,width=0.1)+
coord_cartesian(ylim = quantile(tvabbxp_smax$bias, c(0.05, 0.95)))+
 scale_fill_viridis_d("estimation:",option = "E") +
mytheme+
   facet_wrap(~scenario)+
geom_hline(yintercept=0, color="black",linewidth=1.2) +
 ylab("Smax bias")
 #annotate("text", x = Inf, y = Inf, label = paste("time-varying alpha and Smax"), vjust = 1, hjust = 1.2, size=6)
 pl_tvab_smax_v_scn



tvabbxp_smsy<-df_pbias_mcmc_all %>% filter( paramch=="both"&
                            model %in%c("rwab","hmmab","simple",
                                "rwa","hmma", "rwb","hmmb")&
                            parameter %in% c("smsy"))

tvabbxp_smsy$model<-factor(tvabbxp_smsy$model, levels=c("rwab","hmmab",
        "simple","rwa","hmma", "rwb","hmmb"))



 pl_tvab_smsy_v<-ggplot(tvabbxp_smsy) + 
geom_violin(aes(x=model,y=bias, fill=modtype),scale="width", trim=TRUE,alpha=.7,adjust = 1.8)+
geom_boxplot(aes(x=model,y=bias, fill=modtype),outlier.shape = NA,width=0.1)+
coord_cartesian(ylim = quantile(tvabbxp_smsy$bias, c(0.05, 0.95)))+
 scale_fill_viridis_d("estimation:",option = "E") +
mytheme+
geom_hline(yintercept=0, color="black",linewidth=1.2) +
 ylab("Smsy bias")+
ggtitle("time-varying alpha and Smax")
#annotate("text", x = Inf, y = Inf, label = paste("time-varying alpha and Smax"), vjust = 1, hjust = 1.2, size=6)
 pl_tvab_smsy_v


 pl_tvab_smsy_v_scn<-ggplot(tvabbxp_smsy) + 
geom_violin(aes(x=model,y=bias, fill=modtype),scale="width", trim=TRUE,alpha=.7,adjust = 1.8)+
geom_boxplot(aes(x=model,y=bias, fill=modtype),outlier.shape = NA,width=0.1)+
coord_cartesian(ylim = quantile(tvabbxp_smsy$bias, c(0.05, 0.95)))+
 scale_fill_viridis_d("estimation:",option = "E") +
mytheme+
   facet_wrap(~scenario)+
geom_hline(yintercept=0, color="black",linewidth=1.2) +
 ylab("Smsy bias")
  pl_tvab_smsy_v_scn






pl_tvab_allpars<-plot_grid(pl_tvab_a_v,pl_tvab_smax_v, pl_tvab_smsy_v,ncol=3)

ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/violins/bias_tvab_allpars.png",
 plot=pl_tvab_allpars)




pl_tvab_scn_allpars<-plot_grid(pl_tvab_a_v_scn,pl_tvab_smax_v_scn, pl_tvab_smsy_v_scn,ncol=1)


ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/violins/bias_tvab_scn_allpars.png",
 plot=pl_tvab_scn_allpars)


#===========================================
# stationary scenarios- no time-varying parameters  
unique(df_pbias_mcmc_all$paramch)


notvbxp<-df_pbias_mcmc_all %>% filter( paramch=="stationary"&
                            model %in%c("rwab","hmmab","simple",
                                "rwa","hmma", "rwb","hmmb")&
                            parameter %in%  c("alpha","smax","smsy"))

notvbxp$model<-factor(notvbxp$model, levels=c("rwab","hmmab",
        "simple","rwa","hmma", "rwb","hmmb"))


pl_notv_v<-ggplot(notvbxp) + 
geom_violin(aes(x=model,y=pbias, fill=modtype),scale="width", trim=TRUE,alpha=.7,adjust = 1.8)+
geom_boxplot(aes(x=model,y=pbias, fill=modtype),outlier.shape = NA,width=0.1)+
coord_cartesian(ylim =c(-60,60))+
 scale_fill_viridis_d("estimation:",option = "E") +
mytheme+
facet_wrap(~parameter, ncol=3)+
geom_hline(yintercept=0, color="black",linewidth=1.2) +
 ylab("% bias")
 pl_notv_v


ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/violins/pbias_stn_allpars.png",
 plot= pl_notv_v)



pl_notv_v_scn<-ggplot(notvbxp) + 
geom_violin(aes(x=model,y=pbias, fill=modtype),scale="width", trim=TRUE,alpha=.7,adjust = 1.8)+
geom_boxplot(aes(x=model,y=pbias, fill=modtype),outlier.shape = NA,width=0.1)+
coord_cartesian(ylim =c(-60,60))+
 scale_fill_viridis_d("estimation:",option = "E") +
mytheme+
facet_wrap(scenario~parameter, ncol=3)+
geom_hline(yintercept=0, color="black",linewidth=1.2) +
 ylab("% bias")
 pl_notv_v_scn


ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/violins/pbias_stn_scn_allpars.png",
 plot=  pl_notv_v_scn)


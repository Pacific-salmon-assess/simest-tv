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
simPar <- read.csv("data/sensitivity/SimPars.csv")

## Store relevant object names to help run simulation 
scenNames <- unique(simPar$scenario)

resa1<-readRDS(file = "outs/simest/sensitivity/res_a1.rds")
resa2<-readRDS(file = "outs/simest/sensitivity/res_a2.rds")


restmb<-rbind(resa1,resa2)

head(restmb)

#res<-readRDS(file = "outs/simest/generic/res.rds")
#restmb<-restmb[restmb$convergence==0,]

#

resstana1<-readRDS(file = "outs/simest/sensitivity/resstan_a1.rds")
resstana2<-readRDS(file = "outs/simest/sensitivity/resstan_a2.rds")
resstan<-rbind(resstana1,resstana2)

#resstan<-readRDS(file = "outs/simest/generic/resstan.rds")

res<-rbind(restmb,resstan)

res$parameter[res$parameter=="Smax"]<-"smax"

resparam<-res[res$parameter%in%c("alpha","smax","sigma","smsy","sgen","umsy"),]

head(resparam)

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


allconv<-inner_join(convstatMLE[,-3], convstatMCMC[,-3])

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

df_pbias<-df_pbias[df_pbias$parameter!="sigma"&df_pbias$parameter!="sgen"&!is.na(df_pbias$bias),]


df_pbias_mcmc<-df_pbias[df_pbias$method=="MCMC",]




df_pbias_mcmc$paramch<-"a" 
 
df_pbias_mcmc$type<-"trend"
df_pbias_mcmc$type[df_pbias_mcmc$scenario%in%c("regimeProd1","regimeProd10","regimeProd2","regimeProd5","regimeProd7")]<-"shift"

#time-varying alpha and Smax scenarios
tvbxp<-df_pbias_mcmc %>% filter( 
                            parameter %in% c("alpha","smax","smsy")
                            )

tvbxp$model<-factor(tvbxp$model, levels=c("simple", 
       "autocorr",
       "rwa","hmma",  
       "rwb","hmmb",  
       "rwab","hmmab"))

unique(tvbxp$scenario)

tvbxp$scenario<-factor(tvbxp$scenario, levels=c("regimeProd1",
                                                "regimeProd2",
                                                "regimeProd5",
                                                "regimeProd7",  
                                                "regimeProd10",    
                                                "trendLinearProd1",
                                                "trendLinearProd2",
                                                "trendLinearProd5",
                                                "trendLinearProd7",
                                                "trendLinearProd10" ))


tvbxp$modtype<-case_match(
  tvbxp$model,
  c("simple", "autocorr") ~ "simple",
  c("rwa", "hmma") ~ "tva",
   c("rwb", "hmmb") ~ "tvb", 
    c("rwab", "hmmab") ~ "both")

tvbxp$model<-factor(tvbxp$model, levels=c("simple", 
       "autocorr",
       "rwa","hmma",  
       "rwb","hmmb",  
       "rwab","hmmab"))




#plot bias (instead of pbias)
tvbxp_alpha<-tvbxp[tvbxp$parameter=="alpha",]


bias_alpha_scn<-ggplot(tvbxp_alpha) + 
geom_violin(aes(x=model,y=bias, fill=modtype),scale="width", trim=TRUE, alpha=.7,position = position_dodge(0.9))+
geom_boxplot(aes(x=model,y=bias, fill=modtype),outlier.shape = NA,width=0.1,position = position_dodge(0.9))+
coord_cartesian(ylim =  c(-1.7, 1.5))+
 scale_fill_viridis_d("estimation:",option = "E") +
mytheme+
geom_hline(yintercept=0, linewidth=1.2) +
 ylab("log(alpha) bias")+
 theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
facet_wrap(.~ scenario, scales="free", ncol=5)
 #annotate("text", x = Inf, y = Inf, label = paste("time-varying alpha"), vjust = 1, hjust = 1.2, size=5)
bias_alpha_scn

ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/violins/sens_a/bias_alpha_scn_sensa.png", plot=bias_alpha_scn)
 



tvbxp_smax<-tvbxp[tvbxp$parameter=="smax",]


bias_smax_scn<-ggplot(tvbxp_smax) + 
geom_violin(aes(x=model,y=bias, fill=modtype),scale="width", trim=TRUE, alpha=.7,position = position_dodge(0.9))+
geom_boxplot(aes(x=model,y=bias, fill=modtype),outlier.shape = NA,width=0.1,position = position_dodge(0.9))+
coord_cartesian(ylim =  c(-150000, 250000))+
 scale_fill_viridis_d("estimation:",option = "E") +
mytheme+
geom_hline(yintercept=0, linewidth=1.2) +
 ylab("smax bias")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
facet_wrap(.~ scenario, scales="free", ncol=5)
bias_smax_scn

ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/violins/sens_a/bias_smax_scn_sensa.png", plot=bias_smax_scn)
 







tvbxp_smsy<-tvbxp[tvbxp$parameter=="smsy",]


bias_smsy_scn<-ggplot(tvbxp_smax) + 
geom_violin(aes(x=model,y=bias, fill=modtype),scale="width", trim=TRUE, alpha=.7,position = position_dodge(0.9))+
geom_boxplot(aes(x=model,y=bias, fill=modtype),outlier.shape = NA,width=0.1,position = position_dodge(0.9))+
coord_cartesian(ylim =  c(-150000, 250000))+
 scale_fill_viridis_d("estimation:",option = "E") +
mytheme+
geom_hline(yintercept=0, linewidth=1.2) +
 ylab("smsy bias")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
facet_wrap(.~ scenario, scales="free", ncol=5)
bias_smsy_scn


ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/violins/sens_a/bias_smsy_scn_sensa.png", plot=bias_smsy_scn)
 




init_scales_orig =bias_tva_v_scn$facet$init_scales

tvbxp_a <-tvabxp %>% filter( parameter %in% c("alpha"))
alpha_bias_rng<-as.numeric(quantile(tvbxp_a$bias, c(0.025, 0.975)))


tvbxp_smax <-tvabxp %>% filter( parameter %in% c("smax"))
smax_bias_rng<-as.numeric(quantile(tvbxp_smax$bias, c(0.025, 0.975)))


tvbxp_smsy <-tvabxp %>% filter( parameter %in% c("smsy"))
smsy_bias_rng<-as.numeric(quantile(tvbxp_smsy$bias, c(0.025, 0.975)))

bias_tva_v_ylims = list(alpha_bias_rng, alpha_bias_rng, alpha_bias_rng, alpha_bias_rng, alpha_bias_rng, alpha_bias_rng, alpha_bias_rng,
    smax_bias_rng, smax_bias_rng, smax_bias_rng, smax_bias_rng, smax_bias_rng, smax_bias_rng, smax_bias_rng,
    smsy_bias_rng, smsy_bias_rng, smsy_bias_rng, smsy_bias_rng, smsy_bias_rng, smsy_bias_rng, smsy_bias_rng)


bias_v_scn_ylim<-scale_inidividual_facet_y_axes(bias_tva_v_scn, ylims = bias_tva_v_ylims)
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/bias_v_scn_ylim.png", plot=bias_v_scn_ylim)


ggplot(tvabxp) + 
geom_violin(aes(x=model,y=pbias, fill=modtype), scale="width", trim=TRUE, alpha=.7)+
geom_boxplot(aes(x=model,y=pbias, fill=modtype),outlier.shape = NA,width=0.1)+
coord_cartesian(ylim =  c(-80, 80))+
 scale_fill_viridis_d("estimation:",option = "E") +
mytheme+
facet_grid(parameter2~paramch2)+
geom_hline(yintercept=0, linewidth=1.2) +
 ylab(" % bias")

pbias_tv_v


#alpha
tvabxp_a<-df_pbias_mcmc_all %>% filter( paramch=="a"&
                            model %in%c("rwa","hmma","simple", "rwb","hmmb")&
                            parameter %in% c("alpha"))

tvabxp_a$model<-factor(tvabxp_a$model, levels=c("rwa","hmma",
        "simple", 
       "rwb",
        "hmmb"))

head(tvabxp_a)


pl_tva_a_v<-ggplot(tvabxp_a) + 
geom_violin(aes(x=model,y=bias, fill=modtype),alpha=.7)+
geom_boxplot(aes(x=model,y=bias, fill=modtype),outlier.shape = NA,width=0.1)+
coord_cartesian(ylim = quantile(tvabxp_a$bias, c(0.05, 0.95)))+
 scale_fill_viridis_d("estimation:",option = "E") +
mytheme+
geom_hline(yintercept=0, linewidth=1.2) +
 ylab("log(alpha) bias")+
 annotate("text", x = Inf, y = Inf, label = paste("time-varying alpha"), vjust = 1, hjust = 1.2, size=5)
pl_tva_a_v

head(tvabxp_a)

pl_tva_a_v_scn<-ggplot(tvabxp_a) + 
geom_violin(aes(x=model,y=bias, fill=modtype),alpha=.7)+
geom_boxplot(aes(x=model,y=bias, fill=modtype),outlier.shape = NA,width=0.1)+
coord_cartesian(ylim = quantile(tvabxp_a$bias, c(0.05, 0.95)))+
 scale_fill_viridis_d("estimation:",option = "E") +
mytheme+
geom_hline(yintercept=0, linewidth=1.2) +
 ylab("log(alpha) bias")+
facet_wrap(~scenario)+
 annotate("text", x = Inf, y = Inf, label = paste("time-varying alpha"), vjust = 1, hjust = 1.2, size=5)
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
geom_violin(aes(x=model,y=bias, fill=modtype))+
geom_boxplot(aes(x=model,y=bias, fill=modtype),outlier.shape = NA,width=0.1, color="grey")+
coord_cartesian(ylim = quantile(tvabxp_smax$bias, c(0.05, 0.95)))+
 scale_fill_viridis_d("estimation:",option = "E") +
mytheme+
geom_hline(yintercept=0, color="black",linewidth=1.2) +
 ylab("Smax bias")+
 annotate("text", x = Inf, y = Inf, label = paste("time-varying alpha"), vjust = 1, hjust = 1.2, size=5)
pl_tva_smax_v



pl_tva_smax_v_scn<-ggplot(tvabxp_smax) + 
geom_violin(aes(x=model,y=bias, fill=modtype))+
geom_boxplot(aes(x=model,y=bias, fill=modtype),outlier.shape = NA,width=0.1, color="grey")+
coord_cartesian(ylim = quantile(tvabxp_smax$bias, c(0.05, 0.95)))+
 scale_fill_viridis_d("estimation:",option = "E") +
mytheme+
geom_hline(yintercept=0, color="black",linewidth=1.2) +
 ylab("Smax bias")+
 facet_wrap(~scenario)+
 annotate("text", x = Inf, y = Inf, label = paste("time-varying alpha"), vjust = 1, hjust = 1.2, size=5)
pl_tva_smax_v_scn




tvabxp_smsy<-df_pbias_mcmc_all %>% filter( paramch=="a"&
                            model %in%c("rwa","hmma","simple", "rwb","hmmb")&
                            parameter %in% c("smsy"))

tvabxp_smsy$model<-factor(tvabxp_smsy$model, levels=c("rwa","hmma",
        "simple", 
       "rwb",
        "hmmb"))




pl_tva_smsy_v<-ggplot(tvabxp_smsy) + 
geom_violin(aes(x=model,y=bias, fill=modtype),alpha=.7)+
geom_boxplot(aes(x=model,y=bias, fill=modtype),outlier.shape = NA,width=0.1)+
coord_cartesian(ylim = quantile(tvabxp_smsy$bias, c(0.05, 0.95)))+
 scale_fill_viridis_d("estimation:",option = "E") +
mytheme+
geom_hline(yintercept=0, color="black",linewidth=1.2) +
 ylab("Smsy bias")+
 annotate("text", x = Inf, y = Inf, label = paste("time-varying alpha"), vjust = 1, hjust = 1.2, size=5)
pl_tva_smsy_v





pl_tva_smsy_v_scn<-ggplot(tvabxp_smsy) + 
geom_violin(aes(x=model,y=bias, fill=modtype))+
geom_boxplot(aes(x=model,y=bias, fill=modtype),outlier.shape = NA,width=0.1, color="grey")+
coord_cartesian(ylim = quantile(tvabxp_smsy$bias, c(0.05, 0.95)))+
 scale_fill_viridis_d("estimation:",option = "E") +
mytheme+
geom_hline(yintercept=0, color="black",linewidth=1.2) +
 ylab("Smsy bias")+
  facet_wrap(~scenario)+
 annotate("text", x = Inf, y = Inf, label = paste("time-varying alpha"), vjust = 1, hjust = 1.2, size=5)
pl_tva_smsy_v_scn
ggsave("figures/pl_tva_smsy_v_scn.png", pl_tva_smsy_v_scn)







#Scenario specific para
plot_grid(pl_tva_a_v_scn,pl_tva_smax_v_scn,pl_tva_smsy_v_scn, 
    ncol=3,label_size = 10)



#"Smax" scenarios
unique(df_pbias_mcmc_all$paramch)

tvbbxp_a<-df_pbias_mcmc_all %>% filter( paramch=="b"&
                            model %in%c("rwa","hmma","simple", "rwb","hmmb")&
                            parameter %in% c("alpha"))

tvbbxp_a$model<-factor(tvbbxp_a$model, levels=c("rwa","hmma",
        "simple", 
       "rwb",
        "hmmb"))

pl_tvb_a<-ggplot(tvbbxp_a) + 
geom_boxplot(aes(x=model,y=bias, fill=modtype),outlier.shape = NA)+
coord_cartesian(ylim =  quantile(tvbbxp_a$bias, c(0.05, 0.95)))+
 scale_fill_viridis_d("estimation:",option = "E") +
mytheme+
geom_hline(yintercept=0, color="black",alpha=.7,linewidth=1.2) +
 ylab("log(alpha) bias")+
 annotate("text", x = Inf, y = Inf, label = paste("time-varying Smax"), vjust = 1, hjust = 1.2, size=8)
pl_tvb_a



pl_tvb_a_v<-ggplot(tvbbxp_a) + 
geom_violin(aes(x=model,y=bias, fill=modtype),alpha=.7)+
geom_boxplot(aes(x=model,y=bias, fill=modtype),outlier.shape = NA,width=0.1)+
coord_cartesian(ylim =  quantile(tvbbxp_a$bias, c(0.05, 0.95)))+
 scale_fill_viridis_d("estimation:",option = "E") +
mytheme+
geom_hline(yintercept=0, color="black",alpha=.7,linewidth=1.2) +
 ylab("log(alpha) bias")+
 annotate("text", x = Inf, y = Inf, label = paste("time-varying Smax"), vjust = 1, hjust = 1.2, size=5)
pl_tvb_a_v



pbias_tvb_a_v<-ggplot(tvbbxp_a) + 
geom_violin(aes(x=model,y=pbias, fill=modtype),alpha=.7)+
geom_boxplot(aes(x=model,y=pbias, fill=modtype),outlier.shape = NA,width=0.1)+
coord_cartesian(ylim =  c(-100,100))+
 scale_fill_viridis_d("estimation:",option = "E") +
mytheme+
geom_hline(yintercept=0, color="black",alpha=.7,linewidth=1.2) +
 ylab("log(alpha) % bias")+
 annotate("text", x = Inf, y = Inf, label = paste("time-varying Smax"), vjust = 1, hjust = 1.2, size=5)
pbias_tvb_a_v


pl_tvb_a_v_scn<-ggplot(tvbbxp_a) + 
geom_violin(aes(x=model,y=bias, fill=modtype))+
geom_boxplot(aes(x=model,y=bias, fill=modtype),outlier.shape = NA,width=0.1, color="grey")+
coord_cartesian(ylim =  quantile(tvbbxp_a$bias, c(0.05, 0.95)))+
 scale_fill_viridis_d("estimation:",option = "E") +
mytheme+
geom_hline(yintercept=0, color="black",alpha=.7,linewidth=1.2) +
 ylab("log(alpha) bias")+
   facet_wrap(~scenario)+
 annotate("text", x = Inf, y = Inf, label = paste("time-varying Smax"), vjust = 1, hjust = 1.2, size=5)
pl_tvb_a_v_scn
ggsave("figures/pl_tvb_a_v_scn.png", pl_tvb_a_v_scn)



tvbbxp_smax<-df_pbias_mcmc_all %>% filter( paramch=="b"&
                            model %in%c("rwa","hmma","simple", "rwb","hmmb")&
                            parameter %in% c("smax"))

tvbbxp_smax$model<-factor(tvbbxp_smax$model, levels=c("rwa","hmma",
        "simple", 
       "rwb",
        "hmmb"))


pl_tvb_smax<-ggplot(tvbbxp_smax) + 
geom_boxplot(aes(x=model,y=bias, fill=modtype),outlier.shape = NA)+
coord_cartesian(ylim = quantile(tvbbxp_smax$bias, c(0.05, 0.95)))+
# facet_grid(parameter~., scales="free_y")+
 scale_fill_viridis_d("estimation:",option = "E") +
mytheme+
geom_hline(yintercept=0, color="black",alpha=.7,linewidth=1.2) +
 ylab("Smax bias")+
 annotate("text", x = Inf, y = Inf, label = paste("time-varying Smax"), vjust = 1, hjust = 1.2, size=8)
#ggtitle("Smax varying scenarios")
pl_tvb_smax

pl_tvb_smax_v<-ggplot(tvbbxp_smax) + 
geom_violin(aes(x=model,y=bias, fill=modtype))+
geom_boxplot(aes(x=model,y=bias, fill=modtype),outlier.shape = NA,width=0.1, color="grey")+
coord_cartesian(ylim = quantile(tvbbxp_smax$bias, c(0.05, 0.95)))+
 scale_fill_viridis_d("estimation:",option = "E") +
mytheme+
geom_hline(yintercept=0, color="black",linewidth=1.2) +
 ylab("Smax bias")+
 annotate("text", x = Inf, y = Inf, label = paste("time-varying Smax"), vjust = 1, hjust = 1.2, size=5)
pl_tvb_smax_v


pbias_tvb_smax_v<-ggplot(tvbbxp_smax) + 
geom_violin(aes(x=model,y=pbias, fill=modtype),alpha=.7)+
geom_boxplot(aes(x=model,y=pbias, fill=modtype),outlier.shape = NA,width=0.1)+
coord_cartesian(ylim =  c(-100, 100))+
 scale_fill_viridis_d("estimation:",option = "E") +
mytheme+
geom_hline(yintercept=0, color="black",linewidth=1.2) +
 ylab("Smax % bias")+
 annotate("text", x = Inf, y = Inf, label = paste("time-varying Smax"), vjust = 1, hjust = 1.2, size=5)
pbias_tvb_smax_v

pl_tvb_smax_v_scn<-ggplot(tvbbxp_smax) + 
geom_violin(aes(x=model,y=bias, fill=modtype))+
geom_boxplot(aes(x=model,y=bias, fill=modtype),outlier.shape = NA,width=0.1, color="grey")+
coord_cartesian(ylim = quantile(tvbbxp_smax$bias, c(0.05, 0.95)))+
 scale_fill_viridis_d("estimation:",option = "E") +
mytheme+
geom_hline(yintercept=0, color="black",linewidth=1.2) +
   facet_wrap(~scenario)+
 ylab("Smax bias")+
 annotate("text", x = Inf, y = Inf, label = paste("time-varying Smax"), vjust = 1, hjust = 1.2, size=5)
pl_tvb_smax_v_scn
ggsave("figures/pl_tvb_smax_v_scn.png", pl_tvb_a_v_scn)


tvbbxp_smsy<-df_pbias_mcmc_all %>% filter( paramch=="b"&
                            model %in%c("rwa","hmma","simple", "rwb","hmmb")&
                            parameter %in% c("smsy"))

tvbbxp_smsy$model<-factor(tvbbxp_smsy$model, levels=c("rwa","hmma",
        "simple", 
       "rwb",
        "hmmb"))


pl_tvb_smsy<-ggplot(tvbbxp_smsy) + 
geom_boxplot(aes(x=model,y=bias, fill=modtype),outlier.shape = NA)+
coord_cartesian(ylim = quantile(tvbbxp_smsy$bias, c(0.05, 0.95)))+
 scale_fill_viridis_d("estimation:",option = "E") +
mytheme+
geom_hline(yintercept=0, color="black",alpha=.7,linewidth=1.2) +
 ylab("Smsy bias")+
 annotate("text", x = Inf, y = Inf, label = paste("time-varying Smax"), vjust = 1, hjust = 1.2, size=8)
pl_tvb_smsy




pl_tvb_smsy_v<-ggplot(tvbbxp_smsy) + 
geom_violin(aes(x=model,y=bias, fill=modtype))+
geom_boxplot(aes(x=model,y=bias, fill=modtype),outlier.shape = NA,width=0.1, color="grey")+
coord_cartesian(ylim = quantile(tvbbxp_smsy$bias, c(0.05, 0.95)))+
 scale_fill_viridis_d("estimation:",option = "E") +
mytheme+
geom_hline(yintercept=0, color="black",alpha=.7,linewidth=1.2) +
 ylab("Smsy bias")+
 annotate("text", x = Inf, y = Inf, label = paste("time-varying Smax"), vjust = 1, hjust = 1.2, size=5)
pl_tvb_smsy_v



pbias_tvb_smsy_v<-ggplot(tvbbxp_smsy) + 
geom_violin(aes(x=model,y=pbias, fill=modtype),alpha=.7)+
geom_boxplot(aes(x=model,y=pbias, fill=modtype),outlier.shape = NA,width=0.1)+
coord_cartesian(ylim = c(-100, 100))+
 scale_fill_viridis_d("estimation:",option = "E") +
mytheme+
geom_hline(yintercept=0, color="black",alpha=.7,linewidth=1.2) +
 ylab("Smsy % bias")+
 annotate("text", x = Inf, y = Inf, label = paste("time-varying Smax"), vjust = 1, hjust = 1.2, size=5)
pbias_tvb_smsy_v



pl_tvb_smsy_v_scn<-ggplot(tvbbxp_smsy) + 
geom_violin(aes(x=model,y=bias, fill=modtype))+
geom_boxplot(aes(x=model,y=bias, fill=modtype),outlier.shape = NA,width=0.1, color="grey")+
coord_cartesian(ylim = quantile(tvbbxp_smsy$bias, c(0.05, 0.95)))+
 scale_fill_viridis_d("estimation:",option = "E") +
mytheme+
geom_hline(yintercept=0, color="black",alpha=.7,linewidth=1.2) +
 ylab("Smsy bias")+
   facet_wrap(~scenario)+
 annotate("text", x = Inf, y = Inf, label = paste("time-varying Smax"), vjust = 1, hjust = 1.2, size=5)
pl_tvb_smsy_v_scn
ggsave("figures/pl_tvb_smsy_v_scn.png", pl_tvb_smsy_v_scn)


plot_grid(pl_tva_a,pl_tva_smax,pl_tva_smsy,pl_tvb_a,pl_tvb_smax,pl_tvb_smsy, 
    ncol=3,label_size = 10)


plot_grid(pl_tva_a_v,pl_tva_smax_v,pl_tva_smsy_v,pl_tvb_a_v,pl_tvb_smax_v,pl_tvb_smsy_v, 
    ncol=3,label_size = 8)




plot_grid(pbias_tva_a_v,pbias_tva_smax_v,pbias_tva_smsy_v,pbias_tvb_a_v,pbias_tvb_smax_v,pbias_tvb_smsy_v, 
    ncol=3,label_size = 8)

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
geom_violin(aes(x=model,y=pbias, fill=modtype), scale="width", trim=TRUE, alpha=.7)+
geom_boxplot(aes(x=model,y=pbias, fill=modtype),outlier.shape = NA,width=0.1)+
coord_cartesian(ylim =  c(-80, 80))+
 scale_fill_viridis_d("estimation:",option = "E") +
mytheme+
facet_grid(parameter2~paramch)+
geom_hline(yintercept=0, linewidth=1.2) +
 ylab(" % bias")

pbias_tvab_v





tvabbxp_a<-df_pbias_mcmc_all %>% filter( paramch=="both"&
                            model %in%c("rwab","hmmab","simple",
                                "rwa","hmma", "rwb","hmmb")&
                            parameter %in% c("alpha"))

tvabbxp_a$model<-factor(tvabbxp_a$model, levels=c("rwab","hmmab",
        "simple","rwa","hmma", "rwb","hmmb"))

pl_tvab_a<-ggplot(tvabbxp_a) + 
geom_boxplot(aes(x=model,y=bias, fill=modtype),outlier.shape = NA)+
coord_cartesian(ylim = quantile(tvabbxp_a$bias, c(0.05, 0.95)))+
#facet_grid(rows=vars(parameter), scales="free_y")+
 scale_fill_viridis_d("estimation:",option = "E") +
mytheme+
geom_hline(yintercept=0, color="black",alpha=.7,linewidth=1.2) +
 ylab("log(alpha) bias")+
 annotate("text", x = Inf, y = Inf, label = paste("time-varying alpha and Smax "), vjust = 1, hjust = 1.2, size=6)
 pl_tvab_a

pl_tvab_a_v<-ggplot(tvabbxp_a) + 
geom_violin(aes(x=model,y=bias, fill=modtype),alpha=0.7)+
geom_boxplot(aes(x=model,y=bias, fill=modtype),outlier.shape = NA,width=0.1)+
coord_cartesian(ylim = quantile(tvabbxp_a$bias, c(0.05, 0.95)))+
 scale_fill_viridis_d("estimation:",option = "E") +
mytheme+
geom_hline(yintercept=0, color="black",linewidth=1.2) +
 ylab("log(alpha) bias")+
 annotate("text", x = Inf, y = Inf, label = paste("time-varying alpha and Smax "), vjust = 1, hjust = 1.2, size=6)
 pl_tvab_a_v


pl_tvab_a_v_scn<-ggplot(tvabbxp_a) + 
geom_violin(aes(x=model,y=bias, fill=modtype),alpha=0.7)+
geom_boxplot(aes(x=model,y=bias, fill=modtype),outlier.shape = NA,width=0.1)+
coord_cartesian(ylim = quantile(tvabbxp_a$bias, c(0.05, 0.95)))+
 scale_fill_viridis_d("estimation:",option = "E") +
mytheme+
   facet_wrap(~scenario)+
geom_hline(yintercept=0, color="black",linewidth=1.2) +
 ylab("log(alpha) bias")+
 annotate("text", x = Inf, y = Inf, label = paste("time-varying alpha and Smax "), vjust = 1, hjust = 1.2, size=6)
pl_tvab_a_v_scn


pbias_tvab_a_v<-ggplot(tvabbxp_a) + 
geom_violin(aes(x=model,y=pbias, fill=modtype),alpha=0.7)+
geom_boxplot(aes(x=model,y=pbias, fill=modtype),outlier.shape = NA,width=0.1)+
coord_cartesian(ylim =  c(-100, 100))+
 scale_fill_viridis_d("estimation:",option = "E") +
mytheme+
geom_hline(yintercept=0, color="black",linewidth=1.2) +
 ylab("log(alpha) % bias")+
 annotate("text", x = Inf, y = Inf, label = paste("time-varying alpha and Smax "), vjust = 1, hjust = 1.2, size=6)
 pl_tvab_a_v


tvabbxp_smax<-df_pbias_mcmc_all %>% filter( paramch=="both"&
                            model %in%c("rwab","hmmab","simple",
                                "rwa","hmma", "rwb","hmmb")&
                            parameter %in% c("smax"))

tvabbxp_smax$model<-factor(tvabbxp_smax$model, levels=c("rwab","hmmab",
        "simple","rwa","hmma", "rwb","hmmb"))

pl_tvab_smax<-ggplot(tvabbxp_smax) + 
geom_boxplot(aes(x=model,y=bias, fill=modtype),outlier.shape = NA)+
coord_cartesian(ylim = quantile(tvabbxp_smax$bias, c(0.05, 0.95)))+
#facet_grid(rows=vars(parameter), scales="free_y")+
 scale_fill_viridis_d("estimation:",option = "E") +
mytheme+
geom_hline(yintercept=0, color="black",alpha=.7,linewidth=1.2) +
 ylab("Smax bias")+
 annotate("text", x = Inf, y = Inf, label = paste("time-varying alpha and Smax"), vjust = 1, hjust = 1.2, size=6)
 pl_tvab_smax


pl_tvab_smax_v<-ggplot(tvabbxp_smax) + 
geom_violin(aes(x=model,y=bias, fill=modtype),alpha=0.7)+
geom_boxplot(aes(x=model,y=bias, fill=modtype),outlier.shape = NA,width=0.1)+
coord_cartesian(ylim = quantile(tvabbxp_smax$bias, c(0.05, 0.95)))+
 scale_fill_viridis_d("estimation:",option = "E") +
mytheme+
geom_hline(yintercept=0, color="black",linewidth=1.2) +
 ylab("Smax bias")+
 annotate("text", x = Inf, y = Inf, label = paste("time-varying alpha and Smax"), vjust = 1, hjust = 1.2, size=6)
 pl_tvab_smax_v


pl_tvab_smax_v_scn<-ggplot(tvabbxp_smax) + 
geom_violin(aes(x=model,y=bias, fill=modtype),alpha=0.7)+
geom_boxplot(aes(x=model,y=bias, fill=modtype),outlier.shape = NA,width=0.1)+
coord_cartesian(ylim = quantile(tvabbxp_smax$bias, c(0.05, 0.95)))+
 scale_fill_viridis_d("estimation:",option = "E") +
mytheme+
   facet_wrap(~scenario)+
geom_hline(yintercept=0, color="black",linewidth=1.2) +
 ylab("Smax bias")+
 annotate("text", x = Inf, y = Inf, label = paste("time-varying alpha and Smax"), vjust = 1, hjust = 1.2, size=6)
 pl_tvab_smax_v_scn


pbias_tvab_smax_v<-ggplot(tvabbxp_smax) + 
geom_violin(aes(x=model,y=pbias, fill=modtype),alpha=0.7)+
geom_boxplot(aes(x=model,y=pbias, fill=modtype),outlier.shape = NA,width=0.1)+
coord_cartesian(ylim = c(-100, 100))+
 scale_fill_viridis_d("estimation:",option = "E") +
mytheme+
geom_hline(yintercept=0, color="black",linewidth=1.2) +
 ylab("Smax % bias")+
 annotate("text", x = Inf, y = Inf, label = paste("time-varying alpha and Smax"), vjust = 1, hjust = 1.2, size=6)
 pbias_tvab_smax_v

tvabbxp_smsy<-df_pbias_mcmc_all %>% filter( paramch=="both"&
                            model %in%c("rwab","hmmab","simple",
                                "rwa","hmma", "rwb","hmmb")&
                            parameter %in% c("smsy"))

tvabbxp_smsy$model<-factor(tvabbxp_smsy$model, levels=c("rwab","hmmab",
        "simple","rwa","hmma", "rwb","hmmb"))

pl_tvab_smsy<-ggplot(tvabbxp_smsy) + 
geom_boxplot(aes(x=model,y=bias, fill=modtype),outlier.shape = NA)+
coord_cartesian(ylim = quantile(tvabbxp_smsy$bias, c(0.05, 0.95)))+
#facet_grid(rows=vars(parameter), scales="free_y")+
 scale_fill_viridis_d("estimation:",option = "E") +
mytheme+
geom_hline(yintercept=0, color="black",linewidth=1.2) +
 ylab("Smsy bias")+
 annotate("text", x = Inf, y = Inf, label = paste("time-varying alpha and Smax"), vjust = 1, hjust = 1.2, size=6)
 pl_tvab_smax

 pl_tvab_smsy_v<-ggplot(tvabbxp_smsy) + 
geom_violin(aes(x=model,y=bias, fill=modtype),alpha=0.7)+
geom_boxplot(aes(x=model,y=bias, fill=modtype),outlier.shape = NA,width=0.1)+
coord_cartesian(ylim = quantile(tvabbxp_smsy$bias, c(0.05, 0.95)))+
 scale_fill_viridis_d("estimation:",option = "E") +
mytheme+
geom_hline(yintercept=0, color="black",linewidth=1.2) +
 ylab("Smsy bias")+
 annotate("text", x = Inf, y = Inf, label = paste("time-varying alpha and Smax"), vjust = 1, hjust = 1.2, size=6)
 pl_tvab_smsy_v


 pl_tvab_smsy_v_scn<-ggplot(tvabbxp_smsy) + 
geom_violin(aes(x=model,y=bias, fill=modtype),alpha=0.7)+
geom_boxplot(aes(x=model,y=bias, fill=modtype),outlier.shape = NA,width=0.1)+
coord_cartesian(ylim = quantile(tvabbxp_smsy$bias, c(0.05, 0.95)))+
 scale_fill_viridis_d("estimation:",option = "E") +
mytheme+
   facet_wrap(~scenario)+
geom_hline(yintercept=0, color="black",linewidth=1.2) +
 ylab("Smsy bias")+
 annotate("text", x = Inf, y = Inf, label = paste("time-varying alpha and Smax"), vjust = 1, hjust = 1.2, size=6)
 pl_tvab_smsy_v_scn


 pbias_tvab_smsy_v<-ggplot(tvabbxp_smsy) + 
geom_violin(aes(x=model,y=pbias, fill=modtype),alpha=0.7)+
geom_boxplot(aes(x=model,y=pbias, fill=modtype),outlier.shape = NA,width=0.1)+
coord_cartesian(ylim = c(-100, 100))+
 scale_fill_viridis_d("estimation:",option = "E") +
mytheme+
geom_hline(yintercept=0, color="black",linewidth=1.2) +
 ylab("Smsy bias")+
 annotate("text", x = Inf, y = Inf, label = paste("time-varying alpha and Smax"), vjust = 1, hjust = 1.2, size=6)
 pbias_tvab_smsy_v


plot_grid(pl_tvab_a,pl_tvab_smax,ncol=1)


plot_grid(pl_tvab_a_v,pl_tvab_smax_v, pl_tvab_smsy_v,ncol=3)


plot_grid(pbias_tvab_a_v,pbias_tvab_smax_v, pbias_tvab_smsy_v,ncol=3)


plot_grid(pl_tvab_a_v_scn,pl_tvab_smax_v_scn, pl_tvab_smsy_v_scn,ncol=1)



#===========================================
# stationary scenarios- no time-varying parameters  
unique(df_pbias_mcmc_all$paramch)


notvbxp<-df_pbias_mcmc_all %>% filter( paramch=="stationary"&
                            model %in%c("rwab","hmmab","simple",
                                "rwa","hmma", "rwb","hmmb")&
                            parameter %in%  c("alpha","smax","smsy"))

notvbxp$model<-factor(notvbxp$model, levels=c("rwab","hmmab",
        "simple","rwa","hmma", "rwb","hmmb"))



head(notvbxp)
aggregate(df_pbias_mcmc_all$iteration,
    list(df_pbias_mcmc_all$scenario,
        df_pbias_mcmc_all$model,
        df_pbias_mcmc_all$parameter), function(x){length(unique(x))})




pl_notv_v<-ggplot(notvbxp) + 
geom_violin(aes(x=model,y=pbias, fill=modtype),alpha=0.7, scale="width", trim=TRUE)+
geom_boxplot(aes(x=model,y=pbias, fill=modtype),outlier.shape = NA,width=0.1)+
coord_cartesian(ylim =c(-60,60))+
 scale_fill_viridis_d("estimation:",option = "E") +
mytheme+
facet_wrap(~parameter, ncol=3)+
geom_hline(yintercept=0, color="black",linewidth=1.2) +
 ylab("% bias")
 pl_notv_v



pl_notv_v_scn<-ggplot(notvbxp) + 
geom_violin(aes(x=model,y=pbias, fill=modtype),alpha=0.7, scale="width", trim=TRUE)+
geom_boxplot(aes(x=model,y=pbias, fill=modtype),outlier.shape = NA,width=0.1)+
coord_cartesian(ylim =c(-60,60))+
 scale_fill_viridis_d("estimation:",option = "E") +
mytheme+
facet_wrap(scenario~parameter, ncol=3)+
geom_hline(yintercept=0, color="black",linewidth=1.2) +
 ylab("% bias")
 pl_notv_v_scn


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

#========================================================================================================
#sensitivity smax scenario
#read in data
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



conv_iter<-aggregate(allconv$iteration,
    list(model=allconv$model,scenario=allconv$scenario),
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
unique(df_biasq$scenario)

df_biasq$paramch<-"a"
 
df_biasq$type<-"regime"
df_biasq$type[df_biasq$scenario%in%c("trendLinearSmax025", "trendLinearSmax050",
 "trendLinearSmax150", "trendLinearSmax200", "trendLinearSmax300")]<-"trend"





df_biasq$magnitude_ch<-"smax*0.25"
df_biasq$magnitude_ch[df_biasq$scenario%in%c("regimeSmax050")]<-"smax*0.5"    
df_biasq$magnitude_ch[df_biasq$scenario%in%c("regimeSmax150")]<-"smax*1.5"  
df_biasq$magnitude_ch[df_biasq$scenario%in%c("regimeSmax200")]<-"smax*2.0"  
df_biasq$magnitude_ch[df_biasq$scenario%in%c("regimeSmax300")]<-"smax*3.0"   
df_biasq$magnitude_ch[df_biasq$scenario%in%c("trendLinearSmax025")]<-"smax*0.25"   
df_biasq$magnitude_ch[df_biasq$scenario%in%c("trendLinearSmax050")]<-"smax*0.5"  
df_biasq$magnitude_ch[df_biasq$scenario%in%c("trendLinearSmax150")]<-"smax*1.5" 
df_biasq$magnitude_ch[df_biasq$scenario%in%c("trendLinearSmax200")]<-"smax*2.0" 
df_biasq$magnitude_ch[df_biasq$scenario%in%c("trendLinearSmax300")]<-"smax*3.0"  


df_biasq_MCMC<-df_biasq[df_biasq$method=="MCMC",]


df_alpha_biasq_MCMC<-df_biasq_MCMC[df_biasq_MCMC$parameter=="alpha",]
senssmax_alpha<-ggplot(df_alpha_biasq_MCMC) + 
coord_cartesian(ylim = c(-1,.5))+ 
geom_ribbon(aes(x=by,ymin =  x.10., ymax=x.90.,fill=magnitude_ch,group=scenario),alpha=0.2)+
geom_line(aes(x=by,y= x.50.,color=magnitude_ch,group=scenario),linewidth=1.2)+
       scale_color_viridis_d("magnitude of change:",begin=.1, end=.8) +
       scale_fill_viridis_d("magnitude of change:",begin=.1, end=.8) +
       geom_hline(yintercept=0, color="black",alpha=.7,linewidth=1.2) +
       mytheme+ 
       scale_y_continuous(name = "MCMC bias in log(alpha)", sec.axis = sec_axis(~., name = "varying parameters")) +
       theme(axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank()
       )+
       
       xlab("year") +
       facet_grid(type~model, scales="free_y")
senssmax_alpha  
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/trends/sens_smax/bias_a_MCMC.png", plot=senssmax_alpha  )



df_smax_biasq_MCMC<-df_biasq_MCMC[df_biasq_MCMC$parameter=="smax",]
senssmax_smax<-ggplot(df_smax_biasq_MCMC) + 
coord_cartesian(ylim = c(-180000,180000))+ 
geom_ribbon(aes(x=by,ymin =  x.10., ymax=x.90.,fill=magnitude_ch,group=scenario),alpha=0.2)+
geom_line(aes(x=by,y= x.50.,color=magnitude_ch,group=scenario),linewidth=1.2)+
       scale_color_viridis_d("magnitude of change:",begin=.1, end=.8) +
       scale_fill_viridis_d("magnitude of change:",begin=.1, end=.8) +
       geom_hline(yintercept=0, color="black",alpha=.7,linewidth=1.2) +
       mytheme+ 
       scale_y_continuous(name = "MCMC bias in Smax", sec.axis = sec_axis(~., name = "varying parameters")) +
       theme(axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank()
       )+
       #ylab("MLE bias in log(alpha)") +
       xlab("year") +
       facet_grid(type~model, scales="free_y")
senssmax_smax  
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/trends/sens_smax/bias_smax_MCMC.png",
 plot=senssmax_smax )



df_smsy_biasq_MCMC<-df_biasq_MCMC[df_biasq_MCMC$parameter=="smsy",]
senssmax_smsy<-ggplot(df_smsy_biasq_MCMC) + 
coord_cartesian(ylim = c(-80000,80000))+ 
geom_ribbon(aes(x=by,ymin =  x.10., ymax=x.90.,fill=magnitude_ch,group=scenario),alpha=0.2)+
geom_line(aes(x=by,y= x.50.,color=magnitude_ch,group=scenario),linewidth=1.2)+
       scale_color_viridis_d("magnitude of change:",begin=.1, end=.8) +
       scale_fill_viridis_d("magnitude of change:",begin=.1, end=.8) +
       geom_hline(yintercept=0, color="black",alpha=.7,linewidth=1.2) +
       mytheme+ 
       scale_y_continuous(name = "MCMC bias in Smsy", sec.axis = sec_axis(~., name = "varying parameters")) +
       theme(axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank()
       )+
       #ylab("MLE bias in log(alpha)") +
       xlab("year") +
       facet_grid(type~model, scales="free_y")
senssmax_smsy  
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/trends/sens_smax/bias_smsy_MCMC.png",
 plot=senssmax_smsy   )




#=================================================

df_biasq_MLE<-df_biasq[df_biasq$method=="MLE",]
head(df_biasq_MLE)
df_alpha_biasq_MLE<-df_biasq_MLE[df_biasq_MLE$parameter=="alpha",]
head(df_alpha_biasq_MLE)
senssmax_alpha_MLE<-ggplot(df_alpha_biasq_MLE) + 
coord_cartesian(ylim = c(-1,.5))+ 
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
senssmax_alpha_MLE  
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/trends/sens_smax/bias_a_MLE.png", 
    plot=senssmax_alpha_MLE )



df_smax_biasq_MLE<-df_biasq_MLE[df_biasq_MLE$parameter=="smax",]
senssmax_smax_MLE<- ggplot(df_smax_biasq_MLE) + 
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
senssmax_smax_MLE  
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/trends/sens_smax/bias_smax_MLE.png",
 plot=senssmax_smax_MLE  )


df_smsy_biasq_MLE<-df_biasq_MLE[df_biasq_MLE$parameter=="smsy",]
senssmax_smsy_MLE <-ggplot(df_smsy_biasq_MLE) +  
coord_cartesian(ylim = c(-80000,80000))+ 
geom_ribbon(aes(x=by,ymin =  x.10., ymax=x.90.,fill=magnitude_ch,group=scenario),alpha=0.2)+
geom_line(aes(x=by,y= x.50.,color=magnitude_ch,group=scenario),linewidth=1.2)+
       scale_color_viridis_d(begin=.1, end=.8) +
       scale_fill_viridis_d(begin=.1, end=.8) +
       geom_hline(yintercept=0, color="black",alpha=.7,linewidth=1.2) +
       mytheme+ 
       ylab("MLE bias in Smsy") +
       xlab("year") +
       facet_grid(type~model, scales="free_y")
senssmax_smsy_MLE
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/trends/sens_smax/bias_smsy_MLE.png",
 plot=senssmax_smsy_MLE  )


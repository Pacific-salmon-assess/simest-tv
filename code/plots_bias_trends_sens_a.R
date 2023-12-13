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
convstatMLE<-convstat[convstat$x==0&convstat$method=="MLE",]
convstatMCMC<-convstat[convstat$x==0&convstat$method=="MCMC",]

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
       scale_y_continuous(name = "MCMC bias in log(alpha)", sec.axis = sec_axis(~., name = "varying parameters")) +
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
       scale_y_continuous(name = "MCMC bias in Smax", sec.axis = sec_axis(~., name = "varying parameters")) +
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
       scale_y_continuous(name = "MCMC bias in Smsy", sec.axis = sec_axis(~., name = "varying parameters")) +
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


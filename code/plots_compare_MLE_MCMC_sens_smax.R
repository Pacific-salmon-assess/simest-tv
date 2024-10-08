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
              legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="bold", size=12),
              axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
)

===========================================================================
#=================================================================================================================
#sensitivity smax
simPar <- read.csv("data/Smax_sensitivity/SimPars.csv")

## Store relevant object names to help run simulation 
scenNames <- unique(simPar$scenario)

restmb<-readRDS(file = "outs/simest/Smax_sensitivity/res_smax.rds")

resstansmax1<-readRDS(file = "outs/simest/Smax_sensitivity/resstan_smax1.rds")
resstansmax2<-readRDS(file = "outs/simest/Smax_sensitivity/resstan_smax2.rds")

resstan<-rbind(resstansmax1,resstansmax2)

res<-rbind(restmb,resstan)

res$parameter[res$parameter=="Smax"]<-"smax"
res$parameter[res$parameter=="alpha"]<-"logalpha"

resparam<-res[res$parameter%in%c("logalpha","smax","sigma","smsy","sgen","umsy"),]


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

head(resparam)

df<-reshape2::melt(resparam, id.vars=c("parameter","iteration","scenario","method","model","by", "convergence","pbias","bias"))

#df_alpha<-df[df$parameter%in%c("alpha"),]
df$col<-factor(df$variable,levels=c("median","mode", "sim"))


summarydf<-aggregate(df$value,by=list(scenario=df$scenario,
    parameter=df$parameter,  
    method=df$method, 
    model=df$model,
    by=df$by, 
    variable=df$variable),
    function(x) {quantile(x,probs = c(0.025, .5, 0.975),na.rm=T)})
summarydf<-do.call(data.frame, summarydf)
head(summarydf)
unique(summarydf$scenario)

unique(summarydf$scenario)



summarydf$scenario<-factor(summarydf$scenario,levels=c("regimeSmax025",
                                                       "regimeSmax050",
                                                       "regimeSmax150",     
                                                       "regimeSmax200",
                                                       "regimeSmax300",
                                                       "trendLinearSmax025",
                                                       "trendLinearSmax050",
                                                       "trendLinearSmax150",
                                                       "trendLinearSmax200",
                                                       "trendLinearSmax300"  ))


summarydf$typechange<-dplyr::case_match(summarydf$scenario, 
      "regimeSmax025"~"regime",
        "regimeSmax050"~"regime",
        "regimeSmax150"~"regime",     
        "regimeSmax200"~"regime",
        "regimeSmax300"~"regime",
        "trendLinearSmax025"~"trend",
        "trendLinearSmax050"~"trend",
        "trendLinearSmax150"~"trend",
        "trendLinearSmax200"~"trend",
        "trendLinearSmax300"~"trend"
      )   

unique(summarydf$typechange)

summarydf$magnitude_ch<-0

summarydf$magnitude_ch[summarydf$scenario%in%c("regimeSmax025")]<-"S[max]~x~0.25"  
summarydf$magnitude_ch[summarydf$scenario%in%c("regimeSmax050")]<-"S[max]~x~0.5"    
summarydf$magnitude_ch[summarydf$scenario%in%c("regimeSmax150")]<-"S[max]~x~1.5"  
summarydf$magnitude_ch[summarydf$scenario%in%c("regimeSmax200")]<-"S[max]~x~2.0"  
summarydf$magnitude_ch[summarydf$scenario%in%c("regimeSmax300")]<-"S[max]~x~3.0"   
summarydf$magnitude_ch[summarydf$scenario%in%c("trendLinearSmax025")]<-"S[max]~x~0.25"   
summarydf$magnitude_ch[summarydf$scenario%in%c("trendLinearSmax050")]<-"S[max]~x~0.5"  
summarydf$magnitude_ch[summarydf$scenario%in%c("trendLinearSmax150")]<-"S[max]~x~1.5" 
summarydf$magnitude_ch[summarydf$scenario%in%c("trendLinearSmax200")]<-"S[max]~x~2.0" 
summarydf$magnitude_ch[summarydf$scenario%in%c("trendLinearSmax300")]<-"S[max]~x~3.0"   

summarydf$magnitude_ch<-factor(summarydf$magnitude_ch,levels=c("S[max]~x~0.25",
"S[max]~x~0.5",
"S[max]~x~1.5", 
"S[max]~x~2.0", 
"S[max]~x~3.0"))

unique(summarydf$magnitude_ch)

summarydf$model<-factor(summarydf$model, levels=c("simple", "autocorr", "rwa",  "hmma", "rwb", "hmmb", "rwab", "hmmab"))

summarydf_alpha_sim<- summarydf[summarydf$parameter=="logalpha"&summarydf$variable=="sim",]

summarydf_alpha<- summarydf[summarydf$parameter=="logalpha"&summarydf$variable=="mode",]


p_sens_smax_alpha<-ggplot() + 
geom_pointrange(data=summarydf_alpha,aes(x=by-54,y= x.50.,ymin = x.2.5., ymax = x.97.5., col=method),alpha=.6)+
geom_line(data=summarydf_alpha_sim,aes(x=by-54,y= x.50.),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
coord_cartesian(ylim = c(-0.2,2.0))+ 
mytheme+ 
ylab(expression(log(alpha))) +
xlab("year") +
facet_grid(typechange+magnitude_ch~model, scales="free_y", labeller =  label_parsed)
p_sens_smax_alpha

ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/sens_smax/compareMCMC_MLE_sens_smax_alpha.png",
    plot=p_sens_smax_alpha, width = 15,height = 18 )






summarydf_alpha_sim1<-summarydf_alpha_sim[summarydf_alpha_sim$scenario%in%c( "regimeSmax025",
                                                       "regimeSmax050",
                                                        "trendLinearSmax025",
                                                       "trendLinearSmax050" ),]
summarydf_alpha1<-summarydf_alpha[summarydf_alpha$scenario%in%c( "regimeSmax025",
                                                       "regimeSmax050",
                                                        "trendLinearSmax025",
                                                       "trendLinearSmax050"),]
head(summarydf_alpha1)
head(summarydf_alpha_sim1)

unique(summarydf_alpha1$magnitude_ch)

p_sens_smax_alpha1<-ggplot() + 
geom_pointrange(data=summarydf_alpha1,aes(x=by,y= x.50.,ymin = x.2.5., ymax = x.97.5., col=method),alpha=.6)+
geom_line(data=summarydf_alpha_sim1,aes(x=by,y= x.50.),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
coord_cartesian(ylim = c(-0.2,2.0))+ 
mytheme+ 
ylab(expression("log(alpha)")) +
xlab("year") +
facet_grid(magnitude_ch~model, scales="free_y")
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/sens_smax/compareMCMC_MLE_sens_smax_alpha1.png",
    plot=p_sens_smax_alpha1)





summarydf_alpha_sim2<-summarydf_alpha_sim[summarydf_alpha_sim$scenario%in%c( "regimeSmax150",     
                                                       "regimeSmax200",
                                                       "regimeSmax300",
                                                       "trendLinearSmax150",
                                                       "trendLinearSmax200",
                                                       "trendLinearSmax300"),]
summarydf_alpha2<-summarydf_alpha[summarydf_alpha$scenario%in%c("regimeSmax150",     
                                                       "regimeSmax200",
                                                       "regimeSmax300",
                                                       "trendLinearSmax150",
                                                       "trendLinearSmax200",
                                                       "trendLinearSmax300"),]
  

head(summarydf_alpha2)
head(summarydf_alpha_sim2)  
                                                                
p_sens_smax_alpha2<-ggplot() + 
geom_pointrange(data=summarydf_alpha2,aes(x=by,y= x.50.,ymin = x.2.5., ymax = x.97.5., col=method),alpha=.6)+
geom_line(data=summarydf_alpha_sim2,aes(x=by,y= x.50.),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
coord_cartesian(ylim = c(0.5,1.8))+ 
mytheme + 
ylab("log(alpha)") +
xlab("year") +
facet_grid(magnitude_ch~model, scales="free_y")
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/sens_smax/compareMCMC_MLE_sens_smax_alpha2.png",
    plot=p_sens_smax_alpha2, width = 15,height = 18)
#MLE estimates are less biased and higher than MCMC


#smax



summarydf_smax_sim<- summarydf[summarydf$parameter=="smax"&summarydf$variable=="sim",]

summarydf_smax<- summarydf[summarydf$parameter=="smax"&summarydf$variable=="mode",]


p_sens_smax_smax<-ggplot() + 
geom_pointrange(data=summarydf_smax,aes(x=by-54,y= x.50./1000,ymin = x.2.5./1000,
 ymax = x.97.5./1000, col=method),alpha=.6)+
geom_line(data=summarydf_smax_sim,aes(x=by-54,y= x.50./1000),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
coord_cartesian(ylim = c(10000,670000)/1000)+ #570000
mytheme+ 
ylab(expression(S[max]~"(1000s)")) +
xlab("year") +
facet_grid(typechange+magnitude_ch~model, scales="free_y", labeller =  label_parsed)
p_sens_smax_smax
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/sens_smax/compareMCMC_MLE_sens_smax_smax.png",
    plot=p_sens_smax_smax, width = 15,height = 18)


summarydf_smax_sim1<-summarydf_smax_sim[summarydf_smax_sim$scenario%in%c( "regimeSmax025",
                                                       "regimeSmax050",
                                                        "trendLinearSmax025",
                                                       "trendLinearSmax050"),]
summarydf_smax1<-summarydf_smax[summarydf_smax$scenario%in%c( "regimeSmax025",
                                                       "regimeSmax050",
                                                       "trendLinearSmax025",
                                                       "trendLinearSmax050"),]
head(summarydf_smax1)

p_sens_smax_smax1<-ggplot() + 
geom_pointrange(data=summarydf_smax1,aes(x=by,y= x.50.,ymin = x.2.5., ymax = x.97.5., col=method),alpha=.6)+
geom_line(data=summarydf_smax_sim1,aes(x=by,y= x.50.),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
coord_cartesian(ylim = c(10000,300000))+ #570000
mytheme+ 
ylab("smax") +
xlab("year") +
facet_grid(magnitude_ch~model, scales="free_y")
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/sens_smax/compareMCMC_MLE_sens_smax_smax1.png",
    plot=p_sens_smax_smax1)




summarydf_smax_sim2<-summarydf_smax_sim[summarydf_smax_sim$scenario%in%c("regimeSmax150",     
                                                       "regimeSmax200",
                                                       "regimeSmax300",
                                                       "trendLinearSmax150",
                                                       "trendLinearSmax200",
                                                       "trendLinearSmax300" ),]
summarydf_smax2<-summarydf_smax[summarydf_smax$scenario%in%c( "regimeSmax150",     
                                                       "regimeSmax200",
                                                       "regimeSmax300",
                                                       "trendLinearSmax150",
                                                       "trendLinearSmax200",
                                                       "trendLinearSmax300"),]

head(summarydf_smax2)
head(summarydf_smax_sim2)
p_sens_smax_smax2<-ggplot() + 
geom_pointrange(data=summarydf_smax2,aes(x=by,y= x.50.,ymin = x.2.5., ymax = x.97.5., col=method),alpha=.6)+
geom_line(data=summarydf_smax_sim2,aes(x=by,y= x.50.),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
coord_cartesian(ylim = c(10000,670000))+ 
mytheme + 
ylab("smax") +
xlab("year") +
facet_grid(magnitude_ch~model, scales="free_y")
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/sens_smax/compareMCMC_MLE_sens_smax_smax2.png",
    plot=p_sens_smax_smax2)
#MLE estimates are less biased and higher than MCMC



#=======================================================================
#smsy


summarydf_smsy_sim<- summarydf[summarydf$parameter=="smsy"&summarydf$variable=="sim",]

summarydf_smsy<- summarydf[summarydf$parameter=="smsy"&summarydf$variable=="mode",]


p_sens_smax_smsy<-ggplot() + 
geom_pointrange(data=summarydf_smsy,aes(x=by-54,y= x.50./1000,ymin = x.2.5./1000,
 ymax = x.97.5./1000, col=method),alpha=.6)+
geom_line(data=summarydf_smsy_sim,aes(x=by-54,y= x.50./1000),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
coord_cartesian(ylim = c(0,350000)/1000)+ #570000
mytheme+ 
ylab(expression(S[MSY]~"(1000s)")) +
xlab("year") +
facet_grid(typechange+magnitude_ch~model, scales="free_y", labeller =  label_parsed)
p_sens_smax_smsy
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/sens_smax/compareMCMC_MLE_sens_smax_smsy.png",
    plot=p_sens_smax_smsy, width = 15,height = 18)



summarydf_smsy_sim1<-summarydf_smsy_sim[summarydf_smsy_sim$scenario%in%c( "regimeSmax025",
                                                       "regimeSmax050",
                                                        "trendLinearSmax025",
                                                       "trendLinearSmax050" ),]
summarydf_smsy1<-summarydf_smsy[summarydf_smsy$scenario%in%c( "regimeSmax025",
                                                       "regimeSmax050",
                                                        "trendLinearSmax025",
                                                       "trendLinearSmax050" ),]
head(summarydf_smsy1)


p_sens_smax_smsy1<-ggplot() + 
geom_pointrange(data=summarydf_smsy1,aes(x=by,y= x.50.,ymin = x.2.5., ymax = x.97.5., col=method),alpha=.6)+
geom_line(data=summarydf_smsy_sim1,aes(x=by,y= x.50.),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
coord_cartesian(ylim = c(0,90000))+ 
mytheme+ 
ylab("smsy") +
xlab("year") +
facet_grid(magnitude_ch~model, scales="free_y")
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/sens_smax/compareMCMC_MLE_sens_smax_smsy1.png",
    plot=p_sens_smax_smsy1)



summarydf_smsy_sim2<-summarydf_smsy_sim[summarydf_smsy_sim$scenario%in%c("regimeSmax150",     
                                                       "regimeSmax200",
                                                       "regimeSmax300",
                                                       "trendLinearSmax150",
                                                       "trendLinearSmax200",
                                                       "trendLinearSmax300"),]
summarydf_smsy2<-summarydf_smsy[summarydf_smsy$scenario%in%c("regimeSmax150",     
                                                       "regimeSmax200",
                                                       "regimeSmax300",
                                                       "trendLinearSmax150",
                                                       "trendLinearSmax200",
                                                       "trendLinearSmax300"),]

head(summarydf_smsy2[summarydf_smsy2$scenario=="trendLinearSmax150",])

unique(summarydf_smsy2$scenario)

head(summarydf_smsy2)
head(summarydf_smsy_sim2)
p_sens_smax_smsy2<-ggplot() + 
geom_pointrange(data=summarydf_smsy2,aes(x=by,y= x.50.,ymin = x.2.5., ymax = x.97.5., col=method),alpha=.6)+
geom_line(data=summarydf_smsy_sim2,aes(x=by,y= x.50.),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
coord_cartesian(ylim = c(10000,670000))+ 
mytheme + 
ylab("smsy") +
xlab("year") +
facet_grid(magnitude_ch~model, scales="free_y")
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/sens_smax/compareMCMC_MLE_sens_smax_smsy2.png",
    plot=p_sens_smax_smsy2)
#MLE estimates are less biased and higher than MCMC




#=======================================================================
#umsy -- these are wrong need to be checked
#hmma and hmmab are wrong for stan estimates, fixed function but have not re-run


summarydf_umsy_sim<- summarydf[summarydf$parameter=="umsy"&summarydf$variable=="sim",]





summarydf_umsy<- summarydf[summarydf$parameter=="umsy"&summarydf$variable=="mode",]



summarydf_umsy_sim1<-summarydf_umsy_sim[summarydf_umsy_sim$scenario%in%c( "regimeSmax025",
                                                       "regimeSmax050",
                                                        "trendLinearSmax025",
                                                       "trendLinearSmax050"),]
summarydf_umsy1<-summarydf_umsy[summarydf_umsy$scenario%in%c( "regimeSmax025",
                                                       "regimeSmax050",
                                                        "trendLinearSmax025",
                                                       "trendLinearSmax050"),]

             

ggplot() + 
geom_pointrange(data=summarydf_umsy1,aes(x=by,y= x.50.,ymin = x.2.5., ymax = x.97.5., col=method),alpha=.6)+
geom_line(data=summarydf_umsy_sim1,aes(x=by,y= x.50.),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
coord_cartesian(ylim = c(0.3,.9))+ 
mytheme+ 
ylab("umsy") +
xlab("year") +
facet_grid(scenario~model, scales="free_y")
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/sens_smax/compareMCMC_MLE_sens_smax_umsy1.png")



summarydf_umsy_sim2<-summarydf_umsy_sim[summarydf_umsy_sim$scenario%in%c("regimeSmax150",     
                                                       "regimeSmax200",
                                                       "regimeSmax300",
                                                       "trendLinearSmax150",
                                                       "trendLinearSmax200",
                                                       "trendLinearSmax300"),]
summarydf_umsy2<-summarydf_umsy[summarydf_umsy$scenario%in%c("regimeSmax150",     
                                                       "regimeSmax200",
                                                       "regimeSmax300",
                                                       "trendLinearSmax150",
                                                       "trendLinearSmax200",
                                                       "trendLinearSmax300"),]

head(summarydf_smsy2)
head(summarydf_smsy_sim2)
ggplot() + 
geom_pointrange(data=summarydf_umsy2,aes(x=by,y= x.50.,ymin = x.2.5., ymax = x.97.5., col=method),alpha=.6)+
geom_line(data=summarydf_umsy_sim2,aes(x=by,y= x.50.),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
coord_cartesian(ylim = c(0.3,.9))+ 
mytheme + 
ylab("alpha") +
xlab("year") +
facet_grid(magnitude_ch~model, scales="free_y")
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/sens_smax/compareMCMC_MLE_sens_smax_umsy2.png")
#MLE estimates are less biased and higher than MCMC

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




source("code/read_sensa_data.R")
#=================================================================================================================


df<-reshape2::melt(resparam, id.vars=c("parameter","iteration","scenario","method","model","by", "convergence","conv_warning","pbias","bias"))

#df_alpha<-df[df$parameter%in%c("alpha"),]
df$col<-factor(df$variable,levels=c("median","mode", "sim"))

df$model<-factor(df$model,levels=c("simple",
                                   "autocorr", 
                                   "rwa",
                                   "hmma",
                                   "rwb",
                                   "hmmb",
                                   "rwab",
                                   "hmmab"))


summarydf  <- df %>%
   group_by(scenario,parameter,
    method,model,by,variable, col) %>%
   reframe(qs = quantile(value, c(0.025, .5, 0.975),na.rm=T), prob = c("lower","median", "upper"))



summarydf$scenario<-factor(summarydf$scenario,levels=c( "trendLinearProd1" ,
                                                        "trendLinearProd1.3",
                                                        "trendLinearProd2", 
                                                        "trendLinearProd5", 
                                                        "trendLinearProd7",
                                                        "trendLinearProd10" ,    
                                                        "regimeProd1",
                                                        "regimeProd1.3",  
                                                        "regimeProd2" ,     
                                                        "regimeProd5",     
                                                        "regimeProd7",
                                                        "regimeProd10" ))




summarydf$magnitude_ch<-"a 3.7 -> 1.03"
summarydf$magnitude_ch[summarydf$scenario%in%c("trendLinearProd1")]<-"a 3.7 -> 1.03"

summarydf$magnitude_ch[summarydf$scenario%in%c("regimeProd1.3")]<-"a 3.7 -> 1.35"    
summarydf$magnitude_ch[summarydf$scenario%in%c("trendLinearProd1.3")]<-"a 3.7 -> 1.35"    

summarydf$magnitude_ch[summarydf$scenario%in%c("trendLinearProd10")]<-"a 3.7 -> 10"  
summarydf$magnitude_ch[summarydf$scenario%in%c( "regimeProd10")]<-"a 3.7 -> 10"  

summarydf$magnitude_ch[summarydf$scenario%in%c("trendLinearProd2" )]<-"a 3.7 -> 2"   
summarydf$magnitude_ch[summarydf$scenario%in%c("regimeProd2")]<-"a 3.7 -> 2"   

summarydf$magnitude_ch[summarydf$scenario%in%c("trendLinearProd5")]<-"a 3.7 -> 5"  
summarydf$magnitude_ch[summarydf$scenario%in%c("regimeProd5")]<-"a 3.7 -> 5" 

summarydf$magnitude_ch[summarydf$scenario%in%c("regimeProd7")]<-"a 3.7 -> 7" 
summarydf$magnitude_ch[summarydf$scenario%in%c("trendLinearProd7")]<-"a 3.7 -> 7"               

summarydf$magnitude_ch<-factor(summarydf$magnitude_ch,levels=c("a 3.7 -> 1.03",
       "a 3.7 -> 1.35", 
       "a 3.7 -> 2",
       "a 3.7 -> 5",
       "a 3.7 -> 7", 
       "a 3.7 -> 10"))


summarydf <- reshape2::dcast(data=summarydf,  scenario + parameter + method + model + by + variable  + magnitude_ch ~prob, 
    value.var= "qs",fun.aggregate=mean)

summarydf_alpha_sim<- summarydf[summarydf$parameter=="alpha"&summarydf$variable=="sim",]

summarydf_alpha<- summarydf[summarydf$parameter=="alpha"&summarydf$variable=="mode",]


summarydf_alpha_sim1<-summarydf_alpha_sim[summarydf_alpha_sim$scenario%in%c( "regimeProd1", 
                                                                    "regimeProd1.3", 
                                                                    "regimeProd2" ,     
                                                                    "regimeProd5",     
                                                                    "regimeProd7",
                                                                    "regimeProd10" ),]
summarydf_alpha1<-summarydf_alpha[summarydf_alpha$scenario%in%c("regimeProd1", 
                                                                "regimeProd1.3", 
                                                                "regimeProd2" ,     
                                                                "regimeProd5",     
                                                                "regimeProd7",
                                                                "regimeProd10" ),]

p_sensa_alpha1 <- ggplot() + 
geom_pointrange(data=summarydf_alpha1,aes(x=by,y=median ,ymin =lower , ymax = upper, col=method),alpha=.6)+
geom_line(data=summarydf_alpha_sim1,aes(x=by,y=median),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
coord_cartesian(ylim = c(-3.2,2.5))+ 
mytheme+ 
ylab("alpha") +
xlab("year") +
ggtitle("alpha regime") +
facet_grid(magnitude_ch~model, scales="free_y")
p_sensa_alpha1

ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/sens_a/compareMCMC_MLE_sensa_alpha1.png", plot=p_sensa_alpha1)



                                                       

p_sensa_alpha1_zoom<-ggplot() + 
geom_pointrange(data=summarydf_alpha1,aes(x=by,y=median ,ymin =lower , ymax = upper, col=method),alpha=.6)+
geom_line(data=summarydf_alpha_sim1,aes(x=by, y=median),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
coord_cartesian(ylim = c(-0.1,2.5))+ 
mytheme+ 
ylab("alpha") +
xlab("year") +
ggtitle("alpha regime") +
facet_grid(magnitude_ch~model, scales="free_y")
p_sensa_alpha1_zoom
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/sens_a/compareMCMC_MLE_sensa_alpha1_zoom.png",
 plot=p_sensa_alpha1_zoom)






summarydf_alpha_sim2<-summarydf_alpha_sim[summarydf_alpha_sim$scenario%in%c( "trendLinearProd1" ,
                                                                    "trendLinearProd1.3",
                                                                    "trendLinearProd2",
                                                                    "trendLinearProd5", 
                                                                    "trendLinearProd7",
                                                                    "trendLinearProd10" ),]
summarydf_alpha2<-summarydf_alpha[summarydf_alpha$scenario%in%c( "trendLinearProd1" ,
                                                                 "trendLinearProd1.3",
                                                                    "trendLinearProd2", 
                                                                    "trendLinearProd5", 
                                                                    "trendLinearProd7",
                                                                    "trendLinearProd10" ),]
  

 
                                                                
p_sensa_alpha2<-ggplot() + 
geom_pointrange(data=summarydf_alpha2,aes(x=by,y=median ,ymin =lower , ymax = upper, col=method),alpha=.6)+
geom_line(data=summarydf_alpha_sim2,aes(x=by,y= median),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
coord_cartesian(ylim = c(-0.4,2.7))+ 
mytheme + 
ylab("alpha") +
xlab("year") +
facet_grid(magnitude_ch~model, scales="free_y")
p_sensa_alpha2
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/sens_a/compareMCMC_MLE_sensa_alpha2.png", plot=p_sensa_alpha2)
#MLE estimates are less biased and higher than MCMC


#smax



summarydf_smax_sim<- summarydf[summarydf$parameter=="smax"&summarydf$variable=="sim",]

summarydf_smax<- summarydf[summarydf$parameter=="smax"&summarydf$variable=="mode",]




summarydf_smax_sim1<-summarydf_smax_sim[summarydf_smax_sim$scenario%in%c( "regimeProd1",  
                                                                    "regimeProd2" ,     
                                                                    "regimeProd5",     
                                                                    "regimeProd7",
                                                                    "regimeProd10" ),]
summarydf_smax1<-summarydf_smax[summarydf_smax$scenario%in%c( "regimeProd1",  
                                                                    "regimeProd2" ,     
                                                                    "regimeProd5",     
                                                                    "regimeProd7",
                                                                    "regimeProd10" ),]
head(summarydf_smax1)

p_sensa_smax1<-ggplot() + 
geom_pointrange(data=summarydf_smax1,aes(x=by,y=median ,ymin =lower , ymax = upper, col=method),alpha=.6)+
geom_line(data=summarydf_smax_sim1,aes(x=by,y= median),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
coord_cartesian(ylim = c(10000,370000))+ 
mytheme+ 
ylab("smax") +
xlab("year") +
facet_grid(magnitude_ch~model, scales="free_y")
p_sensa_smax1
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/sens_a/compareMCMC_MLE_sensa_smax1.png", 
    plot=p_sensa_smax1)



summarydf_smax_sim2<-summarydf_smax_sim[summarydf_smax_sim$scenario%in%c( "trendLinearProd1" ,
                                                                    "trendLinearProd1.3",
                                                                    "trendLinearProd2", 
                                                                    "trendLinearProd5", 
                                                                    "trendLinearProd7",
                                                                    "trendLinearProd10" ),]
summarydf_smax2<-summarydf_smax[summarydf_smax$scenario%in%c( "trendLinearProd1" ,
                                                                    "trendLinearProd1.3",
                                                                    "trendLinearProd2", 
                                                                    "trendLinearProd5", 
                                                                    "trendLinearProd7",
                                                                    "trendLinearProd10" ),]


p_sensa_smax2<-ggplot() + 
geom_pointrange(data=summarydf_smax2,aes(x=by,y=median ,ymin =lower , ymax = upper, col=method),alpha=.6)+
geom_line(data=summarydf_smax_sim2,aes(x=by,y= median),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
coord_cartesian(ylim = c(10000,470000))+ 
mytheme + 
ylab("smax") +
xlab("year") +
facet_grid(magnitude_ch~model, scales="free_y")
p_sensa_smax2
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/sens_a/compareMCMC_MLE_sensa_smax2.png",
 plot=p_sensa_smax2)
#MLE estimates are less biased and higher than MCMC



#=======================================================================
#smsy


summarydf_smsy_sim<- summarydf[summarydf$parameter=="smsy"&summarydf$variable=="sim",]

summarydf_smsy<- summarydf[summarydf$parameter=="smsy"&summarydf$variable=="mode",]




summarydf_smsy_sim1<-summarydf_smsy_sim[summarydf_smsy_sim$scenario%in%c( "regimeProd1", 
                                                                    "regimeProd2" ,     
                                                                    "regimeProd5",     
                                                                    "regimeProd7",
                                                                    "regimeProd10" ),]
summarydf_smsy1<-summarydf_smsy[summarydf_smsy$scenario%in%c( "regimeProd1", 
                                                                    "regimeProd2" ,     
                                                                    "regimeProd5",     
                                                                    "regimeProd7",
                                                                    "regimeProd10" ),]
head(summarydf_smsy1)

p_sensa_smsy1<-ggplot() + 
geom_pointrange(data=summarydf_smsy1,aes(x=by,y=median ,ymin =lower , ymax = upper, col=method),alpha=.6)+
geom_line(data=summarydf_smsy_sim1,aes(x=by,y= median),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
coord_cartesian(ylim = c(-100000,170000))+ 
mytheme+ 
ylab("smsy") +
xlab("year") +
facet_grid(magnitude_ch~model, scales="free_y")
p_sensa_smsy1
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/sens_a/compareMCMC_MLE_sensa_smsy1.png", 
    plot=p_sensa_smsy1)



summarydf_smsy_sim2<-summarydf_smax_sim[summarydf_smsy_sim$scenario%in%c( "trendLinearProd1" ,
                                                                    "trendLinearProd1.3", 
                                                                    "trendLinearProd2", 
                                                                    "trendLinearProd5", 
                                                                    "trendLinearProd7",
                                                                    "trendLinearProd10" ),]
summarydf_smsy2<-summarydf_smax[summarydf_smsy$scenario%in%c( "trendLinearProd1" ,
                                                             "trendLinearProd1.3", 
                                                                    "trendLinearProd2", 
                                                                    "trendLinearProd5", 
                                                                    "trendLinearProd7",
                                                                    "trendLinearProd10" ),]

head(summarydf_smsy2)
head(summarydf_smsy_sim2)
p_sensa_smsy2<-ggplot() + 
geom_pointrange(data=summarydf_smsy2,aes(x=by,y=median ,ymin =lower , ymax = upper, col=method),alpha=.6)+
geom_line(data=summarydf_smsy_sim2,aes(x=by,y= median),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
coord_cartesian(ylim = c(10000,470000))+ 
mytheme + 
ylab("smsy") +
xlab("year") +
facet_grid(magnitude_ch~model, scales="free_y")
p_sensa_smsy2
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/sens_a/compareMCMC_MLE_sensa_smsy2.png", 
    plot=p_sensa_smsy2)
#MLE estimates are less biased and higher than MCMC





#=======================================================================
#umsy -- these are wrong need to be checked



summarydf_umsy_sim<- summarydf[summarydf$parameter=="umsy"& 
                              summarydf$variable=="sim"&
                              summarydf$method=="HMC",]

summarydf_umsy<- summarydf[summarydf$parameter=="umsy"&summarydf$variable=="mode",]





summarydf_umsy_sim1<-summarydf_umsy_sim[summarydf_umsy_sim$scenario%in%c( "regimeProd1",  
                                                                    "regimeProd2" ,     
                                                                    "regimeProd5",     
                                                                    "regimeProd7",
                                                                    "regimeProd10" ),]
summarydf_umsy1<-summarydf_umsy[summarydf_umsy$scenario%in%c( "regimeProd1",  
                                                                    "regimeProd2" ,     
                                                                    "regimeProd5",     
                                                                    "regimeProd7",
                                                                    "regimeProd10" ),]




p_sensa_umsy1<-ggplot() + 
geom_pointrange(data=summarydf_umsy1,aes(x=by,y=median ,ymin =lower , ymax = upper, col=method),alpha=.6)+
geom_line(data=summarydf_umsy_sim1,aes(x=by,y= median),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
coord_cartesian(ylim = c(-.3,.9))+ 
mytheme+ 
ylab("umsy") +
xlab("year") +
facet_grid(magnitude_ch~model,, scales="free_y")
p_sensa_umsy1
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/sens_a/compareMCMC_MLE_sensa_umsy1.png",plot=p_sensa_umsy1)



summarydf_umsy_sim2<-summarydf_umsy_sim[summarydf_umsy_sim$scenario%in%c( "trendLinearProd1" ,
                                                                    "trendLinearProd1.3", 
                                                                    "trendLinearProd2", 
                                                                    "trendLinearProd5", 
                                                                    "trendLinearProd7",
                                                                    "trendLinearProd10" ),]
summarydf_umsy2<-summarydf_umsy[summarydf_umsy$scenario%in%c( "trendLinearProd1" ,
                                                              "trendLinearProd1.3", 
                                                                    "trendLinearProd2", 
                                                                    "trendLinearProd5", 
                                                                    "trendLinearProd7",
                                                                    "trendLinearProd10" ),]

p_sensa_umsy2<-ggplot() + 
geom_pointrange(data=summarydf_umsy2,aes(x=by,y=median ,ymin =lower , ymax = upper, col=method),alpha=.6)+
geom_line(data=summarydf_umsy_sim2,aes(x=by,y= median),color="black", alpha=.6,linewidth=1.2)+
scale_color_viridis_d(begin=.1, end=.8) +
scale_fill_viridis_d(begin=.1, end=.8) +
coord_cartesian(ylim = c(-.1,.9))+ 
mytheme + 
ylab("alpha") +
xlab("year") +
facet_grid(magnitude_ch~model, scales="free_y")
p_sensa_umsy2
ggsave("../Best-Practices-time-varying-salmon-SR-models/figures/MCMC_MLE_comp/sens_a/compareMCMC_MLE_sensa_umsy2.png",
 plot=p_sensa_umsy2)
#MLE estimates are less biased and higher than MCMC


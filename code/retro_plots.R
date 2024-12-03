library(ggplot2)
library(dplyr)

mytheme = list(
  theme_classic(16)+
    theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
          legend.position="right", strip.text = element_text(face="bold", size=12),
          axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=15))
)

simPar <- read.csv("data/generic/SimPars.csv")

res<-readRDS(file = "outs/res_retro.rds")
res$parameter[res$parameter=="Smax"]<-"smax"
res<- res[res$conv_warning==0,]

resparam<-res[res$parameter%in%c("logalpha","smax","smsy","sgen","umsy"),]
resparam<-resparam[resparam$model%in%c("autocorr","hmma","hmmb","rwa","rwb"),]
resparam$parameter=factor(resparam$parameter)
resparam$parameter=droplevels(resparam$parameter)
resparam$model=factor(resparam$model)
resparam$model=droplevels(resparam$model)

#resparam<- resparam[resparam$conv_warning==0,]
res2<- resparam %>% group_by(scenario,parameter,by,endyr,model) %>% summarize(median.est=median(mode,na.rm=T),q10.est=quantile(mode,0.1,na.rm=T),q90.est=quantile(mode,0.9,na.rm=T))


#logalpha####
la_retro<-res2[res2$parameter=="logalpha"&res2$model %in%c('autocorr','hmma','rwa'),]
la_retro<- subset(la_retro,endyr %in% seq(min(la_retro$endyr),max(la_retro$endyr),by=2))
la_retro<- la_retro[order(la_retro$model,la_retro$endyr),]

sim=subset(resparam,scenario=='regimeProd'&parameter=='logalpha')
la_simregprod=data.frame(sim=sim$sim[match(seq(min(sim$by),max(sim$by)),sim$by)],by=seq(min(sim$by),max(sim$by))-55)

la_regprod=la_retro[la_retro$scenario=='regimeProd',]
la_regprod$plotMod=dplyr::recode_factor(factor(la_regprod$model),
                                              "autocorr" = "Static",
                                              "rwa"="Random walk",
                                              "hmma"="Hidden Markov model"
                                              )

gregprod<-ggplot(la_regprod) + 
  geom_line(aes(x=by-50,y= median.est,color=factor(endyr),group=factor(endyr)),linewidth=1.2)+
  geom_line(data=la_simregprod,aes(x=by,y=sim),linewidth=1.5)+
  facet_grid(.~plotMod)+
  scale_color_viridis_d('years of data') +
  mytheme+ 
  ylab("log(alpha) estimate") +
  xlab("year")

gregprod  

sim=subset(resparam,scenario=='shiftProd'&parameter=='logalpha')
la_simshiftprod=data.frame(sim=sim$sim[match(seq(min(sim$by),max(sim$by)),sim$by)],by=seq(min(sim$by),max(sim$by))-45)

la_shiftprod=la_retro[la_retro$scenario=='shiftProd',]
la_shiftprod$plotMod=dplyr::recode_factor(factor(la_regprod$model),
                                        "autocorr" = "Static",
                                        "rwa"="Random walk",
                                        "hmma"="Hidden Markov model"
)

gshiftprod<-ggplot(la_shiftprod) + 
  geom_line(aes(x=by-50,y= median.est,color=factor(endyr),group=factor(endyr)),linewidth=1.2)+
  geom_line(data=la_simshiftprod,aes(x=by,y=sim),linewidth=1.5)+
  facet_grid(.~plotMod)+
  scale_color_viridis_d('years of data') +
  mytheme+ 
  ylab("log(alpha) estimate") +
  xlab("year")

gshiftprod

legend <- cowplot::get_legend(
  # create some space to the left of the legend
  gregprod
)
comb<-cowplot::plot_grid(
  gregprod+ theme(legend.position="none"), gshiftprod+ theme(legend.position="none"),
  align = "hv", axis = "bt",
  rel_heights = c(.5,.5),
  nrow=2
)
comb2<-cowplot::plot_grid(comb, legend,
  align = "h", axis = "bt",
  rel_widths = c(.8,.2)
)
comb2
ggsave("outs/figures/retrospective_logalpha.png",plot=comb2,width=14,height=6)


##Guidance doc  version###
sms_retro<-res2[res2$parameter=="smsy"&res2$model %in%c('rwa'),]
sms_retro<- subset(sms_retro,endyr %in% seq(min(sms_retro$endyr),max(sms_retro$endyr),by=2))
sms_retro<- sms_retro[order(sms_retro$model,sms_retro$endyr),]

sim1=subset(resparam,scenario=='regimeProd'&parameter=='smsy')
sms_simregprod=data.frame(sim=c(sim1$sim[match(seq(min(sim1$by),max(sim1$by)),sim1$by)],rep(sim1$sim[1],5)),by=seq(min(sim1$by),max(sim1$by+5))-55)

sim2=subset(resparam,scenario=='shiftProd'&parameter=='smsy')
sms_simshiftprod=data.frame(sim=c(sim2$sim[match(seq(min(sim2$by),max(sim2$by)),sim2$by)],rep(sim2$sim[1],5)),by=seq(min(sim2$by),max(sim2$by+5))-55)

sms_gd1=sms_retro[sms_retro$scenario%in%c('regimeProd'),]

gregprod<-ggplot(sms_gd1) + 
  geom_line(aes(x=by-50,y= median.est,color=factor(endyr),group=factor(endyr)),linewidth=1.2)+
  geom_line(data=sms_simregprod,aes(x=by,y=sim),linewidth=1.5)+
  scale_color_viridis_d('years of data') +
  mytheme+ 
  ylab("Smsy estimate") +
  xlab("year")

gregprod

sms_gd2=sms_retro[sms_retro$scenario%in%c('shiftProd'),]

gshiftprod<-ggplot(sms_gd2) + 
  geom_line(aes(x=by-50,y= median.est,color=factor(endyr),group=factor(endyr)),linewidth=1.2)+
  geom_line(data=sms_simshiftprod,aes(x=by,y=sim),linewidth=1.5)+
  scale_color_viridis_d('years of data') +
  mytheme+ 
  ylab("Smsy estimate") +
  xlab("year")

gshiftprod

legend <- cowplot::get_legend(
  # create some space to the left of the legend
  gregprod
)

comb<-cowplot::plot_grid(
  gregprod+ theme(legend.position="none"), gshiftprod+ theme(legend.position="none"),
  align = "hv", axis = "bt",
  rel_heights = c(.5,.5),
  nrow=2
)
comb2<-cowplot::plot_grid(comb, legend,
                          align = "h", axis = "bt",
                          rel_widths = c(.8,.2)
)
comb2
ggsave("outs/figures/retrospective_smsy_gd.png",plot=comb2,width=8,height=6)


#SMSY####
smsy_retro<-res2[res2$parameter=="smsy"&res2$model %in%c('autocorr','hmma','rwa'),]
smsy_retro<- subset(smsy_retro,endyr %in% seq(min(smsy_retro$endyr),max(smsy_retro$endyr),by=2))
smsy_retro<- smsy_retro[order(smsy_retro$model,smsy_retro$endyr),]

sim=subset(resparam,scenario=='regimeProd'&parameter=='smsy')
smsy_simregprod=data.frame(sim=sim$sim[match(seq(min(sim$by),max(sim$by)),sim$by)],by=seq(min(sim$by),max(sim$by))-50)

smsy_regprod=smsy_retro[smsy_retro$scenario=='regimeProd',]
smsy_regprod$plotMod=dplyr::recode_factor(factor(smsy_regprod$model),
                                        "autocorr" = "Static",
                                        "rwa"="Random walk",
                                        "hmma"="Hidden Markov model"
)

gregprod<-ggplot(smsy_regprod) + 
  geom_line(aes(x=by-50,y= median.est,color=factor(endyr),group=factor(endyr)),linewidth=1.2)+
  geom_line(data=smsy_simregprod,aes(x=by,y=sim),linewidth=1.5)+
  facet_grid(.~plotMod)+
  scale_color_viridis_d('years of data') +
  mytheme+ 
  ylab("Smsy estimate") +
  xlab("year")

gregprod  

sim=subset(resparam,scenario=='shiftProd'&parameter=='smsy')
smsy_simshiftprod=data.frame(sim=sim$sim[match(seq(min(sim$by),max(sim$by)),sim$by)],by=seq(min(sim$by),max(sim$by))-50)

smsy_shiftprod=smsy_retro[smsy_retro$scenario=='shiftProd',]
smsy_shiftprod$plotMod=dplyr::recode_factor(factor(smsy_shiftprod$model),
                                          "autocorr" = "Static",
                                          "rwa"="Random walk",
                                          "hmma"="Hidden Markov model"
)

gshiftprod<-ggplot(smsy_shiftprod) + 
  geom_line(aes(x=by-50,y= median.est,color=factor(endyr),group=factor(endyr)),linewidth=1.2)+
  geom_line(data=smsy_simshiftprod,aes(x=by,y=sim),linewidth=1.5)+
  facet_grid(.~plotMod)+
  scale_color_viridis_d('Years of data') +
  mytheme+ 
  ylab("Smsy estimate") +
  xlab("year")

gshiftprod

legend <- cowplot::get_legend(
  # create some space to the left of the legend
  gregprod
)
comb<-cowplot::plot_grid(
  gregprod+ theme(legend.position="none"), gshiftprod+ theme(legend.position="none"),
  align = "hv", axis = "bt",
  rel_heights = c(.5,.5),
  nrow=2
)
comb2<-cowplot::plot_grid(comb, legend,
                          align = "h", axis = "bt",
                          rel_widths = c(.85,.15)
)
comb2
ggsave("outs/figures/retrospective_smsy.png",plot=comb2,width=14,height=6)

#SMAX####
smax_retro<-res2[res2$parameter=="smax"&res2$model %in%c('autocorr','hmmb','rwb'),]
smax_retro<- subset(smax_retro,endyr %in% seq(min(smax_retro$endyr),max(smax_retro$endyr),by=2))
smax_retro<- smax_retro[order(smax_retro$model,smax_retro$endyr),]

sim=subset(resparam,scenario=='regimeCap'&parameter=='smax')
smax_simregcap=data.frame(sim=sim$sim[match(seq(min(sim$by),max(sim$by)),sim$by)],by=seq(min(sim$by),max(sim$by))-50)

smax_regcap=smax_retro[smax_retro$scenario=='regimeProd',]
smax_regcap$plotMod=dplyr::recode_factor(factor(smax_regcap$model),
                                          "autocorr" = "Static",
                                          "rwb"="Random walk",
                                          "hmmb"="Hidden Markov model"
)

gregcap<-ggplot(smax_regcap) + 
  geom_line(aes(x=by-50,y= median.est,color=factor(endyr),group=factor(endyr)),linewidth=1.2)+
  geom_line(data=smax_simregcap,aes(x=by,y=sim),linewidth=1.5)+
  facet_grid(.~plotMod)+
  scale_color_viridis_d('years of data') +
  mytheme+ 
  ylab("Smax estimate") +
  xlab("year")

gregcap  

sim=subset(resparam,scenario=='shiftCap'&parameter=='smax')
smax_simshiftcap=data.frame(sim=sim$sim[match(seq(min(sim$by),max(sim$by)),sim$by)],by=seq(min(sim$by),max(sim$by))-50)

smax_shiftcap=smax_retro[smax_retro$scenario=='shiftCap',]
smax_shiftcap$plotMod=dplyr::recode_factor(factor(smax_shiftprod$model),
                                            "autocorr" = "Static",
                                            "rwb"="Random walk",
                                            "hmmb"="Hidden Markov model"
)

gshiftcap<-ggplot(smax_shiftcap) + 
  geom_line(aes(x=by-50,y= median.est,color=factor(endyr),group=factor(endyr)),linewidth=1.2)+
  geom_line(data=smax_simshiftcap,aes(x=by,y=sim),linewidth=1.5)+
  facet_grid(.~plotMod)+
  scale_color_viridis_d('Years of data') +
  mytheme+ 
  ylab("Smax estimate") +
  xlab("year")

gshiftcap

legend <- cowplot::get_legend(
  # create some space to the left of the legend
  gregcap
)
comb<-cowplot::plot_grid(
  gregcap+ theme(legend.position="none"), gshiftcap+ theme(legend.position="none"),
  align = "hv", axis = "bt",
  rel_heights = c(.5,.5),
  nrow=2
)
comb2<-cowplot::plot_grid(comb, legend,
                          align = "h", axis = "bt",
                          rel_widths = c(.85,.15)
)
comb2
ggsave("outs/figures/retrospective_smax.png",plot=comb2,width=14,height=6)


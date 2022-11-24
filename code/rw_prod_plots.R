#plot scenarios
library(here)
simPar <- read.csv("data/HarCk/harcnkSimPars.csv")

#Random walk prod
simData=readRDS(paste0("outs/SamSimOutputs/simData/", simPar$nameOM[5],"/",simPar$scenario[5],"/",
               paste(simPar$nameOM[5],"_", simPar$nameMP[5], "_", "CUsrDat.RData",sep="")))$srDatout


#read in estimation 
load("outs/simest/simest_prodcapscenarios5_11.Rdata")
dfpbias<-do.call("rbind",allsimest[[5]])


plot(dat$alpha,ylim=c(min(simData$alpha),max(simData$alpha)),type='n',bty='l',ylab='Ricker alpha',main='Simulated productivity trajectories')
a_slope=numeric(100)

for(u in 1:100){
dat<-simData[simData$iteration==u,]
dat<-dat[dat$year>(max(dat$year)-46),]

lines(dat$alpha,col=adjustcolor('darkgray',alpha=0.3))
a_slope[u]=lm(dat$alpha~dat$year)$coefficients[2]
}

hist(a_slope,breaks=20,main='Slope of alpha over time');abline(v=0,lwd=2,col='red')

#Linear decline scenario
simData=readRDS(paste0("outs/SamSimOutputs/simData/", simPar$nameOM[2],"/",simPar$scenario[2],"/",
                       paste(simPar$nameOM[2],"_", simPar$nameMP[2], "_", "CUsrDat.RData",sep="")))$srDatout

plot(dat$alpha,ylim=c(min(simData$alpha),max(simData$alpha)),type='n',bty='l',ylab='Ricker alpha',main='Simulated productivity trajectories')
a_slope=numeric(100)
for(u in 1:100){
  dat<-simData[simData$iteration==u,]
  dat<-dat[dat$year>(max(dat$year)-46),]
  
  lines(dat$alpha,col=adjustcolor('darkgray',alpha=0.3))
  a_slope[u]=lm(dat$alpha~dat$year)$coefficients[2]
  
  plot(dat$obsRecruits~dat$spawners)
}

hist(a_slope,breaks=20,main='Slope of alpha over time');abline(v=0,lwd=2,col='red')

dfpbias$scenario<-simPar$scenario[5]
dfpbias$slope<-a_slope[dfpbias$iteration]
dfpbias$trend<-"positive"
dfpbias$trend[dfpbias$slope<0]<-"negative"

dff<-dfpbias[dfpbias$parameter=="alpha",]

ggplot(dff,aes(x=model,y=pbias)) +
geom_boxplot(aes(fill=trend)) +
coord_cartesian(ylim = c(-100,100))+
geom_hline(yintercept=0) +
theme_bw(14)+ #theme(legend.position="none")+
facet_wrap(~method, scales="free_y")+
scale_fill_viridis_d(begin=.3, end=.9) +
ggtitle("alpha in rwa scenario")+
#stat_summary(fun.data = give.n, geom = "text", hjust = 0.5,
#    vjust = -2)+
theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
fig
#routines to visualize smulation estimation Results
#Catarina Wor
#November 2022
#============================================



#=================================================================================================================
#sensitivity a
simPar <- read.csv("data/sensitivity/SimPars.csv")

## Store relevant object names to help run simulation 
scenNames <- unique(simPar$scenario)

restmb<-readRDS(file = "outs/simest/sensitivity/res_a.rds")


#restmb<-rbind(resa1,resa2,resa3,resa4)

resstana1<-readRDS(file = "outs/simest/sensitivity/resstan_aq1.rds")
resstana2<-readRDS(file = "outs/simest/sensitivity/resstan_aq2.rds")
resstana3<-readRDS(file = "outs/simest/sensitivity/resstan_aq3.rds")
resstana4<-readRDS(file = "outs/simest/sensitivity/resstan_aq4.rds")

resstan<-rbind(resstana1,resstana2,resstana3,resstana4)

res<-rbind(restmb,resstan)

res$parameter[res$parameter=="Smax"]<-"smax"
res$parameter[res$parameter=="alpha"]<-"logalpha"

resparam<-res[res$parameter%in%c("logalpha","smax","sigma","smsy","sgen","umsy"),]

#exclude outliers
resparam$convergence[resparam$parameter=="logalpha"&resparam$mode>40]<-1
resparam$convergence[resparam$parameter=="smax"&resparam$mode>1e8]<-1



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


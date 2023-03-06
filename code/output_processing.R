library(here);library(ggplot2)
#Output processing:
sc14=readRDS(here('outs','cluster','sc1_4.RDS'))
sc58=readRDS(here('outs','cluster','sc5_8.RDS'))
sc912=readRDS(here('outs','cluster','sc9_12.RDS'))


#set 1: continuous change
out=list()
out[[1]]=subset(sc14,scenario=='stationary')
out[[2]]=subset(sc14,scenario=='autocorr')
out[[3]]=subset(sc14,scenario=='decLinearProd')
out[[4]]=subset(sc58,scenario=='decLinearCap')
out[[5]]=subset(sc912,scenario=='decLinearProdshiftCap')

#set 2: rapid change
out[[6]]=subset(sc912,scenario=='shiftProd')
out[[7]]=subset(sc912,scenario=='shiftCap')
out[[8]]=subset(sc912,scenario=='regimeProdCap')

#summarize model selection
aic_set=list()
aic_set2=list()
bic_set=list()
bic_set2=list()
aic_d90_set=list()
bic_d90_set=list()
aic_d80_set=list()
bic_d80_set=list()
aic_npar2_set=list()
bic_npar2_set=list()
aic_npar2_d90_set=list()
bic_npar2_d90_set=list()
aic_npar2_d80_set=list()
bic_npar2_d80_set=list()
for(i in 1:length(out)){
  aic=subset(out[[i]],parameter=='AIC')
  aic=subset(aic,model %in% c('simple','autocorr','rwa','rwb','rwab'))
  aic_set[[i]]=tidyr::spread(aic[,-9],key=model,value=est)
  aic_set[[i]]=aic_set[[i]][c(12,8,9,10,11)]
  aic2=subset(out[[i]],parameter=='AIC_n')
  aic2=subset(aic2,model %in% c('simple','autocorr','rwa','rwb','rwab'))
  aic_set2[[i]]=tidyr::spread(aic2[,-9],key=model,value=est)
  aic_set2[[i]]=aic_set2[[i]][c(12,8,9,10,11)]
  
  bic=subset(out[[i]],parameter=='BIC')
  bic_set[[i]]=tidyr::spread(bic[,-9],key=model,value=est)
  bic_set[[i]]=bic_set[[i]][c(12,8,9,10,11)]
  bic2=subset(out[[i]],parameter=='bic_n')
  bic_set2[[i]]=tidyr::spread(bic2[,-9],key=model,value=est)
  bic_set2[[i]]=bic_set2[[i]][c(12,8,9,10,11)]
  
  aic_d90=subset(out[[i]],parameter=='AIC_d90')
  aic_d90=subset(aic_d90,model %in% c('simple','autocorr','rwa','rwb','rwab'))
  aic_d90_set[[i]]=tidyr::spread(aic_d90[,-9],key=model,value=est)
  aic_d90_set[[i]]=aic_d90_set[[i]][c(12,8,9,10,11)]
  bic_d90=subset(out[[i]],parameter=='bic_d90')
  bic_d90_set[[i]]=tidyr::spread(bic_d90[,-9],key=model,value=est)
  bic_d90_set[[i]]=bic_d90_set[[i]][c(12,8,9,10,11)]
  
  aic_d80=subset(out[[i]],parameter=='AIC_d80')
  aic_d80=subset(aic_d80,model %in% c('simple','autocorr','rwa','rwb','rwab'))
  aic_d80_set[[i]]=tidyr::spread(aic_d80[,-9],key=model,value=est)
  aic_d80_set[[i]]=aic_d80_set[[i]][c(12,8,9,10,11)]
  bic_d80=subset(out[[i]],parameter=='bic_d80')
  bic_d80_set[[i]]=tidyr::spread(bic_d80[,-9],key=model,value=est)
  bic_d80_set[[i]]=bic_d80_set[[i]][c(12,8,9,10,11)]
  
  aic_npar2=subset(out[[i]],parameter=='AIC_npar2')
  aic_npar2=subset(aic_npar2,model %in% c('simple','autocorr','rwa','rwb','rwab'))
  aic_npar2_set[[i]]=tidyr::spread(aic_npar2[,-9],key=model,value=est)
  aic_npar2_set[[i]]=aic_npar2_set[[i]][c(12,8,9,10,11)]
  bic_npar2=subset(out[[i]],parameter=='bic_npar2')
  bic_npar2_set[[i]]=tidyr::spread(bic_npar2[,-9],key=model,value=est)
  bic_npar2_set[[i]]=bic_npar2_set[[i]][c(12,8,9,10,11)]
  
  aic_npar2_d90=subset(out[[i]],parameter=='AIC_d90npar2')
  aic_npar2_d90=subset(aic_npar2_d90,model %in% c('simple','autocorr','rwa','rwb','rwab'))
  aic_npar2_d90_set[[i]]=tidyr::spread(aic_npar2_d90[,-9],key=model,value=est)
  aic_npar2_d90_set[[i]]=aic_npar2_d90_set[[i]][c(12,8,9,10,11)]
  bic_npar2_d90=subset(out[[i]],parameter=='bic_d90npar2')
  bic_npar2_d90_set[[i]]=tidyr::spread(bic_npar2_d90[,-9],key=model,value=est)
  bic_npar2_d90_set[[i]]=bic_npar2_d90_set[[i]][c(12,8,9,10,11)]
  
  aic_npar2_d80=subset(out[[i]],parameter=='AIC_d80npar2')
  aic_npar2_d80=subset(aic_npar2_d80,model %in% c('simple','autocorr','rwa','rwb','rwab'))
  aic_npar2_d80_set[[i]]=tidyr::spread(aic_npar2_d80[,-9],key=model,value=est)
  aic_npar2_d80_set[[i]]=aic_npar2_d80_set[[i]][c(12,8,9,10,11)]
  bic_npar2_d80=subset(out[[i]],parameter=='bic_d80npar2')
  bic_npar2_d80_set[[i]]=tidyr::spread(bic_npar2_d80[,-9],key=model,value=est)
  bic_npar2_d80_set[[i]]=bic_npar2_d80_set[[i]][c(12,8,9,10,11)]
}



#AIC 1####
#first set of scenarios
sc1=apply(aic_set2[[1]],1,which.min)
cn1=summary(factor(sc1,levels=seq(1:5)))/2000
w1=matrix(ncol=5,nrow=nrow(aic_set2[[1]]))
for(i in 1:nrow(aic_set2[[1]])){
w1[i,]=unlist(samEst::model_weights(aic_set2[[1]][i,],form='AIC'))
}
w_avg1=numeric(5)
for(i in 1:5){
  ws=w1[sc1==i,i]
  w_avg1[i]=mean(ws)
}

sc2=apply(aic_set2[[2]],1,which.min)
cn2=summary(factor(sc2,levels=seq(1:5)))/2000
w2=matrix(ncol=5,nrow=nrow(aic_set2[[2]]))
for(i in 1:nrow(aic_set2[[2]])){
  w2[i,]=unlist(samEst::model_weights(aic_set2[[2]][i,],form='AIC'))
}
w_avg2=numeric(5)
for(i in 1:5){
  ws=w2[sc2==i,i]
  w_avg2[i]=mean(ws)
}

sc3=apply(aic_set2[[3]],1,which.min)
cn3=summary(factor(sc3,levels=seq(1:5)))/2000
w3=matrix(ncol=5,nrow=nrow(aic_set2[[3]]))
for(i in 1:nrow(aic_set2[[3]])){
  w3[i,]=unlist(samEst::model_weights(aic_set2[[3]][i,],form='AIC'))
}
w_avg3=numeric(5)
for(i in 1:5){
  ws=w3[sc3==i,i]
  w_avg3[i]=mean(ws)
}
sc4=apply(aic_set2[[4]],1,which.min)
cn4=summary(factor(sc4,levels=seq(1:5)))/2000
w4=matrix(ncol=5,nrow=nrow(aic_set2[[4]]))
for(i in 1:nrow(aic_set2[[4]])){
  w4[i,]=unlist(samEst::model_weights(aic_set2[[4]][i,],form='AIC'))
}
w_avg4=numeric(5)
for(i in 1:5){
  ws=w4[sc4==i,i]
  w_avg4[i]=mean(ws)
}
sc5=apply(aic_set2[[5]],1,which.min)
cn5=summary(factor(sc5,levels=seq(1:5)))/2000
w5=matrix(ncol=5,nrow=nrow(aic_set2[[5]]))
for(i in 1:nrow(aic_set2[[5]])){
  w5[i,]=unlist(samEst::model_weights(aic_set2[[5]][i,],form='AIC'))
}
w_avg5=numeric(5)
for(i in 1:5){
  ws=w5[sc5==i,i]
  w_avg5[i]=mean(ws)
}
#second set of scenarios
sc6=apply(aic_set2[[6]],1,which.min)
cn6=summary(factor(sc6,levels=seq(1:5)))/2000
w6=matrix(ncol=5,nrow=nrow(aic_set2[[6]]))
for(i in 1:nrow(aic_set2[[6]])){
  w6[i,]=unlist(samEst::model_weights(aic_set2[[6]][i,],form='AIC'))
}
w_avg6=numeric(5)
for(i in 1:5){
  ws=w6[sc6==i,i]
  w_avg6[i]=mean(ws)
}
sc7=apply(aic_set2[[7]],1,which.min)
cn7=summary(factor(sc7,levels=seq(1:5)))/2000
w7=matrix(ncol=5,nrow=nrow(aic_set2[[7]]))
for(i in 1:nrow(aic_set2[[7]])){
  w7[i,]=unlist(samEst::model_weights(aic_set2[[7]][i,],form='AIC'))
}
w_avg7=numeric(5)
for(i in 1:5){
  ws=w7[sc7==i,i]
  w_avg7[i]=mean(ws)
}
sc8=apply(aic_set2[[8]],1,which.min)
cn8=summary(factor(sc8,levels=seq(1:5)))/2000
w8=matrix(ncol=5,nrow=nrow(aic_set2[[8]]))
for(i in 1:nrow(aic_set2[[8]])){
  w8[i,]=unlist(samEst::model_weights(aic_set2[[8]][i,],form='AIC'))
}
w_avg8=numeric(5)
for(i in 1:5){
  ws=w8[sc8==i,i]
  w_avg8[i]=mean(ws)
}

##Confusion matrices
conf_matrix<-expand.grid(EM=c("stationary",
                               "autocorr",
                               "dynamic.a","dynamic.b","dynamic.ab"
                               ),OM=c("stationary",
                                                                        "autocorr",
                                                                        "dec.prod","dec.cap","dec.prodcap"))

conf_matrix$w_AIC=NA
conf_matrix$w_AIC[1:5]=cn1
conf_matrix$w_AIC[6:10]=cn2
conf_matrix$w_AIC[11:15]=cn3
conf_matrix$w_AIC[16:20]=cn4
conf_matrix$w_AIC[21:25]=cn5
conf_matrix$top_w[1:5]=w_avg1
conf_matrix$top_w[6:10]=w_avg2
conf_matrix$top_w[11:15]=w_avg3
conf_matrix$top_w[16:20]=w_avg4
conf_matrix$top_w[21:25]=w_avg5

mytheme = list(
  theme_classic(14)+
    theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
          legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="bold", size=13),
          axis.text=element_text(face="bold"),axis.title = element_text(face="bold",size=16),plot.title = element_text(face = "bold", hjust = 0.5,size=16))
)

library(ggplot2)



p=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = w_AIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_AIC,2)), vjust = 1,size=6) +
  ggtitle("AIC")+
  scale_fill_gradient(low = "white", high = "navy") +
  mytheme + theme(legend.position="none")+xlab("Simulation Scenario")+ylab("Estimation Model")
p

p2=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = top_w), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_AIC,2)), vjust = 1,size=6) +
  ggtitle("AIC")+
  scale_fill_gradient(low = "white", high = "skyblue") +
  mytheme + theme(legend.position="right")+xlab("Simulation Scenario")+ylab("Estimation Model")
p2


conf_matrix2<-expand.grid(EM=c("stationary",
                              "autocorr",
                              "dynamic.a","dynamic.b","dynamic.ab"
),OM=c("stationary",
       "autocorr",
       "regime.prod","regime.cap","regime.prodcap"))

conf_matrix2$w_AIC=NA
conf_matrix2$w_AIC[1:5]=cn1
conf_matrix2$w_AIC[6:10]=cn2
conf_matrix2$w_AIC[11:15]=cn6
conf_matrix2$w_AIC[16:20]=cn7
conf_matrix2$w_AIC[21:25]=cn8
conf_matrix2$top_w[1:5]=w_avg1
conf_matrix2$top_w[6:10]=w_avg2
conf_matrix2$top_w[11:15]=w_avg6
conf_matrix2$top_w[16:20]=w_avg7
conf_matrix2$top_w[21:25]=w_avg8

p=ggplot(data =  conf_matrix2, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = w_AIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_AIC,2)), vjust = 1,size=6) +
  ggtitle("AIC")+
  scale_fill_gradient(low = "white", high = "navy") +
  mytheme + theme(legend.position="none")+xlab("Simulation Scenario")+ylab("Estimation Model")
p

p2=ggplot(data =  conf_matrix2, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = top_w), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_AIC,2)), vjust = 1,size=6) +
  ggtitle("AIC")+
  scale_fill_gradient(low = "white", high = "skyblue") +
  mytheme + theme(legend.position="right")+xlab("Simulation Scenario")+ylab("Estimation Model")
p2


#AIC 2####
#first set of scenarios
sc1=apply(aic_set2[[1]],1,which.min)
cn1=summary(factor(sc1,levels=seq(1:5)))/2000
w1=matrix(ncol=5,nrow=nrow(aic_set2[[1]]))
for(i in 1:nrow(aic_set2[[1]])){
  w1[i,]=unlist(samEst::model_weights(aic_set2[[1]][i,],form='AIC'))
}
w_avg1=numeric(5)
for(i in 1:5){
  ws=w1[sc1==i,i]
  w_avg1[i]=mean(ws)
}

sc2=apply(aic_set2[[2]],1,which.min)
cn2=summary(factor(sc2,levels=seq(1:5)))/2000
w2=matrix(ncol=5,nrow=nrow(aic_set2[[2]]))
for(i in 1:nrow(aic_set2[[2]])){
  w2[i,]=unlist(samEst::model_weights(aic_set2[[2]][i,],form='AIC'))
}
w_avg2=numeric(5)
for(i in 1:5){
  ws=w2[sc2==i,i]
  w_avg2[i]=mean(ws)
}

sc3=apply(aic_set2[[3]],1,which.min)
cn3=summary(factor(sc3,levels=seq(1:5)))/2000
w3=matrix(ncol=5,nrow=nrow(aic_set2[[3]]))
for(i in 1:nrow(aic_set2[[3]])){
  w3[i,]=unlist(samEst::model_weights(aic_set2[[3]][i,],form='AIC'))
}
w_avg3=numeric(5)
for(i in 1:5){
  ws=w3[sc3==i,i]
  w_avg3[i]=mean(ws)
}
sc4=apply(aic_set2[[4]],1,which.min)
cn4=summary(factor(sc4,levels=seq(1:5)))/2000
w4=matrix(ncol=5,nrow=nrow(aic_set2[[4]]))
for(i in 1:nrow(aic_set2[[4]])){
  w4[i,]=unlist(samEst::model_weights(aic_set2[[4]][i,],form='AIC'))
}
w_avg4=numeric(5)
for(i in 1:5){
  ws=w4[sc4==i,i]
  w_avg4[i]=mean(ws)
}
sc5=apply(aic_set2[[5]],1,which.min)
cn5=summary(factor(sc5,levels=seq(1:5)))/2000
w5=matrix(ncol=5,nrow=nrow(aic_set2[[5]]))
for(i in 1:nrow(aic_set2[[5]])){
  w5[i,]=unlist(samEst::model_weights(aic_set2[[5]][i,],form='AIC'))
}
w_avg5=numeric(5)
for(i in 1:5){
  ws=w5[sc5==i,i]
  w_avg5[i]=mean(ws)
}
#second set of scenarios
sc6=apply(aic_set2[[6]],1,which.min)
cn6=summary(factor(sc6,levels=seq(1:5)))/2000
w6=matrix(ncol=5,nrow=nrow(aic_set2[[6]]))
for(i in 1:nrow(aic_set2[[6]])){
  w6[i,]=unlist(samEst::model_weights(aic_set2[[6]][i,],form='AIC'))
}
w_avg6=numeric(5)
for(i in 1:5){
  ws=w6[sc6==i,i]
  w_avg6[i]=mean(ws)
}
sc7=apply(aic_set2[[7]],1,which.min)
cn7=summary(factor(sc7,levels=seq(1:5)))/2000
w7=matrix(ncol=5,nrow=nrow(aic_set2[[7]]))
for(i in 1:nrow(aic_set2[[7]])){
  w7[i,]=unlist(samEst::model_weights(aic_set2[[7]][i,],form='AIC'))
}
w_avg7=numeric(5)
for(i in 1:5){
  ws=w7[sc7==i,i]
  w_avg7[i]=mean(ws)
}
sc8=apply(aic_set2[[8]],1,which.min)
cn8=summary(factor(sc8,levels=seq(1:5)))/2000
w8=matrix(ncol=5,nrow=nrow(aic_set2[[8]]))
for(i in 1:nrow(aic_set2[[8]])){
  w8[i,]=unlist(samEst::model_weights(aic_set2[[8]][i,],form='AIC'))
}
w_avg8=numeric(5)
for(i in 1:5){
  ws=w8[sc8==i,i]
  w_avg8[i]=mean(ws)
}

##Confusion matrices
conf_matrix<-expand.grid(EM=c("stationary",
                              "autocorr",
                              "dynamic.a","dynamic.b","dynamic.ab"
),OM=c("stationary",
       "autocorr",
       "dec.prod","dec.cap","dec.prodcap","regime.prod","regime.cap","regime.prodcap"))

conf_matrix$w_AIC=NA
conf_matrix$w_AIC[1:5]=cn1
conf_matrix$w_AIC[6:10]=cn2
conf_matrix$w_AIC[11:15]=cn3
conf_matrix$w_AIC[16:20]=cn4
conf_matrix$w_AIC[21:25]=cn5
conf_matrix$w_AIC[26:30]=cn6
conf_matrix$w_AIC[31:35]=cn7
conf_matrix$w_AIC[36:40]=cn8

conf_matrix$top_w[1:5]=w_avg1
conf_matrix$top_w[6:10]=w_avg2
conf_matrix$top_w[11:15]=w_avg3
conf_matrix$top_w[16:20]=w_avg4
conf_matrix$top_w[21:25]=w_avg5
conf_matrix$top_w[26:30]=w_avg6
conf_matrix$top_w[31:35]=w_avg7
conf_matrix$top_w[36:40]=w_avg8

library(ggplot2)

mytheme = list(
  theme_classic(14)+
    theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
          legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="bold", size=13),
          axis.text=element_text(face="bold"),axis.title = element_text(face="bold",size=16),plot.title = element_text(face = "bold", hjust = 0.5,size=16))
)



p=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = w_AIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_AIC,2)), vjust = 1,size=6) +
  ggtitle("AIC")+
  scale_fill_gradient(low = "white", high = "navy") +
  mytheme + theme(legend.position="none")+xlab("Simulation Scenario")+ylab("Estimation Model")
p
pdf(file='AIC_conf_matrix.pdf',height=8,width=12)
p
dev.off()

p2=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = top_w), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_AIC,2)), vjust = 1,size=6) +
  ggtitle("AIC")+
  scale_fill_gradient(low = "white", high = "skyblue") +
  mytheme + theme(legend.position="right")+xlab("Simulation Scenario")+ylab("")
p2

pdf(file='AIC_conf_matrix_topw.pdf',height=8,width=12)
p2
dev.off()


conf_matrix2<-expand.grid(EM=c("stationary",
                               "autocorr",
                               "dynamic.a","dynamic.b","dynamic.ab"
),OM=c("stationary",
       "autocorr",
       "regime.prod","regime.cap","regime.prodcap"))

conf_matrix2$w_AIC=NA
conf_matrix2$w_AIC[1:5]=cn1
conf_matrix2$w_AIC[6:10]=cn2
conf_matrix2$w_AIC[11:15]=cn6
conf_matrix2$w_AIC[16:20]=cn7
conf_matrix2$w_AIC[21:25]=cn8
conf_matrix2$top_w[1:5]=w_avg1
conf_matrix2$top_w[6:10]=w_avg2
conf_matrix2$top_w[11:15]=w_avg6
conf_matrix2$top_w[16:20]=w_avg7
conf_matrix2$top_w[21:25]=w_avg8

p=ggplot(data =  conf_matrix2, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = w_AIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_AIC,2)), vjust = 1,size=6) +
  ggtitle("AIC")+
  scale_fill_gradient(low = "white", high = "navy") +
  mytheme + theme(legend.position="none")+xlab("Simulation Scenario")+ylab("Estimation Model")


pdf(file='AIC_conf_matrix_scn5_8.pdf',height=8,width=8)
p
dev.off()


p2=ggplot(data =  conf_matrix2, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = top_w), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_AIC,2)), vjust = 1,size=6) +
  ggtitle("AIC")+
  scale_fill_gradient(low = "white", high = "skyblue") +
  mytheme + theme(axis.text.y=element_blank(),legend.position="right")+xlab("Simulation Scenario")+ylab("")
pdf(file='AIC_conf_matrix_scn5_8_topw.pdf',height=8,width=8)
p2
dev.off()


#By model class - AIC####
sc1m=ifelse(sc1<3,1,sc1)
sc1m=ifelse(sc1>2,2,sc1m)
cn1m=summary(factor(sc1m),levels=seq(1:2))/2000
sc2m=ifelse(sc1<3,1,sc2)
sc2m=ifelse(sc1>2,2,sc2m)
cn2m=summary(factor(sc2m),levels=seq(1:2))/2000
sc3m=ifelse(sc3<3,1,sc3)
sc3m=ifelse(sc3>2,2,sc3m)
cn3m=summary(factor(sc3m),levels=seq(1:2))/2000
sc4m=ifelse(sc4<3,1,sc4)
sc4m=ifelse(sc4>2,2,sc4m)
cn4m=summary(factor(sc4m),levels=seq(1:2))/2000
sc5m=ifelse(sc5<3,1,sc5)
sc5m=ifelse(sc5>2,2,sc5m)
cn5m=summary(factor(sc5m),levels=seq(1:2))/2000
sc6m=ifelse(sc6<3,1,sc6)
sc6m=ifelse(sc6>2,2,sc6m)
cn6m=summary(factor(sc6m),levels=seq(1:2))/2000
sc7m=ifelse(sc7<3,1,sc7)
sc7m=ifelse(sc7>2,2,sc7m)
cn7m=summary(factor(sc7m),levels=seq(1:2))/2000
sc8m=ifelse(sc8<3,1,sc8)
sc8m=ifelse(sc8>2,2,sc8m)
cn8m=summary(factor(sc8m),levels=seq(1:2))/2000


##Confusion matrices
conf_matrix <-expand.grid(EM=c('static','dynamic'),OM=c("stationary",
                                                        "autocorr",
                                                        "dec.prod","dec.cap","dec.prodcap","regime.prod","regime.cap","regime.prodcap"))

conf_matrix$w_AIC=NA
conf_matrix$w_AIC[1:2]=cn1m
conf_matrix$w_AIC[3:4]=cn2m
conf_matrix$w_AIC[5:6]=cn3m
conf_matrix$w_AIC[7:8]=cn4m
conf_matrix$w_AIC[9:10]=cn5m
conf_matrix$w_AIC[11:12]=cn6m
conf_matrix$w_AIC[13:14]=cn7m
conf_matrix$w_AIC[15:16]=cn8m

library(ggplot2)
p=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = w_AIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_AIC,2)), vjust = 1) +
  ggtitle("AIC")+
  scale_fill_gradient(low = "white", high = "navy") +
  mytheme + theme(legend.position="none")+xlab("Operating Model")+ylab("Estimation Model")
pdf(file='model_class_AICconf.pdf',width=10,height=6)
p
dev.off()


#parameter bias####
alpha_pbias=subset(out[[1]],parameter=='alpha')
alpha_est_l=list()
alpha_est_l[[1]]=subset(alpha_pbias,iteration %in% which(sc1==1))
alpha_est_l[[2]]=subset(alpha_pbias,iteration %in% which(sc1==2))
alpha_est_l[[3]]=subset(alpha_pbias,iteration %in% which(sc1==3))

df1=subset(alpha_pbias,iteration==1&model=='simple')
pdf('pbias_sc1.pdf',width=8,height=6)
plot(exp(df1$sim)~c(df1$by-54),type='l',lwd=5,col='darkred',ylab='Maximum Productivity (Recruits/Spawner)',xlab='Year',ylim=c(0,15),main='Stationary Scenario')
for(i in 1:length(unique(alpha_est_l[[3]]$iteration))){
  dfx=subset(alpha_est_l[[3]],iteration==which(sc1==3)[i])
  dfx=subset(dfx,model=='rwa')
  lines(exp(dfx$est),col=adjustcolor('black',alpha.f=0.3))
}
dev.off()

alpha_pbias=subset(out[[2]],parameter=='alpha')
alpha_est_l=list()
alpha_est_l[[1]]=subset(alpha_pbias,iteration %in% which(sc1==1))
alpha_est_l[[2]]=subset(alpha_pbias,iteration %in% which(sc1==2))
alpha_est_l[[3]]=subset(alpha_pbias,iteration %in% which(sc1==3))

df1=subset(alpha_pbias,iteration==1&model=='simple')
pdf('pbias_sc2.pdf',width=8,height=6)
plot(exp(df1$sim)~c(df1$by-54),type='l',lwd=5,col='darkred',ylab='Maximum Productivity (Recruits/Spawner)',xlab='Year',ylim=c(0,15),main='Autocorrelated Scenario')
for(i in 1:length(unique(alpha_est_l[[3]]$iteration))){
  dfx=subset(alpha_est_l[[3]],iteration==which(sc1==3)[i])
  dfx=subset(dfx,model=='rwa')
  lines(exp(dfx$est),col=adjustcolor('black',alpha.f=0.3))
}
dev.off()

#AIC d90####
#first set of scenarios
sc1=apply(aic_d90_set[[1]],1,which.min)
cn1=summary(factor(sc1,levels=seq(1:5)))/2000
w1=matrix(ncol=5,nrow=nrow(aic_d90_set[[1]]))
for(i in 1:nrow(aic_d90_set[[1]])){
  w1[i,]=unlist(samEst::model_weights(aic_d90_set[[1]][i,],form='AIC'))
}
w_avg1=numeric(5)
for(i in 1:5){
  ws=w1[sc1==i,i]
  w_avg1[i]=mean(ws)
}

sc2=apply(aic_d90_set[[2]],1,which.min)
cn2=summary(factor(sc2,levels=seq(1:5)))/2000
w2=matrix(ncol=5,nrow=nrow(aic_d90_set[[2]]))
for(i in 1:nrow(aic_d90_set[[2]])){
  w2[i,]=unlist(samEst::model_weights(aic_d90_set[[2]][i,],form='AIC'))
}
w_avg2=numeric(5)
for(i in 1:5){
  ws=w2[sc2==i,i]
  w_avg2[i]=mean(ws)
}

sc3=apply(aic_d90_set[[3]],1,which.min)
cn3=summary(factor(sc3,levels=seq(1:5)))/2000
w3=matrix(ncol=5,nrow=nrow(aic_d90_set[[3]]))
for(i in 1:nrow(aic_d90_set[[3]])){
  w3[i,]=unlist(samEst::model_weights(aic_d90_set[[3]][i,],form='AIC'))
}
w_avg3=numeric(5)
for(i in 1:5){
  ws=w3[sc3==i,i]
  w_avg3[i]=mean(ws)
}
sc4=apply(aic_d90_set[[4]],1,which.min)
cn4=summary(factor(sc4,levels=seq(1:5)))/2000
w4=matrix(ncol=5,nrow=nrow(aic_d90_set[[4]]))
for(i in 1:nrow(aic_d90_set[[4]])){
  w4[i,]=unlist(samEst::model_weights(aic_d90_set[[4]][i,],form='AIC'))
}
w_avg4=numeric(5)
for(i in 1:5){
  ws=w4[sc4==i,i]
  w_avg4[i]=mean(ws)
}
sc5=apply(aic_d90_set[[5]],1,which.min)
cn5=summary(factor(sc5,levels=seq(1:5)))/2000
w5=matrix(ncol=5,nrow=nrow(aic_d90_set[[5]]))
for(i in 1:nrow(aic_d90_set[[5]])){
  w5[i,]=unlist(samEst::model_weights(aic_d90_set[[5]][i,],form='AIC'))
}
w_avg5=numeric(5)
for(i in 1:5){
  ws=w5[sc5==i,i]
  w_avg5[i]=mean(ws)
}
#second set of scenarios
sc6=apply(aic_d90_set[[6]],1,which.min)
cn6=summary(factor(sc6,levels=seq(1:5)))/2000
w6=matrix(ncol=5,nrow=nrow(aic_d90_set[[6]]))
for(i in 1:nrow(aic_d90_set[[6]])){
  w6[i,]=unlist(samEst::model_weights(aic_d90_set[[6]][i,],form='AIC'))
}
w_avg6=numeric(5)
for(i in 1:5){
  ws=w6[sc6==i,i]
  w_avg6[i]=mean(ws)
}
sc7=apply(aic_d90_set[[7]],1,which.min)
cn7=summary(factor(sc7,levels=seq(1:5)))/2000
w7=matrix(ncol=5,nrow=nrow(aic_d90_set[[7]]))
for(i in 1:nrow(aic_d90_set[[7]])){
  w7[i,]=unlist(samEst::model_weights(aic_d90_set[[7]][i,],form='AIC'))
}
w_avg7=numeric(5)
for(i in 1:5){
  ws=w7[sc7==i,i]
  w_avg7[i]=mean(ws)
}
sc8=apply(aic_d90_set[[8]],1,which.min)
cn8=summary(factor(sc8,levels=seq(1:5)))/2000
w8=matrix(ncol=5,nrow=nrow(aic_d90_set[[8]]))
for(i in 1:nrow(aic_d90_set[[8]])){
  w8[i,]=unlist(samEst::model_weights(aic_d90_set[[8]][i,],form='AIC'))
}
w_avg8=numeric(5)
for(i in 1:5){
  ws=w8[sc8==i,i]
  w_avg8[i]=mean(ws)
}

##Confusion matrices
conf_matrix<-expand.grid(EM=c("stationary",
                              "autocorr",
                              "dynamic.a","dynamic.b","dynamic.ab"
),OM=c("stationary",
       "autocorr",
       "dec.prod","dec.cap","dec.prodcap",
       "regime.prod","regime.cap","regime.prodcap"))

conf_matrix$w_AIC=NA
conf_matrix$w_AIC[1:5]=cn1
conf_matrix$w_AIC[6:10]=cn2
conf_matrix$w_AIC[11:15]=cn3
conf_matrix$w_AIC[16:20]=cn4
conf_matrix$w_AIC[21:25]=cn5
conf_matrix$w_AIC[26:30]=cn6
conf_matrix$w_AIC[31:35]=cn7
conf_matrix$w_AIC[36:40]=cn8
conf_matrix$top_w[1:5]=w_avg1
conf_matrix$top_w[6:10]=w_avg2
conf_matrix$top_w[11:15]=w_avg3
conf_matrix$top_w[16:20]=w_avg4
conf_matrix$top_w[21:25]=w_avg5
conf_matrix$top_w[26:30]=w_avg6
conf_matrix$top_w[31:35]=w_avg7
conf_matrix$top_w[36:40]=w_avg8

mytheme = list(
  theme_classic(14)+
    theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
          legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="bold", size=13),
          axis.text=element_text(face="bold"),axis.title = element_text(face="bold",size=16),plot.title = element_text(face = "bold", hjust = 0.5,size=16))
)

library(ggplot2)



p=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = w_AIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_AIC,2)), vjust = 1,size=6) +
  ggtitle("AIC")+
  scale_fill_gradient(low = "white", high = "navy") +
  mytheme + theme(legend.position="none")+xlab("Simulation Scenario")+ylab("Estimation Model")
pdf(file='AIC_conf_matrix_d90.pdf',height=8,width=12)
p
dev.off()

p2=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = top_w), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_AIC,2)), vjust = 1,size=6) +
  ggtitle("AIC")+
  scale_fill_gradient(low = "white", high = "skyblue") +
  mytheme + theme(legend.position="right")+xlab("Simulation Scenario")+ylab("Estimation Model")
p2


conf_matrix2<-expand.grid(EM=c("stationary",
                               "autocorr",
                               "dynamic.a","dynamic.b","dynamic.ab"
),OM=c("stationary",
       "autocorr",
       "regime.prod","regime.cap","regime.prodcap"))

conf_matrix2$w_AIC=NA
conf_matrix2$w_AIC[1:5]=cn1
conf_matrix2$w_AIC[6:10]=cn2
conf_matrix2$w_AIC[11:15]=cn6
conf_matrix2$w_AIC[16:20]=cn7
conf_matrix2$w_AIC[21:25]=cn8
conf_matrix2$top_w[1:5]=w_avg1
conf_matrix2$top_w[6:10]=w_avg2
conf_matrix2$top_w[11:15]=w_avg6
conf_matrix2$top_w[16:20]=w_avg7
conf_matrix2$top_w[21:25]=w_avg8

p=ggplot(data =  conf_matrix2, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = w_AIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_AIC,2)), vjust = 1,size=6) +
  ggtitle("AIC")+
  scale_fill_gradient(low = "white", high = "navy") +
  mytheme + theme(legend.position="none")+xlab("Simulation Scenario")+ylab("Estimation Model")
p

p2=ggplot(data =  conf_matrix2, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = top_w), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_AIC,2)), vjust = 1,size=6) +
  ggtitle("AIC")+
  scale_fill_gradient(low = "white", high = "skyblue") +
  mytheme + theme(legend.position="right")+xlab("Simulation Scenario")+ylab("Estimation Model")
p2





#AIC d80####
#first set of scenarios
sc1=apply(aic_d80_set[[1]],1,which.min)
cn1=summary(factor(sc1,levels=seq(1:5)))/2000
w1=matrix(ncol=5,nrow=nrow(aic_d80_set[[1]]))
for(i in 1:nrow(aic_d80_set[[1]])){
  w1[i,]=unlist(samEst::model_weights(aic_d80_set[[1]][i,],form='AIC'))
}
w_avg1=numeric(5)
for(i in 1:5){
  ws=w1[sc1==i,i]
  w_avg1[i]=mean(ws)
}

sc2=apply(aic_d80_set[[2]],1,which.min)
cn2=summary(factor(sc2,levels=seq(1:5)))/2000
w2=matrix(ncol=5,nrow=nrow(aic_d80_set[[2]]))
for(i in 1:nrow(aic_d80_set[[2]])){
  w2[i,]=unlist(samEst::model_weights(aic_d80_set[[2]][i,],form='AIC'))
}
w_avg2=numeric(5)
for(i in 1:5){
  ws=w2[sc2==i,i]
  w_avg2[i]=mean(ws)
}

sc3=apply(aic_d80_set[[3]],1,which.min)
cn3=summary(factor(sc3,levels=seq(1:5)))/2000
w3=matrix(ncol=5,nrow=nrow(aic_d80_set[[3]]))
for(i in 1:nrow(aic_d80_set[[3]])){
  w3[i,]=unlist(samEst::model_weights(aic_d80_set[[3]][i,],form='AIC'))
}
w_avg3=numeric(5)
for(i in 1:5){
  ws=w3[sc3==i,i]
  w_avg3[i]=mean(ws)
}
sc4=apply(aic_d80_set[[4]],1,which.min)
cn4=summary(factor(sc4,levels=seq(1:5)))/2000
w4=matrix(ncol=5,nrow=nrow(aic_d80_set[[4]]))
for(i in 1:nrow(aic_d80_set[[4]])){
  w4[i,]=unlist(samEst::model_weights(aic_d80_set[[4]][i,],form='AIC'))
}
w_avg4=numeric(5)
for(i in 1:5){
  ws=w4[sc4==i,i]
  w_avg4[i]=mean(ws)
}
sc5=apply(aic_d80_set[[5]],1,which.min)
cn5=summary(factor(sc5,levels=seq(1:5)))/2000
w5=matrix(ncol=5,nrow=nrow(aic_d80_set[[5]]))
for(i in 1:nrow(aic_d80_set[[5]])){
  w5[i,]=unlist(samEst::model_weights(aic_d80_set[[5]][i,],form='AIC'))
}
w_avg5=numeric(5)
for(i in 1:5){
  ws=w5[sc5==i,i]
  w_avg5[i]=mean(ws)
}
#second set of scenarios
sc6=apply(aic_d80_set[[6]],1,which.min)
cn6=summary(factor(sc6,levels=seq(1:5)))/2000
w6=matrix(ncol=5,nrow=nrow(aic_d80_set[[6]]))
for(i in 1:nrow(aic_d80_set[[6]])){
  w6[i,]=unlist(samEst::model_weights(aic_d80_set[[6]][i,],form='AIC'))
}
w_avg6=numeric(5)
for(i in 1:5){
  ws=w6[sc6==i,i]
  w_avg6[i]=mean(ws)
}
sc7=apply(aic_d80_set[[7]],1,which.min)
cn7=summary(factor(sc7,levels=seq(1:5)))/2000
w7=matrix(ncol=5,nrow=nrow(aic_d80_set[[7]]))
for(i in 1:nrow(aic_d80_set[[7]])){
  w7[i,]=unlist(samEst::model_weights(aic_d80_set[[7]][i,],form='AIC'))
}
w_avg7=numeric(5)
for(i in 1:5){
  ws=w7[sc7==i,i]
  w_avg7[i]=mean(ws)
}
sc8=apply(aic_d80_set[[8]],1,which.min)
cn8=summary(factor(sc8,levels=seq(1:5)))/2000
w8=matrix(ncol=5,nrow=nrow(aic_d80_set[[8]]))
for(i in 1:nrow(aic_d80_set[[8]])){
  w8[i,]=unlist(samEst::model_weights(aic_d80_set[[8]][i,],form='AIC'))
}
w_avg8=numeric(5)
for(i in 1:5){
  ws=w8[sc8==i,i]
  w_avg8[i]=mean(ws)
}

##Confusion matrices
conf_matrix<-expand.grid(EM=c("stationary",
                              "autocorr",
                              "dynamic.a","dynamic.b","dynamic.ab"
),OM=c("stationary",
       "autocorr",
       "dec.prod","dec.cap","dec.prodcap",
       "regime.prod","regime.cap","regime.prodcap"))

conf_matrix$w_AIC=NA
conf_matrix$w_AIC[1:5]=cn1
conf_matrix$w_AIC[6:10]=cn2
conf_matrix$w_AIC[11:15]=cn3
conf_matrix$w_AIC[16:20]=cn4
conf_matrix$w_AIC[21:25]=cn5
conf_matrix$w_AIC[26:30]=cn6
conf_matrix$w_AIC[31:35]=cn7
conf_matrix$w_AIC[36:40]=cn8
conf_matrix$top_w[1:5]=w_avg1
conf_matrix$top_w[6:10]=w_avg2
conf_matrix$top_w[11:15]=w_avg3
conf_matrix$top_w[16:20]=w_avg4
conf_matrix$top_w[21:25]=w_avg5
conf_matrix$top_w[26:30]=w_avg6
conf_matrix$top_w[31:35]=w_avg7
conf_matrix$top_w[36:40]=w_avg8

mytheme = list(
  theme_classic(14)+
    theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
          legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="bold", size=13),
          axis.text=element_text(face="bold"),axis.title = element_text(face="bold",size=16),plot.title = element_text(face = "bold", hjust = 0.5,size=16))
)

library(ggplot2)



p=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = w_AIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_AIC,2)), vjust = 1,size=6) +
  ggtitle("AIC")+
  scale_fill_gradient(low = "white", high = "navy") +
  mytheme + theme(legend.position="none")+xlab("Simulation Scenario")+ylab("Estimation Model")
pdf(file='AIC_conf_matrix_d80.pdf',height=8,width=12)
p
dev.off()

p2=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = top_w), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_AIC,2)), vjust = 1,size=6) +
  ggtitle("AIC")+
  scale_fill_gradient(low = "white", high = "skyblue") +
  mytheme + theme(legend.position="right")+xlab("Simulation Scenario")+ylab("Estimation Model")
p2


conf_matrix2<-expand.grid(EM=c("stationary",
                               "autocorr",
                               "dynamic.a","dynamic.b","dynamic.ab"
),OM=c("stationary",
       "autocorr",
       "regime.prod","regime.cap","regime.prodcap"))

conf_matrix2$w_AIC=NA
conf_matrix2$w_AIC[1:5]=cn1
conf_matrix2$w_AIC[6:10]=cn2
conf_matrix2$w_AIC[11:15]=cn6
conf_matrix2$w_AIC[16:20]=cn7
conf_matrix2$w_AIC[21:25]=cn8
conf_matrix2$top_w[1:5]=w_avg1
conf_matrix2$top_w[6:10]=w_avg2
conf_matrix2$top_w[11:15]=w_avg6
conf_matrix2$top_w[16:20]=w_avg7
conf_matrix2$top_w[21:25]=w_avg8

p=ggplot(data =  conf_matrix2, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = w_AIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_AIC,2)), vjust = 1,size=6) +
  ggtitle("AIC")+
  scale_fill_gradient(low = "white", high = "navy") +
  mytheme + theme(legend.position="none")+xlab("Simulation Scenario")+ylab("Estimation Model")
p

p2=ggplot(data =  conf_matrix2, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = top_w), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_AIC,2)), vjust = 1,size=6) +
  ggtitle("AIC")+
  scale_fill_gradient(low = "white", high = "skyblue") +
  mytheme + theme(legend.position="right")+xlab("Simulation Scenario")+ylab("Estimation Model")
p2







#AIC npar2####
#first set of scenarios
sc1=apply(aic_npar2_set[[1]],1,which.min)
cn1=summary(factor(sc1,levels=seq(1:5)))/2000
w1=matrix(ncol=5,nrow=nrow(aic_npar2_set[[1]]))
for(i in 1:nrow(aic_npar2_set[[1]])){
  w1[i,]=unlist(samEst::model_weights(aic_npar2_set[[1]][i,],form='AIC'))
}
w_avg1=numeric(5)
for(i in 1:5){
  ws=w1[sc1==i,i]
  w_avg1[i]=mean(ws)
}

sc2=apply(aic_npar2_set[[2]],1,which.min)
cn2=summary(factor(sc2,levels=seq(1:5)))/2000
w2=matrix(ncol=5,nrow=nrow(aic_npar2_set[[2]]))
for(i in 1:nrow(aic_npar2_set[[2]])){
  w2[i,]=unlist(samEst::model_weights(aic_npar2_set[[2]][i,],form='AIC'))
}
w_avg2=numeric(5)
for(i in 1:5){
  ws=w2[sc2==i,i]
  w_avg2[i]=mean(ws)
}

sc3=apply(aic_npar2_set[[3]],1,which.min)
cn3=summary(factor(sc3,levels=seq(1:5)))/2000
w3=matrix(ncol=5,nrow=nrow(aic_npar2_set[[3]]))
for(i in 1:nrow(aic_npar2_set[[3]])){
  w3[i,]=unlist(samEst::model_weights(aic_npar2_set[[3]][i,],form='AIC'))
}
w_avg3=numeric(5)
for(i in 1:5){
  ws=w3[sc3==i,i]
  w_avg3[i]=mean(ws)
}
sc4=apply(aic_npar2_set[[4]],1,which.min)
cn4=summary(factor(sc4,levels=seq(1:5)))/2000
w4=matrix(ncol=5,nrow=nrow(aic_npar2_set[[4]]))
for(i in 1:nrow(aic_npar2_set[[4]])){
  w4[i,]=unlist(samEst::model_weights(aic_npar2_set[[4]][i,],form='AIC'))
}
w_avg4=numeric(5)
for(i in 1:5){
  ws=w4[sc4==i,i]
  w_avg4[i]=mean(ws)
}
sc5=apply(aic_npar2_set[[5]],1,which.min)
cn5=summary(factor(sc5,levels=seq(1:5)))/2000
w5=matrix(ncol=5,nrow=nrow(aic_npar2_set[[5]]))
for(i in 1:nrow(aic_npar2_set[[5]])){
  w5[i,]=unlist(samEst::model_weights(aic_npar2_set[[5]][i,],form='AIC'))
}
w_avg5=numeric(5)
for(i in 1:5){
  ws=w5[sc5==i,i]
  w_avg5[i]=mean(ws)
}
#second set of scenarios
sc6=apply(aic_npar2_set[[6]],1,which.min)
cn6=summary(factor(sc6,levels=seq(1:5)))/2000
w6=matrix(ncol=5,nrow=nrow(aic_npar2_set[[6]]))
for(i in 1:nrow(aic_npar2_set[[6]])){
  w6[i,]=unlist(samEst::model_weights(aic_npar2_set[[6]][i,],form='AIC'))
}
w_avg6=numeric(5)
for(i in 1:5){
  ws=w6[sc6==i,i]
  w_avg6[i]=mean(ws)
}
sc7=apply(aic_npar2_set[[7]],1,which.min)
cn7=summary(factor(sc7,levels=seq(1:5)))/2000
w7=matrix(ncol=5,nrow=nrow(aic_npar2_set[[7]]))
for(i in 1:nrow(aic_npar2_set[[7]])){
  w7[i,]=unlist(samEst::model_weights(aic_npar2_set[[7]][i,],form='AIC'))
}
w_avg7=numeric(5)
for(i in 1:5){
  ws=w7[sc7==i,i]
  w_avg7[i]=mean(ws)
}
sc8=apply(aic_npar2_set[[8]],1,which.min)
cn8=summary(factor(sc8,levels=seq(1:5)))/2000
w8=matrix(ncol=5,nrow=nrow(aic_npar2_set[[8]]))
for(i in 1:nrow(aic_npar2_set[[8]])){
  w8[i,]=unlist(samEst::model_weights(aic_npar2_set[[8]][i,],form='AIC'))
}
w_avg8=numeric(5)
for(i in 1:5){
  ws=w8[sc8==i,i]
  w_avg8[i]=mean(ws)
}

##Confusion matrices
conf_matrix<-expand.grid(EM=c("stationary",
                              "autocorr",
                              "dynamic.a","dynamic.b","dynamic.ab"
),OM=c("stationary",
       "autocorr",
       "dec.prod","dec.cap","dec.prodcap",
       "regime.prod","regime.cap","regime.prodcap"))

conf_matrix$w_AIC=NA
conf_matrix$w_AIC[1:5]=cn1
conf_matrix$w_AIC[6:10]=cn2
conf_matrix$w_AIC[11:15]=cn3
conf_matrix$w_AIC[16:20]=cn4
conf_matrix$w_AIC[21:25]=cn5
conf_matrix$w_AIC[26:30]=cn6
conf_matrix$w_AIC[31:35]=cn7
conf_matrix$w_AIC[36:40]=cn8
conf_matrix$top_w[1:5]=w_avg1
conf_matrix$top_w[6:10]=w_avg2
conf_matrix$top_w[11:15]=w_avg3
conf_matrix$top_w[16:20]=w_avg4
conf_matrix$top_w[21:25]=w_avg5
conf_matrix$top_w[26:30]=w_avg6
conf_matrix$top_w[31:35]=w_avg7
conf_matrix$top_w[36:40]=w_avg8


mytheme = list(
  theme_classic(14)+
    theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
          legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="bold", size=13),
          axis.text=element_text(face="bold"),axis.title = element_text(face="bold",size=16),plot.title = element_text(face = "bold", hjust = 0.5,size=16))
)

library(ggplot2)



p=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = w_AIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_AIC,2)), vjust = 1,size=6) +
  ggtitle("AIC")+
  scale_fill_gradient(low = "white", high = "navy") +
  mytheme + theme(legend.position="none")+xlab("Simulation Scenario")+ylab("Estimation Model")

pdf(file='AIC_conf_matrix_k2.pdf',height=8,width=12)
p
dev.off()


p2=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = top_w), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_AIC,2)), vjust = 1,size=6) +
  ggtitle("AIC")+
  scale_fill_gradient(low = "white", high = "skyblue") +
  mytheme + theme(legend.position="right")+xlab("Simulation Scenario")+ylab("Estimation Model")
p2


conf_matrix2<-expand.grid(EM=c("stationary",
                               "autocorr",
                               "dynamic.a","dynamic.b","dynamic.ab"
),OM=c("stationary",
       "autocorr",
       "regime.prod","regime.cap","regime.prodcap"))

conf_matrix2$w_AIC=NA
conf_matrix2$w_AIC[1:5]=cn1
conf_matrix2$w_AIC[6:10]=cn2
conf_matrix2$w_AIC[11:15]=cn6
conf_matrix2$w_AIC[16:20]=cn7
conf_matrix2$w_AIC[21:25]=cn8
conf_matrix2$top_w[1:5]=w_avg1
conf_matrix2$top_w[6:10]=w_avg2
conf_matrix2$top_w[11:15]=w_avg6
conf_matrix2$top_w[16:20]=w_avg7
conf_matrix2$top_w[21:25]=w_avg8

p=ggplot(data =  conf_matrix2, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = w_AIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_AIC,2)), vjust = 1,size=6) +
  ggtitle("AIC")+
  scale_fill_gradient(low = "white", high = "navy") +
  mytheme + theme(legend.position="none")+xlab("Simulation Scenario")+ylab("Estimation Model")
p

p2=ggplot(data =  conf_matrix2, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = top_w), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_AIC,2)), vjust = 1,size=6) +
  ggtitle("AIC")+
  scale_fill_gradient(low = "white", high = "skyblue") +
  mytheme + theme(legend.position="right")+xlab("Simulation Scenario")+ylab("Estimation Model")
p2






#BIC 1####
#first set of scenarios
sc1=apply(bic_set[[1]],1,which.min)
cn1=summary(factor(sc1,levels=seq(1:5)))/2000
w1=matrix(ncol=5,nrow=nrow(bic_set[[1]]))
for(i in 1:nrow(bic_set[[1]])){
  w1[i,]=unlist(samEst::model_weights(bic_set[[1]][i,],form='AIC'))
}
w_avg1=numeric(5)
for(i in 1:5){
  ws=w1[sc1==i,i]
  w_avg1[i]=mean(ws)
}

sc2=apply(bic_set[[2]],1,which.min)
cn2=summary(factor(sc2,levels=seq(1:5)))/2000
w2=matrix(ncol=5,nrow=nrow(bic_set[[2]]))
for(i in 1:nrow(bic_set[[2]])){
  w2[i,]=unlist(samEst::model_weights(bic_set[[2]][i,],form='AIC'))
}
w_avg2=numeric(5)
for(i in 1:5){
  ws=w2[sc2==i,i]
  w_avg2[i]=mean(ws)
}

sc3=apply(bic_set[[3]],1,which.min)
cn3=summary(factor(sc3,levels=seq(1:5)))/2000
w3=matrix(ncol=5,nrow=nrow(bic_set[[3]]))
for(i in 1:nrow(bic_set[[3]])){
  w3[i,]=unlist(samEst::model_weights(bic_set[[3]][i,],form='AIC'))
}
w_avg3=numeric(5)
for(i in 1:5){
  ws=w3[sc3==i,i]
  w_avg3[i]=mean(ws)
}
sc4=apply(bic_set[[4]],1,which.min)
cn4=summary(factor(sc4,levels=seq(1:5)))/2000
w4=matrix(ncol=5,nrow=nrow(bic_set[[4]]))
for(i in 1:nrow(bic_set[[4]])){
  w4[i,]=unlist(samEst::model_weights(bic_set[[4]][i,],form='AIC'))
}
w_avg4=numeric(5)
for(i in 1:5){
  ws=w4[sc4==i,i]
  w_avg4[i]=mean(ws)
}
sc5=apply(bic_set[[5]],1,which.min)
cn5=summary(factor(sc5,levels=seq(1:5)))/2000
w5=matrix(ncol=5,nrow=nrow(bic_set[[5]]))
for(i in 1:nrow(bic_set[[5]])){
  w5[i,]=unlist(samEst::model_weights(bic_set[[5]][i,],form='AIC'))
}
w_avg5=numeric(5)
for(i in 1:5){
  ws=w5[sc5==i,i]
  w_avg5[i]=mean(ws)
}
#second set of scenarios
sc6=apply(bic_set[[6]],1,which.min)
cn6=summary(factor(sc6,levels=seq(1:5)))/2000
w6=matrix(ncol=5,nrow=nrow(bic_set[[6]]))
for(i in 1:nrow(bic_set[[6]])){
  w6[i,]=unlist(samEst::model_weights(bic_set[[6]][i,],form='AIC'))
}
w_avg6=numeric(5)
for(i in 1:5){
  ws=w6[sc6==i,i]
  w_avg6[i]=mean(ws)
}
sc7=apply(bic_set[[7]],1,which.min)
cn7=summary(factor(sc7,levels=seq(1:5)))/2000
w7=matrix(ncol=5,nrow=nrow(bic_set[[7]]))
for(i in 1:nrow(bic_set[[7]])){
  w7[i,]=unlist(samEst::model_weights(bic_set[[7]][i,],form='AIC'))
}
w_avg7=numeric(5)
for(i in 1:5){
  ws=w7[sc7==i,i]
  w_avg7[i]=mean(ws)
}
sc8=apply(bic_set[[8]],1,which.min)
cn8=summary(factor(sc8,levels=seq(1:5)))/2000
w8=matrix(ncol=5,nrow=nrow(bic_set[[8]]))
for(i in 1:nrow(bic_set[[8]])){
  w8[i,]=unlist(samEst::model_weights(bic_set[[8]][i,],form='AIC'))
}
w_avg8=numeric(5)
for(i in 1:5){
  ws=w8[sc8==i,i]
  w_avg8[i]=mean(ws)
}


conf_matrix<-expand.grid(EM=c("stationary",
                              "autocorr",
                              "dynamic.a","dynamic.b","dynamic.ab"
),OM=c("stationary",
       "autocorr",
       "dec.prod","dec.cap","dec.prodcap","regime.prod","regime.cap","regime.prodcap"))

conf_matrix$w_BIC=NA
conf_matrix$w_BIC[1:5]=cn1
conf_matrix$w_BIC[6:10]=cn2
conf_matrix$w_BIC[11:15]=cn3
conf_matrix$w_BIC[16:20]=cn4
conf_matrix$w_BIC[21:25]=cn5
conf_matrix$w_BIC[26:30]=cn6
conf_matrix$w_BIC[31:35]=cn7
conf_matrix$w_BIC[36:40]=cn8

conf_matrix$top_w[1:5]=w_avg1
conf_matrix$top_w[6:10]=w_avg2
conf_matrix$top_w[11:15]=w_avg3
conf_matrix$top_w[16:20]=w_avg4
conf_matrix$top_w[21:25]=w_avg5
conf_matrix$top_w[26:30]=w_avg6
conf_matrix$top_w[31:35]=w_avg7
conf_matrix$top_w[36:40]=w_avg8


p=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = w_BIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_BIC,2)), vjust = 1) +
  ggtitle("BIC")+
  scale_fill_gradient(low = "white", high = "navy") +
  mytheme + theme(legend.position="none")+xlab("Simulation Scenario")+ylab("Estimation Model")

p

p2=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = top_w), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_BIC,2)), vjust = 1,size=6) +
  ggtitle("AIC")+
  scale_fill_gradient(low = "white", high = "skyblue") +
  mytheme + theme(legend.position="right")+xlab("Simulation Scenario")+ylab("Estimation Model")
p2

conf_matrix2<-expand.grid(EM=c("stationary",
                               "autocorr",
                               "dynamic.a","dynamic.b","dynamic.ab"
),OM=c("stationary",
       "autocorr",
       "regime.prod","regime.cap","regime.prodcap"))

conf_matrix2$w_BIC=NA
conf_matrix2$w_BIC[1:5]=cn1
conf_matrix2$w_BIC[6:10]=cn2
conf_matrix2$w_BIC[11:15]=cn6
conf_matrix2$w_BIC[16:20]=cn7
conf_matrix2$w_BIC[21:25]=cn8
conf_matrix2$top_w[1:5]=w_avg1
conf_matrix2$top_w[6:10]=w_avg2
conf_matrix2$top_w[11:15]=w_avg6
conf_matrix2$top_w[16:20]=w_avg7
conf_matrix2$top_w[21:25]=w_avg8

p=ggplot(data =  conf_matrix2, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = w_AIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_AIC,2)), vjust = 1,size=6) +
  ggtitle("AIC")+
  scale_fill_gradient(low = "white", high = "navy") +
  mytheme + theme(legend.position="none")+xlab("Simulation Scenario")+ylab("Estimation Model")
p

p2=ggplot(data =  conf_matrix2, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = top_w), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_AIC,2)), vjust = 1,size=6) +
  ggtitle("AIC")+
  scale_fill_gradient(low = "white", high = "skyblue") +
  mytheme + theme(legend.position="right")+xlab("Simulation Scenario")+ylab("Estimation Model")
p2


#By model class
sc1m=ifelse(sc1<3,1,sc1)
sc1m=ifelse(sc1>2&sc1<6,2,sc1m)
sc1m=ifelse(sc1>5,3,sc1m)
cn1m=summary(factor(sc1m),levels=seq(1:3))/1000
sc2m=ifelse(sc2<3,1,sc2)
sc2m=ifelse(sc2>2&sc2<6,2,sc2m)
sc2m=ifelse(sc2>5,3,sc2m)
cn2m=summary(factor(sc2m),levels=seq(1:3))/1000
sc3m=ifelse(sc3<3,1,sc3)
sc3m=ifelse(sc3>2&sc3<6,2,sc3m)
sc3m=ifelse(sc3>5,3,sc3m)
cn3m=summary(factor(sc3m),levels=seq(1:3))/1000
sc4m=ifelse(sc4<3,1,sc4)
sc4m=ifelse(sc4>2&sc4<6,2,sc4m)
sc4m=ifelse(sc4>5,3,sc4m)
cn4m=summary(factor(sc4m),levels=seq(1:3))/1000
sc5m=ifelse(sc5<3,1,sc5)
sc5m=ifelse(sc5>2&sc5<6,2,sc5m)
sc5m=ifelse(sc5>5,3,sc5m)
cn5m=summary(factor(sc5m),levels=seq(1:3))/1000
sc6m=ifelse(sc6<3,1,sc6)
sc6m=ifelse(sc6>2&sc6<6,2,sc6m)
sc6m=ifelse(sc6>5,3,sc6m)
cn6m=summary(factor(sc6m),levels=seq(1:3))/1000
sc7m=ifelse(sc7<3,1,sc7)
sc7m=ifelse(sc7>2&sc7<6,2,sc7m)
sc7m=ifelse(sc7>5,3,sc7m)
cn7m=summary(factor(sc7m),levels=seq(1:3))/1000
sc8m=ifelse(sc8<3,1,sc8)
sc8m=ifelse(sc8>2&sc8<6,2,sc8m)
sc8m=ifelse(sc8>5,3,sc8m)
cn8m=summary(factor(sc8m),levels=seq(1:3))/1000

sx1m=ifelse(sx1<3,1,sx1)
sx1m=ifelse(sx1>2&sx1<6,2,sx1m)
sx1m=ifelse(sx1>5,3,sx1m)
ck1m=summary(factor(sx1m),levels=seq(1:3))/1000
sx2m=ifelse(sx2<3,1,sx2)
sx2m=ifelse(sx2>2&sx2<6,2,sx2m)
sx2m=ifelse(sx2>5,3,sx2m)
ck2m=summary(factor(sx2m),levels=seq(1:3))/1000
sx3m=ifelse(sx3<3,1,sx3)
sx3m=ifelse(sx3>2&sx3<6,2,sx3m)
sx3m=ifelse(sx3>5,3,sx3m)
ck3m=summary(factor(sx3m,levels=seq(1:3)))/1000
sx4m=ifelse(sx4<3,1,sx4)
sx4m=ifelse(sx4>2&sx4<6,2,sx4m)
sx4m=ifelse(sx4>5,3,sx4m)
ck4m=summary(factor(sx4m),levels=seq(1:3))/1000
sx5m=ifelse(sx5<3,1,sx5)
sx5m=ifelse(sx5>2&sx5<6,2,sx5m)
sx5m=ifelse(sx5>5,3,sx5m)
ck5m=summary(factor(sx5m),levels=seq(1:3))/1000
sx6m=ifelse(sx6<3,1,sx6)
sx6m=ifelse(sx6>2&sx6<6,2,sx6m)
sx6m=ifelse(sx6>5,3,sx6m)
ck6m=summary(factor(sx6m),levels=seq(1:3))/1000
sx7m=ifelse(sx7<3,1,sx7)
sx7m=ifelse(sx7>2&sx7<6,2,sx7m)
sx7m=ifelse(sx7>5,3,sx7m)
ck7m=summary(factor(sx7m),levels=seq(1:3))/1000
sx8m=ifelse(sx8<3,1,sx8)
sx8m=ifelse(sx8>2&sx8<6,2,sx8m)
sx8m=ifelse(sx8>5,3,sx8m)
ck8m=summary(factor(sx8m),levels=seq(1:3))/1000

##Confusion matrices
conf_matrix <-expand.grid(EM=c('static','dynamic','regime'),OM=c("stat",
                                                                 "acorr",
                                                                 "dyn.a","dyn.b","dyn.ab",
                                                                 "reg.a", "reg.b","reg.ab"))

conf_matrix$w_AIC=NA
conf_matrix$w_AIC[1:3]=cn1m
conf_matrix$w_AIC[4:6]=cn2m
conf_matrix$w_AIC[7:9]=cn3m
conf_matrix$w_AIC[10:12]=cn4m
conf_matrix$w_AIC[13:15]=cn5m
conf_matrix$w_AIC[16:18]=cn6m
conf_matrix$w_AIC[19:21]=cn7m
conf_matrix$w_AIC[22:24]=cn8m
conf_matrix$w_BIC=NA
conf_matrix$w_BIC[1:3]=ck1m
conf_matrix$w_BIC[4:6]=ck2m
conf_matrix$w_BIC[7:9]=ck3m
conf_matrix$w_BIC[10:12]=ck4m
conf_matrix$w_BIC[13:15]=ck5m
conf_matrix$w_BIC[16:18]=ck6m
conf_matrix$w_BIC[19:21]=ck7m
conf_matrix$w_BIC[22:24]=ck8m

library(ggplot2)
p=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = w_AIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_AIC,2)), vjust = 1) +
  ggtitle("AIC")+
  scale_fill_gradient(low = "white", high = "navy") +
  mytheme + theme(legend.position="none")+xlab("Operating Model")+ylab("Estimation Model")
p

b=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = w_BIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_BIC,2)), vjust = 1) +
  ggtitle("BIC")+
  scale_fill_gradient(low = "white", high = "navy") +
  mytheme + theme(legend.position="none")+xlab("Operating Model")+ylab("Estimation Model")

b

ps=cowplot::plot_grid(p,b,
                      ncol=2,nrow=1,labels=c("A","B"),label_size=10)

ggsave(filename = "outs/conf_mat_mc_tmb_aicbic.pdf",
       plot=ps,
       width=12,height=4)
ggsave(filename = "outs/conf_mat_mc_tmb_aic.pdf",
       plot=p,
       width=8,height=5)

ggsave(filename = "outs/conf_mat_mc_tmb_bic.pdf",
       plot=b,
       width=12,height=5)




#BIC 2####
#first set of scenarios
sc1=apply(bic_set2[[1]],1,which.min)
cn1=summary(factor(sc1,levels=seq(1:5)))/2000
w1=matrix(ncol=5,nrow=nrow(bic_set2[[1]]))
for(i in 1:nrow(bic_set2[[1]])){
  w1[i,]=unlist(samEst::model_weights(bic_set2[[1]][i,],form='AIC'))
}
w_avg1=numeric(5)
for(i in 1:5){
  ws=w1[sc1==i,i]
  w_avg1[i]=mean(ws)
}

sc2=apply(bic_set2[[2]],1,which.min)
cn2=summary(factor(sc2,levels=seq(1:5)))/2000
w2=matrix(ncol=5,nrow=nrow(bic_set2[[2]]))
for(i in 1:nrow(bic_set2[[2]])){
  w2[i,]=unlist(samEst::model_weights(bic_set2[[2]][i,],form='AIC'))
}
w_avg2=numeric(5)
for(i in 1:5){
  ws=w2[sc2==i,i]
  w_avg2[i]=mean(ws)
}

sc3=apply(bic_set2[[3]],1,which.min)
cn3=summary(factor(sc3,levels=seq(1:5)))/2000
w3=matrix(ncol=5,nrow=nrow(bic_set2[[3]]))
for(i in 1:nrow(bic_set2[[3]])){
  w3[i,]=unlist(samEst::model_weights(bic_set2[[3]][i,],form='AIC'))
}
w_avg3=numeric(5)
for(i in 1:5){
  ws=w3[sc3==i,i]
  w_avg3[i]=mean(ws)
}
sc4=apply(bic_set2[[4]],1,which.min)
cn4=summary(factor(sc4,levels=seq(1:5)))/2000
w4=matrix(ncol=5,nrow=nrow(bic_set2[[4]]))
for(i in 1:nrow(bic_set2[[4]])){
  w4[i,]=unlist(samEst::model_weights(bic_set2[[4]][i,],form='AIC'))
}
w_avg4=numeric(5)
for(i in 1:5){
  ws=w4[sc4==i,i]
  w_avg4[i]=mean(ws)
}
sc5=apply(bic_set2[[5]],1,which.min)
cn5=summary(factor(sc5,levels=seq(1:5)))/2000
w5=matrix(ncol=5,nrow=nrow(bic_set2[[5]]))
for(i in 1:nrow(bic_set2[[5]])){
  w5[i,]=unlist(samEst::model_weights(bic_set2[[5]][i,],form='AIC'))
}
w_avg5=numeric(5)
for(i in 1:5){
  ws=w5[sc5==i,i]
  w_avg5[i]=mean(ws)
}
#second set of scenarios
sc6=apply(bic_set2[[6]],1,which.min)
cn6=summary(factor(sc6,levels=seq(1:5)))/2000
w6=matrix(ncol=5,nrow=nrow(bic_set2[[6]]))
for(i in 1:nrow(bic_set2[[6]])){
  w6[i,]=unlist(samEst::model_weights(bic_set2[[6]][i,],form='AIC'))
}
w_avg6=numeric(5)
for(i in 1:5){
  ws=w6[sc6==i,i]
  w_avg6[i]=mean(ws)
}
sc7=apply(bic_set2[[7]],1,which.min)
cn7=summary(factor(sc7,levels=seq(1:5)))/2000
w7=matrix(ncol=5,nrow=nrow(bic_set2[[7]]))
for(i in 1:nrow(bic_set2[[7]])){
  w7[i,]=unlist(samEst::model_weights(bic_set2[[7]][i,],form='AIC'))
}
w_avg7=numeric(5)
for(i in 1:5){
  ws=w7[sc7==i,i]
  w_avg7[i]=mean(ws)
}
sc8=apply(bic_set2[[8]],1,which.min)
cn8=summary(factor(sc8,levels=seq(1:5)))/2000
w8=matrix(ncol=5,nrow=nrow(bic_set2[[8]]))
for(i in 1:nrow(bic_set2[[8]])){
  w8[i,]=unlist(samEst::model_weights(bic_set2[[8]][i,],form='AIC'))
}
w_avg8=numeric(5)
for(i in 1:5){
  ws=w8[sc8==i,i]
  w_avg8[i]=mean(ws)
}


conf_matrix<-expand.grid(EM=c("stationary",
                              "autocorr",
                              "dynamic.a","dynamic.b","dynamic.ab"
),OM=c("stationary",
       "autocorr",
       "dec.prod","dec.cap","dec.prodcap","regime.prod","regime.cap","regime.prodcap"))

conf_matrix$w_BIC=NA
conf_matrix$w_BIC[1:5]=cn1
conf_matrix$w_BIC[6:10]=cn2
conf_matrix$w_BIC[11:15]=cn3
conf_matrix$w_BIC[16:20]=cn4
conf_matrix$w_BIC[21:25]=cn5
conf_matrix$w_BIC[26:30]=cn6
conf_matrix$w_BIC[31:35]=cn7
conf_matrix$w_BIC[36:40]=cn8

conf_matrix$top_w[1:5]=w_avg1
conf_matrix$top_w[6:10]=w_avg2
conf_matrix$top_w[11:15]=w_avg3
conf_matrix$top_w[16:20]=w_avg4
conf_matrix$top_w[21:25]=w_avg5
conf_matrix$top_w[26:30]=w_avg6
conf_matrix$top_w[31:35]=w_avg7
conf_matrix$top_w[36:40]=w_avg8


p=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = w_BIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_BIC,2)), vjust = 1) +
  ggtitle("BIC")+
  scale_fill_gradient(low = "white", high = "navy") +
  mytheme + theme(legend.position="none")+xlab("Simulation Scenario")+ylab("Estimation Model")

pdf(file='BIC_conf_matrix.pdf',height=8,width=12)
p
dev.off()

p2=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = top_w), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_BIC,2)), vjust = 1,size=6) +
  ggtitle("AIC")+
  scale_fill_gradient(low = "white", high = "skyblue") +
  mytheme + theme(legend.position="right")+xlab("Simulation Scenario")+ylab("Estimation Model")
p2

conf_matrix2<-expand.grid(EM=c("stationary",
                               "autocorr",
                               "dynamic.a","dynamic.b","dynamic.ab"
),OM=c("stationary",
       "autocorr",
       "regime.prod","regime.cap","regime.prodcap"))

conf_matrix2$w_BIC=NA
conf_matrix2$w_BIC[1:5]=cn1
conf_matrix2$w_BIC[6:10]=cn2
conf_matrix2$w_BIC[11:15]=cn6
conf_matrix2$w_BIC[16:20]=cn7
conf_matrix2$w_BIC[21:25]=cn8
conf_matrix2$top_w[1:5]=w_avg1
conf_matrix2$top_w[6:10]=w_avg2
conf_matrix2$top_w[11:15]=w_avg6
conf_matrix2$top_w[16:20]=w_avg7
conf_matrix2$top_w[21:25]=w_avg8

p=ggplot(data =  conf_matrix2, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = w_AIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_BIC,2)), vjust = 1,size=6) +
  ggtitle("AIC")+
  scale_fill_gradient(low = "white", high = "navy") +
  mytheme + theme(legend.position="none")+xlab("Simulation Scenario")+ylab("Estimation Model")
p

p2=ggplot(data =  conf_matrix2, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = top_w), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_AIC,2)), vjust = 1,size=6) +
  ggtitle("AIC")+
  scale_fill_gradient(low = "white", high = "skyblue") +
  mytheme + theme(legend.position="right")+xlab("Simulation Scenario")+ylab("Estimation Model")
p2



sx1m=ifelse(sx1<3,1,sx1)
sx1m=ifelse(sx1>2&sx1<6,2,sx1m)
sx1m=ifelse(sx1>5,3,sx1m)
ck1m=summary(factor(sx1m),levels=seq(1:3))/1000
sx2m=ifelse(sx2<3,1,sx2)
sx2m=ifelse(sx2>2&sx2<6,2,sx2m)
sx2m=ifelse(sx2>5,3,sx2m)
ck2m=summary(factor(sx2m),levels=seq(1:3))/1000
sx3m=ifelse(sx3<3,1,sx3)
sx3m=ifelse(sx3>2&sx3<6,2,sx3m)
sx3m=ifelse(sx3>5,3,sx3m)
ck3m=summary(factor(sx3m,levels=seq(1:3)))/1000
sx4m=ifelse(sx4<3,1,sx4)
sx4m=ifelse(sx4>2&sx4<6,2,sx4m)
sx4m=ifelse(sx4>5,3,sx4m)
ck4m=summary(factor(sx4m),levels=seq(1:3))/1000
sx5m=ifelse(sx5<3,1,sx5)
sx5m=ifelse(sx5>2&sx5<6,2,sx5m)
sx5m=ifelse(sx5>5,3,sx5m)
ck5m=summary(factor(sx5m),levels=seq(1:3))/1000
sx6m=ifelse(sx6<3,1,sx6)
sx6m=ifelse(sx6>2&sx6<6,2,sx6m)
sx6m=ifelse(sx6>5,3,sx6m)
ck6m=summary(factor(sx6m),levels=seq(1:3))/1000
sx7m=ifelse(sx7<3,1,sx7)
sx7m=ifelse(sx7>2&sx7<6,2,sx7m)
sx7m=ifelse(sx7>5,3,sx7m)
ck7m=summary(factor(sx7m),levels=seq(1:3))/1000
sx8m=ifelse(sx8<3,1,sx8)
sx8m=ifelse(sx8>2&sx8<6,2,sx8m)
sx8m=ifelse(sx8>5,3,sx8m)
ck8m=summary(factor(sx8m),levels=seq(1:3))/1000

##Confusion matrices
conf_matrix <-expand.grid(EM=c('static','dynamic','regime'),OM=c("stat",
                                                                 "acorr",
                                                                 "dyn.a","dyn.b","dyn.ab",
                                                                 "reg.a", "reg.b","reg.ab"))

conf_matrix$w_AIC=NA
conf_matrix$w_AIC[1:3]=cn1m
conf_matrix$w_AIC[4:6]=cn2m
conf_matrix$w_AIC[7:9]=cn3m
conf_matrix$w_AIC[10:12]=cn4m
conf_matrix$w_AIC[13:15]=cn5m
conf_matrix$w_AIC[16:18]=cn6m
conf_matrix$w_AIC[19:21]=cn7m
conf_matrix$w_AIC[22:24]=cn8m
conf_matrix$w_BIC=NA
conf_matrix$w_BIC[1:3]=ck1m
conf_matrix$w_BIC[4:6]=ck2m
conf_matrix$w_BIC[7:9]=ck3m
conf_matrix$w_BIC[10:12]=ck4m
conf_matrix$w_BIC[13:15]=ck5m
conf_matrix$w_BIC[16:18]=ck6m
conf_matrix$w_BIC[19:21]=ck7m
conf_matrix$w_BIC[22:24]=ck8m

library(ggplot2)
p=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = w_AIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_AIC,2)), vjust = 1) +
  ggtitle("AIC")+
  scale_fill_gradient(low = "white", high = "navy") +
  mytheme + theme(legend.position="none")+xlab("Operating Model")+ylab("Estimation Model")
p

b=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = w_BIC), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_BIC,2)), vjust = 1) +
  ggtitle("BIC")+
  scale_fill_gradient(low = "white", high = "navy") +
  mytheme + theme(legend.position="none")+xlab("Operating Model")+ylab("Estimation Model")

b

ps=cowplot::plot_grid(p,b,
                      ncol=2,nrow=1,labels=c("A","B"),label_size=10)

ggsave(filename = "outs/conf_mat_mc_tmb_aicbic.pdf",
       plot=ps,
       width=12,height=4)
ggsave(filename = "outs/conf_mat_mc_tmb_aic.pdf",
       plot=p,
       width=8,height=5)

ggsave(filename = "outs/conf_mat_mc_tmb_bic.pdf",
       plot=b,
       width=12,height=5)


#BIC npar2####
#first set of scenarios
sc1=apply(bic_npar2_set[[1]],1,which.min)
cn1=summary(factor(sc1,levels=seq(1:5)))/2000
w1=matrix(ncol=5,nrow=nrow(bic_npar2_set[[1]]))
for(i in 1:nrow(bic_npar2_set[[1]])){
  w1[i,]=unlist(samEst::model_weights(bic_npar2_set[[1]][i,],form='AIC'))
}
w_avg1=numeric(5)
for(i in 1:5){
  ws=w1[sc1==i,i]
  w_avg1[i]=mean(ws)
}

sc2=apply(bic_npar2_set[[2]],1,which.min)
cn2=summary(factor(sc2,levels=seq(1:5)))/2000
w2=matrix(ncol=5,nrow=nrow(bic_npar2_set[[2]]))
for(i in 1:nrow(bic_npar2_set[[2]])){
  w2[i,]=unlist(samEst::model_weights(bic_npar2_set[[2]][i,],form='AIC'))
}
w_avg2=numeric(5)
for(i in 1:5){
  ws=w2[sc2==i,i]
  w_avg2[i]=mean(ws)
}

sc3=apply(bic_npar2_set[[3]],1,which.min)
cn3=summary(factor(sc3,levels=seq(1:5)))/2000
w3=matrix(ncol=5,nrow=nrow(bic_npar2_set[[3]]))
for(i in 1:nrow(bic_npar2_set[[3]])){
  w3[i,]=unlist(samEst::model_weights(bic_npar2_set[[3]][i,],form='AIC'))
}
w_avg3=numeric(5)
for(i in 1:5){
  ws=w3[sc3==i,i]
  w_avg3[i]=mean(ws)
}
sc4=apply(bic_npar2_set[[4]],1,which.min)
cn4=summary(factor(sc4,levels=seq(1:5)))/2000
w4=matrix(ncol=5,nrow=nrow(bic_npar2_set[[4]]))
for(i in 1:nrow(bic_npar2_set[[4]])){
  w4[i,]=unlist(samEst::model_weights(bic_npar2_set[[4]][i,],form='AIC'))
}
w_avg4=numeric(5)
for(i in 1:5){
  ws=w4[sc4==i,i]
  w_avg4[i]=mean(ws)
}
sc5=apply(bic_npar2_set[[5]],1,which.min)
cn5=summary(factor(sc5,levels=seq(1:5)))/2000
w5=matrix(ncol=5,nrow=nrow(bic_npar2_set[[5]]))
for(i in 1:nrow(bic_npar2_set[[5]])){
  w5[i,]=unlist(samEst::model_weights(bic_npar2_set[[5]][i,],form='AIC'))
}
w_avg5=numeric(5)
for(i in 1:5){
  ws=w5[sc5==i,i]
  w_avg5[i]=mean(ws)
}
#second set of scenarios
sc6=apply(bic_npar2_set[[6]],1,which.min)
cn6=summary(factor(sc6,levels=seq(1:5)))/2000
w6=matrix(ncol=5,nrow=nrow(bic_npar2_set[[6]]))
for(i in 1:nrow(bic_npar2_set[[6]])){
  w6[i,]=unlist(samEst::model_weights(bic_npar2_set[[6]][i,],form='AIC'))
}
w_avg6=numeric(5)
for(i in 1:5){
  ws=w6[sc6==i,i]
  w_avg6[i]=mean(ws)
}
sc7=apply(bic_npar2_set[[7]],1,which.min)
cn7=summary(factor(sc7,levels=seq(1:5)))/2000
w7=matrix(ncol=5,nrow=nrow(bic_npar2_set[[7]]))
for(i in 1:nrow(bic_npar2_set[[7]])){
  w7[i,]=unlist(samEst::model_weights(bic_npar2_set[[7]][i,],form='AIC'))
}
w_avg7=numeric(5)
for(i in 1:5){
  ws=w7[sc7==i,i]
  w_avg7[i]=mean(ws)
}
sc8=apply(bic_npar2_set[[8]],1,which.min)
cn8=summary(factor(sc8,levels=seq(1:5)))/2000
w8=matrix(ncol=5,nrow=nrow(bic_npar2_set[[8]]))
for(i in 1:nrow(bic_npar2_set[[8]])){
  w8[i,]=unlist(samEst::model_weights(bic_npar2_set[[8]][i,],form='AIC'))
}
w_avg8=numeric(5)
for(i in 1:5){
  ws=w8[sc8==i,i]
  w_avg8[i]=mean(ws)
}

##Confusion matrices
conf_matrix<-expand.grid(EM=c("stationary",
                              "autocorr",
                              "dynamic.a","dynamic.b","dynamic.ab"
),OM=c("stationary",
       "autocorr",
       "dec.prod","dec.cap","dec.prodcap"))

conf_matrix$w_bic=NA
conf_matrix$w_bic[1:5]=cn1
conf_matrix$w_bic[6:10]=cn2
conf_matrix$w_bic[11:15]=cn3
conf_matrix$w_bic[16:20]=cn4
conf_matrix$w_bic[21:25]=cn5
conf_matrix$top_w[1:5]=w_avg1
conf_matrix$top_w[6:10]=w_avg2
conf_matrix$top_w[11:15]=w_avg3
conf_matrix$top_w[16:20]=w_avg4
conf_matrix$top_w[21:25]=w_avg5

mytheme = list(
  theme_classic(14)+
    theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
          legend.title = element_blank(),legend.position="bottom", strip.text = element_text(face="bold", size=13),
          axis.text=element_text(face="bold"),axis.title = element_text(face="bold",size=16),plot.title = element_text(face = "bold", hjust = 0.5,size=16))
)

library(ggplot2)



p=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = w_bic), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_bic,2)), vjust = 1,size=6) +
  ggtitle("bic")+
  scale_fill_gradient(low = "white", high = "navy") +
  mytheme + theme(legend.position="none")+xlab("Simulation Scenario")+ylab("Estimation Model")
p

p2=ggplot(data =  conf_matrix, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = top_w), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_bic,2)), vjust = 1,size=6) +
  ggtitle("bic")+
  scale_fill_gradient(low = "white", high = "skyblue") +
  mytheme + theme(legend.position="right")+xlab("Simulation Scenario")+ylab("Estimation Model")
p2


conf_matrix2<-expand.grid(EM=c("stationary",
                               "autocorr",
                               "dynamic.a","dynamic.b","dynamic.ab"
),OM=c("stationary",
       "autocorr",
       "regime.prod","regime.cap","regime.prodcap"))

conf_matrix2$w_bic=NA
conf_matrix2$w_bic[1:5]=cn1
conf_matrix2$w_bic[6:10]=cn2
conf_matrix2$w_bic[11:15]=cn6
conf_matrix2$w_bic[16:20]=cn7
conf_matrix2$w_bic[21:25]=cn8
conf_matrix2$top_w[1:5]=w_avg1
conf_matrix2$top_w[6:10]=w_avg2
conf_matrix2$top_w[11:15]=w_avg6
conf_matrix2$top_w[16:20]=w_avg7
conf_matrix2$top_w[21:25]=w_avg8

p=ggplot(data =  conf_matrix2, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = w_bic), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_bic,2)), vjust = 1,size=6) +
  ggtitle("bic")+
  scale_fill_gradient(low = "white", high = "navy") +
  mytheme + theme(legend.position="none")+xlab("Simulation Scenario")+ylab("Estimation Model")
p

p2=ggplot(data =  conf_matrix2, mapping = aes(x = OM, y = EM)) +
  geom_tile(aes(fill = top_w), colour = "white",alpha=0.7) +
  geom_text(aes(label = round(w_bic,2)), vjust = 1,size=6) +
  ggtitle("bic")+
  scale_fill_gradient(low = "white", high = "skyblue") +
  mytheme + theme(legend.position="right")+xlab("Simulation Scenario")+ylab("Estimation Model")
p2



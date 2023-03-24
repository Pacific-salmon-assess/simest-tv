#
library(ggplot2)


#Alpha sensitivity scenarios
log_a<-1.3
#b<-0.00000581
b<-1.176471e-05


yrs<-40

log_af<-log(c(1,2,5,7,10))

stp<-(log_af-log_a)/(yrs-1)
log_a_mat<-matrix(NA,ncol=length(log_af),nrow=yrs)
log_a_mat[1,]<-1.3
for(n in 2:yrs){
   log_a_mat[n,]<-round(log_a_mat[n-1,] +stp,5)
}


S<-seq(0,800000,5000)

R<-matrix(NA,ncol=yrs,nrow=length(S))
R2<-matrix(NA,ncol=yrs,nrow=length(S))
R3<-matrix(NA,ncol=yrs,nrow=length(S))
R4<-matrix(NA,ncol=yrs,nrow=length(S))
R5<-matrix(NA,ncol=yrs,nrow=length(S))



    for(y in seq_len(yrs) ){
      R[,y]<-S*exp(log_a_mat[y,1]-b*S)
      R2[,y]<-S*exp(log_a_mat[y,2]-b*S)
      R3[,y]<-S*exp(log_a_mat[y,3]-b*S)
      R4[,y]<-S*exp(log_a_mat[y,4]-b*S)
      R5[,y]<-S*exp(log_a_mat[y,5]-b*S)
    }

dt<-data.frame(Recruits=c(c(R),c(R2),c(R3),c(R4),c(R5)),
    Spawner=S,
    af=rep(exp(log_af),each=length(c(R))),
    yr=as.factor(rep(seq_len(yrs),each=length(S))))

ggplot(dt)+
geom_line(aes(x=Spawner,y=Recruits,color=yr),linewidth=1.1)+
facet_wrap(~af)+
theme_bw(16)+
scale_colour_viridis_d(end=.9)



#=========================

#beta sensitivity


#Alpha sensitivity scenarios
#log_a<-1.3*2
log_a<-1.3
b<-0.00000581
Smax<-1/b


yrs<-40

bf<-1/c(Smax*.25,Smax*.5,Smax*1.5,Smax*2,Smax*3)

stp<-(bf-b)/(yrs-1)
b_mat<-matrix(NA,ncol=length(bf),nrow=yrs)
b_mat[1,]<-b
for(n in 2:yrs){
   b_mat[n,]<-b_mat[n-1,] +stp
}


S<-seq(0,800000,5000)

R<-matrix(NA,ncol=yrs,nrow=length(S))
R2<-matrix(NA,ncol=yrs,nrow=length(S))
R3<-matrix(NA,ncol=yrs,nrow=length(S))
R4<-matrix(NA,ncol=yrs,nrow=length(S))
R5<-matrix(NA,ncol=yrs,nrow=length(S))



    for(y in seq_len(yrs) ){
      R[,y]<-S*exp(log_a-b_mat[y,1]*S)
      R2[,y]<-S*exp(log_a - b_mat[y,2]*S)
      R3[,y]<-S*exp(log_a - b_mat[y,3]*S)
      R4[,y]<-S*exp(log_a - b_mat[y,4]*S)
      R5[,y]<-S*exp(log_a - b_mat[y,5]*S)
    }

dt<-data.frame(Recruits=c(c(R),c(R2),c(R3),c(R4),c(R5)),
    Spawner=S,
    bf=rep(round(1/bf),each=length(c(R))),
    yr=as.factor(rep(seq_len(yrs),each=length(S))))

ggplot(dt)+
geom_line(aes(x=Spawner,y=Recruits,color=yr),linewidth=1.1)+
facet_wrap(~bf)+
theme_bw(16)+
scale_colour_viridis_d(end=.9)



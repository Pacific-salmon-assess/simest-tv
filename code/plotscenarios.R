#


a<-2
b<-0.0001

at<-rep(c(a-1,a,a-1,a),each=12)[1:40]
at2<-seq(a,a-1,length.out=40)


S<-seq(0,20000,500)

R<-matrix(NA,ncol=40,nrow=length(S))
R2<-matrix(NA,ncol=40,nrow=length(S))


for(i in seq_along(at)){
    R[,i]<-S*exp(at[i]-b*S)
    R2[,i]<-S*exp(at2[i]-b*S)

}

par(mfrow=c(1,2))
matplot(S,R,type="b", ylim=c(0,22000))
matplot(S,R2,type="b", ylim=c(0,22000))

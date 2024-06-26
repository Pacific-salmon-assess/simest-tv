#==========================================
#utility functions
#==========================================




give.n <- function(x){
  return(c(y = median(x,na.rm=T)*1.05, label = length(x))) 
  # experiment with the multiplier to find the perfect position
} 




max.n <- function(x){
  return(c(y = 90*1.05, label = length(x))) 
  # experiment with the multiplier to find the perfect position
} 



sumpair <- function(x,k=2){
  #x<-as.numeric(abs(bhmma$mcmcsummary[grep("^gamma\\[",rownames(bhmma$mcmcsummary)),"Rhat"]-1)>.1)

  y<-matrix(x,ncol=k)
  return(apply(y,1,sum))
  
  # experiment with the multiplier to find the perfect position
} 

mw_tmblfo=function(x){
  w=NA
  for(i in 1:length(x)){w[i]=exp(0.5*x[i])/sum(exp(0.5*x))}
  return(w)
  }


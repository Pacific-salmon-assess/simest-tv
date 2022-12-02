#==========================================
#utility functions
#==========================================




modselecfit <- function(ch,df){
  if(ch=="simple"){
    mod <- ricker_TMB(data=df, priors=1)
    alpha <- rep(mod$alpha,nrow(df))
    smax <- rep(mod$Smax,nrow(df))
    sig <- rep(mod$sig,nrow(df))
    smsy <- rep(mod$Smsy,nrow(df))
    sgen <- rep(sGenCalc(a=mod$alpha,Smsy=mod$Smsy, b=1/mod$Smax)$fit,nrow(df))
    umsy <- rep(mod$umsy,nrow(df))

  }else if(ch=="autocorr"){        
    mod <- ricker_TMB(data=df, AC=TRUE,priors=1)
    alpha <- rep(mod$alpha,nrow(df))
    smax <-rep(mod$Smax,nrow(df))
    sig <- rep(mod$sig,nrow(df))
    smsy<- rep(mod$Smsy,nrow(df))
    sgen <- rep(sGenCalc(a=mod$alpha,Smsy=mod$Smsy, b=1/mod$Smax)$fit,nrow(df))
    umsy <- rep(mod$umsy,nrow(df))

  }else if(ch=="rwa_last"|ch=="rwa_last3"|ch=="rwa_last5"|ch=="rwa"){
    mod <- ricker_rw_TMB(data=df,tv.par='a',priors=1)
    alpha <- mod$alpha
    smax <-rep(mod$Smax,nrow(df))
    sig <- rep(mod$sig,nrow(df))
    smsy <- mod$Smsy
    sgen <- unlist(mapply(sGenCalc,a=alpha,
        Smsy=smsy, 
        b=1/smax))
    umsy <- mod$umsy

  }else if(ch=="rwb_last"|ch=="rwb_last3"|ch=="rwb_last5"|ch=="rwb"){
    mod <- ricker_rw_TMB(data=df, tv.par='b',priors=1)
    alpha <- rep(mod$alpha,nrow(df))
    smax <- mod$Smax
    sig <- rep(mod$sig,nrow(df))
    smsy <- mod$Smsy
    sgen <- unlist(mapply(sGenCalc,a=alpha,
        Smsy=smsy, 
        b=1/smax))
    umsy <- mod$umsy

  }else if(ch=="rwab_last"|ch=="rwab_last3"|"rwab_last5"|"rwab"){
    mod <- ricker_rw_TMB(data=df, tv.par='both',priors=1) 
    alpha <- mod$alpha 
    smax <- mod$Smax
    sig <- rep(mod$sig,nrow(df))
    smsy <- mod$Smsy
    sgen <- unlist(mapply(sGenCalc,a=alpha,
        Smsy=smsy, 
        b=1/smax))
    umsy <- mod$umsy

  }else if(ch=="hmma_last"|ch=="hmma_last3"|ch=="hmma_last5"|ch=="hmma"){
    mod <- ricker_hmm_TMB(data=df, tv.par='a',priors=1)
    alpha <- mod$alpha[mod$regime] 
    smax <-rep(mod$Smax,nrow(df))
    sig <- rep(mod$sigma,nrow(df))
    smsy <- mod$Smsy[mod$regime] 
    sgen <- unlist(mapply(sGenCalc,a=alpha,
        Smsy=smsy, 
        b=1/smax))
    umsy <- mod$umsy

  }else if(ch=="hmmb_last"|ch=="hmmb_last3"|ch=="hmmb_last5"|ch=="hmmb"){
    mod <- ricker_hmm_TMB(data=df, tv.par='b',priors=1)
    alpha <- rep(mod$alpha,nrow(df))
    smax <- mod$Smax[mod$regime] 
    sig <- rep(mod$sigma,nrow(df))
    smsy <- mod$Smsy[mod$regime] 
    sgen <- unlist(mapply(sGenCalc,a=alpha,
        Smsy=smsy, 
        b=1/smax))
    umsy <- mod$umsy

  }else if(ch=="hmm_last"|ch=="hmm_last3"|ch=="hmm_last5"|ch=="hmm"){
    mod <- ricker_hmm_TMB(data=df, tv.par='both',priors=1)
    alpha <- mod$alpha[mod$regime] 
    smax <- mod$Smax[mod$regime]
    sig <- rep(mod$sigma,nrow(df))
    smsy <- mod$Smsy[mod$regime] 
    sgen <- unlist(mapply(sGenCalc,a=alpha,
        Smsy=smsy, 
        b=1/smax))
    umsy <- mod$umsy

  }
  return(list(alpha=alpha,
    smax=smax,
    sig=sig,
    smsy=smsy,
    sgen=sgen,
    umsy=umsy ))

} 




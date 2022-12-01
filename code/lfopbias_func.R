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

  }else if(ch=="rwa_last"|lfoch=="rwa_last3"|lfoch=="rwa_last5"|lfoch=="rwa"){
    mod <- ricker_rw_TMB(data=df,tv.par='a',priors=1)
    alpha <- mod$alpha
    smax <-rep(mod$Smax,nrow(df))
    sig <- rep(mod$sig,nrow(df))
    smsy <- mod$Smsy
    sgen <- unlist(mapply(sGenCalc,a=alpha,
        Smsy=smsy, 
        b=1/smax))
    umsy <- mod$umsy

  }else if(lfoch=="rwb_last"|lfoch=="rwb_last3"|lfoch=="rwb_last5"|lfoch=="rwb"){
    mod <- ricker_rw_TMB(data=df, tv.par='b',priors=1)
    alpha <- rep(mod$alpha,nrow(df))
    smax <- mod$Smax
    sig <- rep(mod$sig,nrow(df))
    smsy <- mod$Smsy
    sgen <- unlist(mapply(sGenCalc,a=alpha,
        Smsy=smsy, 
        b=1/smax))
    umsy <- mod$umsy

  }else if(lfoch=="rwab_last"|lfoch=="rwab_last3"|"rwab_last5"|"rwab"){
    mod <- ricker_rw_TMB(data=df, tv.par='both',priors=1) 
    alpha <- mod$alpha 
    smax <- mod$Smax
    sig <- rep(mod$sig,nrow(df))
    smsy <- mod$Smsy
    sgen <- unlist(mapply(sGenCalc,a=alpha,
        Smsy=smsy, 
        b=1/smax))
    umsy <- mod$umsy

  }else if(lfoch=="hmma_last"|lfoch=="hmma_last3"|lfoch=="hmma_last5"|lfoch=="hmma"){
    mod <- ricker_hmm_TMB(data=df, tv.par='a',priors=1)
    alpha <- mod$alpha[mod$regime] 
    smax <-rep(mod$Smax,nrow(df))
    sig <- rep(mod$sigma,nrow(df))
    smsy <- mod$Smsy[mod$regime] 
    sgen <- unlist(mapply(sGenCalc,a=alpha,
        Smsy=smsy, 
        b=1/smax))
    umsy <- mod$umsy

  }else if(lfoch=="hmmb_last"|lfoch=="hmmb_last3"|lfoch=="hmmb_last5"|lfoch=="hmmb"){
    mod <- ricker_hmm_TMB(data=df, tv.par='b',priors=1)
    alpha <- rep(mod$alpha,nrow(df))
    smax <- mod$Smax[mod$regime] 
    sig <- rep(mod$sigma,nrow(df))
    smsy <- mod$Smsy[mod$regime] 
    sgen <- unlist(mapply(sGenCalc,a=alpha,
        Smsy=smsy, 
        b=1/smax))
    umsy <- mod$umsy

  }else if(lfoch=="hmm_last"|lfoch=="hmm_last3"|lfoch=="hmm_last5"|lfoch=="hmm"){
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




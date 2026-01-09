#' stan_refit function
#'
#' This function refits a stan model (via rstan) to iteratively estimate the out-of-sample likelihood
#' @param sm Compiled stan model, which can be generated with sr_mod function.
#' @param newdata new data to feed to refit models.
#' @param oos What row of the new data should be treated as the out-of-sample test
#' @param regime TRUE or FALSE statement - is this a regime shift model or not
#' @param K Number of potential regime shifts - default 2
#' @param dirichlet_prior k_regime x k_regime matrix. Prior for transition probability matrix, 
#' if NULL (default), prior is set to matrix(1,nrow=k_regime,ncol=k_regime)
#' @param pSmax_mean mean for Smax prior
#' @param pSmax_sig sd for Smax prior
#' @param psig_b sd for sigma_b prior
#' @return returns the model fit
#' @export
#' @examples
#' r=stan_refit(sm=mod3,newdata=df,oos=12)
stan_refit2<- function(mod,newdata,oos,K=2,dirichlet_prior=NULL,pSmax_mean  =NULL,pSmax_sig=NULL,psig_b=NULL){
  #mod = model file name - eg. 'ricker_linear_oos.stan'
  #newdata = data to train model
  #oosdata = data to predict onto
  #regime = TRUE or FALSE for regime shift models (have different data inputs)
  #K = number of potential regimes (2 or 3)
  
  if(is.null(dirichlet_prior)){
    dirichlet_prior<-matrix(1,nrow=K,ncol=K)
  }else if(nrow(dirichlet_prior)!=K |ncol(dirichlet_prior)!=K){
    stop("dirichlet_prior should be a K x K matrix")
  }

  oosdata=newdata[oos,]
  newdata=newdata[-oos,]
  if(is.null(pSmax_mean)){pSmax_mean <- max(newdata$S)/2} 
  if(is.null(pSmax_sig)){pSmax_sig <- max(newdata$S)/2} 
  if(is.null(psig_b)){psig_b <- max(newdata$S)/2} 
  
  dl=list(
    by=newdata$by,
    N=nrow(newdata),
    L=max(newdata$by)-min(newdata$by)+1,
    ii=newdata$by-min(newdata$by)+1,
    R_S =newdata$logRS,
    S=newdata$S,
    y_oos=oosdata$logRS,
    x_oos=oosdata$S,
    K=2,
    alpha_dirichlet= dirichlet_prior,
    pSmax_mean=pSmax_mean,
    pSmax_sig=pSmax_sig,
    psig_b=psig_b
    
  )
  
  r = mod$sample(data=dl,
                  seed=123,
                  chains=6,
                  iter_warmup=200,
                  iter_sampling=500,
                  refresh=0,
                  adapt_delta=0.95,
                  max_treedepth=15)
  return(r)
}


#' stan_lfo_cv function
#'
#' This function implements exact Leave-Future-Out Cross-Validation (LFO-CV) (see https://arxiv.org/abs/1902.06281), 
#' which will iteratively estimate model likelihood in an out-of-sample context using past data as a training set.
#' @param mod Compiled stan model, which can be generated with sr_mod function.
#' @param type The type of stock-recruitment model that is being fit. Either 'static' for the time-invariant model (model 1: static Ricker), 
#' 'tv' for time-varying models (including autocorrelated residual model), or 'regime' for regime shift models
#' @param df Dataframe to use
#' @param L Minimum dataset to retain, default is to 10 years of data
#' @param K Number of potential regime shifts - default 2
#' @param dirichlet_prior k_regime x k_regime matrix. Prior for transition probability matrix, 
#' if NULL prior is set to matrix(1,nrow=k_regime,ncol=k_regime)
#' @param pSmax_mean mean for Smax prior
#' @param pSmax_sig sd for Smax prior
#' @param psig_b sd for sigma_b prior
#' @return returns the pointwise out-of-sample log likelihoods. For 'tv' or 'regime' models returns a vector that uses either 1. the last year's time-varying parameters, 
#' 2. an average of the last 3 years of the time-varying parameters, or 5. an average of the last 5 years of the time-varying parameters. For 'regime' models also includes
#' these estimates that are also probability-weighted (ie. prob regime x regime parameter) in addition to the parameters for the most likely state 1, 3 or 5 years back.
#' @export
#' @examples
#' r=stan_refit(sm=mod3,newdata=df,oos=12)
stan_lfo_cv=function(mod,type=c('static','tv','regime'),df,L=10,K=2,dirichlet_prior=NULL,
  pSmax_mean  =NULL,pSmax_sig=NULL,psig_b=NULL){
  #mod = model to fit (model name for cmdstanr)
  #tv = 0 for static model; 1 for time-varying (for calculating elpds)
  #df = full data frame
  #L = starting point for LFO-CV (default 10)
  # K = number of regimes
  if(is.null(pSmax_mean)){pSmax_mean <- max(df$S)/2} 
  if(is.null(pSmax_sig)){pSmax_sig <- max(df$S)/2} 
  if(is.null(psig_b)){psig_b <- max(df$S)/2} 
  
  loglik_exact <- matrix(nrow = 3000, ncol = length(df$by)) #loglik for static model
  loglik_exact_1b <- matrix(nrow = 3000, ncol = length(df$by)) #loglik for 1-year back estimates of productivity/capacity
  loglik_exact_3b <- matrix(nrow = 3000, ncol = length(df$by)) #loglik for average of last 3-years of productivity/capacity
  loglik_exact_5b <- matrix(nrow = 3000, ncol = length(df$by)) #loglik for average of last 5-years of productivity/capacity
  if(type=='regime'){
    loglik_exact_1bw <- matrix(nrow = 3000, ncol = length(df$by)) #loglik for 1-year back estimates of productivity/capacity
    loglik_exact_3bw <- matrix(nrow = 3000, ncol = length(df$by)) #loglik for average of last 3-years of productivity/capacity
    loglik_exact_5bw <- matrix(nrow = 3000, ncol = length(df$by)) #loglik for average of last 5-years of productivity/capacity
  }
  for (i in L:(nrow(df) - 1)){
    past <- 1:i
    oos <- i + 1
    df_past <- df[past, , drop = FALSE]
    df_oos <- df[c(past, oos), , drop = FALSE]
    if(type=='static'){
      fit_past<- stan_refit2(mod=mod, newdata=df_oos, oos=i+1,pSmax_mean = pSmax_mean, 
                            pSmax_sig = pSmax_sig,psig_b = psig_b)
      ll=as.data.frame(fit_past$draws(variables=c('log_lik_oos'),format='draws_matrix'))
      loglik_exact[,i+1]<- ll$log_lik_oos
    }
    if(type=='tv'){
      fit_past<- stan_refit2(mod=mod,newdata=df_oos,oos=i+1,pSmax_mean = pSmax_mean, 
                            pSmax_sig = pSmax_sig,psig_b = psig_b)
      ll=as.data.frame(fit_past$draws(variables=c('log_lik_oos_1b','log_lik_oos_3b','log_lik_oos_5b'),format='draws_matrix'))
      
      loglik_exact_1b[, i + 1] <-ll$log_lik_oos_1b
      loglik_exact_3b[, i + 1] <-ll$log_lik_oos_3b
      loglik_exact_5b[, i + 1] <-ll$log_lik_oos_5b
    }
    if(type=='regime'){
      fit_past<- stan_refit2(mod=mod,newdata=df_oos,oos=i+1, K=K, dirichlet_prior=dirichlet_prior,pSmax_mean = pSmax_mean, 
                            pSmax_sig = pSmax_sig,psig_b = psig_b)
      ll=as.data.frame(fit_past$draws(variables=c('log_lik_oos_1b','log_lik_oos_3b','log_lik_oos_5b'),format='draws_matrix'))
      loglik_exact_1b[, i + 1] <- ll$log_lik_oos_1b
      loglik_exact_3b[, i + 1] <- ll$log_lik_oos_3b
      loglik_exact_5b[, i + 1] <- ll$log_lik_oos_5b
    }
  }
  
  if(type=='static'){
    exact_elpds<- apply(loglik_exact, 2, log_mean_exp); exact_elpds=exact_elpds[-(1:L)]
    r=exact_elpds
  }
  if(type=='tv'){
    exact_elpds_1b <- apply(loglik_exact_1b, 2, log_mean_exp); exact_elpds_1b=exact_elpds_1b[-(1:L)]
    exact_elpds_3b <- apply(loglik_exact_3b, 2, log_mean_exp); exact_elpds_3b=exact_elpds_3b[-(1:L)]
    exact_elpds_5b <- apply(loglik_exact_5b, 2, log_mean_exp); exact_elpds_5b=exact_elpds_5b[-(1:L)]
    
    r=rbind(exact_elpds_1b,exact_elpds_3b,exact_elpds_5b)
  }
  if(type=='regime'){
    exact_elpds_1b <- apply(loglik_exact_1b, 2, log_mean_exp); exact_elpds_1b=exact_elpds_1b[-(1:L)]
    exact_elpds_3b <- apply(loglik_exact_3b, 2, log_mean_exp); exact_elpds_3b=exact_elpds_3b[-(1:L)]
    exact_elpds_5b <- apply(loglik_exact_5b, 2, log_mean_exp); exact_elpds_5b=exact_elpds_5b[-(1:L)]
   
    r=rbind(exact_elpds_1b,exact_elpds_3b,exact_elpds_5b)
  }
  return(r)
}

#' model_weights function
#'
#' This function implements estimates Pseudo Bayesian Model Averaging (Pseudo-BMA) weights based on Yao et al. 2018  (eq. 8, see http://www.stat.columbia.edu/~gelman/research/published/stacking_paper_discussion_rejoinder.pdf), or AIC/BIC weights.
#' @param x a dataframe of pointwise out-of-sample loglikelihoods, where each row represents model likelihood predictions - these can be estimated with the stan_lfo_cv function
#' @param form option for either pseudo-Bayesian model averaging weight (PBMA) or AIC/BIC weights (AIC)
#' @param type option on what observations to include: full = entire out of sample years, d90 = excluding the 10% of years with the lowest likelihood among all models, d80 = excluding the 20% of years with the lowest likelihood among all models.  
#' @return returns the relative weights for each of the included models
#' @export
#' @examples
#' model_weights(rbind(ll1,ll2))

model_weights<- function(x,form=c('stan','TMB'),type=c('full','d90','d80')){
  if(form=='stan'){
    if(type=='full'){
      elpd_1=apply(x,1,sum) #
      elpd_2=NA
      w=NA
      se_elpd=NA
      for(i in 1:nrow(x)){
        se_elpd[i]=sqrt(sum((x[i,]-(sum(x[i,])/ncol(x)))^2))
        elpd_2[i]=exp(elpd_1[i]-0.5*se_elpd[i])
      }
    }
    if(type=='d90'){
      x2=x[,apply(x,2,sum)>=quantile(apply(x,2,sum),0.1)]
      elpd_1=apply(x2,1,sum) #
      elpd_2=NA
      w=NA
      se_elpd=NA
      for(i in 1:nrow(x2)){
        se_elpd[i]=sqrt(sum((x2[i,]-(sum(x2[i,])/ncol(x2)))^2))
        elpd_2[i]=exp(elpd_1[i]-0.5*se_elpd[i])
      }
      
    }
    if(type=='d80'){
      x2=x[,apply(x,2,sum)>=quantile(apply(x,2,sum),0.2)]
      elpd_1=apply(x2,1,sum) #
      elpd_2=NA
      w=NA
      se_elpd=NA
      for(i in 1:nrow(x2)){
        se_elpd[i]=sqrt(sum((x2[i,]-(sum(x2[i,])/ncol(x2)))^2))
        elpd_2[i]=exp(elpd_1[i]-0.5*se_elpd[i])
      }
      
    }
    for(i in 1:nrow(x)){
      w[i]=elpd_2[i]/sum(elpd_2) 
    }
  }
  if(form=='TMB'){
    w=NA
    for(i in 1:length(x)){w[i]=exp(-0.5*x[i])/sum(exp(-0.5*x))}
  }
  return(round(w,digits=3))
}

#' AIC for stan models
#'
#' This function estimates either AIC/BIC scores based on pointwise log likelihood - taking the average across the posterior.
#' @param x a list of pointwise out-of-sample loglikelihoods from different models, where each row represents model likelihood predictions
#' @param form option for either AIC pr BIC
#' @param type option on what observations to include: full = entire out of sample years, d90 = excluding the 10% of years with the lowest likelihood among all models, d80 = excluding the 20% of years with the lowest likelihood among all models.  
#' @return returns an AIC/BIC score for each model
#' @export
#' 
stan_aic<- function(x,form=c('aic','bic'),type=c('full','d90','d80'),k){
  LL=NA;AIC=NA;BIC=NA
  elpd_1=matrix(nrow=length(x),ncol=ncol(x[[1]]))
    if(type=='full'){
      for(i in 1:length(x)){
        elpd_1[i,]=apply(x[[i]],2,log_mean_exp) #
      }
    LL=apply(elpd_1,1,sum)
    AIC=-2*LL+2*k+(2*k^2+2*k)/(ncol(x[[1]])-k-1)
    BIC=-2*LL+k*log(ncol(x[[1]]))
      
    dAIC=AIC-min(AIC)
    dBIC=BIC-min(BIC)
    w_aic=NA
    w_bic=NA
    for(i in 1:length(x)){w_aic[i]=exp(-0.5*dAIC[i])/sum(exp(-0.5*dBIC))}
    for(i in 1:length(x)){w_bic[i]=exp(-0.5*dBIC[i])/sum(exp(-0.5*dBIC))}
    }
  if(type=='d90'){
    for(i in 1:length(x)){
      elpd_1[i,]=apply(x[[i]],2,log_mean_exp) #
    }
    elpd_1=elpd_1[,apply(elpd_1,2,sum)>=quantile(apply(elpd_1,2,sum),0.1)]
  
    LL=apply(elpd_1,1,sum)
    AIC=-2*LL+2*k+(2*k^2+2*k)/(ncol(x[[1]])-k-1)
    BIC=-2*LL+k*log(ncol(x[[1]]))
    
    dAIC=AIC-min(AIC)
    dBIC=BIC-min(BIC)
    w_aic=NA
    w_bic=NA
    for(i in 1:length(x)){w_aic[i]=exp(-0.5*dAIC[i])/sum(exp(-0.5*dBIC))}
    for(i in 1:length(x)){w_bic[i]=exp(-0.5*dBIC[i])/sum(exp(-0.5*dBIC))}
  }
  if(type=='d80'){
    for(i in 1:length(x)){
      elpd_1[i,]=apply(x[[i]],2,log_mean_exp) #
    }
    elpd_1=elpd_1[,apply(elpd_1,2,sum)>=quantile(apply(elpd_1,2,sum),0.2)]
    
    LL=apply(elpd_1,1,sum)
    AIC=-2*LL+2*k+(2*k^2+2*k)/(ncol(x[[1]])-k-1)
    BIC=-2*LL+k*log(ncol(x[[1]]))
    
    dAIC=AIC-min(AIC)
    dBIC=BIC-min(BIC)
    w_aic=NA
    w_bic=NA
    for(i in 1:length(x)){w_aic[i]=exp(-0.5*dAIC[i])/sum(exp(-0.5*dBIC))}
    for(i in 1:length(x)){w_bic[i]=exp(-0.5*dBIC[i])/sum(exp(-0.5*dBIC))}
  }
  if(form=='aic'){
    return(w_aic)
  }
  if(form=='bic'){
    return(w_bic)
  }
}


#For rstan:::
stan_refit<- function(mod,newdata,oos,K=2,dirichlet_prior=NULL){
  #mod = model file name - eg. 'ricker_linear_oos.stan'
  #newdata = data to train model
  #oosdata = data to predict onto
  #regime = TRUE or FALSE for regime shift models (have different data inputs)
  #K = number of potential regimes (2 or 3)
  
  if(is.null(dirichlet_prior)){
    dirichlet_prior<-matrix(1,nrow=K,ncol=K)
  }else if(nrow(dirichlet_prior)!=K |ncol(dirichlet_prior)!=K){
    stop("dirichlet_prior should be a K x K matrix")
  }
  
  oosdata=newdata[oos,]
  newdata=newdata[-oos,]
  
  dl=list(
    by=newdata$by,
    N=nrow(newdata),
    L=max(newdata$by)-min(newdata$by)+1,
    ii=newdata$by-min(newdata$by)+1,
    R_S =newdata$logRS,
    S=newdata$S,
    y_oos=oosdata$logRS,
    x_oos=oosdata$S,
    K=2,
    alpha_dirichlet= dirichlet_prior
  )
  
  r = rstan::sampling(mod,
                      data=dl,
                      chains=1,
                      warmup=200,
                      iter=700)
  return(r)
}


#' stan_lfo_cv function

stan_lfo_cv2=function(mod,type=c('static','tv','regime'),df,L=10,K=2,dirichlet_prior=NULL){
  #mod = model to fit (model name for cmdstanr)
  #tv = 0 for static model; 1 for time-varying (for calculating elpds)
  #df = full data frame
  #L = starting point for LFO-CV (default 10)
  # K = number of regimes
  
  loglik_exact <- matrix(nrow = 3000, ncol = length(df$by)) #loglik for static model
  loglik_exact_1b <- matrix(nrow = 3000, ncol = length(df$by)) #loglik for 1-year back estimates of productivity/capacity
  loglik_exact_3b <- matrix(nrow = 3000, ncol = length(df$by)) #loglik for average of last 3-years of productivity/capacity
  loglik_exact_5b <- matrix(nrow = 3000, ncol = length(df$by)) #loglik for average of last 5-years of productivity/capacity
  if(type=='regime'){
    loglik_exact_1bw <- matrix(nrow = 3000, ncol = length(df$by)) #loglik for 1-year back estimates of productivity/capacity
    loglik_exact_3bw <- matrix(nrow = 3000, ncol = length(df$by)) #loglik for average of last 3-years of productivity/capacity
    loglik_exact_5bw <- matrix(nrow = 3000, ncol = length(df$by)) #loglik for average of last 5-years of productivity/capacity
  }
  for (i in L:(nrow(df) - 1)){
    past <- 1:i
    oos <- i + 1
    df_past <- df[past, , drop = FALSE]
    df_oos <- df[c(past, oos), , drop = FALSE]
    if(type=='static'){
      fit_past<- stan_refit2(mod=mod, newdata=df_oos, oos=i+1)
      ll=as.data.frame(rstan::extract(fit_past,pars=c('log_lik_oos')))
      loglik_exact[,i+1]<- ll$log_lik_oos
    }
    if(type=='tv'){
      fit_past<- stan_refit2(mod=mod,newdata=df_oos,oos=i+1)
      ll=as.data.frame(rstan::extract(fit_past,pars=c('log_lik_oos_1b','log_lik_oos_3b','log_lik_oos_5b')))
      
      loglik_exact_1b[, i + 1] <-ll$log_lik_oos_1b
      loglik_exact_3b[, i + 1] <-ll$log_lik_oos_3b
      loglik_exact_5b[, i + 1] <-ll$log_lik_oos_5b
    }
    if(type=='regime'){
      fit_past<- stan_refit2(mod=mod,newdata=df_oos,oos=i+1, K=K, dirichlet_prior=dirichlet_prior)
      ll=as.data.frame(rstan::extract(fit_past,pars=c('log_lik_oos_1b','log_lik_oos_3b','log_lik_oos_5b')))
      loglik_exact_1b[, i + 1] <- ll$log_lik_oos_1b
      loglik_exact_3b[, i + 1] <- ll$log_lik_oos_3b
      loglik_exact_5b[, i + 1] <- ll$log_lik_oos_5b
    }
  }
  
  if(type=='static'){
    exact_elpds<- apply(loglik_exact, 2, log_mean_exp); exact_elpds=exact_elpds[-(1:L)]
    r=exact_elpds
  }
  if(type=='tv'){
    exact_elpds_1b <- apply(loglik_exact_1b, 2, log_mean_exp); exact_elpds_1b=exact_elpds_1b[-(1:L)]
    exact_elpds_3b <- apply(loglik_exact_3b, 2, log_mean_exp); exact_elpds_3b=exact_elpds_3b[-(1:L)]
    exact_elpds_5b <- apply(loglik_exact_5b, 2, log_mean_exp); exact_elpds_5b=exact_elpds_5b[-(1:L)]
    
    r=rbind(exact_elpds_1b,exact_elpds_3b,exact_elpds_5b)
  }
  if(type=='regime'){
    exact_elpds_1b <- apply(loglik_exact_1b, 2, log_mean_exp); exact_elpds_1b=exact_elpds_1b[-(1:L)]
    exact_elpds_3b <- apply(loglik_exact_3b, 2, log_mean_exp); exact_elpds_3b=exact_elpds_3b[-(1:L)]
    exact_elpds_5b <- apply(loglik_exact_5b, 2, log_mean_exp); exact_elpds_5b=exact_elpds_5b[-(1:L)]
    
    r=rbind(exact_elpds_1b,exact_elpds_3b,exact_elpds_5b)
  }
  return(r)
}


#' Simple Ricker model estimated with stan
#'
#' @param data A list or data frame containing complete vectors for: brood year (by), Spawners (s) and log(Recruits/Spawners) (logRS) time series.
#' @param ac Logical. Are residuals autocorrelated? Default is FALSE
#' @param smax_priors User-specified smax priors, sampled from a normal distribution, should be formatted as: c(mean, sd). Defaults to mean = half maximum observed spawners, sd = max. observed spawners (from time series). 
#' @param control output of stancontrol
#' @param warmup To be passed to rstan::sampling. A positive integer specifying the number of warmup (aka burnin) iterations per
#'  chain. The default is 200.
#' @param chains To be passed to rstan::sampling. A positive integer specifying the number of Markov chains. The default is 6.
#' @param iter To be passed to rstan::sampling. A positive integer specifying the number of iterations for each chain 
#' (including warmup). The default is 1000.
#' 
#' @param ... Anything else that would be passed to rstan::sampling
#' 
#' 
#' @returns a list containing the 
#' * data - input dataseries
#' * logalpha - median estimates for the log(alpha) parameter 
#' * beta - median estimates for the beta parameter 
#' * Smax - median estimates for the Smax parameter
#' * Smsy - median estimates for the Smsy parameter
#' * Umsy - median estimates for the Umsy parameter
#' * sigma - median estimates for the sigma parameter   
#' * rho - median estimates for the rho parameter 
#' * residuals - median estimates for each epsilon (residual) parameter      
#' * summary - summary of stan model fit
#' * full_posterior - all posterior samples of the above parameters 
#' * stanfit - the stanfit model object
#' 
#' 
#' @importFrom rstan stan extract summary sampling
#' 
#' 
#' @export
#' 
#' @examples
#' data(harck)
#' ricker_stan(data=harck)
#' 
ricker_stan <- function(data,  ac=FALSE, smax_priors=NULL,mod=NULL, control = stancontrol(adapt_delta=0.99), warmup=300, chains = 6, iter = 1000,...) {

    if(is.null(mod)==T){
    sm=sr_mod(type='static',ac=ac,par='n')
  }else{sm=mod}

  if(is.null(smax_priors)==TRUE){
    datm = list(N=nrow(data),
                L=max(data$by)-min(data$by)+1,
                ii=data$by-min(data$by)+1,
                R_S =data$logRS,
                S=data$S,
                pSmax_mean=max(data$S)/2,
                pSmax_sig=max(data$S)*2)
  }else{
      datm = list(N=nrow(data),
                  L=max(data$by)-min(data$by)+1,
                  ii=data$by-min(data$by)+1,
                  R_S =data$logRS,
                  S=data$S,
                  pSmax_mean=smax_priors[1],
                  pSmax_sig=smax_priors[2])
    }

  fit<-rstan::sampling(sm, data=datm,
                      control = control, warmup = warmup, 
                      chains = chains, iter = iter,verbose=FALSE)
 if(ac==FALSE){ 
  mc <- rstan::extract(fit,pars=c('logalpha','beta','Smax','Smsy','Umsy','sigma','mu','epsilon'),permuted=T) 
  mc2=as.data.frame(do.call(cbind,mc))
  colnames(mc2)=c('logalpha','beta','Smax','Smsy','Umsy','sigma',paste('mu[',seq(1:datm$N),']',sep=''),paste('epsilon[',seq(1:datm$N),']',sep=''))
  } 
  if(ac==TRUE){ 
  mc <- rstan::extract(fit,pars=c('logalpha','beta','Smax','Smsy','Umsy','sigma_AR','mu','epsilon','rho'),permuted=T)
  mc2=as.data.frame(do.call(cbind,mc))
  colnames(mc2)=c('logalpha','beta','Smax','Smsy','Umsy','sigma',paste('mu[',seq(1:datm$N),']',sep=''),paste('epsilon[',seq(1:datm$N),']',sep=''),'rho')
  } 
  
  aa <- rstan::summary(fit)
  if(full_posterior==FALSE){
     if(ac==F){
       ans<-list(data=data,
                 logalpha=c(aa$summary["logalpha","50%"]),
                 beta=c(aa$summary["beta","50%"]),
                 Smax=c(aa$summary["Smax","50%"]),
                 Smsy=c(aa$summary["Smsy","50%"]),
                 Umsy=c(aa$summary["Umsy","50%"]),
                 sigma=c(aa$summary["sigma","50%"]),
                 rho=NA,
                 residuals=apply(rstan::extract(fit,pars=c('epsilon'),permuted=T),2,median),
                 summary=aa$summary,
                 full_posterior=mc2,
                 stanfit=fit)
     }
    if(ac==TRUE){
      ans<-list(data=data,
                logalpha=c(aa$summary["logalpha","50%"]),
                beta=c(aa$summary["beta","50%"]),
                Smax=c(aa$summary["Smax","50%"]),
                Smsy=c(aa$summary["Smsy","50%"]),
                Umsy=c(aa$summary["Umsy","50%"]),
                sigma=c(aa$summary["sigma_AR","50%"]),
                rho=c(aa$summary["rho","50%"]),
                residuals=apply(rstan::extract(fit,pars=c('epsilon'),permuted=T),2,median),
                summary=aa$summary,
                full_posterior=mc2,
                stanfit=fit)
    }
    
  }
  return(ans)
  if(any(aa$rhat)>=1.05){print='Warning, some R_hat values are over the threshold of 1.05 - check parameter summary'}
  
}








#' random walks Ricker model estimated with stan
#'
#' @param data A list or data frame containing complete vectors for: brood year (by), Spawners (s) and log(Recruits/Spawners) (logRS) time series. Use sr_format for the correct column names. 
#' @param tv.par Which parameter should vary? Either productivity (intercept, a), capacity (slope, b) or both parameters (both)
#' @param smax_priors User-specified smax priors, sampled from a normal distribution, should be formatted as: c(mean, sd)
#' @param control output of stancontrol
#' @param warmup To be passed to rstan::sampling. A positive integer specifying the number of warmup (aka burnin) iterations per
#'  chain. The default is 200.
#' @param chains To be passed to rstan::sampling. A positive integer specifying the number of Markov chains. The default is 6.
#' @param iter To be passed to rstan::sampling. A positive integer specifying the number of iterations for each chain 
#' (including warmup). The default is 1000.
#' @param ... Anything else that would be passed to rstan::sampling
#' 
#' @returns a list containing the 
#' * data - input dataseries
#' * logalpha - median estimates for the log(alpha) parameter 
#' * beta - median estimates for the beta parameter 
#' * Smax - median estimates for the Smax parameter
#' * Smsy - median estimates for the Smsy parameter
#' * Umsy - median estimates for the Umsy parameter
#' * sigma - median estimates for the sigma parameter   
#' * sigma_a - median estimates for the sigma_a parameter (variance attributed to random walk in log(alpha)) 
#' * sigma_b - median estimates for the sigma_b parameter (variance attributed to random walk in Smax) 
#' * residuals - median estimates for each epsilon (residual) parameter      
#' * summary - summary of stan model fit
#' * full_posterior - all posterior samples of the above parameters 
#' * stanfit - the stanfit model object
#' 
#' 
#' @importFrom rstan stan extract summary sampling
#' 
#' 
#' @export
#' 
#' @examples
#' data(harck)
#' ricker_rw_stan(data=harck)
#' 
ricker_rw_stan <- function(data, tv.par=c('a','b','both'),smax_priors=NULL,control = stancontrol(adapt_delta=0.99), mod=NULL,
  warmup=300,  chains = 6, iter = 1000,...) {
  tv.par=match.arg(tv.par,choices=c('a','b','both'))

  if(is.null(mod)==T){
    sm=sr_mod(type='rw',par=tv.par)
  }else{sm=mod}
  
  if(is.null(smax_priors)==TRUE){
    datm = list(N=nrow(data),
                L=max(data$by)-min(data$by)+1,
                ii=data$by-min(data$by)+1,
                R_S =data$logRS,
                S=data$S,
                pSmax_mean=max(data$S)/2,
                pSmax_sig=max(data$S)*2)
    
  }else{
    datm = list(N=nrow(data),
                L=max(data$by)-min(data$by)+1,
                ii=data$by-min(data$by)+1,
                R_S =data$logRS,
                S=data$S,
                pSmax_mean=smax_priors[1],
                pSmax_sig=smax_priors[2])
  }
  
  fit <- rstan::sampling(sm, data = datm,
                         control = control, warmup = warmup, chains = chains, iter = iter,verbose=FALSE)
  

  aa <- rstan::summary(fit)
  
    if(tv.par=='a'){
      mc <- rstan::extract(fit,pars=c('logalpha','beta','Smax','Smsy','sigma','sigma_a','mu','epsilon'),permuted=T) 
      mc2=as.data.frame(do.call(cbind,mc))
      colnames(mc2)=c(paste('logalpha[',seq(1:datm$L),']',sep=''),'beta','Smax',paste('Smsy[',seq(1:datm$L),']',sep=''),'sigma','sigma_a',paste('mu[',seq(1:datm$N),']',sep=''),paste('epsilon[',seq(1:datm$N),']',sep=''))
      
      ans<-list(data=data,
                logalpha=c(aa$summary[paste('logalpha[',seq(1:datm$L),']',sep=''),"50%"]),
                beta=c(aa$summary["beta","50%"]),
                Smax=c(aa$summary["Smax","50%"]),
                Smsy=c(aa$summary[paste('Smsy[',seq(1:datm$L),']',sep=''),"50%"]),
                Umsy=c(aa$summary[paste('Umsy[',seq(1:datm$L),']',sep=''),"50%"]),
                sigma=c(aa$summary["sigma","50%"]),
                sigma_a=c(aa$summary["sigma_a","50%"]),
                residuals=apply(rstan::extract(fit,pars=c('epsilon'),permuted=T),2,median),
                summary=aa$summary,
                full_posterior=mc2,
                stanfit=fit)
    
    }
    if(tv.par=='b'){
      mc <- rstan::extract(fit,pars=c('logalpha','beta','Smax','Smsy','Umsy','sigma','sigma_b','mu','epsilon'),permuted=T) 
      mc2=as.data.frame(do.call(cbind,mc))
      colnames(mc2)=c('logalpha',paste('beta[',seq(1:datm$L),']',sep=''),paste('Smax[',seq(1:datm$L),']',sep=''),paste('Smsy[',seq(1:datm$L),']',sep=''),'Umsy','sigma','sigma_b',paste('mu[',seq(1:datm$N),']',sep=''),paste('epsilon[',seq(1:datm$N),']',sep=''))
      
      ans<-list(data=data,
                logalpha=c(aa$summary["logalpha","50%"]),
                beta=c(aa$summary[paste('beta[',seq(1:datm$L),']',sep=''),"50%"]),
                Smax=c(aa$summary[paste('Smax[',seq(1:datm$L),']',sep=''),"50%"]),
                Smsy=c(aa$summary[paste('Smsy[',seq(1:datm$L),']',sep=''),"50%"]),
                Umsy=c(aa$summary["Umsy","50%"]),
                sigma=c(aa$summary["sigma","50%"]),
                sigma_b=c(aa$summary["sigma_b","50%"]),
                residuals=apply(rstan::extract(fit,pars=c('epsilon'),permuted=T),2,median),
                summary=aa$summary,
                full_posterior=mc2,
                stanfit=fit)
     }
    if(tv.par=='both'){
      mc <- rstan::extract(fit,pars=c('logalpha','beta','Smax','Smsy','Umsy','sigma','sigma_a','sigma_b','mu','epsilon'),permuted=T) 
      mc2=as.data.frame(do.call(cbind,mc))
      colnames(mc2)=c(paste('logalpha[',seq(1:datm$L),']',sep=''),paste('beta[',seq(1:datm$L),']',sep=''),paste('Smax[',seq(1:datm$L),']',sep=''),paste('Smsy[',seq(1:datm$L),']',sep=''),paste('Umsy[',seq(1:datm$L),']',sep=''),'sigma','sigma_a','sigma_b',paste('mu[',seq(1:datm$N),']',sep=''),paste('epsilon[',seq(1:datm$N),']',sep=''))
      
      
      ans<-list(data=data,
           logalpha=c(aa$summary[paste('logalpha[',seq(1:datm$L),']',sep=''),"50%"]),
           beta=c(aa$summary[paste('beta[',seq(1:datm$L),']',sep=''),"50%"]),
           Smax=c(aa$summary[paste('Smax[',seq(1:datm$L),']',sep=''),"50%"]),
           Smsy=c(aa$summary[paste('Smsy[',seq(1:datm$L),']',sep=''),"50%"]),
           Umsy=c(aa$summary[paste('Umsy[',seq(1:datm$L),']',sep=''),"50%"]),
           sigma=c(aa$summary["sigma","50%"]),
           sigma_a=c(aa$summary["sigma_a","50%"]),
           sigma_b=c(aa$summary["sigma_b","50%"]),
           residuals=apply(rstan::extract(fit,pars=c('epsilon'),permuted=T),2,median),
           summary=aa$summary,
           full_posterior=mc2,
           stanfit=fit)
      
    }
    
  return(ans)
  if(any(aa$rhat)>=1.05){print='Warning, some R_hat values are over the threshold of 1.05 - check parameter summary'}
  
}



#' Hidden markov (regime shift) Ricker model estimated with stan
#'
#' @param data A list or data frame containing Spawners (S) and log(Recruits/Spawners) (logRS) time series. 
#' @param par Which parameter should vary? Either productivity (intercept, a), capacity (slope, b) or both parameters
#' @param k_regime number of regimes in the model, default is 2
#' @param smax_priors Options for custom Smax (capacity) priors - 2 values for the mean and the scale (variance). Defaults to mean = half maximum observed spawners, scale = max. observed spawners (from data) 
#' @param control output of stancontrol
#' @param warmup To be passed to rstan::sampling. A positive integer specifying the number of warmup (aka burnin) iterations per
#'  chain. The default is 200.
#' @param chains To be passed to rstan::sampling. A positive integer specifying the number of Markov chains. The default is 6.
#' @param iter To be passed to rstan::sampling. A positive integer specifying the number of iterations for each chain 
#' (including warmup). The default is 1000.
#' @param sm_ext Default is null, external stan hmm model. Implemented to speed up simulation evaluation, use with caution.
#' @param dirichlet_prior k_regime x k_regime matrix. Prior for transition probability matrix, 
#' if NULL prior is set to matrix(1,nrow=k_regime,ncol=k_regime)
#' @param ... Anything else that would be passed to rstan::sampling
#'  
#'
#' @returns a list containing the 
#' * alpha - median estimates for the alpha parameter vector
#' * beta - median estimates for the beta parameter 
#' * sigobs - median estimates for the observation error sigma         
#' * stanfit - a stanfit model object
#' * mcmcsummary - summary over kept samples
#' * c_mcmcsummary - chain specific summary 
#' * list of samples 
#' 
#' 
#' @importFrom rstan stan extract summary 
#'
#' 
#' @export
#' 
#' @examples
#' data(harck)
#' ricker_hmm_stan(data=harck)
#' 
ricker_hmm_stan <- function(data, par=c('a','b','both'), k_regime=2, smax_priors=NULL,smax_dists=c('normal','lognormal','cauchy'), dirichlet_stasis_prior=2, full_posterior=FALSE,
  control = stancontrol(), warmup=300,  chains = 6, iter = 1000, mod=NULL,...) {

  par=match.arg(tv.par,choices=c('a','b','both'))
  dirichlet_prior<-matrix(c(dirichlet_stasis_prior,1,1,dirichlet_stasis_prior),nrow=k_regime,ncol=k_regime)
  
  if(is.null(mod)){
    sm=samEst::sr_mod(type='hmm',par=par)
  }else{
    sm <-mod
  }
  
  if(is.null(smax_priors)==TRUE){
    datm = list(N=nrow(data),
                L=max(data$by)-min(data$by)+1,
                ii=data$by-min(data$by)+1,
                R_S =data$logRS,
                S=data$S,
                K=k_regime,
                alpha_dirichlet=dirichlet_prior,
                pSmax_mean=max(data$S)/2,
                pSmax_sig=max(data$S)*2)
  }else{
    datm = list(N=nrow(data),
                L=max(data$by)-min(data$by)+1,
                ii=data$by-min(data$by)+1,
                R_S =data$logRS,
                S=data$S,
                K=k_regime,
                alpha_dirichlet=dirichlet_prior,
                pSmax_mean=smax_priors[1],
                pSmax_sig=smax_priors[2])
  }
  
  #if(is.null(sm_ext)){
  #  sm <- sr_mod(type='hmm',ac=FALSE,par=par,loglik=FALSE, modelcode=TRUE)
  #}else{
  #  sm <- sm_ext
  #}
 
  fit <- rstan::sampling(sm, data = datm,
                         control = control, warmup = warmup, chains = chains, iter = iter,verbose=FALSE)
  

  mc <- rstan::extract(fit, 
          inc_warmup=FALSE, permuted=FALSE)
    

  aa <- rstan::summary(fit)
  parts <-stan_regime_rps(m=fit,par=par)
 
  ans<-list(
   logalpha=ifelse(par=="a"|par=="both",list(aa$summary[grep("log_a\\[",row.names(aa$summary)),"50%"]),aa$summary["log_a","50%"])[[1]],
   beta=ifelse(par=="b"|par=="both",list(aa$summary[grep("b\\[",row.names(aa$summary)),"50%"]),aa$summary["b","50%"])[[1]],
   Smax=ifelse(par=="b"|par=="both",list(aa$summary[grep("Smax\\[",row.names(aa$summary)),"50%"]),aa$summary["Smax","50%"])[[1]],
   Smsy=aa$summary[grep("Smsy\\[",row.names(aa$summary)),"50%"],
   Umsy=ifelse(par=="a"|par=="both",list(aa$summary[grep("Umsy\\[",row.names(aa$summary)),"50%"]),aa$summary["Umsy","50%"])[[1]],
   logalpha_regime=ifelse(rep(par=="a"|par=="both",nrow(data)),parts$log_a_t,NA),
   logalpha_wgt=ifelse(rep(par=="a"|par=="both",nrow(data)),parts$log_a_wt,NA),
   Smax_regime=ifelse(rep(par=="b"|par=="both",nrow(data)),parts$Smax_t,NA),
   Smax_wgt=ifelse(rep(par=="b"|par=="both",nrow(data)),parts$Smax_wt,NA),
   Smsy_regime=parts$Smsy_t,
   Smsy_wgt=parts$Smsy_wt,
   Umsy_regime=ifelse(rep(par=="a"|par=="both",nrow(data)),parts$Umsy_t,NA),
   Umsy_wgt=ifelse(rep(par=="a"|par=="both",nrow(data)),parts$Umsy_wt,NA),
   sigma=aa$summary["sigma","50%"],
   pi=aa$summary[grep("pi1",row.names(aa$summary)),"50%"],
   A=aa$summary[grep("A",row.names(aa$summary)),"50%"],
   probregime =matrix(aa$summary[grep("gamma\\[",row.names(aa$summary)),"50%"],ncol=k_regime, byrow=T),
   regime = aa$summary[grep("^zstar",row.names(aa$summary)),"50%"],
   summary=aa$summary,
   full_posterior=mc,
   stanfit=fit) 

  
  if(any(aa$rhat)>=1.05){print='Warning, some R_hat values are over the threshold of 1.05 - check parameter summary'}


}







#' Sampling control options. 
#'
#' Any arguments to pass to [rstan::sampling].
#'
#' @param eval.max Maximum number of evaluations of the objective function
#'   allowed.
#' @param iter.max Maximum number of iterations allowed.
#' @param ... Anything else. See the 'Control parameters' section of
#'   [rstan::stan].
#'
#' @export
stancontrol <- function(adapt_delta = 0.99,max_treedepth = 20,...) {
  list(adapt_delta = adapt_delta, max_treedepth = max_treedepth,...)
}




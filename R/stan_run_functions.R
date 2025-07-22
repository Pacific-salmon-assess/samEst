#' Simple Ricker model estimated with stan
#'
#' @param data A list or data frame containing complete vectors for: brood year (by), Spawners (s) and log(Recruits/Spawners) (logRS) time series.
#' @param AC Logical. Are residuals autocorrelated? Default is FALSE
#' @param smax_priors Options for custom Smax (capacity) priors - 2 values for the mean and the scale (variance). Defaults to mean = half maximum observed spawners, scale = max. observed spawners (from data) 
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
#' rickerstan(data=harck)
#' 
ricker_stan <- function(data,  ac=FALSE, smax_priors=NULL,smax_dists=c('normal','lognormal','cauchy'),mod=NULL,full_posterior=TRUE, control = stancontrol(adapt_delta=0.99), warmup=300, chains = 6, iter = 1000,...) {
  smax_dists=match.arg(smax_dists,choices=c('normal','lognormal','cauchy'))

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
    if(smax_dists=='normal'){
      datm$smax_dist=1
    }
    if(smax_dists=='lognormal'){
      datm$smax_dist=2
    }
    if(smax_dists=='cauchy'){
      datm$smax_dist=1
    }
    
  }else{
      datm = list(N=nrow(data),
                  L=max(data$by)-min(data$by)+1,
                  ii=data$by-min(data$by)+1,
                  R_S =data$logRS,
                  S=data$S,
                  pSmax_mean=smax_priors[1],
                  pSmax_sig=smax_priors[2])
      if(smax_dists=='normal'){
        datm$smax_dist=1
      }
      if(smax_dists=='lognormal'){
        datm$smax_dist=2
      }
      if(smax_dists=='cauchy'){
        datm$smax_dist=3
      }
    }

  fit<-rstan::sampling(sm, data=datm,
                      control = control, warmup = warmup, 
                      chains = chains, iter = iter,verbose=FALSE)
  #fit <- rstan::stan(model_code = sm, 
  #                      data = datm,
  #                      control = control, warmup = warmup, chains = chains, iter = iter)
  #
  if(ac==F){ mc <- rstan::extract(fit,pars=c('log_a','b','Smax','Smsy','Umsy','sigma','mu','epsilon'),permuted=T) 
  mc2=as.data.frame(do.call(cbind,mc))
  colnames(mc2)=c('log_a','b','Smax','Smsy','Umsy','sigma',paste('mu[',seq(1:datm$N),']',sep=''),paste('epsilon[',seq(1:datm$N),']',sep=''))
  } 
  if(ac==T){ mc <- rstan::extract(fit,pars=c('log_a','b','Smax','Smsy','Umsy','sigma','mu','epsilon','rho'),permuted=T)
  mc2=as.data.frame(do.call(cbind,mc))
  colnames(mc2)=c('log_a','b','Smax','Smsy','Umsy','sigma',paste('mu[',seq(1:datm$N),']',sep=''),paste('epsilon[',seq(1:datm$N),']',sep=''),'rho')
  } 
  
  aa <- rstan::summary(fit)
  if(full_posterior==FALSE){
     if(ac==F){
       ans<-list(data=data,
                 alpha=c(median=aa$summary["log_a","50%"],med.cv=aa$summary["log_a","sd"]/aa$summary["log_a","50%"],est2.5=aa$summary["log_a","2.5%"],est97.5=aa$summary["log_a","97.5%"]),
                 beta=c(median=aa$summary["b","50%"],med.cv=aa$summary["b","sd"]/aa$summary["b","50%"],est2.5=aa$summary["b","2.5%"],est97.5=aa$summary["b","97.5%"]),
                 Smax=c(median=aa$summary["Smax","50%"],med.cv=aa$summary["Smax","sd"]/aa$summary["Smax","50%"],est2.5=aa$summary["Smax","2.5%"],est97.5=aa$summary["Smax","97.5%"]),
                 Smsy=c(median=aa$summary["Smsy","50%"],med.cv=aa$summary["Smsy","sd"]/aa$summary["Smsy","50%"],est2.5=aa$summary["Smsy","2.5%"],est97.5=aa$summary["Smsy","97.5%"]),
                 Umsy=c(median=aa$summary["Umsy","50%"],med.cv=aa$summary["Umsy","sd"]/aa$summary["Umsy","50%"],est2.5=aa$summary["Umsy","2.5%"],est97.5=aa$summary["Umsy","97.5%"]),
                 sigma=c(median=aa$summary["sigma","50%"],med.cv=aa$summary["sigma","sd"]/aa$summary["sigma","50%"],est2.5=aa$summary["sigma","2.5%"],est97.5=aa$summary["sigma","97.5%"]))
     }
    if(ac==T){
      ans<-list(data=data,
                alpha=c(median=aa$summary["log_a","50%"],med.cv=aa$summary["log_a","sd"]/aa$summary["log_a","50%"],est2.5=aa$summary["log_a","2.5%"],est97.5=aa$summary["log_a","97.5%"]),
                beta=c(median=aa$summary["b","50%"],med.cv=aa$summary["b","sd"]/aa$summary["b","50%"],est2.5=aa$summary["b","2.5%"],est97.5=aa$summary["b","97.5%"]),
                Smax=c(median=aa$summary["Smax","50%"],med.cv=aa$summary["Smax","sd"]/aa$summary["Smax","50%"],est2.5=aa$summary["Smax","2.5%"],est97.5=aa$summary["Smax","97.5%"]),
                Smsy=c(median=aa$summary["Smsy","50%"],med.cv=aa$summary["Smsy","sd"]/aa$summary["Smsy","50%"],est2.5=aa$summary["Smsy","2.5%"],est97.5=aa$summary["Smsy","97.5%"]),
                Umsy=c(median=aa$summary["Umsy","50%"],med.cv=aa$summary["Umsy","sd"]/aa$summary["Umsy","50%"],est2.5=aa$summary["Umsy","2.5%"],est97.5=aa$summary["Umsy","97.5%"]),
                sigma=c(median=aa$summary["sigma","50%"],med.cv=aa$summary["sigma","sd"]/aa$summary["sigma","50%"],est2.5=aa$summary["sigma","2.5%"],est97.5=aa$summary["sigma","97.5%"]),
                rho=c(median=aa$summary["rho","50%"],med.cv=aa$summary["rho","sd"]/aa$summary["rho","50%"],est2.5=aa$summary["rho","2.5%"],est97.5=aa$summary["rho","97.5%"]))
    }
    
  }else{
    ans<- list(data=data,fit=fit,summary=aa$summary,samples=mc2)
  }
  return(ans)
}








#' random walks Ricker model estimated with stan
#'
#' @param data A data frame containing Spawners (S) and log(Recruits/Spawners) (R_S) time series. Use sr_format for the correct column names. 
#' @param par Which parameter should vary? Either productivity (intercept, a), capacity (slope, b) or both parameters
#' @param control output of stancontrol
#' @param warmup To be passed to rstan::sampling. A positive integer specifying the number of warmup (aka burnin) iterations per
#'  chain. The default is 200.
#' @param chains To be passed to rstan::sampling. A positive integer specifying the number of Markov chains. The default is 6.
#' @param iter To be passed to rstan::sampling. A positive integer specifying the number of iterations for each chain 
#' (including warmup). The default is 1000.
#' @param ... Anything else that would be passed to rstan::sampling
#' @param sm_ext Default is null, external stan rw model. Implemented to speed up simulation evaluation, use with caution.
#' @param lambertW Logical, indicating if lambertW functions should be used in stan code, requires git installation from stan
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
#' ricker_rw_stan(data=harck)
#' 
ricker_rw_stan <- function(data, par=c('a','b','both'),smax_priors=NULL,smax_dists=c('normal','lognormal','cauchy'),full_posterior=TRUE,control = stancontrol(), mod=NULL,
  warmup=300,  chains = 6, iter = 1000,...) {
  par=match.arg(par,choices=c('a','b','both'))
  smax_dists=match.arg(smax_dists,choices=c('normal','lognormal','cauchy'))
  
  if(is.null(mod)==T){
    sm=sr_mod(type='rw',par=par)
  }else{sm=mod}
  
  if(is.null(smax_priors)==TRUE){
    datm = list(N=nrow(data),
                L=max(data$by)-min(data$by)+1,
                ii=data$by-min(data$by)+1,
                R_S =data$logRS,
                S=data$S,
                pSmax_mean=max(data$S)/2,
                pSmax_sig=max(data$S)*2)
    if(smax_dists=='normal'){
      datm$smax_dist=1
    }
    if(smax_dists=='lognormal'){
      datm$smax_dist=2
    }
    if(smax_dists=='cauchy'){
      datm$smax_dist=3
    }
    
  }else{
    datm = list(N=nrow(data),
                L=max(data$by)-min(data$by)+1,
                ii=data$by-min(data$by)+1,
                R_S =data$logRS,
                S=data$S,
                pSmax_mean=smax_priors[1],
                pSmax_sig=smax_priors[2])
    if(smax_dists=='normal'){
      datm$smax_dist=1
    }
    if(smax_dists=='lognormal'){
      datm$smax_dist=2
    }
    if(smax_dists=='cauchy'){
      datm$smax_dist=3
    }
  }
  
  fit <- rstan::sampling(sm, data = datm,
                         control = control, warmup = warmup, chains = chains, iter = iter,verbose=FALSE)
  

  
  aa <- rstan::summary(fit)
  if(full_posterior==FALSE){
    if(par=='a'){
      ans<-list(data=data,
                alpha=data.frame(median=aa$summary[paste('log_a[',seq(1:datm$L),']',sep=''),"50%"],med.cv=aa$summary[paste('log_a[',seq(1:datm$L),']',sep=''),"sd"]/aa$summary[paste('log_a[',seq(1:datm$L),']',sep=''),"50%"],est2.5=aa$summary[paste('log_a[',seq(1:datm$L),']',sep=''),"2.5%"],est97.5=aa$summary[paste('log_a[',seq(1:datm$L),']',sep=''),"97.5%"]),
                beta=c(median=aa$summary["b","50%"],med.cv=aa$summary["b","sd"]/aa$summary["b","50%"],est2.5=aa$summary["b","2.5%"],est97.5=aa$summary["b","97.5%"]),
                Smax=c(median=aa$summary["Smax","50%"],med.cv=aa$summary["Smax","sd"]/aa$summary["Smax","50%"],est2.5=aa$summary["Smax","2.5%"],est97.5=aa$summary["Smax","97.5%"]),
                Smsy=data.frame(median=aa$summary[paste('Smsy[',seq(1:datm$L),']',sep=''),"50%"],med.cv=aa$summary[paste('Smsy[',seq(1:datm$L),']',sep=''),"sd"]/aa$summary[paste('Smsy[',seq(1:datm$L),']',sep=''),"50%"],est2.5=aa$summary[paste('Smsy[',seq(1:datm$L),']',sep=''),"2.5%"],est97.5=aa$summary[paste('Smsy[',seq(1:datm$L),']',sep=''),"97.5%"]),
                Umsy=data.frame(median=aa$summary[paste('Umsy[',seq(1:datm$L),']',sep=''),"50%"],med.cv=aa$summary[paste('Umsy[',seq(1:datm$L),']',sep=''),"sd"]/aa$summary[paste('Umsy[',seq(1:datm$L),']',sep=''),"50%"],est2.5=aa$summary[paste('Umsy[',seq(1:datm$L),']',sep=''),"2.5%"],est97.5=aa$summary[paste('Umsy[',seq(1:datm$L),']',sep=''),"97.5%"]),
                sigma=c(median=aa$summary["sigma","50%"],med.cv=aa$summary["sigma","sd"]/aa$summary["sigma","50%"],est2.5=aa$summary["sigma","2.5%"],est97.5=aa$summary["sigma","97.5%"]),
                sigma_a=c(median=aa$summary["sigma_a","50%"],med.cv=aa$summary["sigma_a","sd"]/aa$summary["sigma_a","50%"],est2.5=aa$summary["sigma_a","2.5%"],est97.5=aa$summary["sigma_a","97.5%"]))
    
    }
    if(par=='b'){
      ans<-list(data=data,
                alpha=c(median=aa$summary["log_a","50%"],med.cv=aa$summary["log_a","sd"]/aa$summary["log_a","50%"],est2.5=aa$summary["log_a","2.5%"],est97.5=aa$summary["log_a","97.5%"]),
                beta=data.frame(median=aa$summary[paste('beta[',seq(1:datm$L),']',sep=''),"50%"],med.cv=aa$summary[paste('beta[',seq(1:datm$L),']',sep=''),"sd"]/aa$summary[paste('beta[',seq(1:datm$L),']',sep=''),"50%"],est2.5=aa$summary[paste('beta[',seq(1:datm$L),']',sep=''),"2.5%"],est97.5=aa$summary[paste('beta[',seq(1:datm$L),']',sep=''),"97.5%"]),
                Smax=data.frame(median=aa$summary[paste('Smax[',seq(1:datm$L),']',sep=''),"50%"],med.cv=aa$summary[paste('Smax[',seq(1:datm$L),']',sep=''),"sd"]/aa$summary[paste('Smax[',seq(1:datm$L),']',sep=''),"50%"],est2.5=aa$summary[paste('Smax[',seq(1:datm$L),']',sep=''),"2.5%"],est97.5=aa$summary[paste('Smax[',seq(1:datm$L),']',sep=''),"97.5%"]),
                Smsy=data.frame(median=aa$summary[paste('Smsy[',seq(1:datm$L),']',sep=''),"50%"],med.cv=aa$summary[paste('Smsy[',seq(1:datm$L),']',sep=''),"sd"]/aa$summary[paste('Smsy[',seq(1:datm$L),']',sep=''),"50%"],est2.5=aa$summary[paste('Smsy[',seq(1:datm$L),']',sep=''),"2.5%"],est97.5=aa$summary[paste('Smsy[',seq(1:datm$L),']',sep=''),"97.5%"]),
                Umsy=c(median=aa$summary["Umsy","50%"],med.cv=aa$summary["Umsy","sd"]/aa$summary["Umsy","50%"],est2.5=aa$summary["Umsy","2.5%"],est97.5=aa$summary["Umsy","97.5%"]),
                sigma=c(median=aa$summary["sigma","50%"],med.cv=aa$summary["sigma","sd"]/aa$summary["sigma","50%"],est2.5=aa$summary["sigma","2.5%"],est97.5=aa$summary["sigma","97.5%"]),
                sigma_b=c(median=aa$summary["sigma_b","50%"],med.cv=aa$summary["sigma_b","sd"]/aa$summary["sigma_b","50%"],est2.5=aa$summary["sigma_b","2.5%"],est97.5=aa$summary["sigma_b","97.5%"]))
      
    }
    if(par=='both'){
      ans<-list(data=data,
                alpha=data.frame(median=aa$summary[paste('log_a[',seq(1:datm$L),']',sep=''),"50%"],med.cv=aa$summary[paste('log_a[',seq(1:datm$L),']',sep=''),"sd"]/aa$summary[paste('log_a[',seq(1:datm$L),']',sep=''),"50%"],est2.5=aa$summary[paste('log_a[',seq(1:datm$L),']',sep=''),"2.5%"],est97.5=aa$summary[paste('log_a[',seq(1:datm$L),']',sep=''),"97.5%"]),
                beta=data.frame(median=aa$summary[paste('beta[',seq(1:datm$L),']',sep=''),"50%"],med.cv=aa$summary[paste('beta[',seq(1:datm$L),']',sep=''),"sd"]/aa$summary[paste('beta[',seq(1:datm$L),']',sep=''),"50%"],est2.5=aa$summary[paste('beta[',seq(1:datm$L),']',sep=''),"2.5%"],est97.5=aa$summary[paste('beta[',seq(1:datm$L),']',sep=''),"97.5%"]),
                Smax=data.frame(median=aa$summary[paste('Smax[',seq(1:datm$L),']',sep=''),"50%"],med.cv=aa$summary[paste('Smax[',seq(1:datm$L),']',sep=''),"sd"]/aa$summary[paste('Smax[',seq(1:datm$L),']',sep=''),"50%"],est2.5=aa$summary[paste('Smax[',seq(1:datm$L),']',sep=''),"2.5%"],est97.5=aa$summary[paste('Smax[',seq(1:datm$L),']',sep=''),"97.5%"]),
                Smsy=data.frame(median=aa$summary[paste('Smsy[',seq(1:datm$L),']',sep=''),"50%"],med.cv=aa$summary[paste('Smsy[',seq(1:datm$L),']',sep=''),"sd"]/aa$summary[paste('Smsy[',seq(1:datm$L),']',sep=''),"50%"],est2.5=aa$summary[paste('Smsy[',seq(1:datm$L),']',sep=''),"2.5%"],est97.5=aa$summary[paste('Smsy[',seq(1:datm$L),']',sep=''),"97.5%"]),
                Umsy=data.frame(median=aa$summary[paste('Umsy[',seq(1:datm$L),']',sep=''),"50%"],med.cv=aa$summary[paste('Umsy[',seq(1:datm$L),']',sep=''),"sd"]/aa$summary[paste('Umsy[',seq(1:datm$L),']',sep=''),"50%"],est2.5=aa$summary[paste('Umsy[',seq(1:datm$L),']',sep=''),"2.5%"],est97.5=aa$summary[paste('Umsy[',seq(1:datm$L),']',sep=''),"97.5%"]),
                sigma=c(median=aa$summary["sigma","50%"],med.cv=aa$summary["sigma","sd"]/aa$summary["sigma","50%"],est2.5=aa$summary["sigma","2.5%"],est97.5=aa$summary["sigma","97.5%"]),
                sigma_a=c(median=aa$summary["sigma_a","50%"],med.cv=aa$summary["sigma_a","sd"]/aa$summary["sigma_a","50%"],est2.5=aa$summary["sigma_a","2.5%"],est97.5=aa$summary["sigma_a","97.5%"]),
                sigma_b=c(median=aa$summary["sigma_b","50%"],med.cv=aa$summary["sigma_b","sd"]/aa$summary["sigma_b","50%"],est2.5=aa$summary["sigma_b","2.5%"],est97.5=aa$summary["sigma_b","97.5%"]))
      
    }
    
  }else{
    if(par=='a'){
      mc <- rstan::extract(fit,pars=c('log_a','b','Smax','Smsy','Umsy','sigma','sigma_a'),permuted=T)
      mc2=as.data.frame(do.call(cbind,mc))
      colnames(mc2)=c(paste('log_a[',seq(1:datm$L),']',sep=''),'b','Smax',paste('Smsy[',seq(1:datm$L),']',sep=''),paste('Umsy[',seq(1:datm$L),']',sep=''),'sigma','sigma_a')
    }
    if(par=='b'){
      mc <- rstan::extract(fit,pars=c('log_a','b','Smax','Smsy','Umsy','sigma','sigma_b'),permuted=T)
      mc2=as.data.frame(do.call(cbind,mc))
      colnames(mc2)=c('log_a',paste('b[',seq(1:datm$L),']',sep=''),paste('Smax[',seq(1:datm$L),']',sep=''),paste('Smsy[',seq(1:datm$L),']',sep=''),'Umsy','sigma','sigma_b')
    }
    if(par=='both'){
      mc <- rstan::extract(fit,pars=c('log_a','b','Smax','Smsy','Umsy','sigma','sigma_a','sigma_b'),permuted=T)
      mc2=as.data.frame(do.call(cbind,mc))
      colnames(mc2)=c(paste('log_a[',seq(1:datm$L),']',sep=''),paste('b[',seq(1:datm$L),']',sep=''),paste('Smax[',seq(1:datm$L),']',sep=''),paste('Smsy[',seq(1:datm$L),']',sep=''),paste('Umsy[',seq(1:datm$L),']',sep=''),'sigma','sigma_a','sigma_b')
    }
   
    ans<- list(data=datm,
               fit=fit,
               summary=aa$summary,
               samples=mc2)
  }
  return(ans)

}



#' random walks Ricker model estimated with stan
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
#' @param lambertW Logical, indicating if lambertW functions should be used in stan code, requires git installation from stan
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
  #par='both'
  par=match.arg(par,choices=c('a','b','both'))
  smax_dists=match.arg(smax_dists,choices=c('normal','lognormal','cauchy'))
  
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
    if(smax_dists=='normal'){
      datm$smax_dist=1
    }
    if(smax_dists=='lognormal'){
      datm$smax_dist=2
    }
    if(smax_dists=='cauchy'){
      datm$smax_dist=3
    }
    
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
    if(smax_dists=='normal'){
      datm$smax_dist=1
    }
    if(smax_dists=='lognormal'){
      datm$smax_dist=2
    }
    if(smax_dists=='cauchy'){
      datm$smax_dist=3
    }
  }
  
  #if(is.null(sm_ext)){
  #  sm <- sr_mod(type='hmm',ac=FALSE,par=par,loglik=FALSE, modelcode=TRUE)
  #}else{
  #  sm <- sm_ext
  #}
 
 
  if(is.null(smax_priors)==TRUE){
  fit <- rstan::sampling(sm, 
                        data = list(N=nrow(data),
                                    R_S=data$logRS,
                                    S=data$S,
                                    K=k_regime,
                                    alpha_dirichlet=dirichlet_prior,
                                    pSmax_mean=max(data$S)/2,
                                    pSmax_sig=max(data$S)),
                        control = control, warmup = warmup, chains = 1, iter = iter,verbose=T)
  }
  if(is.null(smax_priors)==FALSE){
    fit <- rstan::sampling(sm, 
                           data = list(N=nrow(data),
                                       R_S =data$logRS,
                                       S=data$S,
                                       K=k_regime,
                                       alpha_dirichlet=dirichlet_prior,
                                       pSmax_mean=smax_priors[1],
                                       pSmax_sig=smax_priors[2]),
                           control = control, warmup = warmup, chains = chains, iter = iter,verbose=FALSE)
    
  }

  mc <- rstan::extract(fit, 
          inc_warmup=FALSE, permuted=FALSE)
    

  aa <- rstan::summary(fit)
  parts <-stan_regime_rps(m=fit,par=par)
  if(full_posterior==FALSE){
  ans<-list(
   alpha_regime=ifelse(rep(par=="a"|par=="both",nrow(data)),parts$log_a_t,NA),
   alpha_wgt=ifelse(rep(par=="a"|par=="both",nrow(data)),parts$log_a_wt,NA),
   Smax_regime=ifelse(rep(par=="b"|par=="both",nrow(data)),parts$Smax_t,NA),
   Smax_wgt=ifelse(rep(par=="b"|par=="both",nrow(data)),parts$Smax_wt,NA),
   alpha=ifelse(par=="a"|par=="both",list(aa$summary[grep("log_a\\[",row.names(aa$summary)),"50%"]),aa$summary["log_a","50%"])[[1]],
   beta=ifelse(par=="b"|par=="both",list(aa$summary[grep("b\\[",row.names(aa$summary)),"50%"]),aa$summary["b","50%"])[[1]],
   Smax=ifelse(par=="b"|par=="both",list(aa$summary[grep("Smax\\[",row.names(aa$summary)),"50%"]),aa$summary["Smax","50%"])[[1]],
   Smsy_regime=parts$Smsy_t,
   Smsy_wgt=parts$Smsy_wt,
   Umsy_regime=ifelse(rep(par=="a"|par=="both",nrow(data)),parts$Umsy_t,NA),
   Umsy_wgt=ifelse(rep(par=="a"|par=="both",nrow(data)),parts$Umsy_wt,NA),
   Smsy=aa$summary[grep("Smsy\\[",row.names(aa$summary)),"50%"],
   umsy=ifelse(par=="a"|par=="both",list(aa$summary[grep("Umsy\\[",row.names(aa$summary)),"50%"]),aa$summary["Umsy","50%"])[[1]],
   sigobs=aa$summary["sigma","50%"],
   pi=aa$summary[grep("pi1",row.names(aa$summary)),"50%"],
   A=aa$summary[grep("A",row.names(aa$summary)),"50%"],
   probregime =matrix(aa$summary[grep("gamma\\[",row.names(aa$summary)),"50%"],ncol=k_regime, byrow=T),
   regime = aa$summary[grep("^zstar",row.names(aa$summary)),"50%"]) 

  }else{
    ans<- list(fit=fit,summary=aa,samples=mc)
  }
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
  list(adapt_delta = 0.99,max_treedepth = 20,...)
}




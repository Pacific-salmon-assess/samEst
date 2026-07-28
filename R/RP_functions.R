#Function to calculate reference points once SR parameters are computed. 






#sgen functions from samSim



#' compute Sgen likelihood function
#'
#' @param S Spawner numbers set to the interval between 0 and Smsy in sGenCalc.
#' @param loga alpha parameter in Ricker function: R=S*exp(a-b*S)
#' @param b beta parameter in Ricker function: R=S*exp(a-b*S)
#' @param Smsy estimate of Smsy based on the alpha and beta parameters above.
#' 
#' @export
#' 
#' @returns Sgen likelihood 
#' 
#' 
#' 
Sgencompute <- function(S, loga, b, Smsy ) {
	#modified from samsim sGenOptimSmsyum
  
  prt <- S * exp(loga - b * S)
  epsilon <- log(Smsy) - log(prt)
  nLogLike <- sum(dnorm(epsilon, 0, 0.01, log = T))
  return(list(nSS = sum(nLogLike)))
  
}





#' compute Sgen estimate
#'
#' @param loga alpha parameter in Ricker function: R=S*exp(loga-b*S)
#' @param b beta parameter in Ricker function: R=S*exp(loga-b*S)
#' @param Smsy estimate of Smsy based on the alpha and beta parameters above.
#' 
#' @export
#' 
#' @returns Sgen estimate 
#' 
#' @examples
#' m1<- ricker_TMB(data=harck)
#' sgen<-sGenCalc(m1$alpha,m1$beta,m1$Smsy)
#' 
#' 
#' 
#' 
sGenCalc <- function(loga,b, Smsy) {
  #gives the min Ricker log-likelihood
  if(loga>0){
    fnSGen <- function(S, loga, b, Smsy) -1.0 * Sgencompute(S, loga, b, Smsy)$nSS
    fit <- optimize(f = fnSGen, interval = c(0, ((loga / b) * (0.5 - 0.07 * loga))),
                 loga=loga,b=b, Smsy = Smsy)
  }else{
    fit <-list(minimum=NA)
  }

  return(list(fit = fit$minimum))
}






#' compute Smsy estimate based on Scheuerell 2016. An explicit solution for 
#' calculating optimum spawning stock size from Ricker’s stock recruitment model
#'
#' @param loga alpha parameter in Ricker function: R=S*exp(a-b*S)
#' @param b beta parameter in Ricker function: R=S*exp(a-b*S)
#' @param Smsy estimate of Smsy based on the alpha and beta parameters above.
#' 
#' @importFrom gsl lambert_W0
#' 
#' @returns Smsy estimate 
#' 
#' @export
#' 
smsyCalc <- function(loga,b) {
  #gives the min Ricker log-likelihood
  Smsy <- (1 - gsl::lambert_W0(exp(1 - loga))) /b

  return(Smsy)
}



#' compute Umsy estimate based on Scheuerell 2016. An explicit solution for 
#' calculating optimum spawning stock size from Ricker’s stock recruitment model
#'
#' @param loga alpha parameter in Ricker function: R=S*exp(loga-b*S)
#' @param b beta parameter in Ricker function: R=S*exp(loga-b*S)
#' @param Smsy estimate of Smsy based on the alpha and beta parameters above.
#' 
#' @importFrom gsl lambert_W0
#' 
#' @returns Smsy estimate 
#' 
#' @export
#' 
umsyCalc <- function(loga) {
  #gives the min Ricker log-likelihood
  umsy <- (1 - gsl::lambert_W0(exp(1 - loga)))

  return(umsy)
}

#' Compute reference points from parameters of Hidden Markov Model (regime shift) Stan model
#'
#' @param m estimated fit for rstan model
#' @returns Estimates for Smax, Smsy, and Umsy over time, either unweighted (ie. returns parameters for the most likely regime sequence) or weighted (ie. probability of each regime x regime parameters) 
#' 
#' @export
#' 
stan_regime_rps<- function(m,par=c('a','b','both'),lambertW=FALSE){
  d=extract(m)
  if(par=='a'){
    log_a=apply(d$logalpha,2,median)
    beta=median(d$beta)
    Smax=median(d$Smax)
    gamma=cbind(apply(d$gamma[,,1],2,median),apply(d$gamma[,,2],2,median))
    zstar=apply(d$zstar,2,median)
    
    #Most likely sequence of parameters:
    log_a_t=log_a[zstar]
    
    #Weighted parameters (probability of regime state X regime parameter)
    log_a_wt=gamma%*%log_a
    
    if(lambertW){
      Umsy=apply(d$Umsy,2,median)
      Smsy=apply(d$Smsy,2,median)   
      Umsy_t=Umsy[zstar]
      Smsy_t=Smsy[zstar]
      Umsy_wt = gamma%*%Umsy
      Smsy_wt = gamma%*%Smsy
      ans<-data.frame(log_a_t,log_a_wt,Umsy_t,Umsy_wt,Smsy_t,Smsy_wt)
    }else{
      ans<-data.frame(log_a_t,log_a_wt)
    }
  }
  if(par=='b'){
    log_a=median(d$logalpha)
    beta=apply(d$beta,2,median)
    Smax=apply(d$Smax,2,median)
    
    gamma=cbind(apply(d$gamma[,,1],2,median),apply(d$gamma[,,2],2,median))
    zstar=apply(d$zstar,2,median)
    
    #Most likely sequence of parameters:
    beta_t=beta[zstar]
    Smax_t=Smax[zstar]
   
    #Weighted parameters (probability of regime state X regime parameter)
    beta_wt=gamma%*%beta
    Smax_wt = gamma%*%Smax
    
    
    if(lambertW){
      Umsy=median(d$Umsy)
      Smsy=apply(d$Smsy,2,median)
      Smsy_t=S_msy[zstar]
      Smsy_wt = gamma%*%Smsy
      ans<-data.frame(beta_t,beta_wt,Smax_t,Smax_wt,Smsy_t,Smsy_wt)
    }else{
      ans<-data.frame(beta_t,beta_wt,Smax_t,Smax_wt)
    }

  }
  if(par=='both'){
    log_a=apply(d$logalpha,2,median)
    beta=apply(d$beta,2,median)
    Smax=apply(d$Smax,2,median)
    
    gamma=cbind(apply(d$gamma[,,1],2,median),apply(d$gamma[,,2],2,median))
    zstar=apply(d$zstar,2,median)
    
    #Most likely sequence of parameters:
    log_a_t=log_a[zstar]
    beta_t=beta[zstar]
    Smax_t=Smax[zstar]
    
    #Weighted parameters (probability of regime state X regime parameter)
    log_a_wt=gamma%*%log_a
    beta_wt=gamma%*%beta
    Smax_wt = gamma%*%Smax
   

    if(lambertW){
      Umsy=apply(d$Umsy,2,median)
      Smsy=apply(d$Smsy,2,median)
      Umsy_t=Umsy[zstar]
      Smsy_t=Smsy[zstar]
      Umsy_wt = gamma%*%Umsy
      Smsy_wt = gamma%*%Smsy
      ans<-data.frame(log_a_t,log_a_wt,beta_t,beta_wt,Smax_t,Smax_wt,Umsy_t,Umsy_wt,Smsy_t,Smsy_wt)
    }else{
      ans<-data.frame(log_a_t,log_a_wt,beta_t,beta_wt,Smax_t,Smax_wt)
    }
   
  }
  return(ans)  
  
}






#' compute Sgen estimate based on the lambertW function similarly to
#' the approach in Scheuerell 2016. An explicit solution for 
#' calculating optimum spawning stock size from Ricker’s stock recruitment model
#'
#' @param loga alpha parameter in Ricker function: R=S*exp(loga-b*S)
#' @param b beta parameter in Ricker function: R=S*exp(loga-b*S)
#' 
#' @importFrom gsl lambert_W0
#' 
#' @returns DirectSgen estimate 
#' 
#' @export
#' 
sgenCalcDirect <- function(loga, b){
 sMSY <- ( 1 - gsl::lambert_W0 (exp ( 1 - loga) ) ) / b
 a <- exp(loga)
 
 return(-1/b*gsl::lambert_W0(-b*sMSY/a))
}

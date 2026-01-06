# Functions to run the RTMB stock recruitment  models
#===============================================







#' Simple Ricker or AR1 model estimated with RTMB
#'
#' @param data A list or data frame containing Spawners (S), log(Recruits/Spawners) (logRS), and brood year data (by). 
#' @param silent Logical Silent or optimization details? default is FALSE
#' @param control output from TMBcontrol() function, to be passed to nlminb()
#' @param tmb_map optional, mapping list indicating if parameters should be estimated of fixed. 
#' Default is all parameters are estimated
#' @param AC Logical. Are residuals autocorrelated? Default is FALSE
#' @param priors_flag Integer, 1 priors are included in estimation model, 0 priors are not included.
#'  See details for priors documentation. 
#' @param stan_flag Integer, flag indicating wether or not TMB code will be used with TMBstan - Jacobian
#' adjustment implemented. Only works if priors flag is also 1, Default is 0, jacobian adjustment not included. 
#' @param sig_p_sd sd for half normal prior on sigma parameter. default is 1.
#' @param logb_p_mean mean for prior on log b, default is -12.
#' @param logb_p_sd sd for prior on log b, default is 3.#' 
#' @param dautoreg select model code for AR1
#'
#' @details  this function will
#
#' 
#' @returns penalized (if priors_flag= TRUE) negative log likelihood
#' 
#' @importFrom stats nlminb 
#' 
#' 
#' 
#' @examples
#' data(harck)
#' rickerTMB(data=harck)
#' 
ricker_RTMB2 <- function(data, silent = FALSE, control = TMBcontrol(), 
  tmb_map = list(), AC=FALSE, priors_flag=1, stan_flag=0, sig_p_sd=1.0,
  logb_p_mean=-12.0, logb_p_sd=3.0,dautoreg=0){

  
  dat<- list(
    obs_S = as.numeric(data$S),
    obs_logRS = data$logRS,
    by = data$by,
    priors_flag=priors_flag,
    stan_flag=stan_flag,
    sig_p_sd=sig_p_sd,
    logb_p_mean=logb_p_mean,
    logb_p_sd=logb_p_sd
  )
  
  magS <- log10_ceiling(max(dat$obs_S))
  initlm <- lm(obs_logRS~obs_S, data=dat)

  


  if(!AC){
    param <- list(
      logalpha  = initlm$coefficients[[1]],
      logbeta = ifelse(initlm$coefficients[[2]]>0,log(magS),log(-initlm$coefficients[[2]])),
      logsigma = log(1),
      rho=0.0
    )

    tmb_map$rho<-NA
    tmbfn<-function(param){ricker_RTMB_fn(param,dat=dat)}
  
  }else{
    param <- list(
      logalpha  = initlm$coefficients[[1]],
      logbeta = ifelse(initlm$coefficients[[2]]>0,log(magS),log(-initlm$coefficients[[2]])),
      logsigobs = log(1.0),
      ar1_phi=0.0
    )
    if(!dautoreg){
        tmbfn<-function(param){ricker_ac_RTMB_fn(param,dat=dat)}
      }else{
        tmbfn<-function(param){ricker_ac_RTMB_fn_dautoreg(param,dat=dat)}  
      }
   

    
  }
  obj <- RTMB::MakeADFun(tmbfn, param, 
      map = tmb_map,
      silent=TRUE)

  opt <- stats::nlminb(obj$par, 
                       obj$fn, 
                       obj$gr,
                       control = control)

  sd_report <- RTMB::sdreport(obj)
  conv <- get_convergence_diagnostics(sd_report)

  nll <- obj$fn()[1]
  npar <- length(param)
 
  AICc  <- 2*nll + 2*npar +(2*npar*(npar+1)/(nrow(data)-npar-1))
  BIC  <- 2*nll + npar*log(nrow(data))


  
  structure(list(
    logalpha   = obj$report()$logalpha,
    beta       = obj$report()$beta,
    Smax       = obj$report()$Smax,
    sig        = obj$report()$sigobs,
    sigar      = ifelse(AC,obj$report()$sigAR,NA),
    rho        = ifelse(AC,obj$report()$rho,NA),
    Smsy        = obj$report()$Smsy,
    umsy        = obj$report()$umsy,
    Sgen        = obj$report()$Sgen,
    AICc       = AICc,
    BIC        = BIC,
    residuals  = obj$report()$residuals,
    model      = opt,
    data   = data,
    params = param,
    map    = tmb_map,
    obj    = obj,
    gradients  = conv$final_grads,
    bad_eig    = conv$bad_eig,
    conv_problem= conv$conv_problem,
    call       = match.call(expand.dots = TRUE),
    sd_report  = sd_report),
    class      = ifelse(AC,"Ricker_autocorr","Ricker_simple"))


}







#' Simple Ricker model estimated with RTMB
#'
#' @param data A list or data frame containing Spawners (S), log(Recruits/Spawners) (logRS), and brood year data (by). 
#' @param silent Logical Silent or optimization details? default is FALSE
#' @param control output from TMBcontrol() function, to be passed to nlminb()
#' @param tmb_map optional, mapping list indicating if parameters should be estimated of fixed. 
#' Default is all parameters are estimated
#' @param AC Logical. Are residuals autocorrelated? Default is FALSE
#' @param priors_flag Integer, 1 priors are included in estimation model, 0 priors are not included.
#'  See details for priors documentation. 
#' @param stan_flag Integer, flag indicating wether or not TMB code will be used with TMBstan - Jacobian
#' adjustment implemented. Only works if priors flag is also 1, Default is 0, jacobian adjustment not included. 
#' @param sig_p_sd sd for half normal prior on sigma parameter. default is 1.
#' @param logb_p_mean mean for prior on log b, default is -12.
#' @param logb_p_sd sd for prior on log b, default is 3.#' 
#' @param dautoreg select model code for AR1
#'
#' @details  this function will
#
#' 
#' @returns penalized (if priors_flag= TRUE) negative log likelihood
#' 
#' @importFrom stats nlminb 
#' 
#'
#' 
#' @examples
#' data(harck)
#' rickerTMB(data=harck)
#' 
ricker_RTMB <- function(data, silent = FALSE, control = TMBcontrol(), 
  tmb_map = list(), AC=FALSE, priors_flag=1, stan_flag=0, sig_p_sd=1.0,
  logb_p_mean=-12.0, logb_p_sd=3.0,dautoreg=0){

  
  dat<- list(
    obs_S = as.numeric(data$S),
    obs_logRS = data$logRS,
    by = data$by,
    priors_flag=priors_flag,
    stan_flag=stan_flag,
    sig_p_sd=sig_p_sd,
    logb_p_mean=logb_p_mean,
    logb_p_sd=logb_p_sd
  )
  
  magS <- log10_ceiling(max(dat$obs_S))
  initlm <- lm(obs_logRS~obs_S, data=dat)

  


  if(!AC){
    param <- list(
      logalpha  = initlm$coefficients[[1]],
      logbeta = ifelse(initlm$coefficients[[2]]>0,log(magS),log(-initlm$coefficients[[2]])),
      logsigobs = log(1)
    )
    tmbfn<-function(param){ricker_RTMB_fn(param,dat=dat)}
  
  }else{
    param <- list(
      logalpha  = initlm$coefficients[[1]],
      logbeta = ifelse(initlm$coefficients[[2]]>0,log(magS),log(-initlm$coefficients[[2]])),
      logsigobs = log(1.0),
      ar1_phi=0.0
    )
    if(!dautoreg){
        tmbfn<-function(param){ricker_ac_RTMB_fn(param,dat=dat)}
      }else{
        tmbfn<-function(param){ricker_ac_RTMB_fn_dautoreg(param,dat=dat)}  
      }
   

    
  }
  obj <- RTMB::MakeADFun(tmbfn, param, 
      map = tmb_map,
      silent=TRUE)

  opt <- stats::nlminb(obj$par, 
                       obj$fn, 
                       obj$gr,
                       control = control)

  sd_report <- RTMB::sdreport(obj)
  conv <- get_convergence_diagnostics(sd_report)

  nll <- obj$fn()[1]
  npar <- length(param)
 
  AICc  <- 2*nll + 2*npar +(2*npar*(npar+1)/(nrow(data)-npar-1))
  BIC  <- 2*nll + npar*log(nrow(data))


  
  structure(list(
    logalpha   = obj$report()$logalpha,
    beta       = obj$report()$beta,
    Smax       = obj$report()$Smax,
    sig        = obj$report()$sigobs,
    sigar      = ifelse(AC,obj$report()$sigAR,NA),
    rho        = ifelse(AC,obj$report()$rho,NA),
    Smsy        = obj$report()$Smsy,
    umsy        = obj$report()$umsy,
    Sgen        = obj$report()$Sgen,
    AICc       = AICc,
    BIC        = BIC,
    residuals  = obj$report()$residuals,
    model      = opt,
    data   = data,
    params = param,
    map    = tmb_map,
    obj    = obj,
    gradients  = conv$final_grads,
    bad_eig    = conv$bad_eig,
    conv_problem= conv$conv_problem,
    call       = match.call(expand.dots = TRUE),
    sd_report  = sd_report),
    class      = ifelse(AC,"Ricker_autocorr","Ricker_simple"))


}


  









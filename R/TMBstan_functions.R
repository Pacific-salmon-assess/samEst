# Functions to run the TMB stock recruitment  models
#===============================================
#Structure of the code copied from the sdmTMB package :
#https://github.com/pbs-assess/sdmTMB


#' Simple Ricker model estimated with TMBstan
#'
#' @param data A list or data frame containing Spawners (S) and log(Recruits/Spawners) (logRS) time series. 
#' 
#' @param silent Logical Silent or optimization details? default is FALSE
#' @param control output from TMBcontrol() function, to be passed to nlminb()
#' @param tmb_map optional, mapping list indicating if parameters should be estimated of fixed. 
#' Default is all parameters are estimated
#' @param AC Logical. Are residuals autocorrelated? Default is FALSE
#' @param priors_flag Integer, 1 priors are included in estimation model, 0 priors are not included.
#'  See details for priors documentation. See details for priors documentation.
#' @param stan_flag Integer, flag indicating wether or not TMB code will be used with TMBstan - Jacobian
#' adjustment implemented. Default is 0, jacobian adjustment not included.
#' @param sig_p_sd sd for half normal prior on sigma parameter. default is 1.
#' @param chains number of MCMC chains
#' @param iter number of MCMC iterations per chain
#' 
#'
#'
#' @details Priors: Weakly informative priors are included for the main parameterst of the model:
#' alpha ~ gamma(3,1)
#' logbeta ~ N(-12,3)
#' sigobs ~ gamma(2,1/3) 
#' 
#' 
#' @returns a list containing several model outputs - posterior distributions:
#' * alpha - posterior distribution for the alpha parameter vector
#' * beta - posterior distribution for the beta parameter 
#' * sigobs - posterior distribution for the observation error sigma   
#' * Smax - posterior distribution for the Smax parameter vector
#' * umsy - posterior distribution for the umsy parameter 
#' * Smsy - posterior distribution for the Smsy           
#' * AICc - AICc values, given by 2*nll + 2*npar +(2*npar*(npar+1)/(nrow(data)-npar-1)), excluding prior components
#' * BIC - BIC values, excluding prior components
#' 
#' @importFrom stats nlminb 
#' 
#' @export
#' 
#' @examples
#' data(harck)
#' rickerTMB(data=harck)
#' 
ricker_TMBstan <- function(data,  silent = FALSE, control = stancontrol(), 
  tmb_map = list(), AC=FALSE, priors_flag=1, stan_flag=1,sig_p_sd=1,
  Smax_mean=230000,Smax_sd=230000,
  chains=6,iter=5000, warmup = floor(iter/2)) {
 
  priorslogSmax<-log_prior_params(Smax_mean,Smax_sd)
  logsmax_p_mean=priorslogSmax$logsmax_pr_mean
  logsmax_p_sd=priorslogSmax$logsmax_pr_sig

  tmb_data <- list(
    S = data$S,
    logRS = data$logRS,
    priors_flag=priors_flag,
    stan_flag=stan_flag,
    sig_p_sd=sig_p_sd,
    logsmax_p_sd=logsmax_p_sd,
    logsmax_p_mean=logsmax_p_mean,
    y_oos=mean(data$logRS),
    x_oos=mean(data$S)
  )
  
  magS <- log10_ceiling(max(data$S))
  initlm <- lm(logRS~S, data=data)
  tmb_random <- NULL

  tmb_params <- list(
      logalpha   = initlm$coefficients[[1]],
      logSmax = ifelse(initlm$coefficients[[2]]>0,log(magS),log(-1/initlm$coefficients[[2]])),
      logsigma = log(1),
      rho=0
    )
   lowlimit <- c(0.01,4,log(0.01),-100)
   hightlimit <- c(20,20,log(3),100)

  if(!AC){
    tmb_map$rho<-as.factor(NA)
    lowlimit <- c(0.01,4,log(0.01))
    hightlimit <- c(20,20,log(3))
  }  
  
  tmb_obj <- TMB::MakeADFun(data = tmb_data, 
                             parameters = tmb_params, 
                             map = tmb_map,
                             random = tmb_random, 
                             DLL = "Ricker_simple_autocorr", 
                             silent = silent)

  tmb_mcmc <- tmbstan::tmbstan(tmb_obj, chains=chains,
              iter=iter, init="random",
              lower=lowlimit , upper=hightlimit,
              control = control,
              warmup = warmup)

    mc <- extract(tmb_mcmc, pars=names(tmb_obj$par),
              inc_warmup=FALSE, permuted=FALSE)
    

    fit_summary <- summary(tmb_mcmc)

    posterior <- as.matrix(tmb_mcmc)
  
    alpha <- NULL
    sigobs <- NULL  
    beta <- NULL
    Smax <- NULL
    umsy <- NULL 
    Smsy <- NULL  
    rho <- NULL
    ll <- matrix(NA,nrow=nrow(posterior),ncol=length(data$S))
    pred_logR<-matrix(NA,nrow=nrow(posterior),ncol=length(data$S))
    pred_logRS<-matrix(NA,nrow=nrow(posterior),ncol=length(data$S))



    for(i in 1:nrow(posterior)){
      r <- tmb_obj$report(posterior[i,-ncol(posterior)])

      alpha[i] <- r$alpha
      sigobs[i] <- r$sigobs  
      beta[i] <- r$beta
      Smax[i] <- r$Smax      
      umsy[i] <- r$umsy  
      Smsy[i] <- r$Smsy  
      rho[i] <- r$rhoo  
      pred_logR[i,] <- r$pred_logR
      pred_logRS[i,] <- r$pred_logRS
      ll[i,] <-r$ll   
      
    }

    #calculate WAIC
    npar <- length(tmb_params)-length(tmb_map)
 
    elpd_1 <- apply(ll,2,log_mean_exp) #
    AICc  <- -2*sum(elpd_1) + 2*npar +( 2*npar*(npar+1)/(nrow(data)-npar-1))
    BIC  <- -2*sum(elpd_1) + npar*log(nrow(data))
  

  
  
  structure(list(
    alpha = alpha,
    beta = beta,
    Smax = Smax,
    sig = sigobs,  
    umsy = umsy,  
    Smsy = r$Smsy,  
    pred_logR = pred_logR,
    pred_logRS = pred_logRS,
    AICc = AICc,
    BIC=BIC,
    fit_summary =fit_summary,
    posterior=posterior))

}





#' Ricker model with random walk in a, b or both parameters with TMB
#'
#' @param data A list or data frame containing Spawners (S) and log(Recruits/Spawners) (logRS) time series. 
#' @param tv.par Which parameters should vary? Either productivity (intercept, a), capacity (slope, b) or both parameters
#' @param silent Logical Silent or optimization details? default is FALSE
#' @param control list of controls, as it would be passed to stan()
#' @param ini_param Optional. A list with initial parameter guesses. The list should contain: alphao (a number),
#' logbeta (a number), logsigobs (a number), logsiga (a number), and alpha ( a vector with the same length as the data). 
#' @param tmb_map optional, mapping list indicating if parameters should be estimated of fixed. 
#' Default is all parameters are estimated
#' @param priors_flag Integer, 1 priors are included in estimation model, 0 priors are not included.
#'  See details for priors documentation. See details for priors documentation.
#' @param stan_flag Integer, flag indicating wether or not TMB code will be used with TMBstan - Jacobian
#' adjustment implemented. Default is 0, jacobian adjustment not included.
#' @param sig_p_sd sd for half normal prior on sigma parameter. default is 1.
#' @param siga_p_sd sd for half normal prior on sigma for alpha random walk parameter. default is 1.
#' @param sigb_p_sd sd for half normal prior on sigma for beta random walk parameter. default is 1.
#' @param chains number of MCMC chains
#' @param iter number of MCMC iterations per chain
#' @param laplace Use the laplace approximation for random effects, default is FALSE. (TRUE not implemented)
#' 
#' @details Priors: Weakly informative priors are included for the main parameterst of the model:
#' alpha ~ gamma(3,1)
#' logbeta ~ N(-12,3)
#' sigobs ~ gamma(2,1/3)
#' siga ~ gamma(2,1/3) 
#' sigb ~ gamma(2,1/3)  
#' 
#' 
#' @returns a list containing several model outputs:
#' * alpha - MLE estimates for the alpha parameter vector
#' * beta - MLE estimates for the beta parameter 
#' * sig - MLE estimates for the observation error standard deviation     
#' * siga - MLE estimates for the process error (variation in alpha) standard deviation
#' * sigb - MLE estimates for the process error (variation in beta) standard deviation   
#' * AICc - AICc values, given by 2*nll + 2*npar +(2*npar*(npar+1)/(nrow(data)-npar-1)), excluding prior components
#' * BIC - BIC values, excluding prior components
#' * model - opt object, generated by the `stats::nlminb()`  
#' * tmb_data - data provided to TMB,
#' * tmb_params - parameter intial guesses,
#' * tmb_map - mapping indicating which parameters should be fixed or estimated,
#' * tmb_obj - obj object, generated by the `TMB::MakeADFun()` 
#' * gradients - final convergence gradient
#' * bad_eig - eigen value
#' * call - original function call
#' * sd_report - MLE estimates and sdt Error estimates for main parameter estimates
#' * class - name of cpp model
#' 
#' 
#' 
#' @examples 
#' data(harck)
#' ricker_ra_TMB(data=harck,tv.par='a')
#' 
#' 
ricker_rw_TMBstan <- function(data, tv.par=c('a','b','both'), silent = FALSE, 
  control = stancontrol(), ini_param=NULL, tmb_map = list(), priors_flag=1, stan_flag=1,
  sig_p_sd=1, siga_p_sd=1, sigb_p_sd=1, Smax_mean=230000,Smax_sd=230000,
   chains=6,iter=10000 ,laplace=FALSE, warmup = floor(iter/2),...) {

  #===================================
  #prepare TMB input and options
  #===================================
  tmb_data <- list(
    obs_S = data$S,
    obs_logRS = data$logRS,
    priors_flag=priors_flag,
    stan_flag=stan_flag,
    sig_p_sd=sig_p_sd,
    logsmax_p_mean=logsmax_p_mean,
    logsmax_p_sd=logsmax_p_sd,
    siga_p_sd=siga_p_sd,
    sigb_p_sd=sigb_p_sd
  )
  
   ##Smax lognormal prior
  priorslogSmax<-log_prior_params(Smax_mean,Smax_sd)
  logsmax_p_mean=priorslogSmax$logsmax_pr_mean
  logsmax_p_sd=priorslogSmax$logsmax_pr_sig



  if(is.null(ini_param)){
    magS <- log10_ceiling(max(data$S))
    initlm<-lm(logRS~S, data=data)
    
    tmb_params <- list(logalpha   = initlm$coefficients[[1]],
                   logSmax = ifelse(initlm$coefficients[[2]]>0,
                                   log(magS),
                                   log(-1/initlm$coefficients[[2]])),
                   logsigobs = log(.5),
                   logsiga = log(.5),
                   logsigb = log(.5),
                   epslogalpha_t=rep(0.1,length(tmb_data$obs_S)),
                   epslogsmax_t=rep(0.1,length(tmb_data$obs_S)))
  }else{
    tmb_params <-ini_param
  }


  if(tv.par=="a"){

    tmb_data$options_z=c(1,0)

    tmb_random <- "epslogalpha_t"

    tmb_map$logsigb = factor(NA)
    tmb_map$epslogsmax_t = factor( rep(NA,length(tmb_params$epslogsmax_t)) )
     
    npar <- 4
    npar_all <- 4+(length(data$S)-1)

    lowlimit <- c(0.01,4,log(0.01),log(0.01))
    hightlimit <- c(20,20,log(3),log(3))

  }else if(tv.par=="b"){

    tmb_data$options_z=c(0,1)

    tmb_random <- "epslogsmax_t"

    tmb_map$logsiga = factor(NA)
    tmb_map$epslogalpha_t = factor( rep(NA,length(tmb_params$epslogalpha_t)) )

    npar <- 4
    npar_all <- 4+(length(data$S)-1)

    lowlimit <- c(0.01,4,log(0.01),log(0.01))
    hightlimit <- c(20,20,log(3),log(3))

  }else if(tv.par=="both"){

    tmb_data$options_z=c(1,1)
       
    tmb_random <- c("epslogalpha_t", "epslogsmax_t")

    npar <- 5
    npar_all <- 5+(length(data$S)-1)*2

    lowlimit <- c(0.01,4,log(0.01),log(0.01),log(0.01))
    hightlimit <- c(20,20,log(3),log(3),log(3))

    

  }else{
    stop(paste("tv.par",tv.par,"not recognized."))
  }
   
  
  tmb_obj <- TMB::MakeADFun(data = tmb_data, 
                              parameters = tmb_params, 
                              map = tmb_map,
                              random = tmb_random, 
                              DLL = "Ricker_tv_all", 
                              silent = silent)

 

  #===================================
  # TMBstan fit
  #===================================
  #chains=6
  #iter=5000
  tmb_mcmc <- tmbstan::tmbstan(tmb_obj, chains=chains,
              iter=iter, init="random",
              lower=lowlimit , upper=hightlimit,
              control = control,
              laplace=laplace,
              warmup = warmup )

   
  fit_summary <- summary(tmb_mcmc)

  posterior <- as.matrix(tmb_mcmc)
  
  
  postlist<-as.list(as.data.frame(t(posterior[,-ncol(posterior)])))
   

  r <- lapply(postlist,tmb_obj$report)
     
      
  beta <-sapply(r,"[[","beta")
  alpha  <- sapply(r,"[[","alpha")
  sigobs <- sapply(r,"[[","sigobs")
  Smax <- sapply(r,"[[","Smax")
  umsy <-sapply(r,"[[","umsy")
  Smsy <- sapply(r,"[[","Smsy")
  ll <- sapply(r,"[[","ll")
  pred_logR <- sapply(r,"[[","pred_logR")
  pred_logRS <- sapply(r,"[[","pred_logRS")
    
      
  elpd_1 <- apply(ll,1,log_mean_exp) #
  AICc  <- -2*sum(elpd_1) + 2*npar +( 2*npar*(npar+1)/(nrow(data)-npar-1))
  BIC  <- -2*sum(elpd_1) + npar*log(nrow(data))
  
       
      
      
  
  #todo add alpha, beta and sigma parameter esitimates
  structure(list(
    alpha = alpha,
    beta = beta,
    Smax = Smax,
    sig = sigobs,  
    umsy = umsy,  
    Smsy = r$Smsy,  
    pred_logR = pred_logR,
    pred_logRS = pred_logRS,
    AICc = AICc,
    BIC=BIC,
    fit_summary =fit_summary,
    posterior=posterior))

}








#' Ricker hidden markov model with regime shiftts for logalpha and beta. 
#'
#' @param data A list or data frame containing Spawners (S) and log(Recruits/Spawners) (logRS) time series. 
#' @param tv.par Which parameters should vary? Either productivity (intercept, a), capacity (slope, b) or both parameters
#' @param k_regime Number of regimes to be considered
#' @param logalpha_limits vector containing two values: upper and lower limit for logalpha parameters. default is c(0,20)
#' @param beta_limits vector containing two values: upper and lower limit for beta parameters. default is c(1e-10,.1)
#' @param initDist Initial probability of being in each state, default is equal probabilities assigned to all regimes.
#' @param silent Logical Silent or optimization details? default is FALSE
#' @param control output from TMBcontrol() function, to be passed to nlminb()
#' @param ini_param Optional. A list with initial parameter guesses. The list should contain: lalpha (vector of length k_regime),
#' lbeta (vector of length k_regime), logsigma, pi1_tran (vector of length k_regime-1), 
#' and qij_tran (matrix with nrow=k_regime, and ncol=k_regime-1). Keep in mind that the main
#' parameters (logalpha and beta) are transformed with a logistic function to account for the custom upper and lower values.   
#' @param priors_flag Integer, 1 priors are included in estimation model, 0 priors are not included.
#'  See details for priors documentation. See details for priors documentation.
#' @param stan_flag Integer, flag indicating wether or not TMB code will be used with TMBstan - Jacobian
#' adjustment implemented. Default is 0, jacobian adjustment not included.
#' @param sig_p_sd sd for half normal prior on sigma parameter. default is 1.
#' @param dirichlet_prior k_regime x k_regime matrix. Prior for transition probability matrix, 
#' if NULL prior is set to matrix(1,nrow=k_regime,ncol=k_regime) 
#' 
#' @details This model was published in Tang et al. 2021 Identification of recruitment regime shifts with a hidden Markov stock-recruitment model. 
#' The code for this model was a contribution by Xiaozhuo Tang. 
#' 
#' Priors: Weakly informative priors are included for the main parameterst of the model:
#' logalpha ~ N(1.5,2.5)
#' logbeta ~ N(-12,3)
#' sigobs ~ gamma(2,1/3)
#' qi ~ gamma(2,1/3) 
#' sigb ~ gamma(2,1/3)  
#' 
#' @returns a list containing several model outputs:
#' * logalpha - MLE estimates for the logalpha parameter vector
#' * beta - MLE estimates for the beta parameter 
#' * sig  - MLE estimates for the observation error sigma     
#' * pi  - MLE estimates for the initial state probabilities vector of length k_regime
#' * A  - MLE estimates for the transition probabilities, matrix (k_regime x k_regime)   
#' * AICc - AICc values, given by 2*nll + 2*npar +(2*npar*(npar+1)/(nrow(data)-npar-1)), excluding prior components
#' * BIC - BIC values, excluding prior components
#' * model - opt object, generated by the `TMB::MakeADFun()`   
#' * tmb_data - data provided to TMB,
#' * tmb_params = parameter intial guesses,
#' * tmb_map    = mapping indicating which parameters should be fixed or estimated,
#' * tmb_random = tmb_random,
#' * tmb_obj    = tmb_obj,
#' * gradients  = conv$final_grads,
#' * bad_eig    = conv$bad_eig,
#' * call       = match.call(expand.dots = TRUE),
#' * sd_report  = sd_report),
#' * class      = "Ricker_tva"
#' 
#'  @references{
#'   \insertRef{tangIdentificationRecruitmentRegime2021}{samEst}
#' }
#' 
#' 
#' @examples 
#' data(harck)
#' ricker_HMM_TMB(data=harck)
ricker_hmm_TMBstan <- function(data, tv.par=c('a','b','both'),  
                              logalpha_limits = c(0.01,20), 
                              Smax_limits=c(10,100000000), 
                              initDist=NULL,
                              silent = FALSE, 
                              ini_param=NULL, 
                              tmb_map = list(), 
                              priors_flag=1, stan_flag=1,
                              sig_p_sd=1, 
                              Smax_mean=230000,Smax_sd=230000,
                              dirichlet_prior=NULL, 
                              control = stancontrol(),
                              chains=6, iter=10000 , laplace=FALSE, 
                              warmup = floor(iter/2),...) {


  ##Smax lognormal prior
  priorslogSmax<-log_prior_params(Smax_mean,Smax_sd)
  logsmax_p_mean=priorslogSmax$logsmax_pr_mean
  logsmax_p_sd=priorslogSmax$logsmax_pr_sig

  #===================================
  #prepare TMB input and options
  #===================================

  if(is.null(dirichlet_prior)){
    dirichlet_prior<-matrix(1,nrow=k_regime,ncol=k_regime)
  }else if(nrow(dirichlet_prior)!=k_regime |ncol(dirichlet_prior)!=k_regime){
    stop("dirichlet_prior should be a k_regime x k_regime matrix")
  }
  
  if(is.null(initDist)){
    initDist=as.double(rep(1/k_regime,k_regime))
  }else{
    #ensure that initdist sum to 1
    initDist=initDist/sum(initDist)
  }


  
  tmb_data<-list(yt=data$logRS,
    st=data$S, 
    logalpha_u = logalpha_limits[2],
    logalpha_l = logalpha_limits[1],
    Smax_u = Smax_limits[2],
    Smax_l = Smax_limits[1],
    initDist = initDist,
    alpha_dirichlet = dirichlet_prior,
    priors_flag = priors_flag,
    stan_flag = stan_flag,
    sig_p_sd = sig_p_sd,
    logsmax_p_mean = logsmax_p_mean,
    logsmax_p_sd = logsmax_p_sd  
  )

  if(is.null(ini_param)){
    magS <- 1/log10_ceiling(max(data$S))
    initlm<-lm(logRS~S, data=data)
    Smaxguess<-ifelse(initlm$coefficients[[2]]>0,(1/magS),-1/initlm$coefficients[[2]])
  }

  if(is.null(ini_param)){
    tmb_params <- list(        
      vec_logitalpha = rep(find_linit(logalpha_limits[2],logalpha_limits[1],max(initlm$coefficients[[1]],.5)),
            k_regime),
      vec_logitSmax = rep(find_linit(Smax_limits[2],Smax_limits[1],Smaxguess),k_regime),
      logsigma = log(.6),
      scalar_logitalpha = find_linit(logalpha_limits[2],logalpha_limits[1],max(initlm$coefficients[[1]],.5)),
      scalar_logitSmax = find_linit(Smax_limits[2],Smax_limits[1],Smaxguess),      
      qij_tran = matrix(0.1,nrow=k_regime,ncol=k_regime-1)          
    )  
  }else{
      tmb_params <- ini_param
  }

  if(tv.par == "a"){
    
    tmb_map$vec_logitSmax<-as.factor(rep(NA,k_regime))
    tmb_map$scalar_logitalpha<-as.factor(NA)

    tmb_data$options_flag=1
    
  }else if(tv.par=="b"){
  
    
    tmb_map$vec_logitalpha<-as.factor(rep(NA,k_regime))
    tmb_map$scalar_logitSmax<-as.factor(NA)

     tmb_data$options_flag=2
  }else if(tv.par=="both"){

    tmb_map$scalar_logitSmax<-as.factor(NA)
    tmb_map$scalar_logitalpha<-as.factor(NA)

     tmb_data$options_flag=0
    
  }else{
    stop(paste("tv.par",tv.par,"not recognized."))
  } 

  
  tmb_obj <- TMB::MakeADFun(
      data = tmb_data, parameters = tmb_params, map = tmb_map,
      DLL = "SR_HMM_all", silent = silent)

  
  tmb_mcmc <- tmbstan::tmbstan(tmb_obj, chains=chains,
              iter=iter, init="random",
              lower=lowlimit , upper=hightlimit,
              control = control,
              warmup = warmup )

   
  fit_summary <- summary(tmb_mcmc)

  posterior <- as.matrix(tmb_mcmc)
  
  
  postlist<-as.list(as.data.frame(t(posterior[,-ncol(posterior)])))
   

  r <- lapply(postlist,tmb_obj$report)
 
    
  sd_report <- TMB::sdreport(tmb_obj)
  conv <- get_convergence_diagnostics(sd_report)
  
  npar <- length(unlist(tmb_params))-length(unlist(tmb_map))
  nll <- tmb_obj$report()$nll
  
  tmb_obj$report()$pnll
 
  AICc  <- 2*nll + 2*npar +(2*npar*(npar+1)/(nrow(data)-npar-1))
  BIC  <- 2*nll + npar*log(nrow(data))

  structure(list(
    logalpha    = if(tv.par == "b"){tmb_obj$report()$logalpha[1]}else{tmb_obj$report()$logalpha},
    beta     = if(tv.par == "a"){tmb_obj$report()$beta[1]}else{tmb_obj$report()$beta},
    sigma      = tmb_obj$report()$sigma,
    qij      = tmb_obj$report()$qij,
    Smsy      = tmb_obj$report()$Smsy,
    Smax      = if(tv.par == "a"){tmb_obj$report()$Smax[1]}else{tmb_obj$report()$Smax},
    umsy      = if(tv.par == "b"){tmb_obj$report()$umsy[1]}else{tmb_obj$report()$umsy},
    probregime =  tmb_obj$report()$r_pred,
    regime =  apply(tmb_obj$report()$r_pred, 2,which.max),
    AICc       = AICc,
    BIC        = BIC,
    model      = tmb_opt,
    data       = data,
    tmb_data   = tmb_data,
    tmb_params = tmb_params,
    tmb_map    = tmb_map,
    tmb_obj    = tmb_obj,
    gradients  = conv$final_grads,
    bad_eig    = conv$bad_eig,
    conv_problem= conv$conv_problem,
    call       = match.call(expand.dots = TRUE),
    sd_report  = sd_report,
    class      = "SR_HMM_all")
    )

}





# END
#***********************************************************************************

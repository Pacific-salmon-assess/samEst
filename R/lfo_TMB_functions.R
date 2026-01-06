#TMB lfo function
# Dan greenberg and modified by Catarina for the package




#' Simple Ricker model estimated with TMB
#'
#' @param data A list or data frame containing Spawners (S) and log(Recruits/Spawners) (logRS) time series. 
#' @param model string. Which parameters are time-varying? Options are c('static','rw_a','rw_b',
#' 'rw_both', 'HMM', 'HMM_a', 'HMM_b')
#' @param L starting point for LFO-CV (minimum value is 10)
#' @param siglfo string. Indicating whether full variance should be used for lfo of models with random walks in parameters
#' "obs" incates that only observation variance is considered for lfo calculations, "total" indicates that sum of 
#' process and observation variances are used. Option valid only for 'logalpha' tv par. 
#'@param dirichlet_prior Prior for transition probability matrix in HMM models
#' k_regime x k_regime matrix. If NULL prior is set to matrix(1,nrow=k_regime,ncol=k_regime). 
#' @param priors_flag Integer, 1 priors are included in estimation model, 0 priors are not included.
#' @param sig_p_sd sd for half normal prior on sigma parameter. default is 1.
#' @param priorSmax Options are : "default" or "maxobsS", default will rely on fixed values provided by the user, 
#' "maxobsS" will be calculated based on the observed data. 
#' 
#' @returns vector of lfo by year
#' 
#' @importFrom stats nlminb 
#' 
#' @export
#' 
#' @examples
#' data(harck)
#' tmb_mod_lfo_cv(data=harck, model=c('static'))
#' 
tmb_mod_lfo_cv=function(data, model=c('static','staticAC','rw_a','rw_b','rw_both', 'HMM', 'HMM_a','HMM_b'), 
                        L=10, siglfo=c("obs","total")[1],dirichlet_prior=NULL,  priors_flag=1, sig_p_sd=1, 
                        priorSmax=c("default","maxobsS")[1],Smax_mean=exp(12),Smax_sd=exp(12) ){
  
  
  if(L<10){
    warning("L values < 10 are not recommended, results may be unstable")
  }
  if(sum(c("S","logRS") %in% names(data))!=2){
    stop("data must contain 'S' and 'logRS', please revise names(data)") 
  }
  
  conv_problem<-numeric(length(L:(nrow(data) - 1)))
  fail_conv<-numeric(length(L:(nrow(data) - 1)))
  
  if(model=='static'){
    exact_elpds_1b <- numeric(nrow(data)) #loglik for 1-year back estimates of productivity/capacity
    for (i in L:(nrow(data) - 1)) {
      past <- 1:i
      oos <- i + 1
      df_past <- data[past, , drop = FALSE]
      df_oos <- data[c(past, oos), , drop = FALSE]

      if(priorSmax=="maxobsS"){
        Smax_mean<-(max(df_past$S)*.5)
        Smax_sd<-Smax_mean
 
        #logb_p_sd=sqrt(log(1+((1/ Smax_sd)*(1/ Smax_sd))/((1/Smax_mean)*(1/Smax_mean))))
        #logb_p_mean=log(1/(Smax_mean))-0.5*logbeta_pr_sig^2
 
      }
      
      
      fit_past_tmb <-tryCatch({ricker_TMB(data=df_past,silent = TRUE,
                                priors_flag=priors_flag,sig_p_sd=sig_p_sd,
                                Smax_mean=Smax_mean,Smax_sd=Smax_sd)},
                                  error=function(cond){
                                    message(cond)
                                    return(list(fail_conv=1,
                                      conv_problem=1))}
                                  )
      conv_problem[i-(L-1)] <- ifelse(is.null(fit_past_tmb$model$convergence),1,fit_past_tmb$model$convergence)
      fail_conv[i-(L-1)] <- ifelse(is.null(fit_past_tmb$fail_conv),0,fit_past_tmb$fail_conv)

      if(fail_conv[i-(L-1)]==0){
        rs_pred_1b=fit_past_tmb$logalpha-fit_past_tmb$beta*df_oos$S[i + 1]
        exact_elpds_1b[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_1b,sd=fit_past_tmb$sigma))
      }else{
        exact_elpds_1b[i+1] <- NA
      }
    }
    exact_elpds_1b=exact_elpds_1b[-(1:L)]
    return(list(lastparam=exact_elpds_1b,
                conv_problem=conv_problem,
                fail_conv=fail_conv)
    )
  }else if(model=='staticAC'){
    exact_elpds_1b <- numeric(nrow(data)) #loglik for 1-year back estimates of productivity/capacity
    for (i in L:(nrow(data) - 1)) {
      past <- 1:i
      oos <- i + 1
      df_past <- data[past, , drop = FALSE]
      df_oos <- data[c(past, oos), , drop = FALSE]

       if(priorSmax=="maxobsS"){
        Smax_mean<-(max(df_past$S)*.5)
        Smax_sd<-Smax_mean
 
        #logb_p_sd=sqrt(log(1+((1/ Smax_sd)*(1/ Smax_sd))/((1/Smax_mean)*(1/Smax_mean))))
        #logb_p_mean=log(1/(Smax_mean))-0.5*logbeta_pr_sig^2
 
      }
      
      
      #fit_past_tmb<- ricker_TMB(data=df_past,AC=TRUE)
      fit_past_tmb <- tryCatch({ricker_TMB(data=df_past, AC=TRUE,silent = TRUE,
                                priors_flag=priors_flag,sig_p_sd=sig_p_sd,
                                Smax_mean=Smax_mean,Smax_sd=Smax_sd)},
                                  error=function(cond){
                                    message(cond)
                                    return(list(fail_conv=1,
                                      conv_problem=1))}
                                  )

      conv_problem[i-(L-1)] <- ifelse(is.null(fit_past_tmb$model$convergence),1,fit_past_tmb$model$convergence)
      fail_conv[i-(L-1)] <- ifelse(is.null(fit_past_tmb$fail_conv),0,fit_past_tmb$fail_conv)

      if(fail_conv[i-(L-1)]==0){
        rs_pred_1b<-fit_past_tmb$logalpha-fit_past_tmb$beta*df_oos$S[i + 1] + fit_past_tmb$residuals[i] * fit_past_tmb$rho
        exact_elpds_1b[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_1b,sd=fit_past_tmb$sigma))
      }else{
        exact_elpds_1b[i+1] <- NA
      }
    }
    exact_elpds_1b=exact_elpds_1b[-(1:L)]
    return(list(lastparam=exact_elpds_1b,
                conv_problem=conv_problem,
                fail_conv=fail_conv))
  }else if(model=='rw_a'){
    
    exact_elpds_1b <- numeric(nrow(data)) #loglik for 1-year back estimates of productivity/capacity
    exact_elpds_3b <- numeric(nrow(data)) #loglik for average of last 3-years of productivity/capacity
    exact_elpds_5b <- numeric(nrow(data)) #loglik for average of last 5-years of productivity/capacity
    for (i in L:(nrow(data) - 1)) {
      past <- 1:i
      oos <- i + 1
      df_past <- data[past, , drop = FALSE]
      df_oos <- data[c(past, oos), , drop = FALSE]

      if(priorSmax=="maxobsS"){
        Smax_mean<-(max(df_past$S)*.5)
        Smax_sd<-Smax_mean
 
        #logb_p_sd=sqrt(log(1+((1/ Smax_sd)*(1/ Smax_sd))/((1/Smax_mean)*(1/Smax_mean))))
        #logb_p_mean=log(1/(Smax_mean))-0.5*logbeta_pr_sig^2
 
      }
      
      
      fit_past_tv_a_tmb <- tryCatch({ricker_rw_TMB(data=df_past,tv.par='a',silent = TRUE,
                                priors_flag=priors_flag,sig_p_sd=sig_p_sd,
                                Smax_mean=Smax_mean,Smax_sd=Smax_sd)},
                                  error=function(cond){
                                    message(cond)
                                    return(list(fail_conv=1,
                                      conv_problem=1))}
                                  )
      conv_problem[i-(L-1)] <- ifelse(is.null(fit_past_tv_a_tmb$model$convergence),1,fit_past_tv_a_tmb$model$convergence)
      fail_conv[i-(L-1)] <- ifelse(is.null(fit_past_tv_a_tmb$fail_conv),0,fit_past_tv_a_tmb$fail_conv)
      

      if(fail_conv[i-(L-1)]==0){ 
        rs_pred_1b=fit_past_tv_a_tmb$logalpha[i]-fit_past_tv_a_tmb$beta[i]*df_oos$S[i + 1]
        rs_pred_3b=mean(fit_past_tv_a_tmb$logalpha[(i-2):i])-fit_past_tv_a_tmb$beta[i]*df_oos$S[i + 1]
        rs_pred_5b=mean(fit_past_tv_a_tmb$logalpha[(i-4):i])-fit_past_tv_a_tmb$beta[i]*df_oos$S[i + 1]
        
        if(siglfo=="obs"){
          exact_elpds_1b[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_1b,sd=fit_past_tv_a_tmb$sigma))
          exact_elpds_3b[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_3b,sd=fit_past_tv_a_tmb$sigma))
          exact_elpds_5b[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_5b,sd=fit_past_tv_a_tmb$sigma))
        }else if(siglfo=="total"){
          
          sigtot <-sqrt(fit_past_tv_a_tmb$sigma^2+fit_past_tv_a_tmb$sigma_a^2)
          exact_elpds_1b[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_1b,sd=sigtot))
          exact_elpds_3b[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_3b,sd=sigtot))
          exact_elpds_5b[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_5b,sd=sigtot))
        }else{
          stop("siglfo incorrectly defined options are `total` or `obs`")
        }
      }else{
        exact_elpds_1b[i+1] <- NA
        exact_elpds_3b[i+1] <- NA
        exact_elpds_5b[i+1] <- NA
      }
    }
    exact_elpds_1b=exact_elpds_1b[-(1:L)]
    exact_elpds_3b=exact_elpds_3b[-(1:L)]
    exact_elpds_5b=exact_elpds_5b[-(1:L)]
    return(list(lastparam=exact_elpds_1b,
                last3paramavg=exact_elpds_3b,
                last5paramavg=exact_elpds_5b,
                conv_problem=conv_problem,
                fail_conv=fail_conv))
  }else if(model=='rw_b'){
    #stop("not defined")
    exact_elpds_1b <- numeric(nrow(data)) #loglik for 1-year back estimates of productivity/capacity
    exact_elpds_3b <- numeric(nrow(data)) #loglik for average of last 3-years of productivity/capacity
    exact_elpds_5b <- numeric(nrow(data))#loglik for average of last 5-years of productivity/capacity
    for (i in L:(nrow(data) - 1)) {
      past <- 1:i
      oos <- i + 1
      df_past <- data[past, , drop = FALSE]
      df_oos <- data[c(past, oos), , drop = FALSE]
      
      if(priorSmax=="maxobsS"){
        Smax_mean<-(max(df_past$S)*.5)
        Smax_sd<-Smax_mean
 
      }
      
      fit_past_tv_b_tmb <- tryCatch({ricker_rw_TMB(data=df_past,tv.par='b',silent = TRUE,
                                priors_flag=priors_flag,sig_p_sd=sig_p_sd,
                                Smax_mean=Smax_mean,Smax_sd=Smax_sd)},
                                  error=function(cond){
                                    message(cond)
                                    return(list(fail_conv=1,
                                      conv_problem=1))}
                                  )
    
      conv_problem[i-(L-1)] <- ifelse(is.null(fit_past_tv_b_tmb$model$convergence),1,fit_past_tv_b_tmb$model$convergence)
      fail_conv[i-(L-1)] <- ifelse(is.null(fit_past_tv_b_tmb$fail_conv),0,fit_past_tv_b_tmb$fail_conv)
      
      
      if(fail_conv[i-(L-1)]==0){ 
        rs_pred_1b=fit_past_tv_b_tmb$logalpha[i]-fit_past_tv_b_tmb$beta[i]*df_oos$S[i + 1]
        rs_pred_3b=fit_past_tv_b_tmb$logalpha[i]-mean(fit_past_tv_b_tmb$beta[(i-2):i])*df_oos$S[i + 1]
        rs_pred_5b=fit_past_tv_b_tmb$logalpha[i]-mean(fit_past_tv_b_tmb$beta[(i-4):i])*df_oos$S[i + 1]
        
        if(siglfo=="obs"){
          exact_elpds_1b[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_1b,sd=fit_past_tv_b_tmb$sigma))
          exact_elpds_3b[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_3b,sd=fit_past_tv_b_tmb$sigma))
          exact_elpds_5b[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_5b,sd=fit_past_tv_b_tmb$sigma))
        }else if(siglfo=="total"){
          
          sigtot <-sqrt(fit_past_tv_b_tmb$sigma^2+fit_past_tv_b_tmb$sigma_b^2)
          exact_elpds_1b[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_1b,sd=sigtot))
          exact_elpds_3b[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_3b,sd=sigtot))
          exact_elpds_5b[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_5b,sd=sigtot))
        }else{
          stop("siglfo incorrectly defined, options are `total` or `obs`")
        }
      }else{
        exact_elpds_1b[i+1] <- NA
        exact_elpds_3b[i+1] <- NA
        exact_elpds_5b[i+1] <- NA
      }
      
      
    }
    exact_elpds_1b=exact_elpds_1b[-(1:L)]
    exact_elpds_3b=exact_elpds_3b[-(1:L)]
    exact_elpds_5b=exact_elpds_5b[-(1:L)]
    return(list(lastparam=exact_elpds_1b,
                last3paramavg=exact_elpds_3b,
                last5paramavg=exact_elpds_5b,
                conv_problem=conv_problem,
                fail_conv=fail_conv))
  }else if(model=='rw_both'){
    
    exact_elpds_1b <- numeric(nrow(data)) #loglik for 1-year back estimates of productivity/capacity
    exact_elpds_3b <- numeric(nrow(data)) #loglik for average of last 3-years of productivity/capacity
    exact_elpds_5b <- numeric(nrow(data)) #loglik for average of last 5-years of productivity/capacity
    for (i in L:(nrow(data) - 1)) {
      past <- 1:i
      oos <- i + 1
      df_past <- data[past, , drop = FALSE]
      df_oos <- data[c(past, oos), , drop = FALSE]
     
      if(priorSmax=="maxobsS"){
        Smax_mean<-(max(df_past$S)*.5)
        Smax_sd<-Smax_mean
 
        
 
      }
       
     
      fit_past_tv_ab_tmb <- tryCatch({ricker_rw_TMB(data=df_past,tv.par='both',silent = TRUE,
                                priors_flag=priors_flag,sig_p_sd=sig_p_sd,
                                Smax_mean=Smax_mean,Smax_sd=Smax_sd)},
                                  error=function(cond){
                                    message(cond)
                                    return(list(fail_conv=1,
                                      conv_problem=1))}
                                  )

      conv_problem[i-(L-1)] <- ifelse(is.null(fit_past_tv_ab_tmb$model$convergence),1,fit_past_tv_ab_tmb$model$convergence)
      fail_conv[i-(L-1)] <- ifelse(is.null(fit_past_tv_ab_tmb$fail_conv),0,fit_past_tv_ab_tmb$fail_conv)
    
      if(fail_conv[i-(L-1)]==0){ 
        rs_pred_1b=fit_past_tv_ab_tmb$logalpha[i]-fit_past_tv_ab_tmb$beta[i]*df_oos$S[i + 1]
        rs_pred_3b=mean(fit_past_tv_ab_tmb$logalpha[(i-2):i])-mean(fit_past_tv_ab_tmb$beta[(i-2):i])*df_oos$S[i + 1]
        rs_pred_5b=mean(fit_past_tv_ab_tmb$logalpha[(i-4):i])-mean(fit_past_tv_ab_tmb$beta[(i-4):i])*df_oos$S[i + 1]
        
        exact_elpds_1b[i+1] <- log(dnorm(df_oos$logRS[i],mean=rs_pred_1b,sd=exp(fit_past_tv_ab_tmb$sigma)))
        exact_elpds_3b[i+1] <- log(dnorm(df_oos$logRS[i],mean=rs_pred_3b,sd=exp(fit_past_tv_ab_tmb$sigma)))
        exact_elpds_5b[i+1] <- log(dnorm(df_oos$logRS[i],mean=rs_pred_5b,sd=exp(fit_past_tv_ab_tmb$sigma)))
      }else{
        exact_elpds_1b[i+1] <- NA
        exact_elpds_3b[i+1] <- NA
        exact_elpds_5b[i+1] <- NA
      }
    }
    exact_elpds_1b=exact_elpds_1b[-(1:L)]
    exact_elpds_3b=exact_elpds_3b[-(1:L)]
    exact_elpds_5b=exact_elpds_5b[-(1:L)]
    return(list(lastparam=exact_elpds_1b,
                last3paramavg=exact_elpds_3b,
                last5paramavg=exact_elpds_5b,
                conv_problem=conv_problem,
                fail_conv=fail_conv))
    
  }else if(model=='HMM'){
    #stop("not defined")
    exact_elpds_1k <- numeric(nrow(data)) #loglik choosing a specific regime in a given year
    exact_elpds_3k  <- numeric(nrow(data))
    exact_elpds_5k  <- numeric(nrow(data))
    
    for (i in L:(nrow(data) - 1)) {
      past <- 1:i
      oos <- i + 1
      df_past <- data[past, , drop = FALSE]
      df_oos <- data[c(past, oos), , drop = FALSE]
      
      if(priorSmax=="maxobsS"){
        Smax_mean<-(max(df_past$S)*.5)
        Smax_sd<-Smax_mean
 
 
      }
      
      
      fit_past_hmm_tmb <- tryCatch({ricker_hmm_TMB(data=df_past,tv.par='both',silent = TRUE,
                                          dirichlet_prior=dirichlet_prior,
                                priors_flag=priors_flag,sig_p_sd=sig_p_sd,
                                Smax_mean=Smax_mean,Smax_sd=Smax_sd)},
                                  error=function(cond){
                                    message(cond)
                                    return(list(fail_conv=1,
                                      conv_problem=1))}
                                  )

      conv_problem[i-(L-1)] <- ifelse(is.null(fit_past_hmm_tmb$model$convergence),1,fit_past_hmm_tmb$model$convergence)
      fail_conv[i-(L-1)] <- ifelse(is.null(fit_past_hmm_tmb$fail_conv),0,fit_past_hmm_tmb$fail_conv)
      
      if(fail_conv[i-(L-1)]==0){
        logalpha <- fit_past_hmm_tmb$logalpha[fit_past_hmm_tmb$regime]
        beta <- fit_past_hmm_tmb$beta[fit_past_hmm_tmb$regime]
        sigma <- fit_past_hmm_tmb$sigma
      
        rs_pred_1k=logalpha[i]-beta[i]*df_oos$S[i + 1]
        rs_pred_3k=mean(logalpha[(i-2):i])-mean(beta[(i-2):i])*df_oos$S[i + 1]
        rs_pred_5k=mean(logalpha[(i-4):i])-mean(beta[(i-4):i])*df_oos$S[i + 1]
      
      
        exact_elpds_1k[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_1k,sd=sigma))
        exact_elpds_3k[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_3k,sd=sigma))
        exact_elpds_5k[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_5k,sd=sigma))
      }else{
        exact_elpds_1k[i+1] <- NA
        exact_elpds_3k[i+1] <- NA
        exact_elpds_5k[i+1] <- NA
      }
    }
    exact_elpds_1k=exact_elpds_1k[-(1:L)]
    exact_elpds_3k=exact_elpds_3k[-(1:L)]
    exact_elpds_5k=exact_elpds_5k[-(1:L)]
    
    return(list(lastregime_pick=exact_elpds_1k, 
                last3regime_pick=exact_elpds_3k, 
                last5regime_pick=exact_elpds_5k,
                conv_problem=conv_problem,
                fail_conv=fail_conv))
    
  }else if(model=='HMM_a'){
    #stop("not defined")
    exact_elpds_1k <- numeric(nrow(data)) #loglik choosing a specific regime in a given year
    exact_elpds_3k  <- numeric(nrow(data))
    exact_elpds_5k  <- numeric(nrow(data))
    
    for (i in L:(nrow(data) - 1)) {
      past <- 1:i
      oos <- i + 1
      df_past <- data[past, , drop = FALSE]
      df_oos <- data[c(past, oos), , drop = FALSE]
      
      if(priorSmax=="maxobsS"){
        Smax_mean<-(max(df_past$S)*.5)
        Smax_sd<-Smax_mean
 
      
      }
      

      fit_past_hmm_tmb<- tryCatch({ricker_hmm_TMB(data=df_past,tv.par='a',silent = TRUE,
                                            dirichlet_prior=dirichlet_prior,
                                priors_flag=priors_flag,sig_p_sd=sig_p_sd,
                                Smax_mean=Smax_mean,Smax_sd=Smax_sd)},
                                  error=function(cond){
                                    message(cond)
                                    return(list(fail_conv=1,
                                      conv_problem=1))}
                                  )

      conv_problem[i-(L-1)] <- ifelse(is.null(fit_past_hmm_tmb$model$convergence),1,fit_past_hmm_tmb$model$convergence)
      fail_conv[i-(L-1)] <- ifelse(is.null(fit_past_hmm_tmb$fail_conv),0,fit_past_hmm_tmb$fail_conv)
      
      if(fail_conv[i-(L-1)]==0){
        logalpha <- fit_past_hmm_tmb$logalpha[fit_past_hmm_tmb$regime]
        beta <- fit_past_hmm_tmb$beta
        sigma <- fit_past_hmm_tmb$sigma
      
        rs_pred_1k=logalpha[i]-beta*df_oos$S[i + 1]
        rs_pred_3k=mean(logalpha[(i-2):i])-mean(beta)*df_oos$S[i + 1]
        rs_pred_5k=mean(logalpha[(i-4):i])-mean(beta)*df_oos$S[i + 1]
      
        exact_elpds_1k[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_1k,sd=sigma))
        exact_elpds_3k[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_3k,sd=sigma))
        exact_elpds_5k[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_5k,sd=sigma))
      }else{
        exact_elpds_1k[i+1] <- NA
        exact_elpds_3k[i+1] <- NA
        exact_elpds_5k[i+1] <- NA
      }
      
    }
    exact_elpds_1k=exact_elpds_1k[-(1:L)]
    exact_elpds_3k=exact_elpds_3k[-(1:L)]
    exact_elpds_5k=exact_elpds_5k[-(1:L)]
    
    return(list(lastregime_pick=exact_elpds_1k, 
                last3regime_pick=exact_elpds_3k, 
                last5regime_pick=exact_elpds_5k,
                conv_problem=conv_problem,
                fail_conv=fail_conv))
    
  }else if(model=='HMM_b'){
    #stop("not defined")
    exact_elpds_1k <- numeric(nrow(data)) #loglik choosing a specific regime in a given year
    exact_elpds_3k  <- numeric(nrow(data))
    exact_elpds_5k  <- numeric(nrow(data))
    
    for (i in L:(nrow(data) - 1)) {
      past <- 1:i
      oos <- i + 1
      df_past <- data[past, , drop = FALSE]
      df_oos <- data[c(past, oos), , drop = FALSE]
      
      if(priorSmax=="maxobsS"){
        Smax_mean<-(max(df_past$S)*.5)
        Smax_sd<-Smax_mean
 
      }
      

      fit_past_hmm_tmb <- tryCatch({ricker_hmm_TMB(data=df_past,tv.par='b',silent = TRUE,
                                        dirichlet_prior=dirichlet_prior,
                                priors_flag=priors_flag,sig_p_sd=sig_p_sd,
                                Smax_mean=Smax_mean,Smax_sd=Smax_sd)},
                                  error=function(cond){
                                    message(cond)
                                    return(list(fail_conv=1,
                                      conv_problem=1))}
                                  )
      
      conv_problem[i-(L-1)] <- ifelse(is.null(fit_past_hmm_tmb$model$convergence),1,fit_past_hmm_tmb$model$convergence)
      fail_conv[i-(L-1)] <- ifelse(is.null(fit_past_hmm_tmb$fail_conv),0,fit_past_hmm_tmb$fail_conv)
      
      if(fail_conv[i-(L-1)]==0){
        logalpha <- fit_past_hmm_tmb$logalpha
        beta <- fit_past_hmm_tmb$beta[fit_past_hmm_tmb$regime]
        sigma <- fit_past_hmm_tmb$sigma
      
        rs_pred_1k<-logalpha-beta[i]*df_oos$S[i + 1]
        rs_pred_3k=mean(logalpha)-mean(beta[(i-2):i])*df_oos$S[i + 1]
        rs_pred_5k=mean(logalpha)-mean(beta[(i-4):i])*df_oos$S[i + 1]
        
        exact_elpds_1k[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_1k,sd=sigma))
        exact_elpds_3k[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_3k,sd=sigma))
        exact_elpds_5k[i+1] <- log(dnorm(df_oos$logRS[i+1],mean=rs_pred_5k,sd=sigma))
      }else{
        exact_elpds_1k[i+1] <- NA
        exact_elpds_3k[i+1] <- NA
        exact_elpds_5k[i+1] <- NA
      }
    }
    exact_elpds_1k=exact_elpds_1k[-(1:L)]
    exact_elpds_3k=exact_elpds_3k[-(1:L)]
    exact_elpds_5k=exact_elpds_5k[-(1:L)]
    
    return(list(lastregime_pick=exact_elpds_1k, 
                last3regime_pick=exact_elpds_3k, 
                last5regime_pick=exact_elpds_5k,
                conv_problem=conv_problem,
                fail_conv=fail_conv))
    
  }else{
    stop(paste("model", model,"not defined, valid options are c('static','rw_a','rw_b','rw_both', 'HMM', 'HMM_a','HMM_b')"))
  }
}

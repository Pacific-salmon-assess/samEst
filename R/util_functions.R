# util functions
#===============================================
#Structure of the code copied from the sdmTMB package :
#https://github.com/pbs-assess/sdmTMB


#' Find order of magnitude
#'
#' @param x a number. 
#' @export
#' 
#' @returns order of magnitude
#' 
#' 
log10_ceiling <- function(x) {
    10^(ceiling(log10(x)))
}



#' Find lalpha intial value for logistic transformed TMB HMM model parameters
#'
#' @param U a number. upper limit for parameter 
#' @param L a number. lower limit for parameter 
#' @param initpar Inital value of parameter to be transformed 
#' 
#' @export
#' 
#' @returns logistic transformed parameter guess
#' 
#' 
find_linit <- function(U, L, initpar) {
    linitpar <- -log((U-L)/(initpar-L)-1)

    return( linitpar)
}

#'log sum of exponentials
#'
#' @param x a vector 
#' @export
#' 
#' @returns log sum of exponentials
#' 
#' 
log_sum_exp <- function(x) {
  max_x <- max(x)  
  max_x + log(sum(exp(x - max_x)))
}

#' mean of log-sum-exp
#'
#' @param x a vector 
#' @export
#' 
#' @returns mean of exponentials
#' 
#' 
log_mean_exp <- function(x) {
  log_sum_exp(x) - log(length(x))
}

#' Function to format the names for a stock-recruitment dataset to avoid compatibility problems with other samEst functions
#'
#' @param S a vector of spawners
#' @param R a vector of recruits
#' @param R a vector of years for each cohort
#' @export
#' 
#' @returns mean of exponentials
#' 
#' 
sr_format<- function(S,R,by){
  x=data.frame(S=S,R=R,by=by,logRS=log(R/S))
  return(x)
}



#' function to trasnform parameter between minus one an one
#' @param x a value 
#'
#' @export
#' 
#' @returns inverse logit transformed variable 
#' 
#' @importFrom Rdpack reprompt
#' 
#' 
#' 
minus_one_to_one <- function(x) {
  return(2 * (exp(x)/(1+exp(x))) - 1)

}





#' function to generate lognormal prior parameters based on mean and standard deviation in normal 
#' space. Based on Quan and Zhang (2003) Estimate of standard deviation for a log-transformed 
#' variable using arithmetic means and standard deviations. Statist. Med. 22:2723–2736
#' 
#' @param Smax_mean 
#' @param Smax_sd
#'
#' @export
#' 
#' @returns list of mean and sd for lognormal prior on beta 
#' 
#' 
log_prior_params <- function(Smax_mean,Smax_sd) {
  
  logsmax_pr_sig = sqrt(log(1+(Smax_sd^2)/(Smax_mean^2)))
  logsmax_pr = log(Smax_mean)-0.5*logsmax_pr_sig^2

  return(list(logsmax_pr_mean=logsmax_pr,
             logsmax_pr_sig=logsmax_pr_sig))
}






# END
#***********************************************************************************

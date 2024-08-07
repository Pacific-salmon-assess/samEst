% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TMB_functions.R
\name{ricker_rw_TMB_deprecated}
\alias{ricker_rw_TMB_deprecated}
\title{deprecated version based on no centered model Ricker model with random walk in a, b or both parameters with TMB}
\usage{
ricker_rw_TMB_deprecated(
  data,
  tv.par = c("a", "b", "both"),
  silent = FALSE,
  control = TMBcontrol(),
  ini_param = NULL,
  tmb_map = list(),
  priors_flag = 1,
  stan_flag = 0,
  sig_p_sd = 1,
  siga_p_sd = 1,
  sigb_p_sd = 0.3,
  logb_p_mean = -12,
  logb_p_sd = 3,
  AICc_type = c("conditional", "marginal")[1],
  deltaEDF = 1e-04,
  newton_stp = TRUE,
  useEDF = FALSE
)
}
\arguments{
\item{data}{A list or data frame containing Spawners (S) and log(Recruits/Spawners) (logRS) time series.}

\item{tv.par}{Which parameters should vary? Either productivity (intercept, a), capacity (slope, b) or both parameters}

\item{silent}{Logical Silent or optimization details? default is FALSE}

\item{control}{output from TMBcontrol() function, to be passed to nlminb()}

\item{ini_param}{Optional. A list with initial parameter guesses.}

\item{tmb_map}{optional, mapping list indicating if parameters should be estimated of fixed.
Default is all parameters are estimated}

\item{priors_flag}{Integer, 1 priors are included in estimation model, 0 priors are not included.
See details for priors documentation. See details for priors documentation.}

\item{stan_flag}{Integer, flag indicating wether or not TMB code will be used with TMBstan - Jacobian
adjustment implemented. Default is 0, jacobian adjustment not included.}

\item{sig_p_sd}{sd for half normal prior on sigma parameter. default is 1.}

\item{siga_p_sd}{sd for half normal prior on sigma for logalpha random walk parameter. default is 1.}

\item{sigb_p_sd}{sd for half normal prior on sigma for beta random walk parameter. default is 1.}

\item{logb_p_mean}{mean for prior on log b, default is -12.}

\item{logb_p_sd}{sd for prior on log b, default is 3.}

\item{AICc_type}{"conditional" (Thorson 2024) or "marginal"}

\item{newton_stp}{logical. use newton step to improve estimates and gradients}

\item{useEDF}{logical use EDF algorithm describef in Thorson 2024}
}
\value{
a list containing several model outputs:
\itemize{
\item logalpha - MLE estimates for the alpha parameter vector
\item beta - MLE estimates for the beta parameter
\item sig - MLE estimates for the observation error standard deviation
\item siga - MLE estimates for the process error (variation in logalpha) standard deviation
\item sigb - MLE estimates for the process error (variation in beta) standard deviation
\item AICc - AICc values, given by 2\emph{nll + 2}npar +(2\emph{npar}(npar+1)/(nrow(data)-npar-1)), excluding prior components
\item BIC - BIC values, excluding prior components
\item model - opt object, generated by the \code{stats::nlminb()}
\item tmb_data - data provided to TMB,
\item tmb_params - parameter intial guesses,
\item tmb_map - mapping indicating which parameters should be fixed or estimated,
\item tmb_obj - obj object, generated by the \code{TMB::MakeADFun()}
\item gradients - final convergence gradient
\item bad_eig - eigen value
\item call - original function call
\item sd_report - MLE estimates and sdt Error estimates for main parameter estimates
\item class - name of cpp model
}
}
\description{
deprecated version based on no centered model Ricker model with random walk in a, b or both parameters with TMB
}
\details{
Priors: Weakly informative priors are included for the main parameterst of the model:
logalpha ~ N(1.5,2.5)
logbeta ~ N(-12,3)
sigobs ~ gamma(2,1/3)
siga ~ gamma(2,1/3)
sigb ~ gamma(2,1/3)
}
\examples{
data(harck)
ricker_rwa_TMB(data=harck)


}

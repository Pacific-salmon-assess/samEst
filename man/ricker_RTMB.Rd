% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RTMB_run_functions.R
\name{ricker_RTMB}
\alias{ricker_RTMB}
\title{Simple Ricker model estimated with RTMB}
\usage{
ricker_RTMB(
  data,
  silent = FALSE,
  control = TMBcontrol(),
  tmb_map = list(),
  AC = FALSE,
  priors_flag = 1,
  stan_flag = 0,
  sig_p_sd = 1,
  logb_p_mean = -12,
  logb_p_sd = 3,
  dautoreg = 0
)
}
\arguments{
\item{data}{A list or data frame containing Spawners (S), log(Recruits/Spawners) (logRS), and brood year data (by).}

\item{silent}{Logical Silent or optimization details? default is FALSE}

\item{control}{output from TMBcontrol() function, to be passed to nlminb()}

\item{tmb_map}{optional, mapping list indicating if parameters should be estimated of fixed.
Default is all parameters are estimated}

\item{AC}{Logical. Are residuals autocorrelated? Default is FALSE}

\item{priors_flag}{Integer, 1 priors are included in estimation model, 0 priors are not included.
See details for priors documentation.}

\item{stan_flag}{Integer, flag indicating wether or not TMB code will be used with TMBstan - Jacobian
adjustment implemented. Only works if priors flag is also 1, Default is 0, jacobian adjustment not included.}

\item{sig_p_sd}{sd for half normal prior on sigma parameter. default is 1.}

\item{logb_p_mean}{mean for prior on log b, default is -12.}

\item{logb_p_sd}{sd for prior on log b, default is 3.#'}

\item{dautoreg}{select model code for AR1}
}
\value{
penalized (if priors_flag= TRUE) negative log likelihood
}
\description{
Simple Ricker model estimated with RTMB
}
\details{
this function will
}
\examples{
data(harck)
rickerTMB(data=harck)

}

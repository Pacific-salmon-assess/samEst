% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lfo_stan_functions.R
\name{stan_refit2}
\alias{stan_refit2}
\title{stan_refit function}
\usage{
stan_refit2(
  mod,
  newdata,
  oos,
  K = 2,
  dirichlet_prior = NULL,
  pSmax_mean = NULL,
  pSmax_sig = NULL,
  psig_b = NULL
)
}
\arguments{
\item{newdata}{new data to feed to refit models.}

\item{oos}{What row of the new data should be treated as the out-of-sample test}

\item{K}{Number of potential regime shifts - default 2}

\item{dirichlet_prior}{k_regime x k_regime matrix. Prior for transition probability matrix,
if NULL (default), prior is set to matrix(1,nrow=k_regime,ncol=k_regime)}

\item{pSmax_mean}{mean for Smax prior}

\item{pSmax_sig}{sd for Smax prior}

\item{psig_b}{sd for sigma_b prior}

\item{sm}{Compiled stan model, which can be generated with sr_mod function.}

\item{regime}{TRUE or FALSE statement - is this a regime shift model or not}
}
\value{
returns the model fit
}
\description{
This function refits a stan model (via rstan) to iteratively estimate the out-of-sample likelihood
}
\examples{
r=stan_refit(sm=mod3,newdata=df,oos=12)
}

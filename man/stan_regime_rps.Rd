% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RP_functions.R
\name{stan_regime_rps}
\alias{stan_regime_rps}
\title{Compute reference points from parameters of Hidden Markov Model (regime shift) Stan model}
\usage{
stan_regime_rps(m, par = c("a", "b", "both"), lambertW = FALSE)
}
\arguments{
\item{m}{estimated fit for rstan model}
}
\value{
Estimates for Smax, Smsy, and Umsy over time, either unweighted (ie. returns parameters for the most likely regime sequence) or weighted (ie. probability of each regime x regime parameters)
}
\description{
Compute reference points from parameters of Hidden Markov Model (regime shift) Stan model
}

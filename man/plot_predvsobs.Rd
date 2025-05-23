% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_functions.R
\name{plot_predvsobs}
\alias{plot_predvsobs}
\title{Predicted vs observed plot}
\usage{
plot_predvsobs(
  post.obj = NULL,
  main.title = "Stock - Predicted vs. Observed",
  plot.scale = 10^6,
  scale.label = "(Millions)",
  rec.label = "Adult Recruits",
  spn.label = "Spawners",
  file.name = "Plot_FittedvsObs.png"
)
}
\arguments{
\item{post.obj}{posterior samples}

\item{main.title}{plot title}

\item{plot.scale}{plot scale - how to rescale spawners (e.g. thousands, millions, etc.)}

\item{scale.label}{plot scale label}

\item{rec.label}{y-axis label - default 'Adult Recruits'}

\item{spn.label}{x-axis label - default 'Spawners'}

\item{file.name}{filename for output}
}
\description{
Predicted vs observed plot
}

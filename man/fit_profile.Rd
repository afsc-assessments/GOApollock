% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tmb_fns.R
\name{fit_profile}
\alias{fit_profile}
\title{Run and plot likelihood profiles}
\usage{
fit_profile(
  fits,
  xseq,
  par = c("M", "q2", "q6"),
  plot = TRUE,
  maxNLL = 1,
  ylim = c(0, 10),
  estSigR = FALSE,
  return.fits = FALSE
)
}
\arguments{
\item{fits}{A list of fits}

\item{xseq}{The vector over which to run the profile}

\item{par}{Which parameter, currently supports "M", "q2" and "q6".}

\item{plot}{Whether to plot it}

\item{maxNLL}{The maximum NLL below which to filter out individual components, for decluttering the plots.}

\item{ylim}{Y limit for plots}

\item{estSigR}{Whether to estimate sigmaR}

\item{return.fits}{Whether to return the fits (not implemented yet)}
}
\value{
A ggplot object of the profiles
}
\description{
Run and plot likelihood profiles
}

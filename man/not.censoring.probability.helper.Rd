% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper.R
\name{not.censoring.probability.helper}
\alias{not.censoring.probability.helper}
\title{not.censoring.probability.helper}
\usage{
not.censoring.probability.helper(eta, cumulative.hazard, censor.density, time)
}
\arguments{
\item{eta}{vector of linear predictors in a fitted Cox model.}

\item{cumulative.hazard}{vector of baseline cumulative hazard at \code{time}.}

\item{censor.density}{vector of the density of censoring times evaluated at \code{time}.}

\item{time}{vector of survival times.}
}
\value{
vector of the probabilities of not getting censored.
}
\description{
This function is a helper for the function \code{not.censoring.probability}.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DINA.R
\name{DINAEstimator.step2}
\alias{DINAEstimator.step2}
\title{DINAEstimator.step2}
\usage{
DINAEstimator.step2(X, Y, W, family, nuisance.function)
}
\arguments{
\item{X}{matrix of covariates.}

\item{Y}{vector of responses.}

\item{W}{vector of treatment assignment.}

\item{nuisance.function}{nuisance function estimators. If we call nuisance.function.estimator = nuisance.function(X), then nuisance.function.estimator is a list with two members: \code{a} denotes the estimated nuisance function a(x), \code{nu} denotes the estimated nuisance function nu(x).}
}
\description{
This function implements the second step of DINA: estimation of the causal parameter given nuisance function estimators.
}

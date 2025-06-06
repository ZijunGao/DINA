% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper.R
\name{split.weight}
\alias{split.weight}
\title{split.weight}
\usage{
\method{split}{weight}(n, k = 2)
}
\arguments{
\item{n}{sample size.}

\item{k}{number of folds. Default is 2.}
}
\value{
a list indicating the resulting folds.
 \itemize{
 \item{index.resample}{vector of subject indices in the bootstrap resample. For example, assume there are 3 subjects in the dataset and consider a bootstrap resample consisting of two copies of the first subject and one copy of the second subject. In this way, \code{index.resample} is the vector [1,1,2].}
 \item{foldIndex}{list of fold memberships. The j-th member is a vector containing all subjects of the bootstrap resample assigned to the j-th fold.}
 \item{weight}{vector of weights for the original dataset. The j-th coordinate stands for the times that the j-th subject of the original dataset appears in the bootstrap resample.}
}
}
\description{
This function splits a sample with duplicates into non-overlap folds with equal size.
}
\examples{
:
set.seed(318)
test = split.weight(n, k)
test$index.resample[test$foldIndex[[1]]]
test$index.resample[test$foldIndex[[2]]]
test$weight
}

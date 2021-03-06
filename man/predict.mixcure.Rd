% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mixcure.r
\name{predict.mixcure}
\alias{predict.mixcure}
\title{Prediction method for mixture cure models}
\usage{
\method{predict}{mixcure}(object, newdata, times, ...)
}
\arguments{
\item{object}{an object of mixcure}

\item{newdata}{a data frame containing covariate values at which the
survival and cure rate will be estimated. It should contain all the
covariates that are used to build object. Prediction will be made for each
row in newdata.}

\item{times}{a vector of times at which the survival probabilities are
estimated}

\item{...}{for compatibility purpose. Not used.}
}
\value{
a list with the following components:
  \item{cure}{a matrix of 2 columns and the same number of rows as
  newdata. The first column is uncure rates for the rows in newdata and
  the second column is cure rates for the rows in newdata}
  \item{uncuresurv}{a list with the number of components equal to the rows
  of newdata. Each component is a vector of the estimated survival
  probabilities at times for a subject if uncured}
  \item{surv}{similar to uncuresurv except that the survival probabilities
  are the unconditional survival probabilities}
  \item{times}{a vector of times at which the survival probabilities will be
  predicted}
}
\description{
This function computes the estimated survival probabilities and cure rates
for given sets of covariate values using the fitted model from mixcure().
Each set of the covariates values is stored in one row of newdata. newdata
must be a data frame containing all the covariates used in mixcure(). A
used-supplied set of times at which the survival probabilities will be
estimated must be provided in times.
}
\examples{

data(leukaemia)
predict(mixcure(Surv(time, cens) ~ transplant, ~ transplant, data = leukaemia),
newdata = leukaemia[1, ], times = 0:2000)

}
\seealso{
mixcure
}
\author{
Yingwei Peng
}

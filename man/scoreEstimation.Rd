% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gradientscore.R
\name{scoreEstimation}
\alias{scoreEstimation}
\title{Evaluate gradient score function}
\usage{
scoreEstimation(dat, model, weightFun, dWeightFun, alpha, rho, u, p = 1)
}
\arguments{
\item{dat}{\code{n} by \code{d} matrix of data}

\item{model}{string indicating whether the model is \code{dir} or \code{negdir}}

\item{weightFun}{object of class function, must be differentiable}

\item{dWeightFun}{object of class function, gradient of weightfun}

\item{alpha}{numeric vector of positive parameters}

\item{rho}{numeric parameter}

\item{u}{threshold}

\item{p}{norm for the risk function defined in weightFun (default is \eqn{l_p} norm)}
}
\value{
estimated score value
}
\description{
Evaluate gradient score function
}


% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitcomposite.R
\name{margexc}
\alias{margexc}
\title{Negative log-likelihood for marginal threshold exceedances}
\usage{
margexc(margdat, u, scale, shape)
}
\arguments{
\item{margdat}{an \code{m} by \code{d} matrix containing marginal exceedances, one per row}

\item{u}{a \code{d}-vector of thresholds}

\item{scale}{vector of scale parameters}

\item{shape}{vector of shape parameters}
}
\value{
the negative log-likelihood
}
\description{
This function handles missing values in the data matrix that create
marginal only non-missingness. These are not accounted for in the pairwise
likelihood routine and so are handled outside of the C++ code in a uniform
fashion for each of the sub-models
}


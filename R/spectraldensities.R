
#' Spectral density of scaled Dirichlet and negative Dirichlet extreme value distributions
#'
#' @param param vector containing the parameters \eqn{\alpha} and \eqn{\rho}, in this order.
#' @param dat matrix of Fréchet or Pareto observations, of dimension \eqn{n} by \eqn{d}.
#' @param d dimension of model
#' @param model string indicating which model to select from \code{"dir"}, \code{"negdir"} and \code{"ct"}
#' @param transform logical indicating whether parameters are on the log scale. Default to \code{TRUE}
#'
#' @details The function is provided as a wrapper and takes parameters on the log scale for \eqn{\alpha} (\eqn{\rho}).
#' @export
#' @return the log-density for the \code{n} sample
specdens <- function(param, dat, d, model=c("dir","negdir","ct"), transform=TRUE){
  mod <- match.arg(model,choices = c("dir","negdir","ct"))[1]
  switch(mod, dir=.dirspecdens(param=param, dat=dat, d=d, transform = transform),
              negdir=.negdirspecdens(param=param, dat=dat, d=d, transform = transform),
              ct = .ctspecdens(param=param,dat=dat, transform = transform))
}

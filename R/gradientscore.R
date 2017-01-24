#### Gradient score estimation for the scaled Dirichlet and negative Dirichlet extremal models

#lp norm
lpnorm <- function(x, p){
  if(is.null(dim(x)))
    { #vector
  exp(log(sum(x^p))/p)
  } else{ #matrix of observation
    exp(log(rowSums(x^p))/p)
  }
}
#Derivative of lp norm
gradlpnorm <- function(x, p){
  if(is.null(dim(x)))
  { #vector
    exp((1/p-1)*log(sum(x^p)))*x^(p-1)
  } else{ #matrix of observation
    exp((1/p-1)*log(rowSums(x^p)))*x^(p-1)
  }
}

#Weighting function
weightFun <- function(x, u, p){
  x * (1 - exp(-(lpnorm(x, p) / u - 1)))
}

#Partial derivative of weighting function
dWeightFun <- function(x, u, p){
  (1 - exp(-(lpnorm(x,p) / u - 1))) + x*gradlpnorm(x,p)/u * exp( - (lpnorm(x,p) / u - 1))
}

#Gradient of spectral log-likelihood for the scaled negative extremal Dirichlet model
gradient.negdir <- function(x, alpha, rho){
  logA <- lgamma(alpha-rho)-lgamma(alpha)
  (sum(alpha)-rho)*exp(-logA/rho-(1/rho+1)*log(x))/
    (rho*sum(exp(-(logA+log(x))/rho)))-(alpha/rho+1)/x
}

#Gradient of spectral log-likelihood for the scaled extremal Dirichlet model
gradient.dir <- function(x, alpha, rho){
  logC <- lgamma(alpha+rho)-lgamma(alpha)
  -(sum(alpha)+rho)*exp(logC/rho+(1/rho-1)*log(x))/
    (rho*sum(exp((logC+log(x))/rho)))+(alpha/rho-1)/x
}

#' Evaluate gradient score function
#'
#' @param dat \code{n} by \code{d} matrix of data
#' @param model string indicating whether the model is \code{dir} or \code{negdir}
#' @param weightFun object of class function, must be differentiable
#' @param dWeightFun object of class function, gradient of weightfun
#' @param alpha numeric vector of positive parameters
#' @param rho numeric parameter
#' @param u threshold
#' @param p norm for the risk function defined in weightFun (default is \eqn{l_p} norm)
#' @export
#' @return estimated score value
scoreEstimation <- function(dat, model, weightFun, dWeightFun, alpha, rho, u, p=1){
  if(class(weightFun) != "function") {
    stop('weightFun must be a function.')
  }
  if(class(dWeightFun) != "function") {
    stop('dweightFun must be a function.')
  }
  logC <- lgamma(alpha+rho)-lgamma(alpha)
  gradient.dir <- function(x, alpha, rho){
    -(sum(alpha)+rho)*exp(logC/rho+(1/rho-1)*log(x))/
      (rho*sum(exp((logC+log(x))/rho)))+(alpha/rho-1)/x
  }
  laplacian.dir <- function(x, alpha, rho){-(sum(alpha)+rho)*exp(logC/rho+(1/rho-1)*log(x))/
    (rho*sum(exp((logC+log(x))/rho)))*
    ((1/rho-1)/x-exp(logC/rho+(1/rho-1)*log(x))/
       (rho*sum(exp((logC+log(x))/rho))))-
  (alpha/rho-1)*exp(-2*log(x))
  }
  logA <- lgamma(alpha-rho)-lgamma(alpha)
  gradient.negdir <- function(x, alpha, rho){
    (sum(alpha)-rho)*exp(-logA/rho-(1/rho+1)*log(x))/
      (rho*sum(exp(-(logA+log(x))/rho)))-(alpha/rho+1)/x
  }
  laplacian.negdir <- function(x, alpha, rho){(sum(alpha)-rho)*exp(-logA/rho-(1/rho+1)*log(x))/
    (rho*sum(exp(-(logA+log(x))/rho)))*
    (-(1/rho+1)/x+exp(-logA/rho-(1/rho+1)*log(x))/
       (rho*sum(exp(-(logA+log(x))/rho))))+(alpha/rho+1)*exp(-2*log(x))
  }
    weights <- t(apply(dat, 1, weightFun, p=p, u=u))
    dWeights <- t(apply(dat, 1, dWeightFun, p=p, u=u))
    grad <- t(apply(dat, 1, switch(model, dir="gradient.dir", negdir="gradient.negdir"), alpha=alpha, rho=rho))
    lapl <- t(apply(dat, 1, switch(model, dir="laplacian.dir", negdir="laplacian.negdir"), alpha=alpha, rho=rho))
    sum(2 * (weights * dWeights) * grad + weights^2 * lapl + 0.5 * weights^2 * grad^2)
}

#' Fit Dirichlet extremal model using the gradient score
#'
#' Optimization routine to fit a member of the extremal Dirichlet family
#' using the gradient score
#' @export
#' @param start starting values for the parameter (\eqn{alpha} and \eqn{rho}, in this order)
#' @param dat matrix of transformed data on unit Pareto scale lying above the threshold
#' @param model string indicating whether the \code{"dir"} or \code{"negdir"} model should be considered
#' @param p int specifying the degree of the \eqn{l_p} norm used in the weight function \code{weightFun}
#' @param qu probability level for the marginal quantile for the sum
#' @param exceed logical indicating whether the data provided is the exceedances. Default to \code{FALSE}.
#' @return the optimized parameter values
#' @examples
#' samp <- mev::rmev(n=2500,d=3,param=c(1,2,3,0.5),model="negdir")
#' fscore(start=c(1,2,3,0.5), dat=samp, model="negdir")
fscore <- function(start, dat, model, p=1, qu=0.9, exceed=FALSE){
  optim.score <- function(par, dat, model, u, p=1){
      alpha <- exp(par[-length(par)])
      rho <- exp(par[length(par)])
      if(any(alpha>100) || rho < 0 || rho > 5){return(1e30)}
      #Invalid parameter value
      if(min(par)!=par[length(par)]){
        return(1e30)
      }
      scoreEstimation(dat=dat, model = model, weightFun = weightFun, dWeightFun = dWeightFun, p=p, u=u, alpha=alpha, rho=rho)
    }
  #Routine
  if(exceed){
    opt.out <- optim(par=log(start),fn=optim.score, dat=dat, model=model, p=p, u=u, control=list(fnscale=1,maxit=1500))
  } else{
  sums <- apply(dat, 1, lpnorm, p=p)
  u <- quantile(sums, qu)
  exceedances <- dat[sums > u,]
  opt.out <- optim(par=log(start),fn=optim.score, dat=exceedances, model=model, p=p, u=u, control=list(fnscale=1,maxit=1500))
  }
  return(list(par=exp(opt.out$par), convergence=opt.out$convergence, val=opt.out$val))
}

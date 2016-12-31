#' Wrapper for the scaled Dirichlet extremal model
#'
#' This function takes as argument the parameters of the scaled Dirichlet bivariate distribution
#' and passes them to the C++ function.
#'
#' @useDynLib demoroutines
#' @importFrom Rcpp evalCpp
#' @param optpar vector of parameters over which to optimize
#' @param fixedpar vector of parameters to keep constant
#' @param wfixed is vector of int with indices of data points
#' @param numpar vector integer giving the length of each parameter (scale, shape, ...)
#' @param transform boolean indicating whether to exp (expit for rho) the parameters for optimization
#' @param dat data matrix
#' @param thid logical matrix indicating whether entries exceed the marginal threshold
#' @param N total sample size
#' @param lambda vector of marginal rates of exceedance
#' @param u vector of thresholds
#' @return the negative log-likelihood value
#' @keywords internal
nlloptdir <- function(optpar, fixedpar, wfixed, numpar, transform=FALSE, dat, thid, N, lambda, u) {
	if((length(optpar)+length(fixedpar))!=sum(numpar)){
		stop("Invalid arguments passed to function `nllmvctopt'; please check")
	}
	#Transform the parameters if they are to be optimized (otherwise, don't)
	if(1 %in% wfixed){  mscale = fixedpar[1:numpar[1]]; fixedpar <- fixedpar[-(1:numpar[1])]
	} else { mscale = optpar[1:numpar[1]]; optpar <- optpar[-(1:numpar[1])]
	if(transform){mscale <- exp(mscale)}
	}
	if(2 %in% wfixed){  mshape = fixedpar[1:numpar[2]]; fixedpar <- fixedpar[-(1:numpar[2])]
	} else { mshape = optpar[1:numpar[2]]; optpar <- optpar[-(1:numpar[2])] }
	if(3 %in% wfixed){  malpha = fixedpar[1:numpar[3]]; fixedpar <- fixedpar[-(1:numpar[3])]
	} else { malpha = optpar[1:numpar[3]]; optpar <- optpar[-(1:numpar[3])]
	if(transform){malpha <- exp(malpha)}
	}
	if(4 %in% wfixed){  mrho = fixedpar[1:numpar[4]]; fixedpar <- fixedpar[-(1:numpar[4])]
	} else { mrho = optpar[1:numpar[4]]; optpar <- optpar[-(1:numpar[4])]
    	if(transform){
    	  mrho <- exp(mrho)#/(1+exp(mrho)) #@TODO check whether rho must lie in (0,1)
    	}
	}
	#Case of same shape or scale is handled in Marginal transform
  nllmvsdir(y=dat, thid=thid, N=N, lambda=lambda, u=u, alpha=malpha, rho=mrho, scale=mscale, shape=mshape)
}

#' Wrapper for the scaled negative extremal Dirichlet model
#'
#' This function takes as argument the parameters of the scaled Dirichlet bivariate distribution
#' and passes them to the C++ function.
#'
#' @param optpar vector of parameters over which to optimize
#' @param fixedpar vector of parameters to keep constant
#' @param wfixed is vector of int with indices of data points
#' @param numpar vector integer giving the length of each parameter (scale, shape, ...)
#' @param transform boolean indicating whether to exponentiate the parameters for optimization
#' @param dat data matrix
#' @param thid logical matrix indicating whether entries exceed the marginal threshold
#' @param N total sample size
#' @param lambda vector of marginal rates of exceedance
#' @param u vector of thresholds
#' @return the negative log-likelihood value
#' @keywords internal
nlloptnegdir <- function(optpar, fixedpar, wfixed, numpar, transform=FALSE, dat, thid, N, lambda, u) {
  if((length(optpar)+length(fixedpar))!=sum(numpar)){
    stop("Invalid arguments passed to function `nllmvctopt'; please check")
  }
  #Transform the parameters if they are to be optimized (otherwise, don't)
  if(1 %in% wfixed){  mscale = fixedpar[1:numpar[1]]; fixedpar <- fixedpar[-(1:numpar[1])]
  } else { mscale = optpar[1:numpar[1]]; optpar <- optpar[-(1:numpar[1])]
  if(transform){mscale <- exp(mscale)}
  }
  if(2 %in% wfixed){  mshape = fixedpar[1:numpar[2]]; fixedpar <- fixedpar[-(1:numpar[2])]
  } else { mshape = optpar[1:numpar[2]]; optpar <- optpar[-(1:numpar[2])] }
  if(3 %in% wfixed){  malpha = fixedpar[1:numpar[3]]; fixedpar <- fixedpar[-(1:numpar[3])]
  } else { malpha = optpar[1:numpar[3]]; optpar <- optpar[-(1:numpar[3])]
    if(transform){ malpha <- exp(malpha)}
  }
  if(4 %in% wfixed){  mrho = fixedpar[1:numpar[4]]; fixedpar <- fixedpar[-(1:numpar[4])]
  } else { mrho = optpar[1:numpar[4]]; optpar <- optpar[-(1:numpar[4])]
    if(transform){ mrho <- exp(mrho)
    }
  }
  #Case of same shape or scale is handled in Marginal transform
  nllmvsnegdir(y=dat, thid=thid, N=N, lambda=lambda, u=u, alpha=malpha, rho=mrho, scale=mscale, shape=mshape)
}


#' Test for conditional negative semi-definitiveness
#' Is the matrix of rank n-1 satisfying sum-to-zero constraints conditionally negative definite?
#' @param X matrix to test for CNSD
#' @param tol a tolerance factor for the lowest eigenvalue
#' @return a logical indicating whether the matrix \code{X} is CNSD.
.is.CNSD <- function (X, tol = 1e-08){
	n <- nrow(X)
	P <- diag(n)
	if(n>2){
		diag(P[,-1]) <- -1
	} else if(n==2){#error with one dimensional case...
		P[1,-1]	<- -1
	}
	Xhat <- P %*% X %*% t(P)
	eigs <- eigen(Xhat[-n, -n], symmetric=TRUE, only.values=TRUE)$values
	!eigs[1] > tol
}

#' Create a conditionally negative definite matrix from Cholesky decomposition
#'
#' This function takes a vector, fills a Cholesky upper triangular matrix and transforms it to a
#' positive definite matrix, then into a symmetric CNSD matrix with zero diagonal for the HR model
#'
#' @param vec vector of parameters for the Cholesky decomposition, with log-diagonal entries
#' @param dime dimension of the model
#' @return a CNSD matrix
#' @keywords internal
chol2Lambda <- function(vec, dime, EoW=TRUE){#vector of parameter and dimension of the model
	if(!is.integer(dime <- as.integer(dime))){
		stop("Invalid dimension for the Husler-Reiss model; argument must be an integer")
	}
	L = matrix(0, nrow=dime, ncol=dime)
	#diagonal elements are on the log-scale to ensure uniqueness of Cholesky decompo
	diag(L) <- exp(vec[1:dime])
	#other elements stored per row
	if(dime>1){
		L[upper.tri(L, diag=FALSE)] <- vec[-(1:dime)]
	}
	#create Sigma
	Sigma <- t(L)%*%L
	#create Lambda
	Sigma2Lambda(Sigma)
}

#' Transform covariance matrix to conditionally negative definite matrix
#'
#' @param Sigma symmetric positive definite matrix
#' @return a CND matrix
#' @keywords internal
Sigma2Lambda <- function(Sigma){
				rbind(cbind(0, t(diag(Sigma))),cbind(diag(Sigma),
							rep(1,nrow(Sigma))%*%t(diag(Sigma))+diag(Sigma)%*%t(rep(1, nrow(Sigma)))-2*Sigma))/4
}

#' Transform a symmetric conditionally negative definite matrix into a covariance matrix
#' @param Lambda symmetric conditionally negative definite matrix of parameters for the Husler-Reiss model
#' @return a symmetric positive definite matrix
#' @keywords internal
Lambda2Sigma <- function(Lambda){
	2*(matrix(Lambda[,1],nrow(Lambda),nrow(Lambda))+t(matrix(Lambda[1,],nrow(Lambda),nrow(Lambda)))-Lambda)[-1,-1]
}

#' Wrapper for the Husler-Reiss model
#'
#' This function takes as argument the parameters of the Husler-Reiss bivariate distribution
#' and passes them to the C++ function.
#'
#' @inheritParams nlloptdir
#' @return the negative log-likelihood value
#' @keywords internal
nllopthr <- function(optpar, fixedpar, wfixed, numpar, transform=FALSE, dat, thid, N, lambda, u) {
	if((length(optpar)+length(fixedpar))!=sum(numpar)){
		stop("Invalid arguments passed to function `nllmvhropt'; please check")
	}
	#Transform the parameters if they are to be optimized (otherwise, don't)
	if(1 %in% wfixed){  mscale <- fixedpar[1:numpar[1]]; fixedpar <- fixedpar[-(1:numpar[1])]
	} else { mscale <- optpar[1:numpar[1]]; optpar <- optpar[-(1:numpar[1])]
	if(transform){mscale <- exp(mscale)}
	}
	if(2 %in% wfixed){  mshape <- fixedpar[1:numpar[2]]; fixedpar <- fixedpar[-(1:numpar[2])]
	} else { mshape <- optpar[1:numpar[2]]; optpar <- optpar[-(1:numpar[2])] }
	if(3 %in% wfixed){  mLambda = fixedpar;
	#If fixed, there is one check once for the parameters and the latter simply need to be filled
	lambdamat <- matrix(mLambda, ncol(dat),ncol(dat))
	} else { mLambda <- optpar #just whatever is left over for this model
	if(transform){lambdamat <- chol2Lambda(mLambda, ncol(dat)-1)
	} else{#TODO determine whether to fill by row or column #DOES THIS WORK?
		lambdamat <- 	matrix(0, ncol=ncol(dat),nrow=ncol(dat))
		lambdamat[upper.tri(lambdamat)] <- mLambda
		lambdamat <- lambdamat + t(lambdamat)
	}
	}
	#Case of same shape or scale is handled in marginal transform
	nllmvhr(y=dat, thid=thid, N=N, lambda=lambda, u=u, Lambda=lambdamat, scale=mscale, shape=mshape)
}

#' Pairwise composite likelihood estimation for peaks-over-threshold
#'
#' This function finds the maximum composite likelihood estimate
#' for multivariate extreme value models based on pairs of data.
#' It returns a list with the result from \code{\link[stats]{optim}}
#' along with some of the information passed to the function.
#' @param dat \code{n} by \code{d} data matrix on original scale
#' @param u \code{d} vector of threshold parameters on probability scale
#' @param start vector of starting value for the optimization
#' @param N total number of observations if some observations falling under all threshold were removed from \code{x} (can be omitted)
#' @param lambda vector of percentage of threshold exceedances (can be omitted)
#' @param model family of multivariate extreme value distribution. \code{'ct'} and \code{'hr'} are the only family currently supported.
#' @param cscale logical indicating common scale for all variables
#' @param cshape logical indicating common shape for all variables
#' @param sym logical indicating whether equal dependence parameters in the Dirichlet model
#' @param std.err logical indicating whether to calculate standard errors numerically
#' @param method optimization method; see \link[stats]{optim}
#' @param warn.inf whether to check the validity of starting values before optimization
#' @param transform logical indicating whether to perform optimization on unconstrained
#'  space; \code{start} should be adjusted if \code{transform=TRUE}
#' @param ... additional arguments (fixed parameters) passed to the \code{optim}
#' function to be held fixed during the optimization
#' @return a list containing the output of \code{\link[stats]{optim}}, notably
#' parameter values \code{par} and standard errors \code{se}
#' @export
#' @importFrom stats optim var
#' @importFrom numDeriv hessian
#' @importFrom mev gp.fit
#'
#' @examples
#' x <- mev::rmev(n=100000, d=2, param=c(1,2,0.5), model = "negdir")
#' qu <- apply(x, 2, quantile,probs <- 0.99)
#' y <- x[which(rowSums(isAbove(x,threshold = qu))>=1),]
#' fit <- fmvcpot(dat=x, u=qu,model="negdir",lambda=colSums(isAbove(y,qu))/(nrow(x)+1), N=nrow(x),
#'                cscale=TRUE, cshape=TRUE, shape=1, scale=100)
#' fit$par
fmvcpot <- function(dat, u, lambda, N, model=c("ct","dir","negdir","hr"), cscale=FALSE, cshape=FALSE,
                      sym=FALSE, start=NULL, std.err=TRUE, method="BFGS", warn.inf=TRUE,...){
  #@TODO handle missing values in some of the vector (not sure it works) 06.10.2016
  dat <- na.omit(dat)
	#keep for the marginal analysis the cases - so contributions must be handled in the C++ function
  if(length(u)!=ncol(dat)){
    stop("Invalid threshold; must match the dimension of dat");
  }
  if(missing(lambda)){ #sep.bvdata does this automatically
    lambda <- sapply(1:ncol(dat), function(ind){sum(dat[,ind]>u[ind])/(nrow(dat)+1)})
  }
  if(model=="ct"){ model="dir"}
  if(missing(N)){ N <- nrow(dat)}
  thid <- isAbove(dat, u)
  nulls <- which(rowSums(thid)==0)
  #Remove values that are unnecessary (because they are censored)
  #This speeds up the calculations quite drastically, because optim must loop most of the zeros ow.
  if(length(nulls)>0){
    dat <- dat[-nulls,]
    thid <- thid[-nulls,]
  }

  #COMMON ROUTINE

  #Define name of parameters to be optimized - will automatically match `nm',
  #but recall `start' may be provided by user
  param <- c("scale", "shape")
  if(model %in% c("dir","negdir")){ #Currently the only case, but could be changed latter.
    param <- c(param, "alpha","rho")
    numpar <- c(ifelse(cscale,1,ncol(dat)),ifelse(cshape,1,ncol(dat)),ifelse(sym,1,ncol(dat)), 1)
  } else if(model=="hr"){
   param <- c(param, "Lambda")
   numpar <- c(ifelse(cscale,1,ncol(dat)),ifelse(cshape,1,ncol(dat)),ncol(dat)*(ncol(dat)-1)/2)
  }
  nmdots <- names(list(...))
  fixed.param <- list(...)[nmdots %in% param]
  #Define starting values (if not provided by the user)
  if (is.null(start)) {#default
    start <- composite_pot_startvals(x=dat, u=u, model=model, sym=sym, cscale=cscale,
                                     cshape=cshape, transform=TRUE) #TODO check incidence
  } else if (!is.list(start)){
    stop("`start' must be a named list")
  }
  nm <- names(start)
  #For now, just matching the parameters; will deal with length later
  if(isTRUE(any(!(param %in% c(nm,names(fixed.param)))))){
    stop("Unspecified parameters")
  }
  m <- match(nm, param)
  if(any(is.na(m))){
    stop("`start' specifies unknown arguments")
  }
  fixed <- c()
  wfixed <- c()
  if("scale" %in% nmdots){
    fixed <- c(fixed,as.vector(list(...)['scale']))
    start$scale <- NULL
    if(length(unlist(list(...)['scale']))!=numpar[1]){
      stop("Invalid scale parameter: it does not match the constraint provided.")
    }
    wfixed <- c(wfixed, 1)
  }
  if("shape" %in% nmdots){
    fixed <- c(fixed,as.vector(list(...)['shape']))
    start$shape <- NULL
    if(length(unlist(list(...)['shape']))!=numpar[2]){
      stop("Invalid shape parameter: it does not match the constraint provided.")
    }
    wfixed <- c(wfixed, 2)
  }
  #Model specific arguments
  if(model %in% c("dir","negdir")){
    if("alpha" %in% nmdots){
      fixed <- c(fixed,as.vector(list(...)['alpha']))
      start$alpha <- NULL
      if(length(unlist(list(...)['alpha']))!=numpar[3]){
        stop("Invalid `alpha' parameter provided")
      }
      wfixed <- c(wfixed, 3)
    }
    if("rho" %in% nmdots){
      fixed <- c(fixed,as.vector(list(...)['rho']))
      start$rho <- NULL
      if(length(unlist(list(...)['rho']))!=numpar[4]){
        stop("Invalid `rho' parameter provided")
      }
      wfixed <- c(wfixed, 4)
    }
  } else if(model=="hr"){
  	if("Lambda" %in% nmdots){
  		start$Lambda <- NULL
  		numpar[3] <- ncol(dat)^2
  		#User can provide either parameters or the full matrix...
  		lambda_par <- list(...)['Lambda'][[1]]
  		if(is.matrix(lambda_par)){
  			if(all(dim(lambda_par)==ncol(dat))){
  				if(.is.CNSD(lambda_par)){
  					fixed <- c(fixed,as.vector(lambda_par))
  				}
  			}
  		} else if(is.vector(lambda_par)){
  			if(all.equal(length(lambda_par),ncol(dat)*(ncol(dat)-1)/2)){
  				Lambda_mat <- matrix(0,ncol(dat),ncol(dat))
  				Lambda_mat[upper.tri(Lambda_mat,diag=FALSE)] <- lambda_par
  				Lambda_mat <- Lambda_mat + t(Lambda_mat)
  				if(.is.CNSD(Lambda_mat)){
  					fixed <- c(fixed,as.vector(Lambda_mat))
  				}
  			}
  		} else {
  			stop("Invalid `Lambda' parameter provided")
  		}
  		wfixed <- c(wfixed, 3)
  	}
  }
  start.arg <- unlist(start)
  fixed <- unlist(fixed)
  if(is.null(start.arg)){
    #If all parameters are provided, simply evaluate the log-likelihood
  	nll <- switch(model,
  								dir=nlloptdir(start.arg, fixed, wfixed, numpar,transform=FALSE,
                     dat=dat, thid=thid, lambda=lambda, N=N, u=u),
  								negdir=nlloptnegdir(start.arg, fixed, wfixed, numpar,transform=FALSE,
  								              dat=dat, thid=thid, lambda=lambda, N=N, u=u),
  								hr=nllopthr(start.arg, fixed, wfixed, numpar,transform=FALSE,
  														 dat=dat, thid=thid, lambda=lambda, N=N, u=u)
  	)
  	return(nll)
  } else{ #Optimize
  	if(warn.inf && isTRUE(all.equal(do.call(switch(model, dir="nlloptdir", negdir="nlloptnegdir",hr="nllopthr"),#TODO add other models
     list(optpar=start.arg, fixedpar=fixed, wfixed=wfixed, numpar=numpar, transform=TRUE,
                                dat=dat, thid=thid, lambda=lambda, N=N, u=u)),1e10))){
      warning("Negative log-likelihood is infinite at starting values")
    }
    #If start.arg was provided by user, should be inv.transformed
    opt <- optim(par=start.arg, fn=get(switch(model, dir="nlloptdir", negdir="nlloptnegdir", hr="nllopthr")),
    						 hessian = std.err, fixedpar=fixed, wfixed=wfixed,
                 numpar=numpar, transform=TRUE, method = method, control=list(maxit=2500),
                 dat=dat, thid=thid, lambda=lambda, N=N, u=u)
    #print(opt)
    if(opt$convergence!=0){
      warning("Optimization routine did not converge");
    }
	    start <- NULL
	    tpar <- opt$par
	    if(isTRUE(all(opt$convergence==0, !1 %in% wfixed))){
	      if(! (1 %in% wfixed)){
	        start <- c(start,exp(tpar[1:numpar[1]])); tpar <- tpar[-(1:numpar[1])]
	      }
	      if(! (2 %in% wfixed)){
	        start <- c(start,tpar[1:numpar[2]]); tpar <- tpar[-(1:numpar[2])]
	      }
	    }
	    if(model %in% c("dir","negdir")){
	    	if(isTRUE(all(opt$convergence==0, any(! c(3,4)%in% wfixed)))){
		      if(! (3 %in% wfixed)){
		        start <- c(start,exp(tpar[1:numpar[3]])); tpar <- tpar[-(1:numpar[3])]
		      }
		      if(! (4 %in% wfixed)){
		        start <- c(start,exp(tpar[1:numpar[4]]))#/(1+exp(tpar[1:numpar[4]])));
		        tpar <- tpar[-(1:numpar[4])]
		      }
	    	}
	    	opt$par <- start;
	      opt$hessian <- numDeriv::hessian(func=switch(model, dir=nlloptdir, negdir=nlloptnegdir), x=opt$par, fixedpar=fixed, wfixed=wfixed,
  	                                       numpar=numpar, transform=FALSE, dat=dat, thid=thid,
  	                                       lambda=lambda, N=N, u=u)
    	}
				opt$estimate <- list()
				#Put scale and shape in the estimates
				if(1 %in% wfixed){
					opt$estimate[['scale']] <- rep(as.vector(unlist(list(...)['scale'])),length.out = ncol(dat))
				} else{
					opt$estimate[['scale']] <- rep(as.vector(opt$par[grep('scale',names(opt$par))]),length.out = ncol(dat))
				}
					if(2 %in% wfixed){
					opt$estimate[['shape']] <- rep(as.vector(unlist(list(...)['shape'])),length.out = ncol(dat))
				} else{
					opt$estimate[['shape']] <- rep(as.vector(opt$par[grep('shape',names(opt$par))]),length.out = ncol(dat))
				}
				if(model %in% c("dir","negdir")){
					if(3 %in% wfixed){
	      	opt$estimate[['alpha']] <- rep(as.vector(unlist(list(...)['alpha'])),length.out = ncol(dat))
				} else{
					opt$estimate[['alpha']] <- as.vector(opt$par[grep('alpha',names(opt$par))])
				}
	      if(4 %in% wfixed){
	      	opt$estimate['rho'] <- as.vector(unlist(list(...)['rho']))
	      } else{
	      	opt$estimate['rho'] <- opt$par[grep('rho',names(opt$par))]
	     }
	    } else if(model=="hr"){
				#passing on the parameters will fail here because the function nllopthr recognizes fixed arguments
	    	Lambda_mle <- chol2Lambda(tpar, ncol(dat)-1)
	    	start <- c(start,Lambda = c(Lambda_mle[upper.tri(Lambda_mle,diag=FALSE)]))
	    	opt$par <- start;
	    	# opt$hessian <- numDeriv::hessian(func=nllopthr, x=start, fixedpar=fixed, wfixed=wfixed,
	    	# 																 numpar=numpar, transform=FALSE, dat=dat, thid=thid,
	    	# 																 lambda=lambda, N=N, u=u)
				if(3 %in% wfixed){
	      	opt$estimate[['Lambda']] <- list(...)['Lambda']
				} else{
	    		opt$estimate[['Lambda']] <- Lambda_mle
	    	}
	    }


	    #Add the variance matrix and the standard errors to the opt list
	    if(std.err){
	    	if(model %in% c("dir","negdir")){
	         	opt$var <- try(solve(opt$hessian))
		     if(is.character(opt$var)){
		      opt$var <- NULL
		     } else{
		      opt$se <- sqrt(diag(opt$var))
		     }
	    	}
	    }#TODO include the HR model

    #Saving intermediate elements
    opt$model <- model
    opt$u <- as.vector(u)
    opt$cscale <- cscale
    opt$cshape <- cshape
    if(model %in% c("dir","negdir")){opt$sym <- sym}
    opt$N <- N
    opt$lambda <- lambda
    return(opt)
    #print(fixed)
    # 	hin <- function(optpar, fixedpar, wfixed, numpar, transform=FALSE, dat, thid, N, lambda, u){
    # 	  c(optpar[hind],ifelse(4 %in% wfixed, NULL, c(optpar[length(optpar)],1-optpar[length(optpar)])))
    # 	}
    # 	#if parameters were transformed
    # 	#Catch the non-transformed parameters (see which are returned)
    # 	print(hind)
    # 	library(alabama)
    # 	opt <- constrOptim.nl(par=start, fn=nlloptdir, hin = hin,fixedpar=fixed, wfixed=wfixed,
    # 	                      numpar=numpar, transform=FALSE, dat=dat, thid=thid, lambda=lambda, N=N, u=u)
    # 	print(opt)
    # # opt <- optim(par=start, fn=nlloptdir, hessian = std.err, fixedpar=fixed, wfixed=wfixed,
    # #              numpar=numpar, transform=FALSE, method = method, control=list(maxit=2500),
    # # 						x=x, thid=thid, lambda=lambda, N=N, u=u)
    #
  }
  # return(opt)
}



#' Starting values for multivariate peaks-over-threshold via pairwise likelihood methods
#'
#' This function is the analog of \link[evd]{bvstart.vals}, except that it handles higher dimensional models
#'
#' @param x \code{n} by \code{d} data matrix on original scale
#' @param u \code{d} vector of threshold parameters on probability scale
#' @param model string indicating the multivariate extreme value family
#' @param sym logical indicating whether model should be made symmetric
#' @param cscale logical indicating common scale for all variables
#' @param cshape logical indicating common shape for all variables
#' @param transform logical indicating transformation to unconstrained domain for parameters
#' @keywords internal
composite_pot_startvals <- function(x, u=NULL, model, sym=FALSE,
                                    cscale=FALSE, cshape=FALSE, transform=FALSE){
  if(length(u)!=ncol(x)){
    stop("Threshold length must match the provided data")
  }
  #Define containers
  start <- list()
  scale <- rep(1, ncol(x))
  shape <- rep(0.1, ncol(x))
  #Marginal parameters
  if(!cscale || !cshape){
    for(i in 1:ncol(as.matrix(x))){
      out <- mev::gp.fit(xdat=x[,i], threshold=u[i], method="Grimshaw", show=FALSE)
      if(!cscale && !cshape){
        scale[i] <- out$estimate[1][[1]]
        shape[i] <- out$estimate[2][[1]]
        #start[[paste0("scale",i)]] <- out$estimate[1][[1]]
        #start[[paste0("shape",i)]] <- out$estimate[2][[1]]
      } else if(cscale){
        shape[i] <- out$estimate[2][[1]]
      } else if(cshape){
        scale[i] <- out$estimate[1][[1]]
      }
    }
  }
  if(cscale || cshape){
    #If only one parameter is fixed, take average
    out <- mev::gp.fit(xdat=c(t(x)-u), threshold=0, method="Grimshaw",show=FALSE)
    if(cscale){ scale <- out$estimate[1][[1]]	}
    if(cshape){ shape <- out$estimate[2][[1]]	}
  }
  if(transform){
    start[["scale"]] <- log(scale)
  } else{
    start[["scale"]] <- scale
  }
  start[["shape"]] <- shape
  if(model == "log"){
    start["dep"] = 0.75
  } else if(model %in% c("dir","negdir")){
   if(!sym){
      if(transform){
        start[["alpha"]] <- rep(0.1, ncol(x))
      } else{
        start[["alpha"]] <- rep(1.1, ncol(x))
      }
    } else{ #if sym
      if(transform){
        start[["alpha"]] <- 0.1
      } else{
        start[["alpha"]] <- 1.1
      }
    }
      rst = 0.7
      if(transform){
        start["rho"] = log(rst)#-log(1-rst) #logit(irv)
      } else{
        start["rho"] = rst
      }
   } else if(model == "hr"){
   	 sigMatstart <- matrix(0.5,ncol(x)-1,ncol(x)-1)+diag(0.5,ncol(x)-1)
     L <- chol(sigMatstart)
     start[['Lambda']] <- c(log(diag(L)), L[upper.tri(L, diag=FALSE)])
		} else if(model == "xstud"){
		  start["df"] = 1
		  start[[Rho]] <- matrix(0.5,ncol(x)-1,ncol(x)-1)+diag(0.5,ncol(x)-1) #TODO : check dimension
   }
  return(start)
}


#' Estimation of the variability and Godambe information matrix using a
#' nonparametric bootstrap scheme
#'
#' Composite likelihood estimators are consistent and asymptotically normal
#' under mild regularity conditions. The asymptotic covariance matrix is the
#' inverse Godambe information matrix \eqn{G^{-1}}, which is estimated by a
#' nonparametric bootstrap as the empirical covariance of \code{B} replicates.
#' The sensitivity matrix \eqn{H} is calculated from the Hessian matrix at the
#' maximum composite likelihood estimates. Lastly, the variability matrix \eqn{J} is
#' obtained from the relation \eqn{G=HJ^{-1}H}
#'
#' @param dat data matrix
#' @param fitted output from the call to \code{fmvcpot}.
#' @param B number of bootstrap replicates
#' @param ... fixed pararameters to pass to \code{fmvcpot} if any
#' @param use.start logical indicating whether to use MCLE as starting value
#' @return a list with matrices \code{godambe}, \code{sensitivity} and \code{variability}
#' @export
compositemat <- function(dat, fitted, B, use.start=FALSE, ...){
  #Sensitivity matrix H
  Hmat <- fitted$hessian
  bootpar <- matrix(NA,nrow=B, ncol=length(fitted$par))
  for(b in 1:B){
    set.seed(b)
    st <- list()
	  if(fitted$model %in% c("dir","negdir")){
	    parnames <- c("scale","shape","alpha","rho")
	  }
	  #start must be a named list
	  whichpar <- sapply(substr(names(fitted$par),1,3), function(nm){pmatch(nm,parnames)})
	  for(p in parnames){
	    temp  = which(whichpar==switch(p, "scale"=1, "shape"=2,"alpha"=3,"rho"=4))
	    if(length(temp)>0){
	   st[[p]] <- fitted$par[temp]
	    }
	  }
	  if(use.start){
	  optimal <- fmvcpot(dat[sample(1:nrow(dat), size=nrow(dat), replace = TRUE),],
	            u=fitted$u,  model=fitted$model,start=st, std.err=FALSE, method="Nelder-Mead",
	            transform=TRUE,cscale=fitted$cscale, cshape=fitted$cshape, sym=fitted$sym,...)
	  } else{
	  	optimal <- fmvcpot(dat[sample(1:nrow(dat), size=nrow(dat), replace = TRUE),],
	            u=fitted$u,  model=fitted$model, std.err=FALSE, method="Nelder-Mead",
	            transform=TRUE,cscale=fitted$cscale, cshape=fitted$cshape, sym=fitted$sym,...)
	  }
	  if(optimal$convergence==0){
	    bootpar[b,] <- optimal$par
	  }
  }
  print(bootpar)
  covmat <- var(bootpar,na.rm = TRUE)
  godambe <- solve(covmat)
  Jmat <- Hmat%*%covmat%*%Hmat
  return(list(sensitivity=Hmat, godambe=godambe, variability=Jmat))
}

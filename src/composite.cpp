// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
//////////
//Functions by Leonid Rousniak and Adrien de Castro

// Wraps R::pnorm in C++
inline double p_norm(double x) {
  return R::pnorm(x,0,1,true,false);
}

// Wraps R::dnorm in C++
inline double d_norm(double x) {
  return R::dnorm(x,0,1,false);
}

// Wraps R::pnorm in C++
inline double p_t(double x, double nu) {
  return R::pt(x,nu,true,false);
}

// Wraps R::dnorm in C++
inline double d_t(double x, double nu) {
  return R::dt(x,nu,false);
}

///////////////////////////////////////////
// Functions for HR model /////////////////
///////////////////////////////////////////
inline double ztr(double z1, double z2, double a) {
  return a/2.0-log(z1/z2)/a;
}

inline double V_BR(double z1, double z2, double a) {
  return p_norm(ztr(z1,z2,a))/z1 + p_norm(ztr(z2,z1,a))/z2;
}

inline double V_BR1(double z1, double z2, double a) {
  return -(p_norm(ztr(z1,z2,a)) + d_norm(ztr(z1,z2,a))/a)/(z1*z1) +
    d_norm(ztr(z2,z1,a))/(a*z1*z2);
}

inline double V_BR2(double z1, double z2, double a) {
  return V_BR1(z2, z1, a);
}

inline double V_BR12(double z1, double z2, double a) {
  return -(d_norm(ztr(z1,z2,a))/z1*(1-ztr(z1,z2,a)/a) +
           d_norm(ztr(z2,z1,a))/z2*(1-ztr(z2,z1,a)/a)) / (a*z1*z2);
}


///////////////////////////////////////////
// Functions for Extremal-t model /////////
///////////////////////////////////////////
inline double exstudf(double z1, double z2, double b, double rho, double nu) {
  return b*(pow(z2/z1,1/nu) - rho);
}

inline double V_EXST(double z1, double z2, double b, double rho, double nu) {
  return p_t(exstudf(z1,z2,b,rho,nu),nu)/z1 + p_t(exstudf(z2,z1,b,rho,nu),nu)/z2;
}

inline double V_EXST1(double z1, double z2, double b, double rho, double nu) {
  return -p_t(exstudf(z1,z2,b,rho,nu),nu)/(z1*z1);
}

inline double V_EXST2(double z1, double z2, double b, double rho, double nu) {
  return V_EXST1(z2,z1,b,rho,nu);
}

inline double V_EXST12(double z1, double z2, double b, double rho, double nu) {
  return -(b/nu)*d_t(exstudf(z1,z2,b,rho,nu),nu)*pow(z2,1.0/nu-1.0)/pow(z1,1.0/nu+2.0);
}
//////////

//' Marginal transformation of the observations to Frechet or Pareto scale
//'
//' This function uses the probability integral transform to map the observations
//' to unit Frechet scale \code{GEV(1,1,1)} or unit Pareto distribution \code{GP(1,1,1)}.
//' The function returns a list with the transformed observation and the marginal transformation
//' Jacobian.
//'
//' @param data \code{n} by \code{d} matrix of observations on original scale
//' @inheritParams nllmvsdir
//' @param GEVorGP logical indicating whether to transform to unit Frechet or unit Pareto. Default to \code{true} (GEV).
//' @export
// [[Rcpp::export]]
List MarginalTransformation(NumericMatrix data, NumericVector u, NumericVector lambda,
                            NumericVector scale, NumericVector shape, bool GEVorGP = true){
  //Define containers
  NumericMatrix t(data.nrow(), data.ncol());
  NumericMatrix dat(data.nrow(), data.ncol());
  bool outOfDomain = false;
  for(int j=0; j < data.ncol(); j++)  {
    //Scale the data
    dat(_,j) = pmax(data(_,j)-u[j],0.0) / scale[j];
    //Gumbel case if within numerical tolerance
    if(shape[j] < 1e-8 && shape[j] > -1e-8){//if(shape[j] == 0){
      t(_,j) = exp(-dat(_,j));
    } else {
      t(_,j) = 1.0 + shape[j] * dat(_,j);
      //Check first marginal domain restriction
      if(is_true(any(t(_,j) <= 0))) { // TODO O(nd) b/c linear pass in data; is it necessary to keep this here?
        outOfDomain = true;
        return(List::create(_["outOfDomain"]  = outOfDomain));
      }
    }
    t(_,j) = exp(-log(t(_,j)) / shape[j]);
    //Transform to unit Frechet from uniform
    if(GEVorGP){
    dat(_,j) = -1.0/log(1.0 - lambda[j] * t(_,j));
    t(_,j) = lambda[j]/scale[j]*(pow(dat(_,j),2.0)*
      pow(t(_,j),1.0+shape[j])/(1.0-lambda[j]*t(_,j)));
    } else{
      dat(_,j) = 1.0/(lambda[j] * t(_,j));
      t(_,j) = lambda[j]/scale[j]*pow(dat(_,j),2.0)*
        pow(t(_,j),1.0+shape[j]);
    }
    }
  return(List::create(_["outOfDomain"]  = outOfDomain, _["t"] = t, _["dat"] = dat));
}



//' Indicates whether observations lie above the threshold.
//'
//' The function returns a boolean matrix of indicators.
//'
//' @param data matrix of all observations
//' @param threshold vector of thresholds under which observations are censored
//' @return matrix of logical indicating whether or not observations lie above the marginal threshold
//' @keywords internal
//' @export
// [[Rcpp::export]]
LogicalMatrix isAbove(NumericMatrix data, NumericVector threshold){
  if(data.ncol() != threshold.size()){
    Rcpp::stop("Invalid threshold vector");
  }
  LogicalMatrix thid(data.nrow(), data.ncol());
  for(int i=0; i < data.nrow(); i++){
    for(int j=0; j < data.ncol(); j++){
      thid(i,j) = data(i,j)> threshold[j];
    }
  }
  return(thid);
}


//' Negative log-likelihood of scaled Dirichlet extreme value model
//'
//' This is the workhorse with the likelihood, which transforms observations
//' and returns the objective function value.
//'
//' @param y \code{n} by \code{d} matrix of observations on original scale
//' @param u threshold (on \code{y}-scale), to be substracted from y
//' @param thid \code{n} by \code{d} matrix of logical vector indicating whether or not a
//' threshold exceedance was observed
//' @param N total number of observations (not just exceedances)
//' @param scale vector of scale parameters
//' @param shape vector of shape parameters
//' @param lambda vector of percentage of threshold exceedances
//' @param alpha vector of parameters for the Dirichlet model
//' @param rho regular variation parameter
//' @return the negative log-composite likelihood value
//' @section Note: The \code{N} argument is there mostly for cases in which function returns
//' observations stripped from their contribution; in such cases, the count
//' would otherwise be `incorrect'.
//' @export
// [[Rcpp::export]]
NumericVector nllmvsdir(NumericMatrix y, LogicalMatrix thid, int N,
                        NumericVector lambda, NumericVector u, NumericVector alpha,
                        double rho, NumericVector scale, NumericVector shape) {
  if(scale.size()==1){	scale = rep(scale, y.ncol());}
  if(shape.size()==1){	shape = rep(shape, y.ncol());}
  if(alpha.size()==1){	alpha = rep(alpha, y.ncol());} //Symmetric
  if(shape.size()!=y.ncol() || scale.size()!=y.ncol()){
    Rcpp::stop("Invalid marginal parameter vector size for the GPD");
  }
  NumericVector dvec(1);
  NumericVector ll(1);
  //Marginal transformation
  List MT = MarginalTransformation(y, u, lambda, scale, shape);
  if(MT["outOfDomain"]){
    ll[0] = 1e10; //This error code, `Inf' is catched in the optim wrapper
    return ll;
  }
  if(rho < 0 || is_true(any(alpha>100)) || is_true(any(alpha<0))){ //rho > 1 ||
    ll[0] = 1e10; //This `Inf' is catched in the optim wrapper
    return ll;
  }
  NumericVector logC = NumericVector(y.ncol());
  for(int i=0; i<y.ncol(); i++){
    logC[i] = lgamma(alpha[i]+rho)-lgamma(alpha[i]);
  }
  NumericMatrix t = MT["t"];
  NumericMatrix data = MT["dat"];
  NumericVector lambda2 = -1.0/log(1.0-lambda);  //value of threshold on Frechet scale
  //Containers
  NumericVector x(1);   // transformed variable Upsilon (only one point because angle on [0,1])
  NumericVector zdn(1); // contribution for points below both thresholds
  NumericVector v(1);   // exponent measure and derivatives
  NumericVector v1(1);  // these are set to zero by default
  NumericVector v2(1);  // and we overwrite them to save space
  NumericVector v12(1); // joint exceedances
  for(int l=0; l < y.ncol(); l++) { //Loop over variable indices
    checkUserInterrupt();
    for(int k=0; k < l; k++) { // Loop over other pairs (use symmetry)
      // Transformation from Dirichlet mixture space to density
      //Recall k (aka "1") < l (aka "2") - zdn is correct below
      zdn[0] = exp((logC[k] + log(lambda2[k])) / rho) /
        (exp((logC[k] + log(lambda2[k])) / rho) +
          exp(( logC[l] + log(lambda2[l])) / rho));
      //Contribution of values under both threshold
      // V - exponent measure  -- censored contributions, depends only on thresholds
      zdn[0] = -R::pbeta(zdn[0], alpha[k] + rho, alpha[l], 0, 0) / lambda2[k] -
        R::pbeta(zdn[0], alpha[k], alpha[l] + rho, 1, 0) / lambda2[l];
      for(int i=0; i < y.nrow(); i++)  { // Loop over the observations
        if(!R_IsNA(data(i,k)) && !R_IsNA(data(i,l))){
        if(thid(i,k)==true || thid(i,l)==true){
        //If either is marginally exceeding, else skip calculations
        //data is data on unit Frechet scale
        //t contains the marginal censored contributions
        //x are transformed observations
        x[0] = exp((logC[k] + log(data(i,k))) / rho) /
          (exp((logC[k] + log(data(i,k))) / rho) +
            exp((logC[l] + log(data(i,l))) / rho));
        v[0] = exp(R::pbeta(x[0], alpha[k] + rho, alpha[l], 0, 1) - log(data(i,k))) +
          exp(R::pbeta(x[0], alpha[k], alpha[l] + rho, 1, 1)-log(data(i,l)));
        v1[0] = -R::pbeta(x[0], alpha[k] + rho, alpha[l], 0, 0)/ pow(data(i,k),2);
        v2[0] = -R::pbeta(x[0], alpha[k], alpha[l] + rho, 1, 0)/ pow(data(i,l),2);
        v12[0] = - exp((alpha[k]/rho-1)*log(data(i,k))+(alpha[l]/rho-1)*log(data(i,l))+
          (alpha[k]/rho)*logC[k]+(alpha[l]/rho)*logC[l]-log(rho)+lgamma(alpha[k]+alpha[l]+rho)-
          lgamma(alpha[k])-lgamma(alpha[l])-(alpha[k]+alpha[l]+rho)*
          log(exp((logC[k]+log(data(i,k)))/rho)+exp((logC[l]+log(data(i,l)))/rho)));
          // pow(x[0],alpha[k]) * pow(1-x[0], alpha[l])/(data(i,k) * data(i,l) *
          // rho * pow((pow( exp(lgamma(alpha[k] + rho)+lgamma(alpha[l]) -
          // lgamma(alpha[k]+rho+alpha[l])) * data(i,k), 1 / rho) +
          // pow(exp(lgamma(alpha[l] + rho) + lgamma(alpha[k])-lgamma(alpha[k]+rho+alpha[l])) * data(i,l), 1 / rho)), rho))
        }
        //Function from evd: 1 (respectively 2) for marginal exceedance of x_1 (x_2), 3 for joint exceedance

        if(thid(i,k)==true && thid(i,l)==false){
          dvec[0] = log(-v1[0]) + log(t(i,k))- v[0];
        } else if(thid(i,k)==false && thid(i,l)==true){
          dvec[0] = log(-v2[0]) + log(t(i,l))- v[0];
        } else if(thid(i,k)==true && thid(i,l)==true){
          dvec[0] = log(v1[0] * v2[0] - v12[0]) +
            log(t(i,k)) + log(t(i,l)) - v[0];
        } else{
          dvec[0] = zdn[0];
        }
        ll[0] = ll[0] - dvec[0];
        }
      }
      //Extra contribution from observations that are censored, but not in Y.
      if(!Rcpp::is_na(zdn)[0]){ //fixes a bug caused by lambda=1, which creates zdn=NaN
        ll[0] = ll[0]  - (N-y.nrow())*zdn[0];
      }
    }
  }

  return(ll);
}



//' Negative log-likelihood of scaled negative Dirichlet extreme value model
//'
//' This is the workhorse with the likelihood, which transforms observations
//' and returns the objective function value.
//'
//' @param y \code{n} by \code{d} matrix of observations on original scale
//' @param u threshold (on \code{y}-scale), to be substracted from y
//' @param thid \code{n} by \code{d} matrix of logical vector indicating whether or not a
//' threshold exceedance was observed
//' @param N total number of observations (not just exceedances)
//' @param scale vector of scale parameters
//' @param shape vector of shape parameters
//' @param lambda vector of percentage of threshold exceedances
//' @param alpha vector of parameters for the Dirichlet model
//' @param rho index of regular variation parameter, in \eqn{(0, \min(\alpha))}
//' @return the negative log-composite likelihood value
//' @section Note: The \code{N} argument is there mostly for cases in which function returns
//' observations stripped from their contribution; in such cases, the count
//' would otherwise be `incorrect'.
//' @export
// [[Rcpp::export]]
NumericVector nllmvsnegdir(NumericMatrix y, LogicalMatrix thid, int N,
                        NumericVector lambda, NumericVector u, NumericVector alpha,
                        double rho, NumericVector scale, NumericVector shape) {
  if(scale.size()==1){	scale = rep(scale, y.ncol());}
  if(shape.size()==1){	shape = rep(shape, y.ncol());}
  if(alpha.size()==1){	alpha = rep(alpha, y.ncol());} //Symmetric
  if(shape.size()!=y.ncol() || scale.size()!=y.ncol()){
    Rcpp::stop("Invalid marginal parameter vector size for the GPD");
  }
  NumericVector dvec(1);
  NumericVector ll(1);
  //Marginal transformation
  List MT = MarginalTransformation(y, u, lambda, scale, shape);
  if(MT["outOfDomain"] || rho > min(alpha) || rho < 0 || is_true(any(alpha>100)) | is_true(any(alpha<0))){
    ll[0] = 1e10; //This finite approximation of `Inf' is catched in the optim wrapper
    return ll;
  }
  NumericMatrix t = MT["t"];
  NumericMatrix data = MT["dat"];
  NumericVector lambda2 = -1.0/log(1.0-lambda);  //value of threshold on Frechet scale
  //Containers
  NumericVector x(1);   // transformed variable (only one point because angle on [0,1])
  NumericVector zdn(1); // contribution for points below both thresholds
  NumericVector v(1);   // exponent measure and derivatives
  NumericVector v1(1);  // these are set to zero by default
  NumericVector v2(1);  // and we overwrite them to save space
  NumericVector v12(1); // joint exceedances
  NumericVector logA(y.ncol());
  for(int i=0; i < y.ncol(); i++){
    logA[i]  = lgamma(alpha[i] - rho) - lgamma(alpha[i]);
  }
  for(int l=0; l < y.ncol(); l++) { //Loop over variable indices
    checkUserInterrupt();
    for(int k=0; k < l; k++) { // Loop over other pairs (use symmetry)
      //Recall l (aka "1") < k (aka "2") - zdn is correct below
      zdn[0] = exp((logA[k] + log(lambda2[k])) / rho) /
      (exp((logA[k] + log(lambda2[k])) / rho) +
        exp((logA[l] + log(lambda2[l])) / rho));
      //Contribution of values under both threshold
      // V - exponent measure  -- censored contributions, depends only on thresholds
      zdn[0] = -R::pbeta(1-zdn[0], alpha[k] - rho,  alpha[l], 1, 0) / lambda2[k] -
        R::pbeta(zdn[0], alpha[l] - rho, alpha[k], 1, 0) / lambda2[l];
      for(int i=0; i < y.nrow(); i++)  { // Loop over the observations
        if(!R_IsNA(data(i,k)) && !R_IsNA(data(i,l))){
        //If either is marginally exceeding, continue, else skip calculations
        if(thid(i,k)==true || thid(i,l)==true){
        //data is data on unit Frechet scale
        //t contains the marginal censored contributions
        //x are transformed observations
        x[0] = exp( (logA[k] + log(data(i,k))) / rho) /
          (exp((logA[k] + log(data(i,k))) / rho) +
            exp((logA[l] + log(data(i,l))) / rho));
          v[0] = exp(R::pbeta(1-x[0], alpha[k] - rho, alpha[l], 1, 1) - log(data(i,k))) +
            exp(R::pbeta(x[0], alpha[l] - rho, alpha[k], 1, 1) - log(data(i,l)));
          v1[0] = -R::pbeta(1-x[0],  alpha[k] - rho, alpha[l], 1, 0)/ pow(data(i,k),2);
          v2[0] = -R::pbeta(x[0], alpha[l] - rho, alpha[k], 1, 0)/ pow(data(i,l),2);
          v12[0] = - exp(lgamma(alpha[k]+alpha[l]-rho)-lgamma(alpha[k])-lgamma(alpha[l])-log(rho)-
            (alpha[k]+alpha[l]-rho)*log((exp(-(logA[k]+log(data(i,k)))/rho)+exp(-(logA[l]+log(data(i,l)))/rho)))-
            alpha[k]/rho*logA[k]-alpha[l]/rho*logA[l]-
            (alpha[l]/rho + 1.0)*log(data(i,l))-(alpha[k]/rho + 1.0)*log(data(i,k)));
        }
        //Function from evd: 1 (respectively 2) for marginal exceedance of x_1 (x_2), 3 for joint exceedance
        if(thid(i,k)==true && thid(i,l)==false){
          dvec[0] = log(-v1[0]) + log(t(i,k))- v[0];
        } else if(thid(i,k)==false && thid(i,l)==true){
          dvec[0] = log(-v2[0]) + log(t(i,l))- v[0];
        } else if(thid(i,k)==true && thid(i,l)==true){
          dvec[0] = log(v1[0] * v2[0] - v12[0]) +
            log(t(i,k)) + log(t(i,l)) - v[0];
        } else{
          dvec[0] = zdn[0];
        }
        ll[0] = ll[0] - dvec[0];
        }
      }
      //Extra contribution from observations that are censored, but not in Y.
      if(!Rcpp::is_na(zdn)[0]){ //fixes a bug caused by lambda=1, which creates zdn=NaN
        ll[0] = ll[0]  - (N-y.nrow())*zdn[0];
      }
    }
  }

  return(ll);
}



//' Negative log-likelihood of Husler-Reiss model
//'
//' This is the workhorse with the likelihood, which transforms observations
//' and returns the objective function value.
//'
//' @param y \code{n} by \code{d} matrix of observations on original scale
//' @param u threshold (on \code{y}-scale), to be substracted from y
//' @param thid \code{n} by \code{d} matrix of logical vector indicating whether or not a
//' threshold exceedance was observed
//' @param N total number of observations (not just exceedances)
//' @param scale vector of scale parameters
//' @param shape vector of shape parameters
//' @param lambda vector of percentage of threshold exceedances
//' @param Lambda matrix of parameters
//' @return the negative log-composite likelihood value
//' @section Note: The \code{N} argument is there mostly for cases in which function returns
//' observations stripped from their contribution; in such cases, the count
//' would otherwise be `incorrect'.
//' @export
// [[Rcpp::export]]
NumericVector nllmvhr(NumericMatrix y, LogicalMatrix thid, int N,
                      NumericVector lambda, NumericVector u, NumericMatrix Lambda,
                      NumericVector scale, NumericVector shape) {
  if(scale.size()==1){	scale = rep(scale, y.ncol());}
  if(shape.size()==1){	shape = rep(shape, y.ncol());}
  if(shape.size()!=y.ncol() || scale.size()!=y.ncol()){
    Rcpp::stop("Invalid marginal parameter vector size for the GPD");
  }
  NumericVector dvec(1);
  NumericVector ll(1);
  //Marginal transformation
  List MT = MarginalTransformation(y, u, lambda, scale, shape);
  if(MT["outOfDomain"]){
    ll[0] = 1e10; //This error code, `Inf' is catched in the optim wrapper
    return ll;
  }
  if(is_true(any(Lambda>5))){
    ll[0] = 1e10; //This error code, `Inf' is catched in the optim wrapper
    return ll;
  }
  NumericMatrix t = MT["t"];
  NumericMatrix data = MT["dat"];
  NumericVector lambda2 = -1.0/log(1.0-lambda);  //value of threshold on Frechet scale
  //Containers
  NumericVector zdn(1); // contribution for points below both thresholds
  NumericVector v(1);   // exponent measure and derivatives
  NumericVector v1(1);  // these are set to zero by default
  NumericVector v2(1);  // and we overwrite them to save space
  NumericVector v12(1); // joint exceedances
  for(int l=0; l < y.ncol(); l++) { //Loop over variable indices
    checkUserInterrupt();
    for(int k=0; k < l; k++) { // Loop over other pairs (use symmetry)
      // Recall l (aka "1") < k (aka "2") - zdn is correct below
      // Contribution of values under both threshold
      // V - exponent measure  -- censored contributions, depends only on thresholds
      zdn[0] = -V_BR(lambda2[k],lambda2[l],Lambda(k,l)*2.0);
      for(int i=0; i < y.nrow(); i++)  { // Loop over the observations
        if(!R_IsNA(data(i,k)) && !R_IsNA(data(i,l))){
          if(thid(i,k)==true || thid(i,l)==true){
            //data is data on unit Frechet scale
            //t contains the marginal censored contributions
            // v[0] = R::pnorm(Lambda(k,l)-2.0*(log(data(i,k))-log(data(i,l)))/Lambda(k,l),0,1,1,0) / data(i,k) +
            //   R::pnorm(Lambda(k,l)-2.0*(log(data(i,l))-log(data(i,k)]))/Lambda(k,l),0,1,1,0) / data(i,l);
            v[0] = V_BR(data(i,k),data(i,l),Lambda(k,l)*2.0);
            v1[0] = V_BR1(data(i,k),data(i,l),Lambda(k,l)*2.0);
            v2[0] = V_BR2(data(i,k),data(i,l),Lambda(k,l)*2.0);
            v12[0] = V_BR12(data(i,k),data(i,l),Lambda(k,l)*2.0);
          }
        //Function from evd: 1 (respectively 2) for marginal exceedance of x_1 (x_2), 3 for joint exceedance
        if(thid(i,k)==true && thid(i,l)==false){
          dvec[0] = log(-v1[0]) + log(t(i,k))- v[0];
        } else if(thid(i,k)==false && thid(i,l)==true){
          dvec[0] = log(-v2[0]) + log(t(i,l))- v[0];
        } else if(thid(i,k)==true && thid(i,l)==true){
          dvec[0] = log(v1[0] * v2[0] - v12[0]) +
            log(t(i,k)) + log(t(i,l)) - v[0];
        } else{
          dvec[0] = zdn[0];
        }
        ll[0] = ll[0] - dvec[0];
        }
      }
      //Extra contribution from observations that are censored, but not in Y.
      if(!Rcpp::is_na(zdn)[0]){ //fixes a bug caused by lambda=1, which creates zdn=NaN
        ll[0] = ll[0]  - (N-y.nrow())*zdn[0];
      }
    }
  }

  return(ll);
}

//' Negative log-likelihood of extremal Student model
//'
//' This is the workhorse with the likelihood, which transforms observations
//' and returns the objective function value.
//' @param y \code{n} by \code{d} matrix of observations on original scale
//' @param u threshold (on \code{y}-scale), to be substracted from y
//' @param thid \code{n} by \code{d} matrix of logical vector indicating whether or not a
//' threshold exceedance was observed
//' @param N total number of observations (not just exceedances)
//' @param scale vector of scale parameters
//' @param shape vector of shape parameters
//' @param lambda vector of percentage of threshold exceedances
//' @param Rho correlation matrix
//' @param nu degrees of freedom of the model
//' @return the negative log-composite likelihood value
//' @section Note: The \code{N} argument is there mostly for cases in which function returns
//' observations stripped from their contribution; in such cases, the count
//' would otherwise be `incorrect'.
//' @export
// [[Rcpp::export]]
NumericVector nllmvxstud(NumericMatrix y, LogicalMatrix thid, int N,
                          NumericVector lambda, NumericVector u, NumericMatrix Rho,
                          double nu, NumericVector scale, NumericVector shape) {
  if(scale.size()==1){	scale = rep(scale, y.ncol());}
  if(shape.size()==1){	shape = rep(shape, y.ncol());}
  if(shape.size()!=y.ncol() || scale.size()!=y.ncol()){
    Rcpp::stop("Invalid marginal parameter vector size for the GPD");
  }
  NumericVector dvec(1);
  NumericVector ll(1);
  //Marginal transformation
  List MT = MarginalTransformation(y, u, lambda, scale, shape);
  if(MT["outOfDomain"]){
    ll[0] = 1e10; //This error code, `Inf' is catched in the optim wrapper
    return ll;
  }
  if(nu < 0){
    ll[0] = 1e10; //This error code, `Inf' is catched in the optim wrapper
    return ll;
  }
  NumericMatrix t = MT["t"];
  NumericMatrix data = MT["dat"];
  NumericVector lambda2 = -1.0/log(1.0-lambda);  //value of threshold on Frechet scale
  double b;
  //Containers
  NumericVector zdn(1); // contribution for points below both thresholds
  NumericVector v(1);   // exponent measure and derivatives
  NumericVector v1(1);  // these are set to zero by default
  NumericVector v2(1);  // and we overwrite them to save space
  NumericVector v12(1); // joint exceedances
  for(int l=0; l < y.ncol(); l++) { //Loop over variable indices
    checkUserInterrupt();
    for(int k=0; k < l; k++) { // Loop over other pairs (use symmetry)
      b = pow((nu+1.0)/(1.0-pow(Rho(k,l),2)), 0.5);
      // V - exponent measure  -- censored contributions, depends only on thresholds
      zdn[0] = -V_EXST(lambda2(k), lambda2(l), b, Rho(k,l), nu);
      for(int i=0; i < y.nrow(); i++)  { // Loop over the observations
        if(!R_IsNA(data(i,k)) && !R_IsNA(data(i,l))){
          if(thid(i,k)==true || thid(i,l)==true){
            //data is data on unit Frechet scale
            //t contains the marginal censored contributions
            v[0] = V_EXST(data(i,k),data(i,l), b, Rho(k,l), nu);
            v1[0] = V_EXST1(data(i,k),data(i,l), b, Rho(k,l), nu);
            v2[0] = V_EXST2(data(i,k),data(i,l), b, Rho(k,l), nu);
            v12[0] = V_EXST12(data(i,k),data(i,l), b, Rho(k,l), nu);
          }
        //Function from evd: 1 (respectively 2) for marginal exceedance of x_1 (x_2), 3 for joint exceedance
        if(thid(i,k)==true && thid(i,l)==false){
          dvec[0] = log(-v1[0]) + log(t(i,k))- v[0];
        } else if(thid(i,k)==false && thid(i,l)==true){
          dvec[0] = log(-v2[0]) + log(t(i,l))- v[0];
        } else if(thid(i,k)==true && thid(i,l)==true){
          dvec[0] = log(v1[0] * v2[0] - v12[0]) +
            log(t(i,k)) + log(t(i,l)) - v[0];
        } else{
          dvec[0] = zdn[0];
        }
        ll[0] = ll[0] - dvec[0];
        }
      }
      //Extra contribution from observations that are censored, but not in Y.
      if(!Rcpp::is_na(zdn)[0]){ //fixes a bug caused by lambda=1, which creates zdn=NaN
        ll[0] = ll[0]  - (N-y.nrow())*zdn[0];
      }
    }
  }

  return(ll);
}

// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// MarginalTransformation
List MarginalTransformation(NumericMatrix data, NumericVector u, NumericVector lambda, NumericVector scale, NumericVector shape, bool GEVorGP);
RcppExport SEXP ExtLiouv_MarginalTransformation(SEXP dataSEXP, SEXP uSEXP, SEXP lambdaSEXP, SEXP scaleSEXP, SEXP shapeSEXP, SEXP GEVorGPSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type u(uSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type scale(scaleSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type shape(shapeSEXP);
    Rcpp::traits::input_parameter< bool >::type GEVorGP(GEVorGPSEXP);
    rcpp_result_gen = Rcpp::wrap(MarginalTransformation(data, u, lambda, scale, shape, GEVorGP));
    return rcpp_result_gen;
END_RCPP
}
// isAbove
LogicalMatrix isAbove(NumericMatrix data, NumericVector threshold);
RcppExport SEXP ExtLiouv_isAbove(SEXP dataSEXP, SEXP thresholdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type threshold(thresholdSEXP);
    rcpp_result_gen = Rcpp::wrap(isAbove(data, threshold));
    return rcpp_result_gen;
END_RCPP
}
// nllmvsdir
NumericVector nllmvsdir(NumericMatrix y, LogicalMatrix thid, int N, NumericVector lambda, NumericVector u, NumericVector alpha, double rho, NumericVector scale, NumericVector shape);
RcppExport SEXP ExtLiouv_nllmvsdir(SEXP ySEXP, SEXP thidSEXP, SEXP NSEXP, SEXP lambdaSEXP, SEXP uSEXP, SEXP alphaSEXP, SEXP rhoSEXP, SEXP scaleSEXP, SEXP shapeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type y(ySEXP);
    Rcpp::traits::input_parameter< LogicalMatrix >::type thid(thidSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type u(uSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type scale(scaleSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type shape(shapeSEXP);
    rcpp_result_gen = Rcpp::wrap(nllmvsdir(y, thid, N, lambda, u, alpha, rho, scale, shape));
    return rcpp_result_gen;
END_RCPP
}
// nllmvsnegdir
NumericVector nllmvsnegdir(NumericMatrix y, LogicalMatrix thid, int N, NumericVector lambda, NumericVector u, NumericVector alpha, double rho, NumericVector scale, NumericVector shape);
RcppExport SEXP ExtLiouv_nllmvsnegdir(SEXP ySEXP, SEXP thidSEXP, SEXP NSEXP, SEXP lambdaSEXP, SEXP uSEXP, SEXP alphaSEXP, SEXP rhoSEXP, SEXP scaleSEXP, SEXP shapeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type y(ySEXP);
    Rcpp::traits::input_parameter< LogicalMatrix >::type thid(thidSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type u(uSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type scale(scaleSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type shape(shapeSEXP);
    rcpp_result_gen = Rcpp::wrap(nllmvsnegdir(y, thid, N, lambda, u, alpha, rho, scale, shape));
    return rcpp_result_gen;
END_RCPP
}
// nllmvhr
NumericVector nllmvhr(NumericMatrix y, LogicalMatrix thid, int N, NumericVector lambda, NumericVector u, NumericMatrix Lambda, NumericVector scale, NumericVector shape);
RcppExport SEXP ExtLiouv_nllmvhr(SEXP ySEXP, SEXP thidSEXP, SEXP NSEXP, SEXP lambdaSEXP, SEXP uSEXP, SEXP LambdaSEXP, SEXP scaleSEXP, SEXP shapeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type y(ySEXP);
    Rcpp::traits::input_parameter< LogicalMatrix >::type thid(thidSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type u(uSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Lambda(LambdaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type scale(scaleSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type shape(shapeSEXP);
    rcpp_result_gen = Rcpp::wrap(nllmvhr(y, thid, N, lambda, u, Lambda, scale, shape));
    return rcpp_result_gen;
END_RCPP
}
// nllmvxstud
NumericVector nllmvxstud(NumericMatrix y, LogicalMatrix thid, int N, NumericVector lambda, NumericVector u, NumericMatrix Rho, double nu, NumericVector scale, NumericVector shape);
RcppExport SEXP ExtLiouv_nllmvxstud(SEXP ySEXP, SEXP thidSEXP, SEXP NSEXP, SEXP lambdaSEXP, SEXP uSEXP, SEXP RhoSEXP, SEXP nuSEXP, SEXP scaleSEXP, SEXP shapeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type y(ySEXP);
    Rcpp::traits::input_parameter< LogicalMatrix >::type thid(thidSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type u(uSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type Rho(RhoSEXP);
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type scale(scaleSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type shape(shapeSEXP);
    rcpp_result_gen = Rcpp::wrap(nllmvxstud(y, thid, N, lambda, u, Rho, nu, scale, shape));
    return rcpp_result_gen;
END_RCPP
}
// dirspecdens
NumericVector dirspecdens(NumericVector param, NumericMatrix dat, int d, bool transform);
RcppExport SEXP ExtLiouv_dirspecdens(SEXP paramSEXP, SEXP datSEXP, SEXP dSEXP, SEXP transformSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type param(paramSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type dat(datSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< bool >::type transform(transformSEXP);
    rcpp_result_gen = Rcpp::wrap(dirspecdens(param, dat, d, transform));
    return rcpp_result_gen;
END_RCPP
}
// negdirspecdens
NumericVector negdirspecdens(NumericVector param, NumericMatrix dat, int d, bool transform);
RcppExport SEXP ExtLiouv_negdirspecdens(SEXP paramSEXP, SEXP datSEXP, SEXP dSEXP, SEXP transformSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type param(paramSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type dat(datSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< bool >::type transform(transformSEXP);
    rcpp_result_gen = Rcpp::wrap(negdirspecdens(param, dat, d, transform));
    return rcpp_result_gen;
END_RCPP
}
// ctspecdens
NumericVector ctspecdens(NumericVector param, NumericMatrix dat, bool transform);
RcppExport SEXP ExtLiouv_ctspecdens(SEXP paramSEXP, SEXP datSEXP, SEXP transformSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type param(paramSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type dat(datSEXP);
    Rcpp::traits::input_parameter< bool >::type transform(transformSEXP);
    rcpp_result_gen = Rcpp::wrap(ctspecdens(param, dat, transform));
    return rcpp_result_gen;
END_RCPP
}

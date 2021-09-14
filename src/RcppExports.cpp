// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// LinTransf
List LinTransf(arma::mat mean, arma::cube sigma, arma::mat B, arma::mat Z, int G, int d);
RcppExport SEXP _ppgmmga_LinTransf(SEXP meanSEXP, SEXP sigmaSEXP, SEXP BSEXP, SEXP ZSEXP, SEXP GSEXP, SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type B(BSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< int >::type G(GSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(LinTransf(mean, sigma, B, Z, G, d));
    return rcpp_result_gen;
END_RCPP
}
// EntropyGauss
double EntropyGauss(arma::mat S, int d);
RcppExport SEXP _ppgmmga_EntropyGauss(SEXP SSEXP, SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type S(SSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(EntropyGauss(S, d));
    return rcpp_result_gen;
END_RCPP
}
// orth
arma::mat orth(arma::mat A, std::string method);
RcppExport SEXP _ppgmmga_orth(SEXP ASEXP, SEXP methodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< std::string >::type method(methodSEXP);
    rcpp_result_gen = Rcpp::wrap(orth(A, method));
    return rcpp_result_gen;
END_RCPP
}
// encodebasis
NumericMatrix encodebasis(NumericVector par, int d, int p);
RcppExport SEXP _ppgmmga_encodebasis(SEXP parSEXP, SEXP dSEXP, SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type par(parSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(encodebasis(par, d, p));
    return rcpp_result_gen;
END_RCPP
}
// logsumexp
double logsumexp(arma::vec x);
RcppExport SEXP _ppgmmga_logsumexp(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(logsumexp(x));
    return rcpp_result_gen;
END_RCPP
}
// EntropyUT
double EntropyUT(int G, arma::vec pro, arma::mat mean, arma::cube sigma, double d);
RcppExport SEXP _ppgmmga_EntropyUT(SEXP GSEXP, SEXP proSEXP, SEXP meanSEXP, SEXP sigmaSEXP, SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type G(GSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type pro(proSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< double >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(EntropyUT(G, pro, mean, sigma, d));
    return rcpp_result_gen;
END_RCPP
}
// EntropyVAR
double EntropyVAR(int G, NumericVector pro, NumericMatrix mean, arma::cube sigma, int d);
RcppExport SEXP _ppgmmga_EntropyVAR(SEXP GSEXP, SEXP proSEXP, SEXP meanSEXP, SEXP sigmaSEXP, SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type G(GSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type pro(proSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(EntropyVAR(G, pro, mean, sigma, d));
    return rcpp_result_gen;
END_RCPP
}
// EntropySOTE
double EntropySOTE(arma::mat data, int G, arma::vec pro, arma::mat mean, arma::cube sigma);
RcppExport SEXP _ppgmmga_EntropySOTE(SEXP dataSEXP, SEXP GSEXP, SEXP proSEXP, SEXP meanSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type G(GSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type pro(proSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(EntropySOTE(data, G, pro, mean, sigma));
    return rcpp_result_gen;
END_RCPP
}
// EntropyMCapprox
List EntropyMCapprox(arma::mat data, int G, arma::vec pro, arma::mat mean, arma::cube sigma);
RcppExport SEXP _ppgmmga_EntropyMCapprox(SEXP dataSEXP, SEXP GSEXP, SEXP proSEXP, SEXP meanSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type G(GSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type pro(proSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(EntropyMCapprox(data, G, pro, mean, sigma));
    return rcpp_result_gen;
END_RCPP
}
// EntropyGMM
double EntropyGMM(arma::mat data, arma::vec pro, arma::mat mean, arma::cube sigma);
RcppExport SEXP _ppgmmga_EntropyGMM(SEXP dataSEXP, SEXP proSEXP, SEXP meanSEXP, SEXP sigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type data(dataSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type pro(proSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type sigma(sigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(EntropyGMM(data, pro, mean, sigma));
    return rcpp_result_gen;
END_RCPP
}

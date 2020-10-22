// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// sample_pois
arma::umat sample_pois(int n, arma::mat theta, arma::mat gam, int burn_in, int spacing, bool verb);
RcppExport SEXP _pois_sample_pois(SEXP nSEXP, SEXP thetaSEXP, SEXP gamSEXP, SEXP burn_inSEXP, SEXP spacingSEXP, SEXP verbSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type gam(gamSEXP);
    Rcpp::traits::input_parameter< int >::type burn_in(burn_inSEXP);
    Rcpp::traits::input_parameter< int >::type spacing(spacingSEXP);
    Rcpp::traits::input_parameter< bool >::type verb(verbSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_pois(n, theta, gam, burn_in, spacing, verb));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_pois_sample_pois", (DL_FUNC) &_pois_sample_pois, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_pois(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

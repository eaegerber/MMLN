// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]
#include <RcppArmadillo.h>
#include <RcppEigen.h>

using namespace Rcpp; // RcppArmadillo.h includes Rcpp.h interally
using namespace arma;

// Aliases for common matrix/vector types
typedef arma::mat Mat;
typedef arma::vec Vec;
typedef arma::rowvec RowVec;
typedef arma::uvec UVec;

// Verify that the C++ backend is active and callable from R
// [[Rcpp::export]]
void verify_installation()
{
  Rcpp::Rcout << "C++ is working w/ Rcpp" << std::endl;
}
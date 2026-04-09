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

// Some mock utility functions, may be removed in the future or could be used.

// Verify that the C++ backend is active and callable from R
// [[Rcpp::export]]
void verify_installation()
{
  Rcpp::Rcout << "C++ is working w/ Rcpp" << std::endl;
}

// Clamp each element of a vector to [lo, hi], using RcppArmadillo
// [[Rcpp::export]]
arma::vec clamp_vec(arma::vec x, double lo, double hi)
{
  return arma::clamp(x, lo, hi);
}

// Compute column-wise means of a matrix, using RcppArmadillo
// [[Rcpp::export]]
arma::rowvec col_means(arma::mat X)
{
  return arma::mean(X, 0);
}

// Compute eigenvalues of a symmetric matrix, using RcppEigen
// [[Rcpp::export]]
Eigen::VectorXd sym_eigenvalues(Eigen::MatrixXd X)
{
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(X);
  return solver.eigenvalues();
}

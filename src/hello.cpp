#include <Rcpp.h>
using namespace Rcpp;

// DEV NOTE:
// {[[Rcpp::export]]} Makes this function callable from R.
// Remove this annotation for any internal C++ helpers that R doesn't need to call directly.

// [[Rcpp::export]]
void hello_world()
{
  Rcpp::Rcout << "C++ is working w/ Rcpp" << std::endl;
}
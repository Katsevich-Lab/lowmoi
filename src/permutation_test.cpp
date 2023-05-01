#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector low_level_permutation_test(NumericVector y, IntegerMatrix index_mat) {
  int n_tests = index_mat.ncol();
  int n_treatment = index_mat.nrow();
  double sum;
  NumericVector out(n_tests);

  for (int j = 0; j < n_tests; j ++) {
    sum = 0;
    for (int i = 0; i < n_treatment; i ++) {
      sum += y[index_mat(i, j)];
    }
    out[j] = sum;
  }
  return(out);
}

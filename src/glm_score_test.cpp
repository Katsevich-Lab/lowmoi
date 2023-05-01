#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector low_level_score_test_vectorized(NumericVector a, NumericMatrix B, NumericVector w, IntegerMatrix index_mat) {
  int n_tests = index_mat.ncol();
  int n_treatment = index_mat.nrow();
  int dim = B.nrow();
  int i, j, k, idx;
  double lower_right, inner_sum, lower_left, top;
  NumericVector out(n_tests);

  // iterate over the n_tests index vectors
  for (k = 0; k < n_tests; k ++) {
    lower_right = 0;
    inner_sum = 0;
    // iterate over the rows of B
    for (i = 0; i < dim; i ++) {
      inner_sum = 0;
      for (j = 0; j < n_treatment; j ++) {
        inner_sum += B(i, index_mat(j, k));
      }
      lower_right += inner_sum * inner_sum;
    }

    // second, compute the lower-left hand of the denominator; also, compute the top
    lower_left = 0;
    top = 0;
    for (j = 0; j < n_treatment; j ++) {
      idx = index_mat(j, k);
      top += a[idx];
      lower_left += w[idx];
    }

    // finally, compute the z-score
    out[k] = top/sqrt(lower_left - lower_right);
  }

  return out;
}

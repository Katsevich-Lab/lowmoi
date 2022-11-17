#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector run_mw_test_cpp(int n_x, int n_y, NumericVector combined, IntegerMatrix synthetic_treatment_idxs) {
  // define quantities
  int nrow = synthetic_treatment_idxs.nrow(), ncol = synthetic_treatment_idxs.ncol();
  int n_cells = n_x + n_y;
  int trt_pos, ctrl_pos;
  Function rank("rank");
  Function table("table");
  NumericVector r(combined.length());
  NumericVector out(ncol);
  NumericVector curr_combined(combined.length());
  IntegerVector logical_v(n_cells);
  double statistic = 0, z, top_sum, sigma, correction;
  int length_nties, idx;

  for (int j = 0; j < ncol; j ++) {
    for (int i = 0; i < n_cells; i ++) logical_v[i] = 1;
    for (int i = 0; i < n_x; i ++) {
      idx = synthetic_treatment_idxs(i, j) - 1;
      curr_combined[i] = combined[idx];
      logical_v[idx] = 0;
    }
    idx = n_x;
    for (int i = 0; i < n_cells; i ++) {
      if (logical_v[i]) {
        curr_combined[idx ++] = combined[i];
      }
    }

    statistic = 0;
    r = rank(curr_combined);
    for (int i = 0; i < n_x; i ++) {
      statistic += r[i];
    }
    statistic -= n_x * (n_x + 1)/2;
    NumericVector nties = table(r);
    length_nties = nties.length();
    z = statistic - n_x * n_y/2;
    top_sum = 0;
    for (int i = 0; i < length_nties; i ++) {
      top_sum += nties[i] * (nties[i] * nties[i] - 1);
    }
    sigma = sqrt((n_x * n_y/12) * ((n_x + n_y + 1) - top_sum/((n_x + n_y) * (n_x + n_y - 1))));
    if (z < 0) {
      correction = - 0.5;
    } else {
      correction = 0.5;
    }
    z = (z - correction)/sigma;
    out[j] = z;
  }
  return(out);
}

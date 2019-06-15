#include <Rcpp.h>
using namespace Rcpp;

//' Predict honest weights (C++)
//'
//' Computes honest weights from the random forest as in Wager & Athey (2019)
//' for the test sample based on the honest training sample
//'
//' @param x leaf_IDs_test - list of leaf IDs in test data
//' @param y leaf_IDs - list of leaf IDs in honest data
//' @param z leaf_size - list of leaf sizes in honest data
//' @param w binary indicator - equal 1 if marginal effects are being computed, 0 otherwise for normal prediction
// [[Rcpp::export]]
NumericVector pred_weights_C(List x, List y, List z, int w) {

  int nlist = x.size();

  NumericVector f_rows = as<NumericVector>(x[1]);
  NumericVector f_cols = as<NumericVector>(y[1]);

  int nf_rows = f_rows.size();
  int nf_cols = f_cols.size();


  NumericMatrix forest_out(nf_rows, nf_cols);

  for(int l = 0; l < nlist; ++l) {

            NumericVector leaf_IDs_pred = as<NumericVector>(x[l]);
            NumericVector leaf_IDs      = as<NumericVector>(y[l]);
            NumericVector leaf_size     = as<NumericVector>(z[l]);

            int n_leaf_IDs_pred = leaf_IDs_pred.size();
            int n_leaf_IDs      = leaf_IDs.size();

            NumericMatrix tree_out(n_leaf_IDs_pred, n_leaf_IDs);

            for(int i = 0; i < n_leaf_IDs_pred; ++i) {
            for(int j = 0; j < n_leaf_IDs; ++j) {

            tree_out(i,j) = leaf_IDs_pred[i]==leaf_IDs[j];
            tree_out(i,j) = tree_out(i,j)/leaf_size[j];

            }
            }

            for(int i = 0; i < n_leaf_IDs_pred; ++i) {
              for(int j = 0; j < n_leaf_IDs; ++j) {

                forest_out(i,j) = forest_out(i,j) + tree_out(i,j);

              }
            }

  }

  for(int i = 0; i < nf_rows; ++i) {
    for(int j = 0; j < nf_cols; ++j) {

      forest_out(i,j) = forest_out(i,j)/nlist ;

    }
  }

  // take colmeans of weights for saving memory in R if ME are being computed

  if (w == 1) {

    NumericMatrix forest_out_mean(1, nf_cols); // create output matrix for mean weights

    for(int i = 0; i < nf_cols; ++i) {

      forest_out_mean(0, i) = mean(forest_out(_, i));

    }

    return forest_out_mean;

  } else {

    return forest_out;

  }

}

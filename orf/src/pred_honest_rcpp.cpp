#include <Rcpp.h>
using namespace Rcpp;

//' Predict honest predictions (C++)
//'
//' Computes honest predictions from the random forest for a test sample based
//' on the honest training sample
//'
//' @param x unique_leaves (List)[ntree]
//' @param y honest_y (NumericVector)[nrow]
//' @param z honest_leaves (NumericMatrix)[nrow, ntree]
//' @param w test_leaves (NumericMatrix)[nrow, ntree]
//' @keywords internal
// [[Rcpp::export]]
NumericVector pred_honest_C(List x, NumericVector y, NumericMatrix z, NumericMatrix w) {

            // x - unique_leaves (List)[ntree]
            // y - honest_y (NumericVector)[nrow]
            // z - honest_leaves (NumericMatrix)[nrow, ntree]
            // w - test_leaves (NumericMatrix)[nrow, ntree]

            // generate empty variables
            int leaf_ID = 0;
            // how long are unique leaves list
            int n_unique_leaves = x.size();
            // how many rows does data have
            int n_rows = z.nrow();
            int n_rows_test = w.nrow();
            // for which rows are equal leaf_ID
            NumericVector obs_same_all(n_rows);
            NumericVector obs_same_all_test(n_rows_test);
            NumericVector y_same(n_rows);
            double y_mean = 0;

            // output matrix test_pred[nrow_test , ntree]
            NumericMatrix test_pred(n_rows_test, n_unique_leaves);
            NumericVector test_pred_final(n_rows_test);

            // loop over trees
            for(int tree_idx = 0; tree_idx < n_unique_leaves; ++tree_idx) {

              // how long are unique leaves for each tree
              NumericVector leaves = as<NumericVector>(x[tree_idx]);
              int n_leaves = leaves.size();

              // now loop over the unique leaves in the respective tree
              for(int leaf_idx = 0; leaf_idx < n_leaves; ++leaf_idx) {

                // take out leaf ID
                leaf_ID = leaves[leaf_idx];

                // find observations in the same leaf_ID
                for(int row_idx = 0; row_idx < n_rows; ++row_idx) {
                  obs_same_all[row_idx] = z(row_idx, tree_idx) == leaf_ID;
                }
                // for test sample (w)
                for(int row_idx = 0; row_idx < n_rows_test; ++row_idx) {
                  obs_same_all_test[row_idx] = w(row_idx, tree_idx) == leaf_ID;
                }

                // loop over all observations to get the same leaf_IDs
                for(int row_idx = 0; row_idx < n_rows; ++row_idx) {
                  if (obs_same_all[row_idx] == 1) {
                    y_same[row_idx] = y[row_idx];
                  } else {
                    y_same[row_idx] = 0;
                  }
                }

                // set variables
                double y_count = 0;
                double y_sum = 0;

                // loop over all observations to get the sum
                for(int row_idx = 0; row_idx < n_rows; ++row_idx) {
                  if (obs_same_all[row_idx] == 1) {
                    y_count = y_count + obs_same_all[row_idx];
                    y_sum = y_sum + y_same[row_idx];
                  }
                }

                // compute the mean
                y_mean = y_sum/y_count;

                // loop over all observations and feed in mean predictions for test
                for(int row_idx = 0; row_idx < n_rows_test; ++row_idx) {
                  if (obs_same_all_test[row_idx] == 1) {
                    test_pred(row_idx, tree_idx) = y_mean;
                  }
                }


              }

              // check for user interruptions
              Rcpp::checkUserInterrupt();

            }


            // now compute rowmeans for the test_pred matrix
            for (int row_idx = 0; row_idx < n_rows_test; ++row_idx) {
              double total = 0;
              // now go over columns
              for (int col_idx = 0; col_idx < n_unique_leaves; ++col_idx) {
                total += test_pred(row_idx, col_idx);
              }
              test_pred_final[row_idx] = total/n_unique_leaves;
            }

            return test_pred_final;

}

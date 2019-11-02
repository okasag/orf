// -----------------------------------------------------------------------------
// This file is part of orf.
//
// orf is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// orf is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with orf. If not, see <http://www.gnu.org/licenses/>.
//
// Written by:
//
// Gabriel Okasa
// Swiss Institute for Empirical Economic Research
// University of St.Gallen
// Varnbüelstrasse 14
// 9000 St.Gallen
// Switzerland
// -----------------------------------------------------------------------------

#include <Rcpp.h>
using namespace Rcpp;

//' Get honest predictions (C++)
//'
//' Computes honest predictions (fitted values) from the random forest
//' for the train and honest sample based on the honest training sample
//'
//' @param x unique_leaves (List)[ntree]
//' @param y honest_y (NumericVector)[nrow]
//' @param z honest_leaves (NumericMatrix)[nrow, ntree]
//' @param w train_leaves (NumericMatrix)[nrow, ntree]
//' @keywords internal
// [[Rcpp::export]]
NumericVector get_honest_C(List x, NumericVector y, NumericMatrix z, NumericMatrix w) {

            // x - unique_leaves (List)[ntree]
            // y - honest_y (NumericVector)[nrow]
            // z - honest_leaves (NumericMatrix)[nrow, ntree]
            // w - train_leaves (NumericMatrix)[nrow, ntree]

            // generate empty variables
            int leaf_ID = 0;
            // how long are unique leaves list
            int n_unique_leaves = x.size();
            // how many rows does data have
            int n_rows = z.nrow();
            int n_rows_train = w.nrow();
            // for which rows are equal leaf_ID
            NumericVector obs_same_all(n_rows);
            NumericVector obs_same_all_train(n_rows_train);
            NumericVector y_same(n_rows);
            double y_mean = 0;

            // output matrix honest_pred[nrow , ntree]
            NumericMatrix honest_pred(n_rows, n_unique_leaves);
            NumericMatrix train_pred(n_rows_train, n_unique_leaves);
            NumericMatrix all_pred(n_rows+n_rows_train, n_unique_leaves);
            NumericVector all_pred_final(n_rows+n_rows_train);



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

            for(int row_idx = 0; row_idx < n_rows_train; ++row_idx) {
            obs_same_all_train[row_idx] = w(row_idx, tree_idx) == leaf_ID; // for train sample (w)
            }


            for(int row_idx = 0; row_idx < n_rows; ++row_idx) {
            if (obs_same_all[row_idx] == 1) {
            y_same[row_idx] = y[row_idx];
            } else {
            y_same[row_idx] = 0;
            }
            }

            double y_count = 0;
            double y_sum = 0;
            for(int row_idx = 0; row_idx < n_rows; ++row_idx) {
            if (obs_same_all[row_idx] == 1) {
            y_count = y_count + obs_same_all[row_idx];
            y_sum = y_sum + y_same[row_idx];
            }
            }

            y_mean = y_sum/y_count;

            for(int row_idx = 0; row_idx < n_rows; ++row_idx) {
            if (obs_same_all[row_idx] == 1) {
            honest_pred(row_idx, tree_idx) = y_mean;
            }
            }

            // the same for train (lieber n_rows_train!)
            for(int row_idx = 0; row_idx < n_rows_train; ++row_idx) {
            if (obs_same_all_train[row_idx] == 1) {
            train_pred(row_idx, tree_idx) = y_mean;
            }
            }


            }

            }

            // PUT MATRICES TOGETHER
            for (int row_idx = 0; row_idx < n_rows + n_rows_train; ++row_idx) {
              if (row_idx < n_rows) {
                all_pred(row_idx,_) = honest_pred(row_idx,_);
              } else {
                all_pred(row_idx,_) = train_pred(row_idx-n_rows,_);
              }
            }

            // now compute rowmeans for the all_pred matrix
            for (int row_idx = 0; row_idx < n_rows + n_rows_train; ++row_idx) {
                double total = 0;
                // now go over columns
                for (int col_idx = 0; col_idx < n_unique_leaves; ++col_idx) {
                  total += all_pred(row_idx, col_idx);
                }
                all_pred_final[row_idx] = total/n_unique_leaves;
              }

            return all_pred_final;

            }

#include <Rcpp.h>
using namespace Rcpp;

//' Get honest weights (C++)
//'
//' Computes honest weights from the random forest as in Wager & Athey (2019)
//' for the train and honest sample based on the honest training sample
//'
//' @param x leaf_IDs_train - list of leaf IDs in train data
//' @param y leaf_IDs - list of leaf IDs in honest data
//' @param z leaf_size - list of leaf sizes in honest data
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix get_weights_C(List x, List y, List z) {

// leaf_IDs_train - list of leaf IDs in train data
// leaf_IDs - list of leaf IDs in honest data
// leaf_size - list of leaf sizes in honest data

  int nlist = x.size(); // number of trees

  NumericVector f_rows = as<NumericVector>(x[1]);
  NumericVector f_cols = as<NumericVector>(y[1]);

  int nf_rows = f_rows.size(); // rows of train data
  int nf_cols = f_cols.size(); // columns as rows of honest data

  // matrix where rows are rows of train and cols are rows of honest
  NumericMatrix forest_out_train(nf_rows, nf_cols);
  // matrix where rows are rows of honest and cols are rows of honest
  NumericMatrix forest_out_honest(nf_cols, nf_cols);
  // matrix where rows are rows of honest and train and cols are rows of honest
  NumericMatrix forest_out_all(nf_cols+nf_rows, nf_cols);

  // now loop over trees
  for(int l = 0; l < nlist; ++l) {

	          // take elements of lists as vectors
            NumericVector leaf_IDs_train = as<NumericVector>(x[l]);
            NumericVector leaf_IDs       = as<NumericVector>(y[l]);
            NumericVector leaf_size      = as<NumericVector>(z[l]);

            int n_leaf_IDs_train = leaf_IDs_train.size(); // how many different leaves are in train
            int n_leaf_IDs       = leaf_IDs.size(); // how many different leaves are in honest

	          // matrix where rows are obs in train and cols are obs in honest
            NumericMatrix tree_out_train(n_leaf_IDs_train, n_leaf_IDs);
	          // matrix where rows are obs in honest and cols are obs in honest
            NumericMatrix tree_out_honest(n_leaf_IDs, n_leaf_IDs);

	          // loops to go element by element and check the equality of leaves in each tree for train
            for(int i = 0; i < n_leaf_IDs_train; ++i) {
            	for(int j = 0; j < n_leaf_IDs; ++j) {

        		    tree_out_train(i,j) = leaf_IDs_train[i]==leaf_IDs[j]; // are leaves equal
        		    tree_out_train(i,j) = tree_out_train(i,j)/leaf_size[j]; // normalize by leaf size

            	}
            }

	          // loop to add each tree weight to overall forest weight
            for(int i = 0; i < n_leaf_IDs_train; ++i) {
              for(int j = 0; j < n_leaf_IDs; ++j) {

                forest_out_train(i,j) = forest_out_train(i,j) + tree_out_train(i,j);

              }
            }
	          // train sample done

      	    // now do the same weight computation for honest sample
      	    // loops to go element by element and check the equality of leaves in each tree for honest
            for(int i = 0; i < n_leaf_IDs; ++i) {
            	for(int j = 0; j < n_leaf_IDs; ++j) {

        		    tree_out_honest(i,j) = leaf_IDs[i]==leaf_IDs[j]; // are leaves equal
        		    tree_out_honest(i,j) = tree_out_honest(i,j)/leaf_size[j]; // normalize by leaf size

            	}
            }

	          // loop to add each tree weight to overall forest weight
            for(int i = 0; i < n_leaf_IDs; ++i) {
              for(int j = 0; j < n_leaf_IDs; ++j) {

                forest_out_honest(i,j) = forest_out_honest(i,j) + tree_out_honest(i,j);

              }
            }
	          // honest sample done

	          // check for user interruptions
	          Rcpp::checkUserInterrupt();
  }

  // loop to divide each element of a matrix by number of trees to get mean weights
  // train sample
  for(int i = 0; i < nf_rows; ++i) {
    for(int j = 0; j < nf_cols; ++j) {

      forest_out_train(i,j) = forest_out_train(i,j)/nlist ;

    }
  }

  // honest sample (not the same indices as train, be careful)
  for(int i = 0; i < nf_cols; ++i) {
    for(int j = 0; j < nf_cols; ++j) {

      forest_out_honest(i,j) = forest_out_honest(i,j)/nlist ;

    }
  }

  // now put these two matrices together (rbind(honest, train))
  for (int row_idx = 0; row_idx < nf_cols + nf_rows; ++row_idx) {
              if (row_idx < nf_cols) {
                forest_out_all(row_idx,_) = forest_out_honest(row_idx,_);
              } else {
                forest_out_all(row_idx,_) = forest_out_train(row_idx-nf_cols,_);
              }
  }

  return forest_out_all;

}

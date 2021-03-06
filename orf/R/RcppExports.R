# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' Get honest predictions (C++)
#'
#' Computes honest predictions (fitted values) from the random forest
#' for the train and honest sample based on the honest training sample
#'
#' @param x unique_leaves (List)[ntree]
#' @param y honest_y (NumericVector)[nrow]
#' @param z honest_leaves (NumericMatrix)[nrow, ntree]
#' @param w train_leaves (NumericMatrix)[nrow, ntree]
#' @keywords internal
get_honest_C <- function(x, y, z, w) {
    .Call(`_orf_get_honest_C`, x, y, z, w)
}

#' Get honest weights (C++)
#'
#' Computes honest weights from the random forest as in Wager & Athey (2019)
#' for the train and honest sample based on the honest training sample
#'
#' @param x leaf_IDs_train - list of leaf IDs in train data
#' @param y leaf_IDs - list of leaf IDs in honest data
#' @param z leaf_size - list of leaf sizes in honest data
#' @keywords internal
get_weights_C <- function(x, y, z) {
    .Call(`_orf_get_weights_C`, x, y, z)
}

#' Predict honest predictions (C++)
#'
#' Computes honest predictions from the random forest for a test sample based
#' on the honest training sample
#'
#' @param x unique_leaves (List)[ntree]
#' @param y honest_y (NumericVector)[nrow]
#' @param z honest_leaves (NumericMatrix)[nrow, ntree]
#' @param w test_leaves (NumericMatrix)[nrow, ntree]
#' @keywords internal
pred_honest_C <- function(x, y, z, w) {
    .Call(`_orf_pred_honest_C`, x, y, z, w)
}

#' Predict honest weights (C++)
#'
#' Computes honest weights from the random forest as in Wager & Athey (2019)
#' for the test sample based on the honest training sample
#'
#' @param x leaf_IDs_test - list of leaf IDs in test data
#' @param y leaf_IDs - list of leaf IDs in honest data
#' @param z leaf_size - list of leaf sizes in honest data
#' @param w binary indicator - equal 1 if marginal effects are being computed, 0 otherwise for normal prediction
#' @keywords internal
pred_weights_C <- function(x, y, z, w) {
    .Call(`_orf_pred_weights_C`, x, y, z, w)
}


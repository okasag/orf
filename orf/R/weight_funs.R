#' Get Forest Weights
#'
#' get forest weights, i.e. in-sample weights based on honest or train sample
#' produced by the random forest algorithm as defined in Wager & Athey (2018)
#'
#' @param forest estimated forest object of type ranger
#' @param honest_data honest dataframe
#' @param train_data train dataframe
#'
#' @importFrom stats predict ave
#' @import ranger
#'
#' @return matrix of honest forest weights
#'
#' @keywords internal
#'
get_forest_weights <- function(forest, honest_data, train_data) {

  # --------------------------------------------------------------------------------------------------- #

  # first get terminal nodes, i.e get terminal nodes for each obs and each tree
  # run your honest data through the forest structure and look where your observations end up
  leaf_IDs <- predict(forest, honest_data, type = "terminalNodes")$predictions
  # put leaf_IDs into a list (one element for one tree)
  leaf_IDs <- lapply(seq_along(leaf_IDs[1, ]), function(i) leaf_IDs[, i])

  ### formula for single tree: 1(Xi in L(x))/abs(L(x)) - normalized by the size of the leaf
  ## pre-requisites: leaf size and leaf ID
  # get the leaf size as counts of observations in leaves
  leaf_size <- lapply(leaf_IDs, function(x) ave(x, x, FUN=length))

  # now do the same for your training data
  leaf_IDs_train <- predict(forest, train_data, type = "terminalNodes")$predictions
  # put leaf_IDs into a list
  leaf_IDs_train <- lapply(seq_along(leaf_IDs_train[1, ]), function(i) leaf_IDs_train[, i])

  # --------------------------------------------------------------------------------------------------- #

  # get weights for the whole in sample (train + honest) data
  forest_weights <- get_weights_C(leaf_IDs_train, leaf_IDs, leaf_size)

  # --------------------------------------------------------------------------------------------------- #

  ## order everything back according to whole in-sample data
  # combine the two into a complete dataset (first honest rownames, then train rownames)
  rownames(forest_weights) <- c(rownames(honest_data), rownames(train_data)) # put original rownames in
  # order
  forest_weights <- as.matrix(forest_weights[order(as.numeric(row.names(forest_weights))), ])

  # --------------------------------------------------------------------------------------------------- #

  ## return forest weight final matrix
  return(forest_weights)

  # --------------------------------------------------------------------------------------------------- #

}


#' Predict Forest Weights
#'
#' predict forest weights, i.e. out-of-sample weights based on honest or train
#' sample produced by the random forest algorithm as defined in Wager & Athey (2018)
#'
#' @param forest estimated forest object of type ranger
#' @param data train (honest) dataframe
#' @param pred_data prediction dataframe
#'
#' @importFrom stats predict ave
#' @import ranger
#'
#' @return matrix of honest forest weights
#'
#' @keywords internal
#'
predict_forest_weights <- function(forest, data, pred_data) {

  # --------------------------------------------------------------------------------------------------- #

  # first get terminal nodes, i.e get terminal nodes for each obs and each tree
  # run your (new) data through the forest structure and look where your observations end up
  leaf_IDs <- predict(forest, data, type = "terminalNodes")$predictions
  # put leaf_IDs into a list (one element for one tree)
  leaf_IDs <- lapply(seq_along(leaf_IDs[1, ]), function(i) leaf_IDs[, i])
  # get the leaf size as counts of observations in leaves
  leaf_size <- lapply(leaf_IDs, function(x) ave(x, x, FUN=length))

  # now do the same for your prediction data
  leaf_IDs_pred <- predict(forest, pred_data, type = "terminalNodes")$predictions
  # put leaf_IDs into a list
  leaf_IDs_pred <- lapply(seq_along(leaf_IDs_pred[1, ]), function(i) leaf_IDs_pred[, i])

  # --------------------------------------------------------------------------------------------------- #

  # predict weights for the whole prediction sample based on train/honest data
  forest_weights <- pred_weights_C(leaf_IDs_pred, leaf_IDs, leaf_size, 0)

  # --------------------------------------------------------------------------------------------------- #

  # return result
  return(forest_weights)

  # --------------------------------------------------------------------------------------------------- #

}


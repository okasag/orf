#' Get Honest Predictions
#'
#' get honest prediction, i.e. fitted values for in sample data
#' (train and honest sample) based on the honest sample
#'
#' @param forest estimated forest object of type ranger
#' @param honest_data honest dataframe
#' @param train_data train dataframe
#'
#' @importFrom stats predict
#'
#' @return vector of honest forest predictions
get_honest <- function(forest, honest_data, train_data) {

  # needed inputs for the function: forest      - estimated forest object of type ranger
  #                                 honest_data - honest dataframe
  #                                 train_data  - train dataframe

  # ----------------------------------------------------------------------------------- #

  # get terminal nodes for the honest sample, i.e. run honest sample through forest
  honest_leaves <- predict(forest, honest_data, type = "terminalNodes")$predictions
  train_leaves <- predict(forest, train_data, type = "terminalNodes")$predictions
  # filter out unique leaves for each tree
  unique_leaves <- apply(honest_leaves, 2, unique)
  unique_leaves_train <- apply(train_leaves, 2, unique)
  # get honest outcomes
  honest_y <- as.numeric(honest_data[, 1])

  # --------------------------------------------------------------------------------------------------- #

  #### Rcpp implementation ####
  honest_fitted_values <- as.matrix(get_honest_C(unique_leaves, honest_y, honest_leaves, train_leaves))

  # --------------------------------------------------------------------------------------------------- #

  # combine the two into a complete dataset (first honest rownames, then train rownames)
  rownames(honest_fitted_values) <- c(rownames(honest_data), rownames(train_data)) # put original rownames in
  #forest_fitted_values <- rbind(train_fitted_values, honest_fitted_values) # put the sample together
  forest_fitted_values <- honest_fitted_values[order(as.numeric(row.names(honest_fitted_values))), ] # order according to rownames

  # put it into output
  output <- as.numeric(forest_fitted_values)
  #names(output) <- "honestPred"

  # --------------------------------------------------------------------------------------------------- #

  # return the output
  return(output)

  # --------------------------------------------------------------------------------------------------- #

}


#' Predict Honest Predictions
#'
#' predict honest prediction for out of sample data
#' based on the honest sample
#'
#' @param forest estimated forest object of type ranger
#' @param honest_data honest dataframe
#' @param test_data test dataframe
#'
#' @importFrom stats predict
#'
#' @return vector of honest forest predictions
predict_honest <- function(forest, honest_data, test_data) {

  # needed inputs for the function: forest - estimated forest object of type ranger
  #                                 honest_data - honest dataframe
  #                                 test_data - test dataframe

  # ----------------------------------------------------------------------------------- #

  # get terminal nodes for the honest sample, i.e. run honest sample through forest
  honest_leaves <- predict(forest, honest_data, type = "terminalNodes")$predictions
  test_leaves <- predict(forest, test_data, type = "terminalNodes")$predictions
  # filter out unique leaves for each tree
  unique_leaves <- apply(honest_leaves, 2, unique)
  #unique_leaves_test <- apply(test_leaves, 2, unique)
  # take honest outcomes
  honest_y <- as.numeric(honest_data[, 1])

  # ----------------------------------------------------------------------------------- #

  # new Rcpp prediction
  honest_pred <- pred_honest_C(unique_leaves, honest_y, honest_leaves, test_leaves)

  # --------------------------------------------------------------------------------------------------- #

  # put it into output
  output <- as.numeric(honest_pred)

  # return the output
  return(output)

  # ----------------------------------------------------------------------------------- #

}

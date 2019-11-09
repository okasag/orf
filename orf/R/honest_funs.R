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
#' @import ranger
#'
#' @return vector of honest forest predictions
#'
#' @keywords internal
#'
get_honest <- function(forest, honest_data, train_data) {

  # ----------------------------------------------------------------------------------- #

  # get terminal nodes for the honest sample, i.e. run honest sample through forest
  honest_leaves <- predict(forest, honest_data, type = "terminalNodes")$predictions
  train_leaves  <- predict(forest, train_data, type = "terminalNodes")$predictions
  # filter out unique leaves for each tree
  unique_leaves <- apply(honest_leaves, 2, unique)
  unique_leaves_train <- apply(train_leaves, 2, unique)
  # get honest outcomes
  honest_y <- as.numeric(honest_data[, 1])

  # --------------------------------------------------------------------------------------------------- #

  # Rcpp implementation
  honest_fitted_values <- as.matrix(get_honest_C(unique_leaves, honest_y, honest_leaves, train_leaves))

  # --------------------------------------------------------------------------------------------------- #

  # combine the two into a complete dataset (first honest rownames, then train rownames)
  rownames(honest_fitted_values) <- c(rownames(honest_data), rownames(train_data)) # put original rownames in
  # order according to rownames
  forest_fitted_values <- honest_fitted_values[order(as.numeric(row.names(honest_fitted_values))), ]

  # put it into output
  output <- as.numeric(forest_fitted_values)

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
#'
#' @keywords internal
#'
predict_honest <- function(forest, honest_data, test_data) {

  # ----------------------------------------------------------------------------------- #

  # get terminal nodes for the honest sample, i.e. run honest sample through forest
  honest_leaves <- predict(forest, honest_data, type = "terminalNodes")$predictions
  test_leaves   <- predict(forest, test_data, type = "terminalNodes")$predictions
  # filter out unique leaves for each tree
  unique_leaves <- apply(honest_leaves, 2, unique)
  # take honest outcomes
  honest_y <- as.numeric(honest_data[, 1])

  # ----------------------------------------------------------------------------------- #

  # new Rcpp prediction
  honest_pred <- pred_honest_C(unique_leaves, honest_y, honest_leaves, test_leaves)

  # ------------------------------------------------------------------------------------ #

  # put it into output
  output <- as.numeric(honest_pred)

  # return the output
  return(output)

  # ----------------------------------------------------------------------------------- #

}


#' honest sample split
#'
#' Creates honest sample split by randomly selecting prespecified share of observations
#' to belong to honest sample and to training sample
#'
#' @param data dataframe or matrix of features and outcomes to be split honestly
#' @param honesty.fraction share of sample to belong to honesty sample
#' @param orf logical, if honest split should be done for orf or not
#'
#' @return named list of honest and training sample
#'
#' @keywords internal
#'
honest_split <- function(data, honesty.fraction, orf) {

  # ------------------------------------------------------------------------------------ #

  # get number of observations in total
  n <- nrow(data)

  if (orf == TRUE) {

    # initiate repeat indicator
    rep_idx <- 1

    # repeat until optimal honesty split is found (optimal in a sense that all outcome categories are represented in both samples)
    repeat{

      # randomize indices for train and honest sample (take care of uneven numbers with floor and ceiling)
      ind <- sample(c(rep(0, ceiling((1-honesty.fraction)*n)), rep(1, floor(honesty.fraction*n))))
      # indicator for which observations go into train and honest set
      honesty_i <- which(ind == 1)
      # separate training set
      train <- data[-honesty_i, ]
      # separate honest set
      honest <- data[honesty_i, ]

      # check if in both data sets all outcome categories are represented or if too many tries have been done
      if(all(sort(unique(train[, 1])) == sort(unique(honest[, 1]))) | rep_idx == 10){
        break
      }

      # repeat indicator
      rep_idx <- rep_idx + 1

    }

    # check if the above procedure was succesful
    if (all(sort(unique(train[, 1])) != sort(unique(honest[, 1])))) {
      stop("At least one of the categories of the input matrix Y contains too few observations.
           This prevents an optimal honesty split. Consider recoding your outcome into less categories or set honesty = FALSE.")
    }

  } else {

    # for classical regression forest the above condition is not binding
    # randomize indices for train and honest sample (take care of uneven numbers with floor and ceiling)
    ind <- sample(c(rep(0, ceiling((1-honesty.fraction)*n)), rep(1, floor(honesty.fraction*n))))
    # indicator for which observations go into train and honest set
    honesty_i <- which(ind == 1)
    # separate training set
    train <- data[-honesty_i, ]
    # separate honest set
    honest <- data[honesty_i, ]

  }

  # ------------------------------------------------------------------------------------ #

  # put it into output
  output <- list(train, honest)
  # set names
  names(output) <- c("trainData", "honestData")

  # return output
  return(output)

  # ------------------------------------------------------------------------------------ #

}

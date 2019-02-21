#' rrf
#'
#' Regression Random Forests (adaptation of ranger with added honesty sample
#' splitting option, forest weights and inference rpocedure as in Lechner (2018))
#'
#' @param X matrix of features
#' @param Y vector of outcomes (as.matrix acceptable too)
#' @param ntree scalar, number of trees in a forest, i.e. bootstrap replications
#' @param mtry scalar, number of randomly selected features
#' @param nmin scalar, minimum node size
#' @param honesty logical, if TRUE honest forest is built using 50:50 data split
#' @param inference logical, if TRUE the weight based inference is conducted
#'
#' @import ranger
#'
#' @return object of type rrf
#'
#' @export
rrf <- function(X, Y, ntree, mtry, nmin, honesty, inference) {
  # needed inputs for the function: X - matrix of features
  #                                 Y - vector of outcomes (as.matrix acceptable too)
  #                                 ntree - number of trees in a forest
  #                                 mtry - number of randomly selected features
  #                                 nmin - minimum leaf size
  #                                 honesty - logical, if TRUE honest forest is built using 50:50 data split
  #                                 inference - logical, if TRUE the weight based inference is conducted (honesty has to be TRUE)

  # -------------------------------------------------------------------------------- #

  ## check for plausibility of options first:
  if (honesty == FALSE & inference == TRUE) {

    warning("For conducting inference honesty is required. Honesty has been set to TRUE.")
    # set honesty to TRUE
    honesty <- TRUE

  }

  # -------------------------------------------------------------------------------- #

  ## save the inputs:
  inputs <- list(ntree, mtry, nmin, honesty, inference)
  names(inputs) <- c("ntree", "mtry", "nmin", "honesty", "inference")

  ## save colnames
  # Y
  if (is.null(colnames(Y))) { colnames(Y) <- "Y" } # check if Y has name
  Y_name <- colnames(Y) # save the name of Y
  Y <- as.numeric(Y) # numeric response as only regression is supported (so far)
  # X
  if (is.null(colnames(X))) { colnames(X) <- paste0("X", rep(1:ncol(X))) } # check if X has name
  X_name <- colnames(X) # save the name of X

  ## set needed dataframe and local variables
  dat <- as.data.frame(cbind(Y, X)) # dataframe
  colnames(dat) <- c(Y_name, X_name) # column names
  n <- nrow(dat) # number of observations

  # -------------------------------------------------------------------------------- #

  # chekc for honesty and inference
  if (honesty == TRUE & inference == TRUE) {

    ## do honest forest estimation here using 50:50 data split as in Lechner (2018)
    # devide into 50:50 honesty sets
    split_data <- honest_split(dat)
    train_data <- split_data$trainData # take out training data
    honest_data <- split_data$honestData # take out honest data
    rows_train_data <- as.numeric(rownames(train_data)) # take rownames of train data as numeric
    rows_honest_data <- as.numeric(rownames(honest_data)) # take rownames of train data as numeric

    # -------------------------------------------------------------------------------- #

    # built the forest structure using training set, i.e. place splits, use subsampling with 0.5 rate
    forest         <- ranger(dependent.variable.name = paste(Y_name), data = train_data,
                             num.trees = ntree, mtry = mtry, replace = FALSE, sample.fraction = 0.5,
                             min.node.size = nmin, importance = "none")

    # -------------------------------------------------------------------------------- #

    # get honest weights
    forest_weights <- get_forest_weights(forest, honest_data, train_data)
    honest_weights <- forest_weights[rows_honest_data, ] # take out honest sample honest weights
    train_weights <- forest_weights[rows_train_data, ] # take out train sample honest weights

    ## make honest predictions, i.e. fitted values based on honest sample
    # honest sample predictions
    honest_pred <- as.matrix(honest_weights %*% honest_data[, 1]) # honest weights for honest data
    rownames(honest_pred) <- rows_honest_data
    # train sample predictions
    train_pred <- as.matrix(train_weights %*% honest_data[, 1]) # honest weights for train data
    rownames(train_pred) <- rows_train_data
    # put the prediction together for whole sample and order them as original data
    forest_pred <- rbind(honest_pred, train_pred)
    # sort according to rownames
    forest_pred <- as.numeric(forest_pred[order(as.numeric(row.names(forest_pred))), ])

    # compute in-sample honest MSE
    honest_mse <- mean((forest_pred - dat[, 1])^2)

    # -------------------------------------------------------------------------------- #

    ## now do the inference based on the weights
    # honest sample variance
    honest_variance <- as.matrix(get_variance(honest_pred, honest_weights, honest_data[, 1]))
    rownames(honest_variance) <- rows_honest_data
    # train sample variance (also honest_data as predictions are done using honest_data too!)
    train_variance <- as.matrix(get_variance(train_pred, train_weights, honest_data[, 1]))
    rownames(train_variance) <- rows_train_data

    # put it again together and sort according to original data
    forest_variance <- rbind(honest_variance, train_variance)
    # sort according to rownames
    forest_variance <- as.numeric(forest_variance[order(as.numeric(row.names(forest_variance))), ])

    # -------------------------------------------------------------------------------- #

    # save forest information
    forest_info <- list(inputs, train_data, honest_data)
    names(forest_info) <- c("inputs", "trainData", "honestData")

    # pack into output
    output <- list(forest, forest_info, forest_weights, forest_pred, forest_variance, honest_mse)
    names(output) <- c("trainForest", "forestInfo", "honestWeights", "honestPredictions", "honestVariance", "honestMSE")

    # -------------------------------------------------------------------------------- #

  } else if (honesty == TRUE & inference == FALSE) {

    ## do honest forest estimation here using 50:50 data split as in Lechner (2018)
    # devide into 50:50 honesty sets
    split_data <- honest_split(dat)
    train_data <- split_data$trainData # take out training data
    honest_data <- split_data$honestData # take out honest data

    # -------------------------------------------------------------------------------- #

    # built the forest structure using training set, i.e. place splits
    forest         <- ranger(dependent.variable.name = paste(Y_name), data = train_data,
                             num.trees = ntree, mtry = mtry, replace = FALSE, sample.fraction = 0.5,
                             min.node.size = nmin, importance = "none")

    # -------------------------------------------------------------------------------- #

    # compute honest predictions
    honest_pred <- get_honest(forest, honest_data, train_data)

    # compute honest MSE
    honest_mse <- mean((dat[, 1] - honest_pred)^2)

    # -------------------------------------------------------------------------------- #

    # save forest information
    forest_info <- list(inputs, train_data, honest_data)
    names(forest_info) <- c("inputs", "trainData", "honestData")

    # pack into output
    output <- list(forest, forest_info, honest_pred, honest_mse)
    names(output) <- c("trainForest", "forestInfo", "honestPredictions", "honestMSE")

    # -------------------------------------------------------------------------------- #

  } else {

    # no honest splitting, i.e. use all data
    train_data <- dat
    honest_data <- NULL

    # estimate standard random forest as in ranger default
    forest <- ranger(dependent.variable.name = paste(Y_name), data = train_data,
                     num.trees = ntree, mtry = mtry, replace = TRUE,
                     min.node.size = nmin, importance = "none")

    # take OOB predictions based on whole sample
    oob_pred <- forest$predictions

    # compute OOB MSE based on whole sample
    oob_mse <- mean((oob_pred - train_data[, 1])^2)

    # save forest information
    forest_info <- list(inputs, train_data, honest_data)
    names(forest_info) <- c("inputs", "trainData", "honestData")

    # -------------------------------------------------------------------------------- #

    # define output of the function
    output <- list(forest, forest_info, oob_pred, oob_mse)
    names(output) <- c("trainForest",  "forestInfo", "oobPredictions", "oobMSE")

    # -------------------------------------------------------------------------------- #

  }

  # return the output of the function
  return(output)

  # -------------------------------------------------------------------------------- #

}



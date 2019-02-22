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

  ## Set the name for the class
  class(output) <- "rrf"

  # -------------------------------------------------------------------------------- #
  # return the output of the function
  return(output)

  # -------------------------------------------------------------------------------- #

}


#' predict.rrf
#'
#' Prediction for new observations based on estimated forest of type \code{rrf}
#'
#' @param object estimated forest object of type \code{rrf}
#' @param new_data matrix X containing the observations to predict
#' @param ... further arguments (currently ignored)
#'
#' @import ranger
#'
#' @return object of class \code{rrf.prediction} with elements
##'   \tabular{ll}{
##'       \code{forestInfo}    \tab info containing forest inputs and data used \cr
##'       \code{forestPredictions} \tab predicted values \cr
##'       \code{forestVariances} \tab variances of predicted values (only if \code{inference=TRUE} in the passed \code{rrf object}) \cr
##'   }
#'
#' @export
predict.rrf <- function(object, new_data, ...) {
  # needed inputs for the function: forest - forest object coming from random_forest function
  #                                 newdata - matrix X containing the observations to predict

  # -------------------------------------------------------------------------------- #

  ## get forest as na object
  forest <- object
  ## save forest inputs
  inputs <- forest$forestInfo$inputs
  honesty <- inputs$honesty
  inference <- inputs$inference
  honest_data <- forest$forestInfo$honestData
  train_data <- forest$forestInfo$trainData

  # take out ranger object
  forest <- forest$trainForest

  ## get train data names (only X)
  train_data_name <- colnames(train_data)[2:ncol(train_data)]

  ## get X matrix as dataframe and check colnames
  # X
  if (is.null(colnames(new_data))) { colnames(new_data) <- paste0("X", rep(1:ncol(new_data))) } # check if X has name
  new_data_name <- colnames(new_data) # save the name of X
  new_data <- as.data.frame(new_data) # as dataframe

  # -------------------------------------------------------------------------------- #

  # check if its compatible with the data used for training
  if (all(train_data_name != new_data_name) | (ncol(new_data) != ncol(train_data)-1)) {

    stop("New data are not compatible with the training data. Check supplied data. Program terminated.")

  }

  # -------------------------------------------------------------------------------- #

  # check if honest forest was estimated and predict accordingly
  if (honesty == TRUE & inference == TRUE) {

    ## run new Xs through estimated train forest and compute predictions based on honest sample
    # predict weights by using forest train, honest data and new_data
    forest_weights_pred <- predict_forest_weights(forest, honest_data, new_data)

    # get predictions by matrix multiplication of weights with honest responses
    forest_pred <- as.numeric(forest_weights_pred%*%honest_data[, 1])
    # -------------------------------------------------------------------------------- #

    ## now do the inference based on the weights

    forest_pred_variance <- get_variance(forest_pred, forest_weights_pred, honest_data[, 1])

    # -------------------------------------------------------------------------------- #

    # save forest information
    forest_info <- list(inputs, new_data)
    names(forest_info) <- c("inputs", "newData")

    # define output of the function
    output <- list(forest_info, forest_pred, forest_pred_variance)
    names(output) <- c("forestInfo", "forestPredictions", "forestVariances")

    # -------------------------------------------------------------------------------- #

  } else if (honesty == TRUE & inference == FALSE) {

    # -------------------------------------------------------------------------------- #

    ## run new Xs through estimated train forest and compute predictions based on honest sample
    # no need to predict weights, get predictions directly through leaves
    forest_pred <- predict_honest(forest, honest_data, new_data)

    # -------------------------------------------------------------------------------- #

    # save forest information
    forest_info <- list(inputs, new_data)
    names(forest_info) <- c("inputs", "newData")

    # define output of the function
    output <- list(forest_info, forest_pred)
    names(output) <- c("forestInfo", "forestPredictions")

    # -------------------------------------------------------------------------------- #

  } else {

    ## no honest splitting, i.e. use all data

    # predict standard random forest as in ranger default
    forest_pred <- predict(forest, data = new_data)$predictions

    # -------------------------------------------------------------------------------- #

    # save forest information
    forest_info <- list(inputs, new_data)
    names(forest_info) <- c("inputs", "newData")

    # define output of the function
    output <- list(forest_info, forest_pred)
    names(output) <- c("forestInfo", "forestPredictions")

    # -------------------------------------------------------------------------------- #

  }

  ## Set the name for the class
  class(output) <- "rrf.prediction"
  # return the output
  return(output)

  # -------------------------------------------------------------------------------- #

}


#' plot.rrf
#'
#' plot forest object of class \code{rrf}
#'
#' @param x estimated forest object of type \code{rrf}
#' @param ... further arguments (currently ignored)
#'
#' @import ggplot2
#'
#' @export
plot.rrf <- function(x, ...) {

  # needed inputs for the function: forest - forest object coming from rrf function

  # -------------------------------------------------------------------------------- #

  ## get forest as x
  forest <- x
  ## save forest inputs
  inputs <- forest$forestInfo$inputs
  honesty <- inputs$honesty
  inference <- inputs$inference
  honest_data <- forest$forestInfo$honestData
  train_data <- forest$forestInfo$trainData

  # -------------------------------------------------------------------------------- #

  # check if honest forest was estimated and predict accordingly
  if (honesty == TRUE & inference == TRUE) {

    # -------------------------------------------------------------------------------- #

    # put data back together
    all_data <- rbind(honest_data, train_data)
    # combine the two into a complete dataset (first honest rownames, then train rownames)
    rownames(all_data) <- c(rownames(honest_data), rownames(train_data)) # put original rownames in
    #forest_fitted_values <- rbind(train_fitted_values, honest_fitted_values) # put the sample together
    all_data <- all_data[order(as.numeric(row.names(all_data))), ] #

    # -------------------------------------------------------------------------------- #

    # take out predictions
    forest_pred <- forest$honestPredictions
    # take out standard deviations
    forest_sd <- sqrt(forest$honestVariance)
    # take out outcomes
    forest_outcomes <- as.numeric(all_data[, 1])
    # put it into dataframe
    df_plot <- as.data.frame(cbind(forest_outcomes, forest_pred, forest_sd))
    colnames(df_plot) <- c("Observed", "Predicted", "StdDev")

    # -------------------------------------------------------------------------------- #

    ## plot true against predicted values
    forest_plot <- ggplot(df_plot, aes(x = forest_outcomes, y = forest_pred))
    forest_plot + geom_errorbar(aes(ymin=forest_pred-forest_sd, ymax=forest_pred+forest_sd), width=.1, color="grey") +
      geom_point(color="black") +
      geom_abline(intercept=0, slope=1, linetype=2, color="red") +
      xlab("Observed Y") +
      ylab("Predicted Y") +
      ggtitle("Honest RF Predictions with CI")+
      coord_fixed()+
      expand_limits(x = min(df_plot[, 1]), y = min(df_plot[, 2]))+
      theme_bw()+
      theme(plot.title = element_text(hjust = 0.5))

    # -------------------------------------------------------------------------------- #

  } else if (honesty == TRUE & inference == FALSE) {

    # -------------------------------------------------------------------------------- #

    # put data back together
    all_data <- rbind(honest_data, train_data)
    # combine the two into a complete dataset (first honest rownames, then train rownames)
    rownames(all_data) <- c(rownames(honest_data), rownames(train_data)) # put original rownames in
    #forest_fitted_values <- rbind(train_fitted_values, honest_fitted_values) # put the sample together
    all_data <- all_data[order(as.numeric(row.names(all_data))), ] #

    # -------------------------------------------------------------------------------- #

    # take out predictions
    forest_pred <- forest$honestPredictions
    # take out outcomes
    forest_outcomes <- as.numeric(all_data[, 1])
    # put it into dataframe
    df_plot <- as.data.frame(cbind(forest_outcomes, forest_pred))
    #colnames(df_plot) <- c("Observed", "Predicted")

    # -------------------------------------------------------------------------------- #

    ## plot true against predicted values
    forest_plot <- ggplot(df_plot, aes(x = forest_outcomes, y = forest_pred))
    forest_plot + geom_point(color="black") +
      geom_abline(intercept=0, slope=1, linetype=2, color="red") +
      xlab("Observed Y") +
      ylab("Predicted Y") +
      ggtitle("Honest RF Predictions")+
      coord_fixed()+
      expand_limits(x = min(df_plot[, 1]), y = min(df_plot[, 2]))+
      theme_bw()+
      theme(plot.title = element_text(hjust = 0.5))

    # -------------------------------------------------------------------------------- #

  } else {

    # -------------------------------------------------------------------------------- #

    ## no honest splitting, i.e. use all data
    # take out predictions
    forest_pred <- forest$oobPredictions
    # take out outcomes
    forest_outcomes <- as.numeric(train_data[, 1])
    # put it into dataframe
    df_plot <- as.data.frame(cbind(forest_outcomes, forest_pred))
    #colnames(df_plot) <- c("Observed", "Predicted")

    # -------------------------------------------------------------------------------- #

    ## plot true against predicted values
    forest_plot <- ggplot(df_plot, aes(x = forest_outcomes, y = forest_pred))
    forest_plot + geom_point(color="black") +
      geom_abline(intercept=0, slope=1, linetype=2, color="red") +
      xlab("Observed Y") +
      ylab("Predicted Y") +
      ggtitle("RF Predictions")+
      coord_fixed()+
      expand_limits(x = min(df_plot[, 1]), y = min(df_plot[, 2]))+
      theme_bw()+
      theme(plot.title = element_text(hjust = 0.5))

    # -------------------------------------------------------------------------------- #

  }

  # no output to return for plot

}


#' summary.rrf
#'
#' summary of a forest object of class \code{rrf}
#'
#' @param object estimated forest object of type \code{rrf}
#' @param latex logical, TRUE if latex summary should be generated
#' @param ... further arguments (currently ignored)
#'
#' @importFrom xtable xtable
#'
#' @export
summary.rrf <- function(object, latex, ...) {

  # needed inputs for the function: forest - forest object coming from random_forest function
  #                                        - latex : logical if the output should be printed in latex code

  # -------------------------------------------------------------------------------- #

  ## get forest as object
  forest <- object
  ## save forest inputs
  inputs <- forest$forestInfo$inputs
  honesty <- inputs$honesty
  inference <- inputs$inference
  mtry <- inputs$mtry
  ntree <- inputs$ntree
  nmin <- inputs$nmin
  honest_data <- forest$forestInfo$honestData
  train_data <- forest$forestInfo$trainData
  type <- "Regression"

  # -------------------------------------------------------------------------------- #

  # check if honest forest was estimated and predict accordingly
  if (honesty == TRUE & inference == TRUE) {

    # -------------------------------------------------------------------------------- #

    ## honest splitting, i.e. use honest data
    # take out summary statistics
    mse <- round(forest$honestMSE, 5)
    trainsize <- nrow(train_data)
    honestsize <- nrow(honest_data)
    features <- ncol(train_data)-1 # take out the response
    # check if subsampling or bootstrapping was used
    if (forest$trainForest$replace == TRUE) { build <- "Bootstrap" } else { build <- "Subsampling" }

    # -------------------------------------------------------------------------------- #

    # structure summary into a list
    output <- list(type, build, ntree, mtry, nmin, honesty, inference, trainsize, honestsize, features, mse)
    names(output) <- c("type", "build", "ntree", "mtry", "nmin", "honesty", "inference", "trainsize", "honestsize", "features", "mse")

    # output matrix
    output_matrix <- matrix(NA, ncol = 1, nrow = length(output))
    # populate output matrix
    rownames(output_matrix) <- names(output) # rownames are names
    colnames(output_matrix) <- "" # no colname
    output_matrix[, 1] <- unlist(output) # column 2 are values

    # generate latex output if selected
    if (latex == TRUE) { colnames(output_matrix) <- "Attributes"
    output_matrix <- xtable(output_matrix, caption = "Random Forest Summary")
    }

    # -------------------------------------------------------------------------------- #

    # -------------------------------------------------------------------------------- #

  } else if (honesty == TRUE & inference == FALSE) {

    # -------------------------------------------------------------------------------- #

    ## honest splitting, i.e. use honest data
    # take out summary statistics
    mse <- round(forest$honestMSE, 5)
    trainsize <- nrow(train_data)
    honestsize <- nrow(honest_data)
    features <- ncol(train_data)-1 # take out the response
    # check if subsampling or bootstrapping was used
    if (forest$trainForest$replace == TRUE) { build <- "Bootstrap" } else { build <- "Subsampling" }

    # -------------------------------------------------------------------------------- #

    # structure summary into a list
    output <- list(type, build, ntree, mtry, nmin, honesty, inference, trainsize, honestsize, features, mse)
    names(output) <- c("type", "build", "ntree", "mtry", "nmin", "honesty", "inference", "trainsize", "honestsize", "features", "mse")

    # output matrix
    output_matrix <- matrix(NA, ncol = 1, nrow = length(output))
    # populate output matrix
    rownames(output_matrix) <- names(output) # rownames are names
    colnames(output_matrix) <- "" # no colname
    output_matrix[, 1] <- unlist(output) # column 2 are values

    # generate latex output if selected
    if (latex == TRUE) { colnames(output_matrix) <- "Attributes"
    output_matrix <- xtable(output_matrix, caption = "Random Forest Summary")
    }

    # -------------------------------------------------------------------------------- #

  } else {

    # -------------------------------------------------------------------------------- #

    ## no honest splitting, i.e. use all data
    # take out summary statistics
    mse <- round(forest$oobMSE, 5)
    trainsize <- nrow(train_data)
    honestsize <- 0
    features <- ncol(train_data)-1 # take out the response
    # check if subsampling or bootstrapping was used
    if (forest$trainForest$replace == TRUE) { build <- "Bootstrap" } else { build <- "Subsampling" }

    # -------------------------------------------------------------------------------- #

    # structure summary into a list
    output <- list(type, build, ntree, mtry, nmin, honesty, inference, trainsize, honestsize, features, mse)
    names(output) <- c("type", "build", "ntree", "mtry", "nmin", "honesty", "inference", "trainsize", "honestsize", "features", "mse")

    # output matrix
    output_matrix <- matrix(NA, ncol = 1, nrow = length(output))
    # populate output matrix
    rownames(output_matrix) <- names(output) # rownames are names
    colnames(output_matrix) <- "" # no colname
    output_matrix[, 1] <- unlist(output) # column 2 are values

    # generate latex output if selected
    if (latex == TRUE) { colnames(output_matrix) <- "Attributes"
    output_matrix <- xtable(output_matrix, caption = "Random Forest Summary")
    }

    # -------------------------------------------------------------------------------- #

  }

  # return output
  return(output_matrix)

}

#' check X data
#'
#' Checks the input data as a numeric matrix
#'
#' @param X matrix of input features X
#'
check_X <- function(X) {

  if (inherits(X, "matrix") & !is.numeric(X)) {
    stop("The input matrix X must be coded as a numeric matrix. In case of factor variables, encode these as integer values.")
  }

  if (inherits(X, "dgCMatrix")) {
    stop("Currently no sparse matrices of class 'dgCMatrix' are supported. Recode the input data using as.matrix() command.")
  }

  if (any(is.na(X))){
    stop("The input matrix X contains at least one NA value. Check the input data and remove NAs.")
  }

}

#' check Y data
#'
#' Checks the input data as a numeric matrix/vector
#'
#' @param Y matrix of input outcomes Y
#' @param X matrix of input features X
#'
#' @return Y
#'
check_Y <- function(Y, X) {

  if ((inherits(Y, "matrix") | inherits(Y, "numeric")) & !is.numeric(Y)) {
    stop("The input matrix Y must be coded as a numeric matrix or a numeric vector. In case of factor variables, encode these as integer values.")
  }

  if (inherits(Y, "dgCMatrix")) {
    stop("Currently no sparse matrices of class 'dgCMatrix' are supported. Recode the input data using as.matrix() command.")
  }

  if (any(is.na(Y))){
    stop("The input matrix Y contains at least one NA value. Check the input data and remove NAs.")
  }

  if (length(Y) != nrow(X)) {
    stop("The number of outcomes does not equal the number of features. Check your data.")
  }

  if (inherits(Y, "numeric") & is.numeric(Y)) {

    Y <- as.matrix(Y)

  }

  Y

}

#' check if Y is discrete
#'
#' Checks the input data as discrete outcome
#'
#' @param Y matrix of input outcomes Y
#'
#' @return Y
#'
check_discrete_Y <- function(Y) {

  if (all((Y %% 1) != 0)) {
    stop("The input matrix Y must be coded without containing decimals. Recode the inputs. Example: Y = {1,2,3}.")
  }

  if (length(unique(Y)) > 10) {
    warning("The input matrix Y contains more than 10 distinct values. This might be not optimal for an Ordered Choice Model.
            Consider recoding your outcome into less categories.")
  }

  if (any(table(Y)/nrow(Y) < 0.05)) {
    warning("At least one of the categories of the input matrix Y contains less than 5% of observations.
            This might be not optimal for an Ordered/Multinomial Choice Model. Consider recoding your outcome into less categories.")
  }

  if (!all(sort(unique(Y)) == seq_along(unique(Y)))) {
    warning(paste(c("The input matrix Y has been recoded to: ", seq_along(unique(Y))), sep = " ", collapse = " "))
    # recode Y
    for(i in seq_along(unique(Y))) {Y <- replace(Y, Y == sort(unique(Y))[i], seq_along(unique(Y))[i]) }
  }

  Y

}

#' check num.trees
#'
#' Checks the input data of num.trees
#'
#' @param num.trees scalar, number of trees to be estimated
#'
#' @return num.trees
#'
check_num_trees <- function(num.trees) {

  if (is.null(num.trees)) {

    num.trees <- 1000

  } else if (!is.numeric(num.trees) | num.trees <= 0 ) {

    stop("Error: Invalid value for num.trees. num.trees must be a positive number.")
  }

  num.trees

}

#' check mtry
#'
#' Checks the input data of mtry
#'
#' @param mtry scalar, number of randomly selected features
#' @param X matrix of input features X
#'
#' @return mtry
#'
check_mtry <- function(mtry, X) {

  if (is.null(mtry)) {

    mtry <- ceiling(sqrt(ncol(X)))

  } else if (!is.numeric(mtry) | mtry <= 0 | mtry > ncol(X)) {

    stop("Error: Invalid value for mtry. mtry must be a positive number, lower than or equal the maximum number of features.")
  }

  mtry

}

#' check min.node.size
#'
#' Checks the input data of min.node.size
#'
#' @param min.node.size scalar, minimum node size
#' @param X matrix of input features X
#'
#' @return min.node.size
#'
check_min_node_size <- function(min.node.size, X) {

if (is.null(min.node.size)) {

  min.node.size <- 5

} else if (!is.numeric(min.node.size) | min.node.size <= 0 | min.node.size > nrow(X)) {

  stop("Error: Invalid value for min.node.size. min.node.size must be a positive number, lower than number of observations.")
}

  min.node.size

}

#' check replace
#'
#' Checks the input for replace
#'
#' @param replace logical, if TRUE bootstrapping, if FALSE subsampling
#'
#' @return replace
check_replace <- function(replace) {

  if (!(is.logical(replace))) {

    replace <- FALSE
    warning("replace must be logical. replace has been set to FALSE as a default.")

  }

  replace

}

#' check sample.fraction
#'
#' Checks the input data of sample.fraction
#'
#' @param sample.fraction scalar, fraction of data used for subsampling
#'
#' @return sample.fraction
#'
check_sample_fraction <- function(sample.fraction) {

  if (is.null(sample.fraction)) {

    sample.fraction <- 0.5

  } else if (!is.numeric(sample.fraction) | sample.fraction <= 0 | sample.fraction > 1) {

    stop("Error: Invalid input value for sample.fraction. The number must be within the interval of [0,1).")
  }

  sample.fraction

}

#' check honesty
#'
#' Checks the input for honesty
#'
#' @param honesty logical, if TRUE honest forest is built using 50:50 data split
#'
#' @return honesty
check_honesty <- function(honesty) {

  if (!(is.logical(honesty))) {

    honesty <- TRUE
    warning("honesty must be logical. honesty has been set to TRUE as a default.")

  }

  honesty

}

#' check inference
#'
#' Checks the input for inference
#'
#' @param inference ogical, if TRUE the weight based inference is conducted
#'
#' @return inference
check_inference <- function(inference) {

  if (!(is.logical(inference))) {

    inference <- FALSE
    warning("inference must be logical. inference has been set to FALSE as a default.")

  }

  inference

}

#' check newdata
#'
#' Checks the input for newdata for predict.orf/predict.mrf
#'
#' @param new_data matrix X containing the observations to predict
#' @param X matrix of input features X
#'
check_newdata <- function(new_data, X) {

  if (ncol(new_data) != ncol(X)) {

    stop("new_data must have the same column dimension as the training data. Check your data.")

  }

  if (colnames(new_data) != colnames(X)) {

    stop("Column names of new_data differ from the training data. Check your data.")

  }

  check_X(new_data)

}

#' check latex
#'
#' Checks the input for latex
#'
#' @param latex logical, TRUE if latex summary should be generated
#'
#' @return latex
check_latex <- function(latex) {

  if (!(is.logical(latex))) {

    latex <- FALSE
    warning("latex must be logical. latex has been set to FALSE as a default.")

  }

  latex

}

#' check window size
#'
#' Checks the input data of window in margins
#'
#' @param window scalar, share of SD of X used for margins
#'
#' @return window
#'
check_window <- function(window) {

  if (is.null(window)) {

    window <- 0.1

  } else if (!is.numeric(window) | window <= 0 | window > 1) {

    stop("Error: Invalid value for window. window must be within [0,1) interval.")
  }

  window

}

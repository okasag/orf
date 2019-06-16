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

#' check nmin
#'
#' Checks the input data of mtry
#'
#' @param nmin scalar, minimum node size
#' @param X matrix of input features X
#'
#' @return nmin
#'
check_nmin <- function(nmin, X) {

if (is.null(nmin)) {

  nmin <- 5

} else if (!is.numeric(nmin) | nmin <= 0 | nmin > nrow(X)) {

  stop("Error: Invalid value for nmin. nmin must be a positive number, lower than number of observations.")
}

nmin

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


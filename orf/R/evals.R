#' Ranked Probability Score
#'
#' Computes the mean ranked probability score (RPS) for evaluating the accuracy of
#' ordered probability predictions
#'
#' @param predictions matrix of predictions (n x categories)
#' @param observed vector of observed ordered categorical outcomes (n x 1)
#'
#' @return scalar, mean RPS for given predictions
#'
#' @keywords internal
#'
rps <- function(predictions, observed){

  # ------------------------------------------------------------------------------------ #

  # get parameteres
  ncat <- as.numeric(ncol(predictions)) # number of categories
  npred <- as.numeric(nrow(predictions)) # number of observations

  # create probability distribution for observed outcomes
  observed_dist <- matrix(0, nrow = npred, ncol = ncat)
  # populate it
  for (i in 1:npred) {
    observed_dist[i, observed[i]] <- 1
  }

  # ------------------------------------------------------------------------------------ #

  # prepare 0 vectors for rps and cum
  rps <- numeric(npred)
  cum <- numeric(npred)

  # loop over the categories (inspired by and thanks to: https://opisthokonta.net/?p=1333)
  for (i in 1:ncat){

    cum <- cum + (rowSums(matrix(predictions[, 1:i], ncol = i)) - rowSums(matrix(observed_dist[, 1:i], ncol = i)))^2

  }

  # compute rps for each observation (scale it accordingly)
  rps <- (1/(ncat - 1))*cum

  # take mean of rps
  mrps <- mean(rps)

  # ------------------------------------------------------------------------------------ #

  # return the mrps
  return(mrps)

  # ------------------------------------------------------------------------------------ #

}


#' Mean Squared Error
#'
#' computes the mean squared error (MSE) for evaluating the accuracy of
#' ordered/unordered probability predictions
#'
#' @param predictions matrix of predictions (n x categories)
#' @param observed vector of observed ordered categorical outcomes (n x 1)
#'
#' @return scalar, sum MSE for given predictions
#'
#' @keywords internal
#'
mse <- function(predictions, observed){

  # ------------------------------------------------------------------------------------ #

  # get parameteres
  ncat <- as.numeric(ncol(predictions)) # number of categories
  npred <- as.numeric(nrow(predictions)) # number of observations

  # create probability distribution for observed outcomes
  observed_dist <- matrix(0, nrow = npred, ncol = ncat)
  # populate it
  for (i in 1:npred) {
    observed_dist[i, observed[i]] <- 1
  }

  # use apply to calculate mse (multi brier)
  mse <- mean(rowSums(apply((observed_dist - predictions), 2, function(x) {x^2} )))

  # ------------------------------------------------------------------------------------ #

  # return the mse
  return(mse)

  # ------------------------------------------------------------------------------------ #

}

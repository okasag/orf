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

  #               predictions <- matrix of predictions (n x categories)
  #               observed <- vector of observed ordered categorical outcomes (n x 1)

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

  # prepare 0 vector for rps
  rps <- numeric(npred)

  # compute the rps
  for (rr in 1:npred){

    cumulative <- 0
    for (i in 1:ncat){

      cumulative <- cumulative + (sum(predictions[rr, 1:i]) - sum(observed_dist[rr, 1:i]))^2

    }
    rps[rr] <- (1/(ncat - 1))*cumulative
  }

  # take mean of rps
  mrps <- mean(rps)

  # ------------------------------------------------------------------------------------ #

  # return the mrps
  return(mrps)

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

  #               predictions <- matrix of predictions (n x categories)
  #               observed <- vector of observed ordered categorical outcomes (n x 1)

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

}

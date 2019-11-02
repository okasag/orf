# -----------------------------------------------------------------------------
# This file is part of orf.
#
# orf is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# orf is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with orf. If not, see <http://www.gnu.org/licenses/>.
#
# Written by:
#
# Gabriel Okasa
# Swiss Institute for Empirical Economic Research
# University of St.Gallen
# Varnbüelstrasse 14
# 9000 St.Gallen
# Switzerland
# -----------------------------------------------------------------------------

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

  # prepare 0 vector for rps
  rps <- numeric(npred)

  # compute the rps (thanks to: https://opisthokonta.net/?p=1333)
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

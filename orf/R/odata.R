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

#' Simulated Example Dataset
#'
#' @description
#' A simulated example dataset with ordered categorical outcome variable
#' containing different types of covariates for illustration purposes.
#'
#' @format A data frame with 1000 rows and 5 variables
#'
#' @return
#'   \item{Y}{ordered outcome, classes 1, 2, and 3}
#'   \item{X1}{continuous covariate, N(0,1)}
#'   \item{X2}{categorical covariate, values 1, 2, and 3}
#'   \item{X3}{binary covariate, values 0 and 1}
#'   \item{X4}{continuous covariate, N(0,10)}
#'
#' @details
#' For the exact data generating process, see the example below.
#'
#' @examples
#' \dontrun{
#' # generate example data
#'
#' # set seed for replicability
#' set.seed(123)
#'
#' # number of observations
#' n  <- 1000
#'
#' # various covariates
#' X1 <- rnorm(n, 0, 1)    # continuous
#' X2 <- rbinom(n, 2, 0.5) # categorical
#' X3 <- rbinom(n, 1, 0.5) # dummy
#' X4 <- rnorm(n, 0, 10)   # noise
#'
#' # bind into matrix
#' X <- as.matrix(cbind(X1, X2, X3, X4))
#'
#' # deterministic component
#' deterministic <- X1 + X2 + X3
#' # generate continuous outcome with logistic error
#' Y <- deterministic + rlogis(n, 0, 1)
#' # thresholds for continuous outcome
#' cuts <- quantile(Y, c(0, 1/3, 2/3, 1))
#' # discretize outcome into ordered classes 1, 2, 3
#' Y <- as.numeric(cut(Y, breaks = cuts, include.lowest = TRUE))
#'
#' # save data as a dataframe
#' odata <- as.data.frame(cbind(Y, X))
#'
#' # end of data generating
#' }
#'
"odata"

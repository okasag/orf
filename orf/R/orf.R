#' orf
#'
#' An Implementation of the Ordered Random Forest Algorithm
#' as in Lechner & Okasa (2019) and other related estimators for
#' discrete choice models based on the random forest algorithm.
#' These include models with ordered, multinomial as well as binary
#' response. Standard random forest estimator for continuous response
#' is implemented, too. All the forest based algorithms rely on the
#' fast C++ forest implementation from the ranger package. Additionally
#' to common implementations the orf package provides functions for
#' estimating forest weights as well as marginal effects and thus
#' provides similar output as in standard econometric models for
#' ordered choice.
#'
#' @docType package
#'
#' @author Gabriel Okasa \email{gabriel.okasa@@unisg.ch}
#'
#' @name orf
#'
NULL

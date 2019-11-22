#' orf: Ordered Random Forests
#'
#' An implementation of the Ordered Forest estimator as developed
#' in Lechner & Okasa (2019). The Ordered Forest flexibly
#' estimates the conditional probabilities of models with ordered
#' categorical outcomes (so-called ordered choice models).
#' Additionally to common machine learning algorithms the \code{orf}
#' package provides functions for estimating marginal effects as well
#' as statistical inference thereof and thus provides similar output
#' as in standard econometric models for ordered choice. The core
#' forest algorithm relies on the fast C++ forest implementation
#' from the \code{ranger} package (Wright & Ziegler, 2017).
#'
#' @docType package
#'
#' @author Gabriel Okasa, Michael Lechner
#'
#' @examples
#' ## Ordered Forest
#' require(orf)
#'
#' # load example data
#' data(odata)
#'
#' # specify response and covariates
#' Y <- as.numeric(odata[, 1])
#' X <- as.matrix(odata[, -1])
#'
#' # estimate Ordered Forest with default settings
#' orf_fit <- orf(X, Y)
#'
#' # print output of the orf estimation
#' print(orf_fit)
#'
#' # show summary of the orf estimation
#' summary(orf_fit)
#'
#' # plot the estimated probability distributions
#' plot(orf_fit)
#'
#' # predict with the estimated orf
#' predict(orf_fit)
#' \donttest{
#' # estimate marginal effects of the orf
#' margins(orf_fit)
#' }
#'
#' @references
#' \itemize{
#'   \item Lechner, M., & Okasa, G. (2019). Random Forest Estimation of the Ordered Choice Model. arXiv preprint arXiv:1907.02436. \url{https://arxiv.org/abs/1907.02436}
#'   \item Goller, D., Knaus, M. C., Lechner, M., & Okasa, G. (2018). Predicting Match Outcomes in Football by an Ordered Forest Estimator (No. 1811). University of St. Gallen, School of Economics and Political Science. \url{http://ux-tauri.unisg.ch/RePEc/usg/econwp/EWP-1811.pdf}
#'   \item Wright, M. N. & Ziegler, A. (2017). ranger: A fast implementation of random forests for high dimensional data in C++ and R. J Stat Softw 77:1-17. \url{https://doi.org/10.18637/jss.v077.i01}.
#' }
#'
#' @name orf-package
#'
#' @useDynLib orf, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

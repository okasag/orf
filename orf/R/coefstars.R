#' Significance Stars for Estimated Coefficients
#'
#' function for creating significance stars for estimated effects which
#' can be passed into a xtable for printing latex output
#'
#' @param coefs matrix of estimated coefficients or fitted values/predictions
#' @param pvalues matrix of pvalues for the corresponding coefficients
#'
#' @return a character matrix of coefficients with significance stars
#'
#' @examples
#' coefs <- matrix(data = c( 0.25, 0.5, 0.75, 1), ncol = 2, nrow = 2)
#' pvalues <- matrix(data = c( 0.1, 0.04, 0.009, 0.00001), ncol = 2, nrow = 2)
#' coefstars(coefs, pvalues)
#'
#' @export
#'
coefstars <- function(coefs, pvalues) {
#        coefs <- vector of estimated coefficients or fitted values/predictions
#        pvalues <- vector of pvalues for the corresponding coefficients


  # check if coefs and pvalues have the same dimensions
  if (all(dim(coefs))!=all(dim(pvalues))) {
    stop("Dimensions of coefficients and p-values do not match.
         Programme temrinated.")
  }

  # chekc rownames of coefs
  if (is.null(rownames(coefs))) {
    coefs_rownames <- paste0("X", rep(1:nrow(coefs)))
  } else {
    coefs_rownames <- rownames(coefs)
  }

  # chekc colnames of coefs
  if (is.null(rownames(coefs))) {
    coefs_colnames <- paste0("Cat", rep(1:ncol(coefs)))
  } else {
    coefs_colnames <- colnames(coefs)
  }

  # generate stars (thanks to:
  # http://myowelt.blogspot.com/2008/04/beautiful-correlation-tables-in-r.html)
  stars <- ifelse(pvalues < .001, "***", ifelse(pvalues < .01, "** ",
                  ifelse(pvalues < .05, "*  ", "   ")))

  # trunctuate the coefs to three decimals (take care of minus sign as well)
  coefs <- apply(coefs, 2, function(x) format(sprintf("% .3f", round(x, 3))))

  # paste effest with stars together
  star_coefs <- sapply(seq_along(coefs[1, ]), function(i) {

    as.matrix(paste0(coefs[, i], stars[, i]))

    })
  # add colnames and rownames
  colnames(star_coefs) <- coefs_colnames
  rownames(star_coefs) <- coefs_rownames

  # print out table
  output <- star_coefs

  # return  output
  return(output)

}

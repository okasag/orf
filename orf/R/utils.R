#' honest 50/50 sample split
#'
#' Creates honest sample split by randomly selecting 50 percent of observations
#' to belong to honest sample and 50 percent to training sample
#'
#' @param data dataframe or matrix of features and outcomes to be split honestly
#'
#' @return named list of honest and training sample

honest_split <- function(data) {

  # get number of observations in total
  n <- nrow(data)
  # needed inputs for the function: data - dataframe which should be split into 50:50 sets
  ind <- sample(c(rep(0, n/2), rep(1, n/2))) # randomize indicators
  honesty_i <- which(ind == 1) # indicator for whch observations go into train and honest set
  train <- data[-honesty_i, ] # separate training set
  #rownames(train) <- seq(1:nrow(train)) # set rownames
  honest <- data[honesty_i, ] # separate honest set
  #rownames(honest) <- seq(1:nrow(honest)) # set rownames

  # put it into output
  output <- list(train, honest)
  names(output) <- c("trainData", "honestData")

  # return output
  return(output)

}

#' Significance Stars for Estimated Coefficients
#'
#' function for creating significance stars for estimated effects which
#' can be passed into a xtable for printing latex output
#'
#' @param coefs matrix of estimated coefficients or fitted values/predictions
#' @param pvalues matrix of pvalues for the corresponding coefficients
#'
#' @return a character matrix of coefficients with significance stars

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


#' Formatted output for marginal effects with inference
#'
#' function for creating inference table output for estimated effects which
#' can be passed into \code{print.margins.orf}
#'
#' @param x object of type \code{margins.orf}
#'
margins_output <- function(x) {

  output_matrix <- matrix(NA, nrow = 1, ncol = 4)

  cat("ORF Marginal Effects: \n\n")
  cat("---------------------------------------------------------------------------------", "\n")

  for (var_idx in 1:nrow(x$MarginalEffects)) {

    cat(rownames(x$MarginalEffects)[var_idx], "\n")
    cat("                    Cat ", "     Effect", "    StdDev", "    tValue ", "   pValue", "     ", "\n")

    for (cat_idx in 1:ncol(x$MarginalEffects)) {

      # generate stars (thanks to:
      # http://myowelt.blogspot.com/2008/04/beautiful-correlation-tables-in-r.html)
      stars <- ifelse(x$pValues[var_idx, cat_idx] < .001, "***",
                      ifelse(x$pValues[var_idx, cat_idx] < .01, "** ",
                             ifelse(x$pValues[var_idx, cat_idx] < .05, "*  ", "   ")))

      # print estimates for each category iteratively
      output_matrix[1, 1] <- x$MarginalEffects[var_idx, cat_idx]
      output_matrix[1, 2] <- x$StandardErrors[var_idx, cat_idx]
      output_matrix[1, 3] <- x$tValues[var_idx, cat_idx]
      output_matrix[1, 4] <- x$pValues[var_idx, cat_idx]


      cat("                 |  ", cat_idx, "  |  ") # prit out the categories

      cat(format(sprintf("%8.4f", round(output_matrix, 4)), width = 10), stars, "  |  ") # print out the estimates

      cat("\n") # break the line


    }

  }

  cat("---------------------------------------------------------------------------------", "\n")
  cat("Significance levels correspond to: *** .< 0.01, ** .< 0.05, * .< 0.1 \n")
  cat("---------------------------------------------------------------------------------", "\n")

}



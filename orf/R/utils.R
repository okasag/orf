#' honest sample split
#'
#' Creates honest sample split by randomly selecting prespecified share of observations
#' to belong to honest sample and to training sample
#'
#' @param data dataframe or matrix of features and outcomes to be split honestly
#' @param honesty.fraction share of sample to belong to honesty sample
#' @param orf logical, if honest split should be done for orf or not
#'
#' @return named list of honest and training sample

honest_split <- function(data, honesty.fraction, orf) {

  # get number of observations in total
  n <- nrow(data)

  if (orf == TRUE) {

    # initiate repeat indicator
    rep_idx <- 1

    # repeat until optimal honesty split is found (optimal in a sense that all outcome categories are represented in both samples)
    repeat{

      # randomize indices for train and honest sample (take care of uneven numbers with floor and ceiling)
      ind <- sample(c(rep(0, ceiling((1-honesty.fraction)*n)), rep(1, floor(honesty.fraction*n))))
      # indicator for which observations go into train and honest set
      honesty_i <- which(ind == 1)
      # separate training set
      train <- data[-honesty_i, ]
      # separate honest set
      honest <- data[honesty_i, ]

      # check if in both data sets all outcome categories are represented or if too many tries have been done
      if(all(sort(unique(train[, 1])) == sort(unique(honest[, 1]))) | rep_idx == 10){
        break
      }

      # repeat indicator
      rep_idx <- rep_idx + 1

    }

    # check if in the above procedure was succesful
    if (all(sort(unique(train[, 1])) != sort(unique(honest[, 1])))) {
      stop("At least one of the categories of the input matrix Y contains too few observations.
          This prevents an optimal honesty split. Consider recoding your outcome into less categories or set honesty = FALSE.")
    }

  } else {

    # for classical regression forest the above condition is not binding
    # randomize indices for train and honest sample (take care of uneven numbers with floor and ceiling)
    ind <- sample(c(rep(0, ceiling((1-honesty.fraction)*n)), rep(1, floor(honesty.fraction*n))))
    # indicator for which observations go into train and honest set
    honesty_i <- which(ind == 1)
    # separate training set
    train <- data[-honesty_i, ]
    # separate honest set
    honest <- data[honesty_i, ]

  }

  # put it into output
  output <- list(train, honest)
  # set names
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
      stars <- ifelse(x$pValues[var_idx, cat_idx] < .01, "***",
                      ifelse(x$pValues[var_idx, cat_idx] < .05, "** ",
                             ifelse(x$pValues[var_idx, cat_idx] < .1, "*  ", "   ")))

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


#' Formatted latex output for marginal effects with inference
#'
#' function for creating latex inference table output for estimated effects which
#' can be passed into \code{print.margins.orf}
#'
#' @param x object of type \code{margins.orf}
#'
#' @importFrom xtable xtable print.xtable
#'
margins_output_latex <- function(x) {

  # get number of categories
  ncat <- ncol(x$MarginalEffects)
  # get number of variables
  nvar <- nrow(x$MarginalEffects)
  # get variable names
  varnames <- rownames(x$MarginalEffects)

  # create empty output matrix
  output_matrix <- matrix("", nrow = (nvar*ncat), ncol = 7)
  rownames(output_matrix) <- rep("default", (nvar*ncat)) # generate unique identifier

  for (var_idx in 0:(nvar-1)) {


    for (cat_idx in 1:ncat) {

      # generate stars (thanks to:
      # http://myowelt.blogspot.com/2008/04/beautiful-correlation-tables-in-r.html)
      stars <- ifelse(x$pValues[var_idx+1, cat_idx] < .01, "***",
                      ifelse(x$pValues[var_idx+1, cat_idx] < .05, "** ",
                             ifelse(x$pValues[var_idx+1, cat_idx] < .1, "*  ", "   ")))

      # print estimates for each category iteratively
      rownames(output_matrix)[(var_idx*ncat)+cat_idx] <- paste0(varnames[var_idx+1], cat_idx)
      output_matrix[1+(var_idx*ncat), 1] <- varnames[var_idx+1] # fit in variable name
      output_matrix[(var_idx*ncat)+cat_idx, 2] <- cat_idx # fit in category
      output_matrix[(var_idx*ncat)+cat_idx, 3] <- x$MarginalEffects[var_idx+1, cat_idx]
      output_matrix[(var_idx*ncat)+cat_idx, 4] <- x$StandardErrors[var_idx+1, cat_idx]
      output_matrix[(var_idx*ncat)+cat_idx, 5] <- x$tValues[var_idx+1, cat_idx]
      output_matrix[(var_idx*ncat)+cat_idx, 6] <- x$pValues[var_idx+1, cat_idx]
      output_matrix[(var_idx*ncat)+cat_idx, 7] <- stars



    }

  }

  # add colnames
  colnames(output_matrix) <- c("Variable", "Category", "Effect", "Std.Error", "t-Value", "p-Value", " ")
  # define as data.frame
  output_matrix <- as.data.frame(output_matrix)
  # format the output matrix
  for (i in 3:6) {
    output_matrix[, i] <- as.character(output_matrix[, i])
    output_matrix[, i] <- round(as.numeric(output_matrix[, i]), 4)
  }
  # put caption an latex environment
  xoutput <- xtable(output_matrix, digits = 4, caption = "ORF Marginal Effects")
  # put hline after each variable
  print.xtable(xoutput, hline.after = c(0, seq(ncat, ncat*nvar, ncat)), type = "latex", include.rownames = FALSE, comment = FALSE)

}


#' repeat rows of a matrix
#'
#' function for replicating rows of a matrix n number of times
#'
#' @param matrix matrix which rows should be replicated
#' @param n number of times to repeat
#'
rep_row<-function(matrix, n){
  # thanks to: https://www.r-bloggers.com/a-quick-way-to-do-row-repeat-and-col-repeat-rep-row-rep-col/
  matrix(rep(matrix, each = n), nrow = n)
}

### safety edit: 14 June 2019 ###

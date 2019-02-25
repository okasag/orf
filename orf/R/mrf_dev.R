#' Get MRF Variance
#'
#' get variance of multinomial random forest predictions based on honest sample
#' splitting as described in Lechner (2018)
#'
#' @param honest_pred list of vectors of honest forest predictions
#' @param honest_weights list of n x n matrices of honest forest weights
#' @param train_pred list of vectors of honest forets predictions from train sample
#' @param train_weights list of vectors of honest forets predictions from train sample
#' @param Y_ind_honest list of vectors of 0-1 outcomes for the honest sample
#'
#' @return vector of MRF variances
get_mrf_variance <- function(honest_pred, honest_weights, train_pred, train_weights, Y_ind_honest) {

  # needed inputs for the function: honest_pred - list of vectors of honest forest predictions
  #                                 honest_weights - list of n x n matrices of honest forest weights
  #                                 train_pred - list of vectors of honest forets predictions from train sample
  #                                 train_weights - list of vectors of honest forets predictions from train sample
  #                                 Y_ind_honest - list of vectors of 0-1 outcomes for the honest sample

  # ----------------------------------------------------------------------------------- #

  # first get honest and train rownames
  rows_honest_data <- as.numeric(rownames(honest_pred[[1]]))
  rows_train_data <- as.numeric(rownames(train_pred[[1]]))
  # get also categories
  categories <- seq(1:(length(Y_ind_honest)))

  # ----------------------------------------------------------------------------------- #

  ### variance for the orf predictions
  ## compute prerequisities for variance of honest orf predictions
  ## compute the conditional means (predictions): already have this as forest_pred
  # divide it by N to get the "mean"
  honest_pred_mean <- lapply(honest_pred, function(x) x/length(honest_pred[[1]]))
  train_pred_mean <- lapply(train_pred, function(x) x/length(train_pred[[1]]))

  # calculate standard multiplication of weights and outcomes: honest_weights*y_ind_honest
  honest_multi <- mapply(function(x,y) lapply(seq_along(y), function(i) x[i, ] * y), honest_weights, Y_ind_honest, SIMPLIFY = FALSE)
  train_multi <- mapply(function(x,y) lapply(seq_along(y), function(i) x[i, ] * y), train_weights, Y_ind_honest, SIMPLIFY = FALSE)
  # subtract the mean from each obs i

  # subtract the mean from each obs i
  honest_multi_demeaned <- mapply(function(x,y) mapply(function(x,y) x-y, x, y, SIMPLIFY = FALSE), honest_multi, honest_pred_mean, SIMPLIFY = FALSE)
  train_multi_demeaned <- mapply(function(x,y) mapply(function(x,y) x-y, x, y, SIMPLIFY = FALSE), train_multi, train_pred_mean, SIMPLIFY = FALSE)

  # ----------------------------------------------------------------------------------- #

  ## now do the single variances for each category m
  # square the demeaned
  honest_multi_demeaned_sq <- lapply(honest_multi_demeaned, function(x) lapply(x, function(x) x^2))
  train_multi_demeaned_sq <- lapply(train_multi_demeaned, function(x) lapply(x, function(x) x^2))

  # sum all obs i together
  honest_multi_demeaned_sq_sum <- lapply(honest_multi_demeaned_sq, function(x) lapply(x, function(x) sum(x)))
  train_multi_demeaned_sq_sum <- lapply(train_multi_demeaned_sq, function(x) lapply(x, function(x) sum(x)))

  # multiply by N/N-1 (normalize)
  honest_multi_demeaned_sq_sum_norm <- lapply(honest_multi_demeaned_sq_sum, function(x) lapply(x, function(x) x*(length(honest_pred[[1]])/(length(honest_pred[[1]])-1)) ))
  train_multi_demeaned_sq_sum_norm <- lapply(train_multi_demeaned_sq_sum, function(x) lapply(x, function(x) x*(length(train_pred[[1]])/(length(train_pred[[1]])-1)) ))

  # put it into a shorter named object
  honest_variance <- honest_multi_demeaned_sq_sum_norm
  train_variance <- train_multi_demeaned_sq_sum_norm
  ## single variances done (no covariances needed for MRF)

  ## output for final variances of marginal effects
  # coerce to a matrix
  honest_var <- sapply(honest_variance, function(x) sapply(x, function(x) as.matrix(x)))
  train_var <- sapply(train_variance, function(x) sapply(x, function(x) as.matrix(x)))

  ## put it together according to rownames
  rownames(honest_var) <- rows_honest_data # rownames
  rownames(train_var) <- rows_train_data # rownames
  # combine and sort
  forest_variance <- rbind(honest_var, train_var)
  # sort according to rownames
  forest_variance <- forest_variance[order(as.numeric(row.names(forest_variance))), ]

  # add names
  colnames(forest_variance) <- sapply(categories, function(x) paste("Category", x, sep = " "))

  # ----------------------------------------------------------------------------------- #

  ## return the matrix
  output <- forest_variance
  # output
  return(output)

  # ----------------------------------------------------------------------------------- #

}

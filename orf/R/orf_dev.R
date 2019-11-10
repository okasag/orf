#' Get ORF Variance
#'
#' get variance of ordered random forest predictions based on honest sample
#' splitting as described in Lechner (2018)
#'
#' @param honest_pred list of vectors of honest forest predictions
#' @param honest_weights list of n x n matrices of honest forest weights
#' @param train_pred list of vectors of honest forests predictions from train sample
#' @param train_weights list of vectors of honest forests predictions from train sample
#' @param Y_ind_honest list of vectors of 0-1 outcomes for the honest sample
#'
#' @return vector of ORF variances
#'
#' @keywords internal
#'
get_orf_variance <- function(honest_pred, honest_weights, train_pred, train_weights, Y_ind_honest) {

  # ----------------------------------------------------------------------------------- #

  # first get honest and train rownames
  rows_honest_data <- as.numeric(rownames(honest_pred[[1]]))
  rows_train_data <- as.numeric(rownames(train_pred[[1]]))
  # get also categories
  categories <- seq(1:(length(Y_ind_honest)+1))

  # ----------------------------------------------------------------------------------- #

  ## single variances computation
  # compute the conditional means (predictions): already have this as forest_pred
  # divide it by N to get the "mean"
  honest_pred_mean <- lapply(honest_pred, function(x) x/length(honest_pred[[1]]))
  train_pred_mean <- lapply(train_pred, function(x) x/length(train_pred[[1]]))

  # calculate standard multiplication of weights and outcomes: honest_weights*y_ind_honest
  honest_multi <- mapply(function(x,y) lapply(seq_along(y), function(i) x[i, ] * y), honest_weights, Y_ind_honest, SIMPLIFY = FALSE)
  train_multi <- mapply(function(x,y) lapply(seq_along(y), function(i) x[i, ] * y), train_weights, Y_ind_honest, SIMPLIFY = FALSE)

  # subtract the mean from each obs i
  honest_multi_demeaned <- mapply(function(x,y) mapply(function(x,y) x-y, x, y, SIMPLIFY = FALSE), honest_multi, honest_pred_mean, SIMPLIFY = FALSE)
  train_multi_demeaned <- mapply(function(x,y) mapply(function(x,y) x-y, x, y, SIMPLIFY = FALSE), train_multi, train_pred_mean, SIMPLIFY = FALSE)

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

  # ----------------------------------------------------------------------------------- #

  ## covariances computation
  # multiply forest_var_multi_demeaned according to formula for covariance (shifted categories needed for computational convenience)
  # honest sample
  honest_multi_demeaned_0_last <- append(honest_multi_demeaned, list(rep(list(matrix(0, ncol = ncol(honest_multi_demeaned[[1]][[1]]), nrow = nrow(honest_multi_demeaned[[1]][[1]]))), nrow(honest_multi_demeaned[[1]][[1]]) )))
  honest_multi_demeaned_0_first <- append(list(rep(list(matrix(0, ncol = ncol(honest_multi_demeaned[[1]][[1]]), nrow = nrow(honest_multi_demeaned[[1]][[1]]))), nrow(honest_multi_demeaned[[1]][[1]]) )), honest_multi_demeaned)
  # train sample
  train_multi_demeaned_0_last <- append(train_multi_demeaned, list(rep(list(matrix(0, ncol = ncol(train_multi_demeaned[[1]][[1]]), nrow = nrow(train_multi_demeaned[[1]][[1]]))), nrow(train_multi_demeaned[[1]][[1]]) )))
  train_multi_demeaned_0_first <- append(list(rep(list(matrix(0, ncol = ncol(train_multi_demeaned[[1]][[1]]), nrow = nrow(train_multi_demeaned[[1]][[1]]))), nrow(train_multi_demeaned[[1]][[1]]) )), train_multi_demeaned)

  # compute the multiplication of category m with m-1 according to the covariance formula
  honest_multi_demeaned_cov <- mapply(function(x,y) mapply(function(x,y) x*y, x, y, SIMPLIFY = FALSE), honest_multi_demeaned_0_first, honest_multi_demeaned_0_last, SIMPLIFY = FALSE)
  train_multi_demeaned_cov <- mapply(function(x,y) mapply(function(x,y) x*y, x, y, SIMPLIFY = FALSE), train_multi_demeaned_0_first, train_multi_demeaned_0_last, SIMPLIFY = FALSE)

  # sum all obs i together
  honest_multi_demeaned_cov_sum <- lapply(honest_multi_demeaned_cov, function(x) lapply(x, function(x) sum(x)))
  train_multi_demeaned_cov_sum <- lapply(train_multi_demeaned_cov, function(x) lapply(x, function(x) sum(x)))

  # multiply by N/N-1 (normalize)
  honest_multi_demeaned_cov_sum_norm <- lapply(honest_multi_demeaned_cov_sum, function(x) lapply(x, function(x) x*(length(honest_pred[[1]])/(length(honest_pred[[1]])-1)) ))
  train_multi_demeaned_cov_sum_norm <- lapply(train_multi_demeaned_cov_sum, function(x) lapply(x, function(x) x*(length(train_pred[[1]])/(length(train_pred[[1]])-1)) ))

  # multiply by 2
  honest_multi_demeaned_cov_sum_norm_mult2 <- lapply(honest_multi_demeaned_cov_sum_norm, function(x) lapply(x, function(x) x*2 ))
  train_multi_demeaned_cov_sum_norm_mult2 <- lapply(train_multi_demeaned_cov_sum_norm, function(x) lapply(x, function(x) x*2 ))

  # put it into a shorter named object
  honest_covariance <- honest_multi_demeaned_cov_sum_norm_mult2
  train_covariance <- train_multi_demeaned_cov_sum_norm_mult2

  # ----------------------------------------------------------------------------------- #

  ## put everything together according to the whole variance formula
  # shift variances accordingly for ease of next computations (covariance already has the desired format)
  # honest sample
  honest_variance_last <- append(honest_variance, list(rep(list(0), nrow(honest_multi_demeaned[[1]][[1]]) ))) # append zero element list
  honest_variance_first <- append(list(rep(list(0), nrow(honest_multi_demeaned[[1]][[1]]) )), honest_variance) # prepend zero element list
  # train sample
  train_variance_last <- append(train_variance, list(rep(list(0), nrow(train_multi_demeaned[[1]][[1]]) ))) # append zero element list
  train_variance_first <- append(list(rep(list(0), nrow(train_multi_demeaned[[1]][[1]]) )), train_variance) # prepend zero element list

  # put everything together according to formula: var_last + var_first - cov
  honest_variance_final <- mapply(function(x,y,z) mapply(function(x,y,z) x+y-z, x, y, z, SIMPLIFY = FALSE), honest_variance_last, honest_variance_first, honest_covariance, SIMPLIFY = FALSE)
  train_variance_final <- mapply(function(x,y,z) mapply(function(x,y,z) x+y-z, x, y, z, SIMPLIFY = FALSE), train_variance_last, train_variance_first, train_covariance, SIMPLIFY = FALSE)

  ## output for final variances
  # coerce to a matrix
  honest_var <- sapply(honest_variance_final, function(x) sapply(x, function(x) as.matrix(x)))
  train_var <- sapply(train_variance_final, function(x) sapply(x, function(x) as.matrix(x)))

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


#' Predict ORF Variance
#'
#' predict variance of ordered random forest predictions based on honest sample
#' splitting as described in Lechner (2018)
#'
#' @param honest_pred list of vectors of honest forest predictions
#' @param honest_weights list of n x n matrices of honest forest weights
#' @param Y_ind_honest list of vectors of 0-1 outcomes for the honest sample
#'
#' @return vector of ORF variances
#'
#' @keywords internal
#'
pred_orf_variance <- function(honest_pred, honest_weights, Y_ind_honest) {

  # ----------------------------------------------------------------------------------- #

  # get categories
  categories <- seq(1:(length(Y_ind_honest)+1))

  # ----------------------------------------------------------------------------------- #

  ## single variances computation
  # compute the conditional means (predictions): already have this as forest_pred
  # divide it by N to get the "mean"
  honest_pred_mean <- lapply(honest_pred, function(x) x/length(Y_ind_honest[[1]]))

  # calculate standard multiplication of weights and outcomes: honest_weights*y_ind_honest (note with seq_along: as many rows as honest_pred or honest_weights)
  honest_multi <- mapply(function(x,y) lapply(seq_along(x[, 1]), function(i) x[i, ] * y), honest_weights, Y_ind_honest, SIMPLIFY = FALSE)

  # subtract the mean from each obs i
  honest_multi_demeaned <- mapply(function(x,y) mapply(function(x,y) x-y, x, y, SIMPLIFY = FALSE), honest_multi, honest_pred_mean, SIMPLIFY = FALSE)

  ## now do the single variances for each category m
  # square the demeaned
  honest_multi_demeaned_sq <- lapply(honest_multi_demeaned, function(x) lapply(x, function(x) x^2))

  # sum all obs i together
  honest_multi_demeaned_sq_sum <- lapply(honest_multi_demeaned_sq, function(x) lapply(x, function(x) sum(x)))

  # multiply by N/N-1 (normalize)
  honest_multi_demeaned_sq_sum_norm <- lapply(honest_multi_demeaned_sq_sum, function(x) lapply(x, function(x) x*(length(honest_pred[[1]])/(length(honest_pred[[1]])-1)) ))

  # put it into a shorter named object
  honest_variance <- honest_multi_demeaned_sq_sum_norm

  # ----------------------------------------------------------------------------------- #

  ##  covariances computation
  # multiply forest_var_multi_demeaned according to formula for covariance (shifted categories needed for computational convenience)
  # honest sample
  honest_multi_demeaned_0_last <- append(honest_multi_demeaned, list(rep(list(rep(0, length(honest_multi_demeaned[[1]][[1]]))), length(honest_multi_demeaned[[1]]))))
  honest_multi_demeaned_0_first <- append(list(rep(list(rep(0, length(honest_multi_demeaned[[1]][[1]]))), length(honest_multi_demeaned[[1]]))), honest_multi_demeaned)

  # compute the multiplication of category m with m-1 according to the covariance formula
  honest_multi_demeaned_cov <- mapply(function(x,y) mapply(function(x,y) x*y, x, y, SIMPLIFY = FALSE), honest_multi_demeaned_0_first, honest_multi_demeaned_0_last, SIMPLIFY = FALSE)

  # sum all obs i together
  honest_multi_demeaned_cov_sum <- lapply(honest_multi_demeaned_cov, function(x) lapply(x, function(x) sum(x)))

  # multiply by N/N-1 (normalize)
  honest_multi_demeaned_cov_sum_norm <- lapply(honest_multi_demeaned_cov_sum, function(x) lapply(x, function(x) x*(length(honest_pred[[1]])/(length(honest_pred[[1]])-1)) ))

  # multiply by 2
  honest_multi_demeaned_cov_sum_norm_mult2 <- lapply(honest_multi_demeaned_cov_sum_norm, function(x) lapply(x, function(x) x*2 ))

  # put it into a shorter named object
  honest_covariance <- honest_multi_demeaned_cov_sum_norm_mult2

  # ----------------------------------------------------------------------------------- #

  ## put everything together according to the whole variance formula
  # shift variances accordingly for ease of next computations (covariance already has the desired format)
  # honest sample
  honest_variance_last <- append(honest_variance, list(rep(list(0), length(honest_multi_demeaned[[1]])))) # append zero element list
  honest_variance_first <- append(list(rep(list(0), length(honest_multi_demeaned[[1]]))), honest_variance) # prepend zero element list

  # put everything together according to formula: var_last + var_first - cov
  honest_variance_final <- mapply(function(x,y,z) mapply(function(x,y,z) x+y-z, x, y, z, SIMPLIFY = FALSE), honest_variance_last, honest_variance_first, honest_covariance, SIMPLIFY = FALSE)

  ## output for final variances
  # coerce to a matrix
  honest_var <- sapply(honest_variance_final, function(x) sapply(x, function(x) as.matrix(x)))

  # ----------------------------------------------------------------------------------- #

  # save as forest_var
  forest_variance <- honest_var

  # add names
  colnames(forest_variance) <- sapply(categories, function(x) paste("Category", x, sep = " "))

  # ----------------------------------------------------------------------------------- #

  ## return the matrix
  output <- forest_variance
  # output
  return(output)

  # ----------------------------------------------------------------------------------- #

}

#' ORF Predictions for Marginal Effects
#'
#' Fast ORF Predictions for estimation of marginal effects at mean
#'
#' @param forest list of ranger forest objects
#' @param data list of n x n matrices of indicator data for orf
#' @param pred_data list of prediction data (X_mean_up/down)
#'
#' @return list of predictions
#'
#' @keywords internal
#'
predict_forest_preds_for_ME <- function(forest, data, pred_data) {

  # ----------------------------------------------------------------------------------- #

  # get number of observations
  n_col <- ncol(data[[1]])-1 # number of X variables
  n_row <- nrow(pred_data[[1]]) # number of evaluation points rows
  # stack pred_data together and make predictions at once
  pred_data <- do.call(rbind, pred_data)
  # prepare empty list
  forest_preds_up_together <- rep(list(NA), length(forest))
  forest_preds_up <- rep(list(rep(list(NA), n_col)), length(forest))

  # ----------------------------------------------------------------------------------- #

  # predict for desired Xs up
  for (forest_index in seq_along(forest)) {
    # start index for observations of X
    start_index <- 0
    # make honest predictions
    forest_preds_up_together[[forest_index]] <- predict_honest(forest[[forest_index]], data[[forest_index]], pred_data)
    # now assign for each X into respective forest as list entries
    for (X_index in seq_along(pred_data[1,])) {

      stop_index <- (n_row*X_index)
      forest_preds_up[[forest_index]][[X_index]] <- mean(forest_preds_up_together[[forest_index]][(1+start_index):stop_index]) # mean doesnt matter for atmean or atmedian
      start_index <- stop_index

    }
  }

  # ----------------------------------------------------------------------------------- #

  # return predictions
  return(forest_preds_up)

  # ----------------------------------------------------------------------------------- #

}


#' ORF Weight Predictions for Marginal Effects
#'
#' Fast ORF Weight Predictions for estimation of marginal effects at mean
#'
#' @param forest list of ranger forest objects
#' @param data list of n x n matrices of indicator data for orf
#' @param pred_data list of prediction data (X_mean_up/down)
#'
#' @return list of weights
#'
#' @keywords internal
#'
predict_forest_weights_for_ME <- function(forest, data, pred_data) {

  # get leafs only for number of forests and then check the Xs on their own
  # prepare empty list
  forest_weights_up <- rep(list(rep(list(NA), ncol(data))), length(forest))

  # ----------------------------------------------------------------------------------- #

  # extract weights for desired Xs up
  for (forest_index in seq_along(forest)) {

    # first get terminal nodes, i.e get terminal nodes for each obs and each tree
    # run your (new) data through the forest structure and look where your observations end up
    leaf_IDs <- predict(forest[[forest_index]], data, type = "terminalNodes")$predictions
    # put leaf_IDs into a list (one element for one tree)
    leaf_IDs <- lapply(seq_along(leaf_IDs[1,]), function(i) leaf_IDs[,i])
    # get the leaf size as counts of observations in leaves
    leaf_size <- lapply(leaf_IDs, function(x) ave(x, x, FUN=length))

    # start looping over all X_means
    for (X_index in seq_along(pred_data)) {

      # now do the same for your prediction data
      leaf_IDs_pred <- predict(forest[[forest_index]], as.matrix(pred_data[[X_index]]), type = "terminalNodes")$predictions
      # put leaf_IDs into a list
      leaf_IDs_pred <- lapply(seq_along(leaf_IDs_pred[1,]), function(i) leaf_IDs_pred[,i])

      # now average over the bootstraps, i.e. over trees to get final weights
      forest_weights_up[[forest_index]][[X_index]] <- pred_weights_C(leaf_IDs_pred, leaf_IDs, leaf_size, 1)

    }

  }

  # ----------------------------------------------------------------------------------- #

  # return weights
  return(forest_weights_up)

  # ----------------------------------------------------------------------------------- #

}

#' Get ORF Variance
#'
#' get variance of ordered random forest predictions based on honest sample
#' splitting as described in Lechner (2018)
#'
#' @param honest_pred list of vectors of honest forest predictions
#' @param honest_weights list of n x n matrices of honest forest weights
#' @param train_pred list of vectors of honest forets predictions from train sample
#' @param train_weights list of vectors of honest forets predictions from train sample
#' @param Y_ind_honest list of vectors of 0-1 outcomes for the honest sample
#'
#' @return vector of ORF variances
get_orf_variance <- function(honest_pred, honest_weights, train_pred, train_weights, Y_ind_honest) {

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
  categories <- seq(1:(length(Y_ind_honest)+1))

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

  ## now do the single variances for each category m
  # square the demeaned
  honest_multi_demeaned_sq <- lapply(honest_multi_demeaned, function(x) lapply(x, function(x) x^2))
  train_multi_demeaned_sq <- lapply(train_multi_demeaned, function(x) lapply(x, function(x) x^2))

  # sum all obs i together
  honest_multi_demeaned_sq_sum <- lapply(honest_multi_demeaned_sq, function(x) lapply(x, function(x) sum(x)))
  train_multi_demeaned_sq_sum <- lapply(train_multi_demeaned_sq, function(x) lapply(x, function(x) sum(x)))

  # divide by scaling factor
  #forest_multi_demeaned_sq_sum_scaled <- lapply(forest_multi_demeaned_sq_sum, function(x) mapply(function(x,y) x/y, x, scaling_factor_squared, SIMPLIFY = FALSE) )
  # multiply by N/N-1 (normalize)
  honest_multi_demeaned_sq_sum_norm <- lapply(honest_multi_demeaned_sq_sum, function(x) lapply(x, function(x) x*(length(honest_pred[[1]])/(length(honest_pred[[1]])-1)) ))
  train_multi_demeaned_sq_sum_norm <- lapply(train_multi_demeaned_sq_sum, function(x) lapply(x, function(x) x*(length(train_pred[[1]])/(length(train_pred[[1]])-1)) ))

  # put it into a shorter named object
  honest_variance <- honest_multi_demeaned_sq_sum_norm
  train_variance <- train_multi_demeaned_sq_sum_norm
  ## single variances done

  # ----------------------------------------------------------------------------------- #

  ## now compute the covariances
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
  ## covariances done

  # ----------------------------------------------------------------------------------- #

  ## put everything together according to the whole variance formula
  # shift variances accordingly for ease of next computations (covariance already has the desired format)
  # honest sample
  honest_variance_last <- append(honest_variance, list(rep(list(0), nrow(honest_multi_demeaned[[1]][[1]]) ))) # append zero element list
  honest_variance_first <- append(list(rep(list(0), nrow(honest_multi_demeaned[[1]][[1]]) )), honest_variance) # prepend zero element list
  # train sample
  train_variance_last <- append(train_variance, list(rep(list(0), nrow(train_multi_demeaned[[1]][[1]]) ))) # append zero element list
  train_variance_first <- append(list(rep(list(0), nrow(train_multi_demeaned[[1]][[1]]) )), train_variance) # prepend zero element list

  ## put everything together according to formula: var_last + var_first - cov
  honest_variance_final <- mapply(function(x,y,z) mapply(function(x,y,z) x+y-z, x, y, z, SIMPLIFY = FALSE), honest_variance_last, honest_variance_first, honest_covariance, SIMPLIFY = FALSE)
  train_variance_final <- mapply(function(x,y,z) mapply(function(x,y,z) x+y-z, x, y, z, SIMPLIFY = FALSE), train_variance_last, train_variance_first, train_covariance, SIMPLIFY = FALSE)

  ## output for final variances of marginal effects
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
pred_orf_variance <- function(honest_pred, honest_weights, Y_ind_honest) {

  # needed inputs for the function: honest_pred - list of vectors of honest forest predictions
  #                                 honest_weights - list of n x n matrices of honest forest weights
  #                                 Y_ind_honest - list of vectors of 0-1 outcomes for the honest sample

  # ----------------------------------------------------------------------------------- #

  # get categories
  categories <- seq(1:(length(Y_ind_honest)+1))

  # ----------------------------------------------------------------------------------- #

  ### variance for the orf predictions
  ## compute prerequisities for variance of honest orf predictions
  ## compute the conditional means (predictions): already have this as forest_pred
  # divide it by N to get the "mean"
  honest_pred_mean <- lapply(honest_pred, function(x) x/length(Y_ind_honest[[1]]))

  # calculate standard multiplication of weights and outcomes: honest_weights*y_ind_honest (be careful with seq_along!) you have to go as many rows as honest_pred or honest_weights have
  honest_multi <- mapply(function(x,y) lapply(seq_along(x[, 1]), function(i) x[i, ] * y), honest_weights, Y_ind_honest, SIMPLIFY = FALSE)  # subtract the mean from each obs i

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
  ## single variances done

  # ----------------------------------------------------------------------------------- #

  ## now compute the covariances
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
  ## covariances done

  # ----------------------------------------------------------------------------------- #

  ## put everything together according to the whole variance formula
  # shift variances accordingly for ease of next computations (covariance already has the desired format)
  # honest sample
  honest_variance_last <- append(honest_variance, list(rep(list(0), length(honest_multi_demeaned[[1]])))) # append zero element list
  honest_variance_first <- append(list(rep(list(0), length(honest_multi_demeaned[[1]]))), honest_variance) # prepend zero element list

  ## put everything together according to formula: var_last + var_first - cov
  honest_variance_final <- mapply(function(x,y,z) mapply(function(x,y,z) x+y-z, x, y, z, SIMPLIFY = FALSE), honest_variance_last, honest_variance_first, honest_covariance, SIMPLIFY = FALSE)

  ## output for final variances of marginal effects
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


#' ORF margins
#'
#' estimate marginal effects for ordered random forests as defined in
#' Lechner & Okasa (2019) for in-sample
#'
#' @param forest list of ranger forest objects
#' @param data list of n x n matrices of indicator data for orf
#' @param honesty logical if HONESTY should be conducted
#' @param inference logical if INFERENCE should be conducted
#'
#' @importFrom stats predict median pnorm sd
#' @import ranger
#'
#' @return matrix of marginal effects estimates
orf_margins <- function(forest, data, honesty, inference){

  # needed inputs for the function: forest - list of ranger forest objects
  #                                 data - list of n x n matrices of indicator data for orf
  #                                 honesty - logical if HONESTY should be conducted
  #                                 inference - logical if INFERENCE should be conducted

  # ----------------------------------------------------------------------------------- #

  # get number of observations
  n_data <- nrow(data[[1]])
  # get X as matrix
  X <- as.matrix(data[[1]][, -(1)])
  # get categories
  categories <- seq(1:(length(data)+1))
  # evaluate at mean always here in the main orf function
  eval <- "mean"

  # decide if inference should be done or not
  if (inference == FALSE & honesty == FALSE) {

    ## values for evaluation of the marginal effect
    # share of SD to be used
    h_std <- 0.1

    # check if X is continuous or dummy or categorical
    X_type <- apply(X, 2, function(x) length(unique(x)))
    # now determine the type of X
    X_continuous <- which(X_type > 10) # define IDs of continuous Xs
    X_dummy <- which(X_type == 2) # define IDs of dummies
    X_categorical <- which(X_type > 2 & X_type <= 10)
    # additional check for constant variables which are nonsensical
    if (any(X_type == 1) | any(X_type == 0)) {
      stop("Some of the covariates are constant. This is non-sensical for evaluation of marginal effects. Programme terminated.")
    }

    # decide if the marginal effects should be computed at mean or at median
    if (eval=="mean") {
      # variable of interest: X_1 to X_last, ME at mean
      X_mean <- lapply(1:ncol(X), function(x) colMeans(X)) # set all Xs to their mean values (so many times as we have Xs)
    } else if (eval=="median") {
      # variable of interest: X_1 to X_last, ME at median
      X_mean <- lapply(1:ncol(X), function(x) apply(X, 2, median)) # set all Xs to their median values (so many times as we have Xs)
    } else {
      stop("Incorrect evaluation point. Programme terminated.")
    }

    # ----------------------------------------------------------------------------------- #

    # get SD of Xs
    X_sd <- apply(X, 2, sd)
    # create X_up (X_mean + 0.1 * X_sd)
    X_up <- X_mean[[1]] + h_std*X_sd
    # create X_down (X_mean - 0.1 * X_sd)
    X_down <- X_mean[[1]] - h_std*X_sd

    ## now check for the support of X
    # check X_max
    X_max <- apply(X, 2, max)
    # check X_min
    X_min <- apply(X, 2, min)
    # check if X_up is within the range X_min and X_max
    X_up <- (X_up < X_max) * X_up + (X_up >= X_max) * X_max
    X_up <- (X_up > X_min) * X_up + (X_up <= X_min) * (X_min + h_std * X_sd)
    # check if X_down is within the range X_min and X_max
    X_down <- (X_down > X_min) * X_down + (X_down <= X_min) * X_min
    X_down <- (X_down < X_max) * X_down + (X_down >= X_max) * (X_max - h_std * X_sd)
    # some additional checks

    # ----------------------------------------------------------------------------------- #

    ## now we need 2 datasets: one with X_up and second with X_down
    # X_mean_up continous
    X_mean_up <- lapply(seq_along(X_mean), function(i) replace(X_mean[[i]], i, X_up[i]) )
    # X_mean_down continous
    X_mean_down <- lapply(seq_along(X_mean), function(i) replace(X_mean[[i]], i, X_down[i]) )

    # adjust for categorical X (works also for zero categorical)
    for (i in X_categorical) {
      X_mean_up[[i]][[i]] <- ceiling(mean(X[,i]))
      X_mean_down[[i]][[i]] <- floor(mean(X[,i]))
    }
    # adjust for dummies (works also for zero dummies)
    for (i in X_dummy) {
      X_mean_up[[i]][[i]] <- max(X[,i])
      X_mean_down[[i]][[i]] <- min(X[,i])
    }

    # ----------------------------------------------------------------------------------- #

    #### now we do not need weights if we do not need inference (based on out of bag predictions)
    # forest prediction for X_mean_up
    forest_pred_up <- lapply(forest, function(x) lapply(X_mean_up, function(y) predict(x, data = t(as.matrix((y))))$predictions))
    # forest prediction for X_mean_down
    forest_pred_down <- lapply(forest, function(x) lapply(X_mean_down, function(y) predict(x, data = t(as.matrix((y))))$predictions))
    # now subtract the predictions according to the ME formula
    forest_pred_diff_up_down <- mapply(function(x,y) mapply(function(x,y) x-y, x, y,  SIMPLIFY = F), forest_pred_up, forest_pred_down, SIMPLIFY = F)

    # subtract predictions according to formula to isolate categories
    forest_cond_means_0_last <- append(forest_pred_diff_up_down, list(rep(list(0), ncol(X)))) # append zero elemnt list
    forest_cond_means_0_first <- append(list(rep(list(0), ncol(X))), forest_pred_diff_up_down) # prepend zero element list

    # compute the scaling factor: X_up-X_down=2*X_sd
    scaling_factor <- as.list(X_up - X_down)
    # set scaling factor to zero for categorical and dummy variables
    for (i in X_categorical & X_dummy) {
      scaling_factor[[i]] <- 0
    }

    # now compute the differences for marginal effects
    marginal_effects_diff <- mapply(function(x,y) mapply(function(x,y) x-y, x, y, SIMPLIFY = F), forest_cond_means_0_last, forest_cond_means_0_first, SIMPLIFY = F)

    # scale marginal effects
    marginal_effects_scaled <- lapply(marginal_effects_diff, function(x) mapply(function(x,y) x/y, x, scaling_factor, SIMPLIFY = FALSE) )

    ## output for final marginal effects
    # coerce to a matrix
    marginal_effects <- sapply(marginal_effects_scaled, function(x) sapply(x, function(x) as.matrix(x)))
    # add names
    colnames(marginal_effects) <- sapply(categories, function(x) paste("Category", x, sep = " "))
    rownames(marginal_effects) <- colnames(X)

    # put marginal effects into results
    results <- marginal_effects

    # ------------------------------------------------------------------------------------- #

  } else if (inference == FALSE & honesty == TRUE) {

    # ----------------------------------------------------------------------------------- #

    ## values for evaluation of the marginal effect - based on honest sample
    # share of SD to be used
    h_std <- 0.1

    # check if X is continuous or dummy or categorical (now X_honest)
    X_type <- apply(X, 2, function(x) length(unique(x)))
    # now determine the type of X
    X_continuous <- which(X_type > 10) # define IDs of continuous Xs
    X_dummy <- which(X_type == 2) # define IDs of dummies
    X_categorical <- which(X_type > 2 & X_type <= 10)
    # additional check for constant variables which are nonsensical
    if (any(X_type == 1) | any(X_type == 0)) {
      stop("Some of the covariates are constant. This is non-sensical for evaluation of marginal effects. Programme terminated.")
    }

    # decide if the marginal effects should be computed at mean or at median
    if (eval=="mean") {
      # variable of interest: X_1 to X_last, ME at mean
      X_mean <- lapply(1:ncol(X), function(x) colMeans(X)) # set all Xs to their mean values (so many times as we have Xs)
    } else if (eval=="median") {
      # variable of interest: X_1 to X_last, ME at median
      X_mean <- lapply(1:ncol(X), function(x) apply(X, 2, median)) # set all Xs to their median values (so many times as we have Xs)
    } else {
      stop("Incorrect evaluation point. Programme terminated.")
    }

    # ----------------------------------------------------------------------------------- #

    # get SD of Xs
    X_sd <- apply(X, 2, sd)
    # create X_up (X_mean + 0.1 * X_sd)
    X_up <- X_mean[[1]] + h_std*X_sd
    # create X_down (X_mean - 0.1 * X_sd)
    X_down <- X_mean[[1]] - h_std*X_sd

    ## now check for the support of X
    # check X_max
    X_max <- apply(X, 2, max)
    # check X_min
    X_min <- apply(X, 2, min)
    # check if X_up is within the range X_min and X_max
    X_up <- (X_up < X_max) * X_up + (X_up >= X_max) * X_max
    X_up <- (X_up > X_min) * X_up + (X_up <= X_min) * (X_min + h_std * X_sd)
    # check if X_down is within the range X_min and X_max
    X_down <- (X_down > X_min) * X_down + (X_down <= X_min) * X_min
    X_down <- (X_down < X_max) * X_down + (X_down >= X_max) * (X_max - h_std * X_sd)
    # some additional checks

    # ----------------------------------------------------------------------------------- #

    ## now we need 2 datasets: one with X_up and second with X_down
    # X_mean_up
    X_mean_up <- lapply(seq_along(X_mean), function(i) replace(X_mean[[i]], i, X_up[i]) )
    # X_mean_down
    X_mean_down <- lapply(seq_along(X_mean), function(i) replace(X_mean[[i]], i, X_down[i]) )

    # adjust for categorical X (works also for zero categorical)
    for (i in X_categorical) {
      X_mean_up[[i]][[i]] <- ceiling(mean(X[,i]))
      X_mean_down[[i]][[i]] <- floor(mean(X[,i]))
    }
    # adjust for dummies (works also for zero dummies)
    for (i in X_dummy) {
      X_mean_up[[i]][[i]] <- max(X[,i])
      X_mean_down[[i]][[i]] <- min(X[,i])
    }

    # ----------------------------------------------------------------------------------- #

    #### now we do not need weights if we do not need inference (based on out of bag predictions)
    # forest prediction for X_mean_up (use new faster function particularly for ME)
    forest_pred_up <- predict_forest_preds_for_ME(forest, data, X_mean_up)
    # forest prediction for X_mean_down
    forest_pred_down <- predict_forest_preds_for_ME(forest, data, X_mean_down)

    # now subtract the predictions according to the ME formula
    forest_pred_diff_up_down <- mapply(function(x,y) mapply(function(x,y) x-y, x, y,  SIMPLIFY = F), forest_pred_up, forest_pred_down, SIMPLIFY = F)

    # subtract predictions according to formula to isolate categories
    forest_cond_means_0_last <- append(forest_pred_diff_up_down, list(rep(list(0), ncol(X)))) # append zero elemnt list
    forest_cond_means_0_first <- append(list(rep(list(0), ncol(X))), forest_pred_diff_up_down) # prepend zero element list

    # compute the scaling factor: X_up-X_down=2*X_sd
    scaling_factor <- as.list(X_up - X_down)
    # set scaling factor to zero for categorical and dummy variables
    for (i in X_categorical & X_dummy) {
      scaling_factor[[i]] <- 0
    }

    # now compute the differences for marginal effects
    marginal_effects_diff <- mapply(function(x,y) mapply(function(x,y) x-y, x, y, SIMPLIFY = F), forest_cond_means_0_last, forest_cond_means_0_first, SIMPLIFY = F)

    # scale marginal effects
    marginal_effects_scaled <- lapply(marginal_effects_diff, function(x) mapply(function(x,y) x/y, x, scaling_factor, SIMPLIFY = FALSE) )

    ## output for final marginal effects
    # coerce to a matrix
    marginal_effects <- sapply(marginal_effects_scaled, function(x) sapply(x, function(x) as.matrix(x)))
    # add names
    colnames(marginal_effects) <- sapply(categories, function(x) paste("Category", x, sep = " "))
    rownames(marginal_effects) <- colnames(X)

    # put marginal effects into results
    results <- marginal_effects

    # ---------------------------------------------------------------------------------------------------- #

  } else if (inference == TRUE & honesty == TRUE) {

    # ----------------------------------------------------------------------------------- #

    # take out honest indicator Ys
    Y_ind_honest <- lapply(data, function(x) as.matrix(x[, 1]))

    ## values for evaluation of the marginal effect - based on honest sample
    # share of SD to be used
    h_std <- 0.1

    # check if X is continuous or dummy or categorical
    X_type <- apply(X, 2, function(x) length(unique(x)))
    # now determine the type of X
    X_continuous <- which(X_type > 10) # define IDs of continuous Xs
    X_dummy <- which(X_type == 2) # define IDs of dummies
    X_categorical <- which(X_type > 2 & X_type <= 10)
    # additional check for constant variables which are nonsensical
    if (any(X_type == 1) | any(X_type == 0)) {
      stop("Some of the covariates are constant. This is non-sensical for evaluation of marginal effects. Programme terminated.")
    }

    # decide if the marginal effects should be computed at mean or at median
    if (eval=="mean") {
      # variable of interest: X_1 to X_last, ME at mean
      X_mean <- lapply(1:ncol(X), function(x) colMeans(X)) # set all Xs to their mean values (so many times as we have Xs)
    } else if (eval=="median") {
      # variable of interest: X_1 to X_last, ME at median
      X_mean <- lapply(1:ncol(X), function(x) apply(X, 2, median)) # set all Xs to their median values (so many times as we have Xs)
    } else {
      stop("Incorrect evaluation point. Programme terminated.")
    }

    # get SD of Xs
    X_sd <- apply(X, 2, sd)
    # create X_up (X_mean + 0.1 * X_sd)
    X_up <- X_mean[[1]] + h_std*X_sd
    # create X_down (X_mean - 0.1 * X_sd)
    X_down <- X_mean[[1]] - h_std*X_sd

    ## now check for the support of X
    # check X_max
    X_max <- apply(X, 2, max)
    # check X_min
    X_min <- apply(X, 2, min)
    # check if X_up is within the range X_min and X_max
    X_up <- (X_up < X_max) * X_up + (X_up >= X_max) * X_max
    X_up <- (X_up > X_min) * X_up + (X_up <= X_min) * (X_min + h_std * X_sd)
    # check if X_down is within the range X_min and X_max
    X_down <- (X_down > X_min) * X_down + (X_down <= X_min) * X_min
    X_down <- (X_down < X_max) * X_down + (X_down >= X_max) * (X_max - h_std * X_sd)
    # some additional checks

    # ----------------------------------------------------------------------------------- #

    ## now we need 2 datasets: one with X_up and second with X_down
    # X_mean_up
    X_mean_up <- lapply(seq_along(X_mean), function(i) replace(X_mean[[i]], i, X_up[i]) )
    # X_mean_down
    X_mean_down <- lapply(seq_along(X_mean), function(i) replace(X_mean[[i]], i, X_down[i]) )

    # adjust for categorical X (works also for zero categorical)
    for (i in X_categorical) {
      X_mean_up[[i]][[i]] <- ceiling(mean(X[,i]))
      X_mean_down[[i]][[i]] <- floor(mean(X[,i]))
    }
    # adjust for dummies (works also for zero dummies)
    for (i in X_dummy) {
      X_mean_up[[i]][[i]] <- max(X[,i])
      X_mean_down[[i]][[i]] <- min(X[,i])
    }

    # ----------------------------------------------------------------------------------- #

    # extract weights for desired Xs up: get weights from honest sample and predict weights for evaluation points from HONEST sample
    forest_weights_up <- predict_forest_weights_for_ME(forest, X, X_mean_up)
    # extract weights for desired Xs down
    forest_weights_down <- predict_forest_weights_for_ME(forest, X, X_mean_down)
    # now subtract the weights according to the ME formula
    forest_weights_diff_up_down <- mapply(function(x,y) mapply(function(x,y) x-y, x, y,  SIMPLIFY = F), forest_weights_up, forest_weights_down, SIMPLIFY = F)

    ## compute prerequisities for marginal effects
    # compute the conditional means: weights%*%y (predictions are based on honest sample)
    forest_cond_means <- mapply(function(x,y) lapply(x, function(x) x%*%y), forest_weights_diff_up_down, Y_ind_honest, SIMPLIFY = FALSE)
    # subtract conditional means according to formula to isolate categories
    forest_cond_means_0_last <- append(forest_cond_means, list(rep(list(0), ncol(X)))) # append zero elemnt list
    forest_cond_means_0_first <- append(list(rep(list(0), ncol(X))), forest_cond_means) # prepend zero element list
    # compute the scaling factor: X_up-X_down=2*X_sd
    scaling_factor <- as.list(X_up - X_down)
    # set scaling factor to zero for categorical and dummy variables
    for (i in X_categorical & X_dummy) {
      scaling_factor[[i]] <- 0
    }

    # now compute the differences for marginal effects
    marginal_effects_diff <- mapply(function(x,y) mapply(function(x,y) x-y, x, y, SIMPLIFY = F), forest_cond_means_0_last, forest_cond_means_0_first, SIMPLIFY = F)

    # scale marginal effects
    marginal_effects_scaled <- lapply(marginal_effects_diff, function(x) mapply(function(x,y) x/y, x, scaling_factor, SIMPLIFY = FALSE) )

    ## output for final marginal effects
    # coerce to a matrix
    marginal_effects <- sapply(marginal_effects_scaled, function(x) sapply(x, function(x) as.matrix(x)))
    # add names
    colnames(marginal_effects) <- sapply(categories, function(x) paste("Category", x, sep = " "))
    rownames(marginal_effects) <- colnames(X)

    # ----------------------------------------------------------------------------------- #

    ### variance for the marginal effects
    ## compute prerequisities for variance of honest marginal effects
    # scaling factor squared
    scaling_factor_squared <- lapply(scaling_factor, function(x) x^2)
    ## compute the conditional means (predictions): already have this as forest_cond_means
    # divide it by N to get the "mean"
    forest_cond_means_mean <- lapply(forest_cond_means, function(x) lapply(x, function(x) x/n_data))
    # calculate standard multiplication of weights and outcomes: honest_weights*y_ind_honest
    forest_multi <- mapply(function(x,y) lapply(x, function(x) t(x)*y), forest_weights_diff_up_down, Y_ind_honest, SIMPLIFY = FALSE)
    # subtract the mean from each obs i
    forest_multi_demeaned <- mapply(function(x,y) mapply(function(x,y) x-matrix(y, nrow = nrow(x)), x, y, SIMPLIFY = FALSE), forest_multi, forest_cond_means_mean, SIMPLIFY = F)

    ## now do the single variances for each category m
    # square the demeaned
    forest_multi_demeaned_sq <- lapply(forest_multi_demeaned, function(x) lapply(x, function(x) x^2))
    # sum all obs i together
    forest_multi_demeaned_sq_sum <- lapply(forest_multi_demeaned_sq, function(x) lapply(x, function(x) sum(x)))
    # divide by scaling factor
    forest_multi_demeaned_sq_sum_scaled <- lapply(forest_multi_demeaned_sq_sum, function(x) mapply(function(x,y) x/y, x, scaling_factor_squared, SIMPLIFY = FALSE) )
    # multiply by N/N-1 (normalize)
    forest_multi_demeaned_sq_sum_scaled_norm <- lapply(forest_multi_demeaned_sq_sum_scaled, function(x) lapply(x, function(x) x*(n_data/(n_data-1)) ))
    # put it into a shorter named object
    variance <- forest_multi_demeaned_sq_sum_scaled_norm
    ## single variances done

    # ----------------------------------------------------------------------------------- #

    ## now compute the covariances
    # multiply forest_var_multi_demeaned according to formula for covariance (shifted categories needed for computational convenience)
    forest_multi_demeaned_0_last <- append(forest_multi_demeaned, list(rep(list(matrix(0, ncol = ncol(forest_multi_demeaned[[1]][[1]]), nrow = nrow(forest_multi_demeaned[[1]][[1]]))), ncol(X)))) # append zero matrix list
    forest_multi_demeaned_0_first <- append(list(rep(list(matrix(0, ncol = ncol(forest_multi_demeaned[[1]][[1]]), nrow = nrow(forest_multi_demeaned[[1]][[1]]))), ncol(X))), forest_multi_demeaned) # prepend zero matrix list
    # compute the multiplication of category m with m-1 according to the covariance formula
    forest_multi_demeaned_cov <- mapply(function(x,y) mapply(function(x,y) x*y, x, y, SIMPLIFY = FALSE), forest_multi_demeaned_0_first, forest_multi_demeaned_0_last, SIMPLIFY = F)
    # sum all obs i together
    forest_multi_demeaned_cov_sum <- lapply(forest_multi_demeaned_cov, function(x) lapply(x, function(x) sum(x)))
    # divide by scaling factor
    forest_multi_demeaned_cov_sum_scaled <- lapply(forest_multi_demeaned_cov_sum, function(x) mapply(function(x,y) x/y, x, scaling_factor_squared, SIMPLIFY = FALSE) )
    # multiply by N/N-1 (normalize)
    forest_multi_demeaned_cov_sum_scaled_norm <- lapply(forest_multi_demeaned_cov_sum_scaled, function(x) lapply(x, function(x) x*(n_data/(n_data-1)) ))
    # multiply by 2
    forest_multi_demeaned_cov_sum_scaled_norm_mult2 <- lapply(forest_multi_demeaned_cov_sum_scaled_norm, function(x) lapply(x, function(x) x*2 ))
    # put it into a shorter named object
    covariance <- forest_multi_demeaned_cov_sum_scaled_norm_mult2
    ## covariances done

    # ----------------------------------------------------------------------------------- #

    ## put everything together according to the whole variance formula
    # shift variances accordingly for ease of next computations (covariance already has the desired format)
    variance_last <- append(variance, list(rep(list(0), ncol(X)))) # append zero element list
    variance_first <- append(list(rep(list(0), ncol(X))), variance) # prepend zero element list
    # put everything together according to formula: var_last + var_first - cov
    variance_marginal_effects_final <- mapply(function(x,y,z) mapply(function(x,y,z) x+y-z, x, y, z, SIMPLIFY = FALSE), variance_last, variance_first, covariance, SIMPLIFY = F)

    # ----------------------------------------------------------------------------------- #

    ## output for final variances of marginal effects
    # coerce to a matrix
    variance_marginal_effects <- sapply(variance_marginal_effects_final, function(x) sapply(x, function(x) as.matrix(x)))
    # add names
    colnames(variance_marginal_effects) <- sapply(categories, function(x) paste("Category", x, sep = " "))
    rownames(variance_marginal_effects) <- colnames(X)

    # ----------------------------------------------------------------------------------- #

    ## standard deviations
    # take square root of variance
    sd_marginal_effects <- sqrt(variance_marginal_effects)

    #### z scores and p values ####
    z_scores <- (marginal_effects)/(sd_marginal_effects)
    # control for dividing zero by zero
    z_scores[is.nan(z_scores)] = 0
    # p values
    p_values <- 2*pnorm(-abs(z_scores))

    # ----------------------------------------------------------------------------------- #

    # put everzthing into a list of results
    results <- list(marginal_effects, variance_marginal_effects, sd_marginal_effects, p_values)
    names(results) <- c("MarginalEffects", "Variances", "StandardErrors", "pValues")

    # ----------------------------------------------------------------------------------- #

  }

  # return results
  return(results)

  # ----------------------------------------------------------------------------------- #

}


#' Predict ORF Margins
#'
#' estimate marginal effects for ordered random forests as defined in
#' Lechner & Okasa (2019) for prediction sample
#'
#' @param forest list of ranger forest objects
#' @param data list of n x n matrices of indicator data for orf
#' @param newdata matrix of new Xs
#' @param honesty logical if HONESTY should be conducted
#' @param inference logical if INFERENCE should be conducted
#'
#' @importFrom stats predict median pnorm sd
#' @import ranger
#'
#' @return matrix of marginal effects estimates
pred_orf_margins <- function(forest, data, newdata, honesty, inference){

  # needed inputs for the function: forest - list of ranger forest objects
  #                                 data - list of n x n matrices of indicator data for orf
  #                                 newdata - matrix of new Xs
  #                                 honesty - logical if HONESTY should be conducted
  #                                 inference - logical if INFERENCE should be conducted

  # ----------------------------------------------------------------------------------- #

  # get number of observations
  n_data <- nrow(data[[1]])
  # get X as matrix
  X <- as.matrix(newdata)
  # get categories
  categories <- seq(1:(length(data)+1))
  # evaluate at mean always here in the main orf function
  eval <- "mean"

  # ----------------------------------------------------------------------------------- #

  # decide if inference should be done or not
  if (inference == FALSE & honesty == FALSE) {

    ## values for evaluation of the marginal effect
    # share of SD to be used
    h_std <- 0.1

    # check if X is continuous or dummy or categorical
    X_type <- apply(X, 2, function(x) length(unique(x)))
    # now determine the type of X
    X_continuous <- which(X_type > 10) # define IDs of continuous Xs
    X_dummy <- which(X_type == 2) # define IDs of dummies
    X_categorical <- which(X_type > 2 & X_type <= 10)
    # additional check for constant variables which are nonsensical
    if (any(X_type == 1) | any(X_type == 0)) {
      stop("Some of the covariates are constant. This is non-sensical for evaluation of marginal effects. Programme terminated.")
    }

    # decide if the marginal effects should be computed at mean or at median
    if (eval=="mean") {
      # variable of interest: X_1 to X_last, ME at mean
      X_mean <- lapply(1:ncol(X), function(x) colMeans(X)) # set all Xs to their mean values (so many times as we have Xs)
    } else if (eval=="median") {
      # variable of interest: X_1 to X_last, ME at median
      X_mean <- lapply(1:ncol(X), function(x) apply(X, 2, median)) # set all Xs to their median values (so many times as we have Xs)
    } else {
      stop("Incorrect evaluation point. Programme terminated.")
    }

    # ----------------------------------------------------------------------------------- #

    # get SD of Xs
    X_sd <- apply(X, 2, sd)
    # create X_up (X_mean + 0.1 * X_sd)
    X_up <- X_mean[[1]] + h_std*X_sd
    # create X_down (X_mean - 0.1 * X_sd)
    X_down <- X_mean[[1]] - h_std*X_sd

    ## now check for the support of X
    # check X_max
    X_max <- apply(X, 2, max)
    # check X_min
    X_min <- apply(X, 2, min)
    # check if X_up is within the range X_min and X_max
    X_up <- (X_up < X_max) * X_up + (X_up >= X_max) * X_max
    X_up <- (X_up > X_min) * X_up + (X_up <= X_min) * (X_min + h_std * X_sd)
    # check if X_down is within the range X_min and X_max
    X_down <- (X_down > X_min) * X_down + (X_down <= X_min) * X_min
    X_down <- (X_down < X_max) * X_down + (X_down >= X_max) * (X_max - h_std * X_sd)
    # some additional checks

    # ----------------------------------------------------------------------------------- #

    ## now we need 2 datasets: one with X_up and second with X_down
    # X_mean_up continous
    X_mean_up <- lapply(seq_along(X_mean), function(i) replace(X_mean[[i]], i, X_up[i]) )
    # X_mean_down continous
    X_mean_down <- lapply(seq_along(X_mean), function(i) replace(X_mean[[i]], i, X_down[i]) )

    # adjust for categorical X (works also for zero categorical)
    for (i in X_categorical) {
      X_mean_up[[i]][[i]] <- ceiling(mean(X[,i]))
      X_mean_down[[i]][[i]] <- floor(mean(X[,i]))
    }
    # adjust for dummies (works also for zero dummies)
    for (i in X_dummy) {
      X_mean_up[[i]][[i]] <- max(X[,i])
      X_mean_down[[i]][[i]] <- min(X[,i])
    }

    # ----------------------------------------------------------------------------------- #

    #### now we do not need weights if we do not need inference (based on out of bag predictions)
    # forest prediction for X_mean_up
    forest_pred_up <- lapply(forest, function(x) lapply(X_mean_up, function(y) predict(x, data = t(as.matrix((y))))$predictions))
    # forest prediction for X_mean_down
    forest_pred_down <- lapply(forest, function(x) lapply(X_mean_down, function(y) predict(x, data = t(as.matrix((y))))$predictions))
    # now subtract the predictions according to the ME formula
    forest_pred_diff_up_down <- mapply(function(x,y) mapply(function(x,y) x-y, x, y,  SIMPLIFY = F), forest_pred_up, forest_pred_down, SIMPLIFY = F)

    # subtract predictions according to formula to isolate categories
    forest_cond_means_0_last <- append(forest_pred_diff_up_down, list(rep(list(0), ncol(X)))) # append zero elemnt list
    forest_cond_means_0_first <- append(list(rep(list(0), ncol(X))), forest_pred_diff_up_down) # prepend zero element list

    # compute the scaling factor: X_up-X_down=2*X_sd
    scaling_factor <- as.list(X_up - X_down)
    # set scaling factor to zero for categorical and dummy variables
    for (i in X_categorical & X_dummy) {
      scaling_factor[[i]] <- 0
    }

    # now compute the differences for marginal effects
    marginal_effects_diff <- mapply(function(x,y) mapply(function(x,y) x-y, x, y, SIMPLIFY = F), forest_cond_means_0_last, forest_cond_means_0_first, SIMPLIFY = F)

    # scale marginal effects
    marginal_effects_scaled <- lapply(marginal_effects_diff, function(x) mapply(function(x,y) x/y, x, scaling_factor, SIMPLIFY = FALSE) )

    ## output for final marginal effects
    # coerce to a matrix
    marginal_effects <- sapply(marginal_effects_scaled, function(x) sapply(x, function(x) as.matrix(x)))
    # add names
    colnames(marginal_effects) <- sapply(categories, function(x) paste("Category", x, sep = " "))
    rownames(marginal_effects) <- colnames(X)

    # put marginal effects into results
    results <- marginal_effects

    # ------------------------------------------------------------------------------------- #

  } else if (inference == FALSE & honesty == TRUE) {

    # ----------------------------------------------------------------------------------- #

    ## values for evaluation of the marginal effect - based on honest sample
    # share of SD to be used
    h_std <- 0.1

    # check if X is continuous or dummy or categorical (now X_honest)
    X_type <- apply(X, 2, function(x) length(unique(x)))
    # now determine the type of X
    X_continuous <- which(X_type > 10) # define IDs of continuous Xs
    X_dummy <- which(X_type == 2) # define IDs of dummies
    X_categorical <- which(X_type > 2 & X_type <= 10)
    # additional check for constant variables which are nonsensical
    if (any(X_type == 1) | any(X_type == 0)) {
      stop("Some of the covariates are constant. This is non-sensical for evaluation of marginal effects. Programme terminated.")
    }

    # decide if the marginal effects should be computed at mean or at median
    if (eval=="mean") {
      # variable of interest: X_1 to X_last, ME at mean
      X_mean <- lapply(1:ncol(X), function(x) colMeans(X)) # set all Xs to their mean values (so many times as we have Xs)
    } else if (eval=="median") {
      # variable of interest: X_1 to X_last, ME at median
      X_mean <- lapply(1:ncol(X), function(x) apply(X, 2, median)) # set all Xs to their median values (so many times as we have Xs)
    } else {
      stop("Incorrect evaluation point. Programme terminated.")
    }

    # ----------------------------------------------------------------------------------- #

    # get SD of Xs
    X_sd <- apply(X, 2, sd)
    # create X_up (X_mean + 0.1 * X_sd)
    X_up <- X_mean[[1]] + h_std*X_sd
    # create X_down (X_mean - 0.1 * X_sd)
    X_down <- X_mean[[1]] - h_std*X_sd

    ## now check for the support of X
    # check X_max
    X_max <- apply(X, 2, max)
    # check X_min
    X_min <- apply(X, 2, min)
    # check if X_up is within the range X_min and X_max
    X_up <- (X_up < X_max) * X_up + (X_up >= X_max) * X_max
    X_up <- (X_up > X_min) * X_up + (X_up <= X_min) * (X_min + h_std * X_sd)
    # check if X_down is within the range X_min and X_max
    X_down <- (X_down > X_min) * X_down + (X_down <= X_min) * X_min
    X_down <- (X_down < X_max) * X_down + (X_down >= X_max) * (X_max - h_std * X_sd)
    # some additional checks

    # ----------------------------------------------------------------------------------- #

    ## now we need 2 datasets: one with X_up and second with X_down
    # X_mean_up
    X_mean_up <- lapply(seq_along(X_mean), function(i) replace(X_mean[[i]], i, X_up[i]) )
    # X_mean_down
    X_mean_down <- lapply(seq_along(X_mean), function(i) replace(X_mean[[i]], i, X_down[i]) )

    # adjust for categorical X (works also for zero categorical)
    for (i in X_categorical) {
      X_mean_up[[i]][[i]] <- ceiling(mean(X[,i]))
      X_mean_down[[i]][[i]] <- floor(mean(X[,i]))
    }
    # adjust for dummies (works also for zero dummies)
    for (i in X_dummy) {
      X_mean_up[[i]][[i]] <- max(X[,i])
      X_mean_down[[i]][[i]] <- min(X[,i])
    }

    # ----------------------------------------------------------------------------------- #

    #### now we do not need weights if we do not need inference (based on out of bag predictions)
    # forest prediction for X_mean_up (use new faster function particularly for ME)
    forest_pred_up <- predict_forest_preds_for_ME(forest, data, X_mean_up)
    # forest prediction for X_mean_down
    forest_pred_down <- predict_forest_preds_for_ME(forest, data, X_mean_down)

    # now subtract the predictions according to the ME formula
    forest_pred_diff_up_down <- mapply(function(x,y) mapply(function(x,y) x-y, x, y,  SIMPLIFY = F), forest_pred_up, forest_pred_down, SIMPLIFY = F)


    # subtract predictions according to formula to isolate categories
    forest_cond_means_0_last <- append(forest_pred_diff_up_down, list(rep(list(0), ncol(X)))) # append zero elemnt list
    forest_cond_means_0_first <- append(list(rep(list(0), ncol(X))), forest_pred_diff_up_down) # prepend zero element list

    # compute the scaling factor: X_up-X_down=2*X_sd
    scaling_factor <- as.list(X_up - X_down)
    # set scaling factor to zero for categorical and dummy variables
    for (i in X_categorical & X_dummy) {
      scaling_factor[[i]] <- 0
    }

    # now compute the differences for marginal effects
    marginal_effects_diff <- mapply(function(x,y) mapply(function(x,y) x-y, x, y, SIMPLIFY = F), forest_cond_means_0_last, forest_cond_means_0_first, SIMPLIFY = F)

    # scale marginal effects
    marginal_effects_scaled <- lapply(marginal_effects_diff, function(x) mapply(function(x,y) x/y, x, scaling_factor, SIMPLIFY = FALSE) )

    ## output for final marginal effects
    # coerce to a matrix
    marginal_effects <- sapply(marginal_effects_scaled, function(x) sapply(x, function(x) as.matrix(x)))
    # add names
    colnames(marginal_effects) <- sapply(categories, function(x) paste("Category", x, sep = " "))
    rownames(marginal_effects) <- colnames(X)

    # put marginal effects into results
    results <- marginal_effects

    # ---------------------------------------------------------------------------------------------------- #

  } else if (inference == TRUE & honesty == TRUE) {

    # ----------------------------------------------------------------------------------- #

    # take out honest indicator Ys
    Y_ind_honest <- lapply(data, function(x) as.matrix(x[, 1]))
    X_honest <- as.matrix(data[[1]][, -(1)])

    ## values for evaluation of the marginal effect - based on honest sample
    # share of SD to be used
    h_std <- 0.1

    # check if X is continuous or dummy or categorical
    X_type <- apply(X, 2, function(x) length(unique(x)))
    # now determine the type of X
    X_continuous <- which(X_type > 10) # define IDs of continuous Xs
    X_dummy <- which(X_type == 2) # define IDs of dummies
    X_categorical <- which(X_type > 2 & X_type <= 10)
    # additional check for constant variables which are nonsensical
    if (any(X_type == 1) | any(X_type == 0)) {
      stop("Some of the covariates are constant. This is non-sensical for evaluation of marginal effects. Programme terminated.")
    }

    # decide if the marginal effects should be computed at mean or at median
    if (eval=="mean") {
      # variable of interest: X_1 to X_last, ME at mean
      X_mean <- lapply(1:ncol(X), function(x) colMeans(X)) # set all Xs to their mean values (so many times as we have Xs)
    } else if (eval=="median") {
      # variable of interest: X_1 to X_last, ME at median
      X_mean <- lapply(1:ncol(X), function(x) apply(X, 2, median)) # set all Xs to their median values (so many times as we have Xs)
    } else {
      stop("Incorrect evaluation point. Programme terminated.")
    }

    # ----------------------------------------------------------------------------------- #

    # get SD of Xs
    X_sd <- apply(X, 2, sd)
    # create X_up (X_mean + 0.1 * X_sd)
    X_up <- X_mean[[1]] + h_std*X_sd
    # create X_down (X_mean - 0.1 * X_sd)
    X_down <- X_mean[[1]] - h_std*X_sd

    ## now check for the support of X
    # check X_max
    X_max <- apply(X, 2, max)
    # check X_min
    X_min <- apply(X, 2, min)
    # check if X_up is within the range X_min and X_max
    X_up <- (X_up < X_max) * X_up + (X_up >= X_max) * X_max
    X_up <- (X_up > X_min) * X_up + (X_up <= X_min) * (X_min + h_std * X_sd)
    # check if X_down is within the range X_min and X_max
    X_down <- (X_down > X_min) * X_down + (X_down <= X_min) * X_min
    X_down <- (X_down < X_max) * X_down + (X_down >= X_max) * (X_max - h_std * X_sd)
    # some additional checks

    # ----------------------------------------------------------------------------------- #

    ## now we need 2 datasets: one with X_up and second with X_down
    # X_mean_up
    X_mean_up <- lapply(seq_along(X_mean), function(i) replace(X_mean[[i]], i, X_up[i]) )
    # X_mean_down
    X_mean_down <- lapply(seq_along(X_mean), function(i) replace(X_mean[[i]], i, X_down[i]) )

    # adjust for categorical X (works also for zero categorical)
    for (i in X_categorical) {
      X_mean_up[[i]][[i]] <- ceiling(mean(X[,i]))
      X_mean_down[[i]][[i]] <- floor(mean(X[,i]))
    }
    # adjust for dummies (works also for zero dummies)
    for (i in X_dummy) {
      X_mean_up[[i]][[i]] <- max(X[,i])
      X_mean_down[[i]][[i]] <- min(X[,i])
    }

    # ----------------------------------------------------------------------------------- #

    # extract weights for desired Xs up: get weights from honest sample and predict weights for evaluation points from HONEST sample
    forest_weights_up <- predict_forest_weights_for_ME(forest, X_honest, X_mean_up)
    # extract weights for desired Xs down
    forest_weights_down <- predict_forest_weights_for_ME(forest, X_honest, X_mean_down)
    # now subtract the weights according to the ME formula
    forest_weights_diff_up_down <- mapply(function(x,y) mapply(function(x,y) x-y, x, y,  SIMPLIFY = F), forest_weights_up, forest_weights_down, SIMPLIFY = F)

    ## compute prerequisities for marginal effects
    # compute the conditional means: weights%*%y (predictions are based on honest sample)
    forest_cond_means <- mapply(function(x,y) lapply(x, function(x) x%*%y), forest_weights_diff_up_down, Y_ind_honest, SIMPLIFY = FALSE)
    # subtract conditional means according to formula to isolate categories
    forest_cond_means_0_last <- append(forest_cond_means, list(rep(list(0), ncol(X)))) # append zero elemnt list
    forest_cond_means_0_first <- append(list(rep(list(0), ncol(X))), forest_cond_means) # prepend zero element list
    # compute the scaling factor: X_up-X_down=2*X_sd
    scaling_factor <- as.list(X_up - X_down)
    # set scaling factor to zero for categorical and dummy variables
    for (i in X_categorical & X_dummy) {
      scaling_factor[[i]] <- 0
    }

    # now compute the differences for marginal effects
    marginal_effects_diff <- mapply(function(x,y) mapply(function(x,y) x-y, x, y, SIMPLIFY = F), forest_cond_means_0_last, forest_cond_means_0_first, SIMPLIFY = F)

    # scale marginal effects
    marginal_effects_scaled <- lapply(marginal_effects_diff, function(x) mapply(function(x,y) x/y, x, scaling_factor, SIMPLIFY = FALSE) )

    ## output for final marginal effects
    # coerce to a matrix
    marginal_effects <- sapply(marginal_effects_scaled, function(x) sapply(x, function(x) as.matrix(x)))
    # add names
    colnames(marginal_effects) <- sapply(categories, function(x) paste("Category", x, sep = " "))
    rownames(marginal_effects) <- colnames(X)

    # ----------------------------------------------------------------------------------- #

    ### variance for the marginal effects
    ## compute prerequisities for variance of honest marginal effects
    # scaling factor squared
    scaling_factor_squared <- lapply(scaling_factor, function(x) x^2)
    ## compute the conditional means (predictions): already have this as forest_cond_means
    # divide it by N to get the "mean"
    forest_cond_means_mean <- lapply(forest_cond_means, function(x) lapply(x, function(x) x/n_data))
    # calculate standard multiplication of weights and outcomes: honest_weights*y_ind_honest
    forest_multi <- mapply(function(x,y) lapply(x, function(x) t(x)*y), forest_weights_diff_up_down, Y_ind_honest, SIMPLIFY = FALSE)
    # subtract the mean from each obs i
    forest_multi_demeaned <- mapply(function(x,y) mapply(function(x,y) x-matrix(y, nrow = nrow(x)), x, y, SIMPLIFY = FALSE), forest_multi, forest_cond_means_mean, SIMPLIFY = F)

    ## now do the single variances for each category m
    # square the demeaned
    forest_multi_demeaned_sq <- lapply(forest_multi_demeaned, function(x) lapply(x, function(x) x^2))
    # sum all obs i together
    forest_multi_demeaned_sq_sum <- lapply(forest_multi_demeaned_sq, function(x) lapply(x, function(x) sum(x)))
    # divide by scaling factor
    forest_multi_demeaned_sq_sum_scaled <- lapply(forest_multi_demeaned_sq_sum, function(x) mapply(function(x,y) x/y, x, scaling_factor_squared, SIMPLIFY = FALSE) )
    # multiply by N/N-1 (normalize)
    forest_multi_demeaned_sq_sum_scaled_norm <- lapply(forest_multi_demeaned_sq_sum_scaled, function(x) lapply(x, function(x) x*(n_data/(n_data-1)) ))
    # put it into a shorter named object
    variance <- forest_multi_demeaned_sq_sum_scaled_norm
    ## single variances done

    # ----------------------------------------------------------------------------------- #

    ## now compute the covariances
    # multiply forest_var_multi_demeaned according to formula for covariance (shifted categories needed for computational convenience)
    forest_multi_demeaned_0_last <- append(forest_multi_demeaned, list(rep(list(matrix(0, ncol = ncol(forest_multi_demeaned[[1]][[1]]), nrow = nrow(forest_multi_demeaned[[1]][[1]]))), ncol(X)))) # append zero matrix list
    forest_multi_demeaned_0_first <- append(list(rep(list(matrix(0, ncol = ncol(forest_multi_demeaned[[1]][[1]]), nrow = nrow(forest_multi_demeaned[[1]][[1]]))), ncol(X))), forest_multi_demeaned) # prepend zero matrix list
    # compute the multiplication of category m with m-1 according to the covariance formula
    forest_multi_demeaned_cov <- mapply(function(x,y) mapply(function(x,y) x*y, x, y, SIMPLIFY = FALSE), forest_multi_demeaned_0_first, forest_multi_demeaned_0_last, SIMPLIFY = F)
    # sum all obs i together
    forest_multi_demeaned_cov_sum <- lapply(forest_multi_demeaned_cov, function(x) lapply(x, function(x) sum(x)))
    # divide by scaling factor
    forest_multi_demeaned_cov_sum_scaled <- lapply(forest_multi_demeaned_cov_sum, function(x) mapply(function(x,y) x/y, x, scaling_factor_squared, SIMPLIFY = FALSE) )
    # multiply by N/N-1 (normalize)
    forest_multi_demeaned_cov_sum_scaled_norm <- lapply(forest_multi_demeaned_cov_sum_scaled, function(x) lapply(x, function(x) x*(n_data/(n_data-1)) ))
    # multiply by 2
    forest_multi_demeaned_cov_sum_scaled_norm_mult2 <- lapply(forest_multi_demeaned_cov_sum_scaled_norm, function(x) lapply(x, function(x) x*2 ))
    # put it into a shorter named object
    covariance <- forest_multi_demeaned_cov_sum_scaled_norm_mult2
    ## covariances done

    # ----------------------------------------------------------------------------------- #

    ## put everything together according to the whole variance formula
    # shift variances accordingly for ease of next computations (covariance already has the desired format)
    variance_last <- append(variance, list(rep(list(0), ncol(X)))) # append zero element list
    variance_first <- append(list(rep(list(0), ncol(X))), variance) # prepend zero element list
    # put everything together according to formula: var_last + var_first - cov
    variance_marginal_effects_final <- mapply(function(x,y,z) mapply(function(x,y,z) x+y-z, x, y, z, SIMPLIFY = FALSE), variance_last, variance_first, covariance, SIMPLIFY = F)

    # ----------------------------------------------------------------------------------- #

    ## output for final variances of marginal effects
    # coerce to a matrix
    variance_marginal_effects <- sapply(variance_marginal_effects_final, function(x) sapply(x, function(x) as.matrix(x)))
    # add names
    colnames(variance_marginal_effects) <- sapply(categories, function(x) paste("Category", x, sep = " "))
    rownames(variance_marginal_effects) <- colnames(X)

    # ----------------------------------------------------------------------------------- #

    ## standard deviations
    # take square root of variance
    sd_marginal_effects <- sqrt(variance_marginal_effects)

    #### z scores and p values ####
    z_scores <- (marginal_effects)/(sd_marginal_effects)
    # control for dividing zero by zero
    z_scores[is.nan(z_scores)] = 0
    # p values
    p_values <- 2*pnorm(-abs(z_scores))

    # ----------------------------------------------------------------------------------- #

    # put everzthing into a list of results
    results <- list(marginal_effects, variance_marginal_effects, sd_marginal_effects, p_values)
    names(results) <- c("MarginalEffects", "Variances", "StandardErrors", "pValues")

    # ----------------------------------------------------------------------------------- #

  }

  # ----------------------------------------------------------------------------------- #

  # return results
  return(results)

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
predict_forest_preds_for_ME <- function(forest, data, pred_data) {

  # needed inputs for the function: forest - list of ranger forest objects
  #                                 data - list of n x n matrices of indicator data for orf
  #                                 pred_data - list of prediction data (X_mean_up/down)

  # ----------------------------------------------------------------------------------- #

  # get number of observations
  n_col <- ncol(data[[1]])-1 # number of X variables
  # stack pred_data together and make predictions at once
  pred_data <- do.call(rbind, pred_data)
  # prepare empty list
  forest_preds_up_together <- rep(list(NA), length(forest))
  forest_preds_up <- rep(list(rep(list(NA), n_col)), length(forest))

  # ----------------------------------------------------------------------------------- #

  # predict for desired Xs up
  for (forest_index in seq_along(forest)) {
    # make honest predictions
    forest_preds_up_together[[forest_index]] <- predict_honest(forest[[forest_index]], data[[forest_index]], pred_data)
    # now assign for each X into respective forest as list entries
    for (X_index in seq_along(pred_data[1,])) {

      forest_preds_up[[forest_index]][[X_index]] <- forest_preds_up_together[[forest_index]][X_index]

    }
  }

  # ----------------------------------------------------------------------------------- #

  # return predictions
  return(forest_preds_up)

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
#' @return list of weights
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
      # foreach(X_index=seq_along(pred_data)) %do% {

      # now do the same for your prediction data
      leaf_IDs_pred <- predict(forest[[forest_index]], as.data.frame(t(as.matrix(pred_data[[X_index]]))), type = "terminalNodes")$predictions
      # put leaf_IDs into a list
      leaf_IDs_pred <- lapply(seq_along(leaf_IDs_pred[1,]), function(i) leaf_IDs_pred[,i])

      # now average over the bootstraps, i.e. over trees to get final weights
      #forest_weights_pred_final <- Reduce("+", forest_weights_pred) / length(forest_weights_pred)
      forest_weights_pred_final <- pred_weights_C(leaf_IDs_pred, leaf_IDs, leaf_size)

      forest_weights_up[[forest_index]][[X_index]] <- forest_weights_pred_final

    }

  }

  # ----------------------------------------------------------------------------------- #

  # return weights
  return(forest_weights_up)

  # ----------------------------------------------------------------------------------- #

}

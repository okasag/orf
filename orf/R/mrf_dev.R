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


#' Predict MRF Variance
#'
#' predict variance of multinomial random forest predictions based on honest sample
#' splitting as described in Lechner (2018)
#'
#' @param honest_pred list of vectors of honest forest predictions
#' @param honest_weights list of n x n matrices of honest forest weights
#' @param Y_ind_honest list of vectors of 0-1 outcomes for the honest sample
#'
#' @return vector of MRF variances
pred_mrf_variance <- function(honest_pred, honest_weights, Y_ind_honest) {

  # needed inputs for the function: honest_pred - list of vectors of honest forest predictions
  #                                 honest_weights - list of n x n matrices of honest forest weights
  #                                 Y_ind_honest - list of vectors of 0-1 outcomes for the honest sample

  # ----------------------------------------------------------------------------------- #

  # get also categories
  categories <- seq(1:(length(Y_ind_honest)))

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
  # single variances done
  ## no covariances needed for MRF

  # ----------------------------------------------------------------------------------- #

  ## output for final variances
  # coerce to a matrix
  honest_var <- sapply(honest_variance, function(x) sapply(x, function(x) as.matrix(x)))

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


#' MRF margins
#'
#' estimate marginal effects for multinomial random forests as defined in
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
mrf_margins <- function(forest, data, honesty, inference){

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
  categories <- seq(1:(length(data)))
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
    ## now we have to normalize the predictions for MRF
    # get rowsum for each X
    forest_pred_up_sum <- lapply(seq_along(forest_pred_up[[1]]), function(k) sum(unlist(lapply(seq_along(forest_pred_up[[1]]), function(i) lapply(seq_along(forest_pred_up), function(j) unlist(forest_pred_up[[j]][[i]])))[[k]])))
    # normalize, i.e. divide each element by the sum of elements (for each X row)
    forest_pred_up_norm <- lapply(forest_pred_up, function(x) mapply(function(x,y) x/y, x, forest_pred_up_sum, SIMPLIFY = F))

    # forest prediction for X_mean_down
    forest_pred_down <- lapply(forest, function(x) lapply(X_mean_down, function(y) predict(x, data = t(as.matrix((y))))$predictions))
    ## now we have to normalize the predictions for MRF
    # get rowsum for each X
    forest_pred_down_sum <- lapply(seq_along(forest_pred_down[[1]]), function(k) sum(unlist(lapply(seq_along(forest_pred_down[[1]]), function(i) lapply(seq_along(forest_pred_down), function(j) unlist(forest_pred_down[[j]][[i]])))[[k]])))
    # normalize, i.e. divide each element by the sum of elements (for each X row)
    forest_pred_down_norm <- lapply(forest_pred_down, function(x) mapply(function(x,y) x/y, x, forest_pred_down_sum, SIMPLIFY = F))

    # now subtract the predictions according to the ME formula
    forest_pred_diff_up_down <- mapply(function(x,y) mapply(function(x,y) x-y, x, y,  SIMPLIFY = F), forest_pred_up_norm, forest_pred_down_norm, SIMPLIFY = F)

    # compute the scaling factor: X_up-X_down=2*X_sd
    scaling_factor <- as.list(X_up - X_down)
    # set scaling factor to zero for categorical and dummy variables
    for (i in X_categorical & X_dummy) {
      scaling_factor[[i]] <- 0
    }

    # scale marginal effects
    marginal_effects_scaled <- lapply(forest_pred_diff_up_down, function(x) mapply(function(x,y) x/y, x, scaling_factor, SIMPLIFY = FALSE) )

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
    ## now we have to normalize the predictions for MRF
    # get rowsum for each X
    forest_pred_up_sum <- lapply(seq_along(forest_pred_up[[1]]), function(k) sum(unlist(lapply(seq_along(forest_pred_up[[1]]), function(i) lapply(seq_along(forest_pred_up), function(j) unlist(forest_pred_up[[j]][[i]])))[[k]])))
    # normalize, i.e. divide each element by the sum of elements (for each X row)
    forest_pred_up_norm <- lapply(forest_pred_up, function(x) mapply(function(x,y) x/y, x, forest_pred_up_sum, SIMPLIFY = F))

    # forest prediction for X_mean_down
    forest_pred_down <- predict_forest_preds_for_ME(forest, data, X_mean_down)
    ## now we have to normalize the predictions for MRF
    # get rowsum for each X
    forest_pred_down_sum <- lapply(seq_along(forest_pred_down[[1]]), function(k) sum(unlist(lapply(seq_along(forest_pred_down[[1]]), function(i) lapply(seq_along(forest_pred_down), function(j) unlist(forest_pred_down[[j]][[i]])))[[k]])))
    # normalize, i.e. divide each element by the sum of elements (for each X row)
    forest_pred_down_norm <- lapply(forest_pred_down, function(x) mapply(function(x,y) x/y, x, forest_pred_down_sum, SIMPLIFY = F))

    # now subtract the predictions according to the ME formula
    forest_pred_diff_up_down <- mapply(function(x,y) mapply(function(x,y) x-y, x, y,  SIMPLIFY = F), forest_pred_up_norm, forest_pred_down_norm, SIMPLIFY = F)

    # compute the scaling factor: X_up-X_down=2*X_sd
    scaling_factor <- as.list(X_up - X_down)
    # set scaling factor to zero for categorical and dummy variables
    for (i in X_categorical & X_dummy) {
      scaling_factor[[i]] <- 0
    }

    # scale marginal effects
    marginal_effects_scaled <- lapply(forest_pred_diff_up_down, function(x) mapply(function(x,y) x/y, x, scaling_factor, SIMPLIFY = FALSE) )

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
    forest_weights_up <- predict_forest_weights_for_ME(forest, X, X_mean_up)
    # extract weights for desired Xs down
    forest_weights_down <- predict_forest_weights_for_ME(forest, X, X_mean_down)
    # now subtract the weights according to the ME formula
    forest_weights_diff_up_down <- mapply(function(x,y) mapply(function(x,y) x-y, x, y,  SIMPLIFY = F), forest_weights_up, forest_weights_down, SIMPLIFY = F)

    ## compute prerequisities for marginal effects
    # compute the conditional means: i.e. predictions, for MRF nromalize pred_up and pred_down first
    forest_pred_up <- mapply(function(x,y) lapply(x, function(x) x%*%y), forest_weights_up, Y_ind_honest, SIMPLIFY = FALSE)
    # get rowsum for each X
    forest_pred_up_sum <- lapply(seq_along(forest_pred_up[[1]]), function(k) sum(unlist(lapply(seq_along(forest_pred_up[[1]]), function(i) lapply(seq_along(forest_pred_up), function(j) unlist(forest_pred_up[[j]][[i]])))[[k]])))
    # normalize, i.e. divide each element by the sum of elements (for each X row)
    forest_pred_up_norm <- lapply(forest_pred_up, function(x) mapply(function(x,y) x/y, x, forest_pred_up_sum, SIMPLIFY = F))

    # compute the conditional means: i.e. predictions, for MRF nromalize pred_up and pred_down first
    forest_pred_down <- mapply(function(x,y) lapply(x, function(x) x%*%y), forest_weights_down, Y_ind_honest, SIMPLIFY = FALSE)
    # get rowsum for each X
    forest_pred_down_sum <- lapply(seq_along(forest_pred_down[[1]]), function(k) sum(unlist(lapply(seq_along(forest_pred_down[[1]]), function(i) lapply(seq_along(forest_pred_down), function(j) unlist(forest_pred_down[[j]][[i]])))[[k]])))
    # normalize, i.e. divide each element by the sum of elements (for each X row)
    forest_pred_down_norm <- lapply(forest_pred_down, function(x) mapply(function(x,y) x/y, x, forest_pred_down_sum, SIMPLIFY = F))

    # now subtract the predictions according to the ME formula
    forest_pred_diff_up_down <- mapply(function(x,y) mapply(function(x,y) x-y, x, y,  SIMPLIFY = F), forest_pred_up_norm, forest_pred_down_norm, SIMPLIFY = F)

    # compute the scaling factor: X_up-X_down=2*X_sd
    scaling_factor <- as.list(X_up - X_down)
    # set scaling factor to zero for categorical and dummy variables
    for (i in X_categorical & X_dummy) {
      scaling_factor[[i]] <- 0
    }

    # scale marginal effects
    marginal_effects_scaled <- lapply(forest_pred_diff_up_down, function(x) mapply(function(x,y) x/y, x, scaling_factor, SIMPLIFY = FALSE) )

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
    forest_cond_means_mean <- lapply(forest_pred_diff_up_down, function(x) lapply(x, function(x) x/n_data))
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
    # single variances done
    ## covariances not needed for MRF

    # ----------------------------------------------------------------------------------- #

    ## output for final variances of marginal effects
    # coerce to a matrix
    variance_marginal_effects <- sapply(variance, function(x) sapply(x, function(x) as.matrix(x)))
    # add names
    colnames(variance_marginal_effects) <- sapply(categories, function(x) paste("Category", x, sep = " "))
    rownames(variance_marginal_effects) <- colnames(X)

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


#' Predict MRF Margins
#'
#' estimate marginal effects for multinomial random forests as defined in
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
pred_mrf_margins <- function(forest, data, newdata, honesty, inference){

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
  categories <- seq(1:(length(data)))
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
    ## now we have to normalize the predictions for MRF
    # get rowsum for each X
    forest_pred_up_sum <- lapply(seq_along(forest_pred_up[[1]]), function(k) sum(unlist(lapply(seq_along(forest_pred_up[[1]]), function(i) lapply(seq_along(forest_pred_up), function(j) unlist(forest_pred_up[[j]][[i]])))[[k]])))
    # normalize, i.e. divide each element by the sum of elements (for each X row)
    forest_pred_up_norm <- lapply(forest_pred_up, function(x) mapply(function(x,y) x/y, x, forest_pred_up_sum, SIMPLIFY = F))

    # forest prediction for X_mean_down
    forest_pred_down <- lapply(forest, function(x) lapply(X_mean_down, function(y) predict(x, data = t(as.matrix((y))))$predictions))
    ## now we have to normalize the predictions for MRF
    # get rowsum for each X
    forest_pred_down_sum <- lapply(seq_along(forest_pred_down[[1]]), function(k) sum(unlist(lapply(seq_along(forest_pred_down[[1]]), function(i) lapply(seq_along(forest_pred_down), function(j) unlist(forest_pred_down[[j]][[i]])))[[k]])))
    # normalize, i.e. divide each element by the sum of elements (for each X row)
    forest_pred_down_norm <- lapply(forest_pred_down, function(x) mapply(function(x,y) x/y, x, forest_pred_down_sum, SIMPLIFY = F))

    # now subtract the predictions according to the ME formula
    forest_pred_diff_up_down <- mapply(function(x,y) mapply(function(x,y) x-y, x, y,  SIMPLIFY = F), forest_pred_up_norm, forest_pred_down_norm, SIMPLIFY = F)

    # compute the scaling factor: X_up-X_down=2*X_sd
    scaling_factor <- as.list(X_up - X_down)
    # set scaling factor to zero for categorical and dummy variables
    for (i in X_categorical & X_dummy) {
      scaling_factor[[i]] <- 0
    }

    # scale marginal effects
    marginal_effects_scaled <- lapply(forest_pred_diff_up_down, function(x) mapply(function(x,y) x/y, x, scaling_factor, SIMPLIFY = FALSE) )

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
    ## now we have to normalize the predictions for MRF
    # get rowsum for each X
    forest_pred_up_sum <- lapply(seq_along(forest_pred_up[[1]]), function(k) sum(unlist(lapply(seq_along(forest_pred_up[[1]]), function(i) lapply(seq_along(forest_pred_up), function(j) unlist(forest_pred_up[[j]][[i]])))[[k]])))
    # normalize, i.e. divide each element by the sum of elements (for each X row)
    forest_pred_up_norm <- lapply(forest_pred_up, function(x) mapply(function(x,y) x/y, x, forest_pred_up_sum, SIMPLIFY = F))

    # forest prediction for X_mean_down
    forest_pred_down <- predict_forest_preds_for_ME(forest, data, X_mean_down)
    ## now we have to normalize the predictions for MRF
    # get rowsum for each X
    forest_pred_down_sum <- lapply(seq_along(forest_pred_down[[1]]), function(k) sum(unlist(lapply(seq_along(forest_pred_down[[1]]), function(i) lapply(seq_along(forest_pred_down), function(j) unlist(forest_pred_down[[j]][[i]])))[[k]])))
    # normalize, i.e. divide each element by the sum of elements (for each X row)
    forest_pred_down_norm <- lapply(forest_pred_down, function(x) mapply(function(x,y) x/y, x, forest_pred_down_sum, SIMPLIFY = F))

    # now subtract the predictions according to the ME formula
    forest_pred_diff_up_down <- mapply(function(x,y) mapply(function(x,y) x-y, x, y,  SIMPLIFY = F), forest_pred_up_norm, forest_pred_down_norm, SIMPLIFY = F)

    # compute the scaling factor: X_up-X_down=2*X_sd
    scaling_factor <- as.list(X_up - X_down)
    # set scaling factor to zero for categorical and dummy variables
    for (i in X_categorical & X_dummy) {
      scaling_factor[[i]] <- 0
    }

    # scale marginal effects
    marginal_effects_scaled <- lapply(forest_pred_diff_up_down, function(x) mapply(function(x,y) x/y, x, scaling_factor, SIMPLIFY = FALSE) )

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
    # compute the conditional means: i.e. predictions, for MRF nromalize pred_up and pred_down first
    forest_pred_up <- mapply(function(x,y) lapply(x, function(x) x%*%y), forest_weights_up, Y_ind_honest, SIMPLIFY = FALSE)
    # get rowsum for each X
    forest_pred_up_sum <- lapply(seq_along(forest_pred_up[[1]]), function(k) sum(unlist(lapply(seq_along(forest_pred_up[[1]]), function(i) lapply(seq_along(forest_pred_up), function(j) unlist(forest_pred_up[[j]][[i]])))[[k]])))
    # normalize, i.e. divide each element by the sum of elements (for each X row)
    forest_pred_up_norm <- lapply(forest_pred_up, function(x) mapply(function(x,y) x/y, x, forest_pred_up_sum, SIMPLIFY = F))

    # compute the conditional means: i.e. predictions, for MRF nromalize pred_up and pred_down first
    forest_pred_down <- mapply(function(x,y) lapply(x, function(x) x%*%y), forest_weights_down, Y_ind_honest, SIMPLIFY = FALSE)
    # get rowsum for each X
    forest_pred_down_sum <- lapply(seq_along(forest_pred_down[[1]]), function(k) sum(unlist(lapply(seq_along(forest_pred_down[[1]]), function(i) lapply(seq_along(forest_pred_down), function(j) unlist(forest_pred_down[[j]][[i]])))[[k]])))
    # normalize, i.e. divide each element by the sum of elements (for each X row)
    forest_pred_down_norm <- lapply(forest_pred_down, function(x) mapply(function(x,y) x/y, x, forest_pred_down_sum, SIMPLIFY = F))

    # now subtract the predictions according to the ME formula
    forest_pred_diff_up_down <- mapply(function(x,y) mapply(function(x,y) x-y, x, y,  SIMPLIFY = F), forest_pred_up_norm, forest_pred_down_norm, SIMPLIFY = F)

    # compute the scaling factor: X_up-X_down=2*X_sd
    scaling_factor <- as.list(X_up - X_down)
    # set scaling factor to zero for categorical and dummy variables
    for (i in X_categorical & X_dummy) {
      scaling_factor[[i]] <- 0
    }

    # scale marginal effects
    marginal_effects_scaled <- lapply(forest_pred_diff_up_down, function(x) mapply(function(x,y) x/y, x, scaling_factor, SIMPLIFY = FALSE) )

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
    forest_cond_means_mean <- lapply(forest_pred_diff_up_down, function(x) lapply(x, function(x) x/n_data))
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
    # single variances done
    ## covariances not needed for MRF

    # ----------------------------------------------------------------------------------- #

    ## output for final variances of marginal effects
    # coerce to a matrix
    variance_marginal_effects <- sapply(variance, function(x) sapply(x, function(x) as.matrix(x)))
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

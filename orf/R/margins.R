#' margins
#'
#'estimate marginal effects of various discrete choice models based on the
#'random forests algorithm. Applicable classes are \code{orf}, \code{mrf} and
#'\code{brf}.
#'
#' @param forest trained forest object of class \code{orf}/\code{mrf}/\code{brf}
#' @param eval string defining evaluation point for marginal effects. These can be one of "mean", "atmean", or "atmedian"
#' @param newdata matrix of new Xs (currently not supported)
#'
#' @importFrom stats predict median pnorm sd
#' @import ranger
#'
#' @export
margins <- function(forest, eval, newdata) {

  # needed inputs for the function: forest - trained forest object of class orf/mrf/brf
  #                                 eval - string defining evaluation point for marginal effects
  #                                 newdata - matrix of new Xs
  # ----------------------------------------------------------------------------------- #

  ### decide if prediction or in sample marginal effects should be evaluated
  if (is.null(newdata)) {

    # if no newdata supplied, estimate in sample marginal effects
    if (forest$forestInfo$inputs$honesty == FALSE) {

      data <- forest$forestInfo$trainData # take in-sample data

    } else if (forest$forestInfo$inputs$honesty == TRUE) {

      data <- forest$forestInfo$honestData

    }



  } else {

    # check if newdata is compatible with train data
    if (ncol(newdata) != ncol(forest$forestInfo$trainData)) {

      stop("newdata is not compatible with training data. Programme terminated.")

    } else {

      data = newdata

    }

  }
  ### data checks done

  # ----------------------------------------------------------------------------------- #

  ### data preparation and checks
  # get number of observations
  n_data <- as.numeric(nrow(data))
  # get categories
  categories <- forest$forestInfo$categories
  # get X as matrix
  X <- as.matrix(data[, -1])
  # get Y as matrix
  Y <- as.matrix(data[, 1])
  # create indicator variables (outcomes)
  Y_ind <- lapply(categories[1:length(categories)-1], function(x) ifelse((Y <= x), 1, 0))
  # create datasets with indicator outcomes
  data_ind <- lapply(Y_ind, function(x) as.data.frame(cbind(as.matrix(unlist(x)), X)))

  # ----------------------------------------------------------------------------------- #

  ### marginal effects preparation
  # share of SD to be used
  h_std <- 1
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

  # ----------------------------------------------------------------------------------- #

  ### check the evaluation point
  if (eval == "atmean") {
    # variable of interest: X_1 to X_last, ME at mean
    X_mean <- lapply(1:ncol(X), function(x) t(as.matrix(colMeans(X)))) # set all Xs to their mean values (so many times as we have Xs)
  } else if (eval == "atmedian") {
    # variable of interest: X_1 to X_last, ME at median
    X_mean <- lapply(1:ncol(X), function(x) t(as.matrix(apply(X, 2, median)))) # set all Xs to their median values (so many times as we have Xs)
  } else if (eval == "mean") {
    # # variable of interest: X_1 to X_last, mean ME
    X_mean <- lapply(1:ncol(X), function(x) X) # set all Xs to their exact values (so many times as we have Xs)
  } else {
    stop("Incorrect evaluation point. Programme terminated.")
  }

  # ----------------------------------------------------------------------------------- #

  ### get data needed for evaluation of ME
  # get number of evaluation points
  X_rows <- nrow(X_mean[[1]])
  # get number of Xs
  X_cols <- ncol(X_mean[[1]])
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
  # X_mean_up all (continous)
  X_mean_up <- X_mean
  X_mean_down <- X_mean

  # replace values accordingly
  for (i in 1:X_cols) {
    X_mean_up[[i]][, i] <- X_up[, i]
    X_mean_down[[i]][, i] <- X_down[, i]
  }

  # adjust for categorical X (works also for zero categorical)
  for (i in X_categorical) {
    X_mean_up[[i]][, i] <- ceiling(mean(X[, i]))
    X_mean_down[[i]][, i] <- floor(mean(X[, i]))
  }

  # adjust for dummies (works also for zero dummies)
  for (i in X_dummy) {
    X_mean_up[[i]][, i] <- max(X[, i])
    X_mean_down[[i]][, i] <- min(X[, i])
  }

  # ----------------------------------------------------------------------------------- #

  ### check honesty and inference
  if (forest$forestInfo$inputs$honesty == FALSE & forest$forestInfo$inputs$inference == FALSE) {

    #### now we do not need weights if we do not need inference (based on out of bag predictions)
    # forest prediction for X_mean_up
    forest_pred_up <- lapply(forest$trainForests, function(x) lapply(X_mean_up, function(y) predict(x, data = y)$predictions))
    # forest prediction for X_mean_down
    forest_pred_down <- lapply(forest$trainForests, function(x) lapply(X_mean_down, function(y) predict(x, data = y)$predictions))

  } else if (forest$forestInfo$inputs$honesty == TRUE & forest$forestInfo$inputs$inference == FALSE) {

    # do honest predictions
    # forest prediction for X_mean_up (use new faster function particularly for ME)
    forest_pred_up <- predict_forest_preds_for_ME(forest$trainForests, data_ind, X_mean_up)
    # forest prediction for X_mean_down
    forest_pred_down <- predict_forest_preds_for_ME(forest$trainForests, data_ind, X_mean_down)

  } else if (forest$forestInfo$inputs$honesty == TRUE & forest$forestInfo$inputs$inference == TRUE) {

    # do honest predictions with weight based inference
    # extract weights for desired Xs up: get weights from honest sample and predict weights for evaluation points from HONEST sample
    forest_weights_up <- predict_forest_weights_for_ME(forest$trainForests, X, X_mean_up)
    # extract weights for desired Xs down
    forest_weights_down <- predict_forest_weights_for_ME(forest$trainForests, X, X_mean_down)

    ## compute predictions based on weights
    # forest prediction for X_mean_up
    forest_pred_up <- mapply(function(x,y) lapply(x, function(x) as.numeric(x%*%y)), forest_weights_up, Y_ind, SIMPLIFY = FALSE)
    # forest prediction for X_mean_down
    forest_pred_down <- mapply(function(x,y) lapply(x, function(x) as.numeric(x%*%y)), forest_weights_down, Y_ind, SIMPLIFY = FALSE)

  }

  # ----------------------------------------------------------------------------------- #

  ### form ORF predictions
  # prepare up
  forest_pred_up_1 <- append(forest_pred_up, list(rep(list(rep(1, X_rows)), X_cols))) # append 1
  forest_pred_up_0 <- append(list(rep(list(rep(0, X_rows)), X_cols)), forest_pred_up) # prepend 0
  # isolate predictions
  forest_pred_up <- mapply(function(x,y) mapply(function(x,y) x-y, x, y,  SIMPLIFY = F), forest_pred_up_1, forest_pred_up_0, SIMPLIFY = F)
  # avoid negative predictions
  forest_pred_up <- lapply(forest_pred_up, function(x) lapply(x, function(x) ifelse((x < 0), 0, x)))
  # normalize predictions
  forest_pred_up_rowsum <- lapply(seq_along(forest_pred_up[[1]]), function(i) rowSums(matrix(sapply(forest_pred_up, "[[", i), ncol = length(categories), nrow = nrow(X_mean_up[[1]])))) # build rowsums with respect to categories
  forest_pred_up <- lapply(forest_pred_up, function(x) mapply(function(x,y) x/y, x, forest_pred_up_rowsum, SIMPLIFY = FALSE)) # normalize to sum up to 1 (very rare but just to be sure)

  # prepare down
  forest_pred_down_1 <- append(forest_pred_down, list(rep(list(rep(1, X_rows)), X_cols))) # append 1
  forest_pred_down_0 <- append(list(rep(list(rep(0, X_rows)), X_cols)), forest_pred_down) # prepend 0
  # isolate predictions
  forest_pred_down <- mapply(function(x,y) mapply(function(x,y) x-y, x, y,  SIMPLIFY = F), forest_pred_down_1, forest_pred_down_0, SIMPLIFY = F)
  # avoid negative predictions
  forest_pred_down <- lapply(forest_pred_down, function(x) lapply(x, function(x) ifelse((x < 0), 0, x)))
  # normalize predictions
  forest_pred_down_rowsum <- lapply(seq_along(forest_pred_down[[1]]), function(i) rowSums(matrix(sapply(forest_pred_down, "[[", i), ncol = length(categories), nrow = nrow(X_mean_down[[1]])))) # build rowsums with respect to categories
  forest_pred_down <- lapply(forest_pred_down, function(x) mapply(function(x,y) x/y, x, forest_pred_down_rowsum, SIMPLIFY = FALSE)) # normalize to sum up to 1 (very rare but just to be sure)

  # ----------------------------------------------------------------------------------- #

  ### now subtract the predictions according to the ME formula
  forest_pred_diff_up_down <- mapply(function(x,y) mapply(function(x,y) x-y, x, y,  SIMPLIFY = F), forest_pred_up, forest_pred_down, SIMPLIFY = F)
  # compute the scaling factor: X_up-X_down=2*X_sd
  scaling_factor <- split((X_up - X_down), 1:X_cols) # save it as separate list vectors
  # set scaling factor to zero for categorical and dummy variables
  for (i in X_categorical & X_dummy) {
    scaling_factor[[i]] <- rep(0, X_rows)
  }
  # scale the differences to get marginal effects
  marginal_effects_scaled <- lapply(forest_pred_diff_up_down, function(x) mapply(function(x,y) x/y, x, scaling_factor, SIMPLIFY = FALSE) )
  # take means fo marginal effects for each X (for atmean or atmedian this doesnt change anything)
  marginal_effects_mean <- lapply(marginal_effects_scaled, function(x) lapply(x, function(x) mean(x)))

  # ----------------------------------------------------------------------------------- #

  # coerce to a matrix
  marginal_effects <- sapply(marginal_effects_mean, function(x) sapply(x, function(x) as.matrix(x)))
  # add names
  colnames(marginal_effects) <- sapply(categories, function(x) paste("Category", x, sep = " "))
  rownames(marginal_effects) <- colnames(X)

  # ----------------------------------------------------------------------------------- #

  # put marginal effects into results
  results <- marginal_effects

  # ----------------------------------------------------------------------------------- #

  if (forest$forestInfo$inputs$inference == TRUE) {

    ### variance for the marginal effects
    ## compute prerequisities for variance of honest marginal effects
    # scaling factor squared
    scaling_factor_squared <- lapply(scaling_factor, function(x) x^2)

    # now subtract the weights according to the ME formula
    forest_weights_diff_up_down <- mapply(function(x,y) mapply(function(x,y) x-y, x, y,  SIMPLIFY = F), forest_weights_up, forest_weights_down, SIMPLIFY = F)

    # compute the conditional means: 1/N(weights%*%y) (predictions are based on honest sample)
    forest_cond_means <- mapply(function(x,y) lapply(x, function(x) (x%*%y)/nrow(Y_ind[[1]])), forest_weights_diff_up_down, Y_ind, SIMPLIFY = FALSE)
    # calculate standard multiplication of weights and outcomes: honest_weights*y_ind_honest
    if (eval == "mean") {
      forest_multi <- mapply(function(x,y) lapply(x, function(x) apply(t(x), 2, function(x) x*y)), forest_weights_diff_up_down, Y_ind, SIMPLIFY = FALSE)
      # subtract the mean from each obs i
      forest_multi_demeaned <- mapply(function(x,y) mapply(function(x,y) {
        # demean each and every observation
        sapply(seq_along(y[, 1]), function(i) x[, i] - matrix(y[i, 1], nrow = nrow(x))) # sapply as we want to keep matrix structure

        }, x, y, SIMPLIFY = FALSE), forest_multi, forest_cond_means, SIMPLIFY = F)

      ## now do the single variances for each category m
      # square the demeaned and sum it and normalize (colSums as eahc column represents one observations)
      forest_multi_demeaned_sq_sum_norm <- lapply(forest_multi_demeaned, function(x) lapply(x, function(x) (colSums(x^2))*(nrow(Y_ind[[1]])/(nrow(Y_ind[[1]])-1))))

      ## now compute the prerequisites for covariances
      # multiply forest_var_multi_demeaned according to formula for covariance (shifted categories needed for computational convenience)
      forest_multi_demeaned_0_last <- append(forest_multi_demeaned, list(rep(list(matrix(0, ncol = ncol(forest_multi_demeaned[[1]][[1]]), nrow = nrow(forest_multi_demeaned[[1]][[1]]))), ncol(X)))) # append zero matrix list
      forest_multi_demeaned_0_first <- append(list(rep(list(matrix(0, ncol = ncol(forest_multi_demeaned[[1]][[1]]), nrow = nrow(forest_multi_demeaned[[1]][[1]]))), ncol(X))), forest_multi_demeaned) # prepend zero matrix list
      # compute the multiplication of category m with m-1 according to the covariance formula (sum, normalize and multiply by 2)
      forest_multi_demeaned_cov_sum_norm_mult2 <- mapply(function(x,y) mapply(function(x,y) (colSums(x*y))*(nrow(Y_ind[[1]])/(nrow(Y_ind[[1]])-1))*2, x, y, SIMPLIFY = FALSE), forest_multi_demeaned_0_first, forest_multi_demeaned_0_last, SIMPLIFY = F)

    } else {
      forest_multi <- mapply(function(x,y) lapply(x, function(x) t(x)*y), forest_weights_diff_up_down, Y_ind, SIMPLIFY = FALSE)
      # subtract the mean from each obs i
      forest_multi_demeaned <- mapply(function(x,y) mapply(function(x,y) x-matrix(y, nrow = nrow(x)), x, y, SIMPLIFY = FALSE), forest_multi, forest_cond_means, SIMPLIFY = F)

      ## now do the single variances for each category m
      # square the demeaned and sum it and normalize
      forest_multi_demeaned_sq_sum_norm <- lapply(forest_multi_demeaned, function(x) lapply(x, function(x) (sum(x^2))*(nrow(Y_ind[[1]])/(nrow(Y_ind[[1]])-1))))

      ## now compute the prerequisites for covariances
      # multiply forest_var_multi_demeaned according to formula for covariance (shifted categories needed for computational convenience)
      forest_multi_demeaned_0_last <- append(forest_multi_demeaned, list(rep(list(matrix(0, ncol = ncol(forest_multi_demeaned[[1]][[1]]), nrow = nrow(forest_multi_demeaned[[1]][[1]]))), ncol(X)))) # append zero matrix list
      forest_multi_demeaned_0_first <- append(list(rep(list(matrix(0, ncol = ncol(forest_multi_demeaned[[1]][[1]]), nrow = nrow(forest_multi_demeaned[[1]][[1]]))), ncol(X))), forest_multi_demeaned) # prepend zero matrix list
      # compute the multiplication of category m with m-1 according to the covariance formula (sum, normalize and multiply by 2)
      forest_multi_demeaned_cov_sum_norm_mult2 <- mapply(function(x,y) mapply(function(x,y) (sum(x*y))*(nrow(Y_ind[[1]])/(nrow(Y_ind[[1]])-1))*2, x, y, SIMPLIFY = FALSE), forest_multi_demeaned_0_first, forest_multi_demeaned_0_last, SIMPLIFY = F)
    }

    # divide by scaling factor to get the variance
    variance <- lapply(forest_multi_demeaned_sq_sum_norm, function(x) mapply(function(x,y) x/y, x, scaling_factor_squared, SIMPLIFY = FALSE) )
    ## single variances done

    # ----------------------------------------------------------------------------------- #

    # divide by scaling factor
    covariance <- lapply(forest_multi_demeaned_cov_sum_norm_mult2, function(x) mapply(function(x,y) x/y, x, scaling_factor_squared, SIMPLIFY = FALSE) )
    ## covariances done

    # ----------------------------------------------------------------------------------- #

    ## put everything together according to the whole variance formula
    # shift variances accordingly for ease of next computations (covariance already has the desired format)
    variance_last <- append(variance, list(rep(list(rep(0, length(variance[[1]][[1]]))), X_cols))) # append zero element list
    variance_first <- append(list(rep(list(rep(0, length(variance[[1]][[1]]))), X_cols)), variance) # prepend zero element list
    # put everything together according to formula: var_last + var_first - cov (mean doesnt change anything for atmean and atmedian evaluations)
    variance_me <- mapply(function(x,y,z) mapply(function(x,y,z) mean(x+y-z), x, y, z, SIMPLIFY = FALSE), variance_last, variance_first, covariance, SIMPLIFY = F)

    # ----------------------------------------------------------------------------------- #

    ## output for final variances of marginal effects
    # coerce to a matrix
    variance_me <- sapply(variance_me, function(x) sapply(x, function(x) as.matrix(x)))
    # add names
    colnames(variance_me) <- sapply(categories, function(x) paste("Category", x, sep = " "))
    rownames(variance_me) <- colnames(X)

    # ----------------------------------------------------------------------------------- #

    ## standard deviations
    # take square root of variance
    sd_me <- sqrt(variance_me)

    #### z scores and p values ####
    t_value <- (marginal_effects)/(sd_me)
    # control for dividing zero by zero
    t_value[is.nan(t_value)] = 0
    # p values
    p_values <- 2*pnorm(-abs(t_value))

    # ----------------------------------------------------------------------------------- #

    # put everything into a list of results
    results <- list(marginal_effects, variance_me, sd_me, t_value, p_values)
    names(results) <- c("MarginalEffects", "Variances", "StandardErrors", "tValues", "pValues")

    # ----------------------------------------------------------------------------------- #

  }

  # return results
  return(results)

}

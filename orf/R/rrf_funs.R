#' Get Variance
#'
#' get variance of regression random forest predictions based on honest sample
#' splitting as described in Lechner (2018)
#'
#' @param honest_pred vector of honest forest predictions
#' @param honest_weights n x n matrix of honest forest weights (based on train)
#' @param honest_outcomes vector of real outcomes for the honest sample
#'
#' @return vector of variances

get_variance <- function(honest_pred, honest_weights, honest_outcomes) {

  # needed inputs for the function: honest_pred - vector of honest forest predictions
  #                                 honest_weights - n x n matrix of honest forest weights (based on train)
  #                                 outcomes - vector of real outcomes for the honest sample

  # ----------------------------------------------------------------------------------- #

  ## compute the conditional means (predictions): already have this as honest_pred
  # compute mean (divide each predicted value by the number of observations)
  honest_pred_mean <- honest_pred/length(honest_outcomes)
  # calculate standard multiplication of weights and outcomes: honest_weights*y_ind_honest
  forest_multi <- lapply(seq_along(honest_weights[, 1]),  function(i) honest_weights[i, ] * honest_outcomes)
  # subtract the mean from each obs i
  forest_multi_demeaned <- mapply(function(x,y) x - y, forest_multi, honest_pred_mean, SIMPLIFY = FALSE)
  # square the demeaned
  forest_multi_demeaned_sq <- lapply(forest_multi_demeaned, function(x) x^2)
  # sum all obs i together
  forest_multi_demeaned_sq_sum <- lapply(forest_multi_demeaned_sq, function(x) sum(x))
  # multiply by N/N-1 (normalize)
  forest_multi_demeaned_sq_sum_norm <- lapply(forest_multi_demeaned_sq_sum, function(x) x*(length(honest_outcomes)/(length(honest_outcomes)-1)) )
  # put it into a shorter named object
  variance <- unlist(forest_multi_demeaned_sq_sum_norm)
  ## single variances done

  # ----------------------------------------------------------------------------------- #

  # put the variances into output
  output <- as.numeric(variance)
  # return the output
  return(output)

  # ----------------------------------------------------------------------------------- #

}

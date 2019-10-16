#' margins
#'
#' margins is an S3 generic with methods to estimate marginal effects
#' of various discrete choice models based on the random forests algorithm.
#' Applicable classes are \code{orf} and \code{mrf}.
#'
#' @param forest estimated forest object of class \code{orf} or \code{mrf}
#' @param eval string defining evaluation point for marginal effects. These can be one of "mean", "atmean", or "atmedian". (Default is "mean")
#' @param inference logical, if TRUE inference on marginal effects will be conducted (default is inherited from the orf object)
#' @param window numeric, share of standard deviation of X to be used for evaluation of the marginal effect (default is 0.1)
#' @param newdata matrix of new Xs for which marginal effects will be computed
#'
#' @export
margins <- function(forest, eval = NULL, inference = NULL, window = NULL, newdata = NULL) UseMethod("margins")


#' margins.default
#'
#' margins.default is a default for S3 generic with methods to estimate marginal effects
#' of various discrete choice models based on the random forests algorithm.
#' Applicable classes are \code{orf} or \code{mrf} .#'
#' @param forest estimated forest object of class \code{orf} or \code{mrf}
#' @param eval string defining evaluation point for marginal effects. These can be one of "mean", "atmean", or "atmedian". (Default is "mean")
#' @param inference logical, if TRUE inference on marginal effects will be conducted (default is inherited from the orf object)
#' @param window numeric, share of standard deviation of X to be used for evaluation of the marginal effect (default is 0.1)
#' @param newdata matrix of new Xs for which marginal effects will be computed
#'
#' @export
margins.default <- function(forest, eval = NULL, inference = NULL, window = NULL, newdata = NULL) {

  warning(paste("margins does not know how to handle object of class ",
                class(forest),
                ". The supported classes are one of the following: orf/mrf."))

}


#' margins.orf
#'
#' estimate marginal effects of the ordered random forest
#'
#' @param forest trained ordered random forest object of class \code{orf}
#' @param eval string defining evaluation point for marginal effects. These can be one of "mean", "atmean", or "atmedian" (default is "mean")
#' @param inference logical, if TRUE inference on marginal effects will be conducted (default is inherited from the orf object)
#' @param window numeric, share of standard deviation of X to be used for evaluation of the marginal effect (default is 0.1)
#' @param newdata matrix of new Xs for which marginal effects will be computed
#'
#' @importFrom stats predict median pnorm sd
#' @import ranger
#'
#' @return object of type \code{margins.orf}
#'
#' @examples
#' #\dontrun{
#'
#' ## Ordered Forest
#' require(orf)
#'
#' # load example data
#' data(odata)
#'
#' # specify response and covariates
#' Y <- odata[, 1]
#' X <- odata[, -1]
#'
#' # estimate Ordered Forest
#' set.seed(123)
#' orf <- orf(X, Y)
#'
#' # estimate marginal effects of the orf (default)
#' margins(orf)
#'
#' # estimate marginal effects evaluated at the mean
#' margins(orf, eval = "atmean")
#'
#' # estimate marginal effects with inference
#' # (orf object has to be estimated with honesty and subsampling)
#' margins(orf, inference = TRUE)
#'
#' # estimate marginal effects with custom window size
#' margins(orf, window = 0.5)
#'
#' # estimate marginal effects for some new data (within support of X)
#' margins(orf, newdata = X[1:10, ])
#'
#' # estimate marginal effects with all custom settings
#' margins(orf, eval = "atmedian", inference = TRUE, window = 0.5, newdata = X[1:10, ])
#'
#' #}
#'
#' @export
margins.orf <- function(forest, eval = NULL, inference = NULL, window = NULL, newdata = NULL) {

  # needed inputs for the function: forest - trained forest object of class orf/mrf/brf
  #                                 eval - string defining evaluation point for marginal effects
  #                                 newdata - matrix of new Xs for which the marginal effects should be computed
  # ----------------------------------------------------------------------------------- #

  ## save forest inputs
  inputs            <- forest$forestInfo$inputs
  forest_replace    <- inputs$replace
  forest_honesty    <- inputs$honesty
  forest_inference  <- inputs$inference

  # ----------------------------------------------------------------------------------- #

  # check eval input
  eval <- check_eval(eval)

  # ----------------------------------------------------------------------------------- #

  ## check inference possibilities according to previous estimation
  # if inference not specified, take inference argument as it was in the estimation
  if (is.null(inference)) {

    inference <- forest_inference

  }

  # check if inference is logical
  inference <- check_inference(inference)

  # if inference TRUE, but orf was NOT estimated with subsampling AND honesty, no inference possible
  if (inference == TRUE & (forest_replace != FALSE | forest_honesty != TRUE)) {

    warning("Inference is not possible if the orf object was not estimated with both subsampling and honesty.
            For marginal effects with inference, reestimate orf setting replace = FALSE and honesty = TRUE.")
    inference <- FALSE

  }

  # -------------------------------------------------------------------------------- #

  # decide if prediction or in sample marginal effects should be evaluated
  if (is.null(newdata)) {

    # if no newdata supplied, estimate in sample marginal effects
    if (forest_honesty == FALSE) {

      data <- forest$forestInfo$trainData # take in-sample data
      X_eval <- as.matrix(data[, -1])

    } else if (forest_honesty == TRUE) {

      data <- forest$forestInfo$honestData # take honest data
      X_eval <- as.matrix(data[, -1])

    }

  } else {

    # check if X has name
    if (is.null(colnames(newdata))) { colnames(newdata) <- paste0("X", rep(1:ncol(newdata))) }

    # check if newdata is compatible with train data
    if (all(colnames(forest$forestInfo$trainData)[-1] != colnames(newdata)) | ncol(newdata) != (ncol(forest$forestInfo$trainData) - 1)) {

      stop("newdata is not compatible with training data. Programme terminated.")

    } else {

      # check X data
      check_X(newdata)
      # newdata which will be used for evaluating marginal effects
      X_eval <- as.matrix(newdata)

      # get data which will be used for predicting the marginal effect
      if (forest_honesty == FALSE) {

        data <- forest$forestInfo$trainData # take in-sample data

      } else if (forest_honesty == TRUE) {

        data <- forest$forestInfo$honestData # take honest data

      }

    }

  }
  ### data checks done

  # ----------------------------------------------------------------------------------- #

  ### data preparation and checks
  # check the window size
  window <- check_window(window)
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
  h_std <- window
  # check if X is continuous or dummy or categorical
  X_type <- apply(X, 2, function(x) length(unique(x)))
  # now determine the type of X
  X_continuous <- which(X_type > 10) # define IDs of continuous Xs
  X_dummy <- which(X_type == 2) # define IDs of dummies
  X_categorical <- which(X_type > 2 & X_type <= 10)
  # additional check for constant variables which are nonsensical
  if (any(X_type == 1) | any(X_type == 0)) {
    stop("Some of the covariates are constant. This makes no sense for evaluation of marginal effects. Programme terminated.")
  }

  # ----------------------------------------------------------------------------------- #

  ### check the evaluation point
  if (eval == "atmean") {
    # variable of interest: X_1 to X_last, ME at mean
    X_mean <- lapply(1:ncol(X_eval), function(x) t(as.matrix(colMeans(X_eval)))) # set all Xs to their mean values (so many times as we have Xs)
  } else if (eval == "atmedian") {
    # variable of interest: X_1 to X_last, ME at median
    X_mean <- lapply(1:ncol(X_eval), function(x) t(as.matrix(apply(X_eval, 2, median)))) # set all Xs to their median values (so many times as we have Xs)
  } else if (eval == "mean") {
    # # variable of interest: X_1 to X_last, mean ME
    X_mean <- lapply(1:ncol(X_eval), function(x) X_eval) # set all Xs to their exact values (so many times as we have Xs)
  } else {
    stop("Incorrect evaluation point. This must be one of be one of mean, atmean, or atmedian. Programme terminated.")
  }

  # ----------------------------------------------------------------------------------- #

  ### get data needed for evaluation of ME
  # get number of evaluation points
  X_rows <- nrow(X_mean[[1]])
  # get number of Xs
  X_cols <- ncol(X_mean[[1]])
  # get SD of Xs
  X_sd <- rep_row(apply(X, 2, sd), n = X_rows)
  # create X_up (X_mean + 0.1 * X_sd)
  X_up <- X_mean[[1]] + h_std*X_sd
  # create X_down (X_mean - 0.1 * X_sd)
  X_down <- X_mean[[1]] - h_std*X_sd

  ## now check for the support of X
  # check X_max
  X_max <- rep_row(apply(X, 2, max), n = X_rows)
  # check X_min
  X_min <- rep_row(apply(X, 2, min), n = X_rows)
  # check if X_up is within the range X_min and X_max
  X_up <- (X_up < X_max) * X_up + (X_up >= X_max) * X_max
  X_up <- (X_up > X_min) * X_up + (X_up <= X_min) * (X_min + h_std * X_sd)
  # check if X_down is within the range X_min and X_max
  X_down <- (X_down > X_min) * X_down + (X_down <= X_min) * X_min
  X_down <- (X_down < X_max) * X_down + (X_down >= X_max) * (X_max - h_std * X_sd)
  # check if X_up and X_down are same
  if (any(X_up == X_down)) {
    # adjust to higher share of SD
    X_up   <- (X_up > X_down) * X_up   + (X_up == X_down) * (X_up   + 0.5 * h_std * X_sd)
    X_down <- (X_up > X_down) * X_down + (X_up == X_down) * (X_down - 0.5 * h_std * X_sd)
    # check the min max range again
    X_up   <- (X_up < X_max) * X_up + (X_up >= X_max) * X_max
    X_down <- (X_down > X_min) * X_down + (X_down <= X_min) * X_min

  }
  # checks for support of X done

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

  # adjust for categorical X (works also for zero categorical) (adjustment such that the difference is always 1)
  for (i in X_categorical) {
    X_mean_up[[i]][, i] <- ceiling(X_mean_up[[i]][, i])
    X_mean_down[[i]][, i] <- ifelse(ceiling(X_mean_down[[i]][, i]) == ceiling(X_mean_up[[i]][, i]),
                                    floor(X_mean_down[[i]][, i]),
                                    ceiling(X_mean_down[[i]][, i])
                                    )
  }

  # adjust for dummies (works also for zero dummies)
  for (i in X_dummy) {
    X_mean_up[[i]][, i] <- max(X[, i])
    X_mean_down[[i]][, i] <- min(X[, i])
  }

  # ----------------------------------------------------------------------------------- #

  ### check honesty and inference
  if (forest_honesty == FALSE & inference == FALSE) {

    #### now we do not need weights if we do not need inference (based on out of bag predictions)
    # forest prediction for X_mean_up (mean doesnt matter for atmean or atmedian)
    forest_pred_up <- lapply(forest$trainForests, function(x) lapply(X_mean_up, function(y) mean(predict(x, data = y)$predictions)))
    # forest prediction for X_mean_down (mean doesnt matter for atmean or atmedian)
    forest_pred_down <- lapply(forest$trainForests, function(x) lapply(X_mean_down, function(y) mean(predict(x, data = y)$predictions)))

  } else if (forest_honesty == TRUE & inference == FALSE) {

    # do honest predictions
    # forest prediction for X_mean_up (use new faster function particularly for ME)
    forest_pred_up <- predict_forest_preds_for_ME(forest$trainForests, data_ind, X_mean_up)
    # forest prediction for X_mean_down
    forest_pred_down <- predict_forest_preds_for_ME(forest$trainForests, data_ind, X_mean_down)

  } else if (forest_honesty == TRUE & inference == TRUE) {

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
  forest_pred_up_1 <- append(forest_pred_up, list(rep(list(rep(1, 1)), X_cols))) # append 1
  forest_pred_up_0 <- append(list(rep(list(rep(0, 1)), X_cols)), forest_pred_up) # prepend 0
  # isolate predictions
  forest_pred_up <- mapply(function(x,y) mapply(function(x,y) x-y, x, y,  SIMPLIFY = F), forest_pred_up_1, forest_pred_up_0, SIMPLIFY = F)
  # avoid negative predictions
  forest_pred_up <- lapply(forest_pred_up, function(x) lapply(x, function(x) ifelse((x < 0), 0, x)))
  # normalize predictions
  forest_pred_up_rowsum <- lapply(seq_along(forest_pred_up[[1]]), function(i) rowSums(matrix(sapply(forest_pred_up, "[[", i), ncol = length(categories), nrow = 1))) # build rowsums with respect to categories
  forest_pred_up <- lapply(forest_pred_up, function(x) mapply(function(x,y) x/y, x, forest_pred_up_rowsum, SIMPLIFY = FALSE)) # normalize to sum up to 1 (very rare but just to be sure)

  # prepare down
  forest_pred_down_1 <- append(forest_pred_down, list(rep(list(rep(1, 1)), X_cols))) # append 1
  forest_pred_down_0 <- append(list(rep(list(rep(0, 1)), X_cols)), forest_pred_down) # prepend 0
  # isolate predictions
  forest_pred_down <- mapply(function(x,y) mapply(function(x,y) x-y, x, y,  SIMPLIFY = F), forest_pred_down_1, forest_pred_down_0, SIMPLIFY = F)
  # avoid negative predictions
  forest_pred_down <- lapply(forest_pred_down, function(x) lapply(x, function(x) ifelse((x < 0), 0, x)))
  # normalize predictions
  forest_pred_down_rowsum <- lapply(seq_along(forest_pred_down[[1]]), function(i) rowSums(matrix(sapply(forest_pred_down, "[[", i), ncol = length(categories), nrow = 1))) # build rowsums with respect to categories
  forest_pred_down <- lapply(forest_pred_down, function(x) mapply(function(x,y) x/y, x, forest_pred_down_rowsum, SIMPLIFY = FALSE)) # normalize to sum up to 1 (very rare but just to be sure)

  # ----------------------------------------------------------------------------------- #

  ### now subtract the predictions according to the ME formula
  forest_pred_diff_up_down <- mapply(function(x,y) mapply(function(x,y) x-y, x, y,  SIMPLIFY = F), forest_pred_up, forest_pred_down, SIMPLIFY = F)
  # compute the scaling factor: X_up-X_down=2*X_sd
  scaling_factor <- lapply(1:X_cols, function(i) mean(as.numeric((X_up - X_down)[, i]))) # save it as separate list vectors (mean doesnt change anything for "atmean" option)
  # set scaling factor to zero for categorical and dummy variables
  for (i in (union(X_categorical, X_dummy))) {
    scaling_factor[[i]] <- 1
  }
  # scale the differences to get marginal effects
  marginal_effects_scaled <- lapply(forest_pred_diff_up_down, function(x) mapply(function(x,y) x/y, x, scaling_factor, SIMPLIFY = FALSE) )

  # ----------------------------------------------------------------------------------- #

  # coerce to a matrix
  marginal_effects <- sapply(marginal_effects_scaled, function(x) sapply(x, function(x) as.matrix(x)))
  # add names
  colnames(marginal_effects) <- sapply(categories, function(x) paste("Category", x, sep = " "))
  rownames(marginal_effects) <- colnames(X)

  # ----------------------------------------------------------------------------------- #

  if (inference == TRUE) {

    ### variance for the marginal effects
    ## compute prerequisities for variance of honest marginal effects
    # mean of scaling factor and squared afterwards (for atmean and atmedian the averaging doesnt change anything)
    scaling_factor_squared <- lapply(scaling_factor, function(x) (mean(x))^2)

    # now subtract the weights according to the ME formula
    forest_weights_diff_up_down <- mapply(function(x,y) mapply(function(x,y) x-y, x, y,  SIMPLIFY = F), forest_weights_up, forest_weights_down, SIMPLIFY = F)

    # compute the conditional means: 1/N(weights%*%y) (predictions are based on honest sample)
    forest_cond_means <- mapply(function(x,y) lapply(x, function(x) (x%*%y)/nrow(Y_ind[[1]])), forest_weights_diff_up_down, Y_ind, SIMPLIFY = FALSE)
    # compute standard multiplication
    forest_multi <- mapply(function(x,y) lapply(x, function(x) t(x)*y), forest_weights_diff_up_down, Y_ind, SIMPLIFY = FALSE)
    # subtract the mean from each obs i
    forest_multi_demeaned <- mapply(function(x,y) mapply(function(x,y) x-matrix(y, nrow = nrow(x)), x, y, SIMPLIFY = FALSE), forest_multi, forest_cond_means, SIMPLIFY = F)

    ## now do the single variances for each category m
    # square the demeaned and sum it and normalize
    forest_multi_demeaned_sq_sum_norm <- lapply(forest_multi_demeaned, function(x) lapply(x, function(x) (sum(x^2))*(nrow(Y_ind[[1]])/(nrow(Y_ind[[1]])-1))))

    ## now compute the prerequisites for covariances
    # multiply forest_var_multi_demeaned according to formula for covariance (shifted categories needed for computational convenience)
    forest_multi_demeaned_0_last <- append(forest_multi_demeaned, list(rep(list(matrix(0, ncol = ncol(forest_multi_demeaned[[1]][[1]]), nrow = nrow(forest_multi_demeaned[[1]][[1]]))), ncol(X_eval)))) # append zero matrix list
    forest_multi_demeaned_0_first <- append(list(rep(list(matrix(0, ncol = ncol(forest_multi_demeaned[[1]][[1]]), nrow = nrow(forest_multi_demeaned[[1]][[1]]))), ncol(X_eval))), forest_multi_demeaned) # prepend zero matrix list
    # compute the multiplication of category m with m-1 according to the covariance formula (sum, normalize and multiply by 2)
    forest_multi_demeaned_cov_sum_norm_mult2 <- mapply(function(x,y) mapply(function(x,y) (sum(x*y))*(nrow(Y_ind[[1]])/(nrow(Y_ind[[1]])-1))*2, x, y, SIMPLIFY = FALSE), forest_multi_demeaned_0_first, forest_multi_demeaned_0_last, SIMPLIFY = F)

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
    # put everything together according to formula: var_last + var_first - cov
    variance_me <- mapply(function(x,y,z) mapply(function(x,y,z) (x+y-z), x, y, z, SIMPLIFY = FALSE), variance_last, variance_first, covariance, SIMPLIFY = F)

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

  } else {

    # no values for the other parameters if inference is not desired
    variance_me <- NULL
    sd_me       <- NULL
    t_value     <- NULL
    p_values    <- NULL

    # put everything into a list of results
    results <- list(marginal_effects, variance_me, sd_me, t_value, p_values)
    names(results) <- c("MarginalEffects", "Variances", "StandardErrors", "tValues", "pValues")

  }

  # ----------------------------------------------------------------------------------- #

  class(results) <- "margins.orf"

  # return results
  return(results)

  # ----------------------------------------------------------------------------------- #

}


#' print.margins.orf
#'
#' print estimated marginal effects from ordered random forest of class \code{margins.orf}
#'
#' @param x object of type \code{margins.orf}
#' @param latex logical, if latex output should be generated (\code{default = FALSE})
#' @param ... further arguments (currently ignored)
#'
#' @examples
#' #\dontrun{
#'
#' ## Ordered Forest
#' require(orf)
#'
#' # load example data
#' data(odata)
#'
#' # specify response and covariates
#' Y <- odata[, 1]
#' X <- odata[, -1]
#'
#' # estimate Ordered Forest
#' set.seed(123)
#' orf <- orf(X, Y)
#'
#' # estimate marginal effects of the orf
#' orf_margins <- margins(orf)
#'
#' # print marginal effects
#' print(orf_margins)
#'
#' # print marginal effects coded in LaTeX
#' print(orf_margins, latex = TRUE)
#'
#' #}
#'
#' @export
print.margins.orf <- function(x, latex = FALSE, ...) {

  # chekc if inference has been done
  if (!is.null(x$Variances) & latex == FALSE) {

    # print inference output table
    margins_output(x)

   } else if (!is.null(x$Variances) & latex == TRUE) {

    # print inference output table latex
    margins_output_latex(x)

  } else {

    # print just the marginal effects
    print(x$MarginalEffects)

  }


}


#' margins.mrf
#'
#' estimate marginal effects of the multinomial random forest
#'
#' @param forest trained multinomial random forest object of class \code{mrf}
#' @param eval string defining evaluation point for marginal effects. These can be one of "mean", "atmean", or "atmedian". (Default is "mean")
#' @param inference logical, if TRUE inference on marginal effects will be conducted (default is inherited from the mrf object)
#' @param window numeric, share of standard deviation of X to be used for evaluation of the marginal effect (default is 0.1)
#' @param newdata matrix of new Xs for which marginal effects will be computed
#'
#' @importFrom stats predict median pnorm sd
#' @import ranger
#'
#' @return object of type \code{margins.mrf}
#'
#' @export
margins.mrf <- function(forest, eval = NULL, inference = NULL, window = NULL, newdata = NULL) {

  # needed inputs for the function: forest - trained forest object of class orf/mrf/brf
  #                                 eval - string defining evaluation point for marginal effects (default is atmean)
  #                                 newdata - matrix of new Xs
  # ----------------------------------------------------------------------------------- #

  ## save forest inputs
  inputs            <- forest$forestInfo$inputs
  forest_replace    <- inputs$replace
  forest_honesty    <- inputs$honesty
  forest_inference  <- inputs$inference

  # ----------------------------------------------------------------------------------- #

  # check eval input
  eval <- check_eval(eval)

  # ----------------------------------------------------------------------------------- #

  ## check inference possibilities according to previous estimation
  # if inference not specified, take inference argument as it was in the estimation
  if (is.null(inference)) {

    inference <- forest_inference

  }

  # check if inference is logical
  inference <- check_inference(inference)

  # if inference TRUE, but orf was NOT estimated with subsampling AND honesty, no inference possible
  if (inference == TRUE & (forest_replace != FALSE | forest_honesty != TRUE)) {

    warning("Inference is not possible if the orf object was not estimated with both subsampling and honesty.
            For marginal effects with inference, reestimate orf setting replace = FALSE and honesty = TRUE.")
    inference <- FALSE

  }

  # -------------------------------------------------------------------------------- #

  ### decide if prediction or in sample marginal effects should be evaluated
  if (is.null(newdata)) {

    # if no newdata supplied, estimate in sample marginal effects
    if (forest_honesty == FALSE) {

      data <- forest$forestInfo$trainData # take in-sample data
      X_eval <- as.matrix(data[, -1])

    } else if (forest_honesty == TRUE) {

      data <- forest$forestInfo$honestData # take honest data
      X_eval <- as.matrix(data[, -1])

    }

  } else {

    # check if newdata is compatible with train data
    if (ncol(newdata) != (ncol(forest$forestInfo$trainData) - 1)) {

      stop("newdata is not compatible with training data. Programme terminated.")

    } else {

      # check X data
      check_X(newdata)
      # newdata which will be used for evaluating marginal effects
      X_eval <- as.matrix(newdata)

      # get data which will be used for predicting the marginal effect
      if (forest_honesty == FALSE) {

        data <- forest$forestInfo$trainData # take in-sample data

      } else if (forest_honesty == TRUE) {

        data <- forest$forestInfo$honestData # take honest data

      }

    }

  }
  ### data checks done

  # ----------------------------------------------------------------------------------- #

  ### data preparation and checks
  # check the window size
  window <- check_window(window)
  # get number of observations
  n_data <- as.numeric(nrow(data))
  # get categories
  categories <- forest$forestInfo$categories
  # get X as matrix
  X <- as.matrix(data[, -1])
  # get Y as matrix
  Y <- as.matrix(data[, 1])
  # create indicator variables (outcomes) now with equality for each single category
  Y_ind <- lapply(categories, function(x) ifelse((Y == x), 1, 0))
  # create datasets with indicator outcomes
  data_ind <- lapply(Y_ind, function(x) as.data.frame(cbind(as.matrix(unlist(x)), X)))

  # ----------------------------------------------------------------------------------- #

  ### marginal effects preparation
  # share of SD to be used
  h_std <- window
  # check if X is continuous or dummy or categorical
  X_type <- apply(X, 2, function(x) length(unique(x)))
  # now determine the type of X
  X_continuous <- which(X_type > 10) # define IDs of continuous Xs
  X_dummy <- which(X_type == 2) # define IDs of dummies
  X_categorical <- which(X_type > 2 & X_type <= 10)
  # additional check for constant variables which are nonsensical
  if (any(X_type == 1) | any(X_type == 0)) {
    stop("Some of the covariates are constant. This makes no sense for evaluation of marginal effects. Programme terminated.")
  }

  # ----------------------------------------------------------------------------------- #

  ### check the evaluation point
  if (eval == "atmean") {
    # variable of interest: X_1 to X_last, ME at mean
    X_mean <- lapply(1:ncol(X_eval), function(x) t(as.matrix(colMeans(X_eval)))) # set all Xs to their mean values (so many times as we have Xs)
  } else if (eval == "atmedian") {
    # variable of interest: X_1 to X_last, ME at median
    X_mean <- lapply(1:ncol(X_eval), function(x) t(as.matrix(apply(X_eval, 2, median)))) # set all Xs to their median values (so many times as we have Xs)
  } else if (eval == "mean") {
    # # variable of interest: X_1 to X_last, mean ME
    X_mean <- lapply(1:ncol(X_eval), function(x) X_eval) # set all Xs to their exact values (so many times as we have Xs)
  } else {
    stop("Incorrect evaluation point. This must be one of be one of mean, atmean, or atmedian. Programme terminated.")
  }

  # ----------------------------------------------------------------------------------- #

  ### get data needed for evaluation of ME
  # get number of evaluation points
  X_rows <- nrow(X_mean[[1]])
  # get number of Xs
  X_cols <- ncol(X_mean[[1]])
  # get SD of Xs
  X_sd <- rep_row(apply(X, 2, sd), n = X_rows)
  # create X_up (X_mean + 0.1 * X_sd)
  X_up <- X_mean[[1]] + h_std*X_sd
  # create X_down (X_mean - 0.1 * X_sd)
  X_down <- X_mean[[1]] - h_std*X_sd

  ## now check for the support of X
  # check X_max
  X_max <- rep_row(apply(X, 2, max), n = X_rows)
  # check X_min
  X_min <- rep_row(apply(X, 2, min), n = X_rows)
  # check if X_up is within the range X_min and X_max
  X_up <- (X_up < X_max) * X_up + (X_up >= X_max) * X_max
  X_up <- (X_up > X_min) * X_up + (X_up <= X_min) * (X_min + h_std * X_sd)
  # check if X_down is within the range X_min and X_max
  X_down <- (X_down > X_min) * X_down + (X_down <= X_min) * X_min
  X_down <- (X_down < X_max) * X_down + (X_down >= X_max) * (X_max - h_std * X_sd)
  # check if X_up and X_down are same
  if (any(X_up == X_down)) {
    # adjust to higher share of SD
    X_up   <- (X_up > X_down) * X_up   + (X_up == X_down) * (X_up   + 0.5 * h_std * X_sd)
    X_down <- (X_up > X_down) * X_down + (X_up == X_down) * (X_down - 0.5 * h_std * X_sd)
    # check the min max range again
    X_up   <- (X_up < X_max) * X_up + (X_up >= X_max) * X_max
    X_down <- (X_down > X_min) * X_down + (X_down <= X_min) * X_min

  }
  # checks for support of X done

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

  # adjust for categorical X (works also for zero categorical) (adjustment such that the difference is always 1)
  for (i in X_categorical) {
    X_mean_up[[i]][, i] <- ceiling(X_mean_up[[i]][, i])
    X_mean_down[[i]][, i] <- ifelse(ceiling(X_mean_down[[i]][, i]) == ceiling(X_mean_up[[i]][, i]),
                                    floor(X_mean_down[[i]][, i]),
                                    ceiling(X_mean_down[[i]][, i])
    )
  }

  # adjust for dummies (works also for zero dummies)
  for (i in X_dummy) {
    X_mean_up[[i]][, i] <- max(X[, i])
    X_mean_down[[i]][, i] <- min(X[, i])
  }

  # ----------------------------------------------------------------------------------- #

  ### check honesty and inference
  if (forest_honesty == FALSE & inference == FALSE) {

    #### now we do not need weights if we do not need inference (based on out of bag predictions)
    # forest prediction for X_mean_up (mean doesnt matter for atmean or atmedian)
    forest_pred_up <- lapply(forest$trainForests, function(x) lapply(X_mean_up, function(y) mean(predict(x, data = y)$predictions)))
    # forest prediction for X_mean_down (mean doesnt matter for atmean or atmedian)
    forest_pred_down <- lapply(forest$trainForests, function(x) lapply(X_mean_down, function(y) mean(predict(x, data = y)$predictions)))

  } else if (forest_honesty == TRUE & inference == FALSE) {

    # do honest predictions
    # forest prediction for X_mean_up (use new faster function particularly for ME)
    forest_pred_up <- predict_forest_preds_for_ME(forest$trainForests, data_ind, X_mean_up)
    # forest prediction for X_mean_down
    forest_pred_down <- predict_forest_preds_for_ME(forest$trainForests, data_ind, X_mean_down)

  } else if (forest_honesty == TRUE & inference == TRUE) {

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

  ### form MRF predictions
  ## now we have to normalize the predictions for MRF X_up
  # get rowsum for each X
  forest_pred_up_rowsum <- lapply(seq_along(forest_pred_up[[1]]), function(i) rowSums(matrix(sapply(forest_pred_up, "[[", i), ncol = length(categories), nrow = 1))) # build rowsums with respect to categories
  # normalize, i.e. divide each element by the sum of elements (for each X row)
  forest_pred_up_norm <- lapply(forest_pred_up, function(x) mapply(function(x,y) x/y, x, forest_pred_up_rowsum, SIMPLIFY = FALSE)) # normalize to sum up to 1 (very rare but just to be sure)

  ## now we have to normalize the predictions for MRF X_down
  # get rowsum for each X
  forest_pred_down_rowsum <- lapply(seq_along(forest_pred_down[[1]]), function(i) rowSums(matrix(sapply(forest_pred_down, "[[", i), ncol = length(categories), nrow = 1))) # build rowsums with respect to categories
  # normalize, i.e. divide each element by the sum of elements (for each X row)
  forest_pred_down_norm <- lapply(forest_pred_down, function(x) mapply(function(x,y) x/y, x, forest_pred_down_rowsum, SIMPLIFY = FALSE)) # normalize to sum up to 1 (very rare but just to be sure)

  # ----------------------------------------------------------------------------------- #

  ### now subtract the predictions according to the ME formula
  forest_pred_diff_up_down <- mapply(function(x,y) mapply(function(x,y) x-y, x, y,  SIMPLIFY = F), forest_pred_up_norm, forest_pred_down_norm, SIMPLIFY = F)
  # compute the scaling factor: X_up-X_down=2*X_sd
  scaling_factor <- lapply(1:X_cols, function(i) mean(as.numeric((X_up - X_down)[, i]))) # save it as separate list vectors (mean doesnt change anything for "atmean" option)
  # set scaling factor to zero for categorical and dummy variables
  for (i in (union(X_categorical, X_dummy))) {
    scaling_factor[[i]] <- 1
  }
  # scale the differences to get marginal effects
  marginal_effects_scaled <- lapply(forest_pred_diff_up_down, function(x) mapply(function(x,y) x/y, x, scaling_factor, SIMPLIFY = FALSE) )

  # ----------------------------------------------------------------------------------- #

  # coerce to a matrix
  marginal_effects <- sapply(marginal_effects_scaled, function(x) sapply(x, function(x) as.matrix(x)))

  # add names
  colnames(marginal_effects) <- sapply(categories, function(x) paste("Category", x, sep = " "))
  rownames(marginal_effects) <- colnames(X)

  # ----------------------------------------------------------------------------------- #

  if (inference == TRUE) {

    ### variance for the marginal effects
    ## compute prerequisities for variance of honest marginal effects
    # mean of scaling factor and squared afterwards (for atmean and atmedian the averaging doesnt change anything)
    scaling_factor_squared <- lapply(scaling_factor, function(x) (mean(x))^2)

    # now subtract the weights according to the ME formula
    forest_weights_diff_up_down <- mapply(function(x,y) mapply(function(x,y) x-y, x, y,  SIMPLIFY = F), forest_weights_up, forest_weights_down, SIMPLIFY = F)

    # compute the conditional means: 1/N(weights%*%y) (predictions are based on honest sample)
    forest_cond_means <- mapply(function(x,y) lapply(x, function(x) (x%*%y)/nrow(Y_ind[[1]])), forest_weights_diff_up_down, Y_ind, SIMPLIFY = FALSE)

    # compute standard multiplication
    forest_multi <- mapply(function(x,y) lapply(x, function(x) t(x)*y), forest_weights_diff_up_down, Y_ind, SIMPLIFY = FALSE)

    # subtract the mean from each obs i
    forest_multi_demeaned <- mapply(function(x,y) mapply(function(x,y) x-matrix(y, nrow = nrow(x)), x, y, SIMPLIFY = FALSE), forest_multi, forest_cond_means, SIMPLIFY = F)

    ## now do the single variances for each category m
    # square the demeaned and sum it and normalize (# square the demeaned, # sum all obs i together, # multiply by N/N-1 (normalize))
    forest_multi_demeaned_sq_sum_norm <- lapply(forest_multi_demeaned, function(x) lapply(x, function(x) (sum(x^2))*(nrow(Y_ind[[1]])/(nrow(Y_ind[[1]])-1))))

    # divide by scaling factor to get the variance
    variance <- lapply(forest_multi_demeaned_sq_sum_norm, function(x) mapply(function(x,y) x/y, x, scaling_factor_squared, SIMPLIFY = FALSE) )
    ## single variances done

    # ----------------------------------------------------------------------------------- #

    ## covariances not needed for MRF

    # ----------------------------------------------------------------------------------- #

    ## output for final variances of marginal effects
    # coerce to a matrix
    variance_me <- sapply(variance, function(x) sapply(x, function(x) as.matrix(x)))

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

  } else {

    # no values for the other parameters if inference is not desired
    variance_me <- NULL
    sd_me       <- NULL
    t_value     <- NULL
    p_values    <- NULL

    # put everything into a list of results
    results <- list(marginal_effects, variance_me, sd_me, t_value, p_values)
    names(results) <- c("MarginalEffects", "Variances", "StandardErrors", "tValues", "pValues")

  }

  # ----------------------------------------------------------------------------------- #

  class(results) <- "margins.mrf"

  # return results
  return(results)

  # ----------------------------------------------------------------------------------- #

}



#' print.margins.mrf
#'
#' print estimated marginal effects from multinomial random forest of class \code{margins.mrf}
#'
#' @param x object of type \code{margins.mrf}
#' @param latex logical, if latex output should be generated (\code{default = FALSE})
#' @param ... further arguments (currently ignored)
#'
#' @export
print.margins.mrf <- function(x, latex = FALSE, ...) {

  # chekc if inference has been done
  if (!is.null(x$Variances) & latex == FALSE) {

    # print inference output table
    margins_output(x)

  } else if (!is.null(x$Variances) & latex == TRUE) {

    # print inference output table latex
    margins_output_latex(x)

  } else if (is.null(x$Variances) & latex == FALSE){

    # print just the marginal effects
    print(x$MarginalEffects)

  } else {

    # put caption an latex environment
    xoutput <- xtable(x$MarginalEffects, digits = 4, caption = "ORF Marginal Effects")
    # print xtable
    print.xtable(xoutput, type = "latex", include.rownames = TRUE, comment = FALSE)

  }

}

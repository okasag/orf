#' Marginal Effects
#'
#' S3 generic method for estimation of marginal effects
#' of an Ordered Forest objects of class \code{orf}.
#'
#' @seealso \code{\link{margins.orf}}, \code{\link{summary.margins.orf}} and \code{\link{print.margins.orf}}
#'
#' @param forest estimated Ordered Forest object of class \code{orf}
#' @param eval string, defining evaluation point for marginal effects. These can be one of "mean", "atmean", or "atmedian". (Default is "mean")
#' @param inference logical, if TRUE inference on marginal effects will be conducted (default is inherited from the \code{orf} object)
#' @param window numeric, share of standard deviation of X to be used for evaluation of the marginal effect (default is 0.1)
#' @param newdata numeric matrix X containing the new observations for which the marginal effects should be estimated
#'
#' @author Gabriel Okasa
#'
#' @examples
#' ## Ordered Forest
#' require(orf)
#'
#' # load example data
#' data(odata)
#'
#' # specify response and covariates
#' Y <- as.numeric(odata[, 1])
#' X <- as.matrix(odata[, -1])
#'
#' # estimate Ordered Forest
#' orf_fit <- orf(X, Y)
#' \donttest{
#' # estimate default marginal effects of the orf
#' orf_margins <- margins(orf_fit)
#' }
#'
#' @export
margins <- function(forest, eval = NULL, inference = NULL, window = NULL, newdata = NULL) UseMethod("margins")


#' @export
margins.default <- function(forest, eval = NULL, inference = NULL, window = NULL, newdata = NULL) {

  warning(paste("margins does not know how to handle object of class ",
                class(forest),
                ". The supported classes are one of the following: orf."))

}


#' Marginal Effects for the Ordered Forest
#'
#' S3 method for estimation of marginal effects
#' of an Ordered Forest objects of class \code{orf}.
#'
#' \code{margins.orf} estimates marginal effects at the mean, at the median, or
#' the mean marginal effects, depending on the \code{eval} argument. It is advised
#' to increase the number of subsampling replications in the supplied \code{orf}
#' object as the estimation of the marginal effects is a more demanding exercise
#' than a simple Ordered Forest estimation/prediction. Additionally to the estimation
#' of the marginal effects, the weight-based inference for the effects is supported
#' as well. Note, that the inference procedure is much more computationally exhausting
#' exercise due to the computation of the forest weights. Additionally, the evaluation
#' window for the marginal effects can be regulated through the \code{window} argument.
#' Furthermore, new data for which marginal effects should be computed can be supplied
#' as well as long as it lies within the support of \code{X}.
#'
#' @seealso \code{\link{summary.margins.orf}}, \code{\link{print.margins.orf}}
#'
#' @param forest estimated Ordered Forest object of class \code{orf}
#' @param eval string, defining evaluation point for marginal effects. These can be one of "mean", "atmean", or "atmedian". (Default is "mean")
#' @param inference logical, if TRUE inference on marginal effects will be conducted (default is inherited from the \code{orf} object)
#' @param window numeric, share of standard deviation of X to be used for evaluation of the marginal effect (default is 0.1)
#' @param newdata numeric matrix X containing the new observations for which the marginal effects should be estimated
#
#' @importFrom stats predict median pnorm sd
#' @import ranger
#'
#' @return object of type \code{margins.orf} with following elements
#'       \item{info}{info containing forest inputs and data used}
#'       \item{effects}{marginal effects}
#'       \item{variances}{variances of marginal effects}
#'       \item{errors}{standard errors of marginal effects}
#'       \item{tvalues}{t-values of marginal effects}
#'       \item{pvalues}{p-values of marginal effects}
#'
#' @author Gabriel Okasa
#'
#' @examples
#' ## Ordered Forest
#' require(orf)
#'
#' # load example data
#' data(odata)
#'
#' # specify response and covariates
#' Y <- as.numeric(odata[, 1])
#' X <- as.matrix(odata[, -1])
#'
#' # estimate Ordered Forest
#' orf_fit <- orf(X, Y)
#' \donttest{
#' # estimate marginal effects of the orf (default)
#' orf_margins <- margins(orf_fit)
#'
#' # estimate marginal effects evaluated at the mean
#' orf_margins <- margins(orf_fit, eval = "atmean")
#'
#' # estimate marginal effects with inference
#' # (orf object has to be estimated with honesty and subsampling)
#' orf_margins <- margins(orf_fit, inference = TRUE)
#'
#' # estimate marginal effects with custom window size
#' orf_margins <- margins(orf_fit, window = 0.5)
#'
#' # estimate marginal effects for some new data (within support of X)
#' orf_margins <- margins(orf_fit, newdata = X[1:10, ])
#'
#' # estimate marginal effects with all custom settings
#' orf_margins <- margins(orf_fit, eval = "atmedian", inference = TRUE,
#'                                 window = 0.5, newdata = X[1:10, ])
#' }
#'
#' @export
margins.orf <- function(forest, eval = NULL, inference = NULL, window = NULL, newdata = NULL) {

  # ----------------------------------------------------------------------------------- #

  ## save forest inputs
  inputs            <- forest$info$inputs
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

      data   <- forest$info$trainData # take in-sample data
      X_eval <- as.matrix(data[, -1])

    } else if (forest_honesty == TRUE) {

      data   <- forest$info$honestData # take honest data
      X_eval <- as.matrix(data[, -1])

    }

  } else {

    # check if X has name
    if (is.null(colnames(newdata))) { colnames(newdata) <- paste0("X", rep(1:ncol(newdata))) }

    # check if newdata is compatible with train data
    if (all(colnames(forest$info$trainData)[-1] != colnames(newdata)) | ncol(newdata) != (ncol(forest$info$trainData) - 1)) {

      stop("newdata is not compatible with training data. Programme terminated.")

    } else {

      # check X data
      check_X(newdata)
      # newdata which will be used for evaluating marginal effects
      X_eval <- as.matrix(newdata)

      # get data which will be used for predicting the marginal effect
      if (forest_honesty == FALSE) {

        data <- forest$info$trainData # take in-sample data

      } else if (forest_honesty == TRUE) {

        data <- forest$info$honestData # take honest data

      }

    }

  }

  # ----------------------------------------------------------------------------------- #

  ## data preparation and checks
  # check the window size
  window      <- check_window(window)
  # get number of observations
  n_data      <- as.numeric(nrow(data))
  # get categories
  categories  <- forest$info$categories
  # get X as matrix
  X           <- as.matrix(data[, -1])
  # get Y as matrix
  Y           <- as.matrix(data[, 1])
  # create indicator variables (outcomes)
  Y_ind       <- lapply(categories[1:length(categories)-1], function(x) ifelse((Y <= x), 1, 0))
  # create datasets with indicator outcomes
  data_ind    <- lapply(Y_ind, function(x) as.data.frame(cbind(as.matrix(unlist(x)), X)))

  # ----------------------------------------------------------------------------------- #

  ## marginal effects preparation
  # share of SD to be used
  h_std         <- window
  # check if X is continuous or dummy or categorical
  X_type        <- apply(X, 2, function(x) length(unique(x)))
  # now determine the type of X
  X_continuous  <- which(X_type > 10) # define IDs of continuous Xs
  X_dummy       <- which(X_type == 2) # define IDs of dummies
  X_categorical <- which(X_type > 2 & X_type <= 10)
  # additional check for constant variables which are nonsensical
  if (any(X_type == 1) | any(X_type == 0)) {
    stop("Some of the covariates are constant. This makes no sense for evaluation of marginal effects. Programme terminated.")
  }

  # ----------------------------------------------------------------------------------- #

  # check the evaluation point
  if (eval == "atmean") {
    # variable of interest: X_1 to X_last, ME at mean
    X_mean <- lapply(1:ncol(X_eval), function(x) t(as.matrix(colMeans(X_eval)))) # set all Xs to their mean values (so many times as we have Xs)
  } else if (eval == "atmedian") {
    # variable of interest: X_1 to X_last, ME at median
    X_mean <- lapply(1:ncol(X_eval), function(x) t(as.matrix(apply(X_eval, 2, median)))) # set all Xs to their median values (so many times as we have Xs)
  } else if (eval == "mean") {
    # variable of interest: X_1 to X_last, mean ME
    X_mean <- lapply(1:ncol(X_eval), function(x) X_eval) # set all Xs to their exact values (so many times as we have Xs)
  } else {
    stop("Incorrect evaluation point. This must be one of be one of mean, atmean, or atmedian. Programme terminated.")
  }

  # ----------------------------------------------------------------------------------- #

  ## get data needed for evaluation of ME
  # get number of evaluation points
  X_rows <- nrow(X_mean[[1]])
  # get number of Xs
  X_cols <- ncol(X_mean[[1]])
  # get SD of Xs
  X_sd   <- matrix(rep(apply(X, 2, sd), times = 1, each = X_rows), nrow = X_rows)
  # create X_up (X_mean + 0.1 * X_sd)
  X_up   <- X_mean[[1]] + h_std*X_sd
  # create X_down (X_mean - 0.1 * X_sd)
  X_down <- X_mean[[1]] - h_std*X_sd

  ## now check for the support of X
  # check X_max
  X_max   <- matrix(rep(apply(X, 2, max), times = 1, each = X_rows), nrow = X_rows)
  # check X_min
  X_min   <- matrix(rep(apply(X, 2, min), times = 1, each = X_rows), nrow = X_rows)
  # check if X_up is within the range X_min and X_max
  X_up    <- (X_up < X_max) * X_up + (X_up >= X_max) * X_max
  X_up    <- (X_up > X_min) * X_up + (X_up <= X_min) * (X_min + h_std * X_sd)
  # check if X_down is within the range X_min and X_max
  X_down  <- (X_down > X_min) * X_down + (X_down <= X_min) * X_min
  X_down  <- (X_down < X_max) * X_down + (X_down >= X_max) * (X_max - h_std * X_sd)
  # check if X_up and X_down are same
  if (any(X_up == X_down)) {
    # adjust to higher share of SD
    X_up   <- (X_up > X_down) * X_up   + (X_up == X_down) * (X_up   + 0.5 * h_std * X_sd)
    X_down <- (X_up > X_down) * X_down + (X_up == X_down) * (X_down - 0.5 * h_std * X_sd)
    # check the min max range again
    X_up   <- (X_up < X_max) * X_up + (X_up >= X_max) * X_max
    X_down <- (X_down > X_min) * X_down + (X_down <= X_min) * X_min

  }

  # ----------------------------------------------------------------------------------- #

  ## now we need 2 datasets: one with X_up and second with X_down
  # X_mean_up all (continous)
  X_mean_up   <- X_mean
  X_mean_down <- X_mean

  # replace values accordingly
  for (i in 1:X_cols) {
    X_mean_up[[i]][, i]   <- X_up[, i]
    X_mean_down[[i]][, i] <- X_down[, i]
  }

  # adjust for categorical X (works also for zero categorical) (adjustment such that the difference is always 1)
  for (i in X_categorical) {
    X_mean_up[[i]][, i]   <- ceiling(X_mean_up[[i]][, i])
    X_mean_down[[i]][, i] <- ifelse(ceiling(X_mean_down[[i]][, i]) == ceiling(X_mean_up[[i]][, i]),
                                    floor(X_mean_down[[i]][, i]),
                                    ceiling(X_mean_down[[i]][, i])
                                    )
  }

  # adjust for dummies (works also for zero dummies)
  for (i in X_dummy) {
    X_mean_up[[i]][, i]   <- max(X[, i])
    X_mean_down[[i]][, i] <- min(X[, i])
  }

  # ----------------------------------------------------------------------------------- #

  ## check honesty and inference
  if (forest_honesty == FALSE & inference == FALSE) {

    ## now we do not need weights if we do not need inference (based on out of bag predictions)
    # forest prediction for X_mean_up (mean doesnt matter for atmean or atmedian)
    forest_pred_up <- lapply(forest$forests, function(x) lapply(X_mean_up, function(y) mean(predict(x, data = y)$predictions)))
    # forest prediction for X_mean_down (mean doesnt matter for atmean or atmedian)
    forest_pred_down <- lapply(forest$forests, function(x) lapply(X_mean_down, function(y) mean(predict(x, data = y)$predictions)))

  } else if (forest_honesty == TRUE & inference == FALSE) {

    ## do honest predictions
    # forest prediction for X_mean_up
    forest_pred_up <- predict_forest_preds_for_ME(forest$forests, data_ind, X_mean_up)
    # forest prediction for X_mean_down
    forest_pred_down <- predict_forest_preds_for_ME(forest$forests, data_ind, X_mean_down)

  } else if (forest_honesty == TRUE & inference == TRUE) {

    ## do honest predictions with weight based inference
    # extract weights for desired Xs up: get weights from honest sample and predict weights for evaluation points from HONEST sample
    forest_weights_up <- predict_forest_weights_for_ME(forest$forests, X, X_mean_up)
    # extract weights for desired Xs down
    forest_weights_down <- predict_forest_weights_for_ME(forest$forests, X, X_mean_down)

    ## compute predictions based on weights
    # forest prediction for X_mean_up
    forest_pred_up <- mapply(function(x,y) lapply(x, function(x) as.numeric(x%*%y)), forest_weights_up, Y_ind, SIMPLIFY = FALSE)
    # forest prediction for X_mean_down
    forest_pred_down <- mapply(function(x,y) lapply(x, function(x) as.numeric(x%*%y)), forest_weights_down, Y_ind, SIMPLIFY = FALSE)

  }

  # ----------------------------------------------------------------------------------- #

  ## form ORF predictions
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

  ## now subtract the predictions according to the ME formula
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

    ## variance for the marginal effects
    # compute prerequisities for variance of honest marginal effects
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
    # multiply forest_multi_demeaned according to formula for covariance (shifted categories needed for computational convenience)
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

    # t values and p values
    t_value <- (marginal_effects)/(sd_me)
    # control for dividing zero by zero
    t_value[is.nan(t_value)] = 0
    # p values
    p_values <- 2*pnorm(-abs(t_value))

    # ----------------------------------------------------------------------------------- #

  } else {

    # no values for the other parameters if inference is not desired
    variance_me <- NULL
    sd_me       <- NULL
    t_value     <- NULL
    p_values    <- NULL

    # ----------------------------------------------------------------------------------- #

  }

  # ----------------------------------------------------------------------------------- #

  # save forest information
  forest_info <- list(inputs, categories, eval, window, newdata, inference)
  names(forest_info) <- c("inputs", "categories", "eval", "window", "newData", "marginsInference")

  # put everything into a list of results
  results <- list(forest_info, marginal_effects, variance_me, sd_me, t_value, p_values)
  names(results) <- c("info", "effects", "variances", "errors", "tvalues", "pvalues")

  class(results) <- "margins.orf"

  # return results
  return(results)

  # ----------------------------------------------------------------------------------- #

}


#' Summary of the Ordered Forest Marginal Effects
#'
#' summary of estimated marginal effects of the Ordered Forest of class \code{margins.orf}
#'
#' \code{summary.margins.orf} provides estimation results of the Ordered Forest
#' marginal effects. The summary contains the results for the marginal effects
#' for each covariate and each outcome class, optionally with inference as well.
#' Furthermore, summary output as a LaTeX table is supported in order to directly
#' extract the results for the documentation.
#'
#' @param object estimated Ordered Forest Marginal Effect object of type \code{margins.orf}
#' @param latex logical, if TRUE latex coded summary will be generated (default is FALSE)
#' @param ... further arguments (currently ignored)
#'
#' @author Gabriel Okasa
#'
#' @examples
#' ## Ordered Forest
#' require(orf)
#'
#' # load example data
#' data(odata)
#'
#' # specify response and covariates
#' Y <- as.numeric(odata[, 1])
#' X <- as.matrix(odata[, -1])
#'
#' # estimate Ordered Forest
#' orf_fit <- orf(X, Y)
#' \donttest{
#' # estimate marginal effects of the orf
#' orf_margins <- margins(orf_fit)
#'
#' # summary of marginal effects
#' summary(orf_margins)
#'
#' # summary of marginal effects coded in LaTeX
#' summary(orf_margins, latex = TRUE)
#' }
#'
#' @export
summary.margins.orf <- function(object, latex = FALSE, ...) {

  # -------------------------------------------------------------------------------- #

  ## check user inputs
  latex <- check_latex(latex)
  # get object as orf_margins
  orf_margins <- object

  ## save forest margins inputs
  main_class        <- class(orf_margins)[1]
  inputs            <- orf_margins$info$inputs

  honesty           <- inputs$honesty
  honesty.fraction  <- inputs$honesty.fraction
  mtry              <- inputs$mtry
  num.trees         <- inputs$num.trees
  min.node.size     <- inputs$min.node.size
  replace           <- inputs$replace
  sample.fraction   <- inputs$sample.fraction
  inference         <- inputs$inference

  pred_data         <- ifelse(is.null(orf_margins$info$newData), FALSE, TRUE)
  eval_type         <- orf_margins$info$eval
  eval_window       <- orf_margins$info$window
  margins_inference <- orf_margins$info$marginsInference
  categories        <- length(orf_margins$info$categories)
  build             <- ifelse(replace == TRUE, "Bootstrap", "Subsampling")
  type              <- "Ordered Forest Margins"

  # -------------------------------------------------------------------------------- #

  # general output

  # -------------------------------------------------------------------------------- #

  # structure summary into a list
  output        <- list(type, eval_type, eval_window, pred_data, categories, build, num.trees, mtry, min.node.size, replace, sample.fraction, honesty, honesty.fraction, margins_inference)
  names(output) <- c("type", "evaluation.type", "evaluation.window", "new.data", "categories", "build", "num.trees", "mtry", "min.node.size", "replace", "sample.fraction", "honesty", "honesty.fraction", "inference")

  # output matrix
  output_matrix <- matrix(NA, ncol = 1, nrow = length(output))
  # populate output matrix
  rownames(output_matrix) <- names(output) # rownames are names
  colnames(output_matrix) <- "" # no visible colname
  output_matrix[, 1]      <- unlist(output) # column 1 are values

  # generate latex output if selected
  if (latex == TRUE) { colnames(output_matrix) <- "Ordered Forest Margins Summary"
  output_matrix <- xtable(output_matrix, caption = "Summary of the Ordered Forest Margins", align = "ll")
  }

  # pack it into output
  output <- output_matrix

  # -------------------------------------------------------------------------------- #

  cat("Summary of the", type, "\n\n")

  # return output
  print(noquote(output), comment = FALSE)
  cat("\n")

  # -------------------------------------------------------------------------------- #

  # coefficients output

  # -------------------------------------------------------------------------------- #

  # chekc if inference has been done
  if (!is.null(orf_margins$variances) & latex == FALSE) {

     # print inference output table
     margins_output(orf_margins)

   } else if (!is.null(orf_margins$variances) & latex == TRUE) {

     # print inference output table latex
     margins_output_latex(orf_margins)

   } else if (is.null(orf_margins$variances) & latex == TRUE) {

     # put caption an latex environment
     xoutput <- xtable(orf_margins$effects, digits = 4, caption = "ORF Marginal Effects")
     # put hline after each variable
     print.xtable(xoutput, type = "latex", include.rownames = TRUE, comment = FALSE)

   } else {

     # print marginal effects title
     cat("ORF Marginal Effects: \n\n")
     # print just the marginal effects
     print(round(orf_margins$effects, 4))

  }

  # -------------------------------------------------------------------------------- #

}


#' Print of the Ordered Forest Marginal Effects
#'
#' print of estimated marginal effects of the Ordered Forest of class \code{margins.orf}
#'
#' \code{print.margins.orf} provides a first glimpse of the Ordered Forest
#' marginal effects, printed directly to the \code{R} console. The printed information
#' contains the results for the marginal effects for each covariate and each outcome class.
#'
#' @param x estimated Ordered Forest Marginal Effect object of type \code{margins.orf}
#' @param ... further arguments (currently ignored)
#'
#' @author Gabriel Okasa
#'
#' @examples
#' ## Ordered Forest
#' require(orf)
#'
#' # load example data
#' data(odata)
#'
#' # specify response and covariates
#' Y <- as.numeric(odata[, 1])
#' X <- as.matrix(odata[, -1])
#'
#' # estimate Ordered Forest
#' orf_fit <- orf(X, Y)
#' \donttest{
#' # estimate marginal effects of the orf
#' orf_margins <- margins(orf_fit)
#'
#' # print marginal effects
#' print(orf_margins)
#' }
#'
#' @export
print.margins.orf <- function(x, ...) {

  # -------------------------------------------------------------------------------- #

  # save x as orf_margins
  orf_margins       <- x

  ## save forest prediction inputs
  main_class        <- class(orf_margins)[1]
  inputs            <- orf_margins$info$inputs

  honesty           <- inputs$honesty
  mtry              <- inputs$mtry
  num.trees         <- inputs$num.trees
  min.node.size     <- inputs$min.node.size
  replace           <- inputs$replace
  inference         <- inputs$inference

  pred_data         <- orf_margins$info$newData
  eval_type         <- orf_margins$info$eval
  eval_window       <- orf_margins$info$window
  margins_inference <- orf_margins$info$marginsInference
  categories        <- length(orf_margins$info$categories)
  build             <- ifelse(replace == TRUE, "Bootstrap", "Subsampling")
  type              <- "Ordered Forest Margins"

  # -------------------------------------------------------------------------------- #

  cat(type, "object of class", main_class, "\n\n")

  cat("Evaluation Type:                 ", eval_type, "\n")
  cat("Evaluation Window:               ", eval_window, "\n")
  cat("Number of Categories:            ", categories, "\n")
  cat("New Data:                        ", ifelse(is.null(pred_data), FALSE, TRUE), "\n")
  cat("Number of Trees:                 ", num.trees, "\n")
  cat("Build:                           ", build, "\n")
  cat("Mtry:                            ", mtry, "\n")
  cat("Minimum Node Size:               ", min.node.size, "\n")
  cat("Honest Forest:                   ", honesty, "\n")
  cat("Weight-Based Inference:          ", margins_inference, "\n\n")

  # -------------------------------------------------------------------------------- #

  # print marginal effects title
  cat("ORF Marginal Effects: \n\n")
  # print just the marginal effects
  print(round(x$effects, 4))

  # -------------------------------------------------------------------------------- #

}


#' Formatted output for marginal effects with inference
#'
#' function for creating inference table output for estimated effects which
#' can be passed into \code{print.margins.orf}
#'
#' @param x object of type \code{margins.orf}
#'
#' @keywords internal
#'
margins_output <- function(x) {

  output_matrix <- matrix(NA, nrow = 1, ncol = 4)

  cat("ORF Marginal Effects: \n\n")
  cat("---------------------------------------------------------------------------", "\n")

  for (var_idx in 1:nrow(x$effects)) {

    cat(rownames(x$effects)[var_idx], "\n")
    cat("                   Class", "     Effect", "    StdErr", "    tValue ", "   pValue", "     ", "\n")

    for (cat_idx in 1:ncol(x$effects)) {

      # generate stars (inspired by and thanks to: http://myowelt.blogspot.com/2008/04/beautiful-correlation-tables-in-r.html)
      stars <- ifelse(x$pvalues[var_idx, cat_idx] < .01, "***",
                      ifelse(x$pvalues[var_idx, cat_idx] < .05, "** ",
                             ifelse(x$pvalues[var_idx, cat_idx] < .1, "*  ", "   ")))

      # print estimates for each category iteratively
      output_matrix[1, 1] <- x$effects[var_idx, cat_idx]
      output_matrix[1, 2] <- x$errors[var_idx, cat_idx]
      output_matrix[1, 3] <- x$tvalues[var_idx, cat_idx]
      output_matrix[1, 4] <- x$pvalues[var_idx, cat_idx]


      cat("                    ", cat_idx, "     ") # prit out the categories

      cat(format(sprintf("%8.4f", round(output_matrix, 4)), width = 10), stars, "     ") # print out the estimates

      cat("\n") # break the line


    }

  }

  cat("---------------------------------------------------------------------------", "\n")
  cat("Significance levels correspond to: *** .< 0.01, ** .< 0.05, * .< 0.1 \n")
  cat("---------------------------------------------------------------------------", "\n")

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
#' @keywords internal
#'
margins_output_latex <- function(x) {

  # get number of categories
  ncat      <- ncol(x$effects)
  # get number of variables
  nvar      <- nrow(x$effects)
  # get variable names
  varnames  <- rownames(x$effects)

  # create empty output matrix
  output_matrix <- matrix("", nrow = (nvar*ncat), ncol = 7)
  rownames(output_matrix) <- rep("default", (nvar*ncat)) # generate unique identifier

  for (var_idx in 0:(nvar-1)) {


    for (cat_idx in 1:ncat) {

      # generate stars (inspired by and thanks to: http://myowelt.blogspot.com/2008/04/beautiful-correlation-tables-in-r.html)
      stars <- ifelse(x$pvalues[var_idx+1, cat_idx] < .01, "***",
                      ifelse(x$pvalues[var_idx+1, cat_idx] < .05, "** ",
                             ifelse(x$pvalues[var_idx+1, cat_idx] < .1, "*  ", "   ")))

      # print estimates for each category iteratively
      rownames(output_matrix)[(var_idx*ncat)+cat_idx] <- paste0(varnames[var_idx+1], cat_idx)
      output_matrix[1+(var_idx*ncat), 1] <- varnames[var_idx+1] # fit in variable name
      output_matrix[(var_idx*ncat)+cat_idx, 2] <- cat_idx # fit in category
      output_matrix[(var_idx*ncat)+cat_idx, 3] <- x$effects[var_idx+1, cat_idx]
      output_matrix[(var_idx*ncat)+cat_idx, 4] <- x$errors[var_idx+1, cat_idx]
      output_matrix[(var_idx*ncat)+cat_idx, 5] <- x$tvalues[var_idx+1, cat_idx]
      output_matrix[(var_idx*ncat)+cat_idx, 6] <- x$pvalues[var_idx+1, cat_idx]
      output_matrix[(var_idx*ncat)+cat_idx, 7] <- stars



    }

  }

  # add colnames
  colnames(output_matrix) <- c("Variable", "Class", "Effect", "Std.Error", "t-Value", "p-Value", " ")
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

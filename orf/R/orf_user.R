#' orf
#'
#' Ordered Forest for flexible estimation of the ordered choice model as developed in Lechner & Okasa (2019)
#'
#' @param X matrix of features
#' @param Y vector of outcomes (as.matrix acceptable too)
#' @param num.trees scalar, number of trees in a forest, i.e. bootstrap replications (default is 1000 trees)
#' @param mtry scalar, number of randomly selected features (default is the squared root of number of features, rounded up to the nearest integer)
#' @param min.node.size scalar, minimum node size (default is 5 observations)
#' @param replace logical, if TRUE sampling with replacement, i.e. bootstrap is used to grow the trees, otherwise subsampling without replacement is used (default is set to FALSE)
#' @param sample.fraction scalar, subsampling rate (default is 1 for bootstrap and 0.5 for subsampling)
#' @param honesty logical, if TRUE honest forest is built using 50:50 data split (default is set to TRUE)
#' @param honesty.fraction scalar, share of observations belonging to honest sample not used for growing the forest (default is 0.5)
#' @param inference logical, if TRUE the weight based inference is conducted (default is set to FALSE)
#' @param importance logical, if TRUE variable importance measure based on permutation is conducted (default is set to FALSE)
#'
#' @import ranger
#'
#' @return object of type \code{orf} with following elements
#'       \item{trainForests}{saved forests trained for ORF estimations (inherited from \code{ranger})}
#'       \item{forestInfo}{info containing forest inputs and data used}
#'       \item{forestPredictions}{predicted values}
#'       \item{forestVariances}{variances of predicted values}
#'       \item{variableImportance}{weighted measure of permutation based variable importance}
#'       \item{MSE}{in-sample mean squared error}
#'       \item{RPS}{in-sample ranked probability score}
#'
#' @examples
#' \dontrun{
#'
#' # Ordered Forest with default settings
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
#' # print output of the orf estimation
#' print(orf)
#'
#' # show summary of the orf estimation
#' summary(orf)
#'
#' # plot the estimated probability distributions
#' plot(orf)
#'
#' # predict with the estimated orf
#' predict(orf)
#'
#' # estimate marginal effects of the orf
#' margins(orf)
#'
#' }
#'
#' @export
orf <- function(X, Y,
                num.trees = 1000,
                mtry = NULL,
                min.node.size = NULL,
                replace = FALSE,
                sample.fraction = NULL,
                honesty = TRUE,
                honesty.fraction = NULL,
                inference = FALSE,
                importance = FALSE) {

  # needed inputs for the function: X - matrix of features
  #                                 Y - vector of outcomes (as.matrix acceptable too)
  #                                 num.trees - number of trees in a forest
  #                                 mtry - number of randomly selected features
  #                                 min.node.size - minimum node size
  #                                 replace - sampling with or without replacement
  #                                 sample.fraction - subsampling rate
  #                                 honesty - logical, if TRUE honest forest is built using 50:50 data split
  #                                 inference - logical, if TRUE the weight based inference is conducted (honesty has to be TRUE)
  #                                 importance - logical, if TRUE variable importance measure based on permutation is conducted
  # -------------------------------------------------------------------------------- #

  ## standard checks for input data
  check_X(X)
  Y                <- check_Y(Y, X)
  Y                <- check_discrete_Y(Y)
  num.trees        <- check_num_trees(num.trees)
  mtry             <- check_mtry(mtry, X)
  min.node.size    <- check_min_node_size(min.node.size, X)
  replace          <- check_replace(replace)
  sample.fraction  <- check_sample_fraction(sample.fraction, replace)
  honesty          <- check_honesty(honesty)
  honesty.fraction <- check_honesty_fraction(honesty.fraction, honesty)
  inference        <- check_inference(inference)
  importance       <- check_importance(importance)

  # -------------------------------------------------------------------------------- #

  ## check for plausibility of options first:
  if (honesty == FALSE & inference == TRUE) {

    warning("For conducting inference honesty is required. Honesty has been set to TRUE.")
    # set honesty to TRUE
    honesty <- TRUE

  }

  if (replace == TRUE & inference == TRUE) {

    warning("For conducting inference subsampling is required. Replace has been set to FALSE.")
    # set replace to FALSE
    replace <- FALSE

  }

  # --------------------------------------------------------------------------------------- #

  ## save the inputs:
  inputs <- list(num.trees, mtry, min.node.size, replace, sample.fraction, honesty, honesty.fraction, inference, importance)
  names(inputs) <- c("num.trees", "mtry", "min.node.size", "replace", "sample.fraction", "honesty", "honesty.fraction", "inference", "importance")

  ## save colnames
  # Y - numeric response as only regression is supported (so far)
  Y <- as.matrix(as.numeric(Y)) # make sure you can add colname
  if (is.null(colnames(Y))) { colnames(Y) <- "Y" } # check if Y has name
  Y_name <- colnames(Y) # save the name of Y

  # X
  if (is.null(colnames(X))) { colnames(X) <- paste0("X", rep(1:ncol(X))) } # check if X has name
  X_name <- colnames(X) # save the name of X

  ## set needed dataframe and local variables
  dat <- as.data.frame(cbind(Y, X)) # dataframe
  colnames(dat) <- c(Y_name, X_name) # column names
  n <- as.numeric(nrow(dat)) # number of observations

  # parameters (categories)
  categories <- as.numeric(sort(unique(Y))) # sequence of categories
  ncat <- as.numeric(length(categories)) # number of categories
  cat <- categories[1:(ncat-1)] # cat to esitmate / without the last category (not needed cuz P(Y_ind<=last_cat)=1)

  # variable importance definition
  if (importance == TRUE) {

    varimportance <- "permutation"

  } else {

    varimportance <- "none"

  }

  # --------------------------------------------------------------------------------------- #

  ## proceed with estimations
  # decide if honest forest should be estimated
  if (honesty == FALSE & inference == FALSE) {

    # --------------------------------------------------------------------------------------- #

    # no honest splitting, i.e. use all data
    train_data <- dat
    honest_data <- NULL

    ## create variables needed for orf estimations
    # create indicator variables (outcomes)
    Y_ind <- lapply(cat, function(x) ifelse((Y <= x), 1, 0))

    # create dataset for ranger estimation
    data_ind <- lapply(Y_ind, function(x) as.data.frame(cbind(as.matrix(unlist(x)), X)))

    # --------------------------------------------------------------------------------------- #

    # estimate ncat-1 forests (everything on the same data: placing splits and effect estimation), no subsampling
    forest <- lapply(data_ind, function(x) ranger(dependent.variable.name = paste(Y_name), data = x,
                                                  num.trees = num.trees, mtry = mtry, min.node.size = min.node.size,
                                                  replace = replace, sample.fraction = sample.fraction,
                                                  importance = varimportance))

    # --------------------------------------------------------------------------------------- #

    # collect predictions for each forest based on whole sample (oob predictions)
    pred <- lapply(forest, function(x) x$predictions) # collect forest predictions
    # add the probability for the last outcome (always 1)
    pred_1 <- append(pred, list(rep(1, n)))
    # prepend zero vector to predictions for later differencing
    pred_0 <- append(list(rep(0, n)), pred) # append a first 0 elemnt for the list

    # --------------------------------------------------------------------------------------- #

    # total predictions (make sure it returns a list)
    pred_total <- as.list(mapply(function(x,y) x-y, pred_1, pred_0, SIMPLIFY = F))

    # avoid negative predictions
    pred_total <- lapply(pred_total, function(x) ifelse((x < 0), 0, x))

    # coerce to final matrix
    pred_total <- sapply(pred_total, function(x) as.matrix(x))

    # normalize predictions
    pred_final <- matrix(apply(pred_total, 1, function(x) (x)/(sum(x))), ncol = ncat, byrow = T)

    # add names
    colnames(pred_final) <- sapply(categories, function(x) paste("Category", x, sep = " "))

    # compute OOB MSE based on whole sample
    pred_mse <- mse(pred_final, Y)

    # compute OOB RPS based on whole sample
    pred_rps <- rps(pred_final, Y)

    # --------------------------------------------------------------------------------------- #

    ## convert probabilities into class predictions ("classification")
    #pred_class <- as.matrix(apply(pred_final, 1, which.max))
    #colnames(pred_class) <- "Category"

    # --------------------------------------------------------------------------------------- #

    # no inference here
    var_final <- NULL

    # -------------------------------------------------------------------------------- #

  } else if (honesty == TRUE & inference == FALSE) {

    # --------------------------------------------------------------------------------------- #

    ## do honest forest estimation here using 50:50 data split as in Lechner (2018)
    # devide into 50:50 honesty sets
    split_data <- honest_split(dat, honesty.fraction, orf = TRUE)
    # take care of train data
    train_data <- split_data$trainData # take out training data
    rows_train_data <- as.numeric(rownames(train_data)) # take rownames of train data as numeric
    Y_train <- as.matrix(train_data[, 1]) # take out Y train
    colnames(Y_train) <- Y_name # add column name
    X_train <- train_data[, -1] # take out X
    colnames(X_train) <- X_name # add column names
    # take care of honest data
    honest_data <- split_data$honestData # take out honest data
    rows_honest_data <- as.numeric(rownames(honest_data)) # take rownames of train data as numeric
    Y_honest <- as.matrix(honest_data[, 1]) # take out Y train
    colnames(Y_honest) <- Y_name # add column name
    X_honest <- honest_data[, -1] # take out X
    colnames(X_honest) <- X_name # add column names

    # --------------------------------------------------------------------------------------- #

    ## create variables needed for orf estimations
    # create indicator variables (outcomes)
    Y_ind_train <- lapply(cat, function(x) ifelse((Y_train <= x), 1, 0)) # train
    Y_ind_honest <- lapply(cat, function(x) ifelse((Y_honest <= x), 1, 0))

    # create dataset for ranger estimation
    data_ind_train <- lapply(Y_ind_train, function(x) as.data.frame(cbind(as.matrix(unlist(x)), X_train)))
    data_ind_honest <- lapply(Y_ind_honest, function(x) as.data.frame(cbind(as.matrix(unlist(x)), X_honest)))

    # --------------------------------------------------------------------------------------- #

    # estimate ncat-1 forests (everything on the same data: placing splits and effect estimation), no subsampling
    forest <- lapply(data_ind_train, function(x) ranger(dependent.variable.name = paste(Y_name), data = x,
                                                        num.trees = num.trees, mtry = mtry, min.node.size = min.node.size,
                                                        replace = replace, sample.fraction = sample.fraction,
                                                        importance = varimportance))

    # --------------------------------------------------------------------------------------- #

    # compute honest predictions based on honest sample
    pred <- mapply(function(x,y,z) get_honest(x, y, z), forest, data_ind_honest, data_ind_train, SIMPLIFY = F)
    # add the probability for the last outcome (always 1)
    pred_1 <- append(pred, list(rep(1, n)))
    # prepend zero vector to predictions for later differencing
    pred_0 <- append(list(rep(0, n)), pred) # append a first 0 elemnt for the list

    # --------------------------------------------------------------------------------------- #

    # total predictions (make sure it returns a list)
    pred_total <- as.list(mapply(function(x,y) x-y, pred_1, pred_0, SIMPLIFY = F))

    # avoid negative predictions
    pred_total <- lapply(pred_total, function(x) ifelse((x < 0), 0, x))

    # coerce to final matrix
    pred_total <- sapply(pred_total, function(x) as.matrix(x))

    # normalize predictions
    pred_final <- matrix(apply(pred_total, 1, function(x) (x)/(sum(x))), ncol = ncat, byrow = T)

    # add names
    colnames(pred_final) <- sapply(categories, function(x) paste("Category", x, sep = " "))

    # compute honest MSE based on whole sample
    pred_mse <- mse(pred_final, Y)

    # compute honest RPS based on whole sample
    pred_rps <- rps(pred_final, Y)

    # --------------------------------------------------------------------------------------- #

    ## convert probabilities into class predictions ("classification")
    #pred_class <- as.matrix(apply(pred_final, 1, which.max))
    #colnames(pred_class) <- "Category"

    # --------------------------------------------------------------------------------------- #

    # rename for output
    data_ind <- data_ind_honest
    # no inference here
    var_final <- NULL

    # -------------------------------------------------------------------------------- #

  } else if (honesty == TRUE & inference == TRUE) {

    # --------------------------------------------------------------------------------------- #

    ## do honest forest estimation here using 50:50 data split as in Lechner (2018)
    # devide into 50:50 honesty sets
    split_data <- honest_split(dat, honesty.fraction, orf = TRUE)
    # take care of train data
    train_data <- split_data$trainData # take out training data
    rows_train_data <- as.numeric(rownames(train_data)) # take rownames of train data as numeric
    Y_train <- as.matrix(train_data[, 1]) # take out Y train
    colnames(Y_train) <- Y_name # add column name
    X_train <- train_data[, -1] # take out X
    colnames(X_train) <- X_name # add column names
    # take care of honest data
    honest_data <- split_data$honestData # take out honest data
    rows_honest_data <- as.numeric(rownames(honest_data)) # take rownames of train data as numeric
    Y_honest <- as.matrix(honest_data[, 1]) # take out Y train
    colnames(Y_honest) <- Y_name # add column name
    X_honest <- honest_data[, -1] # take out X
    colnames(X_honest) <- X_name # add column names

    # --------------------------------------------------------------------------------------- #

    ## create variables needed for orf estimations
    # create indicator variables (outcomes)
    Y_ind_train <- lapply(cat, function(x) ifelse((Y_train <= x), 1, 0)) # train
    Y_ind_honest <- lapply(cat, function(x) ifelse((Y_honest <= x), 1, 0))

    # create dataset for ranger estimation
    data_ind_train <- lapply(Y_ind_train, function(x) as.data.frame(cbind(as.matrix(unlist(x)), X_train)))
    data_ind_honest <- lapply(Y_ind_honest, function(x) as.data.frame(cbind(as.matrix(unlist(x)), X_honest)))

    # --------------------------------------------------------------------------------------- #

    # estimate ncat-1 forests (everything on the same data: placing splits and effect estimation), no subsampling
    forest <- lapply(data_ind_train, function(x) ranger(dependent.variable.name = paste(Y_name), data = x,
                                                        num.trees = num.trees, mtry = mtry, min.node.size = min.node.size,
                                                        replace = replace, sample.fraction = sample.fraction,
                                                        importance = varimportance))

    # --------------------------------------------------------------------------------------- #

    # get honest weights
    forest_weights <- mapply(function(x,y,z) get_forest_weights(x, y, z), forest, data_ind_honest, data_ind_train, SIMPLIFY = F)
    honest_weights <- lapply(forest_weights, function(x) x[rows_honest_data, ]) # take out honest sample honest weights
    train_weights  <- lapply(forest_weights, function(x) x[rows_train_data, ]) # take out train sample honest weights

    # --------------------------------------------------------------------------------------- #

    ## make honest predictions, i.e. fitted values based on honest sample
    # honest sample predictions
    honest_pred <- mapply(function(x,y) as.matrix(x %*% y[, 1]), honest_weights, data_ind_honest, SIMPLIFY = F) # honest weights for honest data
    # train sample predictions
    train_pred <- mapply(function(x,y) as.matrix(x %*% y[, 1]), train_weights, data_ind_honest, SIMPLIFY = F) # honest weights for train data
    # put the prediction together for whole sample and order them as original data
    forest_pred <- mapply(function(x,y) rbind(x, y), honest_pred, train_pred, SIMPLIFY = F)
    # sort according to rownames
    forest_pred <- lapply(forest_pred, function(x) as.numeric(x[order(as.numeric(row.names(x))), ]))

    # add the probability for the last outcome (always 1)
    pred_1 <- append(forest_pred, list(rep(1, n)))
    # prepend zero vector to predictions for later differencing
    pred_0 <- append(list(rep(0, n)), forest_pred) # append a first 0 elemnt for the list

    # --------------------------------------------------------------------------------------- #

    # total predictions (make sure it returns a list)
    pred_total <- as.list(mapply(function(x,y) x-y, pred_1, pred_0, SIMPLIFY = F))

    # avoid negative predictions
    pred_total <- lapply(pred_total, function(x) ifelse((x < 0), 0, x))

    # coerce to final matrix
    pred_total <- sapply(pred_total, function(x) as.matrix(x))

    # normalize predictions
    pred_final <- matrix(apply(pred_total, 1, function(x) (x)/(sum(x))), ncol = ncat, byrow = T)

    # add names
    colnames(pred_final) <- sapply(categories, function(x) paste("Category", x, sep = " "))

    # compute honest MSE based on whole sample
    pred_mse <- mse(pred_final, Y)

    # compute honest RPS based on whole sample
    pred_rps <- rps(pred_final, Y)

    # --------------------------------------------------------------------------------------- #

    ## convert probabilities into class predictions ("classification")
    #pred_class <- as.matrix(apply(pred_final, 1, which.max))
    #colnames(pred_class) <- "Category"

    # --------------------------------------------------------------------------------------- #

    # compute the variances for the categorical predictions
    var_final <- get_orf_variance(honest_pred, honest_weights, train_pred, train_weights, Y_ind_honest)
    # check for normalization
    #sd_final <- sqrt(var_final)

    # --------------------------------------------------------------------------------------- #

    # rename for output
    data_ind <- data_ind_honest

    # --------------------------------------------------------------------------------------- #

  }

  # -------------------------------------------------------------------------------- #

  # compute simple variable importance if it is desired
  if (importance == TRUE) {

    # pre-requesities for variable importance computation
    var_imp_shares      <- matrix(unlist(lapply(Y_ind_train, function(x) mean(x))), nrow = 1) # matrix of proportions for each but last class
    var_imp_shares_0    <- matrix(c(0, var_imp_shares[1, 1:(ncol(var_imp_shares)-1)]), nrow = 1) # shifted matrix for ease of computation
    var_imp_shares_last <- var_imp_shares[1, ncol(var_imp_shares)] # get the last share of classes
    var_imp_forests     <- matrix(unlist(lapply(forest, function(x) x$variable.importance)), nrow = ncol(var_imp_shares), byrow = T) # get the variable importances for eahc of binary forests

    # compute scaling factor using shares
    var_imp_scaling     <- matrix(rep(t((var_imp_shares - var_imp_shares_0)/var_imp_shares_last), ncol(X_train)), ncol = ncol(X_train))

    # compute var_imp using scaling factor and importance from binary forests
    var_imp             <- as.numeric(round(colSums(var_imp_scaling * var_imp_forests), 4))
    names(var_imp)      <- X_name

  } else {

    var_imp <- NULL

  }

  # -------------------------------------------------------------------------------- #

  # save forest information
  forest_info <- list(inputs, train_data, honest_data, categories, data_ind)
  names(forest_info) <- c("inputs", "trainData", "honestData", "categories", "indicatorData")

  # define output of the function
  output <- list(forest, forest_info, pred_final, var_final, var_imp, pred_mse, pred_rps)
  names(output) <- c("trainForests",  "forestInfo", "forestPredictions", "forestVariances", "variableImportance", "MSE", "RPS")

  # --------------------------------------------------------------------------------------- #

  ## Set the name for the class
  class(output) <- "orf"

  # -------------------------------------------------------------------------------- #

  # return the output of the function
  return(output)

}


#' predict.orf
#'
#' Prediction for new observations based on estimated ordered random forest
#' of type \code{orf}
#'
#' @param object estimated forest object of type \code{orf}
#' @param newdata matrix X containing the observations to predict
#' @param type string specifying the type of the prediction, These can be either "probs" or  "p" for probabilities and "class" or "c" for classes. (Default is "probs").
#' @param inference logical, if TRUE variances for the predictions will be estimated (only feasible for probability predictions).
#' @param ... further arguments (currently ignored)
#'
#' @import ranger
#'
#' @return object of class \code{orf.prediction} with following elements
#'       \item{forestInfo}{info containing forest inputs and data used}
#'       \item{forestPredictions}{predicted values}
#'       \item{forestVariances}{variances of predicted values}
#'
#' @export
predict.orf <- function(object, newdata = NULL, type = NULL, inference = NULL, ...) {
  # needed inputs for the function: forest - orf object coming from orf function
  #                                 newdata - matrix X containing the observations to predict

  # -------------------------------------------------------------------------------- #

  ## standard checks for input data
  if (class(object) != "orf") {
    stop("Forest object is not of class orf. Programme terminated.")
  }

  ## get forest as na object
  forest <- object
  ## save forest inputs
  inputs <- forest$forestInfo$inputs
  categories <- forest$forestInfo$categories
  replace <- inputs$replace
  honesty <- inputs$honesty
  honest_data <- forest$forestInfo$honestData
  train_data <- forest$forestInfo$trainData
  honest_ind_data <- forest$forestInfo$indicatorData # indicator data needed for indicator predictions

  # -------------------------------------------------------------------------------- #

  # if inference not specified, take inference argument as it was in the estimation
  if (is.null(inference)) {

    inference <- inputs$inference

  }

  # check if inference is logical
  inference <- check_inference(inference)

  # check if type is admissible
  type <- check_type(type)

  # check if inference is possible according to prediction type
  if (type == "class" | type == "c") {

   inference <- FALSE

  }

  ## get fitted values if newdata = NULL and inference arguments coincide
  if (is.null(newdata) & (inference == inputs$inference)) {

    # take in sample predictions
    pred_final <- forest$forestPredictions
    var_final  <- forest$forestVariances

    # check the desired type of predictions
    if (type == "class" | type == "c") {

      # convert probabilities into class predictions ("ordered classification")
      pred_final <- as.matrix(apply(pred_final, 1, which.max))
      var_final  <- NULL

    }

  } else if (is.null(newdata) & (inference == FALSE) & (inputs$inference == TRUE)) {

    # then take the estimated values but dont supply the inference results
    pred_final <- forest$forestPredictions
    var_final  <- NULL

    # check the desired type of predictions
    if (type == "class" | type == "c") {

      # convert probabilities into class predictions ("ordered classification")
      pred_final <- as.matrix(apply(pred_final, 1, which.max))

    }

  } else {

    # take out list of ranger objects (be careful, its Forests with S at the end!)
    forest <- forest$trainForests

    ## get train data names (only X)
    train_data_name <- colnames(train_data)[2:ncol(train_data)]

    if (is.null(newdata) & (inference == TRUE) & (inputs$inference == FALSE)) {

      # note that the newdata is the estimation data and inference should reflect that
      flag_newdata <- 1
      # take traindata as newdata and estimate the weights needed for inference
      newdata <- as.data.frame(rbind(train_data, honest_data))
      # sort according to rownames
      newdata <- as.data.frame(newdata[order(as.numeric(row.names(newdata))), ])
      # get further values
      n_newdata <- nrow(newdata) # rows of new data
      n_cat <- as.numeric(length(categories))

    } else if (!is.null(newdata)) {

      ## get X matrix newdata as dataframe and check colnames
      # X
      if (is.null(colnames(newdata))) { colnames(newdata) <- paste0("X", rep(1:ncol(newdata))) } # check if X has name
      newdata_name <- colnames(newdata) # save the name of X
      newdata <- as.data.frame(newdata) # as dataframe
      n_newdata <- nrow(newdata) # rows of new data
      n_cat <- as.numeric(length(categories))

      # -------------------------------------------------------------------------------- #

      # check if its compatible with the data used for training
      if (all(train_data_name != newdata_name) | (ncol(newdata) != ncol(train_data)-1)) {

        stop("Newdata are not compatible with the training data. Check supplied data. Programme terminated.")

      }

    }

    # -------------------------------------------------------------------------------- #

    ## check inference possibilities according to previous estimation
    # if inference TRUE, but orf was NOT estimated with subsampling AND honesty, no inference possible
    if (inference == TRUE & (replace != FALSE | honesty != TRUE)) {

      warning("Inference is not possible if the orf object was not estimated with both subsampling and honesty.
            For predictions with inference, reestimate orf setting replace = FALSE and honesty = TRUE.")
      inference <- FALSE

    }

    # -------------------------------------------------------------------------------- #

    # check if honest forest was estimated and predict accordingly
    if (honesty == FALSE & inference == FALSE) {

      ## no honest splitting, i.e. use all data, no inference

      # predict with ncat-1 forests as in ranger default
      pred <- lapply(forest, function(x) predict(x, data = newdata)$predictions)

      # -------------------------------------------------------------------------------- #

      # add the probability for the last outcome (always 1)
      pred_1 <- append(pred, list(rep(1, n_newdata)))
      # prepend zero vector to predictions for later differencing
      pred_0 <- append(list(rep(0, n_newdata)), pred) # append a first 0 elemnt for the list

      # --------------------------------------------------------------------------------------- #

      # total predictions (make sure it returns a list)
      pred_total <- as.list(mapply(function(x,y) x-y, pred_1, pred_0, SIMPLIFY = F))

      # avoid negative predictions
      pred_total <- lapply(pred_total, function(x) ifelse((x < 0), 0, x))

      # coerce to final matrix
      pred_total <- sapply(pred_total, function(x) as.matrix(x))

      # normalize predictions
      pred_final <- matrix(apply(pred_total, 1, function(x) (x)/(sum(x))), ncol = n_cat, byrow = T)

      # add names
      colnames(pred_final) <- sapply(categories, function(x) paste("Category", x, sep = " "))

      # --------------------------------------------------------------------------------------- #

      # check the desired type of predictions
      if (type == "class" | type == "c") {

        # convert probabilities into class predictions ("ordered classification")
        pred_final <- as.matrix(apply(pred_final, 1, which.max))

      }

      # --------------------------------------------------------------------------------------- #

      # no variances
      var_final <- NULL

      # -------------------------------------------------------------------------------- #

    } else if (honesty == TRUE & inference == FALSE) {

      # -------------------------------------------------------------------------------- #

      ## run new Xs through estimated train forest and compute predictions based on honest sample
      # no need to predict weights, get predictions directly through leaves
      # predict with ncat-1 forests
      pred <- mapply(function(x,y) predict_honest(x, y, newdata), forest, honest_ind_data, SIMPLIFY = FALSE)

      # -------------------------------------------------------------------------------- #

      # add the probability for the last outcome (always 1)
      pred_1 <- append(pred, list(rep(1, n_newdata)))
      # prepend zero vector to predictions for later differencing
      pred_0 <- append(list(rep(0, n_newdata)), pred) # append a first 0 elemnt for the list

      # --------------------------------------------------------------------------------------- #

      # total predictions (make sure it returns a list)
      pred_total <- as.list(mapply(function(x,y) x-y, pred_1, pred_0, SIMPLIFY = F))

      # avoid negative predictions
      pred_total <- lapply(pred_total, function(x) ifelse((x < 0), 0, x))

      # coerce to final matrix
      pred_total <- sapply(pred_total, function(x) as.matrix(x))

      # normalize predictions
      pred_final <- matrix(apply(pred_total, 1, function(x) (x)/(sum(x))), ncol = n_cat, byrow = T)

      # add names
      colnames(pred_final) <- sapply(categories, function(x) paste("Category", x, sep = " "))

      # --------------------------------------------------------------------------------------- #

      # check the desired type of predictions
      if (type == "class" | type == "c") {

        # convert probabilities into class predictions ("ordered classification")
        pred_final <- as.matrix(apply(pred_final, 1, which.max))

      }

      # --------------------------------------------------------------------------------------- #

      # no variances
      var_final <- NULL

      # -------------------------------------------------------------------------------- #

    } else if (honesty == TRUE & inference == TRUE) {

        # -------------------------------------------------------------------------------- #

        # predict weights by using forest train, honest data and newdata (for each category except one)
        forest_weights_pred <- lapply(forest, function(x) predict_forest_weights(x, honest_data, newdata))

        # get predictions by matrix multiplication of weights with honest responses
        forest_pred <- mapply(function(x,y) as.numeric(x%*%y[, 1]), forest_weights_pred, honest_ind_data, SIMPLIFY = FALSE)

        # -------------------------------------------------------------------------------- #

        # add the probability for the last outcome (always 1)
        pred_1 <- append(forest_pred, list(rep(1, n_newdata)))
        # prepend zero vector to predictions for later differencing
        pred_0 <- append(list(rep(0, n_newdata)), forest_pred) # append a first 0 elemnt for the list

        # --------------------------------------------------------------------------------------- #

        # total predictions (make sure it returns a list)
        pred_total <- as.list(mapply(function(x,y) x-y, pred_1, pred_0, SIMPLIFY = F))

        # avoid negative predictions
        pred_total <- lapply(pred_total, function(x) ifelse((x < 0), 0, x))

        # coerce to final matrix
        pred_total <- sapply(pred_total, function(x) as.matrix(x))

        # normalize predictions
        pred_final <- matrix(apply(pred_total, 1, function(x) (x)/(sum(x))), ncol = n_cat, byrow = T)

        # add names
        colnames(pred_final) <- sapply(categories, function(x) paste("Category", x, sep = " "))

        # --------------------------------------------------------------------------------------- #

        # adapt variance computations according to newdata
        if (flag_newdata == 1) {

          # use get_orf_variance
          # compute variances as in within estimation (smaller normalization parameter due to sample size)
          # divide forest_pred and forest_weights into 2 pieces for honest and train (ordering doesnt matter)
          n_halfdata <- floor(length(forest_pred[[1]])/2)
          n_fulldata <- length(forest_pred[[1]])
          # transform to matrix vector
          forest_pred <- lapply(forest_pred, function(x) matrix(x, ncol = 1))
          # take care of rownames
          forest_pred <- lapply(forest_pred, function(x) { rownames(x) <- (1:n_fulldata); return(x) })
          forest_weights_pred <- lapply(forest_weights_pred, function(x) { rownames(x) <- (1:n_fulldata); return(x) })
          # predictions (vectors as matrices)
          honest_pred <- lapply(forest_pred, function(x) as.matrix(x[(1:n_halfdata), ]))
          train_pred  <- lapply(forest_pred, function(x) as.matrix(x[((n_halfdata + 1):n_fulldata), ]))
          # weights (matrices)
          honest_weights_pred <- lapply(forest_weights_pred, function(x) x[(1:n_halfdata), ])
          train_weights_pred  <- lapply(forest_weights_pred, function(x) x[((n_halfdata + 1):n_fulldata), ])

          # get indicator outcomes out of honest indicator data
          Y_ind_honest <- lapply(honest_ind_data, function(x) matrix(x[, 1], ncol = 1))

          # compute the in sample variance
          var_final <- get_orf_variance(honest_pred, honest_weights_pred, train_pred, train_weights_pred, Y_ind_honest)

        } else {

          # use standard pred_orf_variance
          # get indicator outcomes out of honest indicator data
          Y_ind_honest <- lapply(honest_ind_data, function(x) x[, 1])

          # compute the variances for the categorical predictions
          var_final <- pred_orf_variance(forest_pred, forest_weights_pred, Y_ind_honest)

        }

        # --------------------------------------------------------------------------------------- #

    }

    # -------------------------------------------------------------------------------- #

  }

  # -------------------------------------------------------------------------------- #

  # save forest information
  forest_info <- list(inputs, categories, newdata, type, inference)
  names(forest_info) <- c("inputs", "categories", "newData", "predType", "predInference")

  # define output of the function
  output <- list(forest_info, pred_final, var_final)
  names(output) <- c("forestInfo", "forestPredictions", "forestVariances")

  # -------------------------------------------------------------------------------- #

  ## Set the name for the class
  class(output) <- "orf.prediction"

  # return output
  return(output)

  # -------------------------------------------------------------------------------- #

}


#' plot.orf
#'
#' plot ordered random forest object of class \code{orf}
#'
#' @param x estimated ordered random forest object of type \code{orf}
#' @param ... further arguments (currently ignored)
#'
#' @import ggplot2
#' @importFrom utils stack
#'
#' @export
plot.orf <- function(x, ...) {

  # needed inputs for the function: forest - forest object coming from orf function

  # -------------------------------------------------------------------------------- #

  ## get forest as x
  forest <- x
  ## save forest inputs
  inputs <- forest$forestInfo$inputs
  honesty <- inputs$honesty
  categories <- forest$forestInfo$categories
  honest_data <- forest$forestInfo$honestData
  train_data <- forest$forestInfo$trainData

  # -------------------------------------------------------------------------------- #

  # get predictions and estimation data
  probabilities <- forest$forestPredictions # take out honest predictions
  all_data <- rbind(honest_data, train_data) # put data together
  all_data <- all_data[order(as.numeric(row.names(all_data))), ] # sort data as original
  outcomes <- all_data[, 1] # take the observed outcomes

  # -------------------------------------------------------------------------------- #

  ### plot ORF ###

  ## plot realized categories overlayed with predicted category probabilities
  # new colnames
  colnames(probabilities) <- sapply(seq_along(categories), function(i) paste0("P(Y=", i, ")"))

  # cbind together
  df_plot <- as.data.frame(cbind(outcomes, probabilities))

  # subset according to categories
  df_plot_cat <- lapply(seq_along(categories), function(i) as.data.frame(subset(df_plot, outcomes == i)))
  # take colmeans
  df_cat_means <- lapply(df_plot_cat, function(x) t(as.matrix(colMeans(x)[seq_along(categories)+1])))
  # add colmeans to df_plot_cat
  df_plot_cat <- mapply(function(x,y) cbind(x, y), df_plot_cat, df_cat_means, SIMPLIFY = FALSE)

  # reshape data for ggplot
  df_plot_prob <- lapply(df_plot_cat, function(x) stack(x[seq_along(categories)+1]))
  df_plot_mean <- lapply(df_plot_cat, function(x) stack(x[seq_along(categories)+1+length(categories)]))

  # add colnames and column indicating the outcome category
  df_plot_prob <- lapply(seq_along(df_plot_prob), function(i) {
    # add column of outcome category to eahc list entry
    df_plot_prob[[i]]$Outcome <- paste("Class", i, sep = " ")
    # add colnames
    colnames(df_plot_prob[[i]]) <- c("Probability", "Density", "Outcome")
    # return the list
    return(df_plot_prob[[i]])  })

  # stack the dataframes under each other
  df_plot_prob <- as.data.frame(do.call(rbind, df_plot_prob))

  # add colnames and column indicating the outcome category
  df_plot_mean <- lapply(seq_along(df_plot_mean), function(i) {
    # add column of outcome category to eahc list entry
    df_plot_mean[[i]]$Outcome <- paste("Class", i, sep = " ")
    # add colnames
    colnames(df_plot_mean[[i]]) <- c("Probability", "Density", "Outcome")
    # return the list
    return(df_plot_mean[[i]])  })

  # stack the dataframes under each other
  df_plot_mean <- as.data.frame(do.call(rbind, df_plot_mean))

  # -------------------------------------------------------------------------------- #

  # generate ggplot
  ggplot(data = df_plot_prob, aes_string(x = "Probability", fill = "Density"))+
    geom_density(alpha = 0.4, aes_string(y = "..scaled..")) +
    facet_wrap("Outcome", ncol = 1)+
    geom_vline(data = df_plot_mean, aes_string(xintercept = "Probability", color = "Density"), linetype="dashed") +
    ggtitle("Distribution of Ordered Forest Probability Predictions") +
    xlab("Predicted Probability") +
    ylab("Probability Mass") +
    theme_bw() +
    theme(strip.background = element_rect(fill = "gray92")) +
    theme(legend.position = "top") +
    theme(plot.title = element_text(hjust = 0.5))

  # no output to return for plot

  # -------------------------------------------------------------------------------- #

}


#' summary.orf
#'
#' summary of an ordered random forest object of class \code{orf}
#'
#' @param object estimated ordered random forest object of type \code{orf}
#' @param latex logical, TRUE if latex summary should be generated
#' @param ... further arguments (currently ignored)
#'
#' @importFrom xtable xtable
#'
#' @export
summary.orf <- function(object, latex = FALSE, ...) {

  # needed inputs for the function: forest - forest object coming from random_forest function
  #                                        - latex : logical if the output should be printed in latex code

  # -------------------------------------------------------------------------------- #

  ## check user inputs
  latex <- check_latex(latex)

  ## get forest as object
  forest <- object

  # -------------------------------------------------------------------------------- #

  ## save forest inputs
  main_class        <- class(forest)[1]
  inputs            <- forest$forestInfo$inputs

  honesty           <- inputs$honesty
  honesty.fraction  <- inputs$honesty.fraction
  inference         <- inputs$inference
  importance        <- inputs$importance
  mtry              <- inputs$mtry
  num.trees         <- inputs$num.trees
  min.node.size     <- inputs$min.node.size
  replace           <- inputs$replace
  sample.fraction   <- inputs$sample.fraction
  honest_data       <- forest$forestInfo$honestData
  train_data        <- forest$forestInfo$trainData
  categories        <- length(forest$forestInfo$categories)
  type              <- "Ordered Forest"

  # -------------------------------------------------------------------------------- #

  ## honest splitting, i.e. use honest data
  # take out summary statistics
  mse         <- round(forest$MSE, 5)
  rps         <- round(forest$RPS, 5)
  trainsize   <- nrow(train_data)
  honestsize  <- ifelse(is.null(honest_data), 0, nrow(honest_data))
  features    <- ncol(train_data) - 1   # take out the response

  # check if subsampling or bootstrapping was used
  if (forest$trainForests[[1]]$replace == TRUE) { build <- "Bootstrap" } else { build <- "Subsampling" }

  # -------------------------------------------------------------------------------- #

  # structure summary into a list
  output        <- list(type, categories, build, num.trees, mtry, min.node.size, replace, sample.fraction, honesty, honesty.fraction, inference, importance, trainsize, honestsize, features, mse, rps)
  names(output) <- c("type", "categories", "build", "num.trees", "mtry", "min.node.size", "replace", "sample.fraction", "honesty", "honesty.fraction", "inference", "importance", "trainsize", "honestsize", "features", "mse", "rps")

  # output matrix
  output_matrix <- matrix(NA, ncol = 1, nrow = length(output))
  # populate output matrix
  rownames(output_matrix) <- names(output) # rownames are names
  colnames(output_matrix) <- "" # no visible colname
  output_matrix[, 1]      <- unlist(output) # column 1 are values

  # generate latex output if selected
  if (latex == TRUE) { colnames(output_matrix) <- "Ordered Forest Summary"
  output_matrix <- xtable(output_matrix, caption = "Summary of the Ordered Forest Estimation", align = "ll")
  }

  # pack it into output
  output <- output_matrix

  # -------------------------------------------------------------------------------- #

  cat("Summary of the", type, "Estimation \n\n")

  # return output
  print(noquote(output), comment = FALSE)

  # -------------------------------------------------------------------------------- #

}


#' print.orf
#'
#' print of an ordered forest object of class \code{orf}
#'
#' @param x object of type \code{margins.orf}
#' @param ... further arguments (currently ignored)
#'
#' @export
print.orf <- function(x, ...) {

  # needed inputs for the function: forest - forest object coming from orf function
  #                                        - ... additional arguments (currently ignored)

  # -------------------------------------------------------------------------------- #

  ## get forest as x
  forest <- x

  # -------------------------------------------------------------------------------- #

  ## save forest inputs
  main_class        <- class(forest)[1]
  inputs            <- forest$forestInfo$inputs

  honesty           <- inputs$honesty
  mtry              <- inputs$mtry
  num.trees         <- inputs$num.trees
  min.node.size     <- inputs$min.node.size
  replace           <- inputs$replace
  inference         <- inputs$inference
  honest_data       <- forest$forestInfo$honestData
  train_data        <- forest$forestInfo$trainData
  categories        <- length(forest$forestInfo$categories)
  build             <- ifelse(replace == TRUE, "Bootstrap", "Subsampling")
  type              <- "Ordered Forest"

  # -------------------------------------------------------------------------------- #

  cat(type, "object of class", main_class, "\n\n")

  cat("Number of Categories:            ", categories, "\n")
  cat("Sample Size:                     ", sum(nrow(train_data), nrow(honest_data)), "\n")
  cat("Number of Trees:                 ", num.trees, "\n")
  cat("Build:                           ", build, "\n")
  cat("Mtry:                            ", mtry, "\n")
  cat("Minimum Node Size:               ", min.node.size, "\n")
  cat("Honest Forest:                   ", honesty, "\n")
  cat("Weight-Based Inference:          ", inference, "\n")

  # -------------------------------------------------------------------------------- #

}

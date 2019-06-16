#' orf
#'
#' Ordered Random Forests for semi-parametric estimation of ordered choice model
#' as proposed in Lechner & Okasa (2019)
#'
#' @param X matrix of features
#' @param Y vector of outcomes (as.matrix acceptable too)
#' @param ntree scalar, number of trees in a forest, i.e. bootstrap replications
#' @param mtry scalar, number of randomly selected features
#' @param nmin scalar, minimum node size
#' @param honesty logical, if TRUE honest forest is built using 50:50 data split
#' @param inference logical, if TRUE the weight based inference is conducted
#' @param margins logical, if TRUE marginal effects at mean are estimated
#'
#' @import ranger
#'
#' @return object of type orf
#'
#' @export
orf <- function(X, Y, ntree, mtry, nmin, honesty, inference, margins){

  # needed inputs for the function: X - matrix of features
  #                                 Y - vector of outcomes (as.matrix acceptable too)
  #                                 ntree - number of trees in a forest
  #                                 mtry - number of randomly selected features
  #                                 nmin - minimum node size
  #                                 honesty - logical, if TRUE honest forest is built using 50:50 data split
  #                                 inference - logical, if TRUE the weight based inference is conducted (honesty has to be TRUE)
  #                                 margins - logical, if TRUE the marginal effects at mean will be computed

  # -------------------------------------------------------------------------------- #

  ## standard checks for input data
  check_X(X)
  Y <- check_Y(Y, X)
  mtry <- check_mtry(mtry, X)
  nmin <- check_nmin(nmin, X)
  honesty <- check_honesty(honesty)
  inference <- check_inference(inference)

  ## check for plausibility of options first:
  if (honesty == FALSE & inference == TRUE) {

    warning("For conducting inference honesty is required. Honesty has been set to TRUE.")
    # set honesty to TRUE
    honesty <- TRUE

  }

  # --------------------------------------------------------------------------------------- #

  ## save the inputs:
  inputs <- list(ntree, mtry, nmin, honesty, inference, margins)
  names(inputs) <- c("ntree", "mtry", "nmin", "honesty", "inference", "margins")

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

  # --------------------------------------------------------------------------------------- #

  ## proceed with estimations
  # decide if honest forest should be estimated
  if (honesty == FALSE & inference == FALSE & margins == FALSE) {

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
                                                  num.trees = ntree, mtry = mtry, replace = TRUE,
                                                  min.node.size = nmin, importance = "none"))

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
    oob_mse <- mse(pred_final, Y)

    # compute OOB RPS based on whole sample
    oob_rps <- rps(pred_final, Y)

    # --------------------------------------------------------------------------------------- #

    ## convert probabilities into class predictions ("classification")
    pred_class <- as.matrix(apply(pred_final, 1, which.max))
    colnames(pred_class) <- "Category"

    # --------------------------------------------------------------------------------------- #

    # save forest information
    forest_info <- list(inputs, train_data, honest_data, categories, data_ind)
    names(forest_info) <- c("inputs", "trainData", "honestData", "categories", "indicatorData")

    # define output of the function
    output <- list(forest, forest_info, pred_final, pred_class, oob_mse, oob_rps)
    names(output) <- c("trainForests",  "forestInfo", "oobPredictions", "predictedCategories", "oobMSE", "oobRPS")

    # -------------------------------------------------------------------------------- #

  } else if (honesty == TRUE & inference == FALSE & margins == FALSE) {

    # --------------------------------------------------------------------------------------- #

    ## do honest forest estimation here using 50:50 data split as in Lechner (2018)
    # devide into 50:50 honesty sets
    split_data <- honest_split(dat)
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
                                                        num.trees = ntree, mtry = mtry, replace = FALSE,
                                                        sample.fraction = 0.5, min.node.size = nmin, importance = "none"))

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

    # compute OOB MSE based on whole sample
    honest_mse <- mse(pred_final, Y)

    # compute OOB RPS based on whole sample
    honest_rps <- rps(pred_final, Y)

    # --------------------------------------------------------------------------------------- #

    ## convert probabilities into class predictions ("classification")
    pred_class <- as.matrix(apply(pred_final, 1, which.max))
    colnames(pred_class) <- "Category"

    # --------------------------------------------------------------------------------------- #

    # save forest information
    forest_info <- list(inputs, train_data, honest_data, categories, data_ind_honest)
    names(forest_info) <- c("inputs", "trainData", "honestData", "categories", "indicatorData")

    # define output of the function
    output <- list(forest, forest_info, pred_final, pred_class, honest_mse, honest_rps)
    names(output) <- c("trainForests",  "forestInfo", "honestPredictions", "predictedCategories", "honestMSE", "honestRPS")

    # -------------------------------------------------------------------------------- #

  } else if (honesty == TRUE & inference == TRUE & margins == FALSE) {

    # --------------------------------------------------------------------------------------- #

    ## do honest forest estimation here using 50:50 data split as in Lechner (2018)
    # devide into 50:50 honesty sets
    split_data <- honest_split(dat)
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
                                                        num.trees = ntree, mtry = mtry, replace = FALSE,
                                                        sample.fraction = 0.5, min.node.size = nmin, importance = "none"))

    # --------------------------------------------------------------------------------------- #

    # get honest weights
    forest_weights <- mapply(function(x,y,z) get_forest_weights(x, y, z), forest, data_ind_honest, data_ind_train, SIMPLIFY = F)
    honest_weights <- lapply(forest_weights, function(x) x[rows_honest_data, ]) # take out honest sample honest weights
    train_weights <- lapply(forest_weights, function(x) x[rows_train_data, ]) # take out train sample honest weights

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

    # compute OOB MSE based on whole sample
    honest_mse <- mse(pred_final, Y)

    # compute OOB RPS based on whole sample
    honest_rps <- rps(pred_final, Y)

    # --------------------------------------------------------------------------------------- #

    ## convert probabilities into class predictions ("classification")
    pred_class <- as.matrix(apply(pred_final, 1, which.max))
    colnames(pred_class) <- "Category"

    # --------------------------------------------------------------------------------------- #

    # compute the variances for the categorical predictions
    var_final <- get_orf_variance(honest_pred, honest_weights, train_pred, train_weights, Y_ind_honest)
    # check for normalization
    #sd_final <- sqrt(var_final)

    # --------------------------------------------------------------------------------------- #

    # save forest information
    forest_info <- list(inputs, train_data, honest_data, categories, data_ind_honest)
    names(forest_info) <- c("inputs", "trainData", "honestData", "categories", "indicatorData")

    # define output of the function
    output <- list(forest, forest_info, pred_final, var_final, pred_class, honest_mse, honest_rps)
    names(output) <- c("trainForests",  "forestInfo", "honestPredictions", "honestVariance", "predictedCategories", "honestMSE", "honestRPS")

    # --------------------------------------------------------------------------------------- #

  } else if (honesty == FALSE & inference == FALSE & margins == TRUE) {

    # --------------------------------------------------------------------------------------- #

    # put warning that marginal effects are preffered with honesty
    warning("Estimation of marginal effects without honesty might not be optimal.")

    # --------------------------------------------------------------------------------------- #

    # do the estimation of probabilities as usual

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
                                                  num.trees = ntree, mtry = mtry, replace = TRUE,
                                                  min.node.size = nmin, importance = "none"))

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
    oob_mse <- mse(pred_final, Y)

    # compute OOB RPS based on whole sample
    oob_rps <- rps(pred_final, Y)

    # --------------------------------------------------------------------------------------- #

    ## convert probabilities into class predictions ("classification")
    pred_class <- as.matrix(apply(pred_final, 1, which.max))
    colnames(pred_class) <- "Category"

    # --------------------------------------------------------------------------------------- #

    ## now compute marginal effects at mean (without honesty and without inference - not optimal!)

    # --------------------------------------------------------------------------------------- #

    marginal_effects <- orf_margins(forest, data_ind, honesty, inference)

    # --------------------------------------------------------------------------------------- #

    # save forest information
    forest_info <- list(inputs, train_data, honest_data, categories, data_ind)
    names(forest_info) <- c("inputs", "trainData", "honestData", "categories", "indicatorData")

    # define output of the function
    output <- list(forest, forest_info, pred_final, pred_class, oob_mse, oob_rps, marginal_effects)
    names(output) <- c("trainForests",  "forestInfo", "oobPredictions", "predictedCategories", "oobMSE", "oobRPS", "marginalEffects")

    # -------------------------------------------------------------------------------- #

  } else if (honesty == TRUE & inference == FALSE & margins == TRUE) {

    # --------------------------------------------------------------------------------------- #

    # do the estimation of probabilities as usual

    # --------------------------------------------------------------------------------------- #

    ## do honest forest estimation here using 50:50 data split as in Lechner (2018)
    # devide into 50:50 honesty sets
    split_data <- honest_split(dat)
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
                                                        num.trees = ntree, mtry = mtry, replace = FALSE,
                                                        sample.fraction = 0.5, min.node.size = nmin, importance = "none"))

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

    # compute OOB MSE based on whole sample
    honest_mse <- mse(pred_final, Y)

    # compute OOB RPS based on whole sample
    honest_rps <- rps(pred_final, Y)

    # --------------------------------------------------------------------------------------- #

    ## convert probabilities into class predictions ("classification")
    pred_class <- as.matrix(apply(pred_final, 1, which.max))
    colnames(pred_class) <- "Category"

    # --------------------------------------------------------------------------------------- #

    ## now compute marginal effects at mean (with honesty and without inference)

    # --------------------------------------------------------------------------------------- #

    marginal_effects <- orf_margins(forest, data_ind_honest, honesty, inference)

    # --------------------------------------------------------------------------------------- #

    # save forest information
    forest_info <- list(inputs, train_data, honest_data, categories, data_ind_honest)
    names(forest_info) <- c("inputs", "trainData", "honestData", "categories", "indicatorData")

    # define output of the function
    output <- list(forest, forest_info, pred_final, pred_class, honest_mse, honest_rps, marginal_effects)
    names(output) <- c("trainForests",  "forestInfo", "honestPredictions", "predictedCategories", "honestMSE", "honestRPS", "marginalEffects")

    # -------------------------------------------------------------------------------- #

  } else if (honesty == TRUE & inference == TRUE & margins == TRUE) {

    # --------------------------------------------------------------------------------------- #

    ## do honest forest estimation here using 50:50 data split as in Lechner (2018)
    # devide into 50:50 honesty sets
    split_data <- honest_split(dat)
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
                                                        num.trees = ntree, mtry = mtry, replace = FALSE,
                                                        sample.fraction = 0.5, min.node.size = nmin, importance = "none"))

    # --------------------------------------------------------------------------------------- #

    # get honest weights
    forest_weights <- mapply(function(x,y,z) get_forest_weights(x, y, z), forest, data_ind_honest, data_ind_train, SIMPLIFY = F)
    honest_weights <- lapply(forest_weights, function(x) x[rows_honest_data, ]) # take out honest sample honest weights
    train_weights <- lapply(forest_weights, function(x) x[rows_train_data, ]) # take out train sample honest weights

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

    # compute OOB MSE based on whole sample
    honest_mse <- mse(pred_final, Y)

    # compute OOB RPS based on whole sample
    honest_rps <- rps(pred_final, Y)

    # --------------------------------------------------------------------------------------- #

    ## convert probabilities into class predictions ("classification")
    pred_class <- as.matrix(apply(pred_final, 1, which.max))
    colnames(pred_class) <- "Category"

    # --------------------------------------------------------------------------------------- #

    # compute the variances for the categorical predictions
    var_final <- get_orf_variance(honest_pred, honest_weights, train_pred, train_weights, Y_ind_honest)

    # --------------------------------------------------------------------------------------- #

    ## now compute marginal effects at mean (with honesty and with inference - default)

    # --------------------------------------------------------------------------------------- #

    marginal_effects <- orf_margins(forest, data_ind_honest, honesty, inference)

    # --------------------------------------------------------------------------------------- #

    # save forest information
    forest_info <- list(inputs, train_data, honest_data, categories, data_ind_honest)
    names(forest_info) <- c("inputs", "trainData", "honestData", "categories", "indicatorData")

    # define output of the function
    output <- list(forest, forest_info, pred_final, var_final, pred_class, honest_mse, honest_rps, marginal_effects)
    names(output) <- c("trainForests",  "forestInfo", "honestPredictions", "honestVariance", "predictedCategories", "honestMSE", "honestRPS", "marginalEffects")

  }

  # -------------------------------------------------------------------------------- #

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
#' @param new_data matrix X containing the observations to predict
#' @param ... further arguments (currently ignored)
#'
#' @import ranger
#'
#' @return object of class \code{orf.prediction} with elements
##'   \tabular{ll}{
##'       \code{trainForests}  \tab saved forests trained for ORF estimations (inherited from \code{ranger}) \cr
##'       \code{forestInfo}    \tab info containing forest inputs and data used \cr
##'       \code{forestPredictions} \tab predicted values \cr
##'       \code{forestVariances} \tab variances of predicted values (only if \code{inference=TRUE} in the passed \code{orf object}) \cr
##'       \code{predictedCategories} \tab categories predicted with ORF based on maximum probability \cr
##'       \code{marginalEffects} \tab estimated marginal effects from ORF (optionally with inference) \cr
##'   }
#'
#' @export
predict.orf <- function(object, new_data, ...) {
  # needed inputs for the function: forest - orf object coming from orf function
  #                                 new_data - matrix X containing the observations to predict

  # -------------------------------------------------------------------------------- #

  ## standard checks for input data
  if (class(object) != "orf") {
    stop("Forest object is not of class orf. Programme temrinated.")
  }

  ## get forest as na object
  forest <- object
  ## save forest inputs
  inputs <- forest$forestInfo$inputs
  categories <- forest$forestInfo$categories
  honesty <- inputs$honesty
  inference <- inputs$inference
  margins <- inputs$margins
  honest_data <- forest$forestInfo$honestData
  train_data <- forest$forestInfo$trainData
  honest_ind_data <- forest$forestInfo$indicatorData # indicator data needed for indicator predictions

  # take out list of ranger objects (be careful, its Forests with S at the end!)
  forest <- forest$trainForests

  ## get train data names (only X)
  train_data_name <- colnames(train_data)[2:ncol(train_data)]

  ## get X matrix as dataframe and check colnames
  # X
  if (is.null(colnames(new_data))) { colnames(new_data) <- paste0("X", rep(1:ncol(new_data))) } # check if X has name
  new_data_name <- colnames(new_data) # save the name of X
  new_data <- as.data.frame(new_data) # as dataframe
  n_new_data <- nrow(new_data) # rows of new data
  n_cat <- as.numeric(length(categories))

  # -------------------------------------------------------------------------------- #

  # check if its compatible with the data used for training
  if (all(train_data_name != new_data_name) | (ncol(new_data) != ncol(train_data)-1)) {

    stop("New data are not compatible with the training data. Check supplied data. Program terminated.")

  }

  # -------------------------------------------------------------------------------- #

  # check if honest forest was estimated and predict accordingly
  if (honesty == FALSE & inference == FALSE & margins == FALSE) {

    ## no honest splitting, i.e. use all data, no inference and no margins

    # predict with ncat-1 forests as in ranger default
    pred <- lapply(forest, function(x) predict(x, data = new_data)$predictions)

    # -------------------------------------------------------------------------------- #

    # add the probability for the last outcome (always 1)
    pred_1 <- append(pred, list(rep(1, n_new_data)))
    # prepend zero vector to predictions for later differencing
    pred_0 <- append(list(rep(0, n_new_data)), pred) # append a first 0 elemnt for the list

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

    ## convert probabilities into class predictions ("classification")
    pred_class <- as.matrix(apply(pred_final, 1, which.max))
    colnames(pred_class) <- "Category"

    # --------------------------------------------------------------------------------------- #

    # save forest information
    forest_info <- list(inputs, categories, new_data)
    names(forest_info) <- c("inputs", "categories","newData")

    # define output of the function
    output <- list(forest_info, pred_final, pred_class)
    names(output) <- c("forestInfo", "forestPredictions", "predictedCategories")

    # -------------------------------------------------------------------------------- #

  } else if (honesty == TRUE & inference == FALSE & margins == FALSE) {

    # -------------------------------------------------------------------------------- #

    ## run new Xs through estimated train forest and compute predictions based on honest sample
    # no need to predict weights, get predictions directly through leaves
    # predict with ncat-1 forests
    pred <- mapply(function(x,y) predict_honest(x, y, new_data), forest, honest_ind_data, SIMPLIFY = FALSE)

    # -------------------------------------------------------------------------------- #

    # add the probability for the last outcome (always 1)
    pred_1 <- append(pred, list(rep(1, n_new_data)))
    # prepend zero vector to predictions for later differencing
    pred_0 <- append(list(rep(0, n_new_data)), pred) # append a first 0 elemnt for the list

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

    ## convert probabilities into class predictions ("classification")
    pred_class <- as.matrix(apply(pred_final, 1, which.max))
    colnames(pred_class) <- "Category"

    # --------------------------------------------------------------------------------------- #

    # save forest information
    forest_info <- list(inputs, categories, new_data)
    names(forest_info) <- c("inputs", "categories","newData")

    # define output of the function
    output <- list(forest_info, pred_final, pred_class)
    names(output) <- c("forestInfo", "forestPredictions", "predictedCategories")

    # -------------------------------------------------------------------------------- #

  } else if (honesty == TRUE & inference == TRUE & margins == FALSE) {

    # -------------------------------------------------------------------------------- #

    # predict weights by using forest train, honest data and new_data (for each category except one)
    forest_weights_pred <- lapply(forest, function(x) predict_forest_weights(x, honest_data, new_data))

    # get predictions by matrix multiplication of weights with honest responses
    forest_pred <- mapply(function(x,y) as.numeric(x%*%y[, 1]), forest_weights_pred, honest_ind_data, SIMPLIFY = FALSE)

    # -------------------------------------------------------------------------------- #

    # add the probability for the last outcome (always 1)
    pred_1 <- append(forest_pred, list(rep(1, n_new_data)))
    # prepend zero vector to predictions for later differencing
    pred_0 <- append(list(rep(0, n_new_data)), forest_pred) # append a first 0 elemnt for the list

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

    ## convert probabilities into class predictions ("classification")
    pred_class <- as.matrix(apply(pred_final, 1, which.max))
    colnames(pred_class) <- "Category"

    # --------------------------------------------------------------------------------------- #

    # get indicator outcomes out of honest indicator data
    Y_ind_honest <- lapply(honest_ind_data, function(x) x[, 1])

    # compute the variances for the categorical predictions
    var_final <- pred_orf_variance(forest_pred, forest_weights_pred, Y_ind_honest)

    # --------------------------------------------------------------------------------------- #

    # save forest information
    forest_info <- list(inputs, categories, new_data)
    names(forest_info) <- c("inputs", "categories", "newData")

    # define output of the function
    output <- list(forest_info, pred_final, var_final, pred_class)
    names(output) <- c("forestInfo", "forestPredictions", "forestVariances", "predictedCategories")

    # -------------------------------------------------------------------------------- #

  } else if (honesty == FALSE & inference == FALSE & margins == TRUE) {

    # -------------------------------------------------------------------------------- #

    ## no honest splitting, i.e. use all data, no inference and no margins
    # predict with ncat-1 forests as in ranger default
    pred <- lapply(forest, function(x) predict(x, data = new_data)$predictions)

    # -------------------------------------------------------------------------------- #

    # add the probability for the last outcome (always 1)
    pred_1 <- append(pred, list(rep(1, n_new_data)))
    # prepend zero vector to predictions for later differencing
    pred_0 <- append(list(rep(0, n_new_data)), pred) # append a first 0 elemnt for the list

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

    # -------------------------------------------------------------------------------- #

    ## convert probabilities into class predictions ("classification")
    pred_class <- as.matrix(apply(pred_final, 1, which.max))
    colnames(pred_class) <- "Category"

    # --------------------------------------------------------------------------------------- #

    ## now predict marginal effects at mean (without honesty and without inference - not optimal!)

    # --------------------------------------------------------------------------------------- #

    marginal_effects <- pred_orf_margins(forest, honest_ind_data, new_data, honesty, inference)

    # --------------------------------------------------------------------------------------- #

    # save forest information
    forest_info <- list(inputs, train_data, honest_data, categories)
    names(forest_info) <- c("inputs", "trainData", "honestData", "categories")

    # define output of the function
    output <- list(forest, forest_info, pred_final, pred_class, marginal_effects)
    names(output) <- c("trainForests",  "forestInfo", "oobPredictions", "predictedCategories", "marginalEffects")

    # --------------------------------------------------------------------------------------- #

  } else if (honesty == TRUE & inference == FALSE & margins == TRUE) {

    # -------------------------------------------------------------------------------- #

    ## run new Xs through estimated train forest and compute predictions based on honest sample
    # no need to predict weights, get predictions directly through leaves
    # predict with ncat-1 forests
    pred <- mapply(function(x,y) predict_honest(x, y, new_data), forest, honest_ind_data, SIMPLIFY = FALSE)

    # -------------------------------------------------------------------------------- #

    # add the probability for the last outcome (always 1)
    pred_1 <- append(pred, list(rep(1, n_new_data)))
    # prepend zero vector to predictions for later differencing
    pred_0 <- append(list(rep(0, n_new_data)), pred) # append a first 0 elemnt for the list

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

    ## convert probabilities into class predictions ("classification")
    pred_class <- as.matrix(apply(pred_final, 1, which.max))
    colnames(pred_class) <- "Category"

    # --------------------------------------------------------------------------------------- #

    ## now compute marginal effects at mean (with honesty and without inference)

    # --------------------------------------------------------------------------------------- #

    marginal_effects <- pred_orf_margins(forest, honest_ind_data, new_data, honesty, inference)

    # --------------------------------------------------------------------------------------- #

    # save forest information
    forest_info <- list(inputs, train_data, honest_data, categories)
    names(forest_info) <- c("inputs", "trainData", "honestData", "categories")

    # define output of the function
    output <- list(forest, forest_info, pred_final, pred_class, marginal_effects)
    names(output) <- c("trainForests",  "forestInfo", "honestPredictions", "predictedCategories", "marginalEffects")

    # -------------------------------------------------------------------------------- #

  } else if (honesty == TRUE & inference == TRUE & margins == TRUE) {

    # -------------------------------------------------------------------------------- #

    # predict weights by using forest train, honest data and new_data (for each category except one)
    forest_weights_pred <- lapply(forest, function(x) predict_forest_weights(x, honest_data, new_data))

    # get predictions by matrix multiplication of weights with honest responses
    forest_pred <- mapply(function(x,y) as.numeric(x%*%y[, 1]), forest_weights_pred, honest_ind_data, SIMPLIFY = FALSE)

    # -------------------------------------------------------------------------------- #

    # add the probability for the last outcome (always 1)
    pred_1 <- append(forest_pred, list(rep(1, n_new_data)))
    # prepend zero vector to predictions for later differencing
    pred_0 <- append(list(rep(0, n_new_data)), forest_pred) # append a first 0 elemnt for the list

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

    ## convert probabilities into class predictions ("classification")
    pred_class <- as.matrix(apply(pred_final, 1, which.max))
    colnames(pred_class) <- "Category"

    # --------------------------------------------------------------------------------------- #

    # get indicator outcomes out of honest indicator data
    Y_ind_honest <- lapply(honest_ind_data, function(x) x[, 1])

    # compute the variances for the categorical predictions
    var_final <- pred_orf_variance(forest_pred, forest_weights_pred, Y_ind_honest)

    # --------------------------------------------------------------------------------------- #

    ## now compute marginal effects at mean (with honesty and with inference - default)

    # --------------------------------------------------------------------------------------- #

    marginal_effects <- pred_orf_margins(forest, honest_ind_data, new_data, honesty, inference)

    # --------------------------------------------------------------------------------------- #

    # save forest information
    forest_info <- list(inputs, train_data, honest_data, categories)
    names(forest_info) <- c("inputs", "trainData", "honestData", "categories")

    # define output of the function
    output <- list(forest, forest_info, pred_final, var_final, pred_class, marginal_effects)
    names(output) <- c("trainForests",  "forestInfo", "honestPredictions", "honestVariances", "predictedCategories", "marginalEffects")

    # -------------------------------------------------------------------------------- #

  }

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
#' @importFrom gridExtra grid.arrange
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

  # chekc if honesty was used or not and modify objects accordingly
  if (honesty == FALSE) {

    probabilities <- forest$oobPredictions # take out probabilities
    outcomes <- train_data[, 1] # take the observed outcomes


  } else {

    probabilities <- forest$honestPredictions # take out honest predictions
    all_data <- rbind(honest_data, train_data) # put data together
    all_data <- all_data[order(as.numeric(row.names(all_data))), ] # sort data as original
    outcomes <- all_data[, 1] # take the observed outcomes


  }

  # -------------------------------------------------------------------------------- #

  ### plot ORF ###

  ## plot realized categories overlayed with predicted category probabilities
  # new colnames
  colnames(probabilities) <- sapply(seq_along(categories), function(i) paste0("P(Y=", i, ")"))

  # cbind together
  df_plot <- cbind(outcomes, probabilities)

  # subset according to categories
  df_plot_cat <- lapply(seq_along(categories), function(i) as.data.frame(subset(df_plot, outcomes == i)))
  # take colmeans
  df_cat_means <- lapply(df_plot_cat, function(x) t(as.matrix(colMeans(x)[seq_along(categories)+1])))
  # add colmeans to df_plot_cat
  df_plot_cat <- mapply(function(x,y) cbind(x, y), df_plot_cat, df_cat_means, SIMPLIFY = FALSE)

  # reshape data for ggplot
  df_plot_prob <- lapply(df_plot_cat, function(x) stack(x[seq_along(categories)+1]))
  df_plot_mean <- lapply(df_plot_cat, function(x) stack(x[seq_along(categories)+1+length(categories)]))

  # add colnames
  df_plot_prob <- lapply(df_plot_prob, function(x) { colnames(x) <- c("Probability", "Category")
  return(x)}
  )
  # add colnames
  df_plot_mean <- lapply(df_plot_mean, function(x) { colnames(x) <- c("Probability", "Category")
  return(x)}
  )

  # -------------------------------------------------------------------------------- #

  ## generate plot objects
  # Use semi-transparent fill
  plots <- lapply(seq_along(df_plot_prob), function(i) {
    ggplot(df_plot_prob[[i]], aes_string(x="Probability", fill="Category")) +
      geom_density(alpha=0.4) +
      geom_vline(data=df_plot_mean[[i]], aes_string(xintercept="Probability", color="Category"), linetype="dashed") +
      ggtitle(paste("Category", i, sep = " ")) +
      xlab("Predicted Probability") +
      ylab("Probability Mass") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5))
  })

  # print plots
  do.call("grid.arrange", c(plots, ncol=1))

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
  ## save forest inputs
  inputs <- forest$forestInfo$inputs
  honesty <- inputs$honesty
  inference <- inputs$inference
  margins <- inputs$margins
  mtry <- inputs$mtry
  ntree <- inputs$ntree
  nmin <- inputs$nmin
  honest_data <- forest$forestInfo$honestData
  train_data <- forest$forestInfo$trainData
  categories <- length(forest$forestInfo$categories)
  type <- "Ordered Random Forest"

  # -------------------------------------------------------------------------------- #

  # check if honest forest was estimated and predict accordingly
  if (honesty == TRUE & inference == TRUE & margins == TRUE) {

    # -------------------------------------------------------------------------------- #

    ## honest splitting, i.e. use honest data
    # take out summary statistics
    mse <- round(forest$honestMSE, 5)
    rps <- round(forest$honestRPS, 5)
    trainsize <- nrow(train_data)
    honestsize <- nrow(honest_data)
    features <- ncol(train_data)-1 # take out the response
    # check if subsampling or bootstrapping was used
    if (forest$trainForests[[1]]$replace == TRUE) { build <- "Bootstrap" } else { build <- "Subsampling" }

    # -------------------------------------------------------------------------------- #

    # structure summary into a list
    output <- list(type, categories, build, ntree, mtry, nmin, honesty, inference, margins, trainsize, honestsize, features, mse, rps)
    names(output) <- c("type", "categories", "build", "ntree", "mtry", "nmin", "honesty", "inference", "margins", "trainsize", "honestsize", "features", "mse", "rps")

    # output matrix
    output_matrix <- matrix(NA, ncol = 1, nrow = length(output))
    # populate output matrix
    rownames(output_matrix) <- names(output) # rownames are names
    colnames(output_matrix) <- "" # no colname
    output_matrix[, 1] <- unlist(output) # column 2 are values

    # generate latex output if selected
    if (latex == TRUE) { colnames(output_matrix) <- "Attributes"
    output_matrix <- xtable(output_matrix, caption = "Ordered Random Forest Summary")
    }

    # -------------------------------------------------------------------------------- #

    ## summary of marginal effects
    # take out coefficients with their p values
    coefs <- forest$marginalEffects$MarginalEffects
    pvalues <- forest$marginalEffects$pValues

    # produce nice table (latex or no latex)
    coef_table <- coefstars(coefs, pvalues)

    if (latex == TRUE) {
      coef_table <- xtable(coef_table, caption = "ORF Marginal Effects")
    }

    # pack it into output
    output <- list(output_matrix, coef_table)
    names(output) <- c("ForestSummary", "MarginalEffects")

    # -------------------------------------------------------------------------------- #

  } else if (honesty == TRUE & inference == TRUE & margins == FALSE) {

    # -------------------------------------------------------------------------------- #

    ## honest splitting, i.e. use honest data
    # take out summary statistics
    mse <- round(forest$honestMSE, 5)
    rps <- round(forest$honestRPS, 5)
    trainsize <- nrow(train_data)
    honestsize <- nrow(honest_data)
    features <- ncol(train_data)-1 # take out the response
    # check if subsampling or bootstrapping was used
    if (forest$trainForests[[1]]$replace == TRUE) { build <- "Bootstrap" } else { build <- "Subsampling" }

    # -------------------------------------------------------------------------------- #

    # structure summary into a list
    output <- list(type, categories, build, ntree, mtry, nmin, honesty, inference, margins, trainsize, honestsize, features, mse, rps)
    names(output) <- c("type", "categories", "build", "ntree", "mtry", "nmin", "honesty", "inference", "margins", "trainsize", "honestsize", "features", "mse", "rps")

    # output matrix
    output_matrix <- matrix(NA, ncol = 1, nrow = length(output))
    # populate output matrix
    rownames(output_matrix) <- names(output) # rownames are names
    colnames(output_matrix) <- "" # no colname
    output_matrix[, 1] <- unlist(output) # column 2 are values

    # generate latex output if selected
    if (latex == TRUE) { colnames(output_matrix) <- "Attributes"
    output_matrix <- xtable(output_matrix, caption = "Ordered Random Forest Summary")
    }

    # pack it into output
    output <- output_matrix

    # -------------------------------------------------------------------------------- #

  } else if (honesty == FALSE & inference == FALSE & margins == FALSE) {

    # -------------------------------------------------------------------------------- #

    ## no honest splitting, i.e. use all data
    # take out summary statistics
    mse <- round(forest$oobMSE, 5)
    rps <- round(forest$oobRPS, 5)
    trainsize <- nrow(train_data)
    honestsize <- 0
    features <- ncol(train_data)-1 # take out the response
    # check if subsampling or bootstrapping was used
    if (forest$trainForests[[1]]$replace == TRUE) { build <- "Bootstrap" } else { build <- "Subsampling" }

    # -------------------------------------------------------------------------------- #

    # structure summary into a list
    output <- list(type, categories, build, ntree, mtry, nmin, honesty, inference, margins, trainsize, honestsize, features, mse, rps)
    names(output) <- c("type", "categories", "build", "ntree", "mtry", "nmin", "honesty", "inference", "margins", "trainsize", "honestsize", "features", "mse", "rps")

    # output matrix
    output_matrix <- matrix(NA, ncol = 1, nrow = length(output))
    # populate output matrix
    rownames(output_matrix) <- names(output) # rownames are names
    colnames(output_matrix) <- "" # no colname
    output_matrix[, 1] <- unlist(output) # column 2 are values

    # generate latex output if selected
    if (latex == TRUE) { colnames(output_matrix) <- "Attributes"
    output_matrix <- xtable(output_matrix, caption = "Ordered Random Forest Summary")
    }

    # put it into output
    output <- output_matrix

    # -------------------------------------------------------------------------------- #

  } else if (honesty == FALSE & inference == FALSE & margins == TRUE) {

    # -------------------------------------------------------------------------------- #

    ## no honest splitting, i.e. use all data
    # take out summary statistics
    mse <- round(forest$oobMSE, 5)
    rps <- round(forest$oobRPS, 5)
    trainsize <- nrow(train_data)
    honestsize <- 0
    features <- ncol(train_data)-1 # take out the response
    # check if subsampling or bootstrapping was used
    if (forest$trainForests[[1]]$replace == TRUE) { build <- "Bootstrap" } else { build <- "Subsampling" }

    # -------------------------------------------------------------------------------- #

    # structure summary into a list
    output <- list(type, categories, build, ntree, mtry, nmin, honesty, inference, margins, trainsize, honestsize, features, mse, rps)
    names(output) <- c("type", "categories", "build", "ntree", "mtry", "nmin", "honesty", "inference", "margins", "trainsize", "honestsize", "features", "mse", "rps")

    # output matrix
    output_matrix <- matrix(NA, ncol = 1, nrow = length(output))
    # populate output matrix
    rownames(output_matrix) <- names(output) # rownames are names
    colnames(output_matrix) <- "" # no colname
    output_matrix[, 1] <- unlist(output) # column 2 are values

    # generate latex output if selected
    if (latex == TRUE) { colnames(output_matrix) <- "Attributes"
    output_matrix <- xtable(output_matrix, caption = "Ordered Random Forest Summary")
    }

    # -------------------------------------------------------------------------------- #

    ## summary of marginal effects
    # take out coefficients with their p values
    coef_table <- round(forest$marginalEffects, 3)

    # chekc latex
    if (latex == TRUE) {
      coef_table <- xtable(coef_table, caption = "ORF Marginal Effects", digits = 3)
    }

    # pack it into output
    output <- list(output_matrix, coef_table)
    names(output) <- c("ForestSummary", "MarginalEffects")

    # -------------------------------------------------------------------------------- #

  } else if (honesty == TRUE & inference == FALSE & margins == TRUE) {

    # -------------------------------------------------------------------------------- #

    ## honest splitting, i.e. use honest data
    # take out summary statistics
    mse <- round(forest$honestMSE, 5)
    rps <- round(forest$honestRPS, 5)
    trainsize <- nrow(train_data)
    honestsize <- nrow(honest_data)
    features <- ncol(train_data)-1 # take out the response
    # check if subsampling or bootstrapping was used
    if (forest$trainForests[[1]]$replace == TRUE) { build <- "Bootstrap" } else { build <- "Subsampling" }

    # -------------------------------------------------------------------------------- #

    # structure summary into a list
    output <- list(type, categories, build, ntree, mtry, nmin, honesty, inference, margins, trainsize, honestsize, features, mse, rps)
    names(output) <- c("type", "categories", "build", "ntree", "mtry", "nmin", "honesty", "inference", "margins", "trainsize", "honestsize", "features", "mse", "rps")

    # output matrix
    output_matrix <- matrix(NA, ncol = 1, nrow = length(output))
    # populate output matrix
    rownames(output_matrix) <- names(output) # rownames are names
    colnames(output_matrix) <- "" # no colname
    output_matrix[, 1] <- unlist(output) # column 2 are values

    # generate latex output if selected
    if (latex == TRUE) { colnames(output_matrix) <- "Attributes"
    output_matrix <- xtable(output_matrix, caption = "Ordered Random Forest Summary")
    }

    # -------------------------------------------------------------------------------- #

    ## summary of marginal effects
    # take out coefficients without their p values (no pvalues if no inference)
    coef_table <- round(forest$marginalEffects, 3)

    if (latex == TRUE) {
      coef_table <- xtable(coef_table, caption = "ORF Marginal Effects", digits = 3)
    }

    # pack it into output
    output <- list(output_matrix, coef_table)
    names(output) <- c("ForestSummary", "MarginalEffects")

    # -------------------------------------------------------------------------------- #

  } else if (honesty == TRUE & inference == FALSE & margins == FALSE) {

    # -------------------------------------------------------------------------------- #

    ## honest splitting, i.e. use honest data
    # take out summary statistics
    mse <- round(forest$honestMSE, 5)
    rps <- round(forest$honestRPS, 5)
    trainsize <- nrow(train_data)
    honestsize <- nrow(honest_data)
    features <- ncol(train_data)-1 # take out the response
    # check if subsampling or bootstrapping was used
    if (forest$trainForests[[1]]$replace == TRUE) { build <- "Bootstrap" } else { build <- "Subsampling" }

    # -------------------------------------------------------------------------------- #

    # structure summary into a list
    output <- list(type, categories, build, ntree, mtry, nmin, honesty, inference, margins, trainsize, honestsize, features, mse, rps)
    names(output) <- c("type", "categories", "build", "ntree", "mtry", "nmin", "honesty", "inference", "margins", "trainsize", "honestsize", "features", "mse", "rps")

    # output matrix
    output_matrix <- matrix(NA, ncol = 1, nrow = length(output))
    # populate output matrix
    rownames(output_matrix) <- names(output) # rownames are names
    colnames(output_matrix) <- "" # no colname
    output_matrix[, 1] <- unlist(output) # column 2 are values

    # generate latex output if selected
    if (latex == TRUE) { colnames(output_matrix) <- "Attributes"
    output_matrix <- xtable(output_matrix, caption = "Ordered Random Forest Summary")
    }

    # -------------------------------------------------------------------------------- #

    # pack it into output
    output <- output_matrix

    # -------------------------------------------------------------------------------- #

  }

  # return output
  return(output)

}

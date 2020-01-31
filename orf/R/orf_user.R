#' Ordered Forest Estimator
#'
#' An implementation of the Ordered Forest estimator as developed
#' in Lechner & Okasa (2019). The Ordered Forest flexibly
#' estimates the conditional probabilities of models with ordered
#' categorical outcomes (so-called ordered choice models).
#' Additionally to common machine learning algorithms the \code{orf}
#' package provides functions for estimating marginal effects as well
#' as statistical inference thereof and thus provides similar output
#' as in standard econometric models for ordered choice. The core
#' forest algorithm relies on the fast C++ forest implementation
#' from the \code{ranger} package (Wright & Ziegler, 2017).
#'
#' The Ordered Forest function, \code{orf}, estimates the conditional ordered choice
#' probabilities, i.e. P[Y=m|X=x]. Additionally, weight-based inference for
#' the probability predictions can be conducted as well. If inference is desired,
#' the Ordered Forest must be estimated with honesty and subsampling.
#' If prediction only is desired, estimation without honesty and with bootstrapping
#' is recommended for optimal prediction performance.
#'
#' In order to estimate the Ordered Forest user must supply the data in form of
#' matrix of covariates \code{X} and a vector of outcomes 'code{Y} to the \code{orf}
#' function. These data inputs are also the only inputs that must be specified by
#' the user without any defaults. Further optional arguments include the classical forest
#' hyperparameters such as number of trees, \code{num.trees}, number of randomly
#' selected features, \code{mtry}, and the minimum leaf size, \code{min.node.size}.
#' The forest building scheme is regulated by the \code{replace} argument, meaning
#' bootstrapping if \code{replace = TRUE} or subsampling if \code{replace = FALSE}.
#' For the case of subsampling, \code{sample.fraction} argument regulates the subsampling
#' rate. Further, honest forest is estimated if the \code{honesty} argument is set to
#' \code{TRUE}, which is also the default. Similarly, the fraction of the sample used
#' for the honest estimation is regulated by the \code{honesty.fraction} argument.
#' The default setting conducts a 50:50 sample split, which is also generally advised
#' to follow for optimal performance. Inference procedure of the Ordered Forest is based on
#' the forest weights and is controlled by the \code{inference} argument. Note, that
#' such weight-based inference is computationally demanding exercise due to the estimation
#' of the forest weights and as such longer computation time is to be expected. Lastly,
#' the \code{importance} argument turns on and off the permutation based variable
#' importance.
#'
#' \code{orf} is compatible with standard \code{R} commands such as
#' \code{predict}, \code{margins}, \code{plot}, \code{summary} and \code{print}.
#' For further details, see examples below.
#'
#' @seealso \code{\link{summary.orf}}, \code{\link{plot.orf}}
#' \code{\link{predict.orf}}, \code{\link{margins.orf}}
#'
#' @param X numeric matrix of features
#' @param Y numeric vector of outcomes
#' @param num.trees scalar, number of trees in a forest, i.e. bootstrap replications (default is 1000 trees)
#' @param mtry scalar, number of randomly selected features (default is the squared root of number of features, rounded up to the nearest integer)
#' @param min.node.size scalar, minimum node size, i.e. leaf size of a tree (default is 5 observations)
#' @param replace logical, if TRUE sampling with replacement, i.e. bootstrap is used to grow the trees, otherwise subsampling without replacement is used (default is set to FALSE)
#' @param sample.fraction scalar, subsampling rate (default is 1 for bootstrap and 0.5 for subsampling)
#' @param honesty logical, if TRUE honest forest is built using sample splitting (default is set to TRUE)
#' @param honesty.fraction scalar, share of observations belonging to honest sample not used for growing the forest (default is 0.5)
#' @param inference logical, if TRUE the weight based inference is conducted (default is set to FALSE)
#' @param importance logical, if TRUE variable importance measure based on permutation is conducted (default is set to FALSE)
#'
#' @import ranger
#'
#' @return object of type \code{orf} with following elements
#'       \item{forests}{saved forests trained for \code{orf} estimations (inherited from \code{ranger})}
#'       \item{info}{info containing forest inputs and data used}
#'       \item{predictions}{predicted values for class probabilities}
#'       \item{variances}{variances of predicted values}
#'       \item{importance}{weighted measure of permutation based variable importance}
#'       \item{accuracy}{oob measures for mean squared error and ranked probability score}
#'
#' @author Gabriel Okasa
#'
#' @references
#' \itemize{
#'   \item Lechner, M., & Okasa, G. (2019). Random Forest Estimation of the Ordered Choice Model. arXiv preprint arXiv:1907.02436. \url{https://arxiv.org/abs/1907.02436}
#'   \item Goller, D., Knaus, M. C., Lechner, M., & Okasa, G. (2018). Predicting Match Outcomes in Football by an Ordered Forest Estimator (No. 1811). University of St. Gallen, School of Economics and Political Science. \url{http://ux-tauri.unisg.ch/RePEc/usg/econwp/EWP-1811.pdf}
#'   \item Wright, M. N. & Ziegler, A. (2017). ranger: A fast implementation of random forests for high dimensional data in C++ and R. J Stat Softw 77:1-17. \url{https://doi.org/10.18637/jss.v077.i01}.
#' }
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
#' # estimate Ordered Forest with default parameters
#' orf_fit <- orf(X, Y)
#' \donttest{
#' # estimate Ordered Forest with own tuning parameters
#' orf_fit <- orf(X, Y, num.trees = 2000, mtry = 3, min.node.size = 10)
#'
#' # estimate Ordered Forest with bootstrapping and without honesty
#' orf_fit <- orf(X, Y, replace = TRUE, honesty = FALSE)
#'
#' # estimate Ordered Forest with subsampling and with honesty
#' orf_fit <- orf(X, Y, replace = FALSE, honesty = TRUE)
#'
#' # estimate Ordered Forest with subsampling and with honesty
#' # with own tuning for subsample fraction and honesty fraction
#' orf_fit <- orf(X, Y, replace = FALSE, sample.fraction = 0.5,
#'                      honesty = TRUE, honesty.fraction = 0.5)
#'
#' # estimate Ordered Forest with subsampling and with honesty and with inference
#' # (for inference, subsampling and honesty are required)
#' orf_fit <- orf(X, Y, replace = FALSE, honesty = TRUE, inference = TRUE)
#'
#' # estimate Ordered Forest with simple variable importance measure
#' orf_fit <- orf(X, Y, importance = TRUE)
#'
#' # estimate Ordered Forest with all custom settings
#' orf_fit <- orf(X, Y, num.trees = 2000, mtry = 3, min.node.size = 10,
#'                      replace = TRUE, sample.fraction = 1,
#'                      honesty = FALSE, honesty.fraction = 0,
#'                      inference = FALSE, importance = FALSE)
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

    message("For conducting inference honesty is required. Honesty has been set to TRUE.")
    # set honesty to TRUE
    honesty <- TRUE
    # set honesty.fraction to 0.5
    honesty.fraction <- 0.5

  }

  if (replace == TRUE & inference == TRUE) {

    message("For conducting inference subsampling is required. Replace has been set to FALSE.")
    # set replace to FALSE
    replace <- FALSE
    # set sample.fraction to 0.5
    sample.fraction <- 0.5

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
  dat           <- as.data.frame(cbind(Y, X)) # dataframe
  colnames(dat) <- c(Y_name, X_name) # column names
  n             <- as.numeric(nrow(dat)) # number of observations

  # parameters (categories)
  categories <- as.numeric(sort(unique(Y))) # sequence of categories
  ncat       <- as.numeric(length(categories)) # number of categories
  cat        <- categories[1:(ncat-1)] # cat to esitmate / without the last category (not needed cuz P(Y_ind<=last_cat)=1)

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
    train_data  <- dat
    honest_data <- NULL

    ## create variables needed for orf estimations
    X_train     <- X
    # create indicator variables (outcomes)
    Y_ind_train <- lapply(cat, function(x) ifelse((Y <= x), 1, 0))

    # create dataset for ranger estimation
    data_ind    <- lapply(Y_ind_train, function(x) as.data.frame(cbind(as.matrix(unlist(x)), X_train)))

    # --------------------------------------------------------------------------------------- #

    # estimate ncat-1 forests (everything on the same data: placing splits and effect estimation)
    forest <- lapply(data_ind, function(x) ranger(dependent.variable.name = paste(Y_name), data = x,
                                                  num.trees = num.trees, mtry = mtry, min.node.size = min.node.size,
                                                  replace = replace, sample.fraction = sample.fraction,
                                                  importance = varimportance))

    # --------------------------------------------------------------------------------------- #

    # collect predictions for each forest based on whole sample (oob predictions)
    pred   <- lapply(forest, function(x) x$predictions) # collect forest predictions
    # add the probability for the last outcome (always 1)
    pred_1 <- append(pred, list(rep(1, n)))
    # prepend zero vector to predictions for later differencing
    pred_0 <- append(list(rep(0, n)), pred) # append a first 0 element for the list

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

    # no inference here
    var_final <- NULL

    # -------------------------------------------------------------------------------- #

  } else if (honesty == TRUE & inference == FALSE) {

    # --------------------------------------------------------------------------------------- #

    ## do honest forest estimation here using (preferably 50:50) data split as in Lechner (2019)
    # devide into 50:50 honesty sets
    split_data          <- honest_split(dat, honesty.fraction, orf = TRUE)
    # take care of train data
    train_data          <- split_data$trainData # take out training data
    rows_train_data     <- as.numeric(rownames(train_data)) # take rownames of train data as numeric
    Y_train             <- as.matrix(train_data[, 1]) # take out Y train
    colnames(Y_train)   <- Y_name # add column name
    X_train             <- train_data[, -1] # take out X
    colnames(X_train)   <- X_name # add column names
    # take care of honest data
    honest_data         <- split_data$honestData # take out honest data
    rows_honest_data    <- as.numeric(rownames(honest_data)) # take rownames of train data as numeric
    Y_honest            <- as.matrix(honest_data[, 1]) # take out Y train
    colnames(Y_honest)  <- Y_name # add column name
    X_honest            <- honest_data[, -1] # take out X
    colnames(X_honest)  <- X_name # add column names

    # --------------------------------------------------------------------------------------- #

    ## create variables needed for orf estimations
    # create indicator variables (outcomes)
    Y_ind_train <- lapply(cat, function(x) ifelse((Y_train <= x), 1, 0)) # train
    Y_ind_honest <- lapply(cat, function(x) ifelse((Y_honest <= x), 1, 0))

    # create dataset for ranger estimation
    data_ind_train <- lapply(Y_ind_train, function(x) as.data.frame(cbind(as.matrix(unlist(x)), X_train)))
    data_ind_honest <- lapply(Y_ind_honest, function(x) as.data.frame(cbind(as.matrix(unlist(x)), X_honest)))

    # --------------------------------------------------------------------------------------- #

    # estimate ncat-1 forests (on the train data: placing splits)
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

    # rename for output
    data_ind <- data_ind_honest
    # no inference here
    var_final <- NULL

    # -------------------------------------------------------------------------------- #

  } else if (honesty == TRUE & inference == TRUE) {

    # --------------------------------------------------------------------------------------- #

    ## do honest forest estimation here using (preferably 50:50) data split as in Lechner (2019)
    # devide into 50:50 honesty sets
    split_data          <- honest_split(dat, honesty.fraction, orf = TRUE)
    # take care of train data
    train_data          <- split_data$trainData # take out training data
    rows_train_data     <- as.numeric(rownames(train_data)) # take rownames of train data as numeric
    Y_train             <- as.matrix(train_data[, 1]) # take out Y train
    colnames(Y_train)   <- Y_name # add column name
    X_train             <- train_data[, -1] # take out X
    colnames(X_train)   <- X_name # add column names
    # take care of honest data
    honest_data         <- split_data$honestData # take out honest data
    rows_honest_data    <- as.numeric(rownames(honest_data)) # take rownames of train data as numeric
    Y_honest            <- as.matrix(honest_data[, 1]) # take out Y train
    colnames(Y_honest)  <- Y_name # add column name
    X_honest            <- honest_data[, -1] # take out X
    colnames(X_honest)  <- X_name # add column names

    # --------------------------------------------------------------------------------------- #

    ## create variables needed for orf estimations
    # create indicator variables (outcomes)
    Y_ind_train <- lapply(cat, function(x) ifelse((Y_train <= x), 1, 0)) # train
    Y_ind_honest <- lapply(cat, function(x) ifelse((Y_honest <= x), 1, 0))

    # create dataset for ranger estimation
    data_ind_train <- lapply(Y_ind_train, function(x) as.data.frame(cbind(as.matrix(unlist(x)), X_train)))
    data_ind_honest <- lapply(Y_ind_honest, function(x) as.data.frame(cbind(as.matrix(unlist(x)), X_honest)))

    # --------------------------------------------------------------------------------------- #

    # estimate ncat-1 forests (on the train data: placing splits
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

    # compute the variances for the categorical predictions
    var_final <- get_orf_variance(honest_pred, honest_weights, train_pred, train_weights, Y_ind_honest)

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
    var_imp_forests     <- matrix(unlist(lapply(forest, function(x) x$variable.importance)), nrow = ncol(var_imp_shares), byrow = T) # get the variable importances for each of binary forests

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
  forest_info        <- list(inputs, train_data, honest_data, categories, data_ind)
  names(forest_info) <- c("inputs", "trainData", "honestData", "categories", "indicatorData")

  # save forest accuracy measures
  forest_accuracy        <- list(pred_mse, pred_rps)
  names(forest_accuracy) <- c("MSE", "RPS")

  # define output of the function
  output        <- list(forest, forest_info, pred_final, var_final, var_imp, forest_accuracy)
  names(output) <- c("forests",  "info", "predictions", "variances", "importance", "accuracy")

  # --------------------------------------------------------------------------------------- #

  ## Set the name for the class
  class(output) <- "orf"

  # -------------------------------------------------------------------------------- #

  # return the output of the function
  return(output)

  # -------------------------------------------------------------------------------- #

}

#' Plot of the Ordered Forest
#'
#' plot the probability distributions estimated by the Ordered Forest object of class \code{orf}
#'
#' \code{plot.orf} generates probability distributions, i.e. density plots of estimated
#' ordered probabilities by the Ordered Forest for each outcome class considered.
#' The plots effectively visualize the estimated probability density in contrast to
#' a real observed ordered outcome class and as such provide a visual inspection of
#' the overall in-sample estimation accuracy. The dashed lines locate the means of
#' the respective probability distributions.
#'
#' @param x estimated Ordered Forest object of class \code{orf}
#' @param ... further arguments (currently ignored)
#'
#' @import ggplot2
#' @importFrom utils stack
#'
#' @author Gabriel Okasa
#'
#' @examples
#' # Ordered Forest
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
#'
#' # plot the estimated probability distributions
#' plot(orf_fit)
#'
#' @export
plot.orf <- function(x, ...) {

  # -------------------------------------------------------------------------------- #

  ## get forest as x
  forest <- x
  ## save forest inputs
  inputs      <- forest$info$inputs
  honesty     <- inputs$honesty
  categories  <- forest$info$categories
  honest_data <- forest$info$honestData
  train_data  <- forest$info$trainData

  # -------------------------------------------------------------------------------- #

  # get predictions and estimation data
  probabilities <- forest$predictions # take out honest predictions
  all_data      <- rbind(honest_data, train_data) # put data together
  all_data      <- all_data[order(as.numeric(row.names(all_data))), ] # sort data as original
  outcomes      <- all_data[, 1] # take the observed outcomes

  # -------------------------------------------------------------------------------- #

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
    # add column of outcome category to each list entry
    df_plot_prob[[i]]$Outcome <- paste("Class", i, sep = " ")
    # add colnames
    colnames(df_plot_prob[[i]]) <- c("Probability", "Density", "Outcome")
    # return the list
    return(df_plot_prob[[i]])  })

  # stack the dataframes under each other
  df_plot_prob <- as.data.frame(do.call(rbind, df_plot_prob))

  # add colnames and column indicating the outcome category
  df_plot_mean <- lapply(seq_along(df_plot_mean), function(i) {
    # add column of outcome category to each list entry
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


#' Summary of the Ordered Forest
#'
#' summary of an estimated Ordered Forest object of class \code{orf}
#'
#' \code{summary.orf} provides a short summary of the Ordered Forest estimation,
#' including the input information regarding the values of hyperparameters as
#' well as the output information regarding the prediction accuracy.
#'
#' @param object estimated Ordered Forest object of class \code{orf}
#' @param latex logical, if TRUE latex coded summary will be generated (default is FALSE)
#' @param ... further arguments (currently ignored)
#'
#' @importFrom xtable xtable
#'
#' @author Gabriel Okasa
#'
#' @examples
#' # Ordered Forest
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
#'
#' # show summary of the orf estimation
#' summary(orf_fit)
#'
#' # show summary of the orf estimation coded in LaTeX
#' summary(orf_fit, latex = TRUE)
#'
#' @export
summary.orf <- function(object, latex = FALSE, ...) {

  # -------------------------------------------------------------------------------- #

  ## check user inputs
  latex <- check_latex(latex)

  ## get forest as object
  forest <- object

  # -------------------------------------------------------------------------------- #

  ## save forest inputs
  main_class        <- class(forest)[1]
  inputs            <- forest$info$inputs

  honesty           <- inputs$honesty
  honesty.fraction  <- inputs$honesty.fraction
  inference         <- inputs$inference
  importance        <- inputs$importance
  mtry              <- inputs$mtry
  num.trees         <- inputs$num.trees
  min.node.size     <- inputs$min.node.size
  replace           <- inputs$replace
  sample.fraction   <- inputs$sample.fraction
  honest_data       <- forest$info$honestData
  train_data        <- forest$info$trainData
  categories        <- length(forest$info$categories)
  type              <- "Ordered Forest"

  # -------------------------------------------------------------------------------- #

  ## honest splitting, i.e. use honest data
  # take out summary statistics
  mse         <- round(forest$accuracy$MSE, 5)
  rps         <- round(forest$accuracy$RPS, 5)
  trainsize   <- nrow(train_data)
  honestsize  <- ifelse(is.null(honest_data), 0, nrow(honest_data))
  features    <- ncol(train_data) - 1   # take out the response

  # check if subsampling or bootstrapping was used
  if (forest$forests[[1]]$replace == TRUE) { build <- "Bootstrap" } else { build <- "Subsampling" }

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

  cat("Summary of the", type, "Estimation \n")

  # return output
  print(noquote(output), comment = FALSE)

  # -------------------------------------------------------------------------------- #

}


#' Print of the Ordered Forest
#'
#' print of an estimated Ordered Forest object of class \code{orf}
#'
#' \code{print.orf} provides a first glimpse of the Ordered Forest estimation,
#' printed directly to the \code{R} console. The printed information contains
#' the main inputs of the \code{orf} function.
#'
#' @param x estimated Ordered Forest object of class \code{orf}
#' @param ... further arguments (currently ignored)
#'
#' @author Gabriel Okasa
#'
#' @examples
#' # Ordered Forest
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
#'
#' # print output of the orf estimation
#' print(orf_fit)
#'
#' @export
print.orf <- function(x, ...) {

  # -------------------------------------------------------------------------------- #

  ## get forest as x
  forest <- x

  # -------------------------------------------------------------------------------- #

  ## save forest inputs
  main_class        <- class(forest)[1]
  inputs            <- forest$info$inputs

  honesty           <- inputs$honesty
  mtry              <- inputs$mtry
  num.trees         <- inputs$num.trees
  min.node.size     <- inputs$min.node.size
  replace           <- inputs$replace
  inference         <- inputs$inference
  honest_data       <- forest$info$honestData
  train_data        <- forest$info$trainData
  categories        <- length(forest$info$categories)
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


#' Prediction of the Ordered Forest
#'
#' Prediction for new observations based on estimated Ordered Forest of class \code{orf}
#'
#' \code{predict.orf} estimates the conditional ordered choice probabilities,
#' i.e. P[Y=m|X=x] for new data points (matrix X containing new observations
#' of covariates) based on the estimated Ordered Forest object of class \code{orf}.
#' Furthermore, weight-based inference for the probability predictions can be
#' conducted as well. If inference is desired, the supplied Ordered Forest must be
#' estimated with honesty and subsampling. If prediction only is desired, estimation
#' without honesty and with bootstrapping is recommended for optimal prediction
#' performance. Additionally to the probability predictions, class predictions can
#' be estimated as well using the \code{type} argument. In this case, the predicted
#' classes are obtained as classes with the highest predicted probability.
#'
#' @seealso \code{\link{summary.orf.prediction}}, \code{\link{print.orf.prediction}}
#'
#' @param object estimated Ordered Forest object of class \code{orf}
#' @param newdata numeric matrix X containing the observations for which the outcomes should be predicted
#' @param type string, specifying the type of the prediction, These can be either "probs" or  "p" for probabilities and "class" or "c" for classes. (Default is "probs").
#' @param inference logical, if TRUE variances for the predictions will be estimated (only feasible for probability predictions).
#' @param ... further arguments (currently ignored)
#'
#' @import ranger
#'
#' @return object of class \code{orf.prediction} with following elements
#'       \item{info}{info containing forest inputs and data used}
#'       \item{predictions}{predicted values}
#'       \item{variances}{variances of predicted values}
#'
#' @author Gabriel Okasa
#'
#' @examples
#' # Ordered Forest
#' require(orf)
#'
#' # load example data
#' data(odata)
#'
#' # specify response and covariates for train and test
#' idx <- sample(seq(1, nrow(odata), 1), 0.8*nrow(odata))
#'
#' # train set
#' Y_train <- as.numeric(odata[idx, 1])
#' X_train <- as.matrix(odata[idx, -1])
#'
#' # test set
#' Y_test <- as.numeric(odata[-idx, 1])
#' X_test <- as.matrix(odata[-idx, -1])
#'
#' # estimate Ordered Forest
#' orf_fit <- orf(X_train, Y_train)
#'
#' # predict the probabilities with the estimated orf
#' orf_pred <- predict(orf_fit, newdata = X_test)
#' \donttest{
#' # predict the probabilities with estimated orf together with variances
#' orf_pred <- predict(orf_fit, newdata = X_test, inference = TRUE)
#'
#' # predict the classes with estimated orf
#' orf_pred <- predict(orf_fit, newdata = X_test, type = "class")
#' }
#'
#' @export
predict.orf <- function(object, newdata = NULL, type = NULL, inference = NULL, ...) {

  # -------------------------------------------------------------------------------- #

  ## standard checks for input data
  if (class(object) != "orf") {
    stop("Forest object is not of class orf. Programme terminated.")
  }

  ## get forest as an object
  forest <- object
  ## save forest inputs
  inputs          <- forest$info$inputs
  categories      <- forest$info$categories
  replace         <- inputs$replace
  honesty         <- inputs$honesty
  honest_data     <- forest$info$honestData
  train_data      <- forest$info$trainData
  honest_ind_data <- forest$info$indicatorData # indicator data needed for indicator predictions

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
    pred_final <- forest$predictions
    var_final  <- forest$variances

    # check the desired type of predictions
    if (type == "class" | type == "c") {

      # convert probabilities into class predictions ("ordered classification")
      pred_final <- as.matrix(apply(pred_final, 1, which.max))
      var_final  <- NULL

    }

  } else if (is.null(newdata) & (inference == FALSE) & (inputs$inference == TRUE)) {

    # then take the estimated values but dont supply the inference results
    pred_final <- forest$predictions
    var_final  <- NULL

    # check the desired type of predictions
    if (type == "class" | type == "c") {

      # convert probabilities into class predictions ("ordered classification")
      pred_final <- as.matrix(apply(pred_final, 1, which.max))

    }

  } else {

    # take out list of ranger objects
    forest <- forest$forests

    ## get train data names (only X)
    train_data_name <- colnames(train_data)[2:ncol(train_data)]

    if (is.null(newdata) & (inference == TRUE) & (inputs$inference == FALSE)) {

      # note that the newdata is the estimation data and inference should reflect that
      flag_newdata  <- 1
      # take traindata as newdata and estimate the weights needed for inference
      newdata       <- as.data.frame(rbind(train_data, honest_data))
      # sort according to rownames
      newdata       <- as.data.frame(newdata[order(as.numeric(row.names(newdata))), ])
      # get further values
      n_newdata     <- nrow(newdata) # rows of new data
      n_cat         <- as.numeric(length(categories))

    } else if (!is.null(newdata)) {

      ## get X matrix newdata as dataframe and check colnames
      # X
      if (is.null(colnames(newdata))) { colnames(newdata) <- paste0("X", rep(1:ncol(newdata))) } # check if X has name
      newdata_name  <- colnames(newdata) # save the name of X
      newdata       <- as.data.frame(newdata) # as dataframe
      n_newdata     <- nrow(newdata) # rows of new data
      n_cat         <- as.numeric(length(categories))

      # assign flag_newdata to zero
      flag_newdata <- 0

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

        ## use standard pred_orf_variance
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
  names(output) <- c("info", "predictions", "variances")

  # -------------------------------------------------------------------------------- #

  ## Set the name for the class
  class(output) <- "orf.prediction"

  # return output
  return(output)

  # -------------------------------------------------------------------------------- #

}


#' Print of the Ordered Forest Prediction
#'
#' print of Ordered Forest predictions of class \code{orf.prediction}
#'
#' \code{print.orf.prediction} provides a first glimpse of the Ordered Forest
#' prediction, printed directly to the \code{R} console. The printed information
#' contains the main inputs of the \code{predict.orf} function.
#'
#' @param x predicted Ordered Forest object of class \code{orf.prediction}
#' @param ... further arguments (currently ignored)
#'
#' @author Gabriel Okasa
#'
#' @examples
#' # Ordered Forest
#' require(orf)
#'
#' # load example data
#' data(odata)
#'
#' # specify response and covariates for train and test
#' idx <- sample(seq(1, nrow(odata), 1), 0.8*nrow(odata))
#'
#' # train set
#' Y_train <- as.numeric(odata[idx, 1])
#' X_train <- as.matrix(odata[idx, -1])
#'
#' # test set
#' Y_test <- as.numeric(odata[-idx, 1])
#' X_test <- as.matrix(odata[-idx, -1])
#'
#' # estimate Ordered Forest
#' orf_fit <- orf(X_train, Y_train)
#'
#' # predict the probabilities with the estimated orf
#' orf_pred <- predict(orf_fit, newdata = X_test)
#'
#' # print the prediction object
#' print(orf_pred)
#'
#' @export
print.orf.prediction <- function(x, ...) {

  # -------------------------------------------------------------------------------- #

  ## get forest predictions as x
  forest_pred <- x

  # -------------------------------------------------------------------------------- #

  ## save forest prediction inputs
  main_class        <- class(forest_pred)[1]
  inputs            <- forest_pred$info$inputs

  honesty           <- inputs$honesty
  mtry              <- inputs$mtry
  num.trees         <- inputs$num.trees
  min.node.size     <- inputs$min.node.size
  replace           <- inputs$replace
  inference         <- inputs$inference

  pred_data         <- forest_pred$info$newData
  pred_type         <- forest_pred$info$predType
  pred_inference    <- forest_pred$info$predInference
  categories        <- length(forest_pred$info$categories)
  build             <- ifelse(replace == TRUE, "Bootstrap", "Subsampling")
  type              <- "Ordered Forest Prediction"

  # define nice output for pred_type
  if (pred_type == "p" | pred_type == "probs") {

    # probabilities
    pred_type <- "probability"

  } else if (pred_type == "c" | pred_type == "class") {

    # classes
    pred_type <- "class"

  }

  # -------------------------------------------------------------------------------- #

  cat(type, "object of class", main_class, "\n\n")

  cat("Prediction Type:                 ", pred_type, "\n")
  cat("Number of Categories:            ", categories, "\n")
  cat("Sample Size:                     ", nrow(forest_pred$predictions), "\n")
  cat("Number of Trees:                 ", num.trees, "\n")
  cat("Build:                           ", build, "\n")
  cat("Mtry:                            ", mtry, "\n")
  cat("Minimum Node Size:               ", min.node.size, "\n")
  cat("Honest Forest:                   ", honesty, "\n")
  cat("Weight-Based Inference:          ", pred_inference, "\n")

  # -------------------------------------------------------------------------------- #

}


#' Summary of the Ordered Forest Prediction
#'
#' summary of Ordered Forest predictions of class \code{orf.prediction}
#'
#' \code{summary.orf.prediction} provides a main summary of the Ordered Forest
#' prediction, including the input information regarding the values of hyperparameters
#' as well as the inputs of the \code{predict.orf} function.
#'
#' @param object predicted Ordered Forest object of class \code{orf.prediction}
#' @param latex logical, if TRUE latex coded summary will be generated (default is FALSE)
#' @param ... further arguments (currently ignored)
#'
#' @author Gabriel Okasa
#'
#' @examples
#' # Ordered Forest
#' require(orf)
#'
#' # load example data
#' data(odata)
#'
#' # specify response and covariates for train and test
#' idx <- sample(seq(1, nrow(odata), 1), 0.8*nrow(odata))
#'
#' # train set
#' Y_train <- as.numeric(odata[idx, 1])
#' X_train <- as.matrix(odata[idx, -1])
#'
#' # test set
#' Y_test <- as.numeric(odata[-idx, 1])
#' X_test <- as.matrix(odata[-idx, -1])
#'
#' # estimate Ordered Forest
#' orf_fit <- orf(X_train, Y_train)
#'
#' # predict the probabilities with the estimated orf
#' orf_pred <- predict(orf_fit, newdata = X_test)
#'
#' # summary of the prediction object
#' summary(orf_pred)
#'
#' # show summary of the orf prediction coded in LaTeX
#' summary(orf_pred, latex = TRUE)
#'
#' @export
summary.orf.prediction <- function(object, latex = FALSE, ...) {

  # -------------------------------------------------------------------------------- #

  ## check user inputs
  latex <- check_latex(latex)

  ## get forest predictions as object
  forest_pred <- object

  # -------------------------------------------------------------------------------- #

  ## save forest prediction inputs
  main_class        <- class(forest_pred)[1]
  inputs            <- forest_pred$info$inputs

  honesty           <- inputs$honesty
  honesty.fraction  <- inputs$honesty.fraction
  mtry              <- inputs$mtry
  num.trees         <- inputs$num.trees
  min.node.size     <- inputs$min.node.size
  replace           <- inputs$replace
  sample.fraction   <- inputs$sample.fraction
  inference         <- inputs$inference

  pred_data         <- forest_pred$info$newData
  pred_type         <- forest_pred$info$predType
  pred_inference    <- forest_pred$info$predInference
  categories        <- length(forest_pred$info$categories)
  sample_size       <- nrow(forest_pred$predictions)
  build             <- ifelse(replace == TRUE, "Bootstrap", "Subsampling")
  type              <- "Ordered Forest Prediction"

  # define nice output for pred_type
  if (pred_type == "p" | pred_type == "probs") {

    # probabilities
    pred_type <- "probability"

  } else if (pred_type == "c" | pred_type == "class") {

    # classes
    pred_type <- "class"

  }

  # -------------------------------------------------------------------------------- #

  # structure summary into a list
  output        <- list(type, pred_type, categories, build, num.trees, mtry, min.node.size, replace, sample.fraction, honesty, honesty.fraction, pred_inference, sample_size)
  names(output) <- c("type", "prediction.type", "categories", "build", "num.trees", "mtry", "min.node.size", "replace", "sample.fraction", "honesty", "honesty.fraction", "inference", "sample.size")

  # output matrix
  output_matrix <- matrix(NA, ncol = 1, nrow = length(output))
  # populate output matrix
  rownames(output_matrix) <- names(output) # rownames are names
  colnames(output_matrix) <- "" # no visible colname
  output_matrix[, 1]      <- unlist(output) # column 1 are values

  # generate latex output if selected
  if (latex == TRUE) { colnames(output_matrix) <- "Ordered Forest Prediction Summary"
  output_matrix <- xtable(output_matrix, caption = "Summary of the Ordered Forest Prediction", align = "ll")
  }

  # pack it into output
  output <- output_matrix

  # -------------------------------------------------------------------------------- #

  cat("Summary of the", type, "\n")

  # return output
  print(noquote(output), comment = FALSE)

  # -------------------------------------------------------------------------------- #

}

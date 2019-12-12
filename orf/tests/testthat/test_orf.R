library(orf)
context("orf")
set.seed(123)

# load example data
data(odata)

# set X and Y for the estimation
Y <- as.numeric(odata[1:100, 1])
X <- as.matrix(odata[1:100, -1])

# testing the main orf function
test_that("orf object is of class orf", {
  orf <- orf(X, Y)
  expect_s3_class(orf, "orf")
})

test_that("orf object is of type list", {
  orf <- orf(X, Y)
  expect_type(orf, "list")
})

# forests
test_that("orf estimates ncat - 1 forests", {
  orf <- orf(X, Y)
  expect_equal(length(orf$forests), length(unique(Y))-1)
})

# probabilities
test_that("probabilities sum up to 1", {
  orf <- orf(X, Y)
  expect_equal(rowSums(orf$predictions), rep(1, nrow(X)))
})

test_that("probabilities are non-negative", {
  orf <- orf(X, Y)
  expect_true(all(as.numeric(orf$predictions) >= 0))
})

test_that("probabilities reflect all classes", {
  orf <- orf(X, Y)
  expect_equal(ncol(orf$predictions), length(unique(Y)))
})

test_that("probabilities reflect all observations", {
  orf <- orf(X, Y)
  expect_equal(nrow(orf$predictions), length(Y))
})

test_that("fitted values are equal to predicted values in training set", {
  orf <- orf(X, Y)
  orf_fitted <- orf$predictions
  orf_predicted <- predict(orf)$predictions
  expect_equal(orf_fitted, orf_predicted)
})

# variance
test_that("variances of the predictions are positive", {
  orf <- orf(X, Y, inference = TRUE)
  expect_true(all(as.numeric(orf$variances) > 0))
})

test_that("variances of the in sample predictions are the same for predict", {
  orf <- orf(X, Y, inference = TRUE)
  orf_vars <- orf$variances
  orf_preds <- predict(orf, inference = TRUE)$variances
  expect_equal(orf_vars, orf_preds)
})

test_that("variances without honesty and subsampling throw a warning", {
  expect_message(orf(X, Y, honesty = FALSE, replace = TRUE, inference = TRUE))
  expect_message(orf(X, Y, honesty = FALSE, replace = FALSE, inference = TRUE))
  expect_message(orf(X, Y, honesty = TRUE, replace = TRUE, inference = TRUE))
})

test_that("variances are possible with different honesty fractions", {
  orf <- orf(X, Y, inference = FALSE)
  expect_null(orf$variances)
  orf <- orf(X, Y, inference = TRUE)
  expect_vector(orf$variances)
  orf <- orf(X, Y, inference = TRUE, honesty.fraction = 0.6)
  expect_vector(orf$variances)
  orf <- orf(X, Y, inference = TRUE, honesty.fraction = 0.4)
  expect_vector(orf$variances)
})

# prediction error
test_that("prediction errors are positive", {
  orf <- orf(X, Y)
  expect_true(orf$accuracy$MSE > 0)
  expect_true(orf$accuracy$RPS > 0)
})

# variable importance
test_that("variable importance is possible in all settings", {
  orf <- orf(X, Y, importance = FALSE)
  expect_null(orf$importance)
  orf <- orf(X, Y, importance = TRUE, honesty = FALSE, inference = FALSE)
  expect_vector(orf$importance)
  orf <- orf(X, Y, importance = TRUE, honesty = TRUE, inference = FALSE)
  expect_vector(orf$importance)
  orf <- orf(X, Y, importance = TRUE, honesty = TRUE, inference = TRUE)
  expect_vector(orf$importance)
})

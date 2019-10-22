library(orf)
context("orf")

# load example data
data(odata)

# set X and Y for the estimation
Y <- odata[1:100, 1]
X <- odata[1:100, -1]

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
  expect_equal(length(orf$trainForests), length(unique(Y))-1)
})

# probabilities
test_that("probabilities sum up to 1", {
  orf <- orf(X, Y)
  expect_equal(rowSums(orf$forestPredictions), rep(1, 100))
})

test_that("probabilities are non-negative", {
  orf <- orf(X, Y)
  expect_true(all(as.numeric(orf$forestPredictions) >= 0))
})

test_that("probabilities reflect all classes", {
  orf <- orf(X, Y)
  expect_equal(ncol(orf$forestPredictions), length(unique(Y)))
})

test_that("probabilities reflect all observations", {
  orf <- orf(X, Y)
  expect_equal(nrow(orf$forestPredictions), length(Y))
})

test_that("fitted values are equal to predicted values in training set", {
  orf <- orf(X, Y)
  orf_fitted <- orf$forestPredictions
  orf_predicted <- predict(orf)$forestPredictions
  expect_equal(orf_fitted, orf_predicted)
})

# classes
test_that("predicted classes reflect actual classes", {
  orf <- orf(X, Y)
  orf_classes <- predict(orf, type = "class")$forestPredictions
  expect_equal(sort(unique(as.numeric(orf_classes))), sort(unique(Y)))
})

# variance
test_that("variances of the predictions are positive", {
  orf <- orf(X, Y, inference = TRUE)
  expect_true(all(as.numeric(orf$forestVariances) > 0))
})

test_that("variances of the in sample predictions are the same for predict", {
  orf <- orf(X, Y, inference = TRUE)
  orf_vars <- orf$forestVariances
  orf_preds <- predict(orf, inference = TRUE)$forestVariances
  expect_equal(orf_vars, orf_preds)
})

test_that("variances without honesty and subsampling throw a warning", {
  expect_message(orf(X, Y, honesty = FALSE, replace = TRUE, inference = TRUE))
  expect_message(orf(X, Y, honesty = FALSE, replace = FALSE, inference = TRUE))
  expect_message(orf(X, Y, honesty = TRUE, replace = TRUE, inference = TRUE))
})

# prediction error
test_that("prediction errors are positive", {
  orf <- orf(X, Y)
  expect_true(orf$MSE > 0)
  expect_true(orf$RPS > 0)
})
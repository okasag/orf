library(orf)
context("margins")
set.seed(123)

# load example data
data(odata)

# set X and Y for the estimation
Y <- as.numeric(odata[1:100, 1])
X <- as.matrix(odata[1:100, -1])

# testing the main margins function
test_that("margins.orf object is of class margins.orf", {
  orf <- orf(X, Y)
  orf_margins <- margins(orf)
  expect_s3_class(orf_margins, "margins.orf")
})

test_that("margins.orf object is of type list", {
  orf <- orf(X, Y)
  orf_margins <- margins(orf)
  expect_type(orf_margins, "list")
})

# marginal effects
test_that("marginal effects sum up to 0", {
  orf <- orf(X, Y)
  orf_margins <- margins(orf)
  expect_equal(as.numeric(rowSums(orf_margins$effects)), rep(0, ncol(X)))
})

test_that("marginal effects variances etc are positive", {
  orf <- orf(X, Y)
  orf_margins <- margins(orf, inference = TRUE)
  expect_true(all(as.numeric(orf_margins$variances) > 0))
  expect_true(all(as.numeric(orf_margins$errors) > 0))
  expect_true(all(as.numeric(orf_margins$pvalues) >= 0) & all(as.numeric(orf_margins$pvalues) <= 1))
})

test_that("not admissible marginal effects evaluation point throw a warning", {
  orf <- orf(X, Y)
  expect_warning(margins(orf, eval = "something"))
})

test_that("inference cannot be done if orf was not estimated properly", {
  orf <- orf(X, Y, honesty = FALSE, replace = TRUE)
  expect_warning(margins(orf, inference = TRUE))
})

test_that("window must be within 0 and 1", {
  orf <- orf(X, Y)
  expect_error(margins(orf, window = 2))
})

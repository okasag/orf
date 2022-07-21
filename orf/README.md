[![CRAN checks](https://cranchecks.info/badges/summary/orf)](https://cran.r-project.org/web/checks/check_results_orf.html)
[![](https://www.r-pkg.org/badges/version/orf?color=brightgreen)](https://cran.r-project.org/package=orf)
[![](http://cranlogs.r-pkg.org/badges/grand-total/orf?color=blue)](https://cran.r-project.org/package=orf)

# orf: ordered random forests

## Introduction

The `R` package `orf` is an implementation of the Ordered Forest estimator
as developed in Lechner & Okasa (2019). The Ordered Forest flexibly estimates
the conditional probabilities of models with ordered categorical outcomes
(so-called ordered choice models). Additionally to common machine learning
algorithms the `orf` package provides functions for estimating marginal effects
as well as statistical inference thereof and thus provides similar output as in
standard econometric models for ordered choice. The core forest algorithm relies
on the fast `C++` forest implementation from the `ranger` package (Wright & Ziegler, 2017).

## Installation

In order to install the latest `CRAN` released version use:

```r
install.packages("orf", dependencies = c("Imports", "Suggests"))
```

to make sure all the needed packages are installed as well. Note that if you install
the package directly from the source a `C++` compiler is required. For Windows 
users `Rtools` collection is required too.

## Examples

The examples below demonstrate the basic functionality of the `orf` package.

```r
## Ordered Forest
require(orf)

# load example data
data(odata)

# specify response and covariates
Y <- as.numeric(odata[, 1])
X <- as.matrix(odata[, -1])

# estimate Ordered Forest with default settings
orf_fit <- orf(X, Y, num.trees = 1000, mtry = 2, min.node.size = 5,
                     replace = FALSE, sample.fraction = 0.5,
                     honesty = TRUE, honesty.fraction = 0.5,
                     inference = FALSE, importance = FALSE)

# print output of the Ordered Forest estimation
print(orf_fit)

# show summary of the Ordered Forest estimation
summary(orf_fit, latex = FALSE)

# plot the estimated probability distributions
plot(orf_fit)

# predict with the estimated Ordered Forest
predict(orf_fit, newdata = NULL, type = "probs", inference = FALSE)

# estimate marginal effects of the Ordered Forest
margins(orf_fit, newdata = NULL, eval = "mean", window = 0.1, inference = FALSE)
```

For a more detailed examples see the package vignette.

## References

- Lechner, M., & Okasa, G. (2019). Random Forest Estimation of the Ordered Choice Model. arXiv preprint arXiv:1907.02436. <https://arxiv.org/abs/1907.02436>
- Wright, M. N. & Ziegler, A. (2017). ranger: A fast implementation of random forests for high dimensional data in C++ and R. J Stat Softw 77:1-17. <https://doi.org/10.18637/jss.v077.i01>


# orf: ordered random forests

## Introduction

The `R` package `orf` is an implementation of the Ordered Forest estimator
as in Lechner & Okasa (2019). The Ordered Forest flexibly estimates the conditional
probabilities of models with ordered categorical outcomes (so-called ordered
choice models). Additionally to common machine learning algorithms the `orf`
package provides functions for estimating marginal effects as well as
statistical inference thereof and thus provides similar output as in standard
econometric models for ordered choice. The core forest algorithm relies on the
fast `C++` forest implementation from the `ranger` package (Wright & Ziegler, 2017).

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
orf <- orf(X, Y)

# print output of the orf estimation
print(orf)

# show summary of the orf estimation
summary(orf)

# plot the estimated probability distributions
plot(orf)

# predict with the estimated orf
predict(orf)

# estimate marginal effects of the orf
margins(orf)

## end of example
```

For a more detailed examples see the package vignette.

## References

- Lechner, M., & Okasa, G. (2019). Random Forest Estimation of the Ordered Choice Model. arXiv preprint arXiv:1907.02436. <https://arxiv.org/abs/1907.02436>
- Goller, D., Knaus, M. C., Lechner, M., & Okasa, G. (2018). Predicting Match Outcomes in Football by an Ordered Forest Estimator (No. 1811). University of St. Gallen, School of Economics and Political Science. <http://ux-tauri.unisg.ch/RePEc/usg/econwp/EWP-1811.pdf>
- Wright, M. N. & Ziegler, A. (2017). ranger: A fast implementation of random forests for high dimensional data in C++ and R. J Stat Softw 77:1-17. <https://doi.org/10.18637/jss.v077.i01>


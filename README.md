# orf: ordered random forests

An implementation of the Ordered Forest estimator
as in Lechner & Okasa (2019). The Ordered Forest flexibly
estimates the conditional probabilities of models with ordered
categorical outcomes (so-called ordered choice models).
Additionally to common machine learning algorithms the orf
package provides functions for estimating marginal effects as well
as statistical inference thereof and thus provides similar output
as in standard econometric models for ordered choice. The core
forest algorithm relies on the fast C++ forest implementation
from the ranger package (Wright & Ziegler, 2017).

## author

- Gabriel Okasa (gabriel.okasa@unisg.ch)

## dependencies

- main orf functions rely on the fast C++ implementation of the 
random forest algorithm from the ranger package (Wright & Ziegler, 2017)

## to do:

- test checks() on all platforms
- submit to CRAN


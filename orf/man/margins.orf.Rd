% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/margins.R
\name{margins.orf}
\alias{margins.orf}
\title{Marginal Effects for the Ordered Forest}
\usage{
\method{margins}{orf}(forest, eval = NULL, inference = NULL, window = NULL, newdata = NULL)
}
\arguments{
\item{forest}{estimated Ordered Forest object of class \code{orf}}

\item{eval}{string, defining evaluation point for marginal effects. These can be one of "mean", "atmean", or "atmedian". (Default is "mean")}

\item{inference}{logical, if TRUE inference on marginal effects will be conducted (default is inherited from the \code{orf} object)}

\item{window}{numeric, share of standard deviation of X to be used for evaluation of the marginal effect (default is 0.1)}

\item{newdata}{numeric matrix X containing the new observations for which the marginal effects should be estimated}
}
\value{
object of type \code{margins.orf} with following elements
      \item{info}{info containing forest inputs and data used}
      \item{effects}{marginal effects}
      \item{variances}{variances of marginal effects}
      \item{errors}{standard errors of marginal effects}
      \item{tvalues}{t-values of marginal effects}
      \item{pvalues}{p-values of marginal effects}
}
\description{
S3 method for estimation of marginal effects
of an Ordered Forest objects of class \code{orf}.
}
\details{
\code{margins.orf} estimates marginal effects at the mean, at the median, or
the mean marginal effects, depending on the \code{eval} argument. It is advised
to increase the number of subsampling replications in the supplied \code{orf}
object as the estimation of the marginal effects is a more demanding exercise
than a simple Ordered Forest estimation/prediction. Additionally to the estimation
of the marginal effects, the weight-based inference for the effects is supported
as well. Note, that the inference procedure is much more computationally exhausting
exercise due to the computation of the forest weights. Additionally, the evaluation
window for the marginal effects can be regulated through the \code{window} argument.
Furthermore, new data for which marginal effects should be computed can be supplied
as well as long as it lies within the support of \code{X}.
}
\examples{
## Ordered Forest
require(orf)

# load example data
data(odata)

# specify response and covariates
Y <- as.numeric(odata[, 1])
X <- as.matrix(odata[, -1])

# estimate Ordered Forest
orf_fit <- orf(X, Y)
\donttest{
# estimate marginal effects of the orf (default)
orf_margins <- margins(orf_fit)

# estimate marginal effects evaluated at the mean
orf_margins <- margins(orf_fit, eval = "atmean")

# estimate marginal effects with inference
# (orf object has to be estimated with honesty and subsampling)
orf_margins <- margins(orf_fit, inference = TRUE)

# estimate marginal effects with custom window size
orf_margins <- margins(orf_fit, window = 0.5)

# estimate marginal effects for some new data (within support of X)
orf_margins <- margins(orf_fit, newdata = X[1:10, ])

# estimate marginal effects with all custom settings
orf_margins <- margins(orf_fit, eval = "atmedian", inference = TRUE,
                                window = 0.5, newdata = X[1:10, ])
}

}
\seealso{
\code{\link{summary.margins.orf}}, \code{\link{print.margins.orf}}
}
\author{
Gabriel Okasa
}

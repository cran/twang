\name{twang-package}
\alias{twang-package}
\alias{twang}
\docType{package}
\title{
Toolkit for Weighting and Analysis of Nonequivalent Groups
}
\description{
This package offers functions for propensity score
estimating and weighting, nonresponse weighting, and diagnosis
of the weights
}
\details{
\tabular{ll}{
Package: \tab twang\cr
Version: \tab 0.6-3\cr
Date: \tab 2006-2-15\cr
Depends: \tab R (>= 2.2), gbm (>= 1.5-3), survey\cr
License: \tab GPL (version 2 or newer)\cr
Built: \tab R 2.2.1; ; 2006-02-15 10:35:20; windows\cr
}

Index:
\preformatted{
bal.stat                Calculate weighted balance statistics
bal.table               Compute balance table
check.err               Reports on errors and warnings
desc.wts                Diagnosis of weights
diag.plot               Create diagnostic plots
dx.wts                  Propensity score diagnostics
ks.stat                 Functions for evaluating balance
lalonde                 Lalonde's National Supported Work Demonstration
                        data
metric.i                Losses for *.stat functions
print.dxwts             Print a diagnosis of the weights
ps                      Propensity score estimation
ps.summary              Computes balance measures
sensitivity             Sensitivity analysis
stop.methods            Rules for selecting the propensity scores
summary.ps              Summarize a ps object
}

Further information is available in the following vignettes:
\tabular{ll}{
\code{twang} \tab Toolkit for Weighting and Analysis of Nonequivalent Groups: A guide to the twang package (source, pdf)\cr
}
}
\author{
Greg Ridgeway \email{gregr@rand.org}, Dan McCaffrey
\email{danielm@rand.org}, Andrew Morral \email{morral@rand.org}

Maintainer: Greg Ridgeway <gregr@rand.org>
}
\references{
Dan McCaffrey, G. Ridgeway, Andrew Morral (2004). \dQuote{Propensity Score Estimation
with Boosted Regression for Evaluating Adolescent Substance Abuse Treatment,}
\emph{Psychological Methods} 9(4):403-425.

G. Ridgeway (2006). \dQuote{Assessing the effect of race bias in 
post-traffic stop outcomes using propensity scores,} \emph{Journal of 
Quantitative Criminology} 22(1).
}
\keyword{ package }
\seealso{
The \code{\link[gbm:gbm-package]{gbm}} package
}

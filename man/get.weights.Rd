\name{get.weights}
\alias{get.weights}

\title{ Extract propensity score weights }
\description{
  Extracts propensity score weights from a ps object.
}
\usage{
get.weights(ps1,
            type = c("ATT", "ATE")[1],
            stop.method = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ps1}{ a \code{ps} object }
  \item{type}{ indicates whether the weights are for the average treatment effect on the treated (ATT) or the average treatment effect on the population (ATE) }
  \item{stop.method}{ indicates which set of weights to retrieve from the \code{ps} object }
}
\details{
  Weights for ATT are 1 for the treatment cases and p/(1-p) for the control cases.

  Weights for ATE are 1/p for the treatment cases and 1/(1-p) for the control cases.
}

\value{
a vector of weights
}

\seealso{ \code{\link{ps}} }

\keyword{ utilities }

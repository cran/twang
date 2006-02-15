\name{summary.ps}
\alias{summary.ps}

\title{ Summarize a ps object}
\description{
Computes summary information about a stored \code{ps} object
}
\usage{
\method{summary}{ps}(object, ...)
}

\arguments{
  \item{object}{ a \code{\link[twang]{ps}} object }
  \item{\dots}{ additional arguments affecting the summary produced }
}
\details{
Compresses the information in the \code{desc} component of the \code{ps} object
into a short summary table describing the size of the dataset and the quality of
the propensity score weights.
}
\value{
See \code{\link[twang]{ps}} for details on the returned table
}

\seealso{\code{\link[twang]{ps}}}

\keyword{ models }

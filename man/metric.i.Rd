\name{metric.i}
\alias{metric.i}

\title{Losses for *.stat functions}
\description{
Rearranges the arguments of the *.stat functions so that they may be passed to
\code{\link[stats]{optimize}}
}
\usage{
metric.i(i, fun = ks.stat, ...)
}
\arguments{
  \item{i}{ the number of \code{\link[gbm]{gbm}} iterations}
  \item{fun}{ a valid *.stat function }
  \item{\dots}{ other arguments to be passed to \code{fun} }
}

\value{
Evaluates \code{fun} at \code{i}
}

\seealso{ \code{\link{ks.stat}}, \code{\link{es.stat}}, 
          \code{\link{strata.stat}} }

\keyword{ models }% at least one, from doc/KEYWORDS

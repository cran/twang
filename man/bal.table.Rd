\name{bal.table}
\alias{bal.table}
\alias{bal.table.ps}
\alias{bal.table.dxwts}

\title{Compute balance table}
\description{
 Extract the balance table from \code{\link[twang]{ps}} and 
 \code{\link[twang]{dxwts}} objects
}
\usage{
bal.table(x)
}

\arguments{
  \item{x}{a \code{\link[twang]{ps}} or \code{\link[twang]{dxwts}} object}
}
\details{
\code{bal.table} is a generic function that will dispatch the right method for
extracting the balance tables from various \code{twang} objects, currently
\code{\link[twang]{ps}} and \code{\link[twang]{dxwts}} objects
}
\value{
See \code{\link[twang]{ps.summary}} for details on the returned object
}

\seealso{ The example for \code{\link[twang]{ps}} contains an example of the 
use of \code{bal.table}}

\keyword{ models }

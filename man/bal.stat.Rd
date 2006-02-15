\name{bal.stat}
\alias{bal.stat}

\title{ Calculate weighted balance statistics }
\description{
\code{bal.stat} compares the treatment and control subjects by means, standard 
deviations, effect size, and KS statistics
}
\usage{
bal.stat(data, 
         vars = NULL, 
         treat.var, 
         w.all, 
         get.means = TRUE,
         get.ks = TRUE, 
         na.action = "level")
}

\arguments{
  \item{data}{ a data frame containing the data }
  \item{vars}{ a vector of character strings with the names of the variables 
               on which the function will assess the balance }
  \item{treat.var}{ the name of the treatment variable }
  \item{w.all}{ observation weights (e.g. propensity score weights, sampling 
                weights, or both) }
  \item{get.means}{ logical. If \code{TRUE} then \code{bal.stat} will compute means
                 and variances}
  \item{get.ks}{ logical. If \code{TRUE} then \code{bal.stat} will compute KS 
                 statistics}
  \item{na.action}{ a character string indicating how \code{bal.stat} should 
                    handle missing values. Current options are "level", 
                    "exclude", or "lowest" }
}
\details{
    \code{bal.stat} calls \code{\link{ps.summary.f}} and 
    \code{\link{ps.summary.n}} for each variable and assembles the results in 
    a table
}
\value{
See \code{\link[twang]{ps.summary}} for details on the returned object. 
\code{get.means} and \code{get.ks} manipulate the inclusion of certain 
columns in the returned result.
}

\seealso{ The example for \code{\link[twang]{ps}} contains an example of the 
use of \code{bal.table} }

\keyword{ multivariate }

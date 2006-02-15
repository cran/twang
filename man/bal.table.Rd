\name{bal.table}
\alias{bal.table}
\alias{bal.table.ps}
\alias{bal.table.dxwts}

\title{Compute balance table}
\description{
 Extract the balance table from \code{\link{ps}} and 
 \code{\link{dx.wts}} objects
}
\usage{
bal.table(x)
}

\arguments{
  \item{x}{a \code{\link{ps}} or \code{\link{dx.wts}} object}
}

\details{ \code{bal.table} is a generic function for extracting balance 
tables from \code{\link{ps}} and \code{\link{dx.wts}} objects. These objects 
usually have several sets of candidate weights, one for an unweighted 
analysis and perhaps several \code{\link{stop.methods}}. \code{bal.table} 
will return a table for each set of weights combined into a list. Each list 
component will be named as given in the \code{x}, usually the name of the 
\code{stop.method}. The balance table labeled \dQuote{unw} indicates the 
unweighted analysis.}

\value{
    Returns a data frame containing the balance information.            
  \item{tx.mn}{The mean of the treatment group}
  \item{tx.sd}{The standard deviation of the treatment group}
  \item{ct.mn}{The mean of the control group}
  \item{ct.sd}{The standard deviation of the control group}
  \item{std.eff.sz}{The standardized effect size, (tx.mn-ct.mn)/tx.sd}
  \item{stat}{the t-statistic for numeric variables and the chi-square 
              statistic for continuous variables}
  \item{p}{the p-value for the test associated with \code{stat}}
  \item{ks}{the KS statistic}
  \item{ks.pval}{the KS p-value computed using the analytic approximation,
                 which does not necessarily work well with a lot of ties}  
}

\seealso{ The example for \code{\link{ps}} contains an example of the 
use of \code{bal.table}}

\keyword{ models }

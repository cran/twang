\name{check.err}
\alias{check.err}

\title{ Reports on errors and warnings }
\description{
Reports on errors and warnings encountered while running \code{\link{ps}} 
and \code{\link{desc.wts}}. This function is not intended for the user to
call directly
}
\usage{
check.err(cov.table, stage, alerts.stack)
}

\arguments{
  \item{cov.table}{a balance table, intended to be the \code{results} 
                   component of the list that \code{\link{bal.stat}} 
                   returns}
  \item{stage}{a title for the method ``type" used to create the weights, used
               to label the results}
  \item{alerts.stack}{an object for collecting warnings issued during the 
                      analyses}
}
\details{
Checks for treatment standard deviations that are exceedingly small or zero and 
for effect sizes that are unusually large, both indicative of numerical problems
or extreme sample imbalance
}
\value{
\code{check.err} returns no objects but does alter the \code{alerts.stack} 
object
}

\keyword{error}

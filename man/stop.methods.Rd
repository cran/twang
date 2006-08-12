\name{stop.methods}
\alias{stop.methods}

\title{ Rules for selecting the propensity scores }

\description{
  A list of \code{stop.method} objects built into the \code{twang} package 
  that encode rules for selecting propensity score weights
}

\details{
The \code{\link{ps}} function uses a \code{stop.method} object for instructions 
    on how to select the propensity score weights. \code{twang} has some 
    \code{stop.method} objects built in but the user may implement their own if 
    they wish and pass them to \code{\link{ps}} for it to optimize.

A valid \code{stop.method} object is a list that defines the following
\describe{
\item{metric}{a function that evaluates the similarity of the distribution of 
    a variable across the treatment and control groups. Currently, the 
    \code{twang} package has functions \code{\link{es.stat}}, 
    \code{\link{ks.stat}}, and \code{\link{strata.stat}}. The user may 
    implement their own.}
\item{rule.summary}{a function that takes the vector of results from the 
    \code{metric} function and summarizes them into a single number. 
    \code{twang} currently utilizes \code{\link{mean}} and \code{\link{max}} 
    for \code{rule.summary}}
\item{direct}{logical. If \code{TRUE} then \code{\link{ps}} will try to 
    optimize the weights directly rather than utilizing \code{\link[gbm]{gbm}}}
\item{na.action}{a character string indicating how \code{bal.stat} should 
    handle missing values. Current options are "level", 
    "exclude", or "lowest"}
\item{name}{a character string, preferably unique from other 
    \code{stop.method}s for labeling the resulting weights}
}
In addition, the object must have \code{class(mystopmethod) == "stop.method"}
}

\keyword{multivariate}

\name{ks.stat}
\alias{ks.stat}
\alias{es.stat}
\alias{strata.stat}

\title{ Functions for evaluating balance}
\description{
These are a collection of functions that can be used as components of 
\code{\link{stop.methods}} for evaluating the balance of two groups
}
\usage{
ks.stat(logw = NULL, w.ctrl = NULL, 
        gbm1 = NULL, i = 1, data, 
        sampw = rep(1, nrow(data)), 
        rule.summary = mean, na.action = "level", 
        vars, treat.var, collapse.by.var = FALSE, 
        verbose = FALSE)

es.stat(logw = NULL, w.ctrl = NULL, 
        gbm1 = NULL, i = 1, data, 
        sampw = rep(1, nrow(data)), 
        rule.summary = mean, na.action = "level", 
        vars, treat.var, collapse.by.var = FALSE, 
        verbose = FALSE)

strata.stat(logw = NULL, w.ctrl = NULL, 
            gbm1 = NULL, i = 1, data, 
            sampw = rep(1, nrow(data)), 
            rule.summary = mean, na.action = "level", 
            vars, treat.var, collapse.by.var = FALSE, 
            verbose = FALSE)
}
\arguments{
The weights be passed to these functions with any of the first three arguments
  \item{logw}{the logarithm of the weights }
  \item{w.ctrl}{the weights for the control subjects }
  \item{gbm1}{a \code{\link[gbm]{gbm.object}} used for estimating the 
              propensity scores, usually the \code{gbm} component of a \code{ps} 
              object returned from \code{\link{ps}}}
  \item{i}{the iteration of \code{\link[gbm]{gbm}} with which to compute the 
           weights}
  \item{data}{ a data frame with the data }
  \item{sampw}{ optional sampling weights }
  \item{rule.summary}{a function for summarizing the total balance. Used to 
                      collapse statistics across all the covariates. Examples 
                      include \code{mean} and \code{max} }
  \item{na.action}{a string indicating the method for handling missing data}
  \item{vars}{a vector of variable names corresponding to \code{data}}
  \item{treat.var}{the name of the treatment variable}
  \item{collapse.by.var}{if \code{TRUE}, then statistics computed for factors 
                         are collapsed across the levels}
  \item{verbose}{if TRUE, lots of information will be printed to monitor the
                 the progress of the fitting}
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}

\seealso{ \code{\link{stop.methods}} }

\keyword{ models }

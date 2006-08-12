\name{dx.wts}
\alias{dx.wts}

\title{Propensity score diagnostics}
\description{
\code{dx.wts} takes a \code{ps} object or a set of propensity scores and 
computes diagnostics assessing covariates balance.
}
\usage{
dx.wts(x,
       data,
       vars=NULL,
       treat.var,
       x.as.weights=TRUE,
       sampw=NULL,
       perm.test.iters=0,
       plots=TRUE,
       title)
}

\arguments{
  \item{x}{a data frame, matrix, or vector of propensity score weights or a ps 
           object. \code{x} can also be a data frame, matrix, or vector of 
           propensity scores if \code{x.as.weights=FALSE}}
  \item{data}{a data frame}
  \item{vars}{a vector of character strings naming variables in \code{data} on 
              which to assess balance}
  \item{treat.var}{a character string indicating which variable in \code{data}
                   contains the 0/1 treatment group indicator }
  \item{x.as.weights}{\code{TRUE} or \code{FALSE} indicating whether \code{x} 
                      specifies propensity score weights or propensity scores. 
                      Ignored if \code{x} is a ps object}
  \item{sampw}{optional sampling weights. If \code{x} is a ps object then the 
              sampling weights  should have been passed to \code{\link{ps}} and 
              not specified here. \code{dx.wts} will issue a warning if 
              \code{x} is a ps object and \code{sampw} is also specified}
  \item{perm.test.iters}{an non-negative integer giving the number of iterations
        of the permutation test for the KS statistic. If \code{perm.test.iters=0}
        then the function returns an analytic approximation to the p-value. This
        argument is ignored is \code{x} is a \code{ps} object. Setting
        \code{perm.test.iters=200} will yield precision to within 3\% if the true
        p-value is 0.05. Use \code{perm.test.iters=500} to be within 2\%}
  \item{plots}{if \code{plots=TRUE} then \code{dx.wts} will call 
               \code{\link{diag.plot}} generating diagnostic plots}
  \item{title}{a short text title, it will be used in plots and saved files. 
               By default this is set to the current date and time}
}
\details{
Creates a balance table that compares unweighted and weighted means and 
standard deviations, computes effect sizes, and KS statistics to assess the 
ability of the propensity scores to balance the treatment and control groups.
}
\value{
Returns a list containing
  \item{treat}{the vector of 0/1 treatment assignment indicators}
  \item{desc}{a nested list containing detailed diagnostic information on the
       weights. This includes the number of treatment and control subjects, the
       effective sample size, the largest KS statistic, the average absolute
       effect size, and the complete balance table}
  \item{summary.tab}{a data frame showing balance information}
  \item{ps}{the given propensity scores}
  \item{w}{the given weights}
  \item{datestamp}{the date and time of the call to \code{dx.wts}}
  \item{parameters}{the parameters used when calling \code{dx.wts}}
  \item{alerts}{text containing any warnings accumulated during the estimation}
}

\seealso{ The example for \code{\link{ps}} contains an example of the 
use of \code{dx.wts}, \code{\link{diag.plot}}}

\keyword{ models }

\name{ps.summary}
\alias{ps.summary}
\alias{ps.summary.f}
\alias{ps.summary.n}

\title{Computes balance measures}
\description{
Computes balance measures (mean differences and KS statistics) for a particular 
covariate and a set of propensity score weights. This function is not intended
to be called directly by the user but is used by other functions in the package.
}
\usage{
ps.summary(x, t, w, get.means = TRUE, get.ks = TRUE, 
           na.action = c("level", "exclude", "lowest")[1], 
           collapse.by.var = FALSE)
ps.summary.f(x, t, w, get.means = TRUE, get.ks = TRUE, 
             na.action = c("level", "exclude", "lowest")[1], 
             collapse.by.var = TRUE)
ps.summary.n(x, t, w, get.means = TRUE, get.ks = TRUE, 
             na.action = c("level", "exclude", "lowest")[1], 
             collapse.by.var = FALSE)
}


\arguments{
  \item{x}{a vector containing the data for a single covariate}
  \item{t}{a vector of the same length as \code{x} with the 0/1 treatment 
           assigniments}
  \item{w}{a vector of the same length as \code{x} with the weights}
  \item{get.means}{if \code{TRUE}, mean comparisons are computed}
  \item{get.ks}{if \code{TRUE}, the KS statistics are computed}
  \item{na.action}{a string indicating the method for handling missing data}
  \item{collapse.by.var}{if \code{TRUE}, then statistics computed for factors 
                         are collapsed across the levels}
}
\details{
\code{ps.summary} dispatches \code{ps.summary.n} or \code{ps.summary.f} 
depending on whether \code{x} is a numeric vector or a factor.


}
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
                 
  \code{get.means} and \code{get.ks} manipulate the inclusion of certain 
columns in the returned result.
}

\seealso{\code{\link[twang]{bal.stat}}, 
         \code{\link[twang]{ks.stat}}, 
         \code{\link[twang]{es.stat}}}
\examples{
treat <- rbinom(100,1,0.5)
w     <- rexp(100)

# categorical data
x.cat <- factor(sample(letters[1:3],size=100,replace=TRUE))
ps.summary.f(x.cat,treat,w)

# numeric data
x.num <- rnorm(100)
ps.summary.f(x.num,treat,w)

# or let ps.summary figure out which to call
ps.summary(x.num,treat,w)
}
\keyword{ models }

\name{sensitivity}
\alias{sensitivity}

\title{Sensitivity analysis} 

\description{ Produces a table to help the user assess the extent to which a 
hidden bias might remove any differences observed in the propensity score 
analysis.}

\usage{
sensitivity(ps1,
            data,
            outcome,
            order.by.importance = TRUE,
            verbose = TRUE)
}

\arguments{
  \item{ps1}{ a \code{ps} object as returned from \code{\link{ps}} }
  \item{data}{ the data frame used to fit \code{ps1} }
  \item{outcome}{ a character string indicating the name of the variable in
                  \code{data} to use as the outcome }
  \item{order.by.importance}{ if \code{TRUE} then the variables are sorted by
        their relative influence in the \code{\link[gbm]{gbm.object}} used
        to create \code{ps1} }
  \item{verbose}{ if TRUE, lots of information will be printed to monitor the
                 the progress of the fitting }
} 

\details{ This function implements the sensitivity analysis described in 
Ridgeway (2006), Section 5.5. This analysis helps the user assess the extent 
to which a hidden bias might remove any differences observed in the 
propensity score analysis.

If there is an important unobserved factor the odds than the correct 
propensity score weight is not \eqn{w(x_i)}, as the propensity score model 
predicts, but actually \eqn{w(x_i, z_i)} where \eqn{z} represents the 
unobserved factor. Let \eqn{a_i=w(x_i, z_i)/w(x_i)}. These \eqn{a_i}'s give 
an estimate of \eqn{g(a)}, the distribution of the multiplicative errors that 
we observe in the weights when excluding \eqn{z_i}. Changing the values of 
the \eqn{a_i}'s will affect the treatment effect estimate if \eqn{a} is 
correlated with \eqn{y}, the outcome. The stronger the correlation the more 
sensitive the results will be to the hidden bias. \code{sensitivity} computes 
over control group subjects a modified estimate of \eqn{E(Y_0|t=1)}.

\deqn{\frac{\sum_C a_iw_iy_i}{\sum_C a_iw_i}}{(sum a_i*w_i*y_i)/(sum a_i*w_i)}

subject to the constraint that \eqn{a_i \sim g(a)}{a_i ~ g(a)} and 
\eqn{cor(a_i, y_i) = \rho}.

Several \eqn{g(a)}'s are considered by removing each variable from the 
propensity score model in turn and computing the ratio of the original 
weights to the weights with the variable removed. Several choices for 
\eqn{\rho} are also considered, making \eqn{\rho} as large as possible, as 
small as possible, and solving for the ``break even'' \eqn{\rho}, the 
\eqn{\rho} that eliminates any treatment effect. }

\value{ Returns a list where each component contains the sensitivity analysis 
for each \code{stop.method} used in fitting \code{ps1}. Each component 
contains a data frame with a row for each variable in the original propensity 
score model. The columns are
  \item{var}{the name of the variable excluded from the model}
  \item{E0}{the estimated \eqn{E(Y_0|t=1)} with \code{var} excluded from the 
        propensity score model}
  \item{a.min,a.max}{the smallest and largest values of \eqn{a} observed}
  \item{a.cor}{the observed correlation between \eqn{a} and \eqn{\rho}}
  \item{a.mincor,a.maxcor}{the smallest and largest values of \eqn{\rho} 
        possible}
  \item{minE0,maxE0}{the smallest and largest values of estimated 
        \eqn{E(Y_0|t=1)} possible}
  \item{breakeven.cor}{the break even correlation (see Details section)}
} 

\references{ G. Ridgeway (2006). \dQuote{Assessing the effect of race bias in 
post-traffic stop outcomes using propensity scores,} \emph{Journal of 
Quantitative Criminology} 22(1). } 

\author{ Greg Ridgeway \email{gregr@rand.org} }

\seealso{ See \code{\link{ps}} for an example } 

\keyword{ models }

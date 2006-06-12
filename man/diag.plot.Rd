\name{diag.plot}
\alias{diag.plot}
\alias{plot.dxwts}
\alias{plot.ps}

\title{ Create diagnostic plots }
\description{
Creates diagnostic plots of propensity scores including, a side-by-side boxplot
of propensity scores, a histogram of propensity score weights, QQ plots of KS 
and t-statistic p-values, and a plot of absolute effect sizes
}
\usage{
diag.plot(title=NULL,
          treat=NULL,
          p.s=NULL,
          w.ctrl=NULL,
          desc.unw=NULL,
          desc.w=NULL,
          plots="all")

\method{plot}{ps}(x, label = "", ask=FALSE, plots="all", ...)
                  
\method{plot}{dxwts}(x, label = "", ask=FALSE, plots="all", ...)
}

\arguments{
  \item{title}{a title for the plots}
  \item{treat}{a vector of 0/1 treatment indicators}
  \item{p.s}{a vector of propensity scores (optional)}
  \item{w.ctrl}{weights for the control subjects}
  \item{desc.unw}{a list object containing the balance assessment without 
                  weights, usually the result of a call to 
                  \code{\link{desc.wts}} or the \code{desc} component 
                  of a \code{\link{ps}} object}
  \item{desc.w}{a list object containing the weighted balance assessment}
  \item{plots}{a character vector listing the plots to be created. The options 
    are all (the default), optimize, ps boxplot, weight histogram, 
    t pvalues, ks pvalues, es. Any other options (such as "none") will
    produce no plots}
  
  \item{x}{ a \code{ps} object, usually one returned from \code{\link{ps}} }
  \item{label}{ a character string for titling the plots }
  \item{ask}{logical. If \code{TRUE} then the graphics window waits for a 
             response from the user before showing the next graph}
  \item{\dots}{ other arguments passed to the plot function }  
}
\details{
  \code{plot.ps} and \code{plot.dxwts} are wrappers for \code{\link{diag.plot}}
  
The plots include
    \item Boxplot of propensity scores for cases in the treatment and 
comparison conditions
    \item Histogram of comparison condition case weights
    \item P-value plots for unweighted and weighted t statistics, and KS 
statistics
    \item Change in standardized effect size plot. For each model 
covariate, standardized effect sizes before and 
after comparison group case weighting are linked by blue lines if 
weighting reduces the effect size, and by red lines if weights increase 
the effect size.  Significant effect sizes are indicated with a  closed 
red circle.  Standardized effect sizes are defined as the difference 
between the treatment and comparison group means, divided by the treatment 
group standard deviation.  Very large effect sizes are omitted from these 
plots.  When this occurs, a warning is placed at the top of the figure 
}
\value{
No returned objects
}

\seealso{ \code{\link{ps}} }

\keyword{ hplot }

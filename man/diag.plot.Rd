\name{diag.plot}
\alias{diag.plot}
\alias{plot.dxwts}
\alias{plot.ps}

\title{ Create diagnostic plots }
\description{
Creates diagnostic plots of propensity scores including, a side-by-side boxplot
of propensity scores, a histogram of propensity score weights, QQ plots of KS 
and t-statistic p-values, and a ``spaghetti" plot of absolute effect sizes
}
\usage{
diag.plot(til, treat, p.s, w.ctrl, desc.unw, desc.w)

\method{plot}{ps}(x, label = "", ask = FALSE, ...)
                  
\method{plot}{dxwts}(x, label = "", ...)
}

\arguments{
  \item{til}{a title for the plots}
  \item{treat}{a vector of 0/1 treatment indicators}
  \item{p.s}{a vector of propensity scores (optional)}
  \item{w.ctrl}{weights for the control subjects}
  \item{desc.unw}{a list object containing the balance assessment without 
                  weights, usually the result of a call to 
                  \code{\link[twang]{desc.wts}} or the \code{desc} component 
                  of a \code{\link[twang]{ps}} object}
  \item{desc.w}{a list object containing the weighted balance assessment}
  
  \item{x}{ a \code{ps} object, usually one returned from \code{\link[twang]{ps}} }
  \item{label}{ a character string for titling the plots }
  \item{ask}{logical. If \code{TRUE} then the graphics window waits for a 
             response from the user before showing the next graph}
  \item{\dots}{ other arguments passed to the plot function }  
}
\details{
  \code{plot.ps} and \code{plot.dxwts} are wrappers for \code{\link{diag.plot}}
}
\value{
No returned objects
}

\seealso{ \code{\link[twang]{ps}} }

\keyword{ hplot }

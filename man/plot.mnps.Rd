\name{plot.mnps}
\alias{plot.mnps}

\title{Plots for mnps objects}

\description{This function produces a collection of diagnostic plots for ps objects. }

\usage{ 
\method{plot}{mnps}(x, plots = "optimize", figureRows = 1, summaryFcn = max, \dots)
}

\arguments{
  \item{x}{An mnps object.}
  \item{plots}{An indicator of which type of plot is desired.  The options are
\describe{
\item{"optimize" or 1}{A plot of the balance criteria as a function of the GBM iteration}
\item{"boxplot" or 2}{Boxplots of the propensity scores for the treatment and control cases}
\item{"es" or 3}{Plots of the standardized effect size of the pre-treatment variables before and after reweighing}
\item{"t" or 4}{Plots of the p-values from t-statistics comparing means of treated and control subjects for pretreatment variables, before and after weighting.}
\item{"ks" or 5}{Plots of the p-values from Kolmogorov-Smirnov statistics comparing distributions of pretreatment variables of treated and control subjects, before and after weighting.  Not available when the estimand is ATE.}
\item{"histogram" or 6}{Histogram of weights for treated and control subjects.}
}
}
\item{figureRows}{The number of rows of figures that should be used.}
\item{summaryFcn}{Function that summarizes across the multiple \code{ps} fits.  Specifying \code{NULL} forgoes the summarizing process and returns one frame for each treatment category.} 
  \item{...}{Additional arguments that may be passed to the underlying \code{lattice} package plotting functions}

}


\details{
This function produces lattice-style graphics of diagnostic plots.}

\references{
Dan McCaffrey, G. Ridgeway, Andrew Morral (2004). \dQuote{Propensity Score Estimation
with Boosted Regression for Evaluating Adolescent Substance Abuse Treatment,}
\emph{Psychological Methods} 9(4):403-425.
}


\seealso{ \code{\link{mnps}} }


\keyword{models}
\keyword{multivariate}

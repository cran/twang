\name{plot.mnps}
\alias{plot.mnps}

\title{Plots for mnps objects}

\description{This function produces a collection of diagnostic plots for ps objects. }

\usage{ 
\method{plot}{mnps}(x, plots = "optimize", 
	pairwiseMax = TRUE, figureRows = 1,  
	color = TRUE, subset = NULL, \dots)
}

\arguments{
  \item{x}{An mnps object.}
  \item{plots}{An indicator of which type of plot is desired.  The options are
\describe{
\item{"optimize" or 1}{A plot of the balance criteria as a function of the GBM 
	iteration}
\item{"boxplot" or 2}{Boxplots of the propensity scores for the treatment and 
	control cases}
\item{"es" or 3}{Plots of the standardized effect size of the pre-treatment 
	variables before and after reweighing}
\item{"t" or 4}{Plots of the p-values from t-statistics comparing means of 
	treated and control subjects for pretreatment variables, before and after weighting.}
\item{"ks" or 5}{Plots of the p-values from Kolmogorov-Smirnov statistics 
	comparing distributions of pretreatment variables of treated and control 
	subjects, before and after weighting.}
}
}
\item{pairwiseMax}{If \code{FALSE}, the plots for the underlying \code{ps} fits 
	will be returned.  Otherwise, pairwise maxima will be returned.}
\item{figureRows}{The number of rows of figures that should be used.}
\item{color}{If \code{color = FALSE} figures will be gray scale.}
\item{subset}{Used to restrict which of the \code{stop.method}s will be used 
	in the figure. For example \code{subset = c(1,3)} would indicate that the 
	first and third \code{stop.method}s (in alphabetical order of those specified 
	in the original call to \code{mnps}) should be included in the figure.}
  \item{...}{Additional arguments that may be passed to the underlying \code{lattice} 
  	package plotting functions}

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

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.ps.R
\name{plot.ps}
\alias{plot.ps}
\title{Plots for \code{ps} objects}
\usage{
\method{plot}{ps}(x, plots = "optimize", subset = NULL, color = TRUE, ...)
}
\arguments{
\item{x}{A \code{ps} object.}

\item{plots}{An indicator of which type of plot is desired. The options are
\itemize{
\item \verb{"optimize" or 1} A plot of the balance criteria as a function of the GBM
iteration.
\item \verb{"boxplot" or 2} Boxplots of the propensity scores for the treatment and
control cases
\item \verb{"es" or 3} Plots of the standardized effect size of the pre-treatment
variables before and after reweighing
\item \verb{"t" or 4} Plots of the p-values from t-statistics comparing means of
treated and control subjects for pretreatment variables, before and after
weighting.
\item \verb{"ks" or 5} Plots of the p-values from Kolmogorov-Smirnov statistics
comparing distributions of pretreatment variables of treated and control
subjects, before and after weighting.
}}

\item{subset}{If multiple \code{stop.method} rules were used in the \code{ps()} call,
\code{subset} restricts the plots of a subset of the stopping rules that were
employed. This argument expects a subset of the integers from 1 to k,
if k \code{stop.method}s were used.}

\item{color}{If \code{color = FALSE}, figures will be gray scale. Default: \code{TRUE}.}

\item{...}{Additional arguments.}
}
\description{
This function produces a collection of diagnostic plots for \code{ps} objects.
}
\details{
This function produces lattice-style graphics of diagnostic plots.
}
\references{
Dan McCaffrey, G. Ridgeway, Andrew Morral (2004). "Propensity
Score Estimation with Boosted Regression for Evaluating Adolescent
Substance Abuse Treatment", \emph{Psychological Methods} 9(4):403-425.
}
\seealso{
\link{ps}
}

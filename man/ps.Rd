\name{ps}
\alias{ps}

\title{Propensity score estimation}

\description{ \code{ps} calculates propensity scores and diagnoses them using 
a variety of methods, but centered on using boosted logistic regression as 
implemented in \code{\link[gbm]{gbm}} }

\usage{ 
ps(formula = formula(data),
   data,
   sampw = rep(1, nrow(data)),
   title=NULL,
   stop.method = stop.methods[1:2],
   plots="all",
   pdf.plots=FALSE,
   n.trees = 10000,
   interaction.depth = 3,
   shrinkage = 0.01,
   perm.test.iters=0,
   print.level = 2,
   iterlim = 1000,
   verbose = TRUE)
}

\arguments{
  \item{formula}{a formula for the propensity score model with the treatment
                 indicator on the left side of the formula and the potential
                 confounding variables on the right side.}
  \item{title}{a short text title, it will be used in plots and saved files}
  \item{data}{the dataset, includes treatment assignment as well as covariates}
  \item{sampw}{optional sampling weights}
  \item{stop.method}{a \code{\link{stop.methods}} object, or a list of such
                     objects, containing the metrics and rules for evaluating
                     the quality of the propensity scores}
  \item{plots}{a character vector indicating which plots to create. The options 
    are all (the default), optimize, ps boxplot, weight histogram, 
    t pvalues, ks pvalues, es. Any other options (such as "none") will
    produce no plots}
  \item{pdf.plots}{if \code{TRUE} then all plots are dumped to a pdf file with
    the name specified in \code{title}}
  \item{n.trees}{ number of gbm iterations passed on to \code{\link[gbm]{gbm}} }
  \item{interaction.depth}{ \code{interaction.depth} passed on to
        \code{\link[gbm]{gbm}} }
  \item{shrinkage}{ \code{shrinkage} passed on to \code{\link[gbm]{gbm}} }
  \item{perm.test.iters}{a non-negative integer giving the number of iterations
        of the permutation test for the KS statistic. If \code{perm.test.iters=0}
        then the function returns an analytic approximation to the p-value. Setting
        \code{perm.test.iters=200} will yield precision to within 3\% if the true
        p-value is 0.05. Use \code{perm.test.iters=500} to be within 2\%}
  \item{print.level}{ the amount of detail to print to the screen }
  \item{iterlim}{ maximum number of iterations for the direct optimization }
  \item{verbose}{ if TRUE, lots of information will be printed to monitor the
                 the progress of the fitting }
}


\details{
\code{formula} should be something like "treatment ~ X1 + X2 + X3". The
treatment variable should be a 0/1 indicator. There is no need to specify
interaction terms in the formula. \code{interaction.depth} controls the level
of interactions to allow in the propensity score model.

If \code{pdf.plots=TRUE} then \code{ps} causes plots to be saved as a single 
pdf file with the name "[title].pdf" in the working directory. See 
\code{\link{diag.plot}} for details of the plots.}

\value{ 
Returns an object of class \code{ps}, a list containing 
  \item{gbm.obj}{The returned \code{\link[gbm]{gbm}} object}
  \item{ps}{a data frame containing the estimated propensity scores. Each
            column is associated with one of the methods selected in
            \code{stop.methods}}
  \item{w}{a data frame containing the propensity score weights. Each
            column is associated with one of the methods selected in
            \code{stop.methods}. If sampling weights were given then these are
            incorporated into these weights}
  \item{plot.info}{a list containing the raw data used to generate the plots}
  \item{desc}{a list containing balance tables for each method selected in
              \code{stop.methods}. Includes a component for the unweighted
              analysis names \dQuote{unw}. Each \code{desc} component includes
              a list with the following components
  \describe{
     \item{ess}{The effective sample size of the control group}
     \item{n.treat}{The number of subjects in the treatment group}
     \item{n.ctrl}{The number of subjects in the control group}
     \item{max.es}{The largest effect size across the covariates}
     \item{mean.es}{The mean absolute effect size}
     \item{max.ks}{The largest KS statistic across the covariates}
     \item{mean.ks}{The average KS statistic across the covariates}
     \item{bal.tab}{a (potentially large) table summarizing the quality of the 
                    weights for equalizing the distribution of features across 
                    the two groups. This table is best extracted using the
                    \code{\link{bal.table}} method. See the help for 
                    \code{\link{bal.table}} for details on the table's contents}
     \item{n.trees}{The estimated optimal number of \code{\link[gbm]{gbm}} 
                    iterations to optimize the loss function for the associated 
                    \code{\link{stop.methods}}}
     }      
  }
  \item{datestamp}{Records the date of the analysis}
  \item{parameters}{Saves the \code{ps} call}
  \item{alerts}{Text containing any warnings accumulated during the estimation}
 }

\references{
Dan McCaffrey, G. Ridgeway, Andrew Morral (2004). \dQuote{Propensity Score Estimation
with Boosted Regression for Evaluating Adolescent Substance Abuse Treatment,}
\emph{Psychological Methods} 9(4):403-425.
}

\author{
Greg Ridgeway \email{gregr@rand.org},
Dan McCaffrey \email{danielm@rand.org},
Andrew Morral \email{morral@rand.org}
}

\seealso{ \code{\link[gbm]{gbm}} }

\examples{
data(lalonde)
print(nrow(lalonde))

ps.lalonde <- ps(treat ~ age + educ + black + hispan + nodegree + 
                         married + re74 + re75, 
                 data = lalonde,
                 title="Lalonde example",
                 stop.method=stop.methods[c("ks.stat.mean","ks.stat.max")],
                 # generate plots?
                 plots="all",
                 pdf.plots=FALSE,
                 # gbm options
                 n.trees=2000,
                 interaction.depth=3,
                 shrinkage=0.005,
                 perm.test.iters=0,
                 verbose=TRUE)
                 
# get the balance tables
bal.table(ps.lalonde)

# diagnose the weights using a ps object 
a <- dx.wts(ps.lalonde,data=lalonde,treat.var="treat")
print(a)
bal.table(a)

# diagnose the weights as propensity score weights
#    will be the same as before, except for MC variation in the KS p-values
#    when perm.test.iters is greater than 0
w <- with(ps.lalonde, ps/(1-ps))
w[lalonde$treat==1,] <- 1
dx.wts(w,data=lalonde,treat.var="treat",
       perm.test.iters=0)

# diagnose the weights as propensity scores
p <- ps.lalonde$ps
dx.wts(p,data=lalonde,treat.var="treat",x.as.weights=FALSE)

# look at propensity scores
names(ps.lalonde$ps)
hist(ps.lalonde$ps$ks.stat.max)
boxplot(split(ps.lalonde$ps$ks.stat.max,ps.lalonde$treat),
        ylab="estimated propensity scores",
        names=c("control","treatment"))

# check out the balance
names(ps.lalonde$desc)
# unweighted
ps.lalonde$desc$unw
# optimized for ks.stat.max
ps.lalonde$desc$ks.stat.max

# check out the gbm object, indicates which variables are most influential in 
#    estimating the propensity score
summary(ps.lalonde$gbm.obj, n.trees=ps.lalonde$desc$ks.stat.max$n.trees)

# bal.stat() can use an arbitrary set of weights
bal.stat(data=lalonde,
         w.all=w[,1],
         vars=names(lalonde),
         treat.var="treat",
         get.means=TRUE,
         get.ks=TRUE,
         na.action="level")
         
# sensitivity analysis
sensitivity(ps.lalonde,lalonde,"re78")
}
\keyword{models}
\keyword{multivariate}

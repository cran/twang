#' Gradient boosted propensity score estimation
#'
#' \code{ps} calculates propensity scores using gradient boosted logistic
#' regression and diagnoses the resulting propensity scores using a variety of
#' methods
#' 
#' For user more comfortable with the options of \code{\link[xgboost:xgboost]{xgboost::xgboost()}},
#' the options for \code{ps} controlling the behavior of the gradient boosting
#' algorithm can be specified using the \code{\link[xgboost]{xgboost}} naming
#' scheme. This includes \code{nrounds}, \code{max_depth}, \code{eta}, and
#' \code{subsample}. In addition, the list of parameters passed to
#' \code{\link[xgboost]{xgboost}} can be specified with \code{params}.
#' 
#' Note that unlike earlier versions of `twang`, the plotting functions are
#' no longer included in the \code{ps} function. See  \code{\link[twang:plot.ps]{plot}} for
#' details of the plots.
#' 
#' @param formula An object of class \code{\link{formula}}: a symbolic
#'   description of the propensity score model to be fit with the treatment
#'   indicator on the left side of the formula and the potential confounding
#'   variables on the right side.
#' @param data A dataset that includes the treatment indicator as well as the
#'   potential confounding variables.
#' @param n.trees Number of gbm iterations passed on to \code{\link[gbm:gbm]{gbm::gbm()}}. Default: 10000.
#' @param interaction.depth A positive integer denoting the tree depth used in
#'   gradient boosting. Default: 3.
#' @param shrinkage A numeric value between 0 and 1 denoting the learning rate.
#'   See \code{\link[gbm]{gbm}} for more details. Default: 0.01.
#' @param bag.fraction A numeric value between 0 and 1 denoting the fraction of
#'   the observations randomly selected in each iteration of the gradient
#'   boosting algorithm to propose the next tree. See \code{\link[gbm]{gbm}} for
#'   more details. Default: 1.0.
#' @param n.minobsinnode An integer specifying the minimum number of observations 
#'   in the terminal nodes of the trees used in the gradient boosting.  See 
#'   \code{\link[gbm]{gbm}} for
#'   more details. Default: 10.
#' @param perm.test.iters A non-negative integer giving the number of iterations
#'   of the permutation test for the KS statistic. If \code{perm.test.iters=0}
#'   then the function returns an analytic approximation to the p-value. Setting
#'   \code{perm.test.iters=200} will yield precision to within 3\% if the true
#'   p-value is 0.05. Use \code{perm.test.iters=500} to be within 2\%. Default: 0.
#' @param print.level The amount of detail to print to the screen. Default: 2.
#' @param verbose If \code{TRUE}, lots of information will be printed to monitor the
#'   the progress of the fitting. Default: \code{TRUE}.
#' @param estimand \code{"ATE"} (average treatment effect) or \code{"ATT"} (average treatment
#'   effect on the treated) : the causal effect of interest. ATE estimates the
#'   change in the outcome if the treatment were applied to the entire
#'   population versus if the control were applied to the entire population. ATT
#'   estimates the analogous effect, averaging only over the treated population.
#'   Default: \code{"ATE"}.
#' @param stop.method A method or methods of measuring and summarizing balance across pretreatment
#'   variables. Current options are \code{ks.mean}, \code{ks.max}, \code{es.mean}, and \code{es.max}. \code{ks} refers to the
#'   Kolmogorov-Smirnov statistic and es refers to standardized effect size. These are summarized
#'   across the pretreatment variables by either the maximum (\code{.max}) or the mean (\code{.mean}). 
#'   Default: \code{c("ks.mean", "es.mean")}.
#' @param sampw Optional sampling weights.
#' @param version \code{"gbm"}, \code{"xgboost"}, or \code{"legacy"}, indicating which version of the twang package to use.
#'   \describe{
#'     \item{\code{"gbm"}}{ uses gradient boosting from the \code{\link[gbm]{gbm}} package,}
#'     \item{\code{"xgboost"}}{ uses gradient boosting from the \code{\link[xgboost]{xgboost}} package, and}
#'     \item{\code{"legacy"}}{ uses the prior implementation of the \code{ps} function.}
#'   }
#'   Default: \code{"gbm"}.
#' @param ks.exact \code{NULL} or a logical indicating whether the
#'   Kolmogorov-Smirnov p-value should be based on an approximation of exact
#'   distribution from an unweighted two-sample Kolmogorov-Smirnov test. If
#'   \code{NULL}, the approximation based on the exact distribution is computed
#'   if the product of the effective sample sizes is less than 10,000.
#'   Otherwise, an approximation based on the asymptotic distribution is used.
#'   **Warning:** setting \code{ks.exact = TRUE} will add substantial
#'   computation time for larger sample sizes. Default: \code{NULL}.
#' @param n.keep A numeric variable indicating the algorithm should only
#'   consider every \code{n.keep}-th iteration of the propensity score model and
#'   optimize balance over this set instead of all iterations. Default: 1.
#' @param n.grid A numeric variable that sets the grid size for an initial
#'   search of the region most likely to minimize the \code{stop.method}. A
#'   value of \code{n.grid=50} uses a 50 point grid from \code{1:n.trees}. It
#'   finds the minimum, say at grid point 35. It then looks for the actual
#'   minimum between grid points 34 and 36. If specified with \code{n.keep>1}, \code{n.grid} 
#'   corresponds to a grid of points on the kept iterations as defined by \code{n.keep}. Default: 25.
#' @param keep.data A logical variable indicating whether or not the data is saved in 
#'   the resulting \code{ps} object. Default: \code{TRUE}.
#' @param ... Additional arguments that are passed to \code{ps} function.
#' 
#' @return Returns an object of class \code{ps}, a list containing 
#'   \describe{
#'   \item{\code{gbm.obj}}{ The returned \code{\link[gbm]{gbm}} or \code{\link[xgboost]{xgboost}} object.}
#'   \item{\code{treat}}{ The vector of treatment indicators.}
#'   \item{\code{treat.var}}{ The treatment variable.}
#'   \item{\code{desc}}{ A list containing balance tables for each method selected in
#'     \code{stop.methods}. Includes a component for the unweighted
#'     analysis names \dQuote{unw}. Each \code{desc} component includes
#'     a list with the following components
#'     \describe{
#'     \item{\code{ess}}{ The effective sample size of the control group.}
#'     \item{\code{n.treat}}{ The number of subjects in the treatment group.}
#'     \item{\code{n.ctrl}}{ The number of subjects in the control group.}
#'     \item{\code{max.es}}{ The largest effect size across the covariates.}
#'     \item{\code{mean.es}}{ The mean absolute effect size.}
#'     \item{\code{max.ks}}{ The largest KS statistic across the covariates.}
#'     \item{\code{mean.ks}}{ The average KS statistic across the covariates.}
#'     \item{\code{bal.tab}}{ a (potentially large) table summarizing the quality of the 
#'       weights for equalizing the distribution of features across 
#'       the two groups. This table is best extracted using the
#'       \code{\link{bal.table}} method. See the help for \code{\link{bal.table}} for details
#'       on the table's contents.}
#'     \item{\code{n.trees}}{ The estimated optimal number of gradient boosted
#'       iterations to optimize the loss function for the associated 
#'        \code{stop.methods}.}
#'     \item{\code{ps}}{ a data frame containing the estimated propensity scores. Each
#'       column is associated with one of the methods selected in \code{stop.methods}.}
#'     \item{\code{w}}{ a data frame containing the propensity score weights. Each
#'       column is associated with one of the methods selected in \code{stop.methods}.}
#'       If sampling weights are given then these are incorporated into these weights.
#'     \item{\code{estimand}}{ The estimand of interest (ATT or ATE).}
#'     }}
#'  \item{\code{datestamp}}{ Records the date of the analysis.}
#'  \item{\code{parameters}}{ Saves the \code{ps} call.}
#'  \item{\code{alerts}}{ Text containing any warnings accumulated during the estimation.}
#'  \item{\code{iters}}{ A sequence of iterations used in the GBM fits used by \code{\link[twang:plot.ps]{plot}} function.}
#'  \item{\code{balance}}{ The balance measures for the pretreatment covariates used in plotting, with a column for each
#'    \code{stop.method}.}
#'  \item{\code{balance.ks}}{ The KS balance measures for the pretreatment covariates used in plotting, with a column for each
#'    covariate.}
#'  \item{\code{balance.es}}{ The standard differences for the pretreatment covariates used in plotting, with a column for each
#'    covariate.}
#'  \item{\code{ks}}{ The KS balance measures for the pretreatment covariates on a finer grid, with a column for each
#'    covariate.}
#'  \item{\code{es}}{ The standard differences for the pretreatment covariates on a finer grid, with a column for each
#'    covariate.}
#'  \item{\code{n.trees}}{ Maximum number of trees considered in GBM fit.}
#'  \item{\code{data}}{ Data as specified in the \code{data} argument.}
#'  }
#'
#' @seealso \code{\link[gbm]{gbm}}, \code{\link[xgboost]{xgboost}}, \code{\link[twang:plot.ps]{plot}}, \code{\link{bal.table}}
#' @keywords models multivariate
#'
#' @references Dan McCaffrey, G. Ridgeway, Andrew Morral (2004). "Propensity
#'   Score Estimation with Boosted Regression for Evaluating Adolescent
#'   Substance Abuse Treatment", *Psychological Methods* 9(4):403-425.
#'
#' @export
ps<-function(formula = formula(data),
             data,
             # boosting options 
             n.trees = 10000,
             interaction.depth = 3,
             shrinkage = 0.01,
             bag.fraction = 1.0,
             n.minobsinnode =10,
             perm.test.iters = 0,
             print.level = 2,       
             verbose = TRUE,
             estimand = "ATE", 
             stop.method = c("ks.mean", "es.mean"), 
             sampw = NULL, 
             version = "gbm",
             ks.exact = NULL,
             n.keep = 1,
             n.grid = 25,
             keep.data=TRUE,
             ...){

   # collect named arguments from dots
   args         <- list(...)
   args_named   <- names(args)

   # parse hidden options and xgboost parameters
   params       <- args$params
   multinom     <- if (!is.null(args$multinom)) args$multinom else FALSE
   max_depth    <- if (!is.null(args$max_depth)) args$max_depth else interaction.depth
   subsample    <- if (!is.null(args$subsample)) args$subsample else bag.fraction
   nrounds      <- if (!is.null(args$nrounds)) args$nrounds else n.trees
   eta          <- if (!is.null(args$eta)) args$eta else shrinkage
   min_child_weight <- if (!is.null(args$min_child_weight)) args$min_child_weight else n.minobsinnode
   
   # throw some errors if the user specifies two versions of the same option
   if (!missing(n.trees) & ('nrounds' %in% args_named))             stop("Only one of n.trees and nrounds can be specified.")
   if (!missing(interaction.depth) & ('max_depth' %in% args_named)) stop("Only one of interaction.depth and max_depth can be specified.")
   if (!missing(shrinkage) & ('eta' %in% args_named))               stop("Only one of shrinkage and eta can be specified.")
   if (!missing(bag.fraction) & ('subsample' %in% args_named))      stop("Only one of bag.fraction and subsample can be specified.")
   if (!missing(n.minobsinnode) & ('min_child_weight' %in% args_named))      stop("Only one of n.minobsinnode and min_child_weight can be specified.")
   
   # throw error if user specifies params with other options
   if (!missing(interaction.depth) | ('max_depth' %in% args_named) | !missing(shrinkage) | ('eta' %in% args_named) | !missing(bag.fraction) | ('subsample' %in% args_named) ){
      if (!is.null(params)) stop("params cannot be specified with any of interaction.depth, max_depth, shrinkage, eta, bag.fraction, or subsample.")
   }
   
   # check that data is not a tibble or data.table
   if( class(data)[1] %in% c("tbl_df","tbl","data.table")){
      stop("The twang package currently does not support data.table or tibble. Please convert your data object to a data.frame.")
   }else{
      if (class(data)[1] != "data.frame"){
         warning("Data classes other than data.frame may cause errors." , call.=FALSE)
      }
   }

   if (version=="legacy"){
      ## throw some errors if the user specifies an option not allowed in legacy version of ps
      if (!is.null(ks.exact))     stop("Option ks.exact is not allowed with version='legacy'")
      if (!is.null(params))       stop("Option params is not allowed with version='legacy'")
      if (!is.null(args$tree_method))  stop("Option tree_method is not allowed with version='legacy'")
      if (!missing(n.keep))       stop("Option n.keep is not allowed with version='legacy'")
      if (!missing(n.grid))       stop("Option n.grid is not allowed with version='legacy'")
      if (!missing(n.minobsinnode) | !is.null(args$min_child_weight)) stop("Options n.minobsinnode or min_child_weight are not allowed with version='legacy'")
      if (!missing(keep.data)) stop("Option keep.data is not allowed with version='legacy'")
      
      return(ps.old(formula = formula,
                     data = data,                         # data
                     n.trees = nrounds,                   # gbm options
                     interaction.depth = max_depth,
                     shrinkage=eta,
                     bag.fraction = subsample,
                     perm.test.iters = perm.test.iters,
                     print.level = print.level,       
                     verbose = verbose,
                     estimand = estimand, 
                     stop.method = stop.method, 
                     sampw = sampw, 
                     multinom = multinom, 
                     ...))
  }else{
    # xgboost tree method  
    tree_method  <- if (!is.null(args$tree_method)) args$tree_method else "hist"
    # throw error if user specifies params with version=="gbm"
    if ( version=="gbm" & !is.null(params) ) stop("params cannot be specified when version='gbm'.")
    return(ps.fast(formula = formula,
                  data = data,                         # data
                  n.trees = nrounds,                   # gbm options
                  interaction.depth = max_depth,
                  shrinkage = eta,
                  bag.fraction = subsample,
                  n.minobsinnode = min_child_weight,
                  params = params,
                  perm.test.iters = perm.test.iters,
                  print.level = print.level,       
                  verbose = verbose,
                  estimand = estimand, 
                  stop.method = stop.method, 
                  sampw = sampw, 
                  multinom = multinom,
                  ks.exact = ks.exact,
                  version = version,
                  tree_method = tree_method,
                  n.keep = n.keep,
                  n.grid = n.grid,
                  keep.data=keep.data))
  }
}
  

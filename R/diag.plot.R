diag.plot <- function(x,plots, subset, ...)
{
	treat <- x$treat
	propScores <- x$ps
	weights <- x$w
   if (!all(weights[,1]==1)){
      weights   <- cbind(unw=rep(1,nrow(weights)),weights)
      propScores <- cbind(unw=rep(0.5,nrow(propScores)),propScores)
   }
   
   dots <- list(...)
   
   n.tp <- ifelse(class(x) == "dxwts", length(x$desc), ncol(weights))
   hldEffSz <- NULL
   tpHld <- NULL
   hldPVal <- hldksPVal <- hldksPValRanks <- hldTVal <- hldTRank <- hldESBig <- NULL
   desc.unw <- x$desc[[1]]
   for (i in 2:n.tp){
   	desc.temp <- x$desc[[i]]
   	iter <- desc.temp$n.trees
   	tpHld <- c(tpHld, names(x$desc[i]))
   	
   	hldEffSz <- c(hldEffSz, desc.temp$bal.tab$results$std.eff.sz, desc.unw$bal.tab$results$std.eff.sz)
   	hldPVal <- c(hldPVal, desc.temp$bal.tab$results$p, desc.unw$bal.tab$results$p)
    hldksPVal <- c(hldksPVal, desc.temp$bal.tab$results$ks.pval, desc.unw$bal.tab$results$ks.pval)
    hldksPValRanks <- c(hldksPValRanks, rank(desc.temp$bal.tab$results$ks.pval, ties.method = "first"), rank(desc.unw$bal.tab$results$ks.pval, ties.method = "first"))
    hldTVal <- c(hldTVal, desc.temp$bal.tab$results$p, desc.unw$bal.tab$results$p)
    hldTRank <- c(hldTRank, rank(desc.temp$bal.tab$results$p, ties.method = "first"), rank(desc.unw$bal.tab$results$p, ties.method = "first"))
    hldESBig <- c(hldESBig, rep(abs(desc.temp$bal.tab$results$std.eff.sz) > abs(desc.unw$bal.tab$results$std.eff.sz),2))

   	
   }
   
   n.var <- length(hldEffSz)/(2*(n.tp - 1))
   
   	## form data.frame w/: p-val, es, weight/not, whichVar,
   	esDat <- data.frame(effectSize = abs(hldEffSz), 
   	pVal=0+1*(hldPVal > 0.05), 
   	whichComp = as.factor(rep(tpHld, each = 2*n.var)), 
   	weighted = as.factor(rep(rep(c("Weighted", "Unweighted"), each = n.var), (n.tp-1))), 
   	whichVar = factor(rep(1:n.var, 2*(n.tp-1))), 
   	ksPVal = hldksPVal, 
   	ksRank = hldksPValRanks, 
   	tPVal = hldTVal,
   	tRank = hldTRank,
   	esBig = hldESBig)   
   	
#   	return(esDat)


#sbSt <- esDat$stopRule %in% levels(esDat$stopRule)[1:length(levels(esDat$stopRule))]

if (plots == "optimize" || plots == 1) {
	longBal <- matrix(t(x$balance))
	optDat <- data.frame(balance = longBal, iteration = rep(x$iters, each = n.tp-1), stopRule = names(x$desc)[-1])

	if(is.null(subset))
	subset <- 1:length(levels(as.factor(optDat$stopRule)))

	pt1 <- xyplot(balance ~ iteration | stopRule, data = optDat, ylab = "Balance measure", xlab = "Iteration", scales = list(alternating = 1), subset = as.factor(optDat$stopRule) %in% levels(as.factor(optDat$stopRule))[subset], ...)	
	
	
}
   
   if (plots == "es" || plots == 3)	{ ## es plot
   	
   	if(is.null(subset))
   	subset <- 1:length(levels(as.factor(esDat$whichComp)))
   	
   	yMax <- min(3,max(esDat[,1])) + .05	
   	
   	if(max(esDat[,1]) > 3)
   	warning("Some effect sizes are larger than 3 and may not have been plotted.\n")	
   	
   	pt1.1 <- xyplot(effectSize ~ weighted | whichComp, groups = whichVar, data = esDat, scales = list(alternating = 1),
   	ylim = c(-.05, yMax), type = "l", col = "lightblue", 
   	subset = !esBig & (as.factor(esDat$whichComp) %in% levels(as.factor(esDat$whichComp))[subset]), 
   	ylab = "Absolute standard difference", xlab = NULL, ...,
   	panel = function(...){
   		panel.abline(h=c(.2,.5,.8), col="gray80")
   		panel.xyplot(...)
   		
   	})
   	pt1.2 <- xyplot(effectSize ~ weighted | whichComp, groups = whichVar, 
   	data = esDat, ylab = "Absolute standard difference", xlab = NULL,
   	ylim = c(-.05, yMax), type = "l", col = "red", 
   	subset = esBig & (as.factor(esDat$whichComp) %in% levels(as.factor(esDat$whichComp))[subset]), lwd = 2)
   	
   	pt2 <- xyplot(effectSize ~ weighted | whichComp, groups = pVal, data = esDat,
   	ylab = "Absolute standard difference", xlab = NULL, 
   	ylim = c(-.05, yMax), type = "p", col = "red", pch = c(19,1),
   	subset = as.factor(esDat$whichComp) %in% levels(as.factor(esDat$whichComp))[subset])
   	
   	pt1 <- pt1.1 + pt1.2 + pt2
   	
   	}
   						
 			
   
   if (plots == "t" || plots == 4) { ## t p-values plot
   	
   	   	if(is.null(subset))
   	subset <- 1:length(levels(as.factor(esDat$whichComp)))
   	
   	n.var2 <- max(esDat$tRank * (!is.na(esDat$tPVal)), na.rm=TRUE)
   	
   	pt1 <- xyplot(tPVal~tRank|whichComp, groups = weighted, data=esDat, xlab = "Rank of p-value rank for pretreatment variables \n (hollow is weighted, solid is unweighted)", ylab = "T test p-values", pch = c(19,1), col = "black", scales = list(alternating = 1),
   	subset = (as.factor(esDat$whichComp) %in% levels(as.factor(esDat$whichComp))[subset]) & (esDat$tRank <= n.var2), ..., 
   	   	panel = function(...){
   	   		panel.xyplot(x=c(1,n.var2), y=c(0,1), col="lightblue", type="l")
#   		panel.abline(a=0, b=1, col="lightblue")
   		panel.xyplot(...)
   		}
)
   	
   	}
   
   if (plots =="ks" || plots ==5) {  ## ks plot
   	
   	if(is.null(subset))
   	subset <- 1:length(levels(as.factor(esDat$whichComp)))
   	
   	n.var2 <- max(esDat$ksRank*(!is.na(esDat$ksPVal)), na.rm=TRUE)
   	
   	pt1 <- xyplot(ksPVal~ksRank|whichComp, groups=weighted, scales = list(alternating = 1), data = esDat, ..., xlab = "Rank of p-value rank for pretreatment variables \n (hollow is weighted, solid is unweighted)", ylab = "KS test p-values", pch = c(19,1), col="black",
   	subset = (as.factor(esDat$whichComp) %in% levels(as.factor(esDat$whichComp))[subset]) & (esDat$ksRank <= n.var2),
   	panel = function(...){
   		panel.xyplot(x=c(1,n.var2), y=c(0,1), col="lightblue", type="l")
   		panel.xyplot(...)
   	})
   	}
   	
#   	if (plots == "histogram" || plots == 10)
#   	pt1 <- histogram.dxwts(...)
   	
   	if (plots == "boxplot" || plots == 2)
   	pt1 <- boxplot(x, ...)
   	
   	if (plots == "histogram" || plots == 6){
   		
   			treat <- x$treat
	wghts <- x$w
	vars <- x$variableNames
	if (is.null(dots$main)) dots$main <- "Control weights"
#	if(is.null(x$w.ctrl) || is.null(x$treat))
#         stop("For the weight histogram w.ctrl and treat cannot be NULL")	
	if(!all(wghts[,1]==1))   {
		controlWeights <- wghts[x$treat == 0, ]
	}
	else controlWeights <- wghts[x$treat == 0, -1]
	
	controlWeights <- as.matrix(controlWeights)
	nFrame <- ncol(controlWeights)
	longWeights <- matrix(t(controlWeights), ncol=1)
	longWeights <- data.frame(Weights = longWeights, varb = colnames(controlWeights))
	
	if(is.null(subset))
	subset <- 1:length(colnames(controlWeights))
	
	
	pt1 <- histogram(~Weights|varb, data=longWeights, subset = varb %in% colnames(controlWeights)[subset],...)

   		
   		
   	}
   	
   	
   	if(!(plots %in% c(1:6, "boxplot","t","ks","optimize","es", "histogram")))
   	stop("plots must be an integer from 1 to 5, or be one of 'optimize' \n
   	'boxplot','es','t', 'ks', or 'histogram'\n")
   	
   
   return(pt1)


}

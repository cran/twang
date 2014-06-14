plot.mnps.old <- function(x,plots="optimize", figureRows = 1, summaryFcn = max, color = TRUE, pairwise = TRUE, ...)
{
   # Creates diag.plot plots and sends to current device
   # x:     ps object 
   # label: Label added to the plot titles

   # extract the propensity scores and weights from the ps object
   
   if(x$estimand == "ATT") pairwise <- TRUE
   
   	ltBl <- ifelse(color, "lightblue","gray80")
	rdCol <- ifelse(color, "red","black")
	stripBgCol <- ifelse(color, "#ffe5cc", "transparent")
	ptSymCol <- ifelse(color, "#0080ff", "black")

   
   if(!is.null(summaryFcn)) hldFcn <- match.fun(summaryFcn)
   
   ifelse(x$estimand == "ATE", noKS <- TRUE, noKS <- FALSE)
   
   subst <- whichVar <- pVal <- weighted <- NULL 
   
   if(length(plots) > 1) stop("The `plots' argument must be of length 1.")
   if(!(plots %in% c(1,2,3,4,5)) & !(plots %in% c("boxplot","ks","optimize","es","t")))
   stop("Invalid choice of `plots' argument.")
   
   if(plots == 2 | plots == "boxplot"){
   	
   	boxplot(x, color = color, par.settings = list(strip.background = list(col=stripBgCol)),...)
   	
   }
   
#   else if((plots == 5 | plots == "ks") & x$estimand == "ATE"){
#   	
#   	warning("KS p-values are not available for mnps objects with estimand = 'ATE'")
#   }
   
   else{
   
   if(is.null(summaryFcn) | plots == "optimize" | plots == 1){
   
   nPlot <- x$nFits
   ptHld <- vector(mode = "list", length = nPlot)
   for(i in 1:nPlot){
   	if(x$estimand == "ATT") ptNm <- paste("Balance for", x$levExceptTreatATT[i], "versus unweighted", x$treatATT)
   	else ptNm <- paste("Balance for", x$treatLev[i], "against others")
   	ptHld[[i]] <- plot(x$psList[[i]], main = ptNm, plots = plots, noKS = TRUE, color = color, ...)
   }

#pt1 <- twang:::diag.plot(x, plots, subset = subset, ...)


figCol <- ceiling(nPlot/figureRows)

if(dev.cur() == 1) dev.new()

curCol <- curRow <- 1

for(i in 1:(nPlot-1)){
	print(ptHld[[i]], split = c(curCol,curRow,nx = figCol,ny = figureRows), more = TRUE)
	if(curCol < figCol){
		curCol <- curCol + 1
	}
	else {
		curCol <- 1
		curRow <- curRow + 1
	}
}

print(ptHld[[nPlot]], split = c(curCol,curRow,nx = figCol,ny = figureRows), more = FALSE)

}

else{  ## if summaryFcn isn't null and plots == something other than optimize

n.tp <- length(x$psList[[1]]$desc)
n.psFits <- length(x$psList)
   	
#   	return(esDat)   
   	
   if (plots == "es" || plots == 3)	{ ## es plot
   	
   	esDat <- makePlotDat(x$psList[[1]], whichPlot = 3)
   	nVar <- length(makePlotDat(x$psList[[1]], whichPlot = 3, subsetStopMeth = 1, yOnly = TRUE, incUnw = FALSE)$pVal)
   	effSzList <- effSzListUnw <- pValListUnw <- pValList <- matrix(NA, nrow = nVar, ncol = length(x$psList))
   	
   	for(k in 2:n.tp){
   	for(j in 1:length(x$psList)){	
	x2 <- x$psList[[j]]
	effSzList[,j] <- makePlotDat(x2, whichPlot = 3, subsetStopMeth = k, yOnly = TRUE, incUnw = FALSE)$effectSize
	effSzListUnw[,j] <- makePlotDat(x2, whichPlot = 3, subsetStopMeth = 1, yOnly = TRUE, incUnw = FALSE)$effectSize
	pValListUnw[,j] <- makePlotDat(x2, whichPlot = 3, subsetStopMeth = 1, yOnly = TRUE, incUnw = FALSE)$pVal
	pValList[,j] <- makePlotDat(x2, whichPlot = 3, subsetStopMeth = k, yOnly = TRUE, incUnw = FALSE)$pVal	
	}
	effSzList <- abs(effSzList)
	effSzListUnw <- abs(effSzListUnw)
	collapsed <- apply(effSzList, 1, hldFcn)
	collapsedUnw <- apply(effSzListUnw, 1, hldFcn)
	if(pairwise & x$estimand == "ATE"){
		pwc <- pairwiseCompForPlot(x, summaryFcn = summaryFcn, stop.method = (k-1))
		collapsed <- pwc$weightedPairwiseDiff
		collapsedUnw <- pwc$unweightedPairwiseDiff
	}
	collapsedP <- apply(pValList < .05, 1, max)
	collapsedUnwP <- apply(pValListUnw < .05, 1, max)	
	
	esBigHold <- collapsed > collapsedUnw
	esDat$effectSize[(1:(2*nVar)) + (k-2)*2*nVar] <- c(collapsed, collapsedUnw)
	esDat$pVal[(1:(2*nVar)) + (k-2)*2*nVar] <- c(collapsedP, collapsedUnwP)
	esDat$esBig[(1:(2*nVar)) + (k-2)*2*nVar] <- rep(esBigHold, 2)
	}
   	
   	if(is.null(subst)) subst <- 1:length(levels(as.factor(esDat$whichComp)))
   	
   	yMax <- min(3,max(esDat[,1])) + .05	
   	
   	if(max(esDat[,1], na.rm=TRUE) > 3)
   	warning("Some effect sizes are larger than 3 and may not have been plotted.\n")	
   	
    
    nullPlot <- TRUE
       	
   	subsetHold <- !esDat$esBig & (as.factor(esDat$whichComp) %in% levels(as.factor(esDat$whichComp))[subst])
   	
   	if(any(subsetHold)){
   	pt1.1 <- xyplot(effectSize ~ weighted | whichComp, groups = whichVar, data = esDat, scales = list(alternating = 1),
   	ylim = c(-.05, yMax), type = "l", col = ltBl, subset = subsetHold, 
   	ylab = "Absolute standard difference", xlab = NULL, par.settings = list(strip.background = list(col=stripBgCol)),
   	panel = function(...){
   		panel.abline(h=c(.2,.5,.8), col="gray80")
   		panel.xyplot(...)
   		
   	})
   	nullPlot <- FALSE
   	currPt <- pt1.1
   	}
   	
   	subsetHold <- esDat$esBig & (as.factor(esDat$whichComp) %in% levels(as.factor(esDat$whichComp))[subst])
   	
   	if(any(subsetHold)){
   	pt1.2 <- xyplot(effectSize ~ weighted | whichComp, groups = whichVar, 
   	data = esDat, ylab = "Absolute standard difference", xlab = NULL,
   	ylim = c(-.05, yMax), type = "l", col = rdCol, par.settings = list(strip.background = list(col=stripBgCol)),
   	subset = subsetHold, lwd = 2)
   	if(!nullPlot){
   		currPt <- currPt + pt1.2
   	}
   	else{
   		currPt <- pt1.2
   		nullPlot <- FALSE
   		}
   	}
   	
   	subsetHold <- as.factor(esDat$whichComp) %in% levels(as.factor(esDat$whichComp))[subst]
   	
   	if(all(esDat$pVal == 1)) pchHold <- 19
   	else if(all(esDat$pVal == 0)) pchHold <- 1
   	else pchHold <- c(19,1)
   	
   	if(any(subsetHold)){
   	pt2 <- xyplot(effectSize ~ weighted | whichComp, groups = (pVal < 0.05), data = esDat,
   	ylab = "Absolute standard difference", xlab = NULL, 
   	ylim = c(-.05, yMax), type = "p", col = rdCol, pch = pchHold,par.settings = list(strip.background = list(col=stripBgCol)),
   	subset = subsetHold)
   	if(!nullPlot) currPt <- currPt + pt2
   	else currPt <- pt2
   	}
   	
   	return(currPt)
   	
   	}
   						
 			
   
   if (plots == "t" || plots == 4) { ## t p-values plot

   	esDat <- makePlotDat(x$psList[[1]], whichPlot = 4)
   	nVar <- length(makePlotDat(x$psList[[1]], whichPlot = 4, subsetStopMeth = 1, yOnly = TRUE, incUnw = FALSE))
   	effSzList <- effSzListUnw <- matrix(NA, nrow = nVar, ncol = length(x$psList))   	
   	if(x$estimand =="ATE") pwc <- pairwiseComparison(x, collapse.to = "covariate")
  	
   	
   	for(k in 2:n.tp){
   	for(j in 1:length(x$psList)){	
	x2 <- x$psList[[j]]
	effSzList[,j] <- makePlotDat(x2, whichPlot = 4, subsetStopMeth = k, yOnly = TRUE, incUnw = FALSE)
	effSzListUnw[,j] <- makePlotDat(x2, whichPlot = 4, subsetStopMeth = 1, yOnly = TRUE, incUnw = FALSE)
	collapsed <- apply(-1 *effSzList, 1, hldFcn)
	collapsed <- -1 * collapsed
	if(x$estimand == "ATE") collapsed <- pwc$min.p[pwc$stop.method == x$stopMethods[k-1]]	
	collapsedUnw <- apply(-1 * effSzListUnw, 1, hldFcn)
	collapsedUnw <- -1 * collapsedUnw
	if(x$estimand == "ATE") collapsedUnw <- pwc$min.p[pwc$stop.method == 'unw']		
	collRanks <- rank(collapsed, ties.method = "first")
	esBigHold <- collapsed > collapsedUnw
	esDat$tPVal[(1:(2*nVar)) + (k-2)*2*nVar] <- c(collapsed, collapsedUnw)
	esDat$tRank[(1:(2*nVar)) + (k-2)*2*nVar] <- rep(collRanks, 2)
	}
	}
   	
   	   	if(is.null(subst))	subst <- 1:length(levels(as.factor(esDat$whichComp)))
   	
   	n.var2 <- max(esDat$tRank * (!is.na(esDat$tPVal)), na.rm=TRUE)
   	
   	pt1 <- xyplot(tPVal~tRank|whichComp, groups = weighted, data=esDat, xlab = "Rank of p-value rank for pretreatment variables \n (hollow is weighted, solid is unweighted)", ylab = "T test p-values", pch = c(19,1), col = "black", scales = list(alternating = 1), par.settings = list(strip.background = list(col=stripBgCol)),
   	subst = (as.factor(esDat$whichComp) %in% levels(as.factor(esDat$whichComp))[subst]) & (esDat$tRank <= n.var2), ylim = c(-.1, 1.1), ..., 
   	   	panel = function(...){
   	   		panel.xyplot(x=c(1,n.var2), y=c(0,1), col=ltBl, type="l")
#   		panel.abline(a=0, b=1, col="lightblue")
   		panel.xyplot(...)
   		}
)
   	
   	}
   
   if (plots =="ks" || plots ==5) {  ## ks plot
   	if(x$estimand =="ATE") pwc <- pairwiseComparison(x, collapse.to = "covariate")

   	esDat <- makePlotDat(x$psList[[1]], whichPlot = 5)
   	nVar <- length(makePlotDat(x$psList[[1]], whichPlot = 5, subsetStopMeth = 1, yOnly = TRUE, incUnw = FALSE))
   	effSzList <- effSzListUnw <- matrix(NA, nrow = nVar, ncol = length(x$psList))   	
   	
   	for(k in 2:n.tp){
   	for(j in 1:length(x$psList)){	
	x2 <- x$psList[[j]]
	effSzList[,j] <- makePlotDat(x2, whichPlot = 5, subsetStopMeth = k, yOnly = TRUE, incUnw = FALSE)
	effSzListUnw[,j] <- makePlotDat(x2, whichPlot = 5, subsetStopMeth = 1, yOnly = TRUE, incUnw = FALSE)
	collapsed <- apply(effSzList, 1, hldFcn)
	if(x$estimand == "ATE") collapsed <- pwc$min.ks.pval[pwc$stop.method == x$stopMethods[k-1]]
	collapsedUnw <- apply(effSzListUnw, 1, hldFcn)
	if(x$estimand == "ATE") collapsedUnw <- pwc$min.ks.pval[pwc$stop.method == "unw"]	
	collRanks <- rank(collapsed, ties.method = "first")
	esBigHold <- collapsed > collapsedUnw
	esDat$ksPVal[(1:(2*nVar)) + (k-2)*2*nVar] <- c(collapsed, collapsedUnw)
	esDat$ksRank[(1:(2*nVar)) + (k-2)*2*nVar] <- rep(collRanks, 2)
	}
	}


   	
   	if(is.null(subst))	subst <- 1:length(levels(as.factor(esDat$whichComp)))
   	
   	n.var2 <- max(esDat$ksRank*(!is.na(esDat$ksPVal)), na.rm=TRUE)
   	
   	pt1 <- xyplot(ksPVal~ksRank|whichComp, groups=weighted, scales = list(alternating = 1), data = esDat,ylim = c(-.1, 1.1), ..., xlab = "Rank of p-value rank for pretreatment variables \n (hollow is weighted, solid is unweighted)", ylab = "KS test p-values", pch = c(19,1), col="black",par.settings = list(strip.background = list(col=stripBgCol)),
   	subst = (as.factor(esDat$whichComp) %in% levels(as.factor(esDat$whichComp))[subst]) & (esDat$ksRank <= n.var2),
   	panel = function(...){
   		panel.xyplot(x=c(1,n.var2), y=c(0,1), col=ltBl, type="l")
   		panel.xyplot(...)
   	})
   	}
   	
return(pt1)
	
	}

}

}

# Produces a summary table for mnps object 
summary.mnps <- function(object,...){
	nFits <- object$nFits
	summaryList <- vector(mode = "list", length = nFits)
	for(i in 1:nFits){
		summaryList[[i]] <- summary(object$psList[[i]], ...)
	}

	if(object$estimand == "ATT"){
	
	retObj <- list(summaryList = summaryList, nFit = object$nFit, estimand = object$estimand, treatATT = object$treatATT, treatLev = object$treatLev, levExceptTreatATT = object$levExceptTreatATT)
	
	class(retObj) <- "summary.mnps"
	return(retObj)
	
	}
	else{
		hd1 <- pairwiseComparison(object, collapse.to = "stop.method")
		ns <- NULL
		for(i in 1:length(summaryList)) ns <- c(ns, summaryList[[i]][1,1])
		shell <- data.frame(treatment = object$treatLev, n = ns)
		for(i in 1:length(object$stopMethods)){
			hldESS <- NULL
			for(j in 1:length(object$treatLev)){
				hldESS <- c(hldESS, summaryList[[j]][i+1,3])

			}
			shell <- data.frame(shell, currESS = hldESS)
			names(shell)[names(shell) == "currESS"] <- paste("ESS:", object$stopMethods[i], sep = "")
		}
		
		retObj <- list(comp = hd1, ess = shell, estimand = object$estimand)
		class(retObj) <- "summary.mnps"
		return(retObj)
		
		}
	
}


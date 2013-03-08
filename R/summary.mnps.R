# Produces a summary table for ps object 
summary.mnps <- function(object,...)
{
	nFits <- object$nFits
	summaryList <- vector(mode = "list", length = nFits)
	for(i in 1:nFits){
		summaryList[[i]] <- summary(object$psList[[i]], ...)
	}
	
	retObj <- list(summaryList = summaryList, nFit = object$nFit, estimand = object$estimand, treatATT = object$treatATT, treatLev = object$treatLev, levExceptTreatATT = object$levExceptTreatATT)
	
	class(retObj) <- "summary.mnps"
	return(retObj)
	
}


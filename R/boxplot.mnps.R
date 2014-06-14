boxplot.mnps <- function(x, stop.method = NULL, color = TRUE, ...){
	
	ptSymCol <- ifelse(color, "#0080ff", "black")	
	bwCols <- list(col = ptSymCol)
	stripBgCol <- ifelse(color, "#ffe5cc", "transparent")
	
	if(is.null(stop.method)) stop.method <- x$stopMethods
	
	if(length(stop.method) > 1){ 
		if(is.numeric(stop.method))
		warning("Using only the first stop.method, ",  x$stopMethods[1])
		else warning("Using only the first stop.method, ",  stop.method[1])
		stop.method <- stop.method[1]
		}
			
	
	
	if(is.numeric(stop.method)) stop.method = x$stopMethods[stop.method]
	

	
	if(x$estimand == "ATE"){

			bwDat <- whichResp <- NULL
			
			stopMethLong <- paste(stop.method, ".ATE", sep = "")
			
		for(j in 1:x$nFits){
#			bwDat <- rbind(bwDat, data.frame(ps = x$psList[[j]]$ps[,currStopMeth], treat = x$treatVar, whichResp = x$respLev[j]))
			bwDat <- data.frame(ps = x$psList[[j]]$ps[,stopMethLong], treat = x$treatVar, whichResp = x$treatLev[j])

		
		
		
		pt1 <- bwplot(ps ~ treat, groups = whichResp, #layout = c(1,x$nFits), 
		xlab = "Treatment", ylab = "Propensity scores", ylim = c(-.1,1.1), data = bwDat, main = paste(x$treatLev[j], " propensity scores by Tx group"),par.settings = list(strip.background = list(col=stripBgCol), box.rectangle = bwCols, plot.symbol = bwCols, box.umbrella = bwCols), ...)


		print(pt1, split = c(1,j,1,x$nFits), more = (j < x$nFits))

		}
		}
		
		
	else if(x$estimand == "ATT"){
		
		bwDat <- NULL
		
		stopMethLong <- paste(stop.method, ".ATT", sep = "")
		bwDat <- NULL
		
		for(j in 1:x$nFits){
			currCats <- c(x$treatATT, x$levExceptTreatATT[j])
			bwDat <- data.frame(ps = x$psList[[j]]$ps[,stopMethLong], treat = currCats[1 + x$psList[[j]]$data$currResp], respCat = x$levExceptTreatATT[j], attGrp = x$treatATT)
			pt1 <- bwplot(ps ~ treat, data = bwDat, ylim = c(-.1,1.1), ylab = "Propensity scores", xlab = "Treatment", main = paste("Propensity score of ", x$levExceptTreatATT[j], " versus ", x$treatATT, ".", sep = ""),par.settings = list(strip.background = list(col=stripBgCol), box.rectangle = bwCols, plot.symbol = bwCols, box.umbrella = bwCols), ...)
			print(pt1, split = c(1,j,1,x$nFits), more = (j < x$nFits))
			
		}		
	}	


	 

}
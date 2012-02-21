boxplot.ps <- function(x, ...){
	longDat <- matrix(t(as.matrix(x$ps)), ncol = 1)
	nms <- names(x$ps)
	bwDat <- data.frame(ps = longDat, nm = nms, treat = rep(x$treat, each = length(nms)))
	pt1 <- bwplot(treat~ps|nm, data=bwDat, scales = list(alternating = 1), ylab = "Treatment", xlab="Propensity scores")
	
}
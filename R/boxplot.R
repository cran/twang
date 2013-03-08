boxplot <- function(x, stop.method = NULL, ...){
	if(class(x) == "mnps") boxplot.mnps(x=x, stop.method = NULL, ...)
	if(class(x) == "ps") boxplot.ps(x=x, ...)
}
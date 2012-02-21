plot.ps <- function(x,plots="optimize",subset=NULL,...)
{
   # Creates diag.plot plots and sends to current device
   # x:     ps object 
   # label: Label added to the plot titles

   # extract the propensity scores and weights from the ps object

pt1 <- diag.plot(x, plots, subset = subset, ...)

return(pt1)

}

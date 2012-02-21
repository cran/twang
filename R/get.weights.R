get.weights <- function(ps1, stop.method, estimand = NULL)
{
   if(class(ps1)!="ps") stop("ps1 must be a ps object.")
   if(is.null(estimand)) estimand <- ps1$estimand
   
   if(!(estimand %in% c("ATT","ATE"))) stop("estimand must be either \"ATT\" or \"ATE\".")
   if(estimand != ps1$estimand){
   	warning("Estimand specified for get.weights() differs from the estimand used to fit the ps object.")
   }
   if(length(stop.method)>1) stop("More than one stop.method was selected.")
   if(!is.null(stop.method))
   {
   	stop.method.long <- paste(stop.method, ps1$estimand, sep=".")
      i <- match(stop.method.long, names(ps1$w))
      if(is.na(i)) stop("Weights for stop.method=",stop.method, " and estimand=", estimand, " are not available. Please a stop.method and used when fitting the ps object.") 
#      Available options: ",names(ps1$ps),".")
   } else
   {
   	warning("No stop.method specified.  Using ")
      i <- 1
   }
  
   if (estimand == "ATT") return(ps1$w[[i]])
   else if (estimand == "ATE")
   { 
      w <- with(ps1, treat/ps[[i]] + (1-treat)/(1-ps[[i]]))
      return(w)
   }
}

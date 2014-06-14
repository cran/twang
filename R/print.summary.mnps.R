# Produces a summary table for ps object 
print.summary.mnps <- function(x, ...)
{
      dots <- list(...)


      if(x$estimand == "ATE"){
      	cat("Summary of pairwise comparisons:\n")
      	print(x$comp)
      	cat("\nSample sizes and effective sample sizes:\n")
      	print(x$ess)
      		
      }
      if(x$estimand == "ATT"){
      	cat("Summary of mnps object:\n")
      	nSum <- length(x$summaryList)      		
      	for(i in 1:nSum){
      		cat("Summary of observations receiving treatment ", x$levExceptTreatATT[i], " weighted to match the observations receiving treatment ", x$treatATT, ".\n", sep = "")
      	   	if(!is.null(dots$digits)) obj <- round(x$summaryList[[i]][,-c(1,2,3,8)], digits = digits)
      		else obj <- x$summaryList[[i]][,-c(1,2,3,8)]
		    obj <- data.frame(obj)	
#      class(obj) <- "matrix"
	      	names(obj) <- c("ESS", "max.es","mean.es","max.ks","mean.ks","iter")	
      	
      	}

      print(obj)
      cat("\n")
      }
      invisible(x)
      }


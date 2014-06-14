## Generic function for extracting balance tables from ps and other objects
bal.table <- function(x, digits = 3, collapse.to = c("pair","covariate","stop.method")[1]){
	if(class(x) != "mnps"){
   bal.tab <- lapply(x$desc, function(x){return(round(x$bal.tab$results, digits))})
   return(bal.tab)
   }
   else {
   	if(x$estimand == "ATE"){
   		pwc <- pairwiseComparison(x, collapse.to = collapse.to)
   		hldNum <- NULL
   		for(i in 1:length(pwc[1,])) hldNum <- c(hldNum, is.numeric(pwc[1,i]))
   		pwc[,hldNum] <- round(pwc[,hldNum], digits = digits)
   		return(pwc)
   		}
   	nFits <- x$nFits
   	balTabList <- vector(mode = "list", length = nFits)
   	#if(x$estimand == "ATT")
   	cat(paste("Note that `tx` refers to the category specified as the treatATT, ", x$treatATT, ".\n\n", sep = ""))
   	for(i in 1:nFits) balTabList[[i]] <- bal.table(x$psList[[i]], digits = digits)
   	#if(x$estimand == "ATT") 
   	names(balTabList) <- x$levExceptTreatATT
   	for(i in 1:length(balTabList)){
   		for(j in 1:length(balTabList[[i]])){
   			balTabList[[i]][[j]] <- data.frame(var = row.names(balTabList[[i]][[j]]), balTabList[[i]][[j]], control = names(balTabList)[i], stop.method = names(balTabList[[i]])[j])
   		}
#   		balTabList[[i]] <- do.call(rbind, balTabList[[i]])
   	}
#   	balTabList <- do.call(rbind, balTabList)
   	
   	nonTreatATT <- x$levExceptTreatATT
   	newBalTabList <- vector(mode = "list", length = (length(x$stopMethods) + 1) * length(nonTreatATT))
   	cnt <- 1
   	for(i in 1:length(nonTreatATT)){
   		newBalTabList[[cnt]] <- balTabList[[i]][["unw"]]
   		cnt <- cnt + 1
   	}
   	for(i in 1:length(x$stopMethods)){
   		for(j in 1:length(nonTreatATT)){
   			newBalTabList[[cnt]] <- balTabList[[nonTreatATT[j]]][[paste(x$stopMethods[i], ".ATT", sep = "")]]
   			cnt <- cnt + 1
   		}
   	}
   	balTabList <- do.call(rbind, newBalTabList)
   	row.names(balTabList) <- NULL
   	#if(x$estimand == "ATE") names(balTabList) <- x$treatLev
#   	if(pairwise){
#   		allDiffs <- NULL
#   		for(i in 1:length(balTabList[[1]])){
#   			holdDiffs <- rep(0, nrow(balTabList[[1]][[1]]))
#   			for(j in 2:nFits){
#   				for(k in 1:(j-1)){
#   					holdDiffs <- apply(cbind(holdDiffs, abs(balTabList[[j]][[i]]$tx.mn - balTabList[[k]][[i]]$tx.mn)), 1, max, na.rm=TRUE)
#   				}
#   			}
#   			allDiffs <- cbind(allDiffs, holdDiffs/balTabList[[1]][[1]]$ct.sd)
#   		}
#   		row.names(allDiffs) <- row.names(balTabList[[1]][[1]])
#   		dimnames(allDiffs)[[2]] <- names(x$psList[[1]]$desc)
##   		return(list(pairwiseDiff = allDiffs, balanceTable = balTabList))
#   		return(balTabList)
#   		
#   	}
   	#else 
   	return(balTabList)
   }
   	
   	
}



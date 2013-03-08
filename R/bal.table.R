## Generic function for extracting balance tables from ps and other objects
bal.table <- function(x){
	if(class(x) != "mnps"){
   bal.tab <- lapply(x$desc, function(x){return(round(x$bal.tab$results,3))})
   return(bal.tab)
   }
   else {
   	nFits <- x$nFits
   	balTabList <- vector(mode = "list", length = nFits)
   	for(i in 1:nFits) balTabList[[i]] <- bal.table(x$psList[[i]])
   	if(x$estimand == "ATT") names(balTabList) <- x$levExceptTreatATT
   	if(x$estimand == "ATE") names(balTabList) <- x$treatLev
   	return(balTabList)
   }
   	
   	
}



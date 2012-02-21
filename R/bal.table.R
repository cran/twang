## Generic function for extracting balance tables from ps and other objects
bal.table <- function(x){
   bal.tab <- lapply(x$desc, function(x){return(round(x$bal.tab$results,3))})
   return(bal.tab)
}



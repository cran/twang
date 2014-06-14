makeWeightsMNPS <- function(ps, estimand, sampW = NULL, treatATT = NULL, treatATTps = NULL){
	### ps should be a vector of the probability of having received the treatment that was received
	if(estimand == "ATE"){
		wts <- 1/ps
		return(wts)
	}
	if(estimand == "ATT"){
		wts <- sampW
		wts[treatATT] <- sampW[treatATT]
		wts[!treatATT] <- treatATTps[!treatATT]/ps[!treatATT]
		return(wts)
	}
}
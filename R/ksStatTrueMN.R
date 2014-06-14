ksStatTrueMN <- function(gbm1, i, data, sampw, rule.summary = NULL, na.action = "level", vars, treat.var, estimand, treatATT = NULL, collapse.by.var = TRUE){
	if(!(estimand %in% c("ATT","ATE"))) stop("estimand must be either \"ATT\" or \"ATE\".")

	pScores <- predict(gbm1, type = "response")[,,1]
	pMat <- pScores
	treatMat <- NULL
	treatLev <- levels(data[,treat.var])
	nTreat <- length(treatLev)
	for(i in 1:ncol(pScores)) treatMat <- cbind(treatMat, data[,treat.var] == treatLev[i])
	
	if(estimand == "ATE"){
		pMat <- 1/pMat
		pMat <- pMat * treatMat		
		wtVec <- rowSums(pMat)
	}
	
	if(estimand == "ATT"){
		tATT <- which.max(treatLev == treatATT)
		for(i in 1:ncol(pMat)) pMat[,i] <- pMat[,tATT]/pMat[,i]
		pMat[,tATT] <- 1/sum(treatMat[,tATT])
		pMat <- pMat * treatMat
		wtVec <- rowSums(pMat)
	}
	
	nTreat <- ncol(pMat)
	
	nComp <- ifelse(estimand == "ATE", choose(nTreat,2), nTreat-1)
	
	hldKS <- NULL
	
	if(estimand == "ATT"){	
		for(i in 1:(nComp+1)){
			if(i != tATT){
				subDat <- subset(data, data[,treat.var] %in% treatLev[c(i,tATT)])
				subW <- subset(wtVec, data[,treat.var] %in% treatLev[c(i,tATT)])
				subSampW <- subset(sampw, data[,treat.var] %in% treatLev[c(i,tATT)])
				subTreat <- subset(data[,treat.var], data[,treat.var] %in% treatLev[c(i,tATT)])
				subTreat <- subTreat == treatLev[tATT]
				hldKS <- cbind(hldKS, bal.stat(data, vars = vars, treat.var=treat.var, w.all = wtVec, 
					sampw = sampw, get.means = FALSE, na.action = na.action, estimand = estimand, 
					multinom = TRUE)$ks)		
			}
		}
	}
	else{  ## ie, estimand == "ATE"
		for(i in 2:nTreat){
			for(j in 1:(i-1)){
				subDat <- subset(data, data[,treat.var] %in% treatLev[c(i,j)])
				subW <- subset(wtVec, data[,treat.var] %in% treatLev[c(i,j)])
				subSampW <- subset(sampw, data[,treat.var] %in% treatLev[c(i,j)])
				subTreat <- subset(data[,treat.var], data[,treat.var] %in% treatLev[c(i,j)])
				subTreat <- subTreat == treatLev[i]
				hldKS <- cbind(hldKS, bal.stat(data, vars = vars, treat.var=treat.var, w.all = wtVec, 
					sampw = sampw, get.means = FALSE, na.action = na.action, estimand = estimand, 
					multinom = TRUE)$ks)						
			}
		}
	}
	
	
	
	ks <- NULL
	vr <- ifelse(estimand == "ATE", vars, vars[vars != treatATT])
	for(i in 1:length(vr)){
		subDt <- data[data[,treat.var] %in% c(vr[i], treatATT),]
		subW <- wtVec[data[,treat.var] %in% c(vr[i], treatATT),]
		subSampW <- sampw[data[,treat.var] %in% c(vr[i], treatATT),]
		ks <- ps.summary.new2(subDt[,vr], t = as.numeric(subDt[,treat.var]),w=subW, 
			sampw = subSampW, get.means=FALSE, get.ks = TRUE, na.action = na.action,
			collapse.by.var = collapse.by.var, estimand = estimand, multinom = TRUE)
		
	}

	
	
	if(!is.null(rule.summary)) rule.summary <- match.fun(rule.summary)
	
	

	
	
}
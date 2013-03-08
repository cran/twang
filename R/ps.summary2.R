ps.summary <- function(x, t, w,
                       get.means=TRUE, get.ks=TRUE,
                       na.action=c("level","exclude","lowest")[1],
                       collapse.by.var=FALSE,
                       estimand, multinom = FALSE){
   if(!is.factor(x))
      return(ps.summary.n(x=x, t=t, w=w,
                         get.means=get.means, get.ks=get.ks,
                         na.action=na.action,
                         collapse.by.var=collapse.by.var, estimand=estimand, multinom = multinom))
   else
      return(ps.summary.f(x=x, t=t, w=w,
                         get.means=get.means, get.ks=get.ks,
                         na.action=na.action,
                         collapse.by.var=collapse.by.var, estimand=estimand, multinom = multinom))
}


# takes one numeric variable, treatment, and weights and computes balance stats
ps.summary.n <- function(x, t, w,
                         get.means=TRUE, get.ks=TRUE,
                         na.action=c("level","exclude","lowest")[1],
                         collapse.by.var=FALSE, estimand, multinom){

#multinom <- FALSE

if(!(estimand %in% c("ATT","ATE"))) stop("estimand must be either \"ATT\" or \"ATE\".")

if(estimand=="ATT") {

   design <- svydesign(ids=~1, weights=~w, 
                       data=data.frame(x=x,t=t,w=w,miss=is.na(x)))
   if(na.action=="exclude") design <- subset(design, !is.na(x))
   ret <- NULL

   if(get.means)
   {
      design.t <- subset(design, t==1)
      design.c <- subset(design, t==0)
      
      m.t <- svymean(~x, design.t, na.rm=TRUE)
      m.c <- svymean(~x, design.c, na.rm=TRUE)

      sd.t <- sqrt(svyvar(~x, design.t, na.rm=TRUE))
      sd.c <- sqrt(svyvar(~x, design.c, na.rm=TRUE))

      t.n <- summary(svyglm(x~t,design))$coefficients[2,3:4]

      b.n <- ifelse(sd.t==0.0, NA, (m.t-m.c)/sd.t)
      ret <- cbind(m.t, sd.t, m.c, sd.c, b.n, t.n[1], t.n[2])
      colnames(ret) <- c("tx.mn","tx.sd","ct.mn","ct.sd",
                         "std.eff.sz","stat","p")
      if((sum(is.na(x))>0) && (na.action=="level"))
      {
         m.t <- svymean(~is.na(x), design.t, na.rm=TRUE)[2]
         m.c <- svymean(~is.na(x), design.c, na.rm=TRUE)[2]

         sd.t <- sqrt(m.t*(1-m.t))
         sd.c <- sqrt(m.c*(1-m.c))

         test <- try(summary(svyglm(is.na(x)~t,family=quasibinomial,design)), silent=TRUE)
         if(class(test)[1] != "try-error") t.n <- test$coefficients[2,3:4]
         else t.n <- c(NA,NA)

         b.n <- ifelse(sd.t==0.0, NA, (m.t-m.c)/sd.t)
         ret <- rbind(ret, c(m.t, sd.t, m.c, sd.c, b.n, t.n))
      }
   }

   if(get.ks){
      work <- design$variables
      if(na.action=="lowest") work$x[is.na(work$x)] <- min(work$x,na.rm=TRUE)-1
      if(na.action=="level")  work <- subset(work, !is.na(x))

      work$w[work$t==1] <- with(subset(work,t==1),  w/sum(w))
      work$w[work$t==0] <- with(subset(work,t==0), -w/sum(w))
      
      ess <- with(work, sapply(split(w,t), function(w){sum(w)^2/sum(w^2)}))

      ind  <- order(work$x)
      cumv <- abs(cumsum(work$w[ind]))
      cumv <- cumv[diff(work$x[ind])!=0]
      ks   <- ifelse(length(cumv)>0, max(cumv), 0)
      
      # XXX jitter?
#      pval <- 1 - .C("psmirnov2x", p = as.double(ks), 
#              as.integer(ess[2]), as.integer(ess[1]), PACKAGE = "stats")$p
      pval <- 1 - .C("psmirnov2x", p = as.double(ks), 
              as.integer(ess[2]), as.integer(ess[1]), PACKAGE="twang")$p              

      if((sum(is.na(design$variables$x))>0) && (na.action=="level"))
      {
         # for continuous, separate out for NAs
         work <- design$variables
         work$w[work$t==1] <- with(subset(work,t==1),  w/sum(w))
         work$w[work$t==0] <- with(subset(work,t==0), -w/sum(w))
         ks <- c(ks,
                 abs(sum(with(subset(work, is.na(x)), sapply(split(w,t),sum)))))
         pval <- c(pval,
                   as.numeric(svychisq(~miss+t,design=design)$p.value))
      }
      ret <- cbind(ret,ks,ks.pval=pval)
   }

   ret <- data.frame(ret)
   rownames(ret) <- c("","<NA>")[1:nrow(ret)]

   return(ret)
}

if(estimand=="ATE") {
    design <- svydesign(ids = ~1, weights = ~w, data = data.frame(x = x, 
        t = t, w = w, miss = is.na(x)))
    if (na.action == "exclude") 
        design <- subset(design, !is.na(x))
 
    ret <- NULL
    if (get.means) {
        design.t <- subset(design, t == 1)
        if(multinom) design.c <- design
        else 
        design.c <- subset(design, t == 0)
        m.t <- svymean(~x, design.t, na.rm = TRUE)
        m.c <- svymean(~x, design.c, na.rm = TRUE)
        sd.t <- sqrt(svyvar(~x, design.t, na.rm = TRUE))
        sd.c <- sqrt(svyvar(~x, design.c, na.rm = TRUE))
        t.n <- summary(svyglm(x ~ t, design))$coefficients[2, 
            3:4]
            
        if(multinom){
        	
        	m.unwt <- mean(x, na.rm = TRUE)
        	sd.unwt <- poolSD(x, t)
        }    

        m.t.unwt=mean(x[t==1],na.rm=TRUE)
        if(!multinom) m.c.unwt=mean(x[t==0],na.rm=TRUE)
        else m.c.unwt <- mean(x, na.rm = TRUE)
        var.t<-var(x[t==1],na.rm=TRUE)
        var.c<-var(x[t==0],na.rm=TRUE)
        N.t=length(x[t==1])
        N.c=length(x[t==0])
        sd.p=sqrt(((N.t-1)*var.t+(N.c-1)*var.c)/(N.t+N.c-2))
        if(multinom){
        	b.n <- ifelse(sd.unwt == 0, NA, (m.t - m.unwt)/sd.unwt)
        }
        else b.n <- ifelse(sd.t == 0, NA, (m.t - m.c)/sd.p)

#######
        ret <- cbind(m.t, sd.t, m.c, sd.c, b.n, t.n[1], t.n[2])
        colnames(ret) <- c("tx.mn", "tx.sd", "ct.mn", "ct.sd", 
            "std.eff.sz", "stat", "p")
        if ((sum(is.na(x)) > 0) && (na.action == "level")) {
            m.t <- svymean(~is.na(x), design.t, na.rm = TRUE)[2]
            m.c <- svymean(~is.na(x), design.c, na.rm = TRUE)[2]
            sd.t <- sqrt(m.t * (1 - m.t))
            sd.c <- sqrt(m.c * (1 - m.c))
            test <- try(summary(svyglm(is.na(x) ~ t, family = quasibinomial, 
                design)), silent = TRUE)
            if (class(test)[1] != "try-error") 
                t.n <- test$coefficients[2, 3:4]
            else t.n <- c(NA, NA)

          # revised 101310
          # b.n <- ifelse(sd.t == 0, NA, (m.t - m.c)/sd.t)

            p.unw = sum(is.na(x))/length(x)
            sd.p = sqrt(p.unw*(1-p.unw))
            b.n <- ifelse(sd.p == 0, NA, (m.t - m.c)/sd.p)

            ret <- rbind(ret, c(m.t, sd.t, m.c, sd.c, b.n, t.n))
        }
    }
    if (get.ks) {
        work <- design$variables
        if (na.action == "lowest") 
            work$x[is.na(work$x)] <- min(work$x, na.rm = TRUE) - 
                1
        if (na.action == "level") 
            work <- subset(work, !is.na(x))
            

      if(multinom){
      	work$w[work$t == 1] <- with(subset(work, t == 1), w/sum(w))
      	work$w[work$t == 0] <- 0
      	work$w <- work$w - 1/(length(work$w))
      	#### should replace above line with survey weights 11/15/12 #####
      }
      else{ 

      work$w[work$t==1] <- with(subset(work,t==1),  w/sum(w))
      work$w[work$t==0] <- with(subset(work,t==0), -w/sum(w))
      }
            
         ess <- with(work, sapply(split(w, t), function(w) {sum(w)^2/sum(w^2)}))
        ind <- order(work$x)
        cumv <- abs(cumsum(work$w[ind]))
        cumv <- cumv[diff(work$x[ind]) != 0]
        ks <- ifelse(length(cumv) > 0, max(cumv), 0)
        pval <- NULL
        if(!multinom){
#        pval <- 1 - .C("psmirnov2x", p = as.double(ks), as.integer(ess[2]), 
#            as.integer(ess[1]), PACKAGE = "stats")$p
      pval <- 1 - .C("psmirnov2x", p = as.double(ks), 
              as.integer(ess[2]), as.integer(ess[1]), PACKAGE = "twang")$p            
            }
        if ((sum(is.na(design$variables$x)) > 0) && (na.action == 
            "level")) {
            work <- design$variables
            if(multinom){
            	work$w[work$t == 1] <- with(subset(work, t == 1), w/sum(w))
            	work$w[work$t == 0] <- 0
            	work$w <- work$w - 1/(length(work$w))
            	### should use sample weights
            	}
            else{ 
            	work$w[work$t==1] <- with(subset(work,t==1),  w/sum(w))
            	work$w[work$t==0] <- with(subset(work,t==0), -w/sum(w))
            	}
#            work$w[work$t == 1] <- with(subset(work, t == 1), w/sum(w))
#            work$w[work$t == 0] <- with(subset(work, t == 0), -w/sum(w))
            ks <- c(ks, abs(sum(with(subset(work, is.na(x)), 
                sapply(split(w, t), sum)))))
            if(!multinom)    
            pval <- c(pval, as.numeric(svychisq(~miss + t, design = design)$p.value))
        }
        ret <- cbind(ret, ks, ks.pval = pval)
    }
    ret <- data.frame(ret)
    rownames(ret) <- c("", "<NA>")[1:nrow(ret)]
    return(ret)
}

}


# takes one factor variable, treatment, and weights and computes balance stats
ps.summary.f <- function(x, t, w,
                         get.means=TRUE, get.ks=TRUE,
                         na.action=c("level","exclude","lowest")[1],
                         collapse.by.var=TRUE, estimand, multinom)
{

if(!(estimand %in% c("ATT","ATE"))) stop("estimand must be either \"ATT\" or \"ATE\".")

if(estimand=="ATT") {

   if((sum(is.na(x)) > 0) && (na.action %in% c("level","lowest")))
   {
      x <- factor(x,levels=c(levels(x),"<NA>")) # exclude=NULL messes svychisq/ks
      x[is.na(x)] <- "<NA>"
   }
   design <- svydesign(ids=~1, weights=~w, 
                       data=data.frame(x=x,t=t,w=w,miss=is.na(x)))

   if(na.action=="exclude") design <- subset(design, !is.na(x))

   ret <- NULL

   if(get.means){
      p   <- svytable(~x+t,design,Ntotal=1)
      p   <- t(t(p)/colSums(p))
      sd  <- sqrt(p*(1-p))
      
        if(multinom){
        	
        	m.unwt <- mean(x, na.rm = TRUE)
        	sd.unwt <- poolSD(x,t)
        }    
      
      
      b.f <- ifelse(sd[,2]==0, NA, (p[,2]-p[,1])/sd[,2])

      test  <- try( svychisq(~x+t,design), silent=TRUE)
      if(class(test)[1] != "try-error")
      {
         stat  <- c(test$statistic,rep(NA,nrow(p)-1))
         pval <- c(test$p.value,rep(NA,nrow(p)-1))
      } else
      {
         stat <- pval <- NA
      }

      ret  <- cbind(tx.mn=p[,2], tx.sd=sd[,2],
                    ct.mn=p[,1],  ct.sd=sd[,1],
                    std.eff.sz=b.f,
                    stat=stat, p=pval)
   }
   if(get.ks){
      work <- design$variables
      if(multinom){
      	work$w[work$t == 1] <- with(subset(work, t == 1), w/sum(w))
      	work$w[work$t == 0] <- 0
      	work$w <- work$w - 1/(length(work$w))
      }
      else{ 

      work$w[work$t==1] <- with(subset(work,t==1),  w/sum(w))
      work$w[work$t==0] <- with(subset(work,t==0), -w/sum(w))
      }

      
#      work$w[work$t==1] <- with(subset(work,t==1),  w/sum(w))
#      work$w[work$t==0] <- with(subset(work,t==0), -w/sum(w))
      ks <- abs(sapply(split(work$w,work$x),sum))
      if(multinom){
      	pval <- NA
      }
      else{
      	if(sum(ks>0)<=1){ # deal with factors with some empty levels
      		ks[1:length(ks)]   <- 0  # preserves names(ks)
      		pval <- 1
      		} 
      		else{
      			pval <- try(as.numeric(svychisq(~x+t,design=design)$p.value), 
                     silent=TRUE)
                if(class(pval)[1] == "try-error"){
                	pval <- NA
                	}
                }
                }
                ret <- cbind(ret, ks, ks.pval=c(pval,rep(NA,length(ks)-1)))
   			}
   	return(data.frame(ret))
}

if(estimand=="ATE"){
    if ((sum(is.na(x)) > 0) && (na.action %in% c("level", "lowest"))) {
        x <- factor(x, levels = c(levels(x), "<NA>"))
        x[is.na(x)] <- "<NA>"
    }
    design <- svydesign(ids = ~1, weights = ~w, data = data.frame(x = x, 
        t = t, w = w, miss = is.na(x)))
    if (na.action == "exclude") 
        design <- subset(design, !is.na(x))
    ret <- NULL
    if (get.means) {
        p <- svytable(~x + t, design, Ntotal = 1)
        p <- t(t(p)/colSums(p))
        sd <- sqrt(p * (1 - p))
####################
#new denom for ATE - do I need tildas? before the xs
#####################
	N.t=length(x[t==1])
	N.c=length(x[t==0])
	p.unwt <- table(x)
        p.unwt <- t(t(p.unwt))/sum(p.unwt)
        sd.p <- sqrt(p.unwt * (1 - p.unwt))	
 	b.f <- ifelse(sd[, 2] == 0, NA, (p[, 2] - p[, 1])/sd.p)
#######
        test <- try(svychisq(~x + t, design), silent = TRUE)
        if (class(test)[1] != "try-error") {
            stat <- c(test$statistic, rep(NA, nrow(p) - 1))
            p.val <- c(test$p.value, rep(NA, nrow(p) - 1))
        }
        else {
            stat <- p.val <- NA
        }
        ret <- cbind(tx.mn = p[, 2], tx.sd = sd[, 2], ct.mn = p[, 
            1], ct.sd = sd[, 1], std.eff.sz = b.f, stat = stat, 
            p = p.val)
    }
    if (get.ks) {
        work <- design$variables
        work$w[work$t == 1] <- with(subset(work, t == 1), w/sum(w))
        work$w[work$t == 0] <- with(subset(work, t == 0), -w/sum(w))
        ks <- abs(sapply(split(work$w, work$x), sum))
        if (sum(ks > 0) <= 1) {
            ks[1:length(ks)] <- 0
            pval <- 1
        }
        else {
            pval <- try(as.numeric(svychisq(~x + t, design = design)$p.value), 
                silent = TRUE)
            if (class(pval)[1] == "try-error") {
                p.val <- NA
            }
        }
        ret <- cbind(ret, ks, ks.pval = c(pval, rep(NA, length(ks) - 
            1)))
    }
    return(data.frame(ret))
}

}

poolSD <- function(x, t){
	t <- t[!is.na(x)]
	x <- x[!is.na(x)]
	vars <- unlist(lapply(split(x, t), var))
	tTab <- table(t)
	sumVar <- sum(vars * tTab)
	sd <- sqrt(sumVar/(sum(t) - length(tTab)))
	return(sd)
	
}



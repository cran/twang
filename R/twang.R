####### casemix.R
.First.lib <- function(lib, pkg)
{
     require(gbm)
     require(survey)
}

########################################################



### Function returns ps and balance information
ps<-function(formula = formula(data),
             data,                         # data
             sampw=rep(1,nrow(data)),       # sampling weights
             title=NULL,
             stop.method=stop.methods[1:2], # stopping rules
             plots=TRUE,
             n.trees=10000,                 # gbm options
             interaction.depth=3,
             shrinkage=0.01,
             perm.test.iters=0,
             print.level=2,                 # direct optimizer options
             iterlim=1000,
             verbose=TRUE)
{
   terms <- match.call()
   
   # all this is just to extract the variable names
   mf <- match.call(expand.dots = FALSE)
   m <- match(c("formula", "data"), names(mf), 0)
   mf <- mf[c(1, m)]
   mf[[1]] <- as.name("model.frame")
   mf$na.action <- na.pass
   mf$subset <- rep(FALSE, nrow(data)) # drop all the data
   mf <- eval(mf, parent.frame())
   Terms <- attr(mf, "terms")
   var.names <- attributes(Terms)$term.labels
   treat.var <- as.character(formula[[2]])

   # create the desc object. This holds information on variable balance
   if(class(stop.method)=="stop.method") 
   {
      stop.method <- list(stop.method)
   }
   stop.method.names <- sapply(stop.method,function(x){x$name})
   print(stop.method.names)
   i <- which(sapply(stop.method,function(x){x$direct}))
   if(length(i)>0) cat(paste("*** WARNING: Stop method",stop.method.names[i],"involves direct\noptimization, which is very time consuming, especially\nfor datasets with more than a few hundred cases. Consider\nusing another stop.method or be prepared for a very long wait. ***\n\n"))
   desc <- vector("list",1+length(stop.method))
   names(desc) <- c("unw",stop.method.names)

   # create the plot.info object to reproduce the plots later
   plot.info <- vector("list",1+length(stop.method))
   names(plot.info) <- c("unw",stop.method.names)

   # allocate space for the propensity scores and weights
   p.s        <- data.frame(matrix(NA,nrow=nrow(data),
                                      ncol=length(stop.method)))
   names(p.s) <- stop.method.names
   w          <- data.frame(matrix(NA,nrow=nrow(data),
                                      ncol=length(stop.method)))
   names(w)   <- stop.method.names
   
 
   # alert.stack collects all the warnings
   alerts.stack <- textConnection("alert","w")

   # fit the propensity score model
   if(verbose) cat("Fitting gbm model\n")
   # need to reformulate formula to use this environment
   form <- paste(deparse(formula, 500), collapse="") 

   gbm1 <-gbm(formula(form),
              data = data,
              weights=sampw,
              distribution = "bernoulli",
              n.trees = n.trees,
              interaction.depth = interaction.depth,
              n.minobsinnode = 10,
              shrinkage = shrinkage,
              bag.fraction = 0.5,
              train.fraction = 1,
              verbose = verbose,
              keep.data = FALSE)

   if(verbose) cat("Diagnosis of unweighted analysis\n")
   
   desc$unw <- desc.wts(data=data[,c(treat.var,var.names)],
                        treat.var=treat.var,
                        w=sampw,
                        tp="unw",
                        na.action="level",
                        perm.test.iters=perm.test.iters,
                        verbose=verbose,
                        alerts.stack=alerts.stack)
   desc$unw$n.trees <- NA

   if(plots) pdf(file=paste(title,".pdf",sep=""))

   for(i.tp in 1:length(stop.method))
   {
      tp <- stop.method.names[i.tp]
      if(verbose) cat("Optimizing with",tp,"stopping rule\n")

      # get optimal number of iterations
      # Step #1: evaluate at 25 equally spaced points
      iters <- round(seq(1,gbm1$n.trees,length=25))
      bal <- rep(0,length(iters))
      for(j in 1:length(iters))
      {
         bal[j] <- metric.i(
                    iters[j],
                    fun          = match.fun(stop.method[[i.tp]]$metric),
                    vars         = var.names,
                    treat.var    = treat.var,
                    data         = data,
                    sampw        = sampw,
                    rule.summary = match.fun(stop.method[[i.tp]]$rule.summary),
                    na.action    = stop.method[[i.tp]]$na.action,
                    gbm1         = gbm1)
      }
      if(plots) 
      {
         plot(iters,bal,xlab="Iteration",ylab=stop.method.names[i.tp])
      }
      plot.info[[tp]] <- data.frame(iteration=iters,balance=bal)

      # Step #2: find the interval containing the approximate minimum
      interval <- which.min(bal) +c(-1,1)
      interval[1] <- max(1,interval[1])
      interval[2] <- min(length(iters),interval[2])
   
      # Step #3: refine the minimum by searching with the identified interval
      opt<-optimize(metric.i,
                    interval=iters[interval],
                    maximum   = FALSE,
                    tol       = 1,
                    fun       = match.fun(stop.method[[i.tp]]$metric),
                    vars      = var.names,
                    treat.var = treat.var,
                    data      = data,
                    sampw     = sampw,
                    rule.summary = match.fun(stop.method[[i.tp]]$rule.summary),
                    na.action = stop.method[[i.tp]]$na.action,
                    gbm1      = gbm1)
      if(verbose) cat("   Optimized at",round(opt$minimum),"\n")

      # compute propensity score weights
      p.s[,i.tp]  <- predict(gbm1,newdata=data,
                             n.trees=round(opt$minimum),
                             type="response")
      w[,i.tp] <- p.s[,i.tp]/(1-p.s[,i.tp])
      w[data[,treat.var]==1,i.tp] <- 1
      # adjust for sampling weights
      w[,i.tp] <- w[,i.tp]*sampw 

      # Directly optimize the weights if requested
      if(stop.method[[i.tp]]$direct)
      {
         if(verbose) cat("   Proceeding to direct optimization\n")

         obj   <- nlm(match.fun(stop.method[[i.tp]]$metric),
                      p           = log(w[data[,treat.var]==0,i.tp]),
                      fscale      = 0.1,
                      data        = data,
                      vars        = var.names,
                      treat.var   = treat.var,
                      ndigit      = 3,
                      rule.summary = match.fun(stop.method[[i.tp]]$rule.summary),
                      na.action   = stop.method[[i.tp]]$na.action,
                      print.level = print.level,
                      iterlim     = iterlim)
         if(obj$code %in% c(4,5))
            warning("Failed to completely optimize metric directly.")

         w[data[,treat.var]==0,i.tp] <- exp(obj$estimate)
      }

      if(verbose) cat("Diagnosis of",tp,"weights\n")
      desc[[tp]] <- desc.wts(data[,c(treat.var,var.names)],
                             treat.var    = treat.var,
                             w            = w[,i.tp],
                             tp           = type,
                             na.action    = stop.method[[i.tp]]$na.action,
                             perm.test.iters=perm.test.iters,
                             verbose=verbose,
                             alerts.stack = alerts.stack)
      desc[[tp]]$n.trees <- 
         ifelse(stop.method[[i.tp]]$direct, NA, round(opt$minimum))
      
      if(plots)
      {
         diag.plot(til      = paste(title,"-",tp),
                   treat    = data[,treat.var],
                   p.s      = p.s[,i.tp],
                   w.ctrl   = w[data[,treat.var]==0,i.tp],
                   desc.unw = desc$unw,
                   desc.w   = desc[[tp]])
      }
   }

   if(plots) dev.off()

   close(alerts.stack)
   if(verbose) cat(alert,sep="\n")

   result <- list(gbm.obj    = gbm1,
                  treat      = data[,treat.var], # needed for plot.ps
                  treat.var  = treat.var,
                  desc       = desc,
                  ps         = p.s,
                  w          = w,
                  plot.info  = plot.info,
                  datestamp  = date(),
                  parameters = terms,
                  alerts     = alert)
   class(result) <- "ps"
   return(result)
}



########################################################
desc.wts <-function(data,w,
                    vars=NULL,
                    treat.var,
                    tp,
                    na.action="level",
                    perm.test.iters=0,
                    verbose=TRUE,
                    alerts.stack)
{
   if(is.null(vars)) vars <- names(data)[names(data)!=treat.var]
   ess    <- (sum(w[data[,treat.var]==0])^2)/sum(w[data[,treat.var]==0]^2)

   bal.tab   <- bal.stat(data=data,w.all=w,
                         vars=vars,
                         treat.var=treat.var,
                         na.action=na.action)
   pval.maxks <- NA
   # compute permutation p-values for KS statistic
   if(perm.test.iters>0)
   {
      ess.t <- (sum(w[data[,treat.var]==1])^2)/
               sum(w[data[,treat.var]==1]^2)
      ess.c <- (sum(w[data[,treat.var]==0])^2)/
               sum(w[data[,treat.var]==0]^2)
      w.1   <- w
      w.1[data[,treat.var]==1] <- 1/(ess.t+ess.c)
      w.1[data[,treat.var]==0] <- ess.c*w[data[,treat.var]==0]/
                                  sum(w[data[,treat.var]==0])
      w.1[data[,treat.var]==0] <- w.1[data[,treat.var]==0]/
                                  (ess.c+sum(data[,treat.var]==1))

      pval.var   <- rep(0,nrow(bal.tab$results))
      names(pval.var) <- rownames(bal.tab$results)
      pval.maxks <- 0
      if(verbose)
      {
         cat("Permutation test progress:\n")
         progress <- round(perm.test.iters*c(0.01,(1:10)/10))
      }
      for(i.rep in 1:perm.test.iters)
      {
         if(verbose & (i.rep %in% progress))
         {
            cat(round(100*i.rep/perm.test.iters),"%\n",sep="")
         }
         i <- sample(1:nrow(data),ess.t+ess.c,replace=TRUE,prob=w.1)
         temp <- data[i,c(vars,treat.var)]
         temp[,treat.var] <- as.numeric((1:nrow(temp))<=ess.t)
         bal.temp   <- bal.stat(data=temp,
                                w.all=rep(1,length(i)),
                                vars=vars,
                                treat.var=treat.var,
                                na.action=na.action,
                                get.means=FALSE)
         # if a bootstrap sample is missing a level need to make sure
         #   these still align
         j <- match(rownames(bal.temp$results),names(pval.var))
         pval.var[j] <- pval.var[j]   + 
                        (bal.temp$results$ks >= bal.tab$results$ks[j])
         pval.maxks  <- pval.maxks + 
                        (max(bal.temp$results$ks) >= max(bal.tab$results$ks))
      }
      pval.var   <- as.numeric(pval.var/perm.test.iters)
      pval.maxks <- as.numeric(pval.maxks/perm.test.iters)
      
      # replace the analytic p-values computed in bal.tab
      bal.tab$results$ks.pval <- pval.var
   }
   
   check.err(cov.table=bal.tab$results, stage=tp, alerts.stack=alerts.stack)

   max.ks  <- max(bal.tab$results$ks)
   mean.ks  <- mean(bal.tab$results$ks)
   max.es <- with(bal.tab$results, max(abs(std.eff.sz[std.eff.sz<500]),
                 na.rm=TRUE))
   mean.es   <- with(bal.tab$results, mean(abs(std.eff.sz[std.eff.sz<500]),
                  na.rm=TRUE))
   return(list(ess=ess,
               n.treat=sum(data[,treat.var]==1),
               n.ctrl =sum(data[,treat.var]==0),
               max.es=max.es,
               mean.es=mean.es,
               max.ks=max.ks,
               max.ks.p=pval.maxks,
               mean.ks=mean.ks,
               bal.tab=bal.tab))
}


# Rearranges the arguments of the *.stat functions
#    Used in conjunction with optimize
metric.i <- function(i,fun=ks.stat,...)
{
   return(fun(i=round(i),...))
}


# computes KS statistics from weights given in a gbm object, log(w), or control weights
ks.stat<-function(logw=NULL,
                  w.ctrl=NULL,
                  gbm1=NULL,
                  i=1,
                  data,
                  sampw=rep(1,nrow(data)),
                  rule.summary=mean,
                  na.action="level",
                  vars,
                  treat.var,
                  collapse.by.var=FALSE,
                  verbose=FALSE)
{
   if(is.null(gbm1) && is.null(w.ctrl) && is.null(logw))
      stop("No weights given. logw, gbm1, and w.ctrl cannot all be NULL.")
   if(!is.null(rule.summary)) rule.summary <- match.fun(rule.summary)
   w1 <- rep(1/sum(data[,treat.var]==1), nrow(data))
   if(is.null(gbm1))
   {
      if (!is.null(logw))
      {
         w.ctrl<-exp(logw)
      }
      w1[data[,treat.var]==0] <- w.ctrl
   } else
   {
      w <- exp(predict(gbm1,newdata=data[data[,treat.var]==0,],
                           n.trees=i))
      w1[data[,treat.var]==0]<- w
   }

   # compute KS statistics
   w1 <- w1*sampw
   ks <- lapply(data[,vars], ps.summary,
                t=data[,treat.var],
                w=w1,
                get.means=FALSE,
                get.ks=TRUE,
                na.action=na.action)
   ks <- lapply(ks, function(x){x$ks})

   if(collapse.by.var) ks <- sapply(ks,max)
   else ks <- unlist(ks)

   if(!is.null(rule.summary))
   {
      if(verbose) print(rule.summary(ks))
      return(rule.summary(ks))
   } else
   {
      return(ks)
   }
}



es.stat<-function(logw=NULL,
                  w.ctrl=NULL,
                  gbm1=NULL,
                  i=1,
                  data,
                  sampw=rep(1,nrow(data)),
                  rule.summary=mean,
                  na.action="level",
                  vars,
                  treat.var,
                  collapse.by.var=FALSE,
                  verbose=FALSE)
{
   if(is.null(gbm1) && is.null(w.ctrl) && is.null(logw))
      stop("No weights given. logw, gbm1, and w.ctrl cannot all be NULL.")
   if(!is.null(rule.summary)) rule.summary <- match.fun(rule.summary)
   w1 <- rep(1/sum(data[,treat.var]==1), nrow(data))
   if(is.null(gbm1))
   {
      if (!is.null(logw))
      {
         w.ctrl<-exp(logw)
      }
      w1[data[,treat.var]==0] <- w.ctrl
   } else
   {
      w <- exp(predict(gbm1,newdata=data[data[,treat.var]==0,],
                           n.trees=i))
      w1[data[,treat.var]==0] <- w
   }
   w1 <- w1*sampw

   # compute effect sizes
   es <- lapply(data[,vars], ps.summary,
                t=data[,treat.var],
                w=w1,
                get.means=TRUE,
                get.ks=FALSE,
                na.action=na.action)

   if(collapse.by.var)
   {
      es <- sapply(es,function(x,rule.summary){rule.summary(x$std.eff.sz,na.rm=TRUE)},
                   rule.summary=rule.summary)
   } else
   {
      es <- unlist(sapply(es,function(x){x$std.eff.sz}))
   }

   if(!is.null(rule.summary))
   {
      if(verbose) print(rule.summary(es,na.rm=TRUE))
      return(rule.summary(abs(es),na.rm=TRUE))
   } else
   {
      return(es)
   }
}


strata.stat <- function(logw = NULL, 
                        w.ctrl = NULL, 
                        gbm1 = NULL, 
                        i = 1, 
                        data, 
                        sampw = rep(1, nrow(data)), 
                        rule.summary = mean, 
                        na.action = "level", 
                        vars, 
                        treat.var, 
                        collapse.by.var = FALSE, 
                        verbose = FALSE)
{
   if(is.null(gbm1))
         stop("strata.stat() only works with a non-null gbm1.")
   if(!is.null(rule.summary)) rule.summary <- match.fun(rule.summary)
   
   p      <- predict(gbm1, newdata=data, n.trees=i, type="response")

   # generate cut points at quintiles of the treatment propensity scores
   cuts   <- c(0, quantile(p[data[,treat.var]==1], c(0.0,0.2,0.4,0.6,0.8,1.0)))
   strata <- as.numeric(cut(p,cuts,right=FALSE))
   strata[is.na(strata)]<-6 # fix observations at the largest p
   strata.n <- table(strata,data[,treat.var])
   if(any(strata.n[,1]==0))
   {
      return(Inf)
   }else
   {
      strata.w    <- strata.n[,2]/strata.n[,1]
      data$w      <- strata.w[strata]
      data$strata <- strata
      design.t <- svydesign(ids=~1,weights=~w,data=data[data[,treat.var]==1,])
      design.c <- svydesign(ids=~1,weights=~w,data=data[data[,treat.var]==0,])

      # drop those controls below the smallest treatment p
      #   strata 1 are those prop scores between 0 and min treat prop score
      design.t <- subset(design.t, strata!=1)
      design.c <- subset(design.c, strata!=1)

      m.t.n <- sapply(vars,FUN=function(n){
                      svymean(formula(paste("~",n)),design.t,na.rm=TRUE)})
      m.c.n <- sapply(vars,FUN=function(n){
                      svymean(formula(paste("~",n)),design.c,na.rm=TRUE)})

      b <- numeric(length(m.t.n))
      is.fac <- sapply(data[,vars], is.factor)
      for(j in 1:length(m.t.n))
      {
         if(!is.fac[j])
         {
            sd.t <- sqrt(svyvar(formula(paste("~",names(m.t.n)[j])),design.t,na.rm=TRUE))
            if(sd.t!=0) b[j] <- as.numeric( abs(m.t.n[[j]]-m.c.n[[j]])/sd.t )
            else        b[j] <- NA
         } else
         {
            sd.t <- sqrt(m.t.n[[j]]*(1-m.t.n[[j]]))
            b[j] <- as.numeric( rule.summary(abs(m.t.n[[j]]-m.c.n[[j]])/sd.t,na.rm=TRUE) )
         }
      }

      return( rule.summary(abs(b),na.rm=TRUE))
   }
}


stop.methods <- list(ks.stat.mean=       list(metric=ks.stat,
                                              rule.summary=mean,
                                              direct=FALSE,
                                              na.action="level",
                                              name="ks.stat.mean"),
                     es.stat.mean=       list(metric=es.stat,
                                              rule.summary=mean,
                                              direct=FALSE,
                                              na.action="level",
                                              name="es.stat.mean"),
                     ks.stat.max=        list(metric=ks.stat,
                                              rule.summary=max,
                                              direct=FALSE,
                                              na.action="level",
                                              name="ks.stat.max"),
                     es.stat.max=        list(metric=es.stat,
                                              rule.summary=max,
                                              direct=FALSE,
                                              na.action="level",
                                              name="es.stat.max"),
                     ks.stat.max.direct= list(metric=ks.stat,
                                              rule.summary=max,
                                              direct=TRUE,
                                              na.action="level",
                                              name="ks.stat.max.direct"),
                     es.stat.max.direct= list(metric=es.stat,
                                              rule.summary=max,
                                              direct=TRUE,
                                              na.action="level",
                                              name="es.stat.max.direct"),
                     strat.max=          list(metric=strata.stat,
                                              rule.summary=max,
                                              direct=FALSE,
                                              na.action="level",
                                              name="strat.max"))
for(i in 1:length(stop.methods)) class(stop.methods[[i]]) <- "stop.method"
rm(i)

# dispatches ps.summary.n or ps.summary.f, depending on x's type
ps.summary <- function(x, t, w,
                       get.means=TRUE, get.ks=TRUE,
                       na.action=c("level","exclude","lowest")[1],
                       collapse.by.var=FALSE)
{
   if(!is.factor(x))
      return(ps.summary.n(x=x, t=t, w=w,
                         get.means=get.means, get.ks=get.ks,
                         na.action=na.action,
                         collapse.by.var=collapse.by.var))
   else
      return(ps.summary.f(x=x, t=t, w=w,
                         get.means=get.means, get.ks=get.ks,
                         na.action=na.action,
                         collapse.by.var=collapse.by.var))
}


# takes one numeric variable, treatment, and weights and computes balance stats
ps.summary.n <- function(x, t, w,
                         get.means=TRUE, get.ks=TRUE,
                         na.action=c("level","exclude","lowest")[1],
                         collapse.by.var=FALSE)
{
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

   if(get.ks)
   {
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
      pval <- 1 - .C("psmirnov2x", p = as.double(ks), 
              as.integer(ess[2]), as.integer(ess[1]), PACKAGE = "stats")$p

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


# takes one factor variable, treatment, and weights and computes balance stats
ps.summary.f <- function(x, t, w,
                         get.means=TRUE, get.ks=TRUE,
                         na.action=c("level","exclude","lowest")[1],
                         collapse.by.var=TRUE)
{
   if((sum(is.na(x)) > 0) &&
      (na.action %in% c("level","lowest")))
   {
      x <- factor(x,levels=c(levels(x),"<NA>")) # exclude=NULL messes svychisq/ks
      x[is.na(x)] <- "<NA>"
   }
   design <- svydesign(ids=~1, weights=~w, 
                       data=data.frame(x=x,t=t,w=w,miss=is.na(x)))

   if(na.action=="exclude") design <- subset(design, !is.na(x))

   ret <- NULL

   if(get.means)
   {
      p   <- svytable(~x+t,design,Ntotal=1)
      p   <- t(t(p)/colSums(p))
      sd  <- sqrt(p*(1-p))
      b.f <- ifelse(sd[,2]==0, NA, (p[,2]-p[,1])/sd[,2])

      test  <- try( svychisq(~x+t,design), silent=TRUE)
      if(class(test)[1] != "try-error")
      {
         stat  <- c(test$statistic,rep(NA,nrow(p)-1))
         p.val <- c(test$p.value,rep(NA,nrow(p)-1))
      } else
      {
         stat <- p.val <- NA
      }

      ret  <- cbind(tx.mn=p[,2], tx.sd=sd[,2],
                    ct.mn=p[,1],  ct.sd=sd[,1],
                    std.eff.sz=b.f,
                    stat=stat, p=p.val)
   }
   if(get.ks)
   {
      work <- design$variables
      work$w[work$t==1] <- with(subset(work,t==1),  w/sum(w))
      work$w[work$t==0] <- with(subset(work,t==0), -w/sum(w))
      ks <- abs(sapply(split(work$w,work$x),sum))
      if(sum(ks>0)<=1)
      { # deal with factors with some empty levels
         ks[1:length(ks)]   <- 0  # preserves names(ks)
         pval <- 1
      } else
      {
         pval <- as.numeric(svychisq(~x+t,design=design)$p.value)
      }
      ret <- cbind(ret, ks, ks.pval=c(pval,rep(NA,length(ks)-1)))
   }

   return(data.frame(ret))
}



### calculate weighted balance statistic
bal.stat <- function(data,vars=NULL,treat.var,w.all,
                     get.means=TRUE,
                     get.ks=TRUE,
                     na.action="level")
{
   if(is.null(vars)) vars<-names(data)[names(data)!=treat.var]

   is.fac   <- sapply(data[,vars,drop=FALSE],is.factor)
   fac      <- vars[is.fac]
   not.fac  <- vars[!is.fac]

   ret <- vector("list",length(vars))
   names(ret) <- vars

   ##### Calculate stats for numeric variables
   ret[!is.fac] <- lapply(data[,vars[!is.fac],drop=FALSE], ps.summary.n,
                          t=data[,treat.var], w=w.all,
                          get.means=get.means, get.ks=get.ks,
                          na.action=na.action,
                          collapse.by.var=FALSE)

   ##### Calculate stats for factor variables
   ret[is.fac] <- lapply(data[,vars[is.fac],drop=FALSE], ps.summary.f,
                         t=data[,treat.var], w=w.all,
                         get.means=get.means, get.ks=get.ks,
                         na.action=na.action,
                         collapse.by.var=FALSE)
   # this keeps the variables in the same order as vars
   n.rows <- sapply(ret,nrow)
   var.levels <- unlist(sapply(ret, rownames))
   var.names <- rep(names(ret),n.rows)
   var.names[var.levels!=""] <- paste(var.names[var.levels!=""],
                                      var.levels[var.levels!=""],sep=":")

   res <- data.frame(matrix(0,nrow=length(var.names), ncol=ncol(ret[[1]])))
   names(res) <- colnames(ret[[1]])
   rownames(res) <- var.names

   # populate the results table
   i.insert <- 1
   for(i in 1:length(ret))
   {
      res[i.insert:(i.insert+nrow(ret[[i]])-1),] <- ret[[i]]
      i.insert <- i.insert+nrow(ret[[i]])
   }

   res <- list(results=data.frame(res))
   return(res)
}


# need a better description of what this is doing.
# replace function name with something more meaningful
check.err<-function(cov.table, stage, alerts.stack)
{
   ind  <- (cov.table$tx.sd < .0001) | (cov.table$std.ef.sz > 500)
   prob <- cov.table$std.eff.sz[ind]
   if(length(prob)>0)
   {
      sink(alerts.stack, append=TRUE)
      cat("\n problematic standard deviations in stage ",stage,"\n\n")
      print(cov.table[which(ind),c("tx.sd","std.eff.sz")])
      cat("\n\n\n")
      sink()
   }
}


diag.plot <- function(til,treat,p.s,w.ctrl,desc.unw,desc.w)
{
   par(mfrow=c(2,2))

   if(!is.null(p.s))
   {
      boxplot(t(p.s)~t(treat),
              names=c("Control","Treatment"),
              ylab="Propensity Scores",xlab="Treatment assignment",main=til)
   }

   h1<-hist(w.ctrl, xlab="Weight",main="Control Weights")
   h2<-max(h1$counts)-(1:3)*.1*max(h1$counts)

   labels <- c(paste("Tx.n=",sum(treat==1),sep=""),
                 paste("Ct.n=",sum(treat==0),sep=""),
                 paste("ESS=", round(desc.w$ess,1),sep=""))
   text(x=par()$usr[2]-max(strwidth(labels)), 
        y=par()$usr[4],
        labels=paste(labels,collapse="\n"),
        col="dark green",
        adj=c(0,1))        

   stats.unw <- desc.unw$bal.tab$results
   stats.w   <- desc.w$bal.tab$results
   plot(sort(stats.unw$p), col = "red", pch = 20,
        xlab="pretreatment covariates", ylim=c(0,1),
        main="T-test p-values", ylab="red (unw), black (wt)")
   points(sort(stats.w$p))
   lines(c(1,sum(!is.na(stats.w$p))),c(0,1), col="blue")

   # KS p-value plot
   plot(sort(stats.unw$ks.pval), col="red", pch=20,
             xlab="pretreatment covariates", ylim=c(0,1),
             main="KS-test p-values", ylab="red (unw), black (wt)")
   points(sort(stats.w$ks.pval))
   lines(c(1,sum(!is.na(stats.w$ks.pval))), c(0,1), col="blue")

   # ES spaghetti plot
   par(mfrow=c(1,1))
   ases.dat <- data.frame(es.unw=stats.unw$std.eff.sz, es.w=stats.w$std.eff.sz,
                          p.unw =stats.unw$p,          p.w =stats.w$p)
   ases.dat <- abs(subset(ases.dat, !is.na(p.unw)))

   plot(c(0.85,2.15),c(0,min(3,max(unlist(ases.dat[,1:2]), na.rm=TRUE))), type="n",
        xaxt="n",ylab="Absolute Std Difference",xlab="",
        main=til, sub="Std Effect Sizes")
   abline(h=c(.2,.5,.8))
   axis(side=1, at=1:2, labels=c("Unweighted", "Weighted"))
   for(i in 1:nrow(ases.dat))
   {
      points(1:2,abs(ases.dat[i,c("es.unw","es.w")]),type="b",col="skyblue")
   }
   temp1 <- ases.dat[abs(ases.dat$es.unw) < abs(ases.dat$es.w),]
   for(i in 1:nrow(temp1))
   {
         points(1:2,abs(temp1[i,c("es.unw","es.w")]),type="b", col="red",lwd=2)
   }
   ind <- which(ases.dat$p.unw < 0.05)
   points(rep(1,length(ind)), ases.dat$es.unw[ind],pch=19, col="red")
   ind <- which(ases.dat$p.w   < 0.05)
   points(rep(2,length(ind)), ases.dat$es.w[ind],  pch=19, col="red")
   if (max(ases.dat$es.w,na.rm=TRUE)>3) mtext(text="Some weighted effects>3 !!",side=3,col="red")

}


dx.wts <- function(x,
                   data,
                   vars=NULL,
                   treat.var,
                   x.as.weights=TRUE,
                   sampw=NULL,
                   perm.test.iters=0,
                   plots=TRUE,
                   title=format(Sys.time()+10000, "%Y-%m-%d+%H%M"))
{
   if(class(x)!="ps")
   {
      if(!is.vector(x) && !is.matrix(x) && !is.data.frame(x))
         stop("x must be a ps object, a vector, a matrix, or a data frame")
      if(any(x < 0)) stop("x has negative values")
      if(is.null(dim(x))) x <- matrix(x,ncol=1)
      if(nrow(x) != nrow(data)) stop("length(x) != nrow(data)")
      if(x.as.weights)
      {
         w <- x
         p.s <- 1/(1+x)
      } else
      {
         if(any(x > 1) && !x.as.weights) stop("x has values greater than 1. With x.as.weights=FALSE x should be a vector of propensity scores.")
         w <- matrix(1,nrow=nrow(x),ncol=ncol(x))
         i <- data[,treat.var]==0
         w[i,] <- x[i,]/(1-x[i,])
         p.s <- x
      }
      if(any(is.infinite(w))) stop("Some propensity weights are infinite.")
      # add a column for unweighted analysis
   } else
   {
      if(!is.null(sampw)) warning("Sampling weights given when x is a ps object. The sampling weights should be utilized when running ps and are probably not needed here as well and the results may not be correct.")
      # extract the propensity scores and weights from the ps object
      p.s  <- x$ps
      w    <- x$w
      desc <- x$desc
   }
   if(!all(w[,1]==1)) 
   {
      w   <- cbind(unw=rep(1,nrow(w)),w)
      p.s <- cbind(unw=rep(0.5,nrow(p.s)),p.s)
   }
   if(!is.null(sampw)) w <- w*sampw
   if(is.null(vars)) vars <- names(data)[names(data)!=treat.var]

   summary.tab <- NULL   
   zz      <- textConnection("alert","a")
   if(plots) pdf(file=paste(title,".pdf",sep=""))

   n.tp <- ifelse(class(x)=="ps",length(x$desc),ncol(w))
   if(class(x)!="ps")
   { 
     desc<-vector("list",ncol(w))
     names(desc) <- colnames(w)
   }

   for(i.tp in 1:n.tp)
   {
      if(class(x)=="ps")
      {
         desc.temp <- x$desc[[i.tp]]
         iter      <- desc.temp$n.trees
         tp        <- names(x$desc)[i.tp]
      } else
      {
         desc.temp <- desc.wts(data,
                               w=w[,i.tp],
                               vars=vars,
                               treat.var=treat.var,
                               perm.test.iters=perm.test.iters,
                               verbose=TRUE)
         iter <- NA
         tp <- colnames(w)[i.tp]
         desc[[i.tp]] <- desc.temp
      }
      if(is.null(tp)) tp <- "unnamed"

      summary.tab <- rbind(summary.tab,
         with(desc.temp, data.frame(type    = tp,
                                    n.treat = n.treat,
                                    n.ctrl  = n.ctrl,
                                    ess     = ess,
                                    max.es  = max.es,
                                    mean.es = mean.es,
                                    max.ks  = max.ks,
                                    mean.ks = mean.ks,
                                    iter    = iter)))

      if(plots)
      {
         if(i.tp==1)
         {
            desc.unw <- desc.temp
         } else
         {
            diag.plot(til      = paste(title,"-",tp),
                      treat    = data[,treat.var],
                      p.s      = p.s[,i.tp],
                      w.ctrl   = w[data[,treat.var]==0,i.tp],
                      desc.unw = desc.unw,
                      desc.w   = desc.temp)
         }
      }
   } 

   close(zz)
   if(plots) dev.off()

   cat(alert,sep="\n")
   rownames(summary.tab) <- 1:nrow(summary.tab)

   result <- list(treat      = data[,treat.var],
                  desc       = desc,
                  summary.tab = summary.tab,
                  ps         = as.matrix(p.s[,-1]),
                  w          = as.matrix(w[,-1]),
                  datestamp  = date(),
                  parameters = match.call(),
                  alerts     = alert)
   class(result) <- "dxwts"
   return(result)
}

print.dxwts <- function(x,...)
{
   x$summary.tab
}

########################################################
# ps object summary functions added by Dan Mc

plot.ps <- function(x,label="",ask=FALSE, ...)
{
   # Creates diag.plot plots and sends to current device
   # x:     ps object 
   # label: Label added to the plot titles

   # extract the propensity scores and weights from the ps object
   p.s    <- x$ps
   w      <- x$w
   if(!all(w[,1]==1)) 
   {
      w   <- cbind(unw=rep(1,nrow(w)),w)
      p.s <- cbind(unw=rep(0.5,nrow(p.s)),p.s)
   }
   
   par.ask0 <- par()$ask
   par(ask=ask)
   
   n.tp <- ifelse(class(x)=="ps",length(x$desc),ncol(w))
   for(i.tp in 1:n.tp)
   {
      desc.temp <- x$desc[[i.tp]]
      iter      <- desc.temp$n.trees
      tp        <- names(x$desc)[i.tp]

      if(i.tp==1)
      {
         desc.unw <- desc.temp
      } else
      {
         diag.plot(til      = paste(label,tp,sep=""),
                   treat    = x$treat,
                   p.s      = p.s[,i.tp],
                   w.ctrl   = w[x$treat==0,i.tp],
                   desc.unw = desc.unw,
                   desc.w   = desc.temp)
      }
   } 

   par(ask=par.ask0)
   invisible()
}


# Produces a summary table for ps object 
# object: ps object
summary.ps <- function(object,...)
{
  if(class(object)!="ps") stop("object must be a ps object")
  else
  {
      summary.tab <- NULL   
   
      n.tp <- length(object$desc)
      for(i.tp in 1:n.tp)
      {
         desc.temp <- object$desc[[i.tp]]
         iter      <- desc.temp$n.trees
         tp        <- names(object$desc)[i.tp]

         summary.tab <- rbind(summary.tab,
            with(desc.temp, data.frame(type     = tp,
                                       n.treat  = n.treat,
                                       n.ctrl   = n.ctrl,
                                       ess      = ess,
                                       max.es   = max.es,
                                       mean.es  = mean.es,
                                       max.ks   = max.ks,
                                       max.ks.p = max.ks.p,
                                       mean.ks  = mean.ks,
                                       iter     = iter)))
      } 
      return(summary.tab)
   }
}


## Generic function for extracting balance tables from ps and other objects
bal.table <- function(x){
   UseMethod("bal.table")
}


# Table for extracting balance tables from ps objects
# x: ps object
bal.table.ps <- function(x,...)
{  
   bal.tab <- lapply(x$desc, function(x){return(round(x$bal.tab$results,3))})
   return(bal.tab)
}

# Table for extracting balance tables from dxwts objects
# x: dxwts object
bal.table.dxwts <- function(x,...)
{ 
   bal.tab <- lapply(x$desc, function(x){return(round(x$bal.tab$results,3))})
   return(bal.tab)
}

# Creates diag.plot plots and sends to current device
# x:     dxwts object 
# label: Label added to the plot titles
plot.dxwts <- function(x,label="", ...)
{
   # extract the propensity scores and weights from the dxwts object
   p.s    <- x$ps
   w      <- x$w
   if(!all(w[,1]==1)) 
   {
      w   <- cbind(unw=rep(1,nrow(w)),w)
      p.s <- cbind(unw=rep(0.5,nrow(p.s)),p.s)
   }

   n.tp <- ifelse(class(x)=="ps",length(x$desc),ncol(w))
   for(i.tp in 1:n.tp)
   {
      desc.temp <- x$desc[[i.tp]]
      iter      <- desc.temp$n.trees
      tp        <- names(x$desc)[i.tp]

      if(i.tp==1)
      {
         desc.unw <- desc.temp
      } else
      {
         diag.plot(til      = paste(label,tp,sep=""),
                   treat    = x$treat,
                   p.s      = p.s[,i.tp],
                   w.ctrl   = w[x$treat==0,i.tp],
                   desc.unw = desc.unw,
                   desc.w   = desc.temp)
      }
   } 
}

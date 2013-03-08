sensitivity <- function(ps1,data,
                        outcome,
                        order.by.importance=TRUE,
                        verbose=TRUE)
{
   stop.method <- eval(ps1$parameters$stop.method)
   if(class(stop.method)=="stop.method") 
   {
      stop.method <- list(stop.method)
   }
   
   results <- vector("list",length(stop.method))
   for(i.smethod in 1:length(stop.method))
   {
      if(verbose) cat("Sensitivity analysis for",stop.method[[i.smethod]]$name,"\n")
   
      if(order.by.importance)
      {
         best.iter <- ps1$desc[[stop.method[[i.smethod]]$name]]$n.trees
         vars <- as.character(summary(ps1$gbm.obj,n.trees=best.iter,plot=FALSE)$var)
      } else
      {
         vars <- ps1$gbm.obj$var.names
      }
      
      if(is.null(ps1$parameters$sampw))
         ps1$parameters$sampw <- rep(1,nrow(data))
   
      results[[i.smethod]] <- data.frame(var      =vars,
                           E0       =rep(0,length(vars)),
                           a.min    =rep(0,length(vars)),
                           a.max    =rep(0,length(vars)),
                           a.cor    =rep(0,length(vars)),
                           a.mincor =rep(0,length(vars)),
                           a.maxcor =rep(0,length(vars)),
                           minE0  =rep(0,length(vars)),
                           maxE0  =rep(0,length(vars)),
                           breakeven.cor=rep(0,length(vars)))
      if(verbose) cat("Computing sensitivity statistics for:\n")
      for(i.var in 1:length(vars))
      {
         if(verbose) cat(vars[i.var],"\n")
         form <- paste(ps1$treat.var,"~",paste(vars[-i.var],collapse="+"))
      
         ps2 <- ps(formula(form), 
                  data = data,
                  sampw = ps1$parameters$sampw,
                  title = ps1$parameters$title, 
                  stop.method = stop.method[[i.smethod]], 
                  plots = "none", 
                  n.trees = ps1$parameters$n.trees, 
                  interaction.depth = ps1$parameters$interaction.depth, 
                  shrinkage = ps1$parameters$shrinkage, 
                  perm.test.iters = 0,
                  print.level = ps1$parameters$print.level,
                  iterlim = ps1$parameters$iterlim,
                  verbose = verbose)
                  
         # what kind of ai's are there?
         i <- which(data[,ps1$treat.var]==0)
         d0 <- data[i, outcome, drop=FALSE]
         d0$w <- ps1$w[[stop.method[[i.smethod]]$name]][i]
         d0$wa <- ps2$w[[stop.method[[i.smethod]]$name]][i]
         d0$a <- with(d0, w/wa)
   
         design.ps <- svydesign(ids=~1,weights=~wa,data=d0)
         results[[i.smethod]]$E0[i.var] <- 
            as.numeric(svymean(make.formula(outcome),design.ps))
         
         results[[i.smethod]]$a.min[i.var]    <- min(d0$a)
         results[[i.smethod]]$a.max[i.var]    <- max(d0$a)
         results[[i.smethod]]$a.cor[i.var]    <- cor(d0$a,d0[,outcome])
   
         i <- order(d0[,outcome])
         a.sort <- sort(d0$a)
         results[[i.smethod]]$a.mincor[i.var] <- cor(rev(a.sort),d0[i,outcome])
         results[[i.smethod]]$a.maxcor[i.var] <- cor(a.sort,d0[i,outcome])
      
         d0$w1 <- rep(0,nrow(d0))
         # pile the largest weights on the largest outcomes
         d0$w1[i] <- d0$w[i]*a.sort
         design.ps <- svydesign(ids=~1,weights=~w1,data=d0)
         temp <- svymean(make.formula(outcome),design.ps)
         results[[i.smethod]]$maxE0[i.var] <- temp
   
         # pile the largest weights on the shortest durations
         d0$w1[i] <- d0$w[i]*rev(a.sort)
         design.ps <- svydesign(ids=~1,weights=~w1,data=d0)
         temp <- svymean(make.formula(outcome),design.ps)
         results[[i.smethod]]$minE0[i.var] <- temp
   
         # ai partially correlated with duration
         p       <- c(0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.90,0.99)
         eff     <- data.frame(p=c(p,rev(p)))
         eff$rev <- rep(c(TRUE,FALSE),each=nrow(eff)/2)
         eff$rho <- rep(NA,nrow(eff))
         eff$E0  <- rep(NA,nrow(eff))
         b.E0 <- b.cor <- 1:30
         for(i.eff in 1:nrow(eff))
         {
            for(k in 1:30)
            {
                  j <- sample(1:length(a.sort),size=eff$p[i.eff]*length(a.sort))
                  if(eff$rev[i.eff]) a.shuf <- rev(a.sort)
                  else a.shuf <- a.sort
                  a.shuf[j] <- a.shuf[sample(j)]
                  d0$w1[i] <- d0$w[i]*a.shuf
                  
                  design.ps <- svydesign(ids=~1,weights=~w1,data=d0)
                  b.E0[k]  <- as.numeric(svymean(make.formula(outcome),design.ps))
                  b.cor[k] <- cor(a.shuf,d0[i,outcome])
            }
            eff$E0[i.eff]  <- mean(b.E0)
            eff$rho[i.eff] <- mean(b.cor)
         }
         
         design.ps <- svydesign(ids=~1,weights=~w,data=d0)
         E0 <- as.numeric(svymean(make.formula(outcome),design.ps))
         results[[i.smethod]]$breakeven.cor[i.var] <- approx(eff$E0,eff$rho,E0)$y
      }
   }
   
   return(results)
}


### R code from vignette source 'twang.rnw'

###################################################
### code chunk number 1: twang.rnw:105-106
###################################################
options(width=80)


###################################################
### code chunk number 2: twang.rnw:125-127
###################################################
library(twang)
set.seed(1)


###################################################
### code chunk number 3: twang.rnw:135-136
###################################################
data(lalonde)


###################################################
### code chunk number 4: twang.rnw:159-169
###################################################
ps.lalonde <- ps(treat ~ age + educ + black + hispan + nodegree +
                         married + re74 + re75,
                 data = lalonde,
                 n.trees=5000,
                 interaction.depth=2,
                 shrinkage=0.01,
                 perm.test.iters=0,
                 stop.method=c("es.mean","ks.max"),                 
                 estimand = "ATT",
                 verbose=FALSE)


###################################################
### code chunk number 5: iterPt
###################################################
    plot(ps.lalonde)


###################################################
### code chunk number 6: iterPt2
###################################################
    plot(ps.lalonde, subset = 2)


###################################################
### code chunk number 7: twang.rnw:367-370
###################################################
summary(ps.lalonde$gbm.obj,
        n.trees=ps.lalonde$desc$ks.max.ATT$n.trees,
        plot=FALSE)


###################################################
### code chunk number 8: twang.rnw:376-378
###################################################
summary(ps.lalonde$gbm.obj,
        n.trees=ps.lalonde$desc$ks.max.ATT$n.trees)


###################################################
### code chunk number 9: twang.rnw:407-408
###################################################
options(width=85)


###################################################
### code chunk number 10: twang.rnw:411-413
###################################################
lalonde.balance <- bal.table(ps.lalonde)
lalonde.balance


###################################################
### code chunk number 11: twang.rnw:416-417
###################################################
options(width=80)


###################################################
### code chunk number 12: twang.rnw:490-499
###################################################
library(xtable)
pretty.tab <- lalonde.balance$ks.max.ATT[,c("tx.mn","ct.mn","ks")]
pretty.tab <- cbind(pretty.tab, lalonde.balance$unw[,"ct.mn"])
names(pretty.tab) <- c("E(Y1|t=1)","E(Y0|t=1)","KS","E(Y0|t=0)")
xtable(pretty.tab,
       caption = "Balance of the treatment and comparison groups",
       label = "tab:balance",
       digits = c(0, 2, 2, 2, 2),
       align=c("l","r","r","r","r"))


###################################################
### code chunk number 13: twang.rnw:513-514
###################################################
summary(ps.lalonde)


###################################################
### code chunk number 14: twang.rnw:590-591
###################################################
plot(ps.lalonde, plots=2)


###################################################
### code chunk number 15: twang.rnw:636-637
###################################################
plot(ps.lalonde, plots=3)


###################################################
### code chunk number 16: twang.rnw:652-653
###################################################
plot(ps.lalonde, plots = 4)


###################################################
### code chunk number 17: twang.rnw:668-669
###################################################
plot(ps.lalonde, plots = 5)


###################################################
### code chunk number 18: twang.rnw:704-705
###################################################
plot(ps.lalonde, plots = 3, subset = 2)


###################################################
### code chunk number 19: twang.rnw:719-720
###################################################
library(survey)


###################################################
### code chunk number 20: twang.rnw:732-734
###################################################
lalonde$w <- get.weights(ps.lalonde, stop.method="es.mean")
design.ps <- svydesign(ids=~1, weights=~w, data=lalonde)


###################################################
### code chunk number 21: twang.rnw:737-737
###################################################



###################################################
### code chunk number 22: twang.rnw:771-773
###################################################
glm1 <- svyglm(re78 ~ treat, design=design.ps)
summary(glm1)


###################################################
### code chunk number 23: twang.rnw:805-807
###################################################
glm2 <- svyglm(re78 ~ treat + nodegree, design=design.ps)
summary(glm2)


###################################################
### code chunk number 24: twang.rnw:819-823
###################################################
glm3 <- svyglm(re78 ~ treat + age + educ + black + hispan + nodegree +
                      married + re74 + re75,
               design=design.ps)
summary(glm3)


###################################################
### code chunk number 25: twang.rnw:832-836
###################################################
glm4 <- lm(re78 ~ treat + age + educ + black + hispan + nodegree +
                  married + re74 + re75,
           data=lalonde)
summary(glm4)


###################################################
### code chunk number 26: twang.rnw:839-842
###################################################
glm5 <- lm(sqrt(re78) ~ treat + age + educ + black + hispan + nodegree +
                        married + sqrt(re74) + sqrt(re75),
           data=lalonde)


###################################################
### code chunk number 27: twang.rnw:868-874
###################################################
ps.logit <- glm(treat ~ age + educ + black + hispan + nodegree +
                        married + re74 + re75,
                data = lalonde,
                family = binomial)
lalonde$w.logit <- rep(1,nrow(lalonde))
lalonde$w.logit[lalonde$treat==0] <- exp(predict(ps.logit,subset(lalonde,treat==0)))


###################################################
### code chunk number 28: twang.rnw:888-895
###################################################
bal.logit <- dx.wts(x = lalonde$w.logit,
                    data=lalonde,
                    vars=c("age","educ","black","hispan","nodegree",
                      "married","re74","re75"),
                    treat.var="treat",
                    perm.test.iters=0, estimand = "ATT")
bal.logit


###################################################
### code chunk number 29: twang.rnw:903-910
###################################################
bal.logit <- dx.wts(x = lalonde$w.logit,
                    data=lalonde,
                    vars=c("age","educ","black","hispan","nodegree",
                      "married","re74","re75"),
                    treat.var="treat",
                    perm.test.iters=0, estimand = "ATT")
bal.logit


###################################################
### code chunk number 30: twang.rnw:922-930
###################################################
pretty.tab <- bal.table(bal.logit)[[2]][,c("tx.mn","ct.mn","ks")]
pretty.tab <- cbind(pretty.tab, bal.table(bal.logit)[[1]]$ct.mn)
names(pretty.tab) <- c("E(Y1|t=1)","E(Y0|t=1)","KS","E(Y0|t=0)")
xtable(pretty.tab,
       caption = "Logistic regression estimates of the propensity scores",
       label = "tab:balancelogit",
       digits = c(0, 2, 2, 2, 2),
       align=c("l","r","r","r","r"))


###################################################
### code chunk number 31: twang.rnw:940-957
###################################################
bal.gbm <- dx.wts(ps.lalonde,
                  data=lalonde, estimand = "ATE",
                  vars=c("age","educ","black","hispan","nodegree","married","re74","re75"),
                  treat.var="treat",
                  perm.test.iters=0)
pretty.tab <- rbind(bal.logit$summary.tab,
                    bal.gbm$summary.tab[-1,])
rownames(pretty.tab) <- pretty.tab$type
rownames(pretty.tab)[2] <- "logit"
pretty.tab <- pretty.tab[,c("n.treat","ess.ctrl","max.es","mean.es","max.ks","mean.ks")]
xtable(pretty.tab,
       caption = "Summary of the balancing properties of logistic regression and gbm",
       label = "tab:balancecompare",
       digits = c(0, 0, 2, 2, 2, 2, 2),
       #digits = c(0,0, rep(2,9)),
       align=c("l","r","r","r","r","r","r"))
       #align = rep("l", 11))


###################################################
### code chunk number 32: twang.rnw:960-963
###################################################
design.logit <- svydesign(ids=~1, weights=~w.logit, data=lalonde)
glm6 <- svyglm(re78 ~ treat, design=design.logit)
summary(glm6)


###################################################
### code chunk number 33: twang.rnw:976-979
###################################################
glm7 <- svyglm(re78 ~ treat + age + educ + black + hispan + nodegree +
                      married + re74 + re75,
               design=design.logit)


###################################################
### code chunk number 34: twang.rnw:1058-1061
###################################################
data(lindner)
table(lindner$sixMonthSurvive, lindner$abcix)
chisq.test(table(lindner$sixMonthSurvive, lindner$abcix))


###################################################
### code chunk number 35: twang.rnw:1072-1076
###################################################
set.seed(1)
ps.lindner <- ps(abcix ~ stent + height + female + diabetic + 
                 acutemi + ejecfrac + ves1proc, data = lindner,
                 verbose = FALSE, estimand = "ATE")


###################################################
### code chunk number 36: twang.rnw:1088-1089
###################################################
options(width=85)


###################################################
### code chunk number 37: twang.rnw:1092-1093
###################################################
bal.table(ps.lindner)


###################################################
### code chunk number 38: twang.rnw:1096-1097
###################################################
options(width = 80)


###################################################
### code chunk number 39: twang.rnw:1119-1120
###################################################
plot(ps.lindner, plots = 1)


###################################################
### code chunk number 40: twang.rnw:1127-1128
###################################################
plot(ps.lindner, plots = 2)


###################################################
### code chunk number 41: twang.rnw:1134-1135
###################################################
plot(ps.lindner, plots = 3)


###################################################
### code chunk number 42: twang.rnw:1142-1143
###################################################
plot(ps.lindner, plots = 4)


###################################################
### code chunk number 43: twang.rnw:1149-1150
###################################################
plot(ps.lindner, plots = 5)


###################################################
### code chunk number 44: twang.rnw:1163-1164
###################################################
summary(ps.lindner)


###################################################
### code chunk number 45: twang.rnw:1171-1174
###################################################
lindner$w <- get.weights(ps.lindner, stop.method = "es.mean")
design.ps <- svydesign(ids=~1, weights = ~w, data = lindner)
svychisq(~sixMonthSurvive + abcix, design = design.ps)


###################################################
### code chunk number 46: twang.rnw:1234-1235
###################################################
data(egsingle)


###################################################
### code chunk number 47: twang.rnw:1243-1244
###################################################
tmp <- tapply(egsingle$grade, egsingle$childid, unique) 


###################################################
### code chunk number 48: twang.rnw:1252-1253
###################################################
tmp <- lapply(tmp, function(x){return(x %in% 1:4)})


###################################################
### code chunk number 49: twang.rnw:1260-1261
###################################################
tmp <- lapply(tmp, sum)


###################################################
### code chunk number 50: twang.rnw:1266-1267
###################################################
tmp <- sapply(tmp, function(x){as.numeric(x == 4)})


###################################################
### code chunk number 51: twang.rnw:1272-1275
###################################################
tmp <- data.frame(tmp)
names(tmp) <- "resp"
tmp$childid <- row.names(tmp)


###################################################
### code chunk number 52: twang.rnw:1280-1281
###################################################
egsingle <- merge(egsingle, tmp)


###################################################
### code chunk number 53: twang.rnw:1288-1289
###################################################
egsingle.one <-unique(egsingle[,-c(3:6)])


###################################################
### code chunk number 54: twang.rnw:1294-1296
###################################################
egsingle.one$race <- as.factor(race <- ifelse(egsingle.one$black==1, 1,
                                         ifelse(egsingle.one$hispanic==1, 2, 3)))


###################################################
### code chunk number 55: twang.rnw:1308-1314
###################################################
egsingle.ps <-  ps(resp ~ race + female + size + lowinc + mobility,
      data=egsingle.one,
      stop.method=c("es.mean","ks.max"),
      n.trees=5000,
      verbose=FALSE,
      estimand = "ATE")


###################################################
### code chunk number 56: twang.rnw:1325-1326
###################################################
plot(egsingle.ps)


###################################################
### code chunk number 57: twang.rnw:1364-1365
###################################################
egsingle.one$wgt <- get.weights(egsingle.ps, stop.method="ks.max")


###################################################
### code chunk number 58: twang.rnw:1374-1376
###################################################
egtmp <- rbind(data.frame(egsingle.one, nr2=1, wgt2=1), 
                data.frame(egsingle.one, nr2=0, wgt2=egsingle.one$wgt)[egsingle.one$resp==1,])


###################################################
### code chunk number 59: twang.rnw:1382-1387
###################################################
egdxwts <- dx.wts(x=egtmp$wgt2, 
                   data=egtmp, 
                   estimand="ATT", 
                   vars=c("race", "female", "size",  "lowinc",  "mobility"),
                   treat.var="nr2")


###################################################
### code chunk number 60: twang.rnw:1390-1397
###################################################
pretty.tab<-bal.table(egdxwts)[[2]][,c("tx.mn","ct.mn","std.eff.sz","ks")]
names(pretty.tab) <- c("OverallS Sample","Weighted responders","Std ES","KS")
xtable(pretty.tab,
        caption = "Balance of the nonrespondents and respondents",
        label = "tab:balance2",
        digits = c(0, 2, 2, 2, 2),
        align=c("l","r","r","r","r"))


###################################################
### code chunk number 61: twang.rnw:1406-1409
###################################################
egsinge.resp <- merge(subset(egsingle, subset=resp==1),
                        subset(egsingle.one, subset=resp==1,
                               select=c(childid, wgt)) )



### R code from vignette source 'twang.rnw'

###################################################
### code chunk number 1: twang.rnw:100-101
###################################################
options(width=60)


###################################################
### code chunk number 2: twang.rnw:120-122
###################################################
library(twang)
set.seed(1)


###################################################
### code chunk number 3: twang.rnw:130-131
###################################################
data(lalonde)


###################################################
### code chunk number 4: twang.rnw:154-164
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
### code chunk number 7: twang.rnw:361-364
###################################################
summary(ps.lalonde$gbm.obj,
        n.trees=ps.lalonde$desc$ks.max.ATT$n.trees,
        plot=FALSE)


###################################################
### code chunk number 8: twang.rnw:370-372
###################################################
summary(ps.lalonde$gbm.obj,
        n.trees=ps.lalonde$desc$ks.max.ATT$n.trees)


###################################################
### code chunk number 9: twang.rnw:401-402
###################################################
options(width=85)


###################################################
### code chunk number 10: twang.rnw:405-407
###################################################
lalonde.balance <- bal.table(ps.lalonde)
lalonde.balance


###################################################
### code chunk number 11: twang.rnw:410-411
###################################################
options(width=60)


###################################################
### code chunk number 12: twang.rnw:484-493
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
### code chunk number 13: twang.rnw:507-508
###################################################
summary(ps.lalonde)


###################################################
### code chunk number 14: twang.rnw:584-585
###################################################
plot(ps.lalonde, plots=2)


###################################################
### code chunk number 15: twang.rnw:630-631
###################################################
plot(ps.lalonde, plots=3)


###################################################
### code chunk number 16: twang.rnw:645-646
###################################################
plot(ps.lalonde, plots = 4)


###################################################
### code chunk number 17: twang.rnw:661-662
###################################################
plot(ps.lalonde, plots = 5)


###################################################
### code chunk number 18: twang.rnw:697-698
###################################################
plot(ps.lalonde, plots = 3, subset = 2)


###################################################
### code chunk number 19: twang.rnw:714-715
###################################################
library(survey)


###################################################
### code chunk number 20: twang.rnw:727-729
###################################################
lalonde$w <- get.weights(ps.lalonde, stop.method="es.mean")
design.ps <- svydesign(ids=~1, weights=~w, data=lalonde)


###################################################
### code chunk number 21: twang.rnw:732-732
###################################################



###################################################
### code chunk number 22: twang.rnw:765-767
###################################################
glm1 <- svyglm(re78 ~ treat, design=design.ps)
summary(glm1)


###################################################
### code chunk number 23: twang.rnw:798-800
###################################################
glm2 <- svyglm(re78 ~ treat + nodegree, design=design.ps)
summary(glm2)


###################################################
### code chunk number 24: twang.rnw:812-816
###################################################
glm3 <- svyglm(re78 ~ treat + age + educ + black + hispan + nodegree +
                      married + re74 + re75,
               design=design.ps)
summary(glm3)


###################################################
### code chunk number 25: twang.rnw:825-829
###################################################
glm4 <- lm(re78 ~ treat + age + educ + black + hispan + nodegree +
                  married + re74 + re75,
           data=lalonde)
summary(glm4)


###################################################
### code chunk number 26: twang.rnw:832-835
###################################################
glm5 <- lm(sqrt(re78) ~ treat + age + educ + black + hispan + nodegree +
                        married + sqrt(re74) + sqrt(re75),
           data=lalonde)


###################################################
### code chunk number 27: twang.rnw:860-866
###################################################
ps.logit <- glm(treat ~ age + educ + black + hispan + nodegree +
                        married + re74 + re75,
                data = lalonde,
                family = binomial)
lalonde$w.logit <- rep(1,nrow(lalonde))
lalonde$w.logit[lalonde$treat==0] <- exp(predict(ps.logit,subset(lalonde,treat==0)))


###################################################
### code chunk number 28: twang.rnw:880-887
###################################################
bal.logit <- dx.wts(x = lalonde$w.logit,
                    data=lalonde,
                    vars=c("age","educ","black","hispan","nodegree",
                      "married","re74","re75"),
                    treat.var="treat",
                    perm.test.iters=0, estimand = "ATT")
bal.logit


###################################################
### code chunk number 29: twang.rnw:895-898
###################################################

bal.tab.logit <- bal.table(bal.logit)
bal.tab.logit


###################################################
### code chunk number 30: twang.rnw:910-918
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
### code chunk number 31: twang.rnw:928-945
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
### code chunk number 32: twang.rnw:948-951
###################################################
design.logit <- svydesign(ids=~1, weights=~w.logit, data=lalonde)
glm6 <- svyglm(re78 ~ treat, design=design.logit)
summary(glm6)


###################################################
### code chunk number 33: twang.rnw:961-964
###################################################
glm7 <- svyglm(re78 ~ treat + age + educ + black + hispan + nodegree +
                      married + re74 + re75,
               design=design.logit)


###################################################
### code chunk number 34: twang.rnw:1044-1047
###################################################
data(lindner)
table(lindner$sixMonthSurvive, lindner$abcix)
chisq.test(table(lindner$sixMonthSurvive, lindner$abcix))


###################################################
### code chunk number 35: twang.rnw:1057-1061
###################################################
set.seed(1)
ps.lindner <- ps(abcix ~ stent + height + female + diabetic + 
                 acutemi + ejecfrac + ves1proc, data = lindner,
                 verbose = FALSE, estimand = "ATE")


###################################################
### code chunk number 36: twang.rnw:1071-1072
###################################################
options(width=85)


###################################################
### code chunk number 37: twang.rnw:1075-1076
###################################################
bal.table(ps.lindner)


###################################################
### code chunk number 38: twang.rnw:1079-1080
###################################################
options(width = 60)


###################################################
### code chunk number 39: twang.rnw:1101-1102
###################################################
plot(ps.lindner, plots = 1)


###################################################
### code chunk number 40: twang.rnw:1110-1111
###################################################
plot(ps.lindner, plots = 2)


###################################################
### code chunk number 41: twang.rnw:1117-1118
###################################################
plot(ps.lindner, plots = 3)


###################################################
### code chunk number 42: twang.rnw:1125-1126
###################################################
plot(ps.lindner, plots = 4)


###################################################
### code chunk number 43: twang.rnw:1132-1133
###################################################
plot(ps.lindner, plots = 5)


###################################################
### code chunk number 44: twang.rnw:1145-1146
###################################################
summary(ps.lindner)


###################################################
### code chunk number 45: twang.rnw:1152-1155
###################################################
lindner$w <- get.weights(ps.lindner, stop.method = "es.mean")
design.ps <- svydesign(ids=~1, weights = ~w, data = lindner)
svychisq(~sixMonthSurvive + abcix, design = design.ps)


###################################################
### code chunk number 46: twang.rnw:1238-1239
###################################################
data(egsingle)


###################################################
### code chunk number 47: twang.rnw:1245-1247
###################################################
tmp <- sapply(split(egsingle,egsingle$childid),function(x){
              paste(as.character(x$grade),collapse="")})


###################################################
### code chunk number 48: twang.rnw:1251-1254
###################################################
tmp <- data.frame(childid=names(tmp), gpatt=tmp,
                  resp=as.numeric((1:length(tmp)) %in%
                          grep("1234",as.character(tmp))))


###################################################
### code chunk number 49: twang.rnw:1258-1259
###################################################
egsingle <- merge(egsingle, tmp)


###################################################
### code chunk number 50: twang.rnw:1265-1266
###################################################
egsingle.one <-unique(egsingle[,-c(3:6)])


###################################################
### code chunk number 51: twang.rnw:1271-1273
###################################################
egsingle.one$race <- as.factor(race <- ifelse(egsingle.one$black==1, 1,
                                         ifelse(egsingle.one$hispanic==1, 2, 3)))


###################################################
### code chunk number 52: twang.rnw:1285-1292
###################################################
egsingle.ps <-
   ps(resp ~ race + female + size + lowinc + mobility,
      data=egsingle.one,
      stop.method=c("es.mean","ks.max"),
      n.trees=2500,
      verbose=FALSE,
      estimand = "ATE")


###################################################
### code chunk number 53: twang.rnw:1300-1301
###################################################
plot(egsingle.ps)


###################################################
### code chunk number 54: twang.rnw:1325-1332
###################################################
pretty.tab<-bal.table(egsingle.ps)$ks.max.ATE[,c("tx.mn","ct.mn","std.eff.sz","ks")]
names(pretty.tab) <- c("Non-responders","Weighted responders","Std ES","KS")
xtable(pretty.tab,
       caption = "Balance of the nonrespondents and respondents",
       label = "tab:balance2",
       digits = c(0, 2, 2, 2, 2),
       align=c("l","r","r","r","r"))


###################################################
### code chunk number 55: twang.rnw:1337-1338
###################################################
egsingle.one$wgt <- get.weights(egsingle.ps, stop.method="ks.max")


###################################################
### code chunk number 56: twang.rnw:1343-1346
###################################################
 egsinge.resp <- merge(subset(egsingle, subset=resp==1),
                       subset(egsingle.one, subset=resp==1,
                              select=c(childid, wgt)) )



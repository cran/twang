### R code from vignette source 'twang.rnw'

###################################################
### code chunk number 1: twang.rnw:114-115
###################################################
options(width=80)


###################################################
### code chunk number 2: twang.rnw:126-128
###################################################
library(twang)
data(lalonde)


###################################################
### code chunk number 3: twang.rnw:137-150
###################################################
ps.lalonde.gbm = ps(treat ~ age + educ + black + hispan + nodegree +
                         married + re74 + re75,
                 data = lalonde,
                 n.trees=5000,
                 interaction.depth=2,
                 shrinkage=0.01,
                 estimand = "ATT",
                 stop.method=c("es.mean","ks.max"),
                 n.minobsinnode = 10,
                 n.keep = 1,
                 n.grid = 25,
                 ks.exact = NULL,
                 verbose=FALSE)


###################################################
### code chunk number 4: iterPt
###################################################
    plot(ps.lalonde.gbm)


###################################################
### code chunk number 5: twang.rnw:197-198
###################################################
options(width=85)


###################################################
### code chunk number 6: twang.rnw:201-203
###################################################
lalonde.balance <- bal.table(ps.lalonde.gbm)
lalonde.balance


###################################################
### code chunk number 7: twang.rnw:206-207
###################################################
options(width=80)


###################################################
### code chunk number 8: twang.rnw:260-261
###################################################
summary(ps.lalonde.gbm)


###################################################
### code chunk number 9: twang.rnw:311-312
###################################################
plot(ps.lalonde.gbm, plots=2)


###################################################
### code chunk number 10: twang.rnw:319-320
###################################################
plot(ps.lalonde.gbm, plots=3)


###################################################
### code chunk number 11: twang.rnw:330-331
###################################################
plot(ps.lalonde.gbm, plots = 4)


###################################################
### code chunk number 12: twang.rnw:340-341
###################################################
plot(ps.lalonde.gbm, plots = 5)


###################################################
### code chunk number 13: twang.rnw:374-375
###################################################
plot(ps.lalonde.gbm, plots = 3, subset = 2)


###################################################
### code chunk number 14: twang.rnw:383-386
###################################################
summary(ps.lalonde.gbm$gbm.obj,
        n.trees=ps.lalonde.gbm$desc$ks.max.ATT$n.trees,
        plot=FALSE)


###################################################
### code chunk number 15: twang.rnw:391-393
###################################################
summary(ps.lalonde.gbm$gbm.obj,
        n.trees=ps.lalonde.gbm$desc$ks.max.ATT$n.trees)


###################################################
### code chunk number 16: twang.rnw:408-409
###################################################
library(survey)


###################################################
### code chunk number 17: twang.rnw:414-416
###################################################
lalonde$w <- get.weights(ps.lalonde.gbm, stop.method="es.mean")
design.ps <- svydesign(ids=~1, weights=~w, data=lalonde)


###################################################
### code chunk number 18: twang.rnw:419-419
###################################################



###################################################
### code chunk number 19: twang.rnw:437-439
###################################################
glm1 <- svyglm(re78 ~ treat, design=design.ps)
summary(glm1)


###################################################
### code chunk number 20: twang.rnw:449-451
###################################################
glm2 <- svyglm(re78 ~ treat + nodegree, design=design.ps)
summary(glm2)


###################################################
### code chunk number 21: twang.rnw:458-462
###################################################
glm3 <- svyglm(re78 ~ treat + age + educ + black + hispan + nodegree +
                      married + re74 + re75,
               design=design.ps)
summary(glm3)


###################################################
### code chunk number 22: twang.rnw:469-473
###################################################
glm4 <- lm(re78 ~ treat + age + educ + black + hispan + nodegree +
                  married + re74 + re75,
           data=lalonde)
summary(glm4)


###################################################
### code chunk number 23: twang.rnw:476-479
###################################################
glm5 <- lm(sqrt(re78) ~ treat + age + educ + black + hispan + nodegree +
                        married + sqrt(re74) + sqrt(re75),
           data=lalonde)


###################################################
### code chunk number 24: twang.rnw:490-496
###################################################
ps.logit <- glm(treat ~ age + educ + black + hispan + nodegree +
                        married + re74 + re75,
                data = lalonde,
                family = binomial)
lalonde$w.logit <- rep(1,nrow(lalonde))
lalonde$w.logit[lalonde$treat==0] <- exp(predict(ps.logit,subset(lalonde,treat==0)))


###################################################
### code chunk number 25: twang.rnw:502-509
###################################################
bal.logit <- dx.wts(x = lalonde$w.logit,
                    data=lalonde,
                    vars=c("age","educ","black","hispan","nodegree",
                      "married","re74","re75"),
                    treat.var="treat",
                    perm.test.iters=0, estimand = "ATT")
bal.logit


###################################################
### code chunk number 26: twang.rnw:515-516
###################################################
bal.table(bal.logit)


###################################################
### code chunk number 27: twang.rnw:577-580
###################################################
design.logit <- svydesign(ids=~1, weights=~w.logit, data=lalonde)
glm6 <- svyglm(re78 ~ treat, design=design.logit)
summary(glm6)


###################################################
### code chunk number 28: twang.rnw:618-621
###################################################
glm7 <- svyglm(re78 ~ treat + age + educ + black + hispan + nodegree +
                      married + re74 + re75,
               design=design.logit)


###################################################
### code chunk number 29: twang.rnw:661-664
###################################################
data(lindner)
table(lindner$sixMonthSurvive, lindner$abcix)
chisq.test(table(lindner$sixMonthSurvive, lindner$abcix))


###################################################
### code chunk number 30: twang.rnw:670-676
###################################################
set.seed(1)
ps.lindner <- ps(abcix ~ stent + height + female + diabetic + 
                 acutemi + ejecfrac + ves1proc,
                 data = lindner,
                 estimand = "ATE",
                 verbose = FALSE)


###################################################
### code chunk number 31: twang.rnw:682-683
###################################################
options(width=85)


###################################################
### code chunk number 32: twang.rnw:686-687
###################################################
bal.table(ps.lindner)


###################################################
### code chunk number 33: twang.rnw:690-691
###################################################
options(width = 80)


###################################################
### code chunk number 34: twang.rnw:700-701
###################################################
plot(ps.lindner, plots = 1)


###################################################
### code chunk number 35: twang.rnw:706-707
###################################################
plot(ps.lindner, plots = 2)


###################################################
### code chunk number 36: twang.rnw:713-714
###################################################
plot(ps.lindner, plots = 3)


###################################################
### code chunk number 37: twang.rnw:719-720
###################################################
plot(ps.lindner, plots = 4)


###################################################
### code chunk number 38: twang.rnw:726-727
###################################################
plot(ps.lindner, plots = 5)


###################################################
### code chunk number 39: twang.rnw:733-734
###################################################
summary(ps.lindner)


###################################################
### code chunk number 40: twang.rnw:739-742
###################################################
lindner$w <- get.weights(ps.lindner, stop.method = "es.mean")
design.ps <- svydesign(ids=~1, weights = ~w, data = lindner)
svychisq(~sixMonthSurvive + abcix, design = design.ps)


###################################################
### code chunk number 41: twang.rnw:798-799
###################################################
data(egsingle)


###################################################
### code chunk number 42: twang.rnw:807-808
###################################################
tmp <- tapply(egsingle$grade, egsingle$childid, unique)


###################################################
### code chunk number 43: twang.rnw:816-817
###################################################
tmp <- lapply(tmp, function(x){return(x %in% 1:4)})


###################################################
### code chunk number 44: twang.rnw:824-825
###################################################
tmp <- lapply(tmp, sum)


###################################################
### code chunk number 45: twang.rnw:830-831
###################################################
tmp <- sapply(tmp, function(x){as.numeric(x == 4)})


###################################################
### code chunk number 46: twang.rnw:836-839
###################################################
tmp <- data.frame(tmp)
names(tmp) <- "resp"
tmp$childid <- row.names(tmp)


###################################################
### code chunk number 47: twang.rnw:844-845
###################################################
egsingle <- merge(egsingle, tmp)


###################################################
### code chunk number 48: twang.rnw:852-853
###################################################
egsingle.one <-unique(egsingle[,-c(3:6)])


###################################################
### code chunk number 49: twang.rnw:858-860
###################################################
egsingle.one$race <- as.factor(race <- ifelse(egsingle.one$black==1, 1,
                                         ifelse(egsingle.one$hispanic==1, 2, 3)))


###################################################
### code chunk number 50: twang.rnw:872-878
###################################################
egsingle.ps <-  ps(resp ~ race + female + size + lowinc + mobility,
      data=egsingle.one,
      stop.method=c("es.mean","ks.max"),
      n.trees=5000,
      verbose=FALSE,
      estimand = "ATE")


###################################################
### code chunk number 51: twang.rnw:889-890
###################################################
plot(egsingle.ps)


###################################################
### code chunk number 52: twang.rnw:928-929
###################################################
egsingle.one$wgt <- get.weights(egsingle.ps, stop.method="ks.max")


###################################################
### code chunk number 53: twang.rnw:938-940
###################################################
egtmp <- rbind(data.frame(egsingle.one, nr2=1, wgt2=1),
                data.frame(egsingle.one, nr2=0, wgt2=egsingle.one$wgt)[egsingle.one$resp==1,])


###################################################
### code chunk number 54: twang.rnw:946-951
###################################################
egdxwts <- dx.wts(x=egtmp$wgt2,
                   data=egtmp,
                   estimand="ATT",
                   vars=c("race", "female", "size",  "lowinc",  "mobility"),
                   treat.var="nr2")


###################################################
### code chunk number 55: twang.rnw:954-962
###################################################
# pretty.tab<-bal.table(egdxwts)[[2]][,c("tx.mn","ct.mn","std.eff.sz","ks")]
# names(pretty.tab) <- c("OverallS Sample","Weighted responders","Std ES","KS")
# xtable(pretty.tab,
#         caption = "Balance of the nonrespondents and respondents",
#         label = "tab:balance2",
#         digits = c(0, 2, 2, 2, 2),
#         align=c("l","r","r","r","r"))
bal.table(egdxwts)[[2]]


###################################################
### code chunk number 56: twang.rnw:971-974
###################################################
egsinge.resp <- merge(subset(egsingle, subset=resp==1),
                        subset(egsingle.one, subset=resp==1,
                               select=c(childid, wgt)) )



### R code from vignette source 'twang.rnw'

###################################################
### code chunk number 1: twang.rnw:111-112
###################################################
options(width=80)


###################################################
### code chunk number 2: twang.rnw:123-125
###################################################
library(twang)
data(lalonde)


###################################################
### code chunk number 3: twang.rnw:134-147
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
### code chunk number 5: twang.rnw:194-195
###################################################
options(width=85)


###################################################
### code chunk number 6: twang.rnw:198-200
###################################################
lalonde.balance <- bal.table(ps.lalonde.gbm)
lalonde.balance


###################################################
### code chunk number 7: twang.rnw:203-204
###################################################
options(width=80)


###################################################
### code chunk number 8: twang.rnw:257-258
###################################################
summary(ps.lalonde.gbm)


###################################################
### code chunk number 9: twang.rnw:308-309
###################################################
plot(ps.lalonde.gbm, plots=2)


###################################################
### code chunk number 10: twang.rnw:316-317
###################################################
plot(ps.lalonde.gbm, plots=3)


###################################################
### code chunk number 11: twang.rnw:327-328
###################################################
plot(ps.lalonde.gbm, plots = 4)


###################################################
### code chunk number 12: twang.rnw:337-338
###################################################
plot(ps.lalonde.gbm, plots = 5)


###################################################
### code chunk number 13: twang.rnw:371-372
###################################################
plot(ps.lalonde.gbm, plots = 3, subset = 2)


###################################################
### code chunk number 14: twang.rnw:380-383
###################################################
summary(ps.lalonde.gbm$gbm.obj,
        n.trees=ps.lalonde.gbm$desc$ks.max.ATT$n.trees,
        plot=FALSE)


###################################################
### code chunk number 15: twang.rnw:388-390
###################################################
summary(ps.lalonde.gbm$gbm.obj,
        n.trees=ps.lalonde.gbm$desc$ks.max.ATT$n.trees)


###################################################
### code chunk number 16: twang.rnw:405-406
###################################################
library(survey)


###################################################
### code chunk number 17: twang.rnw:411-413
###################################################
lalonde$w <- get.weights(ps.lalonde.gbm, stop.method="es.mean")
design.ps <- svydesign(ids=~1, weights=~w, data=lalonde)


###################################################
### code chunk number 18: twang.rnw:416-416
###################################################



###################################################
### code chunk number 19: twang.rnw:434-436
###################################################
glm1 <- svyglm(re78 ~ treat, design=design.ps)
summary(glm1)


###################################################
### code chunk number 20: twang.rnw:446-448
###################################################
glm2 <- svyglm(re78 ~ treat + nodegree, design=design.ps)
summary(glm2)


###################################################
### code chunk number 21: twang.rnw:455-459
###################################################
glm3 <- svyglm(re78 ~ treat + age + educ + black + hispan + nodegree +
                      married + re74 + re75,
               design=design.ps)
summary(glm3)


###################################################
### code chunk number 22: twang.rnw:466-470
###################################################
glm4 <- lm(re78 ~ treat + age + educ + black + hispan + nodegree +
                  married + re74 + re75,
           data=lalonde)
summary(glm4)


###################################################
### code chunk number 23: twang.rnw:473-476
###################################################
glm5 <- lm(sqrt(re78) ~ treat + age + educ + black + hispan + nodegree +
                        married + sqrt(re74) + sqrt(re75),
           data=lalonde)


###################################################
### code chunk number 24: twang.rnw:486-487
###################################################
library(OVtool)


###################################################
### code chunk number 25: twang.rnw:492-504
###################################################
results = outcome_model(ps_object = NULL,
                        stop.method = NULL, 
                        data = lalonde,
                        weights="w",
                        treatment = "treat",
                        outcome = "re78", 
                        model_covariates = c("age", "educ",
                                             "black", "hispan",
                                             "nodegree", "married",
                                             "re74", "re75"),
                        estimand = "ATT")
summary(results$mod_results)


###################################################
### code chunk number 26: twang.rnw:509-517 (eval = FALSE)
###################################################
## ovtool_att = ov_sim(model_results=results,
##                                   plot_covariates=c("age", "educ", "hispan",
##                                                     "nodegree", "married",
##                                                     "re74", "re75"),
##                                   es_grid = NULL,
##                                   rho_grid = seq(0, 0.40, by = 0.05),
##                                   n_reps = 200,
##                                   progress = TRUE)


###################################################
### code chunk number 27: twang.rnw:520-522
###################################################
# save(ovtool_att, file = "data/ovtool_att.rda", version = 2)
load("ovtool_att.rda")


###################################################
### code chunk number 28: twang.rnw:526-527
###################################################
plot.ov(ovtool_att, col='color', print_graphic = 3)


###################################################
### code chunk number 29: twang.rnw:533-534
###################################################
summary.ov(ovtool_att, model_results = results)


###################################################
### code chunk number 30: twang.rnw:545-551
###################################################
ps.logit <- glm(treat ~ age + educ + black + hispan + nodegree +
                        married + re74 + re75,
                data = lalonde,
                family = binomial)
lalonde$w.logit <- rep(1,nrow(lalonde))
lalonde$w.logit[lalonde$treat==0] <- exp(predict(ps.logit,subset(lalonde,treat==0)))


###################################################
### code chunk number 31: twang.rnw:557-564
###################################################
bal.logit <- dx.wts(x = lalonde$w.logit,
                    data=lalonde,
                    vars=c("age","educ","black","hispan","nodegree",
                      "married","re74","re75"),
                    treat.var="treat",
                    perm.test.iters=0, estimand = "ATT")
bal.logit


###################################################
### code chunk number 32: twang.rnw:570-571
###################################################
bal.table(bal.logit)


###################################################
### code chunk number 33: twang.rnw:632-635
###################################################
design.logit <- svydesign(ids=~1, weights=~w.logit, data=lalonde)
glm6 <- svyglm(re78 ~ treat, design=design.logit)
summary(glm6)


###################################################
### code chunk number 34: twang.rnw:673-676
###################################################
glm7 <- svyglm(re78 ~ treat + age + educ + black + hispan + nodegree +
                      married + re74 + re75,
               design=design.logit)


###################################################
### code chunk number 35: twang.rnw:715-718
###################################################
data(lindner)
table(lindner$sixMonthSurvive, lindner$abcix)
chisq.test(table(lindner$sixMonthSurvive, lindner$abcix))


###################################################
### code chunk number 36: twang.rnw:724-730
###################################################
set.seed(1)
ps.lindner <- ps(abcix ~ stent + height + female + diabetic + 
                 acutemi + ejecfrac + ves1proc,
                 data = lindner,
                 estimand = "ATE",
                 verbose = FALSE)


###################################################
### code chunk number 37: twang.rnw:736-737
###################################################
options(width=85)


###################################################
### code chunk number 38: twang.rnw:740-741
###################################################
bal.table(ps.lindner)


###################################################
### code chunk number 39: twang.rnw:744-745
###################################################
options(width = 80)


###################################################
### code chunk number 40: twang.rnw:754-755
###################################################
plot(ps.lindner, plots = 1)


###################################################
### code chunk number 41: twang.rnw:760-761
###################################################
plot(ps.lindner, plots = 2)


###################################################
### code chunk number 42: twang.rnw:767-768
###################################################
plot(ps.lindner, plots = 3)


###################################################
### code chunk number 43: twang.rnw:773-774
###################################################
plot(ps.lindner, plots = 4)


###################################################
### code chunk number 44: twang.rnw:780-781
###################################################
plot(ps.lindner, plots = 5)


###################################################
### code chunk number 45: twang.rnw:787-788
###################################################
summary(ps.lindner)


###################################################
### code chunk number 46: twang.rnw:793-796
###################################################
lindner$w <- get.weights(ps.lindner, stop.method = "es.mean")
design.ps <- svydesign(ids=~1, weights = ~w, data = lindner)
svychisq(~sixMonthSurvive + abcix, design = design.ps)


###################################################
### code chunk number 47: twang.rnw:803-805
###################################################
glm_ATE = svyglm(sixMonthSurvive ~ abcix, design = design.ps)
summary(glm_ATE)


###################################################
### code chunk number 48: twang.rnw:810-816
###################################################
results_ATE = outcome_model(data = lindner,
                            weights="w",
                            treatment = "abcix",
                            outcome = "sixMonthSurvive",
                            model_covariates = 1,
                            estimand = "ATE")


###################################################
### code chunk number 49: twang.rnw:821-828 (eval = FALSE)
###################################################
## ovtool_ate = ov_sim(model_results=results_ATE,
##                     plot_covariates=c("stent", "height", "female",
##                                       "diabetic", "ejecfrac"),
##                     es_grid = NULL,
##                     rho_grid = NULL,
##                     n_reps = 250,
##                     progress = TRUE)


###################################################
### code chunk number 50: twang.rnw:831-833
###################################################
# save(ovtool_ate, file = "vignettes/ovtool_ate.rda", version = 2)
load("ovtool_ate.rda")


###################################################
### code chunk number 51: twang.rnw:837-838
###################################################
plot.ov(ovtool_ate, col='color', print_graphic = 3, p_contours = c(0.01, 0.05))


###################################################
### code chunk number 52: twang.rnw:844-845
###################################################
summary.ov(ovtool_ate, model_results = results_ATE)


###################################################
### code chunk number 53: twang.rnw:901-902
###################################################
data(egsingle)


###################################################
### code chunk number 54: twang.rnw:910-911
###################################################
tmp <- tapply(egsingle$grade, egsingle$childid, unique)


###################################################
### code chunk number 55: twang.rnw:919-920
###################################################
tmp <- lapply(tmp, function(x){return(x %in% 1:4)})


###################################################
### code chunk number 56: twang.rnw:927-928
###################################################
tmp <- lapply(tmp, sum)


###################################################
### code chunk number 57: twang.rnw:933-934
###################################################
tmp <- sapply(tmp, function(x){as.numeric(x == 4)})


###################################################
### code chunk number 58: twang.rnw:939-942
###################################################
tmp <- data.frame(tmp)
names(tmp) <- "resp"
tmp$childid <- row.names(tmp)


###################################################
### code chunk number 59: twang.rnw:947-948
###################################################
egsingle <- merge(egsingle, tmp)


###################################################
### code chunk number 60: twang.rnw:955-956
###################################################
egsingle.one <-unique(egsingle[,-c(3:6)])


###################################################
### code chunk number 61: twang.rnw:961-963
###################################################
egsingle.one$race <- as.factor(race <- ifelse(egsingle.one$black==1, 1,
                                         ifelse(egsingle.one$hispanic==1, 2, 3)))


###################################################
### code chunk number 62: twang.rnw:975-981
###################################################
egsingle.ps <-  ps(resp ~ race + female + size + lowinc + mobility,
      data=egsingle.one,
      stop.method=c("es.mean","ks.max"),
      n.trees=5000,
      verbose=FALSE,
      estimand = "ATE")


###################################################
### code chunk number 63: twang.rnw:992-993
###################################################
plot(egsingle.ps)


###################################################
### code chunk number 64: twang.rnw:1031-1032
###################################################
egsingle.one$wgt <- get.weights(egsingle.ps, stop.method="ks.max")


###################################################
### code chunk number 65: twang.rnw:1041-1043
###################################################
egtmp <- rbind(data.frame(egsingle.one, nr2=1, wgt2=1),
                data.frame(egsingle.one, nr2=0, wgt2=egsingle.one$wgt)[egsingle.one$resp==1,])


###################################################
### code chunk number 66: twang.rnw:1049-1054
###################################################
egdxwts <- dx.wts(x=egtmp$wgt2,
                   data=egtmp,
                   estimand="ATT",
                   vars=c("race", "female", "size",  "lowinc",  "mobility"),
                   treat.var="nr2")


###################################################
### code chunk number 67: twang.rnw:1057-1065
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
### code chunk number 68: twang.rnw:1074-1077
###################################################
egsinge.resp <- merge(subset(egsingle, subset=resp==1),
                        subset(egsingle.one, subset=resp==1,
                               select=c(childid, wgt)) )



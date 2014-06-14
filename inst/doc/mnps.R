### R code from vignette source 'mnps.rnw'

###################################################
### code chunk number 1: mnps.rnw:67-68
###################################################
options(width=60)


###################################################
### code chunk number 2: mnps.rnw:81-84
###################################################
library(twang)
data(AOD)
set.seed(1)


###################################################
### code chunk number 3: mnps.rnw:127-131
###################################################
mnps.AOD <- mnps(treat ~ illact + crimjust + subprob + subdep + white,
                 data = AOD, estimand = "ATE", verbose = FALSE, 
                 stop.method = c("es.mean", "ks.mean"), 
                 n.trees = 3000)


###################################################
### code chunk number 4: mnps.rnw:179-180
###################################################
    plot(mnps.AOD, plots = 1)


###################################################
### code chunk number 5: mnps.rnw:204-205
###################################################
    plot(mnps.AOD, plots = 2)


###################################################
### code chunk number 6: mnps.rnw:229-230
###################################################
options(width=120)


###################################################
### code chunk number 7: mnps.rnw:234-235
###################################################
    plot(mnps.AOD, plots = 3)


###################################################
### code chunk number 8: mnps.rnw:256-257
###################################################
options(width=120)


###################################################
### code chunk number 9: mnps.rnw:263-264
###################################################
plot(mnps.AOD, plots = 3, pairwiseMax = FALSE, figureRows = 3)


###################################################
### code chunk number 10: mnps.rnw:280-281
###################################################
    plot(mnps.AOD, plots = 4)


###################################################
### code chunk number 11: mnps.rnw:288-289
###################################################
options(width=60)


###################################################
### code chunk number 12: mnps.rnw:308-309
###################################################
options(width=85)


###################################################
### code chunk number 13: mnps.rnw:312-313
###################################################
means.table(mnps.AOD, stop.method = "es.mean", digits = 3)


###################################################
### code chunk number 14: mnps.rnw:336-337
###################################################
bal.table(mnps.AOD, digits = 2)


###################################################
### code chunk number 15: mnps.rnw:349-350
###################################################
bal.table(mnps.AOD, collapse.to = "covariate", digits = 4)


###################################################
### code chunk number 16: mnps.rnw:354-355
###################################################
bal.table(mnps.AOD, collapse.to = "stop.method", digits = 4)


###################################################
### code chunk number 17: mnps.rnw:363-364
###################################################
summary(mnps.AOD)


###################################################
### code chunk number 18: mnps.rnw:367-368
###################################################
options(width=60)


###################################################
### code chunk number 19: mnps.rnw:389-392
###################################################
require(survey)
AOD$w <- get.weights(mnps.AOD, stop.method = "es.mean")
design.mnps <- svydesign(ids=~1, weights=~w, data=AOD)


###################################################
### code chunk number 20: mnps.rnw:397-399
###################################################
glm1 <- svyglm(suf12 ~ as.factor(treat), design = design.mnps)
summary(glm1)


###################################################
### code chunk number 21: mnps.rnw:427-431
###################################################
mnps.AOD.ATT <- mnps(treat ~ illact + crimjust + subprob + subdep + white,
                 data = AOD, estimand = "ATT", treatATT = "community", 
                 verbose = FALSE, n.trees = 3000, 
                 stop.method = c("es.mean", "ks.mean"))


###################################################
### code chunk number 22: mnps.rnw:438-439
###################################################
    plot(mnps.AOD.ATT, plots = 1)


###################################################
### code chunk number 23: mnps.rnw:446-447
###################################################
    plot(mnps.AOD.ATT, plots = 3)


###################################################
### code chunk number 24: mnps.rnw:451-452
###################################################
options(width=85)


###################################################
### code chunk number 25: mnps.rnw:460-461
###################################################
means.table(mnps.AOD.ATT, digits = 3)


###################################################
### code chunk number 26: mnps.rnw:469-470
###################################################
bal.table(mnps.AOD.ATT, digits = 2)


###################################################
### code chunk number 27: mnps.rnw:473-474
###################################################
options(width=60)


###################################################
### code chunk number 28: mnps.rnw:479-482
###################################################
require(survey)
AOD$w.ATT <- get.weights(mnps.AOD.ATT, stop.method = "es.mean")
design.mnps.ATT <- svydesign(ids=~1, weights=~w.ATT, data=AOD)


###################################################
### code chunk number 29: mnps.rnw:486-488
###################################################
glm1 <- svyglm(suf12 ~ as.factor(treat), design = design.mnps.ATT)
summary(glm1)



### R code from vignette source 'mnps.rnw'

###################################################
### code chunk number 1: mnps.rnw:66-67
###################################################
options(width=60)


###################################################
### code chunk number 2: mnps.rnw:80-83
###################################################
library(twang)
data(AOD)
set.seed(1)


###################################################
### code chunk number 3: mnps.rnw:126-130
###################################################
mnps.AOD <- mnps(treat ~ illact + crimjust + subprob + subdep + white,
                 data = AOD, estimand = "ATE", verbose = FALSE, 
                 stop.method = c("es.mean", "ks.mean"), 
                 n.trees = 3000)


###################################################
### code chunk number 4: mnps.rnw:178-179
###################################################
    plot(mnps.AOD, plots = 1)


###################################################
### code chunk number 5: mnps.rnw:203-204
###################################################
    plot(mnps.AOD, plots = 2)


###################################################
### code chunk number 6: mnps.rnw:228-229
###################################################
options(width=120)


###################################################
### code chunk number 7: mnps.rnw:233-234
###################################################
    plot(mnps.AOD, plots = 3)


###################################################
### code chunk number 8: mnps.rnw:253-254
###################################################
options(width=120)


###################################################
### code chunk number 9: mnps.rnw:260-261
###################################################
plot(mnps.AOD, plots = 3, summaryFcn = NULL, figureRows = 3)


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
### code chunk number 14: mnps.rnw:329-330
###################################################
bal.table(mnps.AOD)


###################################################
### code chunk number 15: mnps.rnw:338-339
###################################################
summary(mnps.AOD)


###################################################
### code chunk number 16: mnps.rnw:342-343
###################################################
options(width=60)


###################################################
### code chunk number 17: mnps.rnw:364-367
###################################################
require(survey)
AOD$w <- get.weights(mnps.AOD, stop.method = "es.mean")
design.mnps <- svydesign(ids=~1, weights=~w, data=AOD)


###################################################
### code chunk number 18: mnps.rnw:372-374
###################################################
glm1 <- svyglm(suf12 ~ as.factor(treat), design = design.mnps)
summary(glm1)


###################################################
### code chunk number 19: mnps.rnw:402-405
###################################################
mnps.AOD.ATT <- mnps(treat ~ illact + crimjust + subprob + subdep + white,
                 data = AOD, estimand = "ATT", treatATT = "community", 
                 verbose = FALSE, n.trees = 3000)


###################################################
### code chunk number 20: mnps.rnw:412-413
###################################################
    plot(mnps.AOD.ATT, plots = 1)


###################################################
### code chunk number 21: mnps.rnw:420-421
###################################################
    plot(mnps.AOD.ATT, plots = 3)


###################################################
### code chunk number 22: mnps.rnw:425-426
###################################################
options(width=85)


###################################################
### code chunk number 23: mnps.rnw:434-435
###################################################
means.table(mnps.AOD.ATT, digits = 3)


###################################################
### code chunk number 24: mnps.rnw:443-444
###################################################
bal.table(mnps.AOD.ATT)


###################################################
### code chunk number 25: mnps.rnw:447-448
###################################################
options(width=60)


###################################################
### code chunk number 26: mnps.rnw:453-456
###################################################
require(survey)
AOD$w.ATT <- get.weights(mnps.AOD.ATT, stop.method = "es.mean")
design.mnps.ATT <- svydesign(ids=~1, weights=~w.ATT, data=AOD)


###################################################
### code chunk number 27: mnps.rnw:460-462
###################################################
glm1 <- svyglm(suf12 ~ as.factor(treat), design = design.mnps.ATT)
summary(glm1)



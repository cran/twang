### R code from vignette source 'iptw.rnw'

###################################################
### code chunk number 1: iptw.rnw:81-82
###################################################
options(width=80)


###################################################
### code chunk number 2: iptw.rnw:127-129
###################################################
library(twang)
data(iptwExWide)


###################################################
### code chunk number 3: iptw.rnw:148-158 (eval = FALSE)
###################################################
## iptw.Ex <- iptw(list(tx1 ~ use0 + gender + age, 
##                      tx2 ~ use1 + use0 + tx1 + gender + age, 
##                      tx3 ~ use2 + use1 + use0 + tx2 + tx1 + gender + age),
##                  timeInvariant ~ gender + age,
##                  data = iptwExWide, 
##                  cumulative = FALSE,
##                  priorTreatment = FALSE,
##                  verbose = FALSE, 
##                  stop.method = "es.max", 
##                  n.trees = 5000)


###################################################
### code chunk number 4: iptw.rnw:163-171
###################################################
iptw.Ex <- iptw(list(tx1 ~ use0, tx2 ~ use1, tx3 ~ use2),
                 timeInvariant ~ gender + age,
                 data = iptwExWide, 
                 cumulative = TRUE,
                 priorTreatment = TRUE,
                 verbose = FALSE, 
                 stop.method = "es.max", 
                 n.trees = 5000)


###################################################
### code chunk number 5: iptw.rnw:182-183
###################################################
    plot(iptw.Ex, plots = 1)


###################################################
### code chunk number 6: iptw.rnw:188-189
###################################################
summary(iptw.Ex)


###################################################
### code chunk number 7: iptw.rnw:192-193
###################################################
bal.table(iptw.Ex, timePeriods = 3)


###################################################
### code chunk number 8: iptw.rnw:200-201
###################################################
    plot(iptw.Ex, plots = 2)


###################################################
### code chunk number 9: iptw.rnw:209-210
###################################################
    plot(iptw.Ex, plots = 2, timePeriods = 2:3)


###################################################
### code chunk number 10: iptw.rnw:223-224
###################################################
    plot(iptw.Ex, plots = 3, color = FALSE)


###################################################
### code chunk number 11: iptw.rnw:233-234
###################################################
    plot(iptw.Ex, plots = 4)


###################################################
### code chunk number 12: iptw.rnw:239-240
###################################################
    plot(iptw.Ex, plots = 5)


###################################################
### code chunk number 13: iptw.rnw:250-259 (eval = FALSE)
###################################################
## data(mnIptwExWide)
## mniptw.Ex <- iptw(list(tx1 ~ use0, tx2 ~ use1, tx3 ~ use2),
##                  timeInvariant ~ gender + age,
##                  data = mnIptwExWide, 
##                  cumulative = TRUE,
##                  priorTreatment = TRUE,
##                  verbose = FALSE, 
##                  stop.method = "es.max", 
##                  n.trees = 5000)


###################################################
### code chunk number 14: iptw.rnw:278-279
###################################################
unstabWt1 <- get.weights.unstab(iptw.Ex)


###################################################
### code chunk number 15: iptw.rnw:285-291
###################################################
library(survey)
nTx <- with(iptwExWide, tx1 + tx2 + tx3)
outDatUnstab <- data.frame(outcome = iptwExWide$outcome,
                     nTx, 
                     wt = unstabWt1$es.max.ATE)
sv1unstab <- svydesign(~1, weights = ~wt, data = outDatUnstab)


###################################################
### code chunk number 16: iptw.rnw:296-299
###################################################
fitUnstab <- svyglm(outcome ~ nTx, sv1unstab)
coef(fitUnstab)
confint(fitUnstab)


###################################################
### code chunk number 17: iptw.rnw:305-314
###################################################
fitList <- list(glm(tx1 ~ 1, family = binomial, data = iptwExWide),
                glm(tx2 ~ tx1, family = binomial, data = iptwExWide), 
                glm(tx3 ~ tx1 + tx2, family = binomial, data = iptwExWide))
numWt <- get.weights.num(iptw.Ex, fitList)
stabWt1 <- unstabWt1 * numWt
outDatStab <- data.frame(outcome = iptwExWide$outcome,
                     nTx, 
                     wt = stabWt1$es.max.ATE)
sv1stab <- svydesign(~1, weights = ~wt, data = outDatStab)


###################################################
### code chunk number 18: iptw.rnw:319-322
###################################################
fitStab <- svyglm(outcome ~ nTx, sv1stab)
coef(fitStab)
confint(fitStab)


###################################################
### code chunk number 19: iptw.rnw:331-332
###################################################
confint(lm(iptwExWide$outcome ~ nTx))



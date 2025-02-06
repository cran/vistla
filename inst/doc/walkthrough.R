### R code from vignette source 'walkthrough.Rnw'

###################################################
### code chunk number 1: seed
###################################################
set.seed(1)


###################################################
### code chunk number 2: install (eval = FALSE)
###################################################
## install.packages("vistla")


###################################################
### code chunk number 3: load
###################################################
library(vistla)


###################################################
### code chunk number 4: ju_load
###################################################
data(junction)


###################################################
### code chunk number 5: ju_form
###################################################
head(junction)


###################################################
### code chunk number 6: vistla_run
###################################################
vistla(Y~.,data=junction)->vistla_result


###################################################
### code chunk number 7: vistla_print
###################################################
vistla_result


###################################################
### code chunk number 8: walkthrough.Rnw:102-103
###################################################
plot(vistla_result)


###################################################
### code chunk number 9: walkthrough.Rnw:115-117
###################################################
#Score cut-off
plot(prune(vistla_result,score=.45))


###################################################
### code chunk number 10: walkthrough.Rnw:119-121
###################################################
#Target list limit
plot(prune(vistla_result,targets=c("A2","B3")))


###################################################
### code chunk number 11: path
###################################################
path_to(vistla_result,"B3")
#Also report scores
path_to(vistla_result,"B3",detailed=TRUE)


###################################################
### code chunk number 12: walkthrough.Rnw:146-147
###################################################
data(cchain)


###################################################
### code chunk number 13: walkthrough.Rnw:154-155
###################################################
vistla(Y~.,data=cchain,estimator="kt")->vistla_kt


###################################################
### code chunk number 14: walkthrough.Rnw:160-161
###################################################
plot(vistla_kt)


###################################################
### code chunk number 15: walkthrough.Rnw:174-176
###################################################
vistla(T~.,data=chain,flow="intoup")->vistla_rev
plot(vistla_rev)


###################################################
### code chunk number 16: walkthrough.Rnw:181-182
###################################################
plot(vistla(T~.,data=chain,flow="fromdown"))


###################################################
### code chunk number 17: walkthrough.Rnw:191-193
###################################################
chain_d<-chain
chain_d$M3p<-chain$M3


###################################################
### code chunk number 18: walkthrough.Rnw:198-202
###################################################
set.seed(1)
path_to(vistla(Y~.,data=chain_d),"T")
set.seed(7)
path_to(vistla(Y~.,data=chain_d),"T")


###################################################
### code chunk number 19: walkthrough.Rnw:208-209
###################################################
vistla(Y~.,data=chain_d,ensemble=ensemble(100,resample=FALSE))->vistla_ens


###################################################
### code chunk number 20: walkthrough.Rnw:213-214
###################################################
print(vistla_ens)



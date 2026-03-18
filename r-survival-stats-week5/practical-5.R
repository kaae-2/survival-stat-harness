

# We consider again the Melanoma data, but this time in the version that is in the riskRegression R-package.

# Briefly, in the period 1962-77, 205 patients with malignant melanoma (cancer of the skin) had a radical operation performed at Odense University Hospital, Denmark. All patients were followed until the end of 1977 by which time 134 were still alive while 71 had died (of out whom 57 had died from cancer and 14 from other causes).

# The data has the following variables:

#    time time in days from operation
#    status a numeric with values ‘0=censored’ ‘1=death.malignant.melanoma’ ‘2=death.other.causes’
#    event a factor with levels ‘censored’ ‘death.malignant.melanoma’ ‘death.other.causes’
#    invasion a factor with levels ‘level.0’, ‘level.1’, ‘level.2’
#    ici inflammatory cell infiltration (IFI): 0, 1, 2 or 3
#    epicel a factor with levels ‘not present’ ‘present’
#    ulcer a factor with levels ‘not present’ ‘present’
#    thick tumour thickness (in 1/100 mm)
#    sex a factor with levels ‘Female’ ‘Male’
#    age age at operation (years)
#    logthick tumour thickness on log-scale

# The aim of the exercises are to familiarize you with how to estimate the cumulative incidence and how to regression for the cumulative incidence. We shall thus focus on the competing risks aspect of the data.

library(riskRegression)
library(prodlim)
library(mets)
library(cmprsk)
data(Melanoma,package="riskRegression")
Melanoma$logThick=log(Melanoma$thick)
Melanoma$ulcer <- factor(Melanoma$ulcer)
dtable(Melanoma,~status)

# We note that the status variable has the censorings scored as "0" which is what we want !

# head(Melanoma)
# table(Melanoma$status)
# with(Melanoma,table(status,event))

#  time status                    event invasion ici      epicel       ulcer
# 1   10      2       death.other.causes  level.1   2     present     present
# 2   30      2       death.other.causes  level.0   0 not present not present
# 3   35      0                 censored  level.1   2 not present not present
# 4   99      2       death.other.causes  level.0   2 not present not present
# 5  185      1 death.malignant.melanoma  level.2   2     present     present
# 6  204      1 death.malignant.melanoma  level.2   2 not present     present
#  thick    sex age   logthick   logThick
# 1  6.76   Male  76  1.9110229  1.9110229
# 2  0.65   Male  56 -0.4307829 -0.4307829
# 3  1.34   Male  41  0.2926696  0.2926696
# 4  2.90 Female  71  1.0647107  1.0647107
# 5 12.08   Male  52  2.4915512  2.4915512
# 6  4.84   Male  28  1.5769147  1.5769147

#  0   1   2
# 134  57  14
#      event
# status censored death.malignant.melanoma death.other.causes
#     0      134                        0                  0
#     1        0                       57                  0
#     2        0                        0                 14

# R commands for today

# Some relevant R commands from the cheat-sheet
# 
# # Product limit
# 
# cifu=prodlim(Hist(time,status)~ulcer,data=Melanoma)
# plot(cifu)
# summary(cifu,times=1000)
# 
# # Fine-Gray
# 
# 
# fitFG1 <- cifregFG(Event(time, status) ~ ulcer, data = Melanoma, cause = 1)
# fitFG2 <- cifregFG(Event(time, status) ~ ulcer, data = Melanoma, cause = 2)
# 
# nd <- data.frame(
#   ulcer = factor(c("not present", "present"),
#                  levels = levels(Melanoma$ulcer))
# )
# 
# pfg1 <- predict(fitFG1, newdata = nd)
# pfg2 <- predict(fitFG2, newdata = nd)
# 
# par(mfrow = c(1,2))
# plot(pfg1)
# plot(pfg2)
# # OR-version of CIF regression
# 
# fitOR1=cifreg(Event(time,status)~ulcer,data=Melanoma,cause=1)
# fitOR2=cifreg(Event(time,status)~ulcer,data=Melanoma,cause=2)
# 
# # ## predicting based on data.frame
# nd <- data.frame(ulcer=c(1,2))
# por1 <- predict(fitOR1,nd)
# por2 <- predict(fitOR2,nd)
# 
# par(mfrow=c(1,2))
# plot(por1)
# plot(por2)
# 
# # Based on cause-specific hazards modelling
# 
# library(riskRegression); library(pec)
# fit2 <- CSC(Hist(time,status) ~ ulcer, data=Melanoma)
# 
# nd <- data.frame(
#   ulcer = factor(c("not present","present"),
#                  levels = levels(Melanoma$ulcer))
# )
# 
# times <- seq(10, 4000, by = 100)
# par(mfrow=c(1,1))
# plotPredictEventProb(fit2, newdata = nd, times = times, cause = 1, col = 1:2)
# 
# # If you can install the package on your own, requires some additional work that you may skip, but is my preferred method. GOF (cumulative residuals) for "proportionality"
# 
# # library(crskdiag)
# # out1 <- diag_crr(Crsk(time,status)~ulcer,data=Melanoma,test="prop",seed=1234)
# # print(out1)
# source("https://www.biostat.ku.dk/ts/survival/data/gofFG.R")
# gg <- gofFG(Event(time,status)~ulcer,data=Melanoma)
# summary(gg)

# and for OR-link model

# gofOR1=prop.odds.subdist(Event(time,status)~ulcer,data=Melanoma,cause=1)
# summary(gofOR1)

# We could also do graphical test for proportionality in the two models

# fitFG1=cifreg(Event(time,status)~strata(ulcer),data=Melanoma,cause=1,prop=NULL)
# plot(fitFG1,log="xy")

# Q1. Non-parametric Cumulative Incidence
# Estimate the cumulative incidence for the two causes, do plots and list the estimates with 95% confidence intervals at 2000 days.
fit.ci=prodlim(Hist(time,status)~1,data=Melanoma)

# plots
par(mfrow=c(1,2))
plot(fit.ci, cause = 1, xlab = "Days", ylab = "Cumulative incidence",
     main = "Cause 1: death from melanoma")
plot(fit.ci, cause = 2, xlab = "Days", ylab = "Cumulative incidence",
     main = "Cause 2: death from other causes")

# estimates with 95% CI at 2000 days
ci.cause1 <- summary(fit.ci, times = 2000, cause = 1)
ci.cause2 <- summary(fit.ci, times = 2000, cause = 2)

ci.cause1
ci.cause2

# Q2 Stratified after invasion
# Plot the Aalen Johansen estimates of the cumulative incidence functions 
# stratified by the levels of the covariate invasion. 
# Do a Gray's test for equivalence for the cumulative incidence curves.

# make invasion a factor for grouped curves
Melanoma$invasion <- factor(Melanoma$invasion)

# Aalen-Johansen cumulative incidence curves stratified by invasion
fit.inv <- prodlim(Hist(time, status) ~ invasion, data = Melanoma)

par(mfrow = c(1,2))
plot(fit.inv, cause = 1,
     xlab = "Days", ylab = "Cumulative incidence",
     main = "Cause 1 by invasion")

plot(fit.inv, cause = 2,
     xlab = "Days", ylab = "Cumulative incidence",
     main = "Cause 2 by invasion")

# Gray's test
gray.inv <- cuminc(
  ftime = Melanoma$time,
  fstatus = Melanoma$status,
  group = Melanoma$invasion,
  cencode = 0
)

gray.inv


# Q3 CIF regression, FG-Model

#    Fit a Fine-Gray Model cumulative incidence model for the covariate invasion for the two cumulative incidence functions.
#    Make predictions and compare to the non-parametric estimates.
#    Check the models in terms of the "proportionality"

# R-Hint: to get the data frame with the invasion levels for the predictions

nd=data.frame(invasion=levels(Melanoma$invasion))

fitFG1=cifreg(Event(time,status)~invasion,data=Melanoma,prop=NULL,cause=1)
fitFG2=cifreg(Event(time,status)~invasion,data=Melanoma,prop=NULL,cause=2)

nd=data.frame(invasion=levels(Melanoma$invasion))

par(mfrow=c(1,2))
pfit1=predict(fitFG1,nd)
pfit2=predict(fitFG2,nd)

plot(pfit1,col=1:3)
plot(fit1,add=TRUE,confint=FALSE)
legend("topleft",legend=c("level.0","level.1","level.2"),lty=1:3,col=1:3)
plot(pfit2)
plot(fit1,cause=2,add=TRUE,confint=FALSE)
legend("topleft",legend=c("level.0","level.1","level.2"),lty=1:3,col=1:3)

summary(fitFG1)
summary(fitFG2)

out1 <- gofFG(Event(time,status)~invasion,data=Melanoma)
summary(out1)
# Q4 CIF regression, logistic-link OR-Model
# Fit the logistic-link OR cumulative incidence model for the covariate invasion for the two cumulative incidence functions. 
# Make predictions and compare to the non-parametric estimates.

fitOR1=cifreg(Event(time,status)~invasion,data=Melanoma,cause=1)
fitOR2=cifreg(Event(time,status)~invasion,data=Melanoma,cause=2)

nd=data.frame(invasion=levels(Melanoma$invasion))

par(mfrow=c(1,2))
pfit1=predict(fitOR1,nd)
pfit2=predict(fitOR2,nd)

plot(pfit1,col=1:3)
legend("topleft",legend=c("level.0","level.1","level.2"),lty=1:3,col=1:3)
plot(pfit2)
legend("topleft",legend=c("level.0","level.1","level.2"),lty=1:3,col=1:3)

summary(fitOR1)
summary(fitOR2)

library(timereg)
gofOR1=prop.odds.subdist(Event(time,status)~invasion,data=Melanoma,cause=1)
gofOR2=prop.odds.subdist(Event(time,status)~invasion,data=Melanoma,cause=2)
summary(gofOR1)
summary(gofOR2)


# Q4G
# If you have time you may also do a check that the censoring distribution does not depend on the covariates.
cc <-coxph(Surv(time,status==0)~invasion,data=Melanoma)
summary(cc)

# Q5 Fixed time-horizon IPCW analysis
# Still looking only at invasion-level make a competing risks analysis that describes 
# the risk of dying from malignant melanoma before 2000 days (1000 days) using the binomial regression binreg-function.
cc2000 <-binreg(Event(time,status)~invasion,data=Melanoma,cause=1,time=2000)
cc1000 <-binreg(Event(time,status)~invasion,data=Melanoma,cause=1,time=1000)
summary(cc2000)
summary(cc1000)

# Q6 Cause Specific Hazards
# If time allows you might also compare to the predictions obtained starting instead with the cause-specific hazards. But this may be skipped.
# Fit a univariate model with the covariate invasion for the 2 cause specific hazard functions and plot the cumulative incidence functions for the 3 levels of invasion based on the cause-specific hazards.
#
library(pec)
fit2 <- CSC(Hist(time,status)~invasion,data=Melanoma)
par(mfrow=c(1,2))
nd=data.frame(invasion=levels(Melanoma$invasion))
times=seq(0,5000,by=100)
plotPredictEventProb(fit2,newdata=nd,times,cause=1,col=1:3)
legend("topleft",legend=c("level.0","level.1","level.2"),lty=1,col=1:3)
plotPredictEventProb(fit2,newdata=nd,times,cause=2,col=1:3)
legend("topleft",legend=c("level.0","level.1","level.2"),lty=1,col=1:3)

fit2 <- CSC(Hist(time,status)~invasion,data=Melanoma)
print(fit2)

cause1=predictEventProb(fit2,newdata=data.frame(invasion="level.0"),cause=1,times=fit$time)
cause2=predictEventProb(fit2,newdata=data.frame(invasion="level.0"),cause=2,times=fit$time)
Melanoma$survStatus=ifelse(Melanoma$status>0,1,0)
fitSurv=coxph(Surv(time,survStatus)~invasion,data=Melanoma,x=TRUE)
survProb=predictSurvProb(fitSurv,newdata=data.frame(invasion="level.0"),times=fit$time)
plot(fit$time,cause1,type="s",ylim=c(0,1),col="red",ylab="CIF",xlab="Time In Days")
lines(fit$time,cause1+cause2,col="blue",type="s")
lines(fit$time,survProb,col="magenta",type="s")
legend("topright", legend=c("CIF Cause 1","CIF Cause 2+1","Survival Prob"),
       col=c("red","blue","magenta"),lty=1)

  
## CIF-2-Practicals
# Table of Contents

# We shall now consider multivariate cumulative incidence regression modelling.

library(riskRegression)
library(survival)
library(cmprsk)
library(pec)
library(mets)
data(Melanoma,package="riskRegression")
Melanoma$logThick=log(Melanoma$thick)

# To simplify the R-aspects of the below we code all factors with two levels as a numeric variable and work with some of these in various places below

dnumeric(Melanoma) <- ~.
dtable(Melanoma,~ulcer+ulcer.n)
dtable(Melanoma,~epicel+epicel.n)
dtable(Melanoma,~sex+sex.n)
dtable(Melanoma,~invasion+invasion.n)
dtable(Melanoma,~event+event.n)
Melanoma=transform(Melanoma,ulcer=ulcer.n,epicel=epicel.n,sex=sex.n)

# It is important that the censorings are coded as "0" for the below analyses. 
# This is asssumed by most programs, although a different code can also be specified in the programs.
# Q6 Multivariate Modelling

# Do a Fine-Gray model with the covariates: sex,epicel,ulcer,age,logThick and invasion.

# Investigate the fit of the model

#     Proportionality
#         computer experts may use the diagcrr function in the library(crskdiag)
#         the graphical test is an alternative, here you just look at sex and invasion
#     Interactions (ignore for now)
#     Linearity (ignore for now)
#     Censoring assumption (Q8 below)

# Do predictions for a subject with different invasion levels, and

#     sex=1,epicel=1,ulcer=1,age=50,logThick=1.


fitFG1=cifreg(Event(time,status)~sex+epicel+ulcer+age+logThick+invasion,data=Melanoma,cause=1,prop=NULL)
fitFG2=cifreg(Event(time,status)~sex+epicel+ulcer+age+logThick+invasion,data=Melanoma,cause=2,prop=NULL)

## newdata for predictions
nd=data.frame(invasion=levels(Melanoma$invasion),sex=1,epicel=1,ulcer=1,age=50,logThick=1)

par(mfrow=c(1,2))
pfit1=predict(fitFG1,nd)
pfit2=predict(fitFG2,nd)

plot(pfit1,col=1:3)
legend("topleft",legend=c("ADJ-level.0","ADJ-level.1","ADJ-level.2"),lty=1:3,col=1:3)
plot(pfit2)
legend("topleft",legend=c("ADJ-level.0","ADJ-level.1","ADJ-level.2"),lty=1:3,col=1:3)
summary(fitFG1)
summary(fitFG2)
exp(estimate(coef=fitFG1$coef,vcov=fitFG1$var)$coefmat[,c(1,3,4)])

s1fitFG1=cifreg(Event(time,status)~strata(sex)+epicel+ulcer+age+logThick+invasion,data=Melanoma,cause=1,propodds=NULL)
s2fitFG1=cifreg(Event(time,status)~sex+epicel+ulcer+age+logThick+strata(invasion),data=Melanoma,cause=1,propodds=NULL)

par(mfrow=c(1,2))
plot(s1fitFG1,log="y"); title(main="sex"); 
plot(s2fitFG1,log="y"); title(main="Invasion");
# Q7 Multivariate Modelling OR-modelling
# Do a logistic regression model for the CIF's. Again you may also get around 
# to checking the Goodness of fit of the model, and the predictions.

fitOR1=cifreg(Event(time,status)~sex+epicel+ulcer+age+logThick+invasion,data=Melanoma,cause=1)
fitOR2=cifreg(Event(time,status)~sex+epicel+ulcer+age+logThick+invasion,data=Melanoma,cause=2)

nd=data.frame(invasion=levels(Melanoma$invasion),sex=1,epicel=1,ulcer=1,age=50,logThick=1)

par(mfrow=c(1,2))
pfit1=predict(fitOR1,nd)
pfit2=predict(fitOR2,nd)

plot(pfit1,col=1:3)
legend("topleft",legend=c("ADJ-level.0","ADJ-level.1","ADJ-level.2"),lty=1:3,col=1:3)
plot(pfit2)
legend("topleft",legend=c("ADJ-level.0","ADJ-level.1","ADJ-level.2"),lty=1:3,col=1:3)

summary(fitOR1)
summary(fitOR2)

library(timereg)
gofOR1=prop.odds.subdist(Event(time,status)~sex+epicel+ulcer+age+logThick+invasion,data=Melanoma,cause=1)
gofOR2=prop.odds.subdist(Event(time,status)~sex+epicel+ulcer+age+logThick+invasion,data=Melanoma,cause=2)
summary(gofOR1)
summary(gofOR2)
 # Q8 Multivariate Modelling, censoring assumption
# Check if the censorings are independent of the covariates by doing a Cox model for the censorings. 
# Specifically we fit a Cox regression model for the censorings.

censModel=coxph(Surv(time,status==0)~sex+ulcer+epicel+
                  logThick+age,data=Melanoma)
summary(censModel)

# Q9. Binomial-Regression: Modelling just one time-point
# When there are model fit issues as we have seen above it can be useful to consider 
# just specific time-points such as the cumulative incidence at 2000 days, for example. 
# Then we may fit the binomial regression model with censoring adjustment.
# Do a binomial regression for cause 1 at 2000 days.
# Here we also rely on the assumption that the censorings do not depend on the covariates.

fitOR1t=binreg(Event(time,status)~sex+epicel+ulcer+age+logThick+invasion,data=Melanoma,cause=1,time=2000)
summary(fitOR1t)
dcut(Melanoma,breaks=2) <- ~logThick+age
fitOR1tc=binreg(Event(time,status)~sex+epicel+ulcer+age+logThick+invasion,data=Melanoma,
                cause=1,time=2000,cens.model=~strata(logThickcat.2,agecat.2))
summary(fitOR1tc)


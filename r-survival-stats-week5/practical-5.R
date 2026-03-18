

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

# library(riskRegression)
# library(prodlim)
# library(mets)
# library(cmprsk)
# data(Melanoma,package="riskRegression")
# Melanoma$logThick=log(Melanoma$thick)
# dtable(Melanoma,~status)

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

# Product limit

# cifu=prodlim(Hist(time,status)~factor(ulcer),data=Melanoma)
# plot(cifu)
# summary(cifu,times=1000)

# Fine-Gray

# fitFG1=cifreg(Event(time,status)~factor(ulcer),data=Melanoma,cause=1,prop=NULL)
# fitFG2=cifreg(Event(time,status)~factor(ulcer),data=Melanoma,cause=2,prop=NULL)

# ## predicting based on data.frame
# nd <- data.frame(ulcer=c(1,2))
# pfg1 <- predict(fitFG1,nd)
# pfg2 <- predict(fitFG2,nd)

# ## plotting
# par(mfrow=c(1,2))
# plot(pfg1)
# plot(pfg2)

# OR-version of CIF regression

# fitOR1=cifreg(Event(time,status)~factor(ulcer),data=Melanoma,cause=1)
# fitOR2=cifreg(Event(time,status)~factor(ulcer),data=Melanoma,cause=2)

# ## predicting based on data.frame
# nd <- data.frame(ulcer=c(1,2))
# por1 <- predict(fitOR1,nd)
# por2 <- predict(fitOR2,nd)

# par(mfrow=c(1,2))
# plot(por1)
# plot(por2)

# Based on cause-specific hazards modelling

# library(riskRegression); library(pec)
# fit2 <- CSC(Hist(time,status)~ulcer,data=Melanoma)
# nd=data.frame(ulcer=c(1,2))
# times <- seq(10,4000,by=100)
# plotPredictEventProb(fit2,newdata=nd,times,cause=1,col=1:3)

# If you can install the package on your own, requires some additional work that you may skip, but is my preferred method. GOF (cumulative residuals) for "proportionality"

# #library(crskdiag)
# #out1 <- diag_crr(Crsk(time,status)~ulcer,data=Melanoma,test="prop",seed=1234)
# #print(out1)
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
# Q2 Stratified after invasion
# Plot the Aalen Johansen estimates of the cumulative incidence functions stratified by the levels of the covariate invasion. Do a Gray's test for equivalence for the cumulative incidence curves.
# Q3 CIF regression, FG-Model

#    Fit a Fine-Gray Model cumulative incidence model for the covariate invasion for the two cumulative incidence functions.
#    Make predictions and compare to the non-parametric estimates.
#    Check the models in terms of the "proportionality"

# R-Hint: to get the data frame with the invasion levels for the predictions

# nd=data.frame(invasion=levels(Melanoma$invasion))

# Q4 CIF regression, logistic-link OR-Model
# Fit the logistic-link OR cumulative incidence model for the covariate invasion for the two cumulative incidence functions. Make predictions and compare to the non-parametric estimates.
# Q4G
# If you have time you may also do a check that the censoring distribution does not depend on the covariates.
# Q5 Fixed time-horizon IPCW analysis
# Still looking only at invasion-level make a competing risks analysis that describes the risk of dying from malignant melanoma before 2000 days (1000 days) using the binomial regression binreg-function.
# Q6 Cause Specific Hazards

# If time allows you might also compare to the predictions obtained starting instead with the cause-specific hazards. But this may be skipped.

# Fit a univariate model with the covariate invasion for the 2 cause specific hazard functions and plot the cumulative incidence functions for the 3 levels of invasion based on the cause-specific hazards.
#
# CIF-2-Practicals
# Table of Contents

# We shall now consider multivariate cumulative incidence regression modelling.

# library(riskRegression)
# library(survival)
# library(cmprsk)
# library(pec)
# library(mets)
# data(Melanoma,package="riskRegression")
# Melanoma$logThick=log(Melanoma$thick)

# To simplify the R-aspects of the below we code all factors with two levels as a numeric variable and work with some of these in various places below

# dnumeric(Melanoma) <- ~.
# dtable(Melanoma,~ulcer+ulcer.n)
# dtable(Melanoma,~epicel+epicel.n)
# dtable(Melanoma,~sex+sex.n)
# dtable(Melanoma,~invasion+invasion.n)
# dtable(Melanoma,~event+event.n)
# Melanoma=transform(Melanoma,ulcer=ulcer.n,epicel=epicel.n,sex=sex.n)

# It is important that the censorings are coded as "0" for the below analyses. This is asssumed by most programs, although a different code can also be specified in the programs.
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

# Q7 Multivariate Modelling OR-modelling
# Do a logistic regression model for the CIF's. Again you may also get around to checking the Goodness of fit of the model, and the predictions.
# Q8 Multivariate Modelling, censoring assumption
# Check if the censorings are independent of the covariates by doing a Cox model for the censorings. Specifically we fit a Cox regression model for the censorings.
# Q9. Binomial-Regression: Modelling just one time-point

# When there are model fit issues as we have seen above it can be useful to consider just specific time-points such as the cumulative incidence at 2000 days, for example. Then we may fit the binomial regression model with censoring adjustment.

# Do a binomial regression for cause 1 at 2000 days.

# Here we also rely on the assumption that the censorings do not depend on the covariates.

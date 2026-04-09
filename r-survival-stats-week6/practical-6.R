# Exercise on standardization
# We will look at an extract of data drawn from the Rotterdam tumour bank. The subset
# contains information on lymph node positive breast cancer patients aged 50 or below. We
# will compare relapse-free survival for patients who received chemotherapy to those that did
# not. The data are observational (no randomization). The majority of events are relapses, two
# patients were observed to die without relapse. The variables we have access to are
# Name Content
# time Days from surgery to relapse, death or last follow-up (whichever comes first)
# event 0: alive without relapse, 1: relapse or death
# age age at surgery (years)
# chemo chemotherapy (0: no, 1: yes)
# size tumour size (<=20mm, 20mm-50mm, >50 mm)
# log2nodes number of positive lymph nodes, log2-transformed
# To load the data into a data frame called bca in R write
bca <- read.csv("http://biostat.ku.dk/frank/data/bca-subset.csv")
bca$chemo <- factor(bca$chemo)
# 1. Plot the (event-free) survival curves for patients who did and did not receive chemotherapy,
# respectively. What is the survival difference at 2000 days after surgery?
library(survival)
library(prodlim)
pl <- prodlim(Surv(time,event)~chemo, data=bca)
plot(pl)
# Time
# Survival probability
# 0 1000 2000 3000 4000 5000 6000
# 0 % 50 % 100 %
# 116 65 35 24 18 13 7 3 1 0 0
# 461 382 284 226 167 109 66 28 14 3 1
# chemo
# 0
# 1
summary(pl, times=2000)
# ## chemo time n.risk n.event n.lost surv se.surv lower upper
# ## 1 0 2000 23 0 0 0.224 0.0394 0.147 0.302
# ## 2 1 2000 210 0 0 0.504 0.0235 0.458 0.550
# 1
# 2. Consider the covariates in the data set. Make a table with baseline information, i.e., a
# typical “Table 1”, with the three columns chemotherapy yes and no and total. Reflect on
# any differences between groups. You may want to use the function univariateTable from the
# package Publish
library(Publish)
# univariateTable(chemo~Q(age)+Q(log2nodes)+size, data=bca)
# ## Variable Level chemo = 0 (n=116) chemo = 1 (n=461) Total (n=577)
# ## 1 age median [iqr] 44 [39.8, 48.0] 43 [38, 47] 43 [38, 47]
# ## 2 log2nodes median [iqr] 2.6 [2.0, 3.3] 1 [0.0, 2.3] 1.6 [0.0, 2.6]
# ## 3 size <=20 22 (19.0) 190 (41.2) 212 (36.7)
# ## 4 20-50 74 (63.8) 205 (44.5) 279 (48.4)
# ## 5 >50 20 (17.2) 66 (14.3) 86 (14.9)
# ## p-value
# ## 1 0.07368
# ## 2 < 1e-04
# ## 3
# ## 4
# ## 5 < 1e-04
# 3. Now consider also tumour size. Predict the probability to be event-free after 2000 days for
# all combinations of chemotherapy and size groups (estimate 6 numbers)
# There are multiple ways to do this, for example, you can fit 6 Kaplan-Meier curves with prodlim
pla <- prodlim(Surv(time,event)~chemo+size, data=bca)
# and use summary() to extract the probabilities
summary(pla, times=2000)
# ## chemo size time n.risk n.event n.lost surv se.surv lower upper
# ## 1 0 20-50 2000 13 0 0 0.213 0.0481 0.118 0.307
# ## 2 0 <=20 2000 7 0 0 0.338 0.1035 0.135 0.541
# ## 3 0 >50 2000 3 0 0 0.150 0.0798 0.000 0.306
# ## 4 1 20-50 2000 92 0 0 0.464 0.0350 0.395 0.533
# ## 5 1 <=20 2000 102 0 0 0.619 0.0357 0.549 0.689
# ## 6 1 >50 2000 16 0 0 0.294 0.0573 0.182 0.407
# You can just use the numbers from above and skip to the next question. For those interested in R, another
# way to get the same numbers is to use the predict function
# ## chemo=1
nd1 <- data.frame(expand.grid(chemo=factor(1), size=c("<=20",">50","20-50")))
pred1 <-predict(pla, newdata=nd1, times=2000)
pred1

# ## $`chemo=1, size=<=20`
# ## [1] 0.6191715
# ##
# ## $`chemo=1, size=>50`
# ## [1] 0.2942761
# ##
# ## $`chemo=1, size=20-50`
# ## [1] 0.4641005
# ## chemo=0
nd0 <- data.frame(expand.grid(chemo=factor(0), size=c("<=20",">50","20-50")))
pred0 <- predict(pla, newdata=nd0, times=2000)
pred0
# 2
# ## $`chemo=0, size=<=20`
# ## [1] 0.3380682
# ##
# ## $`chemo=0, size=>50`
# ## [1] 0.15
# ##
# ## $`chemo=0, size=20-50`
# ## [1] 0.2125823
# You can also use a Cox model stratified for both chemo and size and predict the probabilities with the
# function predictRisk from the riskRegression package
library(riskRegression)
# ## riskRegression version 2025.09.17
# ##
# ## Attaching package: 'riskRegression'
# ## The following object is masked from 'package:mets':
# ##
# ## predictRisk
fit.cox <- coxph(Surv(time,event)~ strata(chemo)+strata(size),x=TRUE,data=bca)
# 1-predictRisk(fit.cox, newdata=nd1, time=2000)
# ## [,1]
# ## [1,] 0.6201991
# ## [2,] 0.2998375
# ## [3,] 0.4654216
# 1-predictRisk(fit.cox, newdata=nd0, time=2000)
# ## [,1]
# ## [1,] 0.3537614
# ## [2,] 0.1712884
# ## [3,] 0.2180205
# Since we are stratifying for both (all) covariates, we have six baseline functions. This is essentially the same
# as fitting six Kaplan-Meier estimators.
# 4. How many patients are there in each of the tumour size groups?
# Numbers per group and total sample size
ngr <- table(bca$size)
n <- sum(ngr)
# ngr
# ##
# ## <=20 >50 20-50
# ## 212 86 279
# n
# ## [1] 577
# 5. This step requires calculation, you can use R as a calculator For the chemotherapy group,
# calculate the weighted average of the predicted survival at 2000 days, where each tumour
# size group has weight according to it’s proportion in the data set. This can be achieved by
# multiplying the estimated survival for the group with the number of patients in the group
# and summing over all three chemo=1 estimates and finally dividing by the total sample size.
# You may do this “manually” by copy+paste (it can also be done directly in R). Redo the
# procedure for the no chemotherapy group. Look at the difference between the weighted
# averages. Compare this to the estimates and difference from question 1. Interpret.
weighted.mean(unlist(pred1), ngr)
weighted.mean(unlist(pred0), ngr)
# [1] 0.4957646
# [1] 0.2493604
# 3


# 6. Compare your the standardized estimates to those from the ate function of the riskRegression
# package
fit.cox1 <- coxph(Surv(time,event)~ strata(chemo)+strata(size),x=TRUE,data=bca)
ate1 <- ate(fit.cox1, data = bca, treatment = "chemo", times = 2000, se = TRUE)
# ## Input variables
# ## - Treatment : chemo (2 levels: "0" "1")
# ## - Event : event (cause: TRUE, censoring: 0)
# ## - Time [min;max] : time [38;6270]
# ## - at risk/time : 2000
# ## number in treatment 0 23
# ## number in treatment 1 210
# ##
# ## Estimation procedure
# ## - Estimator : G-formula
# ## - Uncertainty: Gaussian approximation
# ## where the variance is estimated via the influence function
# ##
# ## Processing
# ## - Prepare influence function: outcome done
# ## - Point estimation: done
# ## - Decomposition iid: done
# ## - Confidence intervals: done
summary(ate1)
# ## Average treatment effect
# ##
# ## - Treatment : chemo (2 levels: "0" "1")
# ## - Event : event (cause: TRUE, censoring: 0)
# ## - Time [min;max] : time [38;6270]
# ## - at risk/time : 2000
# ## number in treatment 0 23
# ## number in treatment 1 210
# ##
# ## Estimation procedure
# ## - Estimator : G-formula
# ## - Uncertainty: Gaussian approximation
# ## where the variance is estimated via the influence function
# ##
# ## Testing procedure
# ## - Null hypothesis : given two treatments (A,B) and a specific timepoint, equal risks
# ## - Confidence level : 0.95
# ##
# ## Results:
# ## - Difference in standardized risk (B-A) between time zero and 'time'
# ## reported on the scale [-1;1] (difference between two probabilities)
# ## (difference in average risks when treating all subjects with the experimental treatment (B),
# ## vs. treating all subjects with the reference treatment (A))
# ##
# 4
# ## time chemo=A risk(chemo=A) chemo=B risk(chemo=B) difference ci
# ## 2000 0 0.739 1 0.502 -0.237 [-0.33;-0.14]
# ## p.value
# ## 2.12e-06
# ##
# ## difference : estimated difference in standardized risks
# ## ci : pointwise confidence intervals
# ## p.value : (unadjusted) p-value
# 7. Use either the function ate (riskRegession) or survivalG (mets) to estimate the standardized
# survival difference (ATE) using age, log2nodes and size based on a Cox model. Interpret the
# result. Do you think that the assumption of no unmeasured confounders is fulfilled?
# We can estimate the ATE with the function ate from the package riskRegression
fit.cox2 <- coxph(Surv(time,event)~ chemo+age+size+log2nodes,x=TRUE,data=bca)
summary(ate(fit.cox2, data = bca, treatment = "chemo", times = 2000, se = TRUE))
# ## Input variables
# ## - Treatment : chemo (2 levels: "0" "1")
# ## - Event : event (cause: TRUE, censoring: 0)
# ## - Time [min;max] : time [38;6270]
# ## - at risk/time : 2000
# ## number in treatment 0 23
# ## number in treatment 1 210
# ##
# ## Estimation procedure
# ## - Estimator : G-formula
# ## - Uncertainty: Gaussian approximation
# ## where the variance is estimated via the influence function
# ##
# ## Processing
# ## - Prepare influence function: outcome done
# ## - Point estimation: done
# ## - Decomposition iid: done
# ## - Confidence intervals: done
# ## Average treatment effect
# ##
# ## - Treatment : chemo (2 levels: "0" "1")
# ## - Event : event (cause: TRUE, censoring: 0)
# ## - Time [min;max] : time [38;6270]
# ## - at risk/time : 2000
# ## number in treatment 0 23
# ## number in treatment 1 210
# ##
# ## Estimation procedure
# ## - Estimator : G-formula
# ## - Uncertainty: Gaussian approximation
# ## where the variance is estimated via the influence function
# ##
# ## Testing procedure
# ## - Null hypothesis : given two treatments (A,B) and a specific timepoint, equal risks
# ## - Confidence level : 0.95
# ##
# ## Results:
# ## - Difference in standardized risk (B-A) between time zero and 'time'
# 5
# ## reported on the scale [-1;1] (difference between two probabilities)
# ## (difference in average risks when treating all subjects with the experimental treatment (B),
# ## vs. treating all subjects with the reference treatment (A))
# ##
# ## time chemo=A risk(chemo=A) chemo=B risk(chemo=B) difference ci
# ## 2000 0 0.698 1 0.527 -0.172 [-0.26;-0.09]
# ## p.value
# ## 5.93e-05
# ##
# ## difference : estimated difference in standardized risks
# ## ci : pointwise confidence intervals
# ## p.value : (unadjusted) p-value
# or with survivalG function from the mets package
fit.ph <- phreg(Surv(time,event)~ chemo+age+size+log2nodes,x=TRUE,data=bca)
summary(survivalG(fit.ph, time = 2000, data = bca))
# ## G-estimator :
# ## Estimate Std.Err 2.5% 97.5% P-value
# ## risk0 0.6984 0.04076 0.6185 0.7783 8.208e-66
# ## risk1 0.5266 0.02253 0.4824 0.5707 8.268e-121
# ##
# ## Average Treatment effect: difference (G-estimator) :
# ## Estimate Std.Err 2.5% 97.5% P-value
# ## ps0 -0.1718 0.04404 -0.2581 -0.08549 9.582e-05
# ##
# ## Average Treatment effect: ratio (G-estimator) :
# ## log-ratio:
# ## Estimate Std.Err 2.5% 97.5% P-value
# ## [ps0] -0.2823682 0.06792267 -0.4154942 -0.1492422 3.221699e-05
# ## ratio:
# ## Estimate 2.5% 97.5%
# ## 0.7539960 0.6600140 0.8613604
# ##
# ## Average Treatment effect: survival-difference (G-estimator) :
# ## Estimate Std.Err 2.5% 97.5% P-value
# ## ps0 0.1718064 0.04404226 0.08548519 0.2581277 9.581793e-05
# ##
# ## Average Treatment effect: 1-G (survival)-ratio (G-estimator) :
# ## log-ratio:
# ## Estimate Std.Err 2.5% 97.5% P-value
# ## [ps0] 0.4508392 0.1375572 0.1812321 0.7204463 0.001047428
# ## ratio:
# ## Estimate 2.5% 97.5%
# ## 1.569629 1.198693 2.055350
# 8. Compare this to the non-standardized difference (this corresponds to the difference of the
# unadjusted survival curves from question 1).
fit.cox0 <- coxph(Surv(time,event)~ strata(chemo),x=TRUE,data=bca)
ate0 <- ate(fit.cox0, data = bca, treatment = "chemo", times = 2000, se = TRUE)
# ## Input variables
# ## - Treatment : chemo (2 levels: "0" "1")
# ## - Event : event (cause: TRUE, censoring: 0)
# 6
# ## - Time [min;max] : time [38;6270]
# ## - at risk/time : 2000
# ## number in treatment 0 23
# ## number in treatment 1 210
# ##
# ## Estimation procedure
# ## - Estimator : G-formula
# ## - Uncertainty: Gaussian approximation
# ## where the variance is estimated via the influence function
# ##
# ## Processing
# ## - Prepare influence function: outcome done
# ## - Point estimation: done
# ## - Decomposition iid: done
# ## - Confidence intervals: done
summary(ate0)
# ## Average treatment effect
# ##
# ## - Treatment : chemo (2 levels: "0" "1")
# ## - Event : event (cause: TRUE, censoring: 0)
# ## - Time [min;max] : time [38;6270]
# ## - at risk/time : 2000
# ## number in treatment 0 23
# ## number in treatment 1 210
# ##
# ## Estimation procedure
# ## - Estimator : G-formula
# ## - Uncertainty: Gaussian approximation
# ## where the variance is estimated via the influence function
# ##
# ## Testing procedure
# ## - Null hypothesis : given two treatments (A,B) and a specific timepoint, equal risks
# ## - Confidence level : 0.95
# ##
# ## Results:
# ## - Difference in standardized risk (B-A) between time zero and 'time'
# ## reported on the scale [-1;1] (difference between two probabilities)
# ## (difference in average risks when treating all subjects with the experimental treatment (B),
# ## vs. treating all subjects with the reference treatment (A))
# ##
# ## time chemo=A risk(chemo=A) chemo=B risk(chemo=B) difference ci
# ## 2000 0 0.772 1 0.496 -0.276 [-0.37;-0.19]
# ## p.value
# ## 1.25e-09
# ##
# ## difference : estimated difference in standardized risks
# ## ci : pointwise confidence intervals
# ## p.value : (unadjusted) p-value

library(survival)
library(timereg)
library(prodlim)
GBSG2 <- read.csv2("http://biostat.ku.dk/frank/data/GBSG2.csv")
head(GBSG2)


# Q1-1 Draw survival curves to illustrate the recurrence free survival for patients who did and did not receive hormonal therapy. What do you conclude?

GBSG2$years <- GBSG2$time/365.24
out=survfit(Surv(years, cens) ~horTh, data=GBSG2)
summary(out)


plot(out)
kmplot(out)
title(main="Kaplan-Meyer Plot of GBSG2")

out=survdiff(Surv(years, cens) ~horTh, data=GBSG2)
out

#The conclusion is that with a p-value of 0.0002 it's basically a difference.

# Q1-2 Modify the curves to show the absolute risk of the composite event "recurrence or death" by adding type="risk" to the prodlim plot command.
out=prodlim(Surv(years, cens) ~horTh, data=GBSG2, type="risk")

plot(out)


title(main="Hazard Plot of GBSG2-risk")

out=survdiff(Surv(years,cens) ~horTh, data=GBSG2)
out

# Q1-3A What is the risk of the composite event of recurrence and death in each group within 1 year after mastectomy? Within 5 years?
#   horTh time  surv   risk
# 3   yes    1 0.950 0.0504

out=prodlim(Surv(years, cens) ~horTh, data=GBSG2, type="risk")
summ <- summary(out, times = c(1, 5))
transform(summ[, c("horTh","time","surv")], risk = 1 - surv)


# Q1-4A What is the median time to the event in each group?
# Median with interquartile range (IQR):
#     horTh             Median (IQR)
# 1:     no 1528.00 (629.00;2456.00)
# 2:    yes      2018.00 (859.00;NA)
out = prodlim(Surv(time, cens) ~ horTh, data = GBSG2)
quantile(out, probs = 0.5)
#
#

# Q2-1 Is there evidence for an association between hormonal treatment and the relapse free survival time?
# Test the null hypothesis that the hazards for recurrence or death is the same in the groups by a log-rank test.
#             N Observed Expected (O-E)^2/E (O-E)^2/V
# horTh=no  440      205      180      3.37      8.56
# horTh=yes 246       94      119      5.12      8.56

# Chisq= 8.6  on 1 degrees of freedom, p= 0.003
# Conclusion: rejection with p-value of 0.003
out=survdiff(Surv(time,cens)~horTh, data=GBSG2)
out


# Q2-2 Now assume that the hazards of the composite event are proportional in the two groups. Fit a proportional hazards model.
# Interpret the output.
# What is the value of the estimated hazard ratio comparing the relapse free survival of hormonal treated and untreated patients?
#          exp(coef) exp(-coef) lower .95 upper .95
# horThyes    0.6949      1.439    0.5438    0.8879
# the hazard ratio is about 69,5%
cox_model <- coxph(Surv(years, cens) ~ horTh, data = GBSG2)
summary(cox_model)


# Q2-3 Compare the p-value obtained from the proportional hazards analysis (column Pr(>|z|)) with that from the log rank test.
# Wald test            = 8.47  on 1 df,   p=0.004
# Score (logrank) test = 8.57  on 1 df,   p=0.003
#
# The p values are basically the same.

# Q2-4A Redo the previous analysis only this time estimate the hazard ratio comparing the recurrence free survival of untreated to treated patients.
# Compare the output you get to that from the previous step.
#         exp(coef) exp(-coef) lower .95 upper .95
# horThno     1.439     0.6949     1.126     1.839
# the reverse hazard ratio is 1/0.6949= 1.439
GBSG2$horTh <- factor(GBSG2$horTh, levels = c("yes", "no"))

cox_model <- coxph(Surv(years, cens) ~ horTh, data = GBSG2)
summary(cox_model)
#

#  Now return to question Q1-3A  where you estimated the absolute risks in the groups and use the risks from there.
# You don't have to run any new data analysis for this question.

# Q3-1 What is the ratio of the risk of the event within one year for treated to the risk of untreated patients?
# horTh time  surv   risk
# 1   yes    1 0.950 0.0504
# 3    no    1 0.897 0.1034
0.0504/0.1034
# [1] 0.4874275
# Q3-2A What is the ratio of the risk of the event within five years for treated to the risk of untreated patients?
# horTh time  surv   risk
# 2   yes    5 0.581 0.4188
# 4    no    5 0.437 0.5632
0.4188/ 0.5632
# [1] 0.743608
# Q3-3 Compare the risk ratios above to the hazard ratio. What do you conclude?
# hormone treatment helps, particularly early on.
#
#

quantile(GBSG2$tsize)
hist(GBSG2$tsize, xlab="mm", main="Tumour size")
GBSG2$grouptsize <- cut(GBSG2$tsize, breaks=c(-Inf,20,30,Inf),labels=c("small","medium","large"))
library(mets)
dcut(GBSG2, breaks=c(-Inf,20,30,Inf),labels=c("small","medium","large")) <- grouptsize~tsize
table(GBSG2$grouptsize)
stripchart(tsize ~ grouptsize, method='jitter', vertical=TRUE, data=GBSG2)



# # Q4-1 Make Kaplan-Meier plots for all-cause mortality and tumour thickness in groups. 
# Which of the three groups has the best recurrence free survival prognosis?
# small tumor group has best survival?
out=survfit(Surv(time,cens) ~grouptsize, data=GBSG2)
plot(out)
kmplot(out)
title(main="Kaplan-Meyer of GBSG2 by tumor groups")

# Q4-2A Assuming that the hazards are proportional in the three tumour thickness groups, 
# what is the hazard ratio comparing the hazard in the group with the largest tumours to the hazard in the group with middle sized tumours?
#exp(coef) exp(-coef) lower .95 upper .95
#grouptsizemedium     1.346     0.7429    0.9961     1.819
#grouptsizelarge      1.751     0.5712    1.2861     2.383
cox_model <- coxph(Surv(years, cens) ~ grouptsize, data = GBSG2)
summary(cox_model)
1.7507 /1.3461
# [1] 1.300572

# Q4-3A Make a log-rank test using the tumour size in the three groups. 
# Think about what the null hypothesis states when there are three groups. 
# Is the null hypothesis violated if the survival is the same in the two groups with thicker tumours, but better for those with a smaller tumour?
# bigger than small tumor has a different outcome with a p-value of 0.003
GBSG2$grouptsize_big <- factor(
  ifelse(GBSG2$grouptsize %in% c("medium", "large"), "big", "small"),
  levels = c("small", "big")
)
out=survdiff(Surv(time,cens)~grouptsize_big, data=GBSG2)
out

fit <- survfit(Surv(years, cens) ~ grouptsize, data = GBSG2)

summary(fit, times = 5)
# DONE

# Stanford heart transplant study
#
# Crowley and Hu (J. Amer. Statist. Assoc., 1977) describe a program at Stanford University where 103 patients were registered to recieve a heart transplant and were followed until death or censorship. Of these, 65 received transplants during follow-up, whereas 38 did not
#
# Question: Does getting a heart transplant improve survival?
#
#   The variables in the data set are
# age 	    Age in years at entry into the study
# death 	  0=censored
#           1=dead
# days 	    The follow-up time, days from entry to death/censoring
# trans 	  0=no heart transplant
#           1=heart transplant
# wait 	    Days from entry to transplant
#           If trans=0 then wait is NA
# mismatch 	0=Human leukocyte antigen (HLA) type in donor and patient match
#           1=HLA mismatch
#           If trans=0 then mismatch is NA
#
# Load the data into R and inspect the six first records
library(prodlim);
library(survival);
library(mets);


heart <- read.csv2("http://biostat.ku.dk/frank/data/stanford.csv")
head(heart)

# Q1-1. Plot Kaplan-Meier curves treating transplantation as a time-fixed covariate
# and fit a simple Cox regression model where you “evaluate the effect of transplantation”
# using transplantation as a time-fixed covariate. How should this analysis be interpreted?
# Think carefully in terms of the information provided by the covariate trans.
# Use follow-up time as time axis in the Cox model and adjust for age at entry.
ss=prodlim(Surv(days,death) ~trans+age,data=heart)
plot(ss);

out=coxph(Surv(days,death) ~ trans+age,data=heart)
newdata=data.frame(trans=c(0,1));
pred=survfit(out,newdata);
plot(pred);

plot(pred, col = 1:2, lty = 1:2,
     xlab = "Follow-up time",
     ylab = "Survival probability")
legend("bottomleft", legend = levels(factor(heart$trans)),
       col = 1:2, lty = 1:2)
exp(coef(out))
exp(confint(out))

# This comparison should not be interpreted as the causal effect of transplantation,
# because the transplanted group is defined using future information and necessarily excludes early deaths before transplantation.
# Thus, the analysis is affected by immortal time bias and will typically overstate the benefit of transplantation.
# Transplantation should be analyzed as a time-dependent covariate in a Cox model. Then patients contribute person-time as:
# - non-transplanted before transplant
# - transplanted after transplant


# Q1-2A. Make a landmark analysis at time 30 days on the waiting list.
# That is, consider only patients still alive and under observation after 30 days,
# and use this landmark time as time origin (“time zero”).
# How many patients are there in each group at the (new) baseline?
#   How many were excluded because they were already dead by day 30?
#   How many were excluded because they were censored?

landmark <- 30

lmdata <- subset(heart, days > landmark)

excluded_dead <- subset(heart, days <= landmark & death == 1)
excluded_cens <- subset(heart, days <= landmark & death == 0)

nrow(lmdata)
nrow(excluded_dead)
nrow(excluded_cens)

lmdata$trans30 <- ifelse(lmdata$trans == 1 & lmdata$wait <= landmark, 1, 0)

table(lmdata$trans30, useNA = "ifany")

table(lmdata$trans30)
sum(lmdata$trans30 == 1)

# only patients still alive and under observation at day 30 were included.
# Patients who had died by day 30 were excluded, and patients censored before or at day 30 were also excluded.
# At the new baseline, patients were classified according to whether they had already received transplantation by day 30.

#   Q1-3A Draw Kaplan-Meier curves according to the transplantation status at the landmark time.
# What is the prognosis after two years on the waiting list?
#   Does this analysis suffer from immortal time bias?

lmdata$time_lm <- lmdata$days - landmark
lmdata$event <- ifelse(lmdata$death == 1, 1, 0)


km_lm <- survfit(Surv(time_lm, event) ~ trans30, data = lmdata)
par(mfrow=c(1,1))
plot(km_lm, col = 1:2, lty = 1:2,
    xlab = "Days since landmark (day 30)",
    ylab = "Survival probability")
legend("bottomleft", legend = levels(factor(lmdata$trans30)),
       col = 1:2, lty = 1:2)
ss700 <- summary(km_lm, times = 700)
round(100 * ss700$surv)


data.frame(
  group = ss700$strata,
  surv_percent = round(100 * ss700$surv)
)

# The landmark analysis is specifically designed to avoid the classic immortal time bias
# from classifying patients by future transplantation status.

#   Q1-4. In a landmark analysis after 30 days, calculate the mortality hazard ratio
# comparing transplanted patients according to their transplantation status at the landmark time.

cox_lm <- coxph(Surv(time_lm, event) ~ trans30, data = lmdata)
summary(cox_lm)
exp(cbind(HR = coef(cox_lm), confint(cox_lm)))

# The HR from this model compares mortality after day 30 between:
# - patients transplanted by day 30
# - patients not transplanted by day 30
# among patients alive and still observed at day 30.

# Q1-5A. Can you think about any potential problems in interpreting the
# “effect” of heart transplant on survival in the landmark analysis?
#   Think in terms of selection effects? Who are in the groups at day 30?
#   Who are excluded?
#   Are there other aspects to that might need consideration here?


# The landmark analysis is still subject to selection bias.
# Only patients alive and under follow-up at day 30 are included,
# while those who died or were censored earlier are excluded.
# Thus, the groups compared at day 30 are selected survivors and may differ
# systematically in disease severity, transplant eligibility, and other prognostic factors.
# The estimated association therefore cannot be interpreted straightforwardly
# as a causal effect of heart transplantation.

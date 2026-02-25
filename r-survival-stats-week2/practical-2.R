# Exercise: Vaccinations in Guinea-Bissau
# In rural Guinea-Bissau, 5274 children under 7 months of age were visited two times at home during the years 1990-1996,
# with an interval of six months. Information about vaccination (BCG and DTP) was collected at each visit and at the second visit death during follow-up was registered.
# Some children were censored because they moved away during follow-up or survived until the second visit.

# The variables in the data set are
# id 	Id number
# fuptime 	Follow-up time in days. Time 0 corresponds to first visit.
# dead 	0 = censored, 1 = dead
# bcg 	BCG vaccinated at first visit; 0 = no, 1 = yes
# dtp 	Number of DTP doses at first visit
# age 	Age at first visit in days
# agem 	Age at first visit in months

# Load the data from the file bissau.csv and store it in a data frame with the name bissau by the command
# Remember to load the survival and prodlim packages with library(survival) and library(prodlim).

library(survival)
library(prodlim)
library(timereg)
bissau <- read.csv2("http://biostat.ku.dk/frank/data/bissau.csv")

# 1. Make a Kaplan-Meier plot for children with (DTP>0) and without (DTP=0).
# Use time on study as time scale (time from first visit).
# Calculate the survival probability at time 180 days after first visit.
bissau$dtpcat <- 1* (bissau$dtp>0)
bissau <- transform(bissau, outage = age+fuptime, inage=age)

out=survfit(Surv(fuptime, dead) ~dtpcat, data=bissau)
plot(out, mark.time=TRUE)
kmplot(out, ylim=c(0.85,1.0))
title(main="Kaplan-Meyer Plot of DTP>1 vs DTP=0 for bissau by time of inclusion")

out=prodlim(Surv(inage, outage, dead) ~dtpcat, data=bissau)
ylim <- c(0.85,1.0)
plot(out, ylim=ylim, marktime=TRUE, legend="bottom-right")
title(main="Kaplan-Meyer Plot of DTP>1 vs DTP=0 for bissau by age of child")

summary(out, times=180)
# 2.Test if having received a DTP dose is good by a Cox model. Write a conclusion sentence.
# Also test the association between DTP any dose and mortality with a log-rank test.
cox <- coxph(Surv(fuptime, dead) ~dtpcat, data=bissau)
summary(cox)
survdiff(Surv(fuptime,dead)~dtpcat,data=bissau)
# Conclusion, we cannot say, since p Hazard ratio is small, confidence interval crosses 1.0 and p-value is high
# We cannot say there is an association


# 3. Are children vaccinated with DTP also vaccinated with BCG?
# Make a 2x2 table of BCG and DTP (any dose) to check this.
table(dtp=bissau$dtpcat, bcg=bissau$bcg, useNA="ifany") 
# It looks like there is a strong association between getting dtp and having bcg

# 4. Perform a Cox regression to assess the association between DTP any dose and mortality adjusted for BCG.
# Why is this a good idea? Write a conclusion sentence.
cox <- coxph(Surv(fuptime, dead) ~dtpcat+bcg, data=bissau)
summary(cox)
# given the strong coupling of receiving a dtp if already having bcg, we want to check the dtp impact alone (and it looks bad)


# 5. Redo the analysis from the previous question adjusting for the age of the children by using age as time-scale.
cox <- coxph(Surv(inage, outage, dead) ~dtpcat+bcg, data=bissau)
summary(cox)
# with age as time-scale, the p-value increases significantly to >10%, so we cannot disregard the null hypothesis when comparing by age


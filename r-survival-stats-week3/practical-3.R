# We consider the TRACE data that was collected by the TRACE study group, see for example Jensen, G.V., Torp-Pedersen, C., Hildebrandt, P., Kober, L., F. E. Nielsen, Melchior, T., Joen, T. and P. K. Andersen (1997), Does in-hospital ventricular fibrillation affect prognosis after myocardial infarction?, European Heart Journal 18, 919-924.

# The data contains 1877 patients and is a subset of a data set consisting of approximately 6000 patients. It contains data relating survival of patients after myocardial infarction to various risk factors.

#    time : observation time (years)
#    status : status at exit (0: censored, 9: dead from myorcaridal, 8,7:other death )
#    chf : (1: present, 0: not present )
#    vf : (1: present, 0: not present )
#    diabetes : (1: present, 0: not present )
#    sex : sex of patient
#    age a numeric vector code. Age of patient.
#    wmi a numeric vector. Measure of heart pumping effect based on ultrasound measurements where 2 is normal and 0 is worst.

# To load the data

# library(mets);
# library(prodlim);
# library(survival);
# library(timereg)
# data(TRACE)

# head(TRACE)

#        id wmi status chf    age sex diabetes  time vf
# X6440 6440 1.9      9   0 66.114   0        0 6.237  0
# X5201 5201 2.0      9   1 67.668   1        0 0.858  0
# X2643 2643 1.8      8   0 42.243   1        1 6.154  0
# X1667 1667 1.4      9   1 77.579   0        1 0.036  1
# X79     79 1.6      9   0 51.863   1        0 4.609  0
# X3334 3334 1.2      0   1 63.559   0        0 7.647  0

# For the non-R users

# ttt <- read.table("https://www.biostat.ku.dk/ts/survival/data/TRACE.txt")

# R commands for GOF

# Some relevant R commands from the cheat-sheet

# Cumulative residuals for proportionality

# library(mets)
# data(pbc); pbc=transform(pbc,lbili=log(bili),lpro=log(protime))
# cox1 <- phreg(Surv(time,status==2)~age+edema+lbili+lpro+albumin,data=pbc)
# gof1 <- gof(cox1)

# par(mfrow=c(2,3));
# plot(gof1);

# Cumulative residuals for linearity

# pbc <- na.omit(pbc)
# cox2 <- gofZ.phreg(Surv(time,status==2)~age+edema+lbili+lpro+albumin,data=pbc)
# summary(cox2)

# Timedependent Cox to describe effects

# coxt <- coxph(Surv(time,status==2)~age+edema+lbili+lpro+albumin+tt(edema),data=pbc,
#              tt=function(x,t,...) x*(t>1))
# summary(coxt)

# The graphical old-school test

# cox1 <- phreg(Surv(time,status==2)~age+strata(edema)+lbili+lpro+albumin,data=pbc)
# bplot(cox1,log="y")

# we can also make the time-axis on log-scale, but we are still cheking if lines are parallel l

# cox1 <- phreg(Surv(time,status==2)~age+strata(edema)+lbili+lpro+albumin,data=pbc)
# bplot(cox1,log="xy")

# Same plot (up to handling of ties) are done by

# cox1 <- coxph(Surv(time,status==2)~age+strata(edema)+lbili+lpro+albumin,data=pbc)
# ps <- survfit(cox1,data.frame(age=40,edema=c(0,0.5,1),lbili=0,lpro=0,albumin=0))
# kmplot(ps,fun="cloglog",col=1:3)

# Q1. GOF for Cox model

# We look a cox model with vf, diabetes and age, to keep it simple at a first go.

# You might start by thinking about what we expect considering the predictors in terms of how the predict the risk of death.

#    check for proportionality using cumulative residuals
#        what can we do to remedy problems
#    check for linearity using cumulative residuals
#    check for proportionality using the graphical approach

# What do you conclude about the suitability of the Cox model.

# We shall not explore any possible correlations, again to keep life simple.
# Q1. Improving the model
# We look a cox model with vf, age and diabetes as predcitors.

#    We can stratify on vf and examine the fit
#        proportionality
#        linearity
#    Make interactions with time
#        still need linearity
#        we would cut time early
#    Settle for conclusions at specific time points
#        binomial regression (rmst analysis later)

# Q3. Describing the effect of vf using Cox models.

# We shall now fit a Cox model with an early and a late effect. Looking at the survival plots or the cumulative residual plot you may guess on a useful cut-point for the early and late effects.

# Try to express the results in a useful way.
# Q4. Survival Predictions
# Using the simple cox model.

#    Make predictions using the simple Cox model
#    possibly also compare with binomial regression at specfic time points (not now).

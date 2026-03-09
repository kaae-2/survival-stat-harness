# PBC3

# PBC3 was a multi-centre randomized clinical trial conducted in six European hospitals (Lombard, et al., Gastroenterology, 1993, and Andersen & Skovgaard, Regression with linear predictors, 2010). We consider  the subset of patients with information on biochemical markers for liver function. Between 1 Jan. 1983 and 1 Jan. 1987, 341 patients with the liver disease primary biliary cirrhosis (PBC) were randomized to either treatment with Cyclosporin A (CyA, 173 patients) or placebo (168 patients).

# Patients were then followed from randomization until treatment failure, drop-out or 1 Jan, 1989

# • 59 patients died (CyA: 28, placebo: 31)
# • 28 were transplanted (CyA: 14, placebo: 14)
# • 4 patients were lost to follow-up before 1 January 1989.

# At entry a number of clinical, biochemical and histological variables were recorded.
# days 	observation time (days)
# status 	status at exit (0: censored, 1: liver transplantation, 2 : dead)
# tment 	treatment (0: placebo, 1: CyA)
# sex 	(1: males, 0: females)
# age 	years
# alb 	albumin (g/L)
# bili 	bilirubin (micromoles/L)
# asptr 	aspartate transaminase (IU/L)


library(prodlim);
library(survival);
library(mets);
# The purpose of the trial was to study the effect of treatment on the survival time. However, during the course of the trial an increased use of liver transplantation for patients with this disease made the investigators redefine the main response variable to be time to “failure of medical treatment” defined as either death or liver transplantation.


# Load the data from the file pbc3subset.csv and store it in a data frame with the name pbc3 by the command

pbc3 <- read.csv2("http://biostat.ku.dk/frank/data/pbc3subset.csv")

# We will study “treatment failure” defined as either death or liver transplantation.
# For this purpose make a new status variable which is 0 if the patient was still alive
# and without liver transplant when the study closed in 1989 and 1 if the patient died or had a liver transplant.
# Q1-1 Draw Kaplan-Meier curves to compare survival in the treatment group.

pbc3$dead <- pbc3$status!=0
outkm <- survfit(Surv(days, dead) ~tment, data= pbc3);

kmplot(outkm)

# look the same??

# Q1-2A Fit a simple Cox model with treatment as explanatory variable to test if the CyA has
# an effect on the outcome “treatment failure”. Is the CyA treatment effective?
cox1 <- coxph(Surv(days, dead) ~tment, data=pbc3);
summary(cox1)
(1-0.91287)*100
#based on this, can't say yes

# Q1-3 Assess the proportional hazards assumption by the test and graphical procedure based
# on the score process (cumulative martingale residual procedure for proportionality).

prop <- phreg(Surv(days, dead) ~tment, data=pbc3);
gofprop <- gof(prop);
plot(gofprop)

#it is proportional, because the p-value is above 0.05?
# the plot also looks fine

# You can get the test for proportional hazards by fitting the Cox model with
# the phreg command from the mets package and then applying the gof function to the phreg object.
# You get the plots with the generic plot command, e.g., plot(gof.fit)
# where gof.fit is the result from the gof command.
# PBC is a slowly progressing liver disease with patients diagnosed at varying disease stages,
# thereby making populations of patients with PBC rather heterogeneous.
# The liver function is not measured directly but only via biochemical markers such as bilirubin,
# albumin and aspartate transminase. From previous studies
# it is known that low albumin values and high levels of the other two markers indicate poor liver function.

# Q2-1A Consider how these variables affect survival by fitting simple unadjusted Cox models for each of them separately.
# Assess the functional form of the covariates (log-hazard linearity) by the method based on cumulative martingale residuals.
# Try also to log-transform each variable and assess the fit for the transformed variable.

pbc3 = transform(pbc3, lalb=log(alb), lbili=log(bili), lasptr=log(asptr));
cox1 <- phreg(Surv(days, dead) ~tment+alb+lalb+bili+lbili+asptr+lasptr, data=pbc3);
gvf <- gof(cox1);
par(mfrow=c(3,3))
plot(gvf)
summary(gvf);

vars <- c("tment","alb","lalb","bili","lbili","asptr","lasptr")

par(mfrow = c(3,3))   # one plot at a time

for(v in vars){

  form <- as.formula(paste("Surv(days, dead) ~", v))

  fit <- phreg(form, data = pbc3)

  gvf <- gof(fit)

  plot(gvf, main = paste("Cumulative score process:", v))

  print(summary(gvf))
}

# just bili?

# To test for linearity (on the log-hazard scale) you can use command gofZ.phreg from the mets package.
# To get the plot, use the plot function with the extra argument type='z', i.e., use plot(gofZ.fit, type='z')
# where gofZ.fit is the result from the gofZ.phreg command.

vars <- c("alb","lalb","bili","lbili","asptr","lasptr")

par(mfrow = c(3,3))   # one plot at a time

for(v in vars){

  form <- as.formula(paste("Surv(days, dead) ~", v))

  fit <- tryCatch(
    gofZ.phreg(form, data = pbc3),
    error = function(e) {
      message("Skipping ", v, ": ", e$message)
      return(NULL)
    }
  )
  if (is.null(fit)) next


  plot(fit, type='z')

  print(summary(fit))
}

# Based on the p-values of the tests found under “Cumulative residuals versus model matrix”, which variables (raw and log-transformed) are rejected at the 5% significance level?

# bili and asptr
# apparently its only bili


# Q2-2 Assess the proportional hazards assumption for the functional form you find suitable in question Q2-1A

par(mfrow=c(2,3))
cox_final <- phreg(Surv(days, dead) ~ tment +lalb + lbili + lasptr, data = pbc3)
gof_final <- gof(cox_final)
plot(gof_final)
summary(gof_final)

# the proportional assumption holds

# Q2-3A Do your results verify that a poor liver function at randomization yields a higher hazard for “treatment failure”?

# fit <- phreg(Surv(days, dead) ~ lalb + lbili + lasptr, data = pbc3)
fit <- coxph(Surv(days, dead) ~ lalb + lbili + lasptr, data = pbc3)
s <- summary(fit)
s
# apparently all of them have hazard 

# hazard ratios with CI
HR  <- exp(coef(fit))
CI  <- exp(confint(fit))
cbind(HR, CI)

# the albumin and bilirubin are significant.

# Because this trial was randomized, the liver function should be balanced in the CyA and placebo groups.
# Thus, the simple analysis in Q1 including only treatment could be considered the “correct” analysis and end
# of the story.

# It should however be kept in mind that randomization only ensures complete balancing of treatment groups for very large trials.
# You will investigate if patients in one of the treatment groups is more severely ill than those in the other group, despite the randomization.

# Q3-1 For each of the three biomarkers bilirubin, albumin and aspartate transminase
# draw box plots for the CyA treatment and placebo treated groups.
# Do the distributions appear to be similar? Which of the treatment groups has the worst liver function?

par(mfrow = c(2,3))

boxplot(bili ~ tment, data = pbc3,
        main = "Bilirubin by treatment",
        xlab = "Treatment",
        ylab = "Bilirubin")

boxplot(alb ~ tment, data = pbc3,
        main = "Albumin by treatment",
        xlab = "Treatment",
        ylab = "Albumin")

boxplot(asptr ~ tment, data = pbc3,
        main = "AST by treatment",
        xlab = "Treatment",
        ylab = "AST")

boxplot(lbili ~ tment, data = pbc3, main="log(Bilirubin)")
boxplot(lalb ~ tment, data = pbc3, main="log(Albumin)")
boxplot(lasptr ~ tment, data = pbc3, main="log(AS)")

# i would say they are close, maybe log albumin and bilirubin

# Q3-2A Calculate the mean and standard deviation for albumin, log2-transformed bilirubin
# and log2-transformed aspartate transminase for the CyA treatment and placebo treated groups.
# Which of the treatment groups has the worst liver function?

tapply(pbc3$alb, pbc3$tment,
       function(x) c(mean=mean(x,na.rm=TRUE), sd=sd(x,na.rm=TRUE)))

tapply(pbc3$lbili, pbc3$tment,
       function(x) c(mean=mean(x,na.rm=TRUE), sd=sd(x,na.rm=TRUE)))

tapply(pbc3$lasptr, pbc3$tment,
       function(x) c(mean=mean(x,na.rm=TRUE), sd=sd(x,na.rm=TRUE)))

# the treatment group with the low albumin has the worst function?

# Q3-3 It is the size of the difference in the markers (in absolute terms) between the groups that matters.
# Because the study was randomized, it does not make sense to make statistical significance tests to compare the groups. Why?
# With the findings from the previous steps we fit a multiple Cox model for treatment
# adjusting for the three liver function markers.
# No interactions were considered of special interest in this study,
# so we don’t pursue a study of interactions between treatment and other variables.

coxtrt <- coxph(Surv(days, dead) ~ tment + alb + bili + asptr, data = pbc3)
summary(coxtrt)

# it does not make sense to compare the groups because they are the same? that's the whole point of randomizing?
#
#


# Q4-1. Fit a Cox model with treatment, albumin, log2(bilirubin) and log2(aspartate transaminase)
# as explanatory variables. Check the model assumptions. Interpret the results.

coxltrt <- coxph(Surv(days, dead) ~ tment + alb + lbili + lasptr, data = pbc3)
summary(coxltrt)

# treatment halves hazard significantly, doubling the bili 2.5x the risk, higher albumin lowers the hazard by 10%
# lasptr is not significant

# Q4-2A Compare your adjusted model to the unadjusted model from Q1. Is the “message” on
# the potential effect of the treatment different? Given the randomization, which of the models
# do you prefer?

# with the simple model, treatment seems to have no significant effect.
# Adjusting for the different proteins, we do find a significant value of the treatment
# I would prefer the adjusted model, as we're trying to isolate the effect of treatment

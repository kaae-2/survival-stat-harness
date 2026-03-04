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



# The purpose of the trial was to study the effect of treatment on the survival time. However, during the course of the trial an increased use of liver transplantation for patients with this disease made the investigators redefine the main response variable to be time to “failure of medical treatment” defined as either death or liver transplantation.


# Load the data from the file pbc3subset.csv and store it in a data frame with the name pbc3 by the command

# pbc3 <- read.csv2("http://biostat.ku.dk/frank/data/pbc3subset.csv")

# We will study “treatment failure” defined as either death or liver transplantation. For this purpose make a new status variable which is 0 if the patient was still alive and without liver transplant when the study closed in 1989 and 1 if the patient died or had a liver transplant.
# Q1-1 Draw Kaplan-Meier curves to compare survival in the treatment group.

# Q1-2A Fit a simple Cox model with treatment as explanatory variable to test if the CyA has
# an effect on the outcome “treatment failure”. Is the CyA treatment effective?

# Q1-3 Assess the proportional hazards assumption by the test and graphical procedure based
# on the score process (cumulative martingale residual procedure for proportionality).

# You can get the test for proportional hazards by fitting the Cox model with the phreg command from the mets package and then applying the gof function to the phreg object. You get the plots with the generic plot command, e.g., plot(gof.fit) where gof.fit is the result from the gof command.
# PBC is a slowly progressing liver disease with patients diagnosed at varying disease stages, thereby making populations of patients with PBC rather heterogeneous. The liver function is not measured directly but only via biochemical markers such as bilirubin, albumin and aspartate transminase. From previous studies it is known that low albumin values and high levels of the other two markers indicate poor liver function.

# Q2-1A Consider how these variables affect survival by fitting simple unadjusted Cox models for each of them separately. Assess the functional form of the covariates (log-hazard linearity) by the method based on cumulative martingale residuals. Try also to log-transform each variable and assess the fit for the transformed variable.

# To test for linearity (on the log-hazard scale) you can use command gofZ.phreg from the mets package. To get the plot, use the plot function with the extra argument type='z', i.e., use plot(gofZ.fit, type='z') where gofZ.fit is the result from the gofZ.phreg command.

# Based on the p-values of the tests found under “Cumulative residuals versus model matrix”, which variables (raw and log-transformed) are rejected at the 5% significance level?

# Q2-2 Assess the proportional hazards assumption for the functional form you find suitable in question Q2-1A

# Q2-3A Do your results verify that a poor liver function at randomization yields a higher hazard for “treatment failure”?

# Because this trial was randomized, the liver function should be balanced in the CyA and placebo groups. Thus, the simple analysis in Q1 including only treatment could be considered the “correct” analysis and end
# of the story.

# It should however be kept in mind that randomization only ensures complete balancing of treatment groups for very large trials. You will investigate if patients in one of the treatment groups is more severely ill than those in the other group, despite the randomization.

# Q3-1 For each of the three biomarkers bilirubin, albumin and aspartate transminase draw box plots for the CyA treatment and placebo treated groups. Do the distributions appear to be similar? Which of the treatment groups has the worst liver function?

# Q3-2A Calculate the mean and standard deviation for albumin, log2-transformed bilirubin and log2-transformed aspartate transminase for the CyA treatment and placebo treated groups. Which of the treatment groups has the worst liver function?

# Q3-3 It is the size of the difference in the markers (in absolute terms) between the groups that matters. Because the study was randomized, it does not make sense to make statistical significance tests to compare the groups. Why?


# With the findings from the previous steps we fit a multiple Cox model for treatment adjusting for the three liver function markers. No interactions were considered of special interest in this study, so we don’t pursue a study of interactions between treatment and other variables.

# Q4-1. Fit a Cox model with treatment, albumin, log2(bilirubin) and log2(aspartate transaminase) as explanatory variables. Check the model assumptions. Interpret the results.

# Q4-2A Compare your adjusted model to the unadjusted model from Q1. Is the “message” on
# the potential effect of the treatment different? Given the randomization, which of the models
# do you prefer?

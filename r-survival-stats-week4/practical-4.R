

# Practical Competing risks I
# Pneumonia on admission to intensive care unit

# The dataset sir.adm from the compeir package contains data on a random subsample of 747 
# intensive care unit patients from the SIR 3 (Spread of nonsocomial Infections and Resistant pathogens) 
# cohort collected to examine the effect of hospital-acquired infections in intensive care (Wolkewitz et al., 2008). 
# Competing endpoints are discharge from the unit and death on the unit. 
#The data was analyzed in the book Competing risks and multistate models with R (Beyersmann, Schumacher, Allignol, Springer, 2012).

# The variables
# Name 	Content
# time 	length of stay (days)
# status 	0, censored before end of unit stay; 1, discharged (alive); 2 death on unit
# pneu 	pneumonia status on admission; 0 no pneumonia, 1 pneumonia

library(prodlim);
library(survival);
library(mets);


# Load the data from the file sir.adm.csv and store it in a data frame with the name sir.adm by the command

sir.adm <- read.csv2("http://biostat.ku.dk/frank/data/sir.adm.csv")
table(pneumonia=sir.adm$pneu,useNA="ifany")
# We will study the impact of pneumonia on on-unit mortality. 
# Pneumonia is a severe infection, suspected to increase mortality. 
# Death is the event of interest and discharge is the competing event.

# 1. Estimate the cause-specific cumulative hazards for both unit death and discharge by the Nelson-Aalen estimator. 
# Plot the cumulative hazards for patients with and without pneumomia at admission.
par(mfrow=c(1,2))
# cc=phreg(Surv(time,status==1)~strata(pneu),data=sir.adm); 
# plot(cc);
# cc=phreg(Surv(time,status==2)~strata(pneu),data=sir.adm); 
# plot(cc);


ss=survfit(Surv(time,status==1) ~pneu,data=sir.adm);
# plot(ss,fun="cumhaz");
kmplot(ss,fun="cumhaz", main="discharged");

ss=survfit(Surv(time,status==2) ~pneu,data=sir.adm);
# plot(ss,fun="cumhaz"); 
kmplot(ss,fun="cumhaz", main="death");

par(mfrow=c(1,2))
aj <- prodlim(Hist(time, status)~pneu, data=sir.adm)
plot(aj, cause=1, plot.main="Discharge", col=c("black","red"),
     legend.legend=c("no pneu","yes pneu"))
plot(aj, cause=2, plot.main="DEATH", col=c("black","red"),
     legend.legend=c("no pneu","yes pneu"))

# Looking at the plots (ignore sampling variability and confidence intervals for now), consider.

#     Does pneumonia have an impact on the mortality rate? 
#       cirka no for DEATH
#     What about discharge? Do patients with pneumonia at admission stay longer at the unit?
#       chance of being discharged is much higher if you don't have pneumonia. 
#       it would appear that they stay longer, since rate of DEATH is same, annd
#       rate of discharge is lower => you must still be at the ward
#     Based on your answers to a-b above, would you expect to see more or fewer patients with 
#     pneumonia die on unit than patients without pneumonia?
#           More, since the absolute risk is higher

# 2. How many percent were not still at the intensive care unit after 50 days? 
# Draw the curve for the overall probability of leaving the unit (composite event death or discharge).
par(mfrow=c(1,1))
ahh <- prodlim(Hist(time, status!=0)~pneu, data=sir.adm, type="risk")
plot(ahh)
summary(ahh, time=50)

# 3. Draw the 1-KM (Kaplan-Meier) cuves for on-unit death censoring for discharge. 
#How many percent died on-unit within 50 days according to this method?
par1 <- prodlim(Surv(time, status==2)~pneu, data=sir.adm, type="risk")
summary(par1, time=50)
plot(par1)

# 4. Redo the analysis from question 3 but for live discharge from unit. 
# What is the sum of the two curves (for example at 50 days) compared to the overall 
# probability of leaving the unit in question 2?
par2 <- prodlim(Surv(time, status==1)~pneu, data=sir.adm, type="risk")
summary(par2, time=50)
plot(par2)

# 5. Estimate the cumulative incidence of death and discharge for patients with 
# and without pneumonia by the Aalen-Johansen product limit estimator. 
# Plot the cumulative incidence curves. Do they agree with your expectations from question 1?
aj1 <- prodlim(Hist(time, status)~factor(pneu), data=sir.adm)
par(mfrow=c(1,2))
plot(aj1, plot.main="Discharge", cause=1)
plot(aj1, plot.main="Death", cause=2)


# 6. What is the sum of the two curves (for example at 50 days) compared to 
# the overall probability of leaving the unit in question 2 and question 4?
summary(aj1, cause=1, times=50)
summary(aj1, cause=2, times=50)

# 7. Discuss the merits of the two types of plots, the cumulative hazards and cumulative incidence. 
# Which one do you prefer? Which gives a better description of the data?
# i prefer the cumulative incident plots, because it takes into account the total population
# however, it doesn't say anything about the rates of events

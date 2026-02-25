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

library(survival)
library(prodlim)

# Load the data from the file bissau.csv and store it in a data frame with the name bissau by the command
bissau <- read.csv2("http://biostat.ku.dk/frank/data/bissau.csv")

# Remember to load the survival and prodlim packages with library(survival) and library(prodlim).
# 1. Make a Kaplan-Meier plot for children with (DTP>0) and without (DTP=0).
# Use time on study as time scale (time from first visit).
# Calculate the survival probability at time 180 days after first visit.


# 2.Test if having received a DTP dose is good by a Cox model. Write a conclusion sentence.
# Also test the association between DTP any dose and mortalty with a log-rank test.


# 3. Are children vaccinated with DTP also vaccinated with BCG?
# Make a 2x2 table of BCG and DTP (any dose) to check this.



# 4. Perform a Cox regression to assess the association between DTP any dose and mortality adjusted for BCG.
# Why is this a good idea? Write a conclusion sentence.


# 5. Redo the analysis from the previous question adjusting for the age of the children by using age as time-scale.

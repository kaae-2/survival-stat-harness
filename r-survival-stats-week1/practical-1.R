options(repos = c(CRAN = "https://cloud.r-project.org"))

update.packages(ask = FALSE, checkBuilt = TRUE)

required_packages <- c("survival", "prodlim", "timereg")
missing_packages <- setdiff(required_packages, rownames(installed.packages()))

if (length(missing_packages) > 0) {
  install.packages(missing_packages)
}

invisible(lapply(required_packages, library, character.only = TRUE))

melanoma <- read.csv2("http://biostat.ku.dk/frank/data/melanoma.csv")
melanoma$years <- melanoma$days/365.24
melanoma$sex <- factor(melanoma$sex, labels=c("female","male"))

# A:
# A1. What is the time origin (“time zero”)?
# - Time of operation?
# A2. Are the survival times censored (if yes, by what and what kind of censoring?)
# End of 1977 (right censoring?) Administrative censoring
# A3. Are the survival times truncated (if yes, by what and what kind of truncaton?)
# patients are only included when they get an operation?
# consider inclusion bias
# A4. Are there competing risks (if yes, what are the possible endpoints)?
# death by other causes? relapse? 


melanoma$dead <- (melanoma$status != 2) # status not equal to 2 means dead from any cause


# B1. Limit the follow-up period to five years after surgery 
# and define the binary outcome “died before 5 years” (yes/no). 
# Make the 2x2 table relating this outcome to sex. 
# How do the censored patients appear in this table?
# Try to understand the following code snippet. 
# factor is R-language for a categorical variable (see the R help ?factor).
# Try removing useNA="ifany", what happens?
# - they are ignored - use t

melanoma$dead5yrs <- (melanoma$dead & melanoma$years<=5)
melanoma$dead5yrs[!melanoma$dead & melanoma$years<=5] <- NA
# use useNA="ifany" first,even if not used later to ensure we know how many are excluded
table(sex=melanoma$sex, dead5yrs=melanoma$dead5yrs, useNA="ifany") 
table(sex=melanoma$sex, dead5yrs=melanoma$dead5yrs)


# B2. In order to distinguish between survivors and censored patients, 
# make instead an outcome variable with three categories: 
# - “Dead before 5 years”, 
# - “Censored before 5 years” and 
# - “Seen alive after 5 years”, 
# again for the first five years after surgery. 
# Make the corresponding 2x3 table and comment.
# - done:
melanoma$status5yrs <- melanoma$dead
melanoma$status5yrs[melanoma$years>5] <- 2
melanoma$status5yrs <- factor(melanoma$status5yrs,labels=c("Censored before 5 years",
                              "Dead before 5 years","Seen alive after 5 years"))
table(sex=melanoma$sex, status5yrs=melanoma$status5yrs)

# B3. Would it be sensible to make any formal comparisons between treatment 
# based on any one of these variables from the two questions above?
# - 
# B4. Why cannot the mean survival time be calculated from the current data set 
# (at least not without imposing assumptions not verifiable from the data)?
# - 
# Do you expect the mean survival time among those who were observed to die,
# - we cannot calculate mean, because people are not all dead?
mean(melanoma$years[melanoma$dead])
# B5. What if we (somehow) knew that censoring was independent of both sex 
# (i.e. distributed equally in for males and females) and death
#  - we would have to look at the data, but then we might be able to conclude that 
#    there is a difference in outcomes.

# C1. Make Kaplan-Meier plots for male and female patients. 
# How can you see the timing of the deaths in the curves? 
# What about the censored individuals? 
###  smart survival function with automatic labels and risk set
out=prodlim(Surv(days/365,status!=2) ~sex ,data=melanoma)
summary(out)
plot(out)
title(main="Melanoma data")

out=survfit(Surv(days/365,status!=2) ~sex, data=melanoma)
summary(out)
plot(out)
kmplot(out) ### adds labels to plot with timereg package
title(main="Melanoma data")

# For the latter you may want to add the argument marktime=TRUE.
# drops in curve is death, marks on line is censoring
# C2. Do males or females have a better survival prognosis?
# - females
out=survdiff(Surv(days/365,status!=2) ~sex ,data=melanoma)
out
# C3. What is the probability (with 95% confidence interval) to still be alive 5 years after surgery?
# 0.79 for females
# 0.63 for males
# - assuming independence in censoring

out=survdiff(Surv(days/365,status!=2) ~sex+strata(ulc), data=melanoma)
out
# C4. Use a log-rank test to formally test the hypothesis that the survival is the same for males and females.
# - it is higher than "very low" so we shouldn't discard the null hypothesis
# C5. Plot the cumulative hazard for males and females, respectively. 
# Interpret and conclude based on the plots. 
# Relate the plot to the corresponding survival curves.

outna=survfit(Surv(days/365,status!=2)~sex, data=melanoma)
kmplot(outna,fun="cumhaz")


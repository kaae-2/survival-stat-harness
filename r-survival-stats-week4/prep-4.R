# Prep exercises for this week
# TODO: add your prep tasks here.
#Malignant melanoma

# In this exercise we will look at the malignant melanoma data from Odense University Hospital again. In the
# period 1962-77, 205 patients had their tumour removed and were followed until 1977. At the end of 1977:

# • 57 died of malignant melanoma (status=1)
# • 134 were still alive (status=2)
# • 14 died from other causes (status=3)

# The variables in the data set are
# no 	Patient code
# status 	Survival status. 1: dead from melanoma, 2: alive, 3: dead from other cause
# days 	Survival time
# ulc 	Ulceration, 1: present, 0: absent
# thick 	Tumour thickness (1/100 mm)
# sex 	0: female, 1: male

# Purpose: Study effect of sex, age, thickness of tumor and ulceration on survival.

# Load the data from the R package timereg and start by recode the censorings as 0
library(prodlim);
library(survival);
library(mets);
library(cmprsk);
library(riskRegression);


data(melanoma);
melanoma=dtransform(melanoma,status=0,status==2)
#
# Q1-1  With the melanoma data we shall consider the probability of dying from malignant melanoma (mm) or dying from other causes.
# Start by doing a simple cumulative incidence estimator (prodlim) of the 4 groups by sex and ulc.
# Plot and do summaries. In particular what are the probabilities of dying from mm at 2000 days.
prod <- prodlim(Hist(days, status) ~ulc+sex, data=melanoma)
plot(prod)
summary(prod, times=2000, cause=1)

# the risk of dying from mm is  
# 9.5% if female-no_ulc
# 11.9% if male-no_ulc
# 31.9% if female-ulc
# 46.6% if male-ulc


# Q1-2 Do a Grays test for the equivalence of the cumulative incidence for mm (cheat sheet section V section 3).
 cifa=with(melanoma,cuminc(days,status,interaction(ulc,sex)));
 print(cifa); 
 par(mfrow=c(1,1))
 plot(cifa);
 
  # There is very strong evidence that the cumulative incidence of melanoma death differs between the 4 ulc × sex group

# Q1-3 Take out summaries at 2000 days for the non-parametric estimator.
tp <- timepoints(cifa, times=2000)
 tp
 # 0.0 = 9.58%
 # 1.0 = 31.9%
 # 0.1 = 11.8%
 # 1.1 = 46.6%
 
# Q1-4 The extra challenge is to estimate and compare the cumulative incidence via cause-specific cox models
# with sex+ulc and compare to the estimates from the simple non-parametric estimator from prodlim.
# Specifically, you may do cox-regressions for each cause with for example the covariates sex, ulc, and their interaction,
# and combine these estimates into cumulative incidence predictions. See the cheat-sheet with the CSC function sextion IX.
# Plot the predictions  for the four groups based based on CSC
# (Hint to get a new data.frame with the four groups expand.grid(ulc=0:1,sex=0:1) will do the thing).
 nd <- expand.grid(ulc = 0:1, sex = 0:1)
 
 fit <- CSC(Hist(days, status) ~ ulc * sex, data = melanoma)
 
 tt <- seq(0, max(melanoma$days), by = 1)
 
 p1 <- predictRisk(fit, newdata = nd, cause = 1, times = tt)

 matplot(tt, t(p1), type = "s")
 
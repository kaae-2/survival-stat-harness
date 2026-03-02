#In this exercise we will look at the malignant melanoma data from Odense University Hospital again. In the
# period 1962-77, 205 patients had their tumour removed and were followed until 1977. At the end of 1977:

# • 57 died of malignant melanoma (status=1)
# • 134 were still alive (status=2)
# • 14 died from other causes (status=3)

# The variables in the data set are
# no 	Patient code
# status 	Survival status. 1: dead from melanoma, 2: alive, 3: dead from other cause
# days 	Survival time
# ulc 	Ulceration, 1: present, 2: absent
# thick 	Tumour thickness (1/100 mm)
# sex 	0: female, 1: male

# Purpose: Study effect of sex, age, thickness of tumor and ulceration on survival.

# Load the data from the file melanoma.csv and store it in a data frame with the name melanoma by the
# command
library(survival)
library(prodlim)
library(mets)

melanoma <- read.csv2("http://biostat.ku.dk/frank/data/melanoma.csv")

# We will study death due to malignant melanoma or other cause (overall survival). For this purpose make a
# new status variable which is 0 if the patient was still alive when the study closed in 1977 and 1 if the patient
# died from any cause.
melanoma$died <- 1*(melanoma$status!=2)
melanoma$ulcpres <- 1*(melanoma$ulc==1)

# Home exercise 2: Cox II Q1

# Q1-1 Draw a complementary log log (cloglog) survival curve to assess the proportional hazards assumption for the ulceration covariate.
# Interpret the plot, Is there a problem with proportionality?
cox1=phreg(Surv(days, died)~strata(ulcpres), data=melanoma)
plot(cox1, log='xy')

# I think so ? they cross twice

# Q1-2A-1 Test the proportional hazards assumption by testing for an interaction ulceration by f(time), where f(t)=I(t>2000).
# That is can you reject that the hazard ratio for ulceration differs for time t<= 2000 days and t>2000 days?
cox=coxph(Surv(days, died)~ulcpres + tt(ulcpres), data=melanoma, tt = function(x,t, ...) x * (t > 2000))
summary(cox, times=2000)
# p-value for tt(ulc) is lower than 0.033, the hazard ratio changes over time

# Q1-2A-2 Using the above model, what is the effect of ulceration for t<=2000, and what is the effect when t>2000.
summary(cox)
# the effect before 2000 is 4.8486
# the effect after is exp(1.5787-1.3034) = 1.316926
exp(1.5787-1.3034)

# Q1-3A Test the proportional hazards assumption by testing for an interaction ulceration by f(time), where f(t)=log(t).
# Can you reject that the hazard ratio for ulceration is constant over time?
coxlog=coxph(Surv(days, died) ~ ulc + tt(ulc), data=melanoma, tt = function(x, t, ...) x* log(t))
summary(coxlog)

# given p-value = 0.172 so we cannot reject the null hypothesis

# Q2-1 Fit a Cox model to evaluate the effect of tumour thickness of the hazard. Fit a model
# with both thick and log2(thick). Do you judge that the model be reduced to only include one
# of them without loss of fit? Or are both necessary? If, yes, which one?

melanoma$logthick <- log2(melanoma$thick)
coxboththick=coxph(Surv(days, died) ~ thick + logthick, data=melanoma)
coxthick=coxph(Surv(days, died) ~ thick, data=melanoma)
coxlogthick=coxph(Surv(days, died) ~ logthick, data=melanoma)
summary(coxboththick)
summary(coxthick)
summary(coxlogthick)
# we can ignore thick given the high p-value?


# Q2-2A What is the hazard ratio comparing patients where one group has tumours of twice the
# size of the other group? Compare thicker tumours to thinner tumours.

#about 1.628 ?


# Q2-3A Fit a Cox models for ulceration and base 2 logarithm of tumour thickness, i.e. log2(thick).
# Compare a simple unadjusted model for ulceration to an additive model with ulceration and
# log2(thick). Does controlling for log2(thick) have an impact on the conclusion for the importance of ulceration?
coxulclogthick = coxph(Surv(days, died) ~ulc + logthick, data=melanoma)
coxulc = coxph(Surv(days, died) ~ ulc, data=melanoma)
summary(coxulclogthick)
summary(coxulc)

#yes, the ulceration is less important when controlling for the thickness?

# Q2-4 Fit a model with interaction between ulceration and log2-transformed tumour thickness.
# Test if you can reject that there is an effect modification (interaction).
fit1=coxph(Surv(days, died) ~ulc+logthick + ulc*logthick, data=melanoma)
fit2 = coxph(Surv(days, died) ~ulc+logthick, data=melanoma)
anova(fit1,fit2)

#Given the p-value is over 5% it does not suggest we can reject the null-hypothesis

# Q2-5A Still in the interaction model, what is the hazard ratio comparing patients with ulceration
# where one group has tumours of twice the size of the other group?

summary(fit1)
# then it is exp(0.32606 + 0.005242)
exp(0.32606 + 0.05242)
# I would say it's  1.460064

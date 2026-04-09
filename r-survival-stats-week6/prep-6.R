# Malignant melanoma
# In this exercise we will look at the malignant melanoma data from Odense University Hospital again. In the period 1962-77, 205 patients had their tumour removed and were followed until 1977. At the end of 1977:
# • 57 died of malignant melanoma (status=1) • 134 were still alive (status=2) • 14 died from other causes (status=3)
# The variables in the data set are
#     no Patient code
#     status Survival status. 1: dead from melanoma, 2: alive, 3: dead from other cause
#     days Survival time
#     ulc Ulceration, 1: present, 0: absent
#     thick Tumour thickness (1/100 mm)
#     sex 0: female, 1: male
# Purpose: Study effect of sex, age, thickness of tumor and ulceration on survival.
# Reading data and recoding the censorings as 0 and making dead variable.
library(mets)
library(survival)

data(melanoma)

melanoma <- dtransform(melanoma,status=0,status==2)
melanoma <- dtransform(melanoma,dead=status!=0)
head(melanoma)

melanoma  <- dfactor(melanoma,fulc~ulc,labels=c("ulc-not-present","ulc-present"))
dtable(melanoma,~fulc+ulc)


# Answer the following questions and then answer the questions in the quiz

# 1) Considering overall death, compute how many days you loose up to 5 years, by computing the restricted residual mean at 5 years in the two groups. How many days do you loose more if you have ulceration.
# 5 years in days
tau <- 5 * 365

# Kaplan-Meier by ulceration group
fit <- survfit(Surv(days, dead) ~ fulc, data = melanoma)

# Restricted mean survival time up to 5 years
tab <- summary(fit, rmean = tau)$table

# Results: RMST and days lost up to 5 years
res <- data.frame(
  group = rownames(tab),
  rmst_5y_days = tab[, "rmean"],
  days_lost_5y = tau - tab[, "rmean"]
)

print(res)
460-87
# Extra days lost if ulceration is present
extra_days_lost <- with(
  res,
  days_lost_5y[group == "fulc=ulc-present"] -
    days_lost_5y[group == "fulc=ulc-not-present"]
)

cat("Extra days lost with ulceration up to 5 years:",
    round(extra_days_lost, 1), "\n")


# 2) We now wish to see how the "time lost" decompose on the two causes, malignant-melanoma death and other deaths. Compute how much of the time lost is due to malignant-melanoma deaths. What percentage of the time lost for the ulceration group is due to malignant melanoma deaths ?
# Restricted mean time lost due to malignant melanoma death
yl_melanoma <- resmeanIPCW(
  Event(days, status) ~ -1 + fulc,
  data = melanoma,
  time = tau,
  cause = 1,
  cens.model = ~ strata(fulc),
  model = "lin"
)

# Restricted mean time lost due to other deaths
yl_other <- resmeanIPCW(
  Event(days, status) ~ -1 + fulc,
  data = melanoma,
  time = tau,
  cause = 3,
  cens.model = ~ strata(fulc),
  model = "lin"
)

# Extract estimates
e1 <- estimate(yl_melanoma)
e2 <- estimate(yl_other)

res <- data.frame(
  group = names(coef(e1)),
  lost_due_to_melanoma = as.numeric(coef(e1)),
  lost_due_to_other = as.numeric(coef(e2))
)


res$total_time_lost <- res$lost_due_to_melanoma + res$lost_due_to_other
res$pct_due_to_melanoma <- 100 * res$lost_due_to_melanoma / res$total_time_lost

print(res)

# Answer specifically for the ulceration group
ulc_row <- grep("ulc-present", res$group)

cat(
  "Time lost due to malignant melanoma death in the ulceration group:",
  round(res$lost_due_to_melanoma[ulc_row], 1), "days\n"
)

cat(
  "Percentage of total time lost in the ulceration group due to malignant melanoma death:",
  round(res$pct_due_to_melanoma[ulc_row], 1), "%\n"
)

# 3) We finally wish to express how the ulceration changes the survival at 5 years. We will do this by computing a standardized survival estimate for ulceration. Do this by using a Cox model adjusting for thick(ness) and sex. You may use the surivalG function for this. What is the standardized risk difference due to ulceration ?

melanoma$sex <- factor(
  melanoma$sex,
  levels = c(0, 1),
  labels = c("female", "male")
)
# Cox model adjusted for thickness and sex
fit <- phreg(Surv(days, dead) ~ fulc + thick + sex, data = melanoma)

# Standardized 5-year survival by ulceration
g5 <- survivalG(fit, data = melanoma, time = tau)

# This prints the answer table
summary(g5)

# the standardized risk difference due to ulceration at 5 years is:
#   
# 0.2053, or 20.5 percentage points.
# 
# Interpretation:
#   
# Standardized 5-year risk with no ulceration: 0.1727 = 17.3%
# Standardized 5-year risk with ulceration: 0.3780 = 37.8%
# Risk difference: 0.3780 - 0.1727 = 0.2053

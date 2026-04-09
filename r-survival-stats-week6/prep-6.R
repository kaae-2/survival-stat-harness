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

data(melanoma)

melanoma <- dtransform(melanoma,status=0,status==2)
melanoma <- dtransform(melanoma,dead=status!=0)
head(melanoma)

melanoma  <- dfactor(melanoma,fulc~ulc,labels=c("ulc-not-present","ulc-present"))
dtable(melanoma,~fulc+ulc)
# Answer the following questions and then answer the questions in the quiz

# 1) Considering overall death, compute how many days you loose up to 5 years, by computing the restricted residual mean at 5 years in the two groups. How many days do you loose more if you have ulceration.

# 2) We now wish to see how the "time lost" decompose on the two causes, malignant-melanoma death and other deaths. Compute how much of the time lost is due to malignant-melanoma deaths. What percentage of the time lost for the ulceration group is due to malignant melanoma deaths ?

# 3) We finally wish to express how the ulceration changes the survival at 5 years. We will do this by computing a standardized survival estimate for ulceration. Do this by using a Cox model adjusting for thick(ness) and sex. You may use the surivalG function for this. What is the standardized risk difference due to ulceration ?

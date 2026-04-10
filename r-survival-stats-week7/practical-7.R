
# Setting up the bladder cancer data

#    Data on recurrences of bladder cancer, used by many people to demonstrate methodology for recurrent event modelling.
#    Bladder1 is the full data set from the study. It contains all three treatment arms and all recurrences for 118 subjects; the maximum observed number of recurrences is 9.

# The data (bladder1)

#    id: Patient id
#    treatment: Placebo, pyridoxine (vitamin B6), or thiotepa
#    number: Initial number of tumours (8=8 or more)
#    size: Size (cm) of largest initial tumour
#    recur: Number of recurrences
#    start,stop: The start and end time of each time interval
#    status: End of interval code, 0=censored, 1=recurrence, 2=death from bladder disease, 3=death other/unknown cause
#    rtumor: Number of tumors found at the time of a recurrence
#    rsize: Size of largest tumor at a recurrence
#    enum: Event number (observation number within patient)

library(mets); library(survival);
head(bladder)
dtable(bladder1,~status+treatment,level=1)
bladder1$event <- 1*(bladder1$status!=0)
# ## break ties
blad <- tie.breaker(bladder1,start="start",stop="stop",status="event",id="id")
# ## avoid "0" risk data point
blad <- dtransform(blad,stop=0.5,id==49)
## make lag-count
blad <- count.history(blad,status="event",types=1)
dlist(blad,.~id| id %in% c(26,27,91))

# R commands for today

# Some relevant R commands from the cheat-sheet

# Marginal mean:

c1 <- phreg(Surv(start,stop,status==1)~+strata(treatment)+cluster(id),blad)
c2 <- phreg(Surv(start,stop,status>=2)~+strata(treatment)+cluster(id),blad)
plot(c2,se=1); title(main="Hazard of Death")
m1 <- recurrentMarginal(c1,c2)
bplot(m1,se=1);  title(main="Average events")

# Ghosh-Lin model (ipcw estimation)

gl <- recreg(
  Event(start,stop,status)~+treatment+cluster(id),
  blad,
  death.code=c(2,3),
  cause=1,
  cens.code=0
)
summary(gl)

# Andersen-Gill model

c1 <- phreg(Surv(start,stop,status==1)~+treatment+Count1+cluster(id),blad)

# Q1. The intensities Andersen-Gill

#    First we study the effects of treatment on the hazard of event among those 
#    still at risk while correcting for size and number of initial tumor(s).
#    Also to have a look at how these covariates affects who are alive.
#    Skip goodness of fit considerations.

# Q2 Marginal Mean number of events

#    Now estimate the marginal mean given treatment groups. 
#    Remember carefully what we are estimating. Make the plot. 
#    Check when the mice are dying in the different groups.
#    Estimate now the proportional effect of treatment (Ghosh-Lin) 
#    on the marginal mean adjusted and un-adjusted. 
#    Adjusting for initial size and number.

# Q4 Cox power calculation

#    I would suggest to do this via the homepages for power-calculations, 
#    even though I here show an R based solution.
#    With hazard ratio 0.7, power=0.9 and the significance level 5%. 
#    What is the number of event needed.
#    If 0.3 of the subjects in the baseline group gets the event in 3 years, 
#    what is the number of subjects needed in two equally sized groups.

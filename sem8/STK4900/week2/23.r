# The data are from a old study of children with leukemia.
# The children were treated so that they had no symptoms (one say that they are in remission).
# They were then given an active treatment ("Drug") to keep the symptoms from coming back or placebo ("Control").
# Time until the symptoms came back (called relapse) were measured.
# For some children the symptoms had not come back at the end of the study. These give rise to censored observations.
 
# Read the data into a dataframe and inspect the data:
gehan=read.table("http://www.uio.no/studier/emner/matnat/math/STK4900/v11/gehan.dat", header=T)
gehan
 
# Check that the data correspond to those given in the exercise.
# "time" is time in weeks to relapse or censoring
# "cens" is 1 for relapse and 0 for censoring
# "treat" is 1 for "Control" and 2 for "Drug"

#         It is correct

# Attach the R-library for survival analysis:
library(survival)
 
 
# QUESTION 1
 
# Compute Kaplan-Meier estimates for the two groups (without confidence intervals)
fit.1=survfit(Surv(time,cens)~treat, conf.type="none", data=gehan)
summary(fit.1)
 
# Make sure you understand what the output tells!
 
 
# Plot the Kaplan-Meier estimates:
plot(fit.1,lty=1:2)
 
# Interpret the plots.
#- The plots show the survival chance for treated and untreated patients. 
#- It is clear that the treated patients have a longer survival rate than the untreated group.

# Read from the plots (approximately) what is the median time to relapse for the two groups
#- Untreated patients survival chance drop below .5 between 8 and 10, so we say 8.
#- Treated patients drop below .5 between 22 and 23, so we say 22.

# Check the results by giving the command print(fit.1)
print(fit.1)
#-Corresponds to what I read from the plot
 
# QUESTION 2
 
# Compute Kaplan-Meier estimates for the two groups with confidence intervals
# (the default confidence interval in R is not a good choice, so we make an explicit choice of the type of confidence interval)
fit.2=survfit(Surv(time,cens)~treat, conf.type="plain",data=gehan)
summary(fit.2)
plot(fit.2, conf.int=T, lty=1:2)
 
# Interpret the output and the plot.
#-The plot shows the confidence interval around the model estimates.
#-We note that there is some overlap between the treated and untreated patients in the plot.
 
# QUESTION 3
 
# Log-rank test for difference between the groups:
survdiff(Surv(time,cens)~treat, data=gehan)

# What does the output tell you?
#-This is basically a deviance test which ends up with a p-value=4.17e-5 indicating a
#-significant effect of the treatment.
 
# QUESTION 4
 
# Cox regresion with treatment group as covariate:
fit.4=coxph(Surv(time,cens)~factor(treat), data=gehan)
summary(fit.4)
 
# Interpret the results!
#-We get output indicating a huge probability that the treatment is significant to the survival time.

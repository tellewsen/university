# Read the data into a dataframe:
melanom=read.table("http://www.uio.no/studier/emner/matnat/math/STK4900/v11/melanoma.dat",header=T)
 
# Attach the R-library for survival analysis (if you have not already done so):
library(survival)
 
 
# QUESTION a
 
# Commands:
fit.sex=survfit(Surv(lifetime,status==1)~sex, data=melanom)
plot(fit.sex, lty=1:2, mark.time=F)
#-Plot: Women seem to surive longer than men on average
survdiff(Surv(lifetime,status==1)~sex, data=melanom)
#-There is a significance in survival from gender.
 
# QUESTIONS b and c are done by similar commands
#b
fit.grthick=survfit(Surv(lifetime,status==1)~grthick, data=melanom)
plot(fit.grthick, lty=1:2, mark.time=F)
#-Plot: Thickness group seems to have an effect
survdiff(Surv(lifetime,status==1)~grthick, data=melanom)
#-Confirmed by a very low p-value

#c
fit.ulcer=survfit(Surv(lifetime,status==1)~ulcer, data=melanom)
plot(fit.ulcer, lty=1:2, mark.time=F)
#-Plot: Seems to be a difference whether or not there is ulcerations as well
survdiff(Surv(lifetime,status==1)~ulcer, data=melanom)
#-low p-value confirms the suspission

#D:
# Commands for sex
coxfit.sex=coxph(Surv(lifetime,status==1)~factor(sex),data=melanom)
summary(coxfit.sex)
#Seems to be a small effect
coxfit.thickn=coxph(Surv(lifetime,status==1)~thickn,data=melanom)
summary(coxfit.thickn)
#Seem to be a strong significance
coxfit.ulcer=coxph(Surv(lifetime,status==1)~factor(ulcer),data=melanom)
summary(coxfit.ulcer)
#Strong effect here as well
 
#Test whether to use thickn,logthickn, or grthick
coxfit.grthick=coxph(Surv(lifetime,status==1)~factor(grthick),data=melanom)
summary(coxfit.grthick)

coxfit.logthick=coxph(Surv(lifetime,status==1)~logthick,data=melanom)
summary(coxfit.logthick)

#It seems that the thickness variable gets the lowest log rank test p value so I choose to go with that one.

#E:
coxfit.many=coxph(Surv(lifetime,status==1)~logthick+thickn+factor(grthick)+sex+ulcer,data=melanom)
summary(coxfit.many)
#We get a model with failry low significance to most variables.
#There are two ways to remove terms here. The worst one is one of the grthick factors, but the other of them is much more significant, so we could remove logthick instead since that is the next worst one. 
#Depeneding on which one we remove we eliminate down to the following two models.

coxfit.many1=coxph(Surv(lifetime,status==1)~factor(grthick)+ulcer,data=melanom)
coxfit.many2=coxph(Surv(lifetime,status==1)~logthick+ulcer,data=melanom)
summary(coxfit.many1)
summary(coxfit.many2)
#These have very close values for all three tests, but since the wald test returns one order of magnitude lower p value for the last model I choose to make this the final model.


cirrhosis <- read.table("cirrhosis.txt",sep="\t",header=TRUE)
cirrhosis
summary(cirrhosis)
png('boxplot.png')
boxplot(cirrhosis)
dev.off()
# Attach the R-library for survival analysis
library(survival)

fit.1=survfit(Surv(time)~treat, conf.type="none", data=cirrhosis)
png(filename='treat.png')
plot(fit.1,lty=1:2,xlab='time[days]')
dev.off()
#dev.copy(png,'treat.png')
#dev.off()
#No noticable effect from the plot

fit.2=survfit(Surv(time)~sex, conf.type="none", data=cirrhosis)
png('sex.png')
plot(fit.2,lty=1:2,xlab='time[days]')
dev.off()
#no noticable effect

fit.3=survfit(Surv(time)~asc, conf.type="none", data=cirrhosis)
png('asc.png')
plot(fit.3,lty=1:3,xlab='time[days]')
dev.off()
#noticable effect

fit.4=survfit(Surv(time)~factor(agegr), conf.type="none", data=cirrhosis)
png('agegr.png')
plot(fit.4),lty=1:3,xlab='time[days]')
dev.off()
#noticable effect

survdiff(Surv(time)~treat, data=cirrhosis)
survdiff(Surv(time)~sex, data=cirrhosis)
survdiff(Surv(time)~asc, data=cirrhosis)
survdiff(Surv(time)~agegr, data=cirrhosis)

#Treatment not significant!
coxfit.1=coxph(Surv(time,status==1)~treat+sex+asc+age,data=cirrhosis)
summary(coxfit.1)
coxfit.2=coxph(Surv(time,status==1)~sex+asc+age,data=cirrhosis)
summary(coxfit.2)
#Use this one further and include interactions
coxfit.3=coxph(Surv(time,status==1)~sex+asc+age+sex:age+sex:asc+asc:age,data=cirrhosis)
summary(coxfit.3)
#remove interaction age:asc
coxfit.4=coxph(Surv(time,status==1)~sex+asc+age+sex:age+sex:asc,data=cirrhosis)
summary(coxfit.4)
#remove sex:asc
coxfit.5=coxph(Surv(time,status==1)~sex+asc+age+sex:age,data=cirrhosis)
summary(coxfit.5)
#remove sex
coxfit.6=coxph(Surv(time,status==1)~asc+age+sex:age,data=cirrhosis)
summary(coxfit.6)

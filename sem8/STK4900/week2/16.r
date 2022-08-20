#Test change in support for Hoyre
n1=969
y1=310
p1=y1/n1
se1=sqrt(p1*(1-p1)/n1)
n2=968
y2=303
p2=y2/n2
se2=sqrt(p2*(1-p2)/n2)
se=sqrt(se1^2+se2^2)
change=p1-p2
margin=1.96*se
lower=change-margin
upper=change+margin
cbind(change,margin,lower,upper)

# change     margin       lower      upper
#[1,] 0.006900912 0.04142406 -0.03452315 0.04832497

#Note that the margin is larger than the change giving us a confidence interval(CI) for the change that includes zero. Thus one would say that the support has not changed.

#Next we test the null hypothesis that the support has not changed
p=(y1+y2)/(n1+n2)
se0=sqrt(p*(1-p)/n1+p*(1-p)/n2)
z=(p1-p2)/se0
pval=2*(1-pnorm(abs(z)))
cbind(z,pval)
#We get a large p-value for this indicating that there is a good chance that the support has not changed, supporting the conclusion from the CI.

#Test matrix function

hoyre=matrix(c(y1,y2,n1-y1,n2-y2),nrow=2)   # give the data for Høyre in a 2x2 table (cf. slide 10)
prop.test(hoyre,correct=F)
#z*z = chi-squared like it should be since a chi test of 2 proportions is the same as doing a z-test

#The next step is to repeat this for AP
n1=969
y1=263
p1=y1/n1
se1=sqrt(p1*(1-p1)/n1)
n2=968
y2=238
p2=y2/n2
se2=sqrt(p2*(1-p2)/n2)
se=sqrt(se1^2+se2^2)
change=p1-p2
margin=1.96*se
lower=change-margin
upper=change+margin
cbind(change,margin,lower,upper)

#Same conclusion as above. CI includes zero, proposing no change.

p=(y1+y2)/(n1+n2)
se0=sqrt(p*(1-p)/n1+p*(1-p)/n2)
z=(p1-p2)/se0
pval=2*(1-pnorm(abs(z)))
cbind(z,pval)

#Get a lower p-value now but it is still large enough to not be negligible. Probability of no change is good.

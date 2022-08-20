#load data
wcgs=read.table("http://www.uio.no/studier/emner/matnat/math/STK4900/v11/wcgs.txt",sep="\t", header=T,na.strings=".")

#Make model and check corrolation
fit.smoke=glm(chd69~smoke,data=wcgs,family=binomial)
summary(fit.smoke)
#Signifcant corrolation between smoking and chd!

#Calc CI and oddsratio and CI for oddsratio
#CI smoking reg cooeff:
#Output from summary says
#smoke         0.6299     0.1337    4.71 2.47e-06 ***
lower = .6299 - 1.96*0.1337
#lower = 0.367848
upper = .6299 + 1.96*0.1337
#upper = 0.891952	

expcoef=function(glmobj) 
{
regtab=summary(glmobj)$coef
expcoef=exp(regtab[,1])
lower=expcoef*exp(-1.96*regtab[,2])
upper=expcoef*exp(1.96*regtab[,2])
cbind(expcoef,lower,upper) 
}
 
expcoef(fit.smoke) 
#Output:
#               expcoef      lower      upper
#(Intercept) 0.06306306 0.05141861 0.07734457
#smoke       1.87735347 1.44451020 2.43989698


 # We then use logistic regression to study the effect of age (at entry to the study) for the risk of developing CHD:

fit.age=glm(chd69~age,data=wcgs,family=binomial)
summary(fit.age)

expcoef(fit.age)
#Output:
#                expcoef        lower       upper
#(Intercept) 0.002633304 0.0008972465 0.007728412
#age         1.077261913 1.0536601655 1.101392334

#Use 10 years instead
fit.age10=glm(chd69~I(age/10),data=wcgs,family=binomial)
summary(fit.age10)
#Everything looks the same
expcoef(fit.age10)
#We get a difference in the CI and expcoef for age and age/10

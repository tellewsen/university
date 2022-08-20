cancerdata=hers=read.table("http://www.uio.no/studier/emner/matnat/math/STK4900/v11/cancer.txt")
names(cancerdata)=c("age","cig","pyr","cancer")
cancerdata

# We first consider the model  E(Y) = n*exp{b0+b1*s+b2*a}, where
#      Y=number of cancer cases (=cancer),
#       n=number of person years (= pyr),
#       s=number of cigarettes smoked per day (=cig)
#       a = age in years (=age)
# We may write the model on the form  E(Y)= exp{1*log(n)+b0+b1*s+b2*a}.
# Note that  log(n)  appears as a sort of "covariate" where we know that the regression coefficient takes the value 1. This is called an OFFSET.
 
# We fit the model and look at the result::
cancerfit.1=glm(cancer~offset(log(pyr))+age+cig, data=cancerdata, family=poisson)
summary(cancerfit.1)
#Significant effect of age and cigarettes smoked.

#make excoef function
expcoef=function(glmobj) 
{
regtab=summary(glmobj)$coef
expcoef=exp(regtab[,1])
lower=expcoef*exp(-1.96*regtab[,2])
upper=expcoef*exp(1.96*regtab[,2])
cbind(expcoef,lower,upper) 
}

expcoef(cancerfit.1)
#We see a 7% increase in odds for each cigarette smoked, and 11% increase in odds for each year we age.


# We then look at a model with second order terms and interaction:
cancerfit.2=glm(cancer~offset(log(pyr))+ age+I(age^2)+cig+I(cig^2)+age:cig, data=cancerdata, family=poisson)
summary(cancerfit.3)

#Interaction term not significant, rest are significant, see output:

#Coefficients:
#              Estimate Std. Error z value Pr(>|z|)    
#(Intercept) -2.316e+01  3.291e+00  -7.039 1.94e-12 ***
#age          3.587e-01  9.952e-02   3.604 0.000313 ***
#I(age^2)    -1.880e-03  7.797e-04  -2.411 0.015918 *  
#cig          2.041e-01  6.521e-02   3.130 0.001745 ** 
#I(cig^2)    -1.998e-03  5.656e-04  -3.532 0.000412 ***
#age:cig     -7.090e-04  8.513e-04  -0.833 0.404985   

#G approx 22 with 3 df between this and the last model. This one much better. Still I dont want the interaction term.
#remove it and check then.
anova(cancerfit.2)
#Note the deviance(G)=0.703 with 1 degree of freedom. This corresponds to a p-value between 0.95 and 0.05 making the interaction term not significant. Thus we remove it.
cancerfit.3=glm(cancer~offset(log(pyr))+ age+I(age^2)+cig+I(cig^2), data=cancerdata, family=poisson)
#The rest of the terms have G >3.841 with 1df giving them all significance. THus these are kept making this the final model. 
#The estimates are the Beta values in or model P=pyr*exp(cancer+age+age^2+cig+cig^2).
#To find the increase by smoking for example one cigarette one calculates exp(cig) = exp(1.561e-01)=1.168.
#Thus for each cigarette you smoke the odds of getting cancer increase by 16%.
#Conclusion, either people do something that gives them cancer every time they smoke, or smoking causes cancer by itself. 

#Alternatively we could use intervals of age and number of cigarettes smoked.
#This could be given by the following model:
cancerfit.a=glm(cancer~offset(log(pyr))+factor(age)+factor(cig), data=cancerdata, family=poisson)
#From this model we see that age is not significant until one reaches the 52year group. One can also see that the significance of cigarettes is only weekly significant for less than 5.2 cigarettes per day. We should also note that the residual deviance is lower for this model.

#Doing a deviance test between cancerfit.3 and cancerfit.a gives G = approx 10 with 10 degrees of freedom. This corresponds to p-value between 0.95 and 0.05 indicating no significance in changing model.

#load data
wcgs=read.table("http://www.uio.no/studier/emner/matnat/math/STK4900/v11/wcgs.txt", sep="\t",header=T,na.strings=".")

#make table with smoke and chd
table(wcgs$smoke,wcgs$chd69)

#Output:
#       0    1
#  0 1554   98
#  1 1343  159

#use values 
p1=159/(1343+159)
p0=98/(1554+98)
RR=p1/p0
OR=(p1/(1-p1))/(p0/(1-p0))
cbind(p1,p0,RR,OR)

#Output:
#            p1         p0       RR       OR
#[1,] 0.1058589 0.05932203 1.784478 1.877353

#p1,p0 proportions of groups getting CHD.
#RR Increase in chance of getting CHD if smoker over non smoker

#OR Odds ratio. Ratio between the odds of getting, and not getting, chd for smoker and non smoker.

#reg model for smoke
fit.smoke=glm(chd69~smoke,data=wcgs,family=binomial)
print(fit.smoke)

#OUtput:
#Call:  glm(formula = chd69 ~ smoke, family = binomial, data = wcgs)
#
#Coefficients:
#(Intercept)        smoke  
#    -2.7636       0.6299  

#Degrees of Freedom: 3153 Total (i.e. Null);  3152 Residual
#Null Deviance:	    1781 
#Residual Deviance: 1758 	AIC: 1762

#We can use that the odds ratio corresponding to one increase in smoke corresponds to exp(0.6299) = 1.877423 fitting almost perfectly to the odds ratio calculated earlier

#Then age
fit.age=glm(chd69~age,data=wcgs,family=binomial)
print(fit.age)
#Output:
#Call:  glm(formula = chd69 ~ age, family = binomial, data = wcgs)
#
#Coefficients:
#(Intercept)          age  
#   -5.93952      0.07442  
#
#Degrees of Freedom: 3153 Total (i.e. Null);  3152 Residual
#Null Deviance:	    1781 
#Residual Deviance: 1738 	AIC: 1742

#Try the same as before to calculate odds ratio: 
exp(0.07442)# = 1.077259
#and 
exp(0.07442*10)#=2.104757

#Here we see that the odds are almost equal after 1 year, but the odds have increased to twice as much for the smokers after 10 years.


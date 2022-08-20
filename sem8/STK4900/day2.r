#ex 5

#a)

pef=c(494,395,516,434,476,413,442,433)
minipef =c(512,430,520,428,500,364,380,445)
boxplot(pef,minipef)


#b)

cor(pef,minipef)
# 0.8154886


#c)
cor.test(minipef,pef)

#	Pearson's product-moment correlation


#data:  minipef and pef
#t = 3.4513, df = 6, p-value = 0.01361
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# 0.2605298 0.9653948
#sample estimates:
#      cor 
#0.8154886 


#d)
fit =lm(minipef~pef)
summary(fit)

#Call:
#lm(formula = minipef ~ pef)
#
#Residuals:
#    Min      1Q  Median      3Q     Max 
#-57.625 -12.796   6.763  19.088  47.090 
#
#Coefficients:
#            Estimate Std. Error t value Pr(>|t|)  
#(Intercept) -76.9345   152.4750  -0.505   0.6319  
#pef           1.1642     0.3373   3.451   0.0136 *
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 36.9 on 6 degrees of freedom
#Multiple R-squared:  0.665,	Adjusted R-squared:  0.6092 
#F-statistic: 11.91 on 1 and 6 DF,  p-value: 0.01361

#We get a slope of 1.1642 indicating that the two measurements rise with almost the same amount for each step
#except that one of them starts with a value approx 77 lower than the other.


#e)
1.1642*sd(pef)/sd(minipef)
# 0.8155151

#The number we get out is almost identical to the pearson correlation. 



#########
#ex6
#a)
hers.sample=read.table("http://www.uio.no/studier/emner/matnat/math/STK4900/v11/hers.sample.txt",header=T)
plot(hers.sample$age,hers.sample$sbp)
#no clear visual correlation
#b)
hers.fit.b=lm(sbp~age,data=hers.sample)
summary(hers.fit.b)

#Call:
#lm(formula = sbp ~ age, data = hers.sample)
#
#Residuals:
#    Min      1Q  Median      3Q     Max 
#-43.550 -13.520  -2.431  12.578  62.736 
#
#Coefficients:
#            Estimate Std. Error t value Pr(>|t|)    
#(Intercept) 105.7130    12.4024   8.524 1.05e-15 ***
age           0.4405     0.1865   2.363   0.0188 *  
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 19.76 on 274 degrees of freedom
#Multiple R-squared:  0.01997,	Adjusted R-squared:  0.01639 
#F-statistic: 5.582 on 1 and 274 DF,  p-value: 0.01884

#There seems to be an increase in sbp as age increase. 
#It is not very large, but it is there.


#c)
hers.fit.c=lm(sbp~I(age-67),data=hers.sample)
summary(hers.fit.c)
#Call:
#lm(formula = sbp ~ I(age - 67), data = hers.sample)
#
#Residuals:
#    Min      1Q  Median      3Q     Max 
#-43.550 -13.520  -2.431  12.578  62.736 
#
#Coefficients:
#            Estimate Std. Error t value Pr(>|t|)    
#(Intercept) 135.2284     1.1985 112.829   <2e-16 ***
#I(age - 67)   0.4405     0.1865   2.363   0.0188 *  
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Residual standard error: 19.76 on 274 degrees of freedom
#Multiple R-squared:  0.01997,	Adjusted R-squared:  0.01639 
#F-statistic: 5.582 on 1 and 274 DF,  p-value: 0.01884


#Comparing the least square estimates from c and b reveal that the standard error of the intercept is much lower in the c case

#d)
#It is the same thing, only the \beta_0 value is multiplied by 10 in the last one.


#Ex6
#a)
n=25
rho=0.30
m=matrix(c(0,0),nrow=2)
S=matrix(c(1,rho,rho,1),nrow=2)
obs=mvrnorm(n,m,S)
x=obs[,1]
y=obs[,2]
cor(x,y)
plot(x,y)

#the pearson correlation varies. This is caused by the low number of values that are generated.

#b)
n=25
rho=0.60
m=matrix(c(0,0),nrow=2)
S=matrix(c(1,rho,rho,1),nrow=2)
obs=mvrnorm(n,m,S)
x=obs[,1]
y=obs[,2]
cor(x,y)
plot(x,y)

#The variation increases a lot , but there is a clear indication of a correlation.

n=25
rho=0.90
m=matrix(c(0,0),nrow=2)
S=matrix(c(1,rho,rho,1),nrow=2)
obs=mvrnorm(n,m,S)
x=obs[,1]
y=obs[,2]
cor(x,y)
plot(x,y)


#c)

n=100
rho=0.30
m=matrix(c(0,0),nrow=2)
S=matrix(c(1,rho,rho,1),nrow=2)
obs=mvrnorm(n,m,S)
x=obs[,1]
y=obs[,2]
cor(x,y)
plot(x,y)


n=100
rho=0.60
m=matrix(c(0,0),nrow=2)
S=matrix(c(1,rho,rho,1),nrow=2)
obs=mvrnorm(n,m,S)
x=obs[,1]
y=obs[,2]
cor(x,y)
plot(x,y)


n=100
rho=0.90
m=matrix(c(0,0),nrow=2)
S=matrix(c(1,rho,rho,1),nrow=2)
obs=mvrnorm(n,m,S)
x=obs[,1]
y=obs[,2]
cor(x,y)
plot(x,y)


n=400
rho=0.30
m=matrix(c(0,0),nrow=2)
S=matrix(c(1,rho,rho,1),nrow=2)
obs=mvrnorm(n,m,S)
x=obs[,1]
y=obs[,2]
cor(x,y)
plot(x,y)


n=400
rho=0.60
m=matrix(c(0,0),nrow=2)
S=matrix(c(1,rho,rho,1),nrow=2)
obs=mvrnorm(n,m,S)
x=obs[,1]
y=obs[,2]
cor(x,y)
plot(x,y)


n=400
rho=0.90
m=matrix(c(0,0),nrow=2)
S=matrix(c(1,rho,rho,1),nrow=2)
obs=mvrnorm(n,m,S)
x=obs[,1]
y=obs[,2]
cor(x,y)
plot(x,y)

#As expected increasing the number of samples brings the pearson correlation closer to the real rho
#As rho increases the correlation becomes more and more apparent visually as it should.

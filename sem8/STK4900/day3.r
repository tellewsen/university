#1)
cafe=read.table("http://www.uio.no/studier/emner/matnat/math/STK4900/v11/exer3_1.dat")
names(cafe)=c("no","sale")
cafe
attach(cafe) #makes no and sale available outside cafe
plot(no,sale)

#looks fairly linear

linfit=lm(sale~no)
linfit
abline(linfit)

#fits fairly well.

#Try second order
x=seq(0,7,0.1)
koef=lm(sale~no+I(no^2))$coef
koef
lines(x,koef[1]+koef[2]*x+koef[3]*x^2,lty=2)

#fits much better.
#I would choose the second order one

#2)
insurance=read.table("http://www.uio.no/studier/emner/matnat/math/STK4900/v11/exer3_2.dat")
names(insurance)=c("income","riskave","amount")
insurance
attach(insurance)
summary(insurance)
#summary tells us the distribution of the different variables

par(mfrow=c(1,2))
plot(income,amount)
plot(riskave,amount)
par(mfrow=c(1,1))
#plots reveal that there seams to be a linear correlation between income and amount
#I can not see any clear indication of a correlation between risk aversion and amount

cor(insurance)
#          income   riskave    amount
#income  1.0000000 0.2544991 0.9058552
#riskave 0.2544991 1.0000000 0.3988377
#amount  0.9058552 0.3988377 1.0000000

#the results of the correlation agree with the visual conclusion.
#there is however a 25% correlation between riskaversion and income.  

#income is clearly the major predictor in this case.

fit3=lm(amount~income+riskave)
summary(fit3)
#together incoem and riskave can be said to explain 85% of the behaviour of insurance
#with linear regression



#3)
hers=read.table("http://www.uio.no/studier/emner/matnat/math/STK4900/v11/hers.txt",sep="\t",header=T,na.strings=".")
hers.no=hers[hers$diabetes==0, ]

summary(hers.no$glucose[hers.no$exercise==0])
summary(hers.no$glucose[hers.no$exercise==1])
boxplot(hers.no$glucose~hers.no$exercise)

#The minimum is higher for those who don't exercise, both the median #and mean are higher for then non exercising people. The maxes are the #same. The boxplot agrees with this.

t.test(glucose~exercise, var.equal=T,data=hers.no)
#a t.test also agrees, giving a p-value 0.000113
#with a CI (.83, 2.55)
#exercise has an effect on glucose



insects=read.table("http://www.uio.no/studier/emner/matnat/math/STK4900/v11/insects.txt",header=T)

proportion<-insects$DEAD/insects$NUMBER
plot(insects$LOGDOSE, proportion, ylim=c(0,1),pch=16,ylab="Proportion dead",xlab="Log conc. rotenone")
#Get a nice plot that starts off linear but looks to be topping out at the end.


# Question b)

# Fit a logistic regression model with logdose as covariate and look at the result:
# (Note that since we have grouped data, not binary, the response has to be specified as
# cbind(y,n-y) where n=number of individuals in a group and y=number of "successes" in the group.)
fit<-glm(cbind(DEAD,NUMBER-DEAD)~LOGDOSE, data=insects,family=binomial)
summary(fit)
  
# Question c)
# Compute the probabilities obtained from the fitted model and include them in the plot from question a
pred.prop=predict(fit,type="response")
points(insects$LOGDOSE ,pred.prop)
# We may also draw a curve that describes the fitted logistic model:
logdose=seq(0.4,1,0.01)
new.doses=data.frame(LOGDOSE=logdose)
lines(logdose,predict(fit,new.doses,type="response"))

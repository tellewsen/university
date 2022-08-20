#load in data
no2data <- read.table("no2.txt",sep="\t",header=TRUE)

#make data available outside no2data and cirrhosis
attach(no2data)
plot(log.cars,no2)

#make a linear fit to the plot
linfit=lm(no2~log.cars)

#make the line appear on the plot
abline(linfit)

#study the linear fit
summary(linfit)

#make a boxplot
boxplot(no2,log.cars)

#Test the correlation between no2 and log.cars
cor.test(no2,log.cars)

#Plot all the data
plot(no2data)

#Correlations between all variables
cor(no2data)

#Make a terrible multiple regression model
multireg3 = lm(no2~log.cars+temp+wind.speed+hour.of.day+log(wind.speed)+log(hour.of.day) +I(log.cars^2)+I(temp^2)+I(wind.speed^2)+I(hour.of.day^2)+I(log(wind.speed)^2)+I(log(hour.of.day)^2))
#Has adjusted R^2=0.4983

#Remove obviously useless terms
multireg4 = lm(no2~log.cars+temp+wind.speed+hour.of.day+log(wind.speed)+log(hour.of.day) +I(log.cars^2)+I(temp^2)+I(log(wind.speed)^2)+I(log(hour.of.day)^2))
#Has adjusted R^2=0.5003

#Remove more useless terms
multireg5 = lm(no2~temp+wind.speed+hour.of.day+log(wind.speed)+log(hour.of.day) +I(log.cars^2)+I(log(wind.speed)^2)+I(log(hour.of.day)^2))
#Has adjusted R^2=0.5004

#Try removing even more terms
multireg6 = lm(no2~temp+wind.speed+log(wind.speed)+I(log.cars^2)+I(log(wind.speed)^2)+I(log(hour.of.day)^2))
#Has adjusted R^2=0.4968

#R has started going down. Can still remove more terms
multireg7 = lm(no2~temp+wind.speed+log(wind.speed)+I(log.cars^2)+I(log(wind.speed)^2))
#Has adjusted R^2=0.4899

#Keep removing
multireg8 = lm(no2~temp+log(wind.speed)+I(log.cars^2)+I(log(wind.speed)^2))
#Adj R^2=0.4804

#Continue
multireg9 = lm(no2~temp+log(wind.speed)+I(log.cars^2))
#Adj R^2=0.4791


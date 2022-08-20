light = c(28,-44,29,30,24,28,37,32,36,27,26,28
	   ,29, 26,27,22,23,20,25,25,36,23,31,32
         ,27, 33,16,24,29,36,21,28,26,27,27,25
         ,28, 24,40,21,31,32,28,26,30,27,26,24
         ,32, 29,34,-2,25,19,36,29,30,22,28,33
         ,39, 25,16,23)
hist(light)
boxplot(light)
plot(ecdf(light),verticals=T, do.points=F)
quantile(light)

#upper and lower of 95% confidence interval
mean(light) - qt(0.95,63)*(sd(age)/sqrt(64))    
# 20.50491
mean(light) + qt(0.95,63)*(sd(age)/sqrt(64))    
# 31.80759

light2 = c(28,29,30,24,28,37,32,36,27,26,28
	   ,29, 26,27,22,23,20,25,25,36,23,31,32
         ,27, 33,16,24,29,36,21,28,26,27,27,25
         ,28, 24,40,21,31,32,28,26,30,27,26,24
         ,32, 29,34,25,19,36,29,30,22,28,33
         ,39, 25,16,23)
mean(light2) - qt(0.95,61)*(sd(age)/sqrt(62))
# 21.99736
mean(light2) + qt(0.95,61)*(sd(age)/sqrt(62))
# 33.48651



#the right value is shown to be 33.02
#We see that the first confidence interval does not contain this value.
#That is fine, we did choose the 95% interval so there was a 5% probability
#of the real value being outside


#The version without the outliers does contain it so that's nice


#Start of ex 1b
rock.age=c(249,254,243,268,253,269,287,241,273,306,303,280,260,256,278,344,304,283,310)
bootsamp=sample(rock.age,replace=T)
sort(bootsamp); sort(rock.age)

bootagemean<-numeric(0)
for (i in 1:1000) bootagemean[i]<-mean(sample(rock.age,replace=T))
sort(bootagemean)[c(25,975)]
# 265.0526 289.9474

bootagemean[1:10]
#increase to 10000 samples
bootagemean<-numeric(0)
for (i in 1:10000) bootagemean[i]<-mean(sample(rock.age,replace=T))
sort(bootagemean)[c(25,975)]
# 261.4737 269.1053

#increase to 50000 samples
bootagemean<-numeric(0)
for (i in 1:50000) bootagemean[i]<-mean(sample(rock.age,replace=T))
sort(bootagemean)[c(25,975)]
# 259.0000 265.1053


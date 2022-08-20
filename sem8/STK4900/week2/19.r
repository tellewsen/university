#load data
wcgs=read.table("http://www.uio.no/studier/emner/matnat/math/STK4900/v11/wcgs.txt",sep="\t", header=T,na.strings=".")

# Our starting point in this exercise is the model for CHD risk using the predictors age (per 10 years), cholesterol (per 50 mg/dL), systolic blood pressure (per 50 mmHg), body mass index (per 10 kg/m2), smoking (yes, no), and behavioral pattern (with the four groups 1=A1, 2=A2, 3=B3, 4=B4).

#Fit the model
wcgs$behcat=factor(wcgs$behpat)
wcgs.beh=glm(chd69~age_10+chol_50+sbp_50+bmi_10+smoke+behcat, data=wcgs, family=binomial, subset=(chol<600))
summary(wcgs.beh)

#age,chol,sbp, smoking strongly correlated. BMI not as strong but correlated. Behcat 2 not strong. 3, partially like bmi, 4 weakly.

# In question a we consider a model with four behavioral groups. One may wonder if it is sufficient to consider only two behavioral groups (A-behavior and B-behavior). To this end we fit a model with the binary covariate dibpat (which is coded as 0 for B3 and B4, and as 1 for A1 and A2):

wcgs.beh2=glm(chd69~age_10+chol_50+sbp_50+bmi_10+smoke+dibpat, 
data=wcgs, family=binomial, subset=(chol<600))

# To test the null hypothesis that the effects behavioral patterns A1 and A2 are the same and  the effects behavioral patterns B1 and B2 are the same, we compare the (residual) deviances for the model in a and the model considered here:

anova(wcgs.beh2,wcgs.beh, test="Chisq")
#We get a small deviance between the two models, meaning that we can easily use this model instead of the old one.

 # We then fit a model without behavioral pattern and compare the three models in an analysis of deviance table:

wcgs.resc=glm(chd69~age_10+chol_50+bmi_10+sbp_50+smoke, data=wcgs, family=binomial, subset=(chol<600))

anova(wcgs.resc,wcgs.beh2,wcgs.beh, test="Chisq")
#lage deviance between model 1 and 2. Not as large between 2 and 3.

#Still model 2 that is the significant one. Meaning tha we do not need to distinguish between the two B behaviors and the two A behaviours. Only between A and B.
#I would keep model 2.

expcoef(wcgs.beh2)
#Output:
#               expcoef      lower      upper
#(Intercept) 0.03359428 0.02512283 0.04492233
#age_10      1.83025048 1.44753341 2.31415509
#chol_50     1.70240514 1.46583042 1.97716134
#sbp_50      2.46791666 1.64802505 3.69570393
#bmi_10      1.73236105 1.02991824 2.91389615
#smoke       1.82916180 1.38725756 2.41183250
#dibpat      2.00685487 1.51225096 2.66322626

#We note that inceasing sbp by 50 increases odds of CHD by 2.47. The same things can be said about other things as well, only with different degrees of raised odds. Note that dibpat and bmi are scary factors as well.

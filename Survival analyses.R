############################################
############################################
######       Survival Analyses        ######
############################################
############################################
############################################
# Outline
############################################
# Part 0: summary statistics
# Part 1: Fit CPH models to the different datasets
############################################
############################################
############################################
rm(list = ls())
library('survival')
setwd("~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/final_datasets_copied_from_phdfolder/survival")
data<-read.csv("~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/final_datasets_copied_from_phdfolder/survival/brucsurvival_TB3controls_longresidnomissing_noerrors_season2.csv")
data$age_yr2<-floor(data$age_yr)

data2<-read.csv("~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/final_datasets_copied_from_phdfolder/survival/brucsurvival_TB3controls_longresidnomissing_noerrors_season2_final_fixed.csv")
data2$age_yr2<-floor(data2$age_yr)

bolus <- c("R22", "R24", "R27", "R35", "R45", "R45b", "R50", "R6", "Y20", "Y31c", "Y31d", "Y39", "Y43", "Y44")
#"W1" is ok.
data3<- data2[!(data2$animal %in% bolus),]
length(unique(data2$animal)); length(unique(data3$animal))  # 142 to 128 = 14

data3<- data3[data3$animal!= "O10",]
data<- data3
data3$animal<-as.character(data3$animal)
############################################
############################################
# Part I: Fit CHP models to different datasets
# shows logical anlayses works great!
############################################
############################################

############################################
# initial plots & descriptive statistics
############################################
# plot time to death, doesn't tell you about time between infection and death, just overall time
d<-data[data$death.time==1,]
par(mfrow=c(2,2))
hist(d$stop2, ylab="Number of animals (all)") 
hist(d$stop2[d$TB_3==1 & d$brucella==0], ylab="Number of TB+ buffalo")
hist(d$stop2[d$TB_3==0 & d$brucella==1], ylab="Number of Bruc+ buffalo")
hist(d$stop2[d$TB_3==1 & d$brucella==1], ylab="Number of Coinfected buffalo")

data<- data3
# want a plot from time bruc+ to death in TB+ and TB- buffalo. 
mort<-data[data$death.time==1,]
mort2<-data2[data2$death.time==1,]
length(mort[,1])  # 38 total deaths
length(mort[mort$TB_3==1,1])  # 15 with bTB (6 bTB only)
length(mort[mort$brucella==1,1])  # 21 with brucellosis (12)
length(mort[mort$brucella==1& mort$TB_3==1,1]) # 9 of which are co-infected
length(mort[mort$brucella==0& mort$TB_3==0,1]) # 11 uninfected


hist(mort$start2, xlab="Months (from June 2008)", col="light gray", main="")
hist(mort2$start2, xlab="Months (from June 2008)", col="light gray", main="")
length(unique(mort$animal)); length(unique(mort2$animal)) # 42 vs. 49

mort<-data[data$death.time==1 & data$herd=="LS",]  # 18 total; 9 TB totoal (3 TB only), 10 BR only (4 Br only); 6 co.  

# Need a sum of time observing uninfected, TB+, Bruc+, Coinfected animals
id<-as.character(unique(data3$animal))
data3$animal<-as.character(data3$animal)

timedf<-data.frame(id=id, min_neg=NA, max_neg=NA, min_TB=NA, 
	max_TB=NA, min_B=NA, max_B=NA, min_C=NA, max_C=NA)
for (i in 1:length(id)){
	newdf<- data3[data3$animal==id[i],]
	# need if statments added in case of NA vlaues
	if(length(newdf$animal[newdf$TB_3==0 & newdf$brucella==0])>0){
	timedf$min_neg[i]<-min(newdf$start2[newdf$TB_3==0 & newdf$brucella==0])
	timedf$max_neg[i]<-max(newdf$stop2[newdf$TB_3==0 & newdf$brucella==0])
	} else {
		timedf$min_neg[i]<-0
		timedf$max_neg[i]<-0
	}
	
	if(length(newdf$animal[newdf$TB_3==1 & newdf$brucella==0]>0)){
		timedf$min_TB[i]<-min(newdf$start2[newdf$TB_3==1 & newdf$brucella==0])
		timedf$max_TB[i]<-max(newdf$stop2[newdf$TB_3==1 & newdf$brucella==0])
	} else {
		timedf$min_TB[i]<-0
		timedf$max_TB[i]<-0
	}

	if(length(newdf$animal[newdf$TB_3==0 & newdf$brucella==1]>0)){
		timedf$min_B[i]<-min(newdf$start2[newdf$TB_3==0 & newdf$brucella==1])
		timedf$max_B[i]<-max(newdf$stop2[newdf$TB_3==0 & newdf$brucella==1])
	} else {
		timedf$min_B[i]<-0
		timedf$max_B[i]<-0
	}

	if(length(newdf$animal[newdf$TB_3==1 & newdf$brucella==1]>0)){
	timedf$min_C[i]<-min(newdf$start2[newdf$TB_3==1 & newdf$brucella==1])
	timedf$max_C[i]<-max(newdf$stop2[newdf$TB_3==1 & newdf$brucella==1])
	} else {
		timedf$min_C[i]<-0
		timedf$max_C[i]<-0
	}
}
timeneg<-timedf$max_neg-timedf$min_neg; sum(timeneg) # 2031 all; 
timeTB<-timedf$max_TB-timedf$min_TB #600; 
timeB<-timedf$max_B-timedf$min_B; sum(timeB) # 1164; 
timeC<-timedf$max_C-timedf$min_C # 579; 
sum(timeB)+sum(timeneg)+sum(timeTB)+sum(timeC) # 4386 total months
6/600
12/1182
8/546
# mortality rates
10/(2109/12) #neg 
6/(666/12) # TB
13/(1086/12) # Bruc
9/(525/12) # co

length(mort$TB_3[mort$brucella==1& mort$TB_3==1])
length(mort$TB_3[mort$TB_3==1])
length(mort$TB_3[mort$brucella==1])

##############
# how many animals did we observe Br+ for > 2years
length(timeB[timeB/12 >2])
length(timeB[timeB/12 >2])
for (i in 1:length(timedf[,1])){
	timedf$maxBC[i] <- max(timedf$max_B[i], timedf$max_C[i])
	temp <- c(timedf$min_B[i], timedf$min_C[i])
	Bruc <- c(timedf$min_B[i], timedf$max_B[i])
	Co <- c(timedf$min_C[i], timedf$max_C[i])
	if(sum(Bruc)>0 & sum(Co)>0){
		timedf$minBC[i] <- min(temp)
	}
	if(sum(Bruc)>0 & sum(Co)==0){
		timedf$minBC[i] <- timedf$min_B[i]
	} 
	if(sum(Bruc)==0 & sum(Co)<0){
		timedf$minBC[i] <- timedf$min_C[i]
	} 
	if (sum(Bruc) == 0 & sum(Co)==0){
		timedf$minBC[i] <- 0
	}
	rm(temp, Bruc, Co)
	timedf$t[i] <- timedf$maxBC[i] - timedf$minBC[i]
}
length(timedf$t[timedf$t/12 >2])	
length(timedf$t[timedf$t/12 >1])	
	
###############
# plot of when animals died
par(mfrow=c(2,2))
hist(d$start2); hist(d$stop2) # spread throughout study period, maybe more in beginning/younger?
hist((d$start2 %% 12)/3); hist((d$stop2 %% 12)/3)

a<- hist((d$start2 %% 12)/3)


############################################
############################################
# New data, no bolus- age category
############################################
############################################
# = AGE1 in text
# 0-2.9, 3-5.9, 6+
data3$age2.2 <- data3$age2
data3$age2.2[data3$age_yr2 == 4] <- "subadult"
data3$age2.2[data3$age_yr2 == 5] <- "subadult"
data3$age2.2[data3$age_yr2 == 6] <- "subadult"
data3$age2.2[data3$age_yr2 == 7] <- "matureadult"
data3$age2.2 <- as.factor(as.character(data3$age2.2))
data3$age2.3 <- data3$age2.2
data3$age2.3[data3$age_yr2 == 3] <- "juvenile"
data3$age2.4 <- data3$age2.3
data3$age2.4[data3$age_yr2 == 6] <- "matureadult"
data3$age2.4[data3$age_yr2 == 3] <- "subadult"
data3$age2.4 <- as.character(data3$age2.4)
table(data3$age2.4, data3$age_yr)
test.mod<-coxph(Surv(start, stop, death.time)~brucella*age2.4+TB_3+herd2+ cluster(animal), data=data3); summary(test.mod)
test.mod<-coxph(Surv(start, stop, death.time)~age2.4+ cluster(animal), data=data3); summary(test.mod)
logLik(test.mod); extractAIC(test.mod)

# = AGE3 in text
# 0-1.9, 2-4.9, 5+ 
data3$age3.2 <- data3$age3
data3$age3.2[data3$age_yr2 == 4] <- "subadult"
data3$age3.2 <- as.factor(as.character(data3$age3.2))
table(data3$age3.2, data3$age_yr)
#test.mod<-coxph(Surv(start, stop, death.time)~brucella*age3.2+TB_3+herd2+ cluster(animal), data=data3); summary(test.mod)  # AIC = 316.9
test.mod<-coxph(Surv(start, stop, death.time)~age3.2 + cluster(animal), data=data3); summary(test.mod)
logLik(test.mod); extractAIC(test.mod)

# AGE 4 in TEXT!!
# 0-1.9, 2+
table(data3$age5, data3$age_yr)
test.mod<-coxph(Surv(start, stop, death.time)~brucella+TB_3+herd2+ age5+ cluster(animal), data=data3); summary(test.mod)
test.mod<-coxph(Surv(start, stop, death.time)~ age5+ cluster(animal), data=data3); summary(test.mod)
logLik(test.mod); extractAIC(test.mod)


# Additional Age Categories
# 1.4-1.9, 2-3.9, 4-8, 8+
# juveniles (0-1) and subadults (2-4) both have higher mortality than adults (4-8); subadults not significantly different than juveniles (p=0.11)
test.mod<-coxph(Surv(start, stop, death.time)~brucella*age1+TB_3+herd2+ age1+ cluster(animal), data=data3); summary(test.mod) # AIC = 317.686
test.mod<-coxph(Surv(start, stop, death.time)~age1+ cluster(animal), data=data3); test.mod$loglik # AIC = 336.93

# 1.4-2.9, 3-4.9, 5-7.9, 8+
# AIC = 314.7509; LogLike =  -167.0254 -151.3755
test.mod<-coxph(Surv(start, stop, death.time)~brucella*age2+TB_3+herd2+ cluster(animal), data=data3); summary(test.mod)  #LL error
test.mod<-coxph(Surv(start, stop, death.time)~age2+ cluster(animal), data=data3); summary(test.mod)  # AIC = 333.33

# 1-1.9, 2-3.9, 4+ 
test.mod<-coxph(Surv(start, stop, death.time)~brucella*age3+TB_3+herd2+ age3+ cluster(animal), data=data3); summary(test.mod)  # AIC = 314.51
test.mod<-coxph(Surv(start, stop, death.time)~age3+ cluster(animal), data=data3); summary(test.mod)# AIC = 335.5

# 0-2, 3-4, 5+ -> 0-3 yr olds have highest mortality. 3-5 vs 5+ similar
test.mod<-coxph(Surv(start, stop, death.time)~brucella+TB_3+herd2+ age4+ cluster(animal), data=data3); summary(test.mod)
# AIC = 313.242; LogLike = -167.0254 -151.6210




# AGE 2 IN TEXT
#0-2.9, 3+
data3$age6[data3$age_yr == 3] <- "adult"
test.mod<-coxph(Surv(start2, stop2, death.time)~ age6 + cluster(animal), data=data3); summary(test.mod)

final.mod<-coxph(Surv(start2, stop2, death.time)~brucella+TB_3+herd2+ age6+ cluster(animal), data=data3); summary(final.mod)
# AIC = 312.7035; LogLike = -167.0254 -152.3518
#coxph(formula = Surv(start, stop, death.time) ~ brucella + TB_3 + 
#    herd2 + age6 + cluster(animal), data = data3)
# n= 1462, number of events= 38 

#               coef exp(coef) se(coef) robust se     z Pr(>|z|)    
#brucella     1.1060    3.0224   0.3365    0.3505 3.155 0.001602 ** 
#TB_3         1.0370    2.8207   0.3609    0.3483 2.977 0.002907 ** 
#herd2CB      0.7351    2.0858   0.3426    0.3236 2.272 0.023087 *  
#age6juvenile 1.1825    3.2625   0.3408    0.3342 3.539 0.000402 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#             exp(coef) exp(-coef) lower .95 upper .95
#brucella         3.022     0.3309     1.521     6.008
#TB_3             2.821     0.3545     1.425     5.582
#herd2CB          2.086     0.4794     1.106     3.933
#age6juvenile     3.263     0.3065     1.695     6.280

#Concordance= 0.717  (se = 0.049 )
#Rsquare= 0.02   (max possible= 0.204 )
#Likelihood ratio test= 29.35  on 4 df,   p=6.645e-06
#Wald test            = 34.51  on 4 df,   p=5.86e-07
#Score (logrank) test = 32.82  on 4 df,   p=1.298e-06,   Robust = 15.46  p=0.00383

vcov(final.mod)
# var for bruc + bTB = (1)*V(bruc) + (1)*V(bTB) + 2*(1)*(1)*Cov(bruc,bTB)
# se = squrt(var) = 0.5004329
(0.121814 + 0.1186884 + 0.00993067)^0.5
# CI: 
exp(1.106 + 1.0370)
exp(1.106+1.0370 + 1.96*0.500)
exp(1.106+1.0370 - 1.96*0.500)

###################################################
data3$infection <- NA
data3$infection[data3$TB_3 == 0 & data3$brucella == 0] <- "AS"
data3$infection[data3$TB_3 == 1 & data3$brucella == 0] <- "TB"
data3$infection[data3$TB_3 == 0 & data3$brucella == 1] <- "BR"
data3$infection[data3$TB_3 == 1 & data3$brucella == 1] <- "C"
test.mod <- coxph(Surv(start2, stop2, death.time)~infection+herd2+ age6+ cluster(animal), data=data3)


######!!!!!!!!  # SEE AGE-SPECIFIC PATTERN WITH BRUCELLOSIS!
data3$age2.4 <- as.factor(data3$age2.4)
data3$age2.4 <- relevel(data3$age2.4, "matureadult")
data3$age2.4 <- relevel(data3$age2.4, "subadult")
final.mod<-coxph(Surv(start2, stop2, death.time)~brucella+TB_3+herd2+ age2.4*brucella+ cluster(animal), data=data3); summary(final.mod)
predRes <- predict(final.mod, type="risk") # residuals
head(predRes, n=10)
Shat2 <- survexp(~ age2.4 + TB_3, ratetable=final.mod, data=data3[data3$brucella == 0,])
with(Shat2, head(data.frame(time, surv), n=4))
Shat2 <- survexp(~ age2.4 + TB_3, ratetable=final.mod, data=data3[data3$brucella == 1,])
with(Shat2, head(data.frame(time, surv), n=4))


# Juvenile: 
#                               coef exp(coef) se(coef) robust se      z Pr(>|z|)   
#brucella                    1.58217   4.86552  0.50558   0.51212  3.089  0.00201 **
#TB_3                        1.13865   3.12255  0.36471   0.34905  3.262  0.00111 **
#herd2CB                     0.75853   2.13514  0.34770   0.32821  2.311  0.02083 * 
#age2.4matureadult          -1.19107   0.30390  1.10847   1.16339 -1.024  0.30594   
#age2.4subadult             -0.61855   0.53872  0.55199   0.57149 -1.082  0.27910   
#brucella:age2.4matureadult  0.08662   1.09048  1.22112   1.26645  0.068  0.94547   
#brucella:age2.4subadult    -1.57650   0.20670  0.78290   0.77776 -2.027  0.04267 * 


# Subadult: 
#                                coef exp(coef)  se(coef) robust se      z Pr(>|z|)   
#brucella                    0.005675  1.005691  0.609981  0.608143  0.009  0.99255   
#TB_3                        1.138651  3.122554  0.364711  0.349049  3.262  0.00111 **
#herd2CB                     0.758531  2.135138  0.347704  0.328212  2.311  0.02083 * 
#age2.4matureadult          -0.572513  0.564106  1.057936  1.096687 -0.522  0.60164   
#age2.4juvenile              0.618554  1.856242  0.551986  0.571493  1.082  0.27910   
#brucella:age2.4matureadult  1.663119  5.275739  1.263485  1.294638  1.285  0.19892   
#brucella:age2.4juvenile     1.576498  4.837983  0.782904  0.777763  2.027  0.04267 * 




data3$age7 <- data3$age3							 
data3$age7[data3$age7 == "juvenile"] <- "subadult"   # pre-repro (0-3)
data3$age7 <- as.character(data3$age7)				# post-repro (4+)

data3$age8 <- data3$age7							 
data3$age8[data3$age_yr2 == 4] <- "subadult"         # pre-repro (0-4)
data3$age8 <- as.character(data3$age8)				 # post-repro (5+)

data3$age9 <- data3$age2							 
data3$age9[data3$age_yr2 == 3] <- "juvenile"         # pre-repro (0-3)
data3$age9[data3$age_yr2 == 5] <- "subadult" # pre-repro (4-6)
data3$age9[data3$age_yr2 == 6] <- "subadult" # pre-repro (4-6)
data3$age9[data3$age_yr2 == 7] <- "matureadult"        #(7+)
data3$age9 <- as.factor(as.character(data3$age9))				 # post-repro (6+)
data3$age9 <- relevel(data3$age9, "matureadult")
########################################################
# CPH model diagnostics
# no change with start/stop designations: exact same covariates
# model with tb*herd prefered by AIC; but parameter value= 0.1
test.mod<- final.mod

# test proportional hazards assumption
temp<-cox.zph(test.mod)
print(temp);
par(mfrow=c(2,3))
 plot(temp)
#P. Grambsch and T. Therneau (1994), Proportional hazards tests and diagnostics based on weighted residuals. Biometrika, 81, 515-26.

# residuals
res<-resid(test.mod)

# get predicted
predRes <- predict(test.mod, type="risk")
head(predRes, n=10)
Shat2 <- survexp(~ TB_3+ brucella, ratetable=test.mod, data=data3)
with(Shat2, head(data.frame(time, surv), n=4))

test.mod<- final.mod
predRes <- predict(test.mod, type="risk")
head(predRes, n=10)
Shat2 <- survexp(~ age6 + herd, ratetable=test.mod, data=data3)
with(Shat2, head(data.frame(time, surv), n=4))

test.mod<- final.mod
predRes <- predict(test.mod, type="risk")
head(predRes, n=10)
Shat2 <- survexp(~ age6, ratetable=test.mod, data=data3)
with(Shat2, head(data.frame(time, surv), n=8))

plot(survfit(Surv(start, stop, death.time)~herd2+ age6+ cluster(animal), data=data3), conf.int=FALSE)

plot(survfit(Surv(start, stop, death.time)~herd2+ age6+ cluster(animal), data=data3), mark.time=FALSE)
lines(survexp(~ age6 + herd, ratetable=test.mod, data=data3), col='purple')

# predicted proportional increase at average buffalo (herd in between the LS & CB)
#data$herd2<-NA
#for(i in 1:length(data$herd)){
#ifelse(data$herd[i]=="LS", data$herd2[i]<-1, data$herd2[i]<-0)
#}
#summary(data$herd2) #0.5438
#data$herd3<-data$herd2-0.5438
#full.mod<-coxph(Surv(start, stop, death.time)~brucella*herd3+ TB_3*herd3+age_yr2+ I(age_yr2^2), data=data)


# test proportional hazards assumption
temp<-cox.zph(final.mod)
print(temp);
par(mfrow=c(2,3))
 plot(temp)
#P. Grambsch and T. Therneau (1994), Proportional hazards tests and diagnostics based on weighted residuals. Biometrika, 81, 515-26.

# residuals
res<-resid(test.mod)

# get predicted
predRes <- predict(final.mod, type="risk")
head(predRes, n=10)
Shat2 <- survexp(~ age6 + TB_3+ brucella, ratetable=final.mod, data=data3)
with(Shat2, head(data.frame(time, surv), n=4))


# CUT - Delete later

############################
# New data, no bolus- age continuous
############################
# just brucellosis
test.mod<-coxph(Surv(start, stop, death.time)~brucella+ cluster(animal), data=data3) # 0.0119
test.mod<-coxph(Surv(start, stop, death.time)~brucella+herd2+ cluster(animal), data=data3); summary(test.mod) #0.0134 
test.mod<-coxph(Surv(start, stop, death.time)~brucella+herd2+age_yr2+ cluster(animal), data=data3); summary(test.mod) #0.00457
test.mod<-coxph(Surv(start, stop, death.time)~brucella+herd2+age_yr2+ I(age_yr2^2)+ cluster(animal), data=data3); summary(test.mod) # 0.00533
test.mod<-coxph(Surv(start, stop, death.time)~brucella*herd2+age_yr2+ I(age_yr2^2)+ cluster(animal), data=data3); summary(test.mod) # 0.01, int n.s. 

# just tb
test.mod<-coxph(Surv(start, stop, death.time)~ TB_3+ cluster(animal), data=data3) ; summary(test.mod) #0.0172 
test.mod<-coxph(Surv(start, stop, death.time)~ TB_3 +herd2+ cluster(animal), data=data3); summary(test.mod)  # 0.0125
test.mod<-coxph(Surv(start, stop, death.time)~ TB_3 +herd2+age_yr2+ cluster(animal), data=data3) ; summary(test.mod) #0.00934
test.mod<-coxph(Surv(start, stop, death.time)~ TB_3 +herd2+age_yr2+ I(age_yr2^2)+ cluster(animal), data=data3); # 0.00216 summary(test.mod)
test.mod<-coxph(Surv(start, stop, death.time)~ TB_3*herd2+age_yr2+ I(age_yr2^2)+ cluster(animal), data=data3); summary(test.mod)
#0.00110  

# both
test.mod<-coxph(Surv(start, stop, death.time)~brucella+ TB_3+ cluster(animal), data=data3); summary(test.mod) # Br- 0.01,TB-0.02
test.mod<-coxph(Surv(start, stop, death.time)~brucella*TB_3+ cluster(animal), data=data3); summary(test.mod) #ns
test.mod<-coxph(Surv(start, stop, death.time)~brucella+ TB_3+herd2+ cluster(animal), data=data3); summary(test.mod)#.01,0.01
#*******test.mod<-coxph(Surv(start, stop, death.time)~brucella+TB_3+herd2 + age_yr2 +I(age_yr2^2)+ cluster(animal), data=data3); summary(test.mod)


test.mod<-coxph(Surv(start, stop, death.time)~brucella*herd2+TB_3*herd2 + age_yr2 +I(age_yr2^2)+ cluster(animal), data=data3); summary(test.mod)

final.mod<-coxph(Surv(start, stop, death.time)~brucella+TB_3+herd2 + age_yr2 +I(age_yr2^2)+ cluster(animal), data=data3); summary(final.mod)

standardize = function(datavec){
    temp<- NA
    for (i in 1:length(datavec)){
        temp[i]<- (datavec[i]-mean(datavec)) / (2*sd(datavec))
    }
    return(temp)
}
data3$age_yrsd<- NA; data3$age_yrsq<- NA; data3$herdsd<-NA; data3$tbsd<- NA; data3$brucsd<-NA
data3$age_yrsd<-standardize(data3$age_yr)
t<- data3$age_yr*data3$age_yr
data3$age_yrsq<-standardize(t)
data3$therd<-NA
data3$therd[data3$herd=="CB"]<- 1
data3$therd[data3$herd=="LS"]<- 0
data3$herdsd<- standardize(data3$therd)
data3$tbsd<- standardize(data3$TB_3)
data3$brucsd<- standardize(data3$brucella)

# standardized output
test.mod<-coxph(Surv(start2, stop2, death.time)~brucsd+tbsd+ herdsd + age_yrsd +I(age_yrsd^2)+ cluster(animal), data=data3)  
#                coef exp(coef) se(coef) robust se      z Pr(>|z|)   
#brucsd         0.8826    2.4172   0.3365    0.3427  2.576  0.01001 * 
#tbsd           0.8983    2.4555   0.3251    0.3058  2.937  0.00331 **
#herdsd         0.7981    2.2214   0.3628    0.3466  2.303  0.02129 * 
#age_yrsd      -1.5095    0.2210   0.5051    0.5078 -2.973  0.00295 **
#I(age_yrsd^2)  0.8707    2.3887   0.3199    0.2981  2.921  0.00349 **
#              exp(coef) exp(-coef) lower .95 upper .95
#brucsd            2.417     0.4137    1.2348    4.7316
#tbsd              2.456     0.4072    1.3484    4.4716
#herdsd            2.221     0.4502    1.1261    4.3818
#age_yrsd          0.221     4.5243    0.0817    0.5979
#(age_yrsd^2)     2.389     0.4186    1.3316    4.2849

library('JMbayes')
library('JM')
#http://www.r-bloggers.com/joint-models-for-longitudinal-and-survival-data/
#http://www.r-bloggers.com/dynamic-predictions-using-joint-models/

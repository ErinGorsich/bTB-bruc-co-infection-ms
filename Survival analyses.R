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
# Part II: Fit best (no bolus, add data) with Bayesain model
# Part III: Joint models
############################################
############################################
############################################
# Started: 1 April 2016
library('JMbayes')
library('JM')
#http://www.r-bloggers.com/joint-models-for-longitudinal-and-survival-data/
#http://www.r-bloggers.com/dynamic-predictions-using-joint-models/
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

###############
# plot of when animals died
par(mfrow=c(2,2))
hist(d$start2); hist(d$stop2) # spread throughout study period, maybe more in beginning/younger?
hist((d$start2 %% 12)/3); hist((d$stop2 %% 12)/3)

a<- hist((d$start2 %% 12)/3)




############################################
# run cph model
############################################
############################
# Part 1: Old data, from thesis
############################
# just brucellosis
test.mod<-coxph(Surv(start, stop, death.time)~brucella, data=data)  # 0.11
test.mod<-coxph(Surv(start, stop, death.time)~brucella+herd2, data=data) # 0.12
test.mod<-coxph(Surv(start, stop, death.time)~brucella+herd2+age_yr2, data=data) # 0.07
test.mod<-coxph(Surv(start, stop, death.time)~brucella+herd2+age_yr2+ I(age_yr2^2), data=data) #ns.
test.mod<-coxph(Surv(start, stop, death.time)~brucella*herd2+age_yr2+ I(age_yr2^2), data=data) #in LS, notCB

# just tb
test.mod<-coxph(Surv(start, stop, death.time)~ TB_3, data=data)  # 0.11
test.mod<-coxph(Surv(start, stop, death.time)~ TB_3 +herd2, data=data)
test.mod<-coxph(Surv(start, stop, death.time)~ TB_3 +herd2+age_yr2, data=data)
test.mod<-coxph(Surv(start, stop, death.time)~ TB_3 +herd2+age_yr2+ I(age_yr2^2), data=data)
test.mod<-coxph(Surv(start, stop, death.time)~ TB_3*herd2+age_yr2+ I(age_yr2^2), data=data)

# both
test.mod<-coxph(Surv(start, stop, death.time)~brucella+ TB_3, data=data)  # TB=0.06, bruc=0.11
test.mod<-coxph(Surv(start, stop, death.time)~brucella*TB_3, data=data)  #  all n.s.
test.mod<-coxph(Surv(start, stop, death.time)~brucella+ TB_3+herd2, data=data) #TB=0.05,Bruc=0.13,herd=0.04
test.mod<-coxph(Surv(start, stop, death.time)~brucella*TB_3+herd2, data=data) # all n.s. 

test.mod<-coxph(Surv(start, stop, death.time)~brucella+ TB_3+herd2+age_yr2, data=data) # TB=0.04, Br=0.08
test.mod<-coxph(Surv(start, stop, death.time)~brucella*TB_3+herd2+age_yr2+I(age_yr2^2),data=data)#TB=0.09,ns
test.mod<-coxph(Surv(start, stop, death.time)~brucella*herd2+ TB_3+age_yr2+ I(age_yr2^2), data=data)#all sig
test.mod<-coxph(Surv(start, stop, death.time)~brucella*herd2+ TB_3*herd2+age_yr2+ I(age_yr2^2), data=data)
test.mod<-coxph(Surv(start, stop, death.time)~brucella*herd2*TB_3+age_yr2+ I(age_yr2^2), data=data) #nope

#Selection
test.mod<-coxph(Surv(start, stop, death.time)~brucella*herd2+ TB_3*herd2+brucella*age_yr2+ TB_3*age_yr2+ brucella*I(age_yr2^2)+ TB_3*I(age_yr2^2), data=data)
test.mod<-coxph(Surv(start, stop, death.time)~brucella*herd2+ TB_3*herd2+brucella*age_yr2+ TB_3*age_yr2+ I(age_yr2^2), data=data)
# top two
full.mod<-coxph(Surv(start, stop, death.time)~brucella*herd2+ TB_3*herd2+age_yr2+ I(age_yr2^2), data=data)
red.mod<-coxph(Surv(start, stop, death.time)~brucella+ TB_3*herd2+age_yr2+ I(age_yr2^2), data=data)
# choose model with lowest AIC= full.mod
# choose based on drop in deviance= red.mod (chi.sq= 2.55, df=1, p=0.11)


############################
# Part 2: New data
############################
# just brucellosis
test.mod<-coxph(Surv(start, stop, death.time)~brucella, data=data2)  # 0.0264
test.mod<-coxph(Surv(start, stop, death.time)~brucella+herd2, data=data2); summary(test.mod) #0.0352  
test.mod<-coxph(Surv(start, stop, death.time)~brucella+herd2+age_yr2, data=data2); summary(test.mod) #0.02
test.mod<-coxph(Surv(start, stop, death.time)~brucella+herd2+age_yr2+ I(age_yr2^2), data=data2); summary(test.mod) # 0.03417 *
test.mod<-coxph(Surv(start, stop, death.time)~brucella*herd2+age_yr2+ I(age_yr2^2), data=data2); summary(test.mod) # 0.01769 

# just tb
test.mod<-coxph(Surv(start, stop, death.time)~ TB_3, data=data2) ; summary(test.mod) #0.0776 
test.mod<-coxph(Surv(start, stop, death.time)~ TB_3 +herd2, data=data2); summary(test.mod) #0.04835 * 
test.mod<-coxph(Surv(start, stop, death.time)~ TB_3 +herd2+age_yr2, data=data2) 
test.mod<-coxph(Surv(start, stop, death.time)~ TB_3 +herd2+age_yr2+ I(age_yr2^2), data=data2); # 0.01074 summary(test.mod)
test.mod<-coxph(Surv(start, stop, death.time)~ TB_3*herd2+age_yr2+ I(age_yr2^2), data=data2); summary(test.mod)
#0.001185 

# both
test.mod<-coxph(Surv(start, stop, death.time)~brucella+ TB_3, data=data2); summary(test.mod) # Br- 0.02,TB-0.08
test.mod<-coxph(Surv(start, stop, death.time)~brucella*TB_3, data=data2); summary(test.mod) #ns
test.mod<-coxph(Surv(start, stop, death.time)~brucella+ TB_3+herd2, data=data2); summary(test.mod)#.04,.06,.01
test.mod<-coxph(Surv(start, stop, death.time)~brucella*TB_3+herd2, data=data2); summary(test.mod) #ns

test.mod<-coxph(Surv(start, stop, death.time)~brucella+ TB_3+herd2+age_yr2, data=data2); summary(test.mod)
******test.mod<-coxph(Surv(start, stop, death.time)~brucella+TB_3+herd2+age_yr2+I(age_yr2^2),data=data2); summary(test.mod) # 0.04, 0.01!
test.mod<-coxph(Surv(start, stop, death.time)~brucella*herd2+ TB_3+age_yr2+ I(age_yr2^2), data=data2); summary(test.mod) # 0.02, 0.01, 0.13
test.mod<-coxph(Surv(start, stop, death.time)~brucella*herd2+ TB_3*herd2+age_yr2+ I(age_yr2^2), data=data2); summary(test.mod) # 0.02, 0.001, int suggestive
*****test.mod<-coxph(Surv(start, stop, death.time)~brucella+ TB_3*herd2+age_yr2+ I(age_yr2^2), data=data2); summary(test.mod) # 0.05, 0.001


############################
# Part 3: New data, no bolus
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
*******test.mod<-coxph(Surv(start, stop, death.time)~brucella+TB_3+herd2 + age_yr2 +I(age_yr2^2)+ cluster(animal), data=data3); summary(test.mod)


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
                coef exp(coef) se(coef) robust se      z Pr(>|z|)   
brucsd         0.8826    2.4172   0.3365    0.3427  2.576  0.01001 * 
tbsd           0.8983    2.4555   0.3251    0.3058  2.937  0.00331 **
herdsd         0.7981    2.2214   0.3628    0.3466  2.303  0.02129 * 
age_yrsd      -1.5095    0.2210   0.5051    0.5078 -2.973  0.00295 **
I(age_yrsd^2)  0.8707    2.3887   0.3199    0.2981  2.921  0.00349 **
              exp(coef) exp(-coef) lower .95 upper .95
brucsd            2.417     0.4137    1.2348    4.7316
tbsd              2.456     0.4072    1.3484    4.4716
herdsd            2.221     0.4502    1.1261    4.3818
age_yrsd          0.221     4.5243    0.0817    0.5979
I(age_yrsd^2)     2.389     0.4186    1.3316    4.2849

# CPH model diagnostics
# no change with start/stop designations: exact same covariates
# model with tb*herd prefered by AIC; but parameter value= 0.1
test.mod<-coxph(Surv(start2, stop2, death.time)~brucella+TB_3+herd2 + age_yr2 +I(age_yr2^2)+ cluster(animal), data=data3)  

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
Shat2 <- survexp(~ TB_3+ brucella, ratetable=test.mod, data=data)
with(Shat2, head(data.frame(time, surv), n=4))
#time TB_3.0..brucella.0 TB_3.0..brucella.1 TB_3.1..brucella.0 TB_3.1..brucella.1
#1    3          0.9156157          0.8279407          0.8255238          0.6367819
#2    6          0.9006887          0.8003969          0.7969077          0.5878509
#3   12          0.8948555          0.7898502          0.7858978          0.5697411
#4   15          0.8628729          0.7341025          0.7272129          0.4796214

# predicted proportional increase at average buffalo (herd in between the LS & CB)
#data$herd2<-NA
#for(i in 1:length(data$herd)){
#ifelse(data$herd[i]=="LS", data$herd2[i]<-1, data$herd2[i]<-0)
#}
#summary(data$herd2) #0.5438
#data$herd3<-data$herd2-0.5438
#full.mod<-coxph(Surv(start, stop, death.time)~brucella*herd3+ TB_3*herd3+age_yr2+ I(age_yr2^2), data=data)
full.mod<-test.mod

# plot predicted survival time (same as before but good): 
df<-data.frame(start=data3$start2, stop=data3$stop2, death.time=data3$death.time, TB=data3$TB_3, Br=data3$brucella, age=data3$age_yr2, herd=data3$herd)
# set at average herd value (0.53) so not LS or CB. 
mort<-with(df, data.frame(TB=0, Br=0, age=rep(mean(age),2), herd=rep(mean(herd=="LS"), 2)))
mortTB<-with(df, data.frame(TB=1, Br=0, age=rep(mean(age),2), herd=rep(mean(herd=="LS"), 2)))
mortBr<-with(df, data.frame(Br=1, TB=0, age=rep(mean(age),2), herd=rep(mean(herd=="LS"), 2)))
mortco<-with(df, data.frame(Br=1, TB=1, age=rep(mean(age),2), herd=rep(mean(herd=="LS"), 2)))
plot_add.mod<-coxph(Surv(start, stop, death.time)~Br+herd+TB+herd+age+I(age^2), data=df)
m<-survfit(plot_add.mod, newdata=mort)
mt<-survfit(plot_add.mod, newdata=mortTB)
mb<-survfit(plot_add.mod, newdata=mortBr)
mco<-survfit(plot_add.mod, newdata=mortco)

plot(m, conf.int=FALSE, ylab="Survival", xlab="Time (months)", lty=c(1, 2),
	ylim=c(0.1, 1), cex.lab=1.4, col="dark blue")
lines(mt, lty=c(5, 5), conf.int= FALSE, col="dark green")
lines(mb, lty=c(4, 4), conf.int= FALSE, col="purple")
lines(mco, lty=c(3,3), conf.int= FALSE, col="dark red")
legend("bottomleft", legend=c("Uninfected", "bTB+", "Brucellosis +", "Coinfected"),
	lty=c(1 ,5, 4, 3), inset=0.02, bty="n", col=c("dark blue", "dark green", "purple", "dark red"))



############################################
############################################
# Repeat above with age category
############################################
############################################
# 0-2, 2-4, 4-8, 8+
# juveniles (0-2) and subadults (2-4) both have higher mortality than adults (4-8); subadults not significantly different than juveniles (p=0.11)
test.mod<-coxph(Surv(start, stop, death.time)~brucella+TB_3+herd2+ age1+ cluster(animal), data=data3); summary(test.mod) 
# AIC = 320.9481; LogLike =  -167.0254 -154.4741
test.mod<-coxph(Surv(start, stop, death.time)~age1+ cluster(animal), data=data3); test.mod$loglik 


# 0-3, 3-5, 5-8, 8+
test.mod<-coxph(Surv(start, stop, death.time)~brucella+TB_3+herd2+ age2+ cluster(animal), data=data3); summary(test.mod) 
# AIC = 314.7509; LogLike =  -167.0254 -151.3755

# 0-2, 2-4, 4+ -> 0-2 yr olds have highest mortality. 2-4 vs 4+ similar 
test.mod<-coxph(Surv(start, stop, death.time)~brucella+TB_3+herd2+ age3+ cluster(animal), data=data3); summary(test.mod)
# AIC = 319.3569; LogLike = -167.0254  -154.6784

# 0-3, 3-5, 5+ -> 0-3 yr olds have highest mortality. 3-5 vs 5+ similar
test.mod<-coxph(Surv(start, stop, death.time)~brucella+TB_3+herd2+ age4+ cluster(animal), data=data3); summary(test.mod)
# AIC = 313.242; LogLike = -167.0254 -151.6210

# 0-2, 2+
test.mod<-coxph(Surv(start, stop, death.time)~brucella+TB_3+herd2+ age5+ cluster(animal), data=data3); summary(test.mod)
test.mod<-coxph(Surv(start, stop, death.time)~ age5+ cluster(animal), data=data3); summary(test.mod)
# AIC = 323.0442; LogLike = -167.0254 -157.5221

#0-3, 3+
final.mod<-coxph(Surv(start, stop, death.time)~brucella+TB_3+herd2+ age6+ cluster(animal), data=data3); summary(final.mod)
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








# age as a continuous predictor: 
test.mod<-coxph(Surv(start2, stop2, death.time)~brucella+TB_3+ herd2 + age_yrsd +I(age_yrsd^2)+ cluster(animal), data=data3)  
# AIC = 317.2427, LogLike = -167.0254 -153.6213


# still no herd effects
test.mod<-coxph(Surv(start, stop, death.time)~brucella+TB_3*herd2+ age4+ cluster(animal), data=data3); summary(test.mod)
test.mod<-coxph(Surv(start, stop, death.time)~brucella*herd2+TB_3*herd2+ age4+ cluster(animal), data=data3); summary(test.mod)
test.mod<-coxph(Surv(start, stop, death.time)~brucella*herd2+TB_3+herd2+ age4+ cluster(animal), data=data3); summary(test.mod)

# nor age... but worry about age*bruc colinearity
test.mod<-coxph(Surv(start, stop, death.time)~brucella*age4+TB_3+herd2+ age4+ cluster(animal), data=data3); summary(test.mod)
test.mod<-coxph(Surv(start, stop, death.time)~brucella+age4*TB_3+herd2+ age4+ cluster(animal), data=data3); summary(test.mod)

# nor co-infection
test.mod<-coxph(Surv(start, stop, death.time)~brucella*TB_3+herd2+ age4+ cluster(animal), data=data3); summary(test.mod)




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





############################################
############################################
# Part II: Fit best (no bolus, add data) with Bayesain model.... no luck!
############################################
############################################
#!!!!!!!!!!!!!!
#library("spBayesSurv") nope
# survregbayes nope
#https://cran.r-project.org/web/packages/spBayesSurv/spBayesSurv.pdf


fit = coxph(Surv(start2, stop2, death.time)~brucella+TB_3+herd2 + age_yr2 +I(age_yr2^2), data=data3)

# can't fit time dependent covariates with survreg... 
# fit0 = survreg(formula = Surv(start2, stop2, death.time)~brucella+
#	TB_3+herd2 + age_yr2 +I(age_yr2^2), data=data3, dist = "exponential")

test = survregbayes(formula = Surv(start2, stop2, death.time)~brucella+TB_3+
	herd2 + age_yr2 +I(age_yr2^2), data= data3, survmodel="PH")

coxph(Surv(start2, stop2, death.time)~brucella+TB_3+herd2 + age_yr2 +I(age_yr2^2),, data=data3)



############################################
############################################
# Part III: Joint models
############################################
############################################
final.mod<-coxph(Surv(start, stop, death.time)~brucella+TB_3+herd2 + age_yr2 +I(age_yr2^2), data=data3); summary(final.mod)

library("nlme")
library("lme4")
immune<-read.csv("~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/final_datasets_copied_from_phdfolder/cross_sectional_data_withdz_cleandisease_nofinal_Feb2016_capturetime_forsurv.csv")  # 5 animals added from last time
# added these five to cross_sectional_data_withdz_cleandisease_nofinal_Feb2016_capturetime.csv because in survival analyses.
immune2<-immune[!is.na(immune$ifng),]
immune3<-immune2[!is.na(immune2$ifng_plate),]  # 1 NA plate value!
immune3<-immune3[!is.na(immune3$ifng_delay),]
immune3<-immune3[immune3$ifng_CV<15,]
immune3$ifng_plate<-as.factor(immune3$ifng_plate)
immune3$herd<-relevel(immune3$herdorig, "LS")
immuneg<-immune3

# datasets need same numbers of individuals
im_id<-unique(immuneg$id); mort_id<-unique(data3$animal)
exclude_from_surv<- mort_id[!(mort_id %in% im_id)] # All buffs in the mort dataset are in immune data!
exclude_from_immune<- im_id[!(im_id %in% mort_id)] # 24 removed from gamma dataset!
immuneg2<-immuneg[!(immuneg$id %in% exclude_from_immune),]
immuneg2$id<-as.factor(as.character(immuneg2$id))
data3$animal<-as.factor(as.character(data3$animal))
unique(immuneg2$id)
unique(data3$animal)


# is there an association with current gamma levels and mortality: 
############################################
sub<-as.character(unique(data3$animal))
df <- data.frame(id= sub, eventtime= NA, maxstart= NA, maxstop= NA, maxtime= NA)
for (i in 1:length(sub)){
	temp = data3[data3$animal==sub[i],]
	if (temp$death[1] > 0){
	df$eventtime[i]<-temp$start2[temp$death.time==1]}
	df$maxstart[i]<- max(temp$start2)
	df$maxstop[i]<- max(temp$stop2)
	temp = immuneg[immuneg$id==sub[i],]
	df$maxtime[i] <- max(temp$capturetime)
	######
}
df$id[df$maxtime > df$maxstart+4]  # only O10, typeo in gamma dat... no an animal censored early when it went missing

immuneg2$problem<-NA
#for (i in 1:length(df[,1])){
#	buff<- as.character(df$id[i])
#	timeremove<- df$maxtime[i]
#	immuneg2$problem[immuneg2$id==buff & immuneg2$capturetime==timeremove]<-1
#}

#immuneg2$problem[immuneg2$id=="O10" & immuneg2$capturetime > 10]<- 1
#immuneg2$problem[is.na(immuneg2$problem)]<-0
#immuneg3<- immuneg2[immuneg2$problem==0,]

#immuneg3<- immuneg2[immuneg2$id != "O10",]

#im_id<-unique(immuneg3$id); mort_id<-unique(data3$animal)
#exclude_from_surv<- mort_id[!(mort_id %in% im_id)] # All buffs in the mort dataset are in immune data!
#exclude_from_immune<- im_id[!(im_id %in% mort_id)]


#sub<- c("O1b", "O26b", "O30", "O31", "O42", "O47b", "R20", "R28", "R34b", "R38", "R39b", "R43", "R46b", "R7b", "Y10b", "Y16b", "Y33b", "Y36", "Y38b", "Y4", "Y40", "Y42", "Y46b", "Y47")
#immune4<- immune3[!(immune3$id %in% sub),]
#data4<- data3[!(data3$animal %in% sub),]

data3$start3<-data3$start2+3
data3$stop3<-data3$stop2+3
mort.fit = coxph(Surv(start3, stop3, death.time)~brucella+TB_3+herd2 + age_yr2 
	+I(age_yr2^2) + cluster(animal), data=data3, model=TRUE, x=TRUE)
#gamma.fit0 = lmer(log(ifng)~tb+ bruc+ herd+ bruc:herd+ (capturetime|id) + (1| ifng_plate), data=immuneg2)
#immuneg2$gamma<- log(immuneg2$ifng)
#gamma.fit = lme(log(ifng) ~ tb+ bruc+ herd+ bruc:herd+ capturetime, random=~capturetime|id, data=immuneg2)
immuneg2$time<- immuneg2$capturetime+3.1
gamma.fit = lme(log(ifng) ~ tb + bruc+ herd+ bruc:herd+ time, random= ~1|id, data=immuneg2)

test<-jointModelBayes(gamma.fit, mort.fit, timeVar="time", method="spline-PH-aGH")






# remove those with 0 event times: B42, O29, O32, O47, R46, 
rm<-c("B42", "O29", "O32", "O47", "R46")
im<-immuneg2[!(immuneg2$id %in% rm), ]
mort<- data3[!(data3$animal %in% rm),]

gamma.fit <- lme(log(ifng) ~ tb+ bruc+ herd+ bruc:herd+ time, random=~1|id, data=im)
mort.fit <- coxph(Surv(start3, stop3, death.time)~brucella+TB_3+herd2 + age_yr2 
                  +I(age_yr2^2) + cluster(animal), data=mort, model=TRUE, x=TRUE)
test<-jointModelBayes(gamma.fit, mort.fit, timeVar="time", method="spline-PH-aGH")



###########################################
# Hapto
###########################################
immune<-read.csv("~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/final_datasets_copied_from_phdfolder/cross_sectional_data_withdz_cleandisease_nofinal_Feb2016_capturetime_forsurv.csv")  # 5 animals added from last time
# added these five to cross_sectional_data_withdz_cleandisease_nofinal_Feb2016_capturetime.csv because in survival analyses.
immune2<-immune[!is.na(immune$hapto),]
immune3<-immune2[!is.na(immune2$hapto_plate),]  # 1 NA plate value!
immune3$hapto_plate<-as.factor(immune3$hapto_plate)
immuneg<-immune3

# datasets need same numbers of individuals
im_id<-unique(immuneg$id); mort_id<-unique(data3$animal)
exclude_from_surv<- mort_id[!(mort_id %in% im_id)] # 12 buffs in the mort dataset are not in immune data
exclude_from_immune<- im_id[!(im_id %in% mort_id)] # 16 removed from gamma dataset!
immuneg2<-immuneg[!(immuneg$id %in% exclude_from_immune),]
immuneg<-immuneg2
data<- data3[!(data3$animal %in% exclude_from_surv),]
length(unique(immuneg$id)); length(unique(data$animal)) # should be 115 each
table(immuneg$id); table(data$animal)
immuneg$id<-as.factor(as.character(immuneg$id))
data$animal<-as.factor(as.character(data$animal))
colnames(data)[3]<-"start2"


# remove animals that died in the first capture period, none. 
# remove animals that have only one measurement. 
table(data$animal)
rm<-c("R15c", "B13", "B2", "O23", "O51", "O8", "R44", "R31", "R5", "Y1", "Y26")
data<-data[!(data$animal %in% rm),]
data$animal<- as.factor(as.character(data$animal))
immuneg<- immuneg[!(immuneg$id %in% rm),]
immuneg$id<-as.factor(as.character(immuneg$id))
unique(immuneg$id)
unique(data$animal)

maxtime<-tapply(data$stop2, data$animal, max)
maxitime<-tapply(immuneg$capturetime, immuneg$id, max)
df<-NA
for (i in 1:length(maxtime)){
     df[i]<-maxtime[[1]]>maxitime[[1]]
 }
immuneg$capturetime2<- immuneg$capturetime+
hapto.fit<-lme(log(hapto)~tb+age_yr+capturetime2, random = ~capturetime2|id, data=immuneg)
mort.fit = coxph(Surv(start2, stop2, death.time)~brucella+TB_3+herd2 + age_yr2 
                 +I(age_yr2^2) + cluster(animal), data=data,  x=TRUE, model=TRUE)
test<-jointModel(hapto.fit, mort.fit, timeVar='capturetime2')


# Does the association betwee current gamma levels and mortality vary by infection status
############################################




# Is there an association between change in gamma leves and mortality? or lagged gamma levels? 
############################################


############################################
############################################
############################################



############################################
############################################
############################################


############################################
############################################
######       Incidence Analyses        ######
############################################
############################################
# time until conversion: 
br<-read.csv("")
tb<-read.csv("")

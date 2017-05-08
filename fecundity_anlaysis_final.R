##############################################
##############################################
# Fecundity analyses- Calving
##############################################
##############################################
library(lme4)
library(arm)
library(ggplot2)
library(MASS) #glmmpql
library(MuMIn)
########################################################################
# Read in data, calculate approximation for calving, in calving period
########################################################################
data<- read.csv("~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/final_datasets_copied_from_phdfolder/cross_sectional_data_withdz_cleandisease_nofinal_Feb2016_capturetime_forsurv.csv")
data2<- data[data$herdorig=="LS",]
calftime <- c(0, 12, 24, 36)
data4<- data2[data2$capturetime %in% calftime,]
data4$fec<- NA
table(data4$milk , data4$calf)
# With milk one 10 year old, two 12 year olds, no 9, 11, 0, or 1 yr olds (good)
# one 11 year old without milk
# called age <=10 yes
# called 10 unknowns with milk yes; 3 unknowns with no/unknown milks no
include = c( 10, 2, 3,4, 5, 6, 7, 8, 9, "yes")
data4$fec[data4$calf %in% include] <- 1 
data4$fec[data4$calf %in% c(11, 12, "no")] <- 0
data4$fec[data4$calf=="unknown" & data4$milk=="yes"] <- 1
data4$fec[data4$calf=="unknown" & data4$milk=="unknown"] <- 0
data4$fec[data4$calf=="unknown" & data4$milk=="no"] <- 0

# overall rates
table(data4$fec, data4$tb, data4$bruc)
# 19.4 in uninfected
# 28% in bruc -; bTB + 
# 20.7% in bruc +; bTB - 
# 33% in co


pregdata<- data2[!(data2$capturetime %in% calftime),]
pregdata$pregmilk <- paste(pregdata$preg, pregdata$milk)
table(pregdata$pregmilk, pregdata$tb, pregdata$bruc)  
(23+10+4)/(23+10+58) # four not pregnants had milk, so assume they just gave birth...  made yes 

# plot acf plot
acf(data4$fec)
pacf(data4$fec)

########################################################################
# NEED TO ACCOUNT FOR AGE- 
# Figures show no calves in buffalo < 4yrs and no TB in buffalo 9+ years, so analyses subset to capture these age categories. 
########################################################################
par(mfrow=c(2,2))
hist(as.numeric(data4$age_yr[data4$tb==0 & data4$bruc=="negative"]), main="Uninfected", xlab= "Age")
hist(data4$age_yr[data4$tb==0 & data4$bruc=="positive"], main = "Brucellosis + ", xlab= "Age")
hist(data4$age_yr[data4$tb==1 & data4$bruc=="negative"], main = "TB +", xlab= "Age")
hist(data4$age_yr[data4$tb==1 & data4$bruc=="positive"], main = "Coinfected", xlab= "Age")

d<- data4
d$tb[d$tb==0] <- "negative"
d$tb[d$tb==1] <- "positive"
d$tb <- as.factor(d$tb)
a<- data.frame(table(d$fec, d$tb, d$bruc,  d$age_yr)); colnames(a)<- c("Fec", "TB", "Brucellosis", "Age", "Freq")
a$col <- c("blue4", "brown4", "dodgerblue", "brown1", "darkslategray3", "lightcoral", "lightskyblue", "mistyrose1")
# preg= red, not= blue

# make a stacked barplot
make_stacked_barplot = function(age){
	temp<- a[a$Age==age,]
	temp$Brucellosis<- relevel(temp$Brucellosis, "negative")
	temp$TB<- relevel(temp$TB, "negative")
	a<- temp[with(temp, order(temp$Fec, temp$Brucellosis, temp$TB)),]
	mat<- matrix(a$Freq, nrow=2, ncol=4, byrow=TRUE, 
	dimnames=list(c("no calf", "calf"), c("Uninfected", "TB+", "Brucellosis+", "Co-infected")))
	prop <- prop.table(mat, margin=2)
	par(mar=c(5.1, 4.1, 4.1, 7.1), xpd=TRUE)
	barplot(mat, col=heat.colors(length(rownames(mat))), width=2, main=paste("Age= ", age), las=1, 
	ylab= "Number of buffalo")
	legend("topright", fill=heat.colors(length(rownames(mat))), legend=rownames(mat),bty="n")
}	

#ages<- c(2, 3,4,5,6,7,8, 9, 10)
ages<- c(4,5,6,7,8, 9)
par(mfrow=c(2,3))
for (age in ages){
	make_stacked_barplot(age)
}

########################################################################
#Choose age category that is most appropriate
########################################################################
d<- data4[data4$age_yr > 3 & data4$age_yr <= 9,]
table(d$fec, d$tb, d$bruc)
# uninfected: 31.4%
# bTB only: 0.3
# bruc : 24.4
# co: 33%

# with > 4
# uninfected: 56%
# bTB only: 30.4%
# bruc only: 30%
# co: 42.8%


d$age5 <- "adult"; 
d$age5[d$age_yr %in% c(4) ] <- "juvenile"
d$age5[d$age_yr %in% c(8, 9) ] <- "mature"

d$age6 <- "adult"; 
d$age6[d$age_yr %in% c(4, 5) ] <- "juvenile"
d$age6[d$age_yr %in% c(8, 9) ] <- "mature"

d$age1 <- "adult" 
d$age1[d$age_yr %in% c(4) ] <- "juvenile"

d$age3 <- "adult" 
d$age3[d$age_yr %in% c(4, 5) ] <- "juvenile"

d$age7<- "adult"
d$age7[d$age_sel/12 < 5] <- "subadult"
d$age7[d$age_sel/12 > 7] <- "mature"

# FIX ME- THINK ABOUT BETTER PACKAGE... 
t1<-glmmPQL(fec~ age1, correlation= corAR1(form=~ capturetime |id), random= ~ 1|id, data=d, family=binomial); summary(t1); r.squaredGLMM(t1)
t1<-glmmPQL(fec~ age3, correlation= corAR1(form=~ capturetime |id), random= ~ 1|id, data=d, family=binomial); summary(t1); r.squaredGLMM(t1)
t1<-glmmPQL(fec~ age5, correlation= corAR1(form=~ capturetime |id), random= ~ 1|id, data=d, family=binomial); summary(t1) #mature didn't improve
t1<-glmmPQL(fec~ age6, correlation= corAR1(form=~ capturetime |id), random= ~ 1|id, data=d, family=binomial); summary(t1); r.squaredGLMM(t1)
t1<-glmmPQL(fec~ age7, correlation= corAR1(form=~ capturetime |id), random= ~ 1|id, data=d, family=binomial); summary(t1); r.squaredGLMM(t1)

# age 1 and age 3 prefered. Going with age1 = 4 and younger
# age 1 = c(4, rest) higher values than age 3 = c(4,5), rest
# age 6 = equivalent of age3 with a mature stage. 
age1: 0.2866831 0.6014588  **
age5: 0.3147492 0.5651035  # age 1 with matrue stage

age3: 0.2176504 0.5318982 
age6: 0.2256599 0.5171490 # age 3 with a mature stage. 
age7: 0.4188140 0.5353242 

t1<-glmmPQL(fec~ age1+tb + bruc, correlation= corAR1(form=~ capturetime |id), random= ~ 1|id, data=d, family=binomial); summary(t1); r.squaredGLMM(t1)
t1<-glmmPQL(fec~ age3+tb + bruc, correlation= corAR1(form=~ capturetime |id), random= ~ 1|id, data=d, family=binomial); summary(t1); r.squaredGLMM(t1)
t1<-glmmPQL(fec~ age5+tb + bruc, correlation= corAR1(form=~ capturetime |id), random= ~ 1|id, data=d, family=binomial); summary(t1) #mature didn't imporve
t1<-glmmPQL(fec~ age6+tb + bruc, correlation= corAR1(form=~ capturetime |id), random= ~ 1|id, data=d, family=binomial); summary(t1); r.squaredGLMM(t1)
t1<-glmmPQL(fec~ age7+tb + bruc, correlation= corAR1(form=~ capturetime |id), random= ~ 1|id, data=d, family=binomial); summary(t1); r.squaredGLMM(t1)

age1: 0.2907100 0.6482049 **
age5: 0.3252764 0.6080539

age3: 0.2286035 0.5511505
age6: 0.2425155 0.5270254 
age7: 0.4362000 0.5852685 

table(d$fec, d$tb, d$bruc, d$age1)


t1<-glmmPQL(fec~ age1*tb + age1*bruc, correlation= corAR1(form=~ capturetime |id), random= ~ 1|id, data=d, family=binomial); summary(t1)
t1<-glmmPQL(fec~ age1*tb + bruc, correlation= corAR1(form=~ capturetime |id), random= ~ 1|id, data=d, family=binomial); summary(t1)
t1<-glmmPQL(fec~ age1*tb + tb*bruc, correlation= corAR1(form=~ capturetime |id), random= ~ 1|id, data=d, family=binomial); summary(t1)
t1<-glmmPQL(fec~ age1 + tb*bruc, correlation= corAR1(form=~ capturetime |id), random= ~ 1|id, data=d, family=binomial); summary(t1)
t1<-glmmPQL(fec~ age1*bruc + tb*bruc, correlation= corAR1(form=~ capturetime |id), random= ~ 1|id, data=d, family=binomial); summary(t1)
###############
t1<-glmmPQL(fec~ age3*tb + age3*bruc, correlation= corAR1(form=~ capturetime |id), random= ~ 1|id, data=d, family=binomial); summary(t1)
t1<-glmmPQL(fec~ age3*tb + bruc, correlation= corAR1(form=~ capturetime |id), random= ~ 1|id, data=d, family=binomial); summary(t1)
t1<-glmmPQL(fec~ age3*tb + tb*bruc, correlation= corAR1(form=~ capturetime |id), random= ~ 1|id, data=d, family=binomial); summary(t1)
t1<-glmmPQL(fec~ age3 + tb*bruc, correlation= corAR1(form=~ capturetime |id), random= ~ 1|id, data=d, family=binomial); summary(t1)
t1<-glmmPQL(fec~ age3*bruc + tb*bruc, correlation= corAR1(form=~ capturetime |id), random= ~ 1|id, data=d, family=binomial); summary(t1)

t1<-glmmPQL(fec~ age1*bruc + tb*bruc + age1*tb, correlation= corAR1(form=~ capturetime |id), random= ~ 1|id, data=d, family=binomial); summary(t1)

# TO REPORT:                     
#Fixed effects: fec ~ age1 + tb * bruc + tb * age1 + bruc * age1 
#                              Value Std.Error DF   t-value p-value
#(Intercept)                0.343779 0.4832828 75  0.711341  0.4791
#age1juvenile              -2.846057 0.7897330 75 -3.603822  0.0006
#tb                        -1.619684 0.7148484 75 -2.265773  0.0264
#brucpositive              -1.504051 0.6718430 75 -2.238694  0.0281
#tb:brucpositive            2.283633 0.9866149 75  2.314614  0.0234
#age1juvenile:tb            3.615669 1.2631557 75  2.862409  0.0054
#age1juvenile:brucpositive -2.623061 1.4908507 75 -1.759439  0.0826

d$age1<- as.factor(d$age1)
d$age1j<-relevel(d$age1, "juvenile")
t1<-glmmPQL(fec~ age1j*bruc + tb*bruc + age1j*tb, correlation= corAR1(form=~ capturetime |id), random= ~ 1|id, data=d, family=binomial); summary(t1)






# CI, use t-distribution on 75 df, 
# brucellosis effect
exp(-1.504)
exp (-1.504 - 1.98*0.6718430); exp (-1.504 + 1.98*0.6718430)
# tb
exp(-1.6197)
exp(-1.6197- 1.98* 0.7148484); exp(-1.6197+ 1.98* 0.7148484)
# co
exp(-1.504- 1.6197 + 2.2836)
d$tb2<- as.factor(d$tb)
d$tb2<- relevel(d$tb2, "1")
d$bruc2 <- as.factor(d$bruc)
d$bruc2 <- relevel(d$bruc2, "positive")
t1<-glmmPQL(fec~ age1*bruc2 + tb2*bruc2 + age1*tb2, correlation= corAR1(form=~ capturetime |id), random= ~ 1|id, data=d, family=binomial); summary(t1)

########################################################################
Consider brucellosis time since infection 
########################################################################
length(d[,1])
dsave<- d[d$bruc_beforeafter != "pfc",]
length(dsave[,1])
dsave$tb[dsave$tb==0] <- "negative"
dsave$tb[dsave$tb==1] <- "positive"
dsave$tb <- as.factor(dsave$tb)

a<- data.frame(table(dsave$fec, dsave$tb, dsave$bruc,  dsave$age_yr)); colnames(a)<- c("Fec", "TB", "Brucellosis", "Age", "Freq")
a$col <- c("blue4", "brown4", "dodgerblue", "brown1", "darkslategray3", "lightcoral", "lightskyblue", "mistyrose1")
ages<- c(4,5,6,7,8, 9)
par(mfrow=c(2,3))
for (age in ages){
	make_stacked_barplot(age)
}


#write.csv(dsave, "~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/final_datasets_copied_from_phdfolder/fec_timesincebruc.csv")
# don't really have more than two captures past incidence, so can't get at what happens over two years post conversion.


t1<-glmmPQL(fec~ age1*bruc + tb*bruc + age1*tb, correlation= corAR1(form=~ capturetime |id), random= ~ 1|id, data=dsave, family=binomial); summary(t1)
t1<-glmmPQL(fec~ age1 + bruc + tb*bruc + age1*tb, correlation= corAR1(form=~ capturetime |id), random= ~ 1|id, data=dsave, family=binomial); summary(t1)
t1<-glmmPQL(fec~ age1 + bruc + tb + age1*tb, correlation= corAR1(form=~ capturetime |id), random= ~ 1|id, data=dsave, family=binomial); summary(t1)
#on only  those where we know time since first infection
t1<-glmmPQL(fec~ age1 + bruc + tb + age1*tb, correlation= corAR1(form=~ capturetime |id), random= ~ 1|id, data=dsave, family=binomial); summary(t1)


########################################################################
Above non-conclusive... now consider model with brucellosis as pfc vs recent infection
########################################################################

# Including pfc as a covariate improves R2 slightly.  In this simple model, converters have reduced fecundity and pfc animals do not.  
d$bruc_beforeafter[d$bruc_beforeafter=="nc"] <- 0
t1<-glmmPQL(fec~ age1 + bruc + tb + age1*tb, correlation= corAR1(form=~ capturetime |id), random= ~ 1|id, data=d, family=binomial); summary(t1)
r.squaredGLMM(t1)  #0.3816404 0.7114846 
t1<-glmmPQL(fec~ age1 + bruc_beforeafter + tb + age1*tb, correlation= corAR1(form=~ capturetime |id), random= ~ 1|id, data=d, family=binomial); summary(t1)
r.squaredGLMM(t1)  #0.4317796 0.7539974 


# In full model, R2 is improved much more by including pfcs.
# Now pfc's are negatively associated with fecundity unless co-infected; converters n.s. 
# Model not fitting well because pfc animals are older... so can't really say much at all. 
t1<-glmmPQL(fec~ age1*bruc + tb*bruc + age1*tb, correlation= corAR1(form=~ capturetime |id), random= ~ 1|id, data=d, family=binomial); summary(t1)
r.squaredGLMM(t1) # 0.519, 0.8666
t1<-glmmPQL(fec~ age1*bruc_beforeafter + tb*bruc_beforeafter + age1*tb, correlation= corAR1(form=~ capturetime |id), random= ~ 1|id, data=d, family=binomial); summary(t1)
r.squaredGLMM(t1)  #0.9145931 0.9875749 
 t1<-glmmPQL(fec~ age1+bruc_beforeafter + tb*bruc_beforeafter + age1*tb, correlation= corAR1(form=~ capturetime |id), random= ~ 1|id, data=d, family=binomial); summary(t1)


a<- data.frame(table(dsave$fec, dsave$tb, dsave$bruc_beforeafter,  dsave$age_yr)); colnames(a)<- c("Fec", "TB", "Brucellosis", "Age", "Freq")
a$Brucellosis<- as.character(a$Brucellosis)
a$Brucellosis[a$Brucellosis=="nc"] <- "negative"
a$Brucellosis[a$Brucellosis =="0"] <- "negative"
a$Brucellosis[a$Brucellosis =="1"] <- "converter"
a$Brucellosis<-as.factor(a$Brucellosis)

b<- a[a$Age != 4,]

# make a new data frame to hold updated frequenceis wihtout the ages...
c<- b[b$Age == 5,]
c <- c[c(1:8, 13:16),]
c$Freq <- 0
for (i in 1:length(c$Freq)){
	tbval <- c$TB[i]
	brucval <- c$Brucellosis[i]
	ageval <- c$Age[i]
	fecval <- c$Fec[i]
	temp <- b[b$TB == tbval & b$Brucellosis == brucval & b$Age == ageval & b$Fec == fecval,]	
	#if (length(temp$Freq) >1) { 
	#	print ("subset of length longer than 1, ")
	#	print ( temp)}
	#else {
		c$Freq[i] <- c$Freq[i] + sum(temp$Freq)
	#}
}


make_stacked_barplot = function(age){
	temp<- b
	temp$Brucellosis<- relevel(temp$Brucellosis, "negative")
	temp$TB<- relevel(temp$TB, "negative")
	a<- temp[with(temp, order(temp$Fec, temp$Brucellosis, temp$TB)),]


	mat<- matrix(a$Freq, nrow=2, ncol=6, byrow=TRUE,dimnames=list(c("no calf", "calf"), c("Uninfected", "TB+", "Brucellosis+", "BrucellosisConverter", "Co-infected+", "Co-Converter")))
	prop <- prop.table(mat, margin=2)
	par(mar=c(5.1, 4.1, 4.1, 7.1), xpd=TRUE)
	barplot(mat, col=heat.colors(length(rownames(mat))), width=2, main=paste("Age= ", age), las=1, 
	ylab= "Number of buffalo")
	legend("topright", fill=heat.colors(length(rownames(mat))), legend=rownames(mat),bty="n")
}	

age<- c(5,6,7,8, 9)
par(mfrow=c(2,3))
for (age in ages){
	make_stacked_barplot(age)
}




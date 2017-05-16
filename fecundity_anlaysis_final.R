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
# Figures show no calves in buffalo < 4yrs and no TB in buffalo 9+ years,
# so analyses subset to capture these age categories. 
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

# Define age categories
########################################################################
d$age1 <- "adult"; 
d$age1[d$age_yr < 6.55] <- "juvenile"

d$age3 <- "adult"; 
d$age3[d$age_yr < 5.55] <- "juvenile"

# Do we need random intercept and random time component or AR1
########################################################################
t1<-glmmPQL(fec~ age1+ tb + bruc, correlation= corAR1(form=~ capturetime |id), random= ~ 1|id, data=d, family=binomial); summary(t1); r.squaredGLMM(t1) # phi estimated at 0...

t <- glmer(d$fec ~ d$age1 + d$tb + d$bruc + d$capturetime+ (1|d$id), family = binomial(link = "logit")); summary(t)
r <- resid(t)
par(mfrow = c(1,2))
acf(r)
pacf(r)

# not much variaince from random slopes, no linear trend with time
t <- glmer(d$fec ~ d$age1 + d$tb + d$bruc + d$capturetime+ (d$capturetime|d$id), family = binomial(link = "logit")); summary(t)

# Model selection:
########################################################################
t <- glmer(d$fec ~ d$age1 + (1|d$id), family = binomial(link = "logit")); summary(t); extractAIC(t); logLik(t) #166.4655; -80.23276
t <- glmer(d$fec ~ d$age3 + (1|d$id), family = binomial(link = "logit")); summary(t); extractAIC(t); logLik(t) #164.6394; -79.31971

# for figure B2
t<- d[d$tb == 1,]
table(t$fec, t$bruc, t$age3)
t<- d[d$tb == 0,]
table(t$fec, t$bruc, t$age3)

#d$age3 <- as.factor(d$age3)
#d$age3 <- relevel(d$age3, "juvenile")

full.mod <- glmer(d$fec ~ d$age3 + d$tb + d$bruc + d$tb*d$age3 + d$bruc*d$age3 + d$tb*d$bruc + (1|d$id), family = binomial(link = "logit")); summary(full.mod) #164.6394; -79.31971
lwm <- glmer(d$fec ~ d$age3 + d$tb + d$bruc + d$tb*d$age3 + d$tb*d$bruc + (1|d$id), family = binomial(link = "logit")); summary(lwm) 
lwm <- glmer(d$fec ~ d$age3 + d$tb + d$bruc + d$tb*d$bruc + (1|d$id), family = binomial(link = "logit")); summary(lwm)
lwm <- glmer(d$fec ~ d$age3 + d$tb + d$bruc + (1|d$id), family = binomial(link = "logit")); summary(lwm)
lwm <- glmer(d$fec ~ d$age3 + d$bruc + (1|d$id), family = binomial(link = "logit")); summary(lwm)
lwm <- glmer(d$fec ~ d$age3 + (1|d$id), family = binomial(link = "logit")); summary(lwm)

t <- d[d$age3 == "adult",]
lwm <- glmer(t$fec ~ t$age3 + t$tb*t$bruc + (1|t$id), family = binomial(link = "logit")); summary(lwm)


# Play with age categories ... not used!
########################################################################
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

##############################################
##############################################
# Fecundity analyses- Calving
##############################################
##############################################
library(lme4)
library(arm)
library(ggplot2)

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

# overall rates: 
table(data4$fec, data4$tb, data4$bruc)

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
	barplot(mat, col=heat.colors(length(rownames(mat))), width=2, main=paste("Age= ", age))
	legend("topright", fill=heat.colors(length(rownames(mat))), legend=rownames(mat),bty="n")
}	

ages<- c(2, 3,4,5,6,7,8, 9, 10)
par(mfrow=c(3,3))
for (age in ages){
	make_stacked_barplot(age)
}

########################################################################
Choose age category that is most appropriate
########################################################################
d<- data4[data4$age_yr > 3 & data4$age_yr <= 9,]
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

t1<-glmmPQL(fec~ age1, correlation= corAR1(form=~ capturetime |id), random= ~ 1|id, data=d, family=binomial); summary(t1) 
t1<-glmmPQL(fec~ age3, correlation= corAR1(form=~ capturetime |id), random= ~ 1|id, data=d, family=binomial); summary(t1)
t1<-glmmPQL(fec~ age5, correlation= corAR1(form=~ capturetime |id), random= ~ 1|id, data=d, family=binomial); summary(t1) #mature didn't improve
t1<-glmmPQL(fec~ age6, correlation= corAR1(form=~ capturetime |id), random= ~ 1|id, data=d, family=binomial); summary(t1)
t1<-glmmPQL(fec~ age7, correlation= corAR1(form=~ capturetime |id), random= ~ 1|id, data=d, family=binomial); summary(t1)
# age 1 and age 3 prefered. Going with age1 = 4 and younger

t1<-glmmPQL(fec~ age1+tb + bruc, correlation= corAR1(form=~ capturetime |id), random= ~ 1|id, data=d, family=binomial); summary(t1)
t1<-glmmPQL(fec~ age3+tb + bruc, correlation= corAR1(form=~ capturetime |id), random= ~ 1|id, data=d, family=binomial); summary(t1)
t1<-glmmPQL(fec~ age5+tb + bruc, correlation= corAR1(form=~ capturetime |id), random= ~ 1|id, data=d, family=binomial); summary(t1) #mature didn't imporve
t1<-glmmPQL(fec~ age6+tb + bruc, correlation= corAR1(form=~ capturetime |id), random= ~ 1|id, data=d, family=binomial); summary(t1)
t1<-glmmPQL(fec~ age7+tb + bruc, correlation= corAR1(form=~ capturetime |id), random= ~ 1|id, data=d, family=binomial); summary(t1)


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

# TO REPORT:                     Value Std.Error DF   t-value p-value
(Intercept)      0.336268 0.4444201 76  0.756645  0.4516
age1juvenile    -2.880588 0.7645851 76 -3.767518  0.0003
tb              -1.234166 0.6440562 76 -1.916239  0.0591
brucpositive    -1.410167 0.6099359 76 -2.311992  0.0235
tb:brucpositive  1.588655 0.8786644 76  1.808034  0.0746
age1juvenile:tb  2.366394 1.0623551 76  2.227498  0.0289

Fixed effects: fec ~ age1 + tb * bruc + tb * age1 + bruc * age1 
                              Value Std.Error DF   t-value p-value
(Intercept)                0.343779 0.4832828 75  0.711341  0.4791
age1juvenile              -2.846057 0.7897330 75 -3.603822  0.0006
tb                        -1.619684 0.7148484 75 -2.265773  0.0264
brucpositive              -1.504051 0.6718430 75 -2.238694  0.0281
tb:brucpositive            2.283633 0.9866149 75  2.314614  0.0234
age1juvenile:tb            3.615669 1.2631557 75  2.862409  0.0054
age1juvenile:brucpositive -2.623061 1.4908507 75 -1.759439  0.0826







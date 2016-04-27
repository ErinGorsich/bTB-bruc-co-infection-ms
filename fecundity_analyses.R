##############################################
##############################################
# Fecundity analyses
##############################################
##############################################
data<- read.csv("~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/final_datasets_copied_from_phdfolder/cross_sectional_data_withdz_cleandisease_nofinal_Feb2016_capturetime_forsurv.csv")
library(lme4)
library(arm)
library(ggplot2)

######################################################
# 1 Proportion of pregnant buffalo in LS or with <=2mo calves 
#   During October-January, the capture right before the birthing season
######################################################
data2<- data[data$herdorig=="LS",]
table(data2$capturetime)
# 0  6  9 12 15 18 24 30 33 36 39 42 48 
# 50 49  2 52  1 50 51 52  1 51  1 47  1 
# October -January cooresponds to capturetime= 
pregtime <- c(6, 18, 30, 42)
data3<- data2[data2$capturetime %in% pregtime,]

table(data3$preg, data3$calf)

# fecundity==1 if mid, late mid/late pregnant. 
# two early pregnancies, both not counted because early... will be for next year. 
data3$fec<- NA
data3$fec[data3$preg == "late" || data3$preg=="mid"|| data3$preg=="mid/late"]<- 1
data3$fec[data3$preg == "early" || data3$preg == "unknown"] <- 0
data3$fec[data3$preg == "no" & data3$calf %in% c(0, 1, 2)] <- 1
data3$fec[data3$preg == "no" & data3$calf %in% c(8, 10, 11, 12, "no", "yes", "unknown")] <- 0

# overall rates: 
table(data3$fec, data3$tb, data3$bruc)
# 90 TB - Br -
# 30 TB + Br-
# 51 TB - BR+
# 24 TB + BR +

# BR - TB -
34/(56+34)   # 0.37777
# Br - TB + 
17/(17+16)   # 0.515

# BR + TB -
23/(28+23)  # 0.4509
# BR+ TB +
5/(5+19)    # 0.208

# NEED TO ACCOUNT FOR AGE!
par(mfrow=c(2,2))
hist(as.numeric(data3$age_yr[data3$tb==0 & data3$bruc=="negative"]), main="Uninfected", xlab= "Age")
hist(data3$age_yr[data3$tb==0 & data3$bruc=="positive"], main = "Brucellosis + ", xlab= "Age")
hist(data3$age_yr[data3$tb==1 & data3$bruc=="negative"], main = "TB +", xlab= "Age")
hist(data3$age_yr[data3$tb==1 & data3$bruc=="positive"], main = "Coinfected", xlab= "Age")

# youngest pregnant buffalo was 2.6yr, it is the only pregnant buffalo less than 3. 
par(mfrow=c(2,2))
hist(as.numeric(data3$age_yr[data3$tb==0 & data3$bruc=="negative" & data3$fec==1]), main="Uninfected & Pregnant", xlab= "Age")
hist(data3$age_yr[data3$tb==0 & data3$bruc=="positive" & data3$fec==1], main = "Brucellosis +  & Pregnant", xlab= "Age")
hist(data3$age_yr[data3$tb==1 & data3$bruc=="negative" & data3$fec==1], main = "TB + & Pregnant", xlab= "Age")
hist(data3$age_yr[data3$tb==1 & data3$bruc=="positive" & data3$fec==1], main = "Coinfected & Pregnant", xlab= "Age")



# age is sig, neither infection alone (but neg), sig positve on interaction term 
test.mod<-glmer(data3$fec ~ age_yr + tb + bruc+ (1|id), data=data3, family=binomial(link="logit"))
test.mod<-glmer(data3$fec ~ age_yr+ I(age_yr^2) + tb + bruc+ (1|id), data=data3, family=binomial(link="logit")); summary(test.mod)
test.mod<-glmer(data3$fec ~ age_yr+ I(age_yr^2) + tb *bruc+ (1|id), data=data3, family=binomial(link="logit")); summary(test.mod)  # 234.2

# age*tb (0.05); age*bruc n.s.
test.mod<-glmer(data3$fec ~ age_yr*bruc+ I(age_yr^2) + tb + bruc+ (1|id), data=data3, family=binomial(link="logit")); summary(test.mod)
test.mod<-glmer(data3$fec ~ age_yr*tb+ I(age_yr^2) + tb + bruc+ (1|id), data=data3, family=binomial(link="logit")); summary(test.mod)

# no shift in both terms
test.mod<-glmer(data3$fec ~ age_yr*bruc+ I(age_yr^2)*bruc + tb + bruc+ (1|id), data=data3, family=binomial(link="logit")); summary(test.mod) 
test.mod<-glmer(data3$fec ~ age_yr*tb+ I(age_yr^2)*tb + tb + bruc+ (1|id), data=data3, family=binomial(link="logit")); summary(test.mod)

# age^2 by tb significant (0.06)
test.mod<-glmer(data3$fec ~ age_yr+ I(age_yr^2)*bruc + tb + bruc+ (1|id), data=data3, family=binomial(link="logit")); summary(test.mod) # AIC = 174.3
test.mod<-glmer(data3$fec ~ age_yr+ I(age_yr^2)*tb + tb + bruc+ (1|id), data=data3, family=binomial(link="logit")); summary(test.mod)

# no age*bruc patterns, just tb*bruc in all
test.mod<-glmer(data3$fec ~ age_yr*bruc+ I(age_yr^2) + tb *bruc+ (1|id), data=data3, family=binomial(link="logit")); summary(test.mod)  # 172.7
test.mod<-glmer(data3$fec ~ age_yr*bruc+ I(age_yr^2)*bruc + tb *bruc+ (1|id), data=data3, family=binomial(link="logit")); summary(test.mod)
test.mod<-glmer(data3$fec ~ age_yr+ I(age_yr^2)*bruc + tb *bruc+ (1|id), data=data3, family=binomial(link="logit")); summary(test.mod)

# tb*age terms not sig, but suggests a shit interaction with both most supported.
****test.mod<-glmer(data3$fec ~ age_yr*tb+ I(age_yr^2) + tb *bruc+ (1|id), data=data3, family=binomial(link="logit")); summary(test.mod) # 230.3
test.mod<-glmer(data3$fec ~ age_yr*tb+ I(age_yr^2)*tb + tb *bruc+ (1|id), data=data3, family=binomial(link="logit")); summary(test.mod)
test.mod<-glmer(data3$fec ~ age_yr+ I(age_yr^2)*tb + tb *bruc+ (1|id), data=data3, family=binomial(link="logit")); summary(test.mod)  #  230.6


final.mod<-glmer(data3$fec ~ age_yr*tb+ I(age_yr^2) + tb *bruc+ (capturetime|id), data=data3, family=binomial(link="logit")); summary(test.mod)
#Fixed effects:
#                Estimate Std. Error z value Pr(>|z|)    
#(Intercept)     -8.50652    1.62722  -5.228 1.72e-07 ***
#age_yr           2.49164    0.51788   4.811 1.50e-06 ***
#tb               2.96936    1.37984   2.152  0.03140 *  
#I(age_yr^2)     -0.15091    0.03596  -4.197 2.71e-05 ***
#brucpositive    -0.89174    0.46157  -1.932  0.05336 .  
#age_yr:tb       -0.59721    0.24206  -2.467  0.01362 *  
#tb:brucpositive  2.17197    0.78507   2.767  0.00566 ** 


data3$age_yr2<-NA
for (i in 1:length(data3$age_yr)){
data3$age_yr2[i]<- (data3$age_yr[i]- mean(data3$age_yr))/sd(data3$age_yr)}

final.mod<-glmer(data3$fec ~ age_yr2*tb+ I(age_yr2^2) + tb *bruc+ (capturetime|id), data=data3, family=binomial(link="logit")); summary(final.mod)

#Fixed effects:
#                Estimate Std. Error z value Pr(>|z|)    
#(Intercept)       0.5815     0.3208   1.813  0.06989 .  
#age_yr2           1.7516     0.3196   5.480 4.25e-08 ***
#tb               -0.2791     0.4630  -0.603  0.54658    
#I(age_yr2^2)     -0.6410     0.1528  -4.196 2.72e-05 ***
#brucpositive     -0.8917     0.4616  -1.932  0.05337 .  
#age_yr2:tb       -1.2308     0.4989  -2.467  0.01362 *  
#tb:brucpositive   2.1720     0.7851   2.767  0.00566 **


# BUT...youngest pregnant buffalo was 2.6yr, it is the only pregnant buffalo less than 3. 
# no TB + buffalo older than 9... 
d <- data3[data3$age_yr >= 3,]
d$age_yr2 <- d$age_yr- 3  # puts zero point on buffalo that are age 3

final.mod<-glmer(fec ~ age_yr2*tb+ I(age_yr2^2) + tb *bruc+ (capturetime|id), data=d, family=binomial(link="logit")); summary(final.mod)

#                Estimate Std. Error z value Pr(>|z|)    
#(Intercept)      0.17836    0.29003   0.615  0.53856    
#age_yr2          0.96966    0.18321   5.293 1.21e-07 ***
#tb              -0.01828    0.47053  -0.039  0.96901    
#I(age_yr2^2)    -0.14859    0.03650  -4.071 4.69e-05 ***
#brucpositive    -0.88487    0.46067  -1.921  0.05475 .  
#age_yr2:tb      -0.58896    0.24298  -2.424  0.01536 *  
#tb:brucpositive  2.16352    0.78432   2.758  0.00581 ** 
# says odds of pregnancy in 3 year old uninfected buffalo is 0.178 +/- 1.96*0.29
#0.1511036

# refit at other ages... 
ages <- seq(3, 8, 1)
fitsummary <- data.frame(Age.ind=rep(ages, 4), Coef= NA, SE= NA, LB=NA, UB=NA, Z= NA, p= NA, TB= NA, Bruc= NA)
#d <- data3[data3$age_yr >= 3,]

d<- data3[data3$age_yr >= 3 & data3$age_yr <= 8,]
d$age_yr2<- d$age_yr - 4
#final.mod<-glmer(fec ~ age_yr2*tb+ I(age_yr2^2) + tb *bruc+ (capturetime|id), data=d, family=binomial(link="logit")); summary(final.mod)
#final.mod<-glmer(fec ~ age_yr2+ I(age_yr2^2) + tb *bruc+ (capturetime|id), data=d, family=binomial(link="logit")); summary(final.mod)

# this is why: 
#table(d$fec, d$tb, d$age_yr)
d$tb[d$tb==0] <- "negative"
d$tb[d$tb==1] <- "positive"
d$tb <- as.factor(d$tb)

get_model_statistics= function(data, age){
	data$age_yr2 <- data$age_yr - age
	#temp.mod<- glmer(fec ~ age_yr*tb+ I(age_yr^2) + tb *bruc+ (capturetime|id), data=data, family=binomial(link="logit"))
	temp.mod<- glmer(fec ~ age_yr2 + I(age_yr2^2) + tb *bruc+ (capturetime|id), data=data, family=binomial(link="logit"))
	Coef <- temp.mod@beta[1]
	SE <- se.fixef(temp.mod)[[1]]
	LB <- Coef - 1.96*SE
	UB <- Coef + 1.96*SE 
	Z <- coef(summary(temp.mod))[1,3]
	p <- coef(summary(temp.mod))[1,4]
	rowvals <- c(Coef, SE, LB, UB, Z, p)
	return(rowvals)
	rm(temp.mod, Coef, SE, LB, UB, Z, p)
}

for (i in 1:length(ages)){
	# TB-, Bruc - 
	d$tb<- relevel(d$tb, "negative")
	d$bruc<- relevel(d$bruc, "negative")
	fitsummary[i,2:7] <- get_model_statistics(d, ages[i])
	fitsummary$TB[i] <- "negative"
	fitsummary$Bruc[i] <- "negative"

	# TB + Bruc -
	d$tb<- relevel(d$tb, "positive")
	d$bruc<- relevel(d$bruc, "negative")
	fitsummary[i+length(ages),2:7] <- get_model_statistics(d, ages[i])
	fitsummary$TB[i+length(ages)] <- "positive"
	fitsummary$Bruc[i+length(ages)] <- "negative"

	# TB + Bruc -
	d$tb<- relevel(d$tb, "negative")
	d$bruc<- relevel(d$bruc, "positive")
	fitsummary[i+2*length(ages),2:7] <- get_model_statistics(d, ages[i])
	fitsummary$TB[i+2*length(ages)] <- "negative"
	fitsummary$Bruc[i+2*length(ages)] <- "positive"

	# TB + Bruc +
	d$tb<- relevel(d$tb, "positive")
	d$bruc<- relevel(d$bruc, "positive")
	fitsummary[i+3*length(ages),2:7] <- get_model_statistics(d, ages[i])
	fitsummary$TB[i+3*length(ages)] <- "positive"
	fitsummary$Bruc[i+3*length(ages)] <- "positive"
}



# plot predicted values
fitsummary$prob<- NA; fitsummary$plb<- NA; fitsummary$pub <- NA
fitsummary$prob<-exp(fitsummary$Coef)/ (1+ exp(fitsummary$Coef))
fitsummary$plb<-exp(fitsummary$LB)/ (1+ exp(fitsummary$LB))
fitsummary$pub<-exp(fitsummary$UB)/ (1+ exp(fitsummary$UB))

#write.csv(fitsummary, "~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/final_datasets_copied_from_phdfolder/fecundity/prob_pregnant.csv")
fitsummary$dz<- NA
fitsummary$dz[fitsummary$TB=="negative" & fitsummary$Bruc=="negative"] <- "Uninfected"
fitsummary$dz[fitsummary$TB=="negative" & fitsummary$Bruc=="positive"] <- "Brucellosis"
fitsummary$dz[fitsummary$TB=="positive" & fitsummary$Bruc=="negative"] <- "TB"
fitsummary$dz[fitsummary$TB=="positive" & fitsummary$Bruc=="positive"] <- "Co-infected"

#ggplot(fitsummary, aes(x= fitsummary$Age.ind, y= fitsummary$Coef, colour = fitsummary$dz)) + geom_point()

#fs<- fitsummary[fitsummary$Age.ind < 9 & fitsummary$Age.ind > 3,]
fs <- fitsummary
ggplot(fs, aes(x= fs$Age.ind, y= fs$prob, colour = fs$dz)) + geom_line()
ggplot(fs, aes(x= fs$Age.ind, y= fs$prob, colour = fs$dz)) + geom_point() + geom_line()
#+ geom_errorbar(aes(ymin= fs$plb, ymax=fs$pub, width=0.1))
write.csv(fitsummary, "~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/final_datasets_copied_from_phdfolder/fecundity/prob_pregnant.csv")




d<- data3[data3$age_yr >= 3 & data3$age_yr <= 8,]
d<- data3
d$tb[d$tb==0] <- "negative"
d$tb[d$tb==1] <- "positive"
d$tb <- as.factor(d$tb)
ages<- c(3,4,5,6,7,8)
a<- data.frame(table(d$fec, d$tb, d$bruc,  d$age_yr)); colnames(a)<- c("Fec", "TB", "Brucellosis", "Age", "Freq")
a$col <- c("blue4", "brown4", "dodgerblue", "brown1", "darkslategray3", "lightcoral", "lightskyblue", "mistyrose1")
# preg= red, not= blue
a$concat<- paste(a$Fec, a$TB, a$Brucellosis, sep="_")
par(mfrow=c(2,3))
for (age in ages){
	temp<- a[a$Age==age,]
	barplot(temp$Freq, col=temp$col, main=paste("Age=", age, sep=" "))
}


# make a stacked barplot
make_stacked_barplot = function(age){
	temp<- a[a$Age==age,]
	temp$Brucellosis<- relevel(temp$Brucellosis, "negative")
	temp$TB<- relevel(temp$TB, "negative")
	a<- temp[with(temp, order(temp$Fec, temp$Brucellosis, temp$TB)),]
	mat<- matrix(a$Freq, nrow=2, ncol=4, byrow=TRUE, 
	dimnames=list(c("non-pregnant", "pregnant"), c("Uninfected", "TB+", "Brucellosis+", "Co-infected")))
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

######################################################
######################################################
# 1 Proportion of buffalo in LS or calves calves 
#   During May-August, the capture right after the birthing season
######################################################
######################################################
pregtime <- c(0, 12, 24, 36)
data4<- data2[data2$capturetime %in% pregtime,]
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

# BR - TB -
19/ (79+19)  # 19.39
# Br - TB + 
9/(9+23)  # 28.125

# BR + TB -
11/ (11+42)  # 20.75
# BR+ TB +
7/(7+14)   # 33.333

# NEED TO ACCOUNT FOR AGE!
par(mfrow=c(2,2))
hist(as.numeric(data4$age_yr[data4$tb==0 & data4$bruc=="negative"]), main="Uninfected", xlab= "Age")
hist(data4$age_yr[data4$tb==0 & data4$bruc=="positive"], main = "Brucellosis + ", xlab= "Age")
hist(data4$age_yr[data4$tb==1 & data4$bruc=="negative"], main = "TB +", xlab= "Age")
hist(data4$age_yr[data4$tb==1 & data4$bruc=="positive"], main = "Coinfected", xlab= "Age")






test.mod<-glmer(data4$fec ~ age_yr + tb + bruc+ (1|id), data=data4, family=binomial(link="logit"))
test.mod<-glmer(data4$fec ~ age_yr+ I(age_yr^2) + tb + bruc+ (1|id), data=data4, family=binomial(link="logit")); summary(test.mod)
test.mod<-glmer(data4$fec ~ age_yr+ I(age_yr^2) + tb *bruc+ (1|id), data=data4, family=binomial(link="logit")); summary(test.mod)  # 173.2

# no age*tb; suggestive age*bruc
test.mod<-glmer(data4$fec ~ age_yr*bruc+ I(age_yr^2) + tb + bruc+ (1|id), data=data4, family=binomial(link="logit")); summary(test.mod)
test.mod<-glmer(data4$fec ~ age_yr*tb+ I(age_yr^2) + tb + bruc+ (1|id), data=data4, family=binomial(link="logit")); summary(test.mod)

# no shift in both terms
test.mod<-glmer(data4$fec ~ age_yr*bruc+ I(age_yr^2)*bruc + tb + bruc+ (1|id), data=data4, family=binomial(link="logit")); summary(test.mod) 
test.mod<-glmer(data4$fec ~ age_yr*tb+ I(age_yr^2)*tb + tb + bruc+ (1|id), data=data4, family=binomial(link="logit")); summary(test.mod)

# some evidence for the bruc term (0.09)
**test.mod<-glmer(data4$fec ~ age_yr+ I(age_yr^2)*bruc + tb + bruc+ (1|id), data=data4, family=binomial(link="logit")); summary(test.mod) # AIC = 174.3
test.mod<-glmer(data4$fec ~ age_yr+ I(age_yr^2)*tb + tb + bruc+ (1|id), data=data4, family=binomial(link="logit")); summary(test.mod)

# best is int with age not age2:  age*bruc (0.125) alone with TB*bruc(0.068)
**test.mod<-glmer(data4$fec ~ age_yr*bruc+ I(age_yr^2) + tb *bruc+ (1|id), data=data4, family=binomial(link="logit")); summary(test.mod)  # 172.7
test.mod<-glmer(data4$fec ~ age_yr*bruc+ I(age_yr^2)*bruc + tb *bruc+ (1|id), data=data4, family=binomial(link="logit")); summary(test.mod)
test.mod<-glmer(data4$fec ~ age_yr+ I(age_yr^2)*bruc + tb *bruc+ (1|id), data=data4, family=binomial(link="logit")); summary(test.mod)

# tb*age terms not sig, but suggests a shit interaction with both most supported.
test.mod<-glmer(data4$fec ~ age_yr*tb+ I(age_yr^2) + tb *bruc+ (1|id), data=data4, family=binomial(link="logit")); summary(test.mod)
test.mod<-glmer(data4$fec ~ age_yr*tb+ I(age_yr^2)*tb + tb *bruc+ (1|id), data=data4, family=binomial(link="logit")); summary(test.mod)
test.mod<-glmer(data4$fec ~ age_yr+ I(age_yr^2)*tb + tb *bruc+ (1|id), data=data4, family=binomial(link="logit")); summary(test.mod)

final.mod<-glmer(data4$fec ~ age_yr+ I(age_yr^2) + tb *bruc+ (1|id), data=data4, family=binomial(link="logit")); summary(test.mod)
# Fixed effects:
#                Estimate Std. Error z value Pr(>|z|)    
#(Intercept)     -9.31088    2.23218  -4.171 3.03e-05 ***
#age_yr           2.34775    0.71609   3.279  0.00104 ** 
#I(age_yr^2)     -0.12902    0.05289  -2.440  0.01470 *  
#tb              -0.70186    0.59720  -1.175  0.23989    
#brucpositive    -1.37787    0.57802  -2.384  0.01714 *  
#tb:brucpositive  1.87521    0.92290   2.032  0.04217 *  

data4$age_yr2<-NA
for (i in 1:length(data4$age_yr)){
data4$age_yr2[i]<- (data4$age_yr[i]- mean(data4$age_yr))/sd(data4$age_yr)}

final.mod<-glmer(data4$fec ~ age_yr2+ I(age_yr2^2) + tb *bruc+ (1|id), data=data4, family=binomial(link="logit")); summary(final.mod)
Fixed effects:
                Estimate Std. Error z value Pr(>|z|)    
(Intercept)      -0.9990     0.3532  -2.828  0.00468 ** 
age_yr2           2.2636     0.4783   4.733 2.21e-06 ***
I(age_yr2^2)     -0.5408     0.2217  -2.439  0.01472 *  
tb               -0.7018     0.5972  -1.175  0.23991    
brucpositive     -1.3779     0.5780  -2.384  0.01714 *  
tb:brucpositive   1.8752     0.9229   2.032  0.04217 * 

d<- data4[data4$age_yr >= 3 & data4$age_yr <= 9,]
d$age_yr3 <- d$age_yr - 3
final.mod<-glmer(fec ~ age_yr2+ I(age_yr2^2) + tb *bruc+ (1|id), data=d, family=binomial(link="logit")); summary(final.mod)
final.mod<-glmer(fec ~ age_yr2+ I(age_yr2^2) + tb +bruc+ (1|id), data=d, family=binomial(link="logit")); summary(final.mod)
for (i in 1:length(data4$age_yr)){
data4$age_yr2[i]<- (data4$age_yr[i]- mean(data4$age_yr))/sd(data4$age_yr)}
final.mod<-glmer(fec ~ age_yr2+ I(age_yr2^2) + tb*bruc+ (1|id), data=d, family=binomial(link="logit")); summary(final.mod)

d<- data4

d$tb[d$tb==0] <- "negative"
d$tb[d$tb==1] <- "positive"
d$tb <- as.factor(d$tb)
ages<- c(3,4,5,6,7,8)
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




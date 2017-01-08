#########################################
#########################################
# Incidence analyses
#########################################
#########################################
library(ggplots2)
library(lme4)
library(pbkrtest)
library("JM")
library("lattice")
library("survival")

############################################################
############################################################
data<-read.csv("~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/final_datasets_copied_from_phdfolder/cross_sectional_data_withdz_cleandisease_nofinal_Feb2016_capturetime_forsurv.csv")  # 5 animals added from last time

# make datasets for Courtney
convert<-data[data $brucconvert==1 & data$convert==1,]
write.csv(convert, "~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/final_datasets_copied_from_phdfolder/convertersonly_Feb2016.csv")
length(convert$id); length(unique(convert$id))  # 14 both, 98 time observations

brconvert<- data[data$brucconvert==1,]
write.csv(brconvert, "~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/final_datasets_copied_from_phdfolder/brucconverters_Feb2016.csv")
length(brconvert$id); length(unique(brconvert$id))  # 29 became infected with brucellosis, 190 time observations

tbconvert<- data[data$convert==1,]
write.csv(tbconvert, "~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/final_datasets_copied_from_phdfolder/tbconverters_Feb2016.csv")
length(tbconvert$id); length(unique(tbconvert$id))  # 47 became infected with bTB; 288 time observations


############################################################
############################################################
# Summary statistics: 
############################################################
############################################################
# age of first infection
##TB
quantile(data$age_sel[data$incid==1], c(0.05, 0.95))/12

## Bruc
incidbr<- data[data$incid_bruc_jo ==1,c(2:5,9, 15, 30)]
quantile(incidbr$age_sel, c(0.05, 0.95))/12
summary(incidbr)

# who first
table(incidbr$tb)  # brucellosis convertes with and without bTB first
table(data$bruc[data$incid==1])  # bTB converters with and without brucellosis


# month and year histograms for suppliment
incidbr<- data[data$incid_bruc_jo ==1,c(2:5,9)]
incidtb<- data[data$incid==1, c(2:5,9)]
incidbr$month<- NA; incidtb$month<- NA
get_month = function(c){
	val <- as.character(c)
	t<- strsplit(val, '-')[[1]][2]
	ifelse (as.numeric((strsplit(t, "")[[1]][1])) > 0, 
		new <- paste(strsplit(t, "")[[1]][1], strsplit(t, "")[[1]][2], sep=""), 
		new <- strsplit(t, "")[[1]][2] )
	r <- as.numeric(new)
	return(r)
}
for (i in 1:length(incidbr$capid)){
	incidbr$month[i] <- get_month(incidbr$capid[i])
}
for (i in 1:length(incidtb$capid)){
	incidtb$month[i] <- get_month(incidtb$capid[i])
}
table(incidtb$month, incidtb$herdorig)
table(incidbr$month, incidbr$herdorig)

# initial prevalence, bTB 
table(data$tb[data$capturetime==0])
table(data$tb[data$capturetime==3])

# initial prevalence, brucellosis
table(data$bruc[data$capturetime==0])
table(data$bruc[data$capturetime==3])


############################################################
############################################################
# Models- BRUCELLOSIS INCIDENCE
############################################################
############################################################
data<- read.csv("~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/final_datasets_copied_from_phdfolder/brucellosis_incidence.csv")
data$animal<-as.character(data$animal)
data$age_yr2<- floor(data$age_yr)

# just tb
test.mod<-coxph(Surv(start, stop, convert.time)~ TB_3+ cluster(animal), data=data)  # 0.11
test.mod<-coxph(Surv(start, stop, convert.time)~ TB_3 +herd2+ cluster(animal), data=data); summary(test.mod)
test.mod<-coxph(Surv(start, stop, convert.time)~ TB_3 +herd2+age_yr2+ cluster(animal), data=data); summary(test.mod)
test.mod<-coxph(Surv(start, stop, convert.time)~ TB_3 +herd2+age_yr2+ I(age_yr2^2)+ cluster(animal), data=data); summary(test.mod) # no age^2
test.mod<-coxph(Surv(start, stop, convert.time)~ TB_3*herd2+age_yr2+ cluster(animal), data=data); summary(test.mod)
test.mod<-coxph(Surv(start, stop, convert.time)~ TB_3*age_yr2+herd2+ I(age_yr2^2)+ cluster(animal), data=data); summary(test.mod)
test.mod<-coxph(Surv(start, stop, convert.time)~ TB_3*I(age_yr2^2)+age_yr2+herd2+  cluster(animal), data=data); summary(test.mod) # no

test.mod<-coxph(Surv(start, stop, convert.time)~ TB_3*herd2+age_yr2+ cluster(animal), data=data); summary(test.mod)
test.mod<-coxph(Surv(start, stop, convert.time)~ TB_3*age_yr2+herd2+ cluster(animal), data=data); summary(test.mod)
test.mod<-coxph(Surv(start, stop, convert.time)~ TB_3*age_yr2+TB_3*herd2+ cluster(animal), data=data); summary(test.mod)


#Selection (age, continuous)
test.mod<-coxph(Surv(start, stop, convert.time)~ TB_3*herd2+age_yr2+ cluster(animal), data=data); extractAIC(test.mod) # 244.6, 
test.mod<-coxph(Surv(start, stop, convert.time)~ TB_3*age_yr2+herd2+ cluster(animal), data=data); extractAIC(test.mod)
temp<- data[data$convert.time==1,]

table(temp$age_yr2, temp$herd2)
table(temp$age_yr2, temp$TB_3)

par(mfrow=c(1,2))
hist(data$age_yr2, xlab="Age of first capture", col="lightgray")
hist(temp$age_yr2, xlab="Age of first capture of converters", col="lightgray")


data2<- data[data$age_yr2 <6,]
length(data$convert.time[data$convert.time==1])  # 30 converters. 
length(data2$convert.time[data2$convert.time==1])  # 30 converters. 

# TB*age effect remains significant even if subset by buffalo < 5yrs
test.mod<-coxph(Surv(start, stop, convert.time)~ TB_3+herd2+age_yr2+ cluster(animal), data=data2); summary(test.mod)
test.mod<-coxph(Surv(start, stop, convert.time)~ TB_3*herd2+age_yr2+ cluster(animal), data=data2); summary(test.mod)
test.mod<-coxph(Surv(start, stop, convert.time)~ TB_3*age_yr2+herd2+ cluster(animal), data=data2); summary(test.mod)

data2<- data[data$age_yr2 <14,] # yes, still holds


# And if use categories...above (REPORTED WITH CONTINUOUS AGE...)
data$testage<- NA
data$testage[data$age_yr <3] <- "young"
data$testage[data$age_yr >= 3] <- "old"
test.mod<-coxph(Surv(start, stop, convert.time)~ TB_3*herd2+testage+ cluster(animal), data=data); summary(test.mod) # 244.4092
test.mod<-coxph(Surv(start, stop, convert.time)~ TB_3*herd2+TB_3*testage+ cluster(animal), data=data); summary(test.mod) #2 46.9492
test.mod<-coxph(Surv(start, stop, convert.time)~ herd2+TB_3*testage+ cluster(animal), data=data); summary(test.mod) # 248.7528
test.mod<-coxph(Surv(start, stop, convert.time)~ TB_3+herd2+testage+ cluster(animal), data=data); summary(test.mod) # 248.75

data$age4<- NA
data$age4[data$age_yr <4] <- "young"
data$age4[data$age_yr >= 4] <- "old"
data$age5<- NA
data$age5[data$age_yr < 5] <- "young"
data$age5[data$age_yr >= 5] <- "old"
data$age2 <- NA
data$age2[data$age_yr < 2] <- "young"
data$age2[data$age_yr >= 2] <- "old"
data$age3p <- NA
data$age3p[data$age_yr < 3] <- "young"
data$age3p[data$age_yr >= 3 & data$age_yr <5] <- "adult"
data$age3p[data$age_yr >=5] <- "old"
data$age2p <- NA
data$age2p[data$age_yr < 2] <- "young"
data$age2p[data$age_yr >= 2 & data$age_yr <4] <- "adult"
data$age2p[data$age_yr >=4] <- "old"


data$agelate <- NA
data$agelate[data$age_yr < 2.5] <- "young"
data$agelate[data$age_yr < 5.5 & data$age_yr >= 2.5] <- "subadult"
data$agelate[data$age_yr >= 5.5] <- "mature"

data$ageearly <- NA
data$ageearly[data$age_yr <1.5] <- "young"
data$ageearly[data$age_yr <5.5 & data$age_yr >=1.5] <- "subadult"
data$ageearly[data$age_yr >= 5.5] <- "mature"

test.mod<-coxph(Surv(start, stop, convert.time)~ TB_3+age4+TB_3*herd2+ cluster(animal), data=data); summary(test.mod) # 247.3
test.mod<-coxph(Surv(start, stop, convert.time)~ TB_3+age5+TB_3*herd2+ cluster(animal), data=data); summary(test.mod)  # 248.9
test.mod<-coxph(Surv(start, stop, convert.time)~ TB_3*testage+herd2+ cluster(animal), data=data); summary(test.mod) #244.4

# age category alone!
test.mod<-coxph(Surv(start, stop, convert.time)~ testage + cluster(animal), data=data); summary(test.mod) # 248.0299
test.mod<-coxph(Surv(start, stop, convert.time)~ age3p + cluster(animal), data=data); summary(test.mod) # 249.862
test.mod<-coxph(Surv(start, stop, convert.time)~ age2 + cluster(animal), data=data); summary(test.mod) # 251.63
test.mod<-coxph(Surv(start, stop, convert.time)~ age2p + cluster(animal), data=data); summary(test.mod) # 251.53
test.mod<-coxph(Surv(start, stop, convert.time)~ age4 + cluster(animal), data=data); summary(test.mod) # 249.5333
test.mod<-coxph(Surv(start, stop, convert.time)~ age5 + cluster(animal), data=data); summary(test.mod) # 250.5

# Effect size for bTB: 
1) In LS, continuous age, 3.9
2) In CB, continuous age, 0.39 (TB*herd+age)
3) In LS, >3, 4.32
4) In CB, >3 0.39
5) overall, continuous age (TB + herd + age), 1.9624 (p = 0.123, age is significant and negative)
6) overall, categorical age, 2.129

# EVALUATION: Do we trust that age term- or is this a
test.mod<-coxph(Surv(start, stop, convert.time)~ TB_3*herd2+testage+ cluster(animal), data=data); summary(test.mod) # 244.4092
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
Shat2 <- survexp(~ TB_3, ratetable=test.mod, data=data)
with(Shat2, head(data.frame(time, surv), n=4))
#                 rho  chisq      p
#TB_3         -0.0289 0.0388 0.8439
#herd2CB      -0.3780 5.1152 0.0237
#testageyoung -0.1638 1.2866 0.2567
#GLOBAL            NA 5.1168 0.1634

test.mod<-coxph(Surv(start, stop, convert.time)~ TB_3*herd2+testage+ herd2:stop+ herd2:TB_3:stop + cluster(animal), data=data); summary(test.mod) # 244.4092
test.mod<-coxph(Surv(start, stop, convert.time)~ TB_3*herd2+testage+ herd2:stop+ cluster(animal), data=data); summary(test.mod) # 244.4092

test.mod<-coxph(Surv(start, stop, convert.time)~ TB_3*herd2+testage+ cluster(animal), data=data[data$start > 6,]); summary(test.mod) # 244.4092
test.mod<-coxph(Surv(start, stop, convert.time)~ TB_3*herd2+testage+ cluster(animal), data=data[data$start < 6,]); summary(test.mod) # 244.4092

# There is a postiive association with bTB in the later half of the study in Lower Sabie; earlier half of the study in CB. 
# Age effect is independent of all this cluster!
# Global constant hazard is ok; but some evidence that proporional hazards assumptions is violated for herd parameter... 
# Argue that model inference overall is valid... 
# coxph(formula = Surv(start, stop, convert.time) ~ TB_3 * herd2 + 
#    testage + herd2:stop + herd2:stop + cluster(animal), data = data)

#  n= 1049, number of events= 30 

#                 coef exp(coef) se(coef) robust se      z Pr(>|z|)   
#TB_3           1.3996    4.0538   0.5098    0.5430  2.578  0.00994 **
#herd2CB        2.7468   15.5934   1.3391    1.0903  2.519  0.01176 * 
#testageyoung   0.8466    2.3317   0.3891    0.4208  2.012  0.04422 * 
#TB_3:herd2CB  -2.2221    0.1084   1.1568    1.2207 -1.820  0.06870 . 
#herd2aLS:stop  0.2133    1.2378   0.1294    0.1117  1.910  0.05610 . 
#herd2CB:stop       NA        NA   0.0000    0.0000     NA       NA   

test.mod<-coxph(Surv(start, stop, convert.time)~ TB_3*herd2+testage+ herd2:stop + cluster(animal), data=data); summary(test.mod) # 244.4092

# Average effect of bTB, taken accross herds: 




############################################################
############################################################
# Models- BTB INCIDENCE
############################################################
############################################################
data<- read.csv("~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/final_datasets_copied_from_phdfolder/tb_incidence.csv")
data$animal<-as.character(data$animal)
data$age_yr2<- floor(data$age_yr)

length(data$convert.time[data$convert.time==1])  # 41 converters. (131 buffalo)

# just brucella additive models, age continuous (none significant)
test.mod<-coxph(Surv(start, stop, convert.time)~ brucella+ cluster(animal), data=data)  # 0.23
test.mod<-coxph(Surv(start, stop, convert.time)~ brucella +herd2+ cluster(animal), data=data); summary(test.mod)
test.mod<-coxph(Surv(start, stop, convert.time)~ brucella +herd2+age_yr2+ cluster(animal), data=data); summary(test.mod)
test.mod<-coxph(Surv(start, stop, convert.time)~ brucella +herd2+age_yr2+ I(age_yr2^2)+ cluster(animal), data=data); summary(test.mod) # no age^2

# 2 way
test.mod<-coxph(Surv(start, stop, convert.time)~ brucella*herd2+age_yr2+ cluster(animal), data=data); summary(test.mod)
test.mod<-coxph(Surv(start, stop, convert.time)~ brucella*age_yr2+herd2+ I(age_yr2^2)+ cluster(animal), data=data); summary(test.mod) # sig
test.mod<-coxph(Surv(start, stop, convert.time)~ brucella*I(age_yr2^2)+age_yr2+herd2+  cluster(animal), data=data); summary(test.mod) # sig
test.mod<-coxph(Surv(start, stop, convert.time)~ brucella*herd2+age_yr2+ cluster(animal), data=data); summary(test.mod)
test.mod<-coxph(Surv(start, stop, convert.time)~ brucella*age_yr2+herd2+ cluster(animal), data=data); summary(test.mod) # sig
test.mod<-coxph(Surv(start, stop, convert.time)~ brucella*age_yr2+brucella*herd2+ cluster(animal), data=data); summary(test.mod)


#Selection (age, continuous)-  models not significant after excluding the one 14 year old buffalo
test.mod<-coxph(Surv(start, stop, convert.time)~ brucella*I(age_yr2^2)+age_yr2+herd2+  cluster(animal), data=data); extractAIC(test.mod) # 350.76
test.mod<-coxph(Surv(start, stop, convert.time)~ brucella*age_yr2+herd2+ I(age_yr2^2)+ cluster(animal), data=data); extractAIC(test.mod) # 352.12
test.mod<-coxph(Surv(start, stop, convert.time)~ brucella*age_yr2+herd2+ cluster(animal), data=data); extractAIC(test.mod)  # 350.18
temp<- data[data$convert.time==1,]

# set up age categories
par(mfrow=c(1,2))
hist(data$age_yr2, xlab="Age of first capture", col="lightgray")
temp<- data[data$convert.time==1,]
hist(temp$age_yr2, xlab="Age of first capture of converters", col="lightgray")

data$testage<- NA
data$testage[data$age_yr <3] <- "young"
data$testage[data$age_yr >= 3] <- "old"


data$age4<- NA
data$age4[data$age_yr <4] <- "young"
data$age4[data$age_yr >= 4] <- "old"
data$age5<- NA
data$age5[data$age_yr < 5] <- "young"
data$age5[data$age_yr >= 5] <- "old"

# just brucella additive models, age categorical (none significant)
test.mod<-coxph(Surv(start, stop, convert.time)~ brucella +herd2+testage+ cluster(animal), data=data); summary(test.mod) # 1.29 times higher... but n.s.
test.mod<-coxph(Surv(start, stop, convert.time)~ brucella +herd2+age4+ cluster(animal), data=data); summary(test.mod)
test.mod<-coxph(Surv(start, stop, convert.time)~ brucella +herd2+age5+ cluster(animal), data=data); summary(test.mod)

# two-way interactions
test.mod<-coxph(Surv(start, stop, convert.time)~ brucella*testage +herd2+ cluster(animal), data=data); summary(test.mod) #n.s.
test.mod<-coxph(Surv(start, stop, convert.time)~ brucella* age4 +herd2+ cluster(animal), data=data); summary(test.mod)
test.mod<-coxph(Surv(start, stop, convert.time)~ brucella* age5 +herd2+ cluster(animal), data=data); summary(test.mod)

# Selection (age, categorical)
test.mod<-coxph(Surv(start, stop, convert.time)~ herd2 + testage+ cluster(animal), data=data); extractAIC(test.mod) 	# 347.86
test.mod<-coxph(Surv(start, stop, convert.time)~ herd2 + age4+ cluster(animal), data=data); extractAIC(test.mod) 		# 350.79
test.mod<-coxph(Surv(start, stop, convert.time)~ herd2 + age5+ cluster(animal), data=data); extractAIC(test.mod) 		# 349.6

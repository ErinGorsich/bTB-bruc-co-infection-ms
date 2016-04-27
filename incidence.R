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
quantile(data$age_sel[data$incid==1], c(0.5, 0.95))/12

brpos<- brconvert[brconvert$bruc_beforeafter==1,]
b<- NA
a<- tapply(brpos$age_sel, brpos$id, min)
for (i in 1:length(a)){
	if (!is.na(a[[i]])){
		b[i] <- a[[i]]
	}}

b[is.na(b)]<-0
c<- b[b>0]
quantile(c, c(0.5, 0.95))/12


# who first


# month and year histograms for suppliment





############################################################
# make incidence dataset
surv<-read.csv("~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/final_datasets_copied_from_phdfolder/survival/brucsurvival_TB3controls_longresidnomissing_noerrors_season2_final.csv")

# get IDs of brucellosis converters; tb converters
incidb<-data[data$bruc_beforeafter!="pfc",]
ids <- as.character(unique(incidb$id)); ids
convertids<- as.character(unique(brconvert$id))
incidtb<-data[data$tb_beforeafter!="pfc",]
idst <- as.character(unique(incidtb$id))

ids[ids %in% surv$animal]
ids[!(ids %in% surv$animal)]  
# 19 buffalo missing from survival dataset: "B22b" "B28"  "B37b" "B42b" "O35"  "R15"  "R34"  "R40"  "R7"   "R9"   "Y2"   "Y23b" "Y24"  "Y26b" "Y33"  "Y45"  "Y7"   "Y7b"  "Y8" 
# added all but R9 (not cleear dz status) and Y2 bTB status is fucked up
# check not bolus

incidb<-surv
incidb$convert.time<-NA

incidt<-surv
incidt$convert.time<-NA

# remove all but the smallest time associated with being br positive.
newdf<-NA
temp <- incidb[incidb$animal==ids[1],]
time <- min(temp$start[temp$brucella==1])
temp2 <- temp[temp$start <= time,]
temp2$convert.time[temp2$start==time]<-1
temp2$convert.time[temp2$start!=time]<-0
newdf<-temp2

for (i in 2:length(ids)){
	temp <- incidb[incidb$animal==ids[i],]
	if (length(temp$brucella[temp$brucella == 1])>0) {
		time <- min(temp$start[temp$brucella == 1])
		temp2 <- temp[temp$start <= time,]
		temp2$convert.time[temp2$start==time]<-1
		temp2$convert.time[temp2$start!=time]<-0
	newdf<- rbind(newdf, temp2)
	}
	else{
		temp2 <- temp
		temp2$convert.time<-0
		newdf<- rbind(newdf, temp2)
	}
}
write.csv(newdf, "~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/final_datasets_copied_from_phdfolder/brucellosis_incidence.csv")
# will have to add the remaining buffalo in by hand...
br<-read.csv("~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/final_datasets_copied_from_phdfolder/brucellosis_incidence.csv")
ids <- as.character(unique(br$animal))

incidb<-as.character(unique(data$id[data$bruc_beforeafter!="pfc"]))
toadd<- incidb[!(incidb %in% ids)]

dftoadd<- surv[surv$animal %in% toadd,]
write.csv(dftoadd, "test.csv")
# [1] "B25"  "B26b" "B28b" "B2b"  "B31"  "B32"  "B37"  "B41"  "B45"  "O1"   "O12"  "O14"  "O21"  "O23"  "O23b"
#[16] "O26"  "O27"  "O29"  "O29b" "O32"  "O32b" "O36"  "O47"  "O47b" "O51b" "O8"   "O8b"  "R13"  "R15b" "R15c"
#[31] "R19"  "R21"  "R28"  "R3"   "R32b" "R33"  "R34b" "R34c" "R44"  "R44b" "R46b" "R52"  "R7b"  "R9"   "W1b" 
#[46] "W6"   "Y10"  "Y16"  "Y16b" "Y17"  "Y2"   "Y21b" "Y24b" "Y26"  "Y2b"  "Y3"   "Y33b" "Y34"  "Y34b" "Y4"  
# [61] "Y42"  "Y46b" "Y47" 
# added to spreadsheet



############################################################
############################################################
# Models: 
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


#Selection
test.mod<-coxph(Surv(start, stop, convert.time)~ TB_3*herd2+age_yr2+ cluster(animal), data=data); extractAIC(test.mod)
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


######################################
######################################
# Prep cross sectional data
# this script differes from the previous in that I added different bTB statuses
######################################
######################################

# Capture data from All data in excel_June2013
cap<-read.csv("~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/final_datasets_copied_from_phdfolder/capturedf.csv")
bolus<-read.csv("~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/final_datasets_copied_from_phdfolder/bolusdf.csv")
d<-cap[,c(1,2,3,5,6, 17, 18, 20:22, 28, 29)]
colnames(d)<-c("capid", "id", "n", "yr", "date", "cond", "bolusgiven", "preg", "milk", "calf", "age_sel", "age_unsel")
d$shouldbebolus<-bolus$bolus.assignment_AnID.table[match(d$id, bolus$Animal.ID)]

########################################
# remove data taken from old Paul Cross study and 2 calves
idremove<-c("PC1", "PC2", "PC13", "Calf1", "Calf2", "S1", "D1")
d<-d[!(d$id %in% idremove),]

########################################
# remove bolus animals, and flag those with accidentally given bolus
# Summary: Removed all bolused animals; Y2, O47b, O10 all were given one bolus 
# but were primarily unbolused, these were labeled in the bolus_problem column.
d[d$bolusgiven=="TRUE" & d$shouldbebolus=="control",]   
# O47b-1109 and Y2-0609 were given a bolus but remainder of sequence control
temp<-d[d$bolusgiven=="FALSE" & d$shouldbebolus=="bolus",]
length(temp[,1])
# 132 not given bolus but labelled so!    

# remove those animal IDS who were bolus + AND always received their bolus
id<-unique(d$id[d$bolus=="TRUE" & d$shouldbebolus=="bolus"])  # 155
idproblem<-unique(d$id[d$bolusgiven=="FALSE" & d$shouldbebolus=="bolus"]) # 98, so can at least remove 58
length(unique(d$id[d$id %in% id & !(d$id %in% idproblem)]))
id2<-setdiff(id, idproblem)  # remove these, those that were always bolused AND should be bolused. 
d2<- d[!(d$id %in% id2),]

# many of the ones with one bolus are from cull, now remove those who should be bolus and only one false negative
d3<- d2[d2$id %in% id,]
table(d3$id, d3$bolusgiven)
# rm= those with one or two FALSE but still mostly bolused.
rm<- c("B12", "B15", "B17", "B18", "B23", "B24", "B27", "B3", "B34", "B35", "B36b", "B38", "B39c", "B4", "B40", "B40b", "B44", 
"B48b", "B50b", "B53", "B7", "B8", "B8b", "B9", "O11", "O15", "O16", "O17b", "O18", "O19", "O19b", "O20", "O22", 
"O24", "O25", "O28", "O28b", "O37", "O38", "O39c", "O3b", "O44", "O46", "O48", "O49", "O50", "O52", "O52c", "O6",
 "O7", "R1", "R10", "R14b", "R16", "R17", "R18", "R2", "R22", "R23b", "R24b", "R25b", "R26b", "R27b", "R29b", "R30", 
 "R36", "R37b", "R38", "R42b", "R45c", "R47", "R4b", "R4c", "R50b", "R51", "W1", "W2", "W4", "W5b", "Y11", "Y12", "Y14", "Y15c", 
 "Y18", "Y20b", "Y27", "Y28", "Y29", "Y31d", "Y35d", "Y37b", "Y39b", "Y43b", "Y44c", "Y5c", "Y6", "Y9")
length(d2[,1])
d3<-d2[!(d2$id %in% rm),]
#length(d3[,1])
#d3[d3$bolusgiven==TRUE,]   #O47b, O10, Y2

# Y2, O47b, O10 all were given one bolus but were primarily unbolused.  I keep these for now and add a column
d3$bolus_problem<-NA
prob<-c("Y2", "O47b", "O10")
for (i in 1:length(d3[,1])){
ifelse(d3$id[i] %in% prob, d3$bolus_problem[i]<-TRUE, d3$bolus_problem[i]<-FALSE)
}
data<-d3
#write.csv(d3, "~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/final_datasets_copied_from_phdfolder/cross_sectional_data_Feb2016.csv")



########################################
# add brucellosis data
# summary: will need to exclude 4 with NA status and final capture
# (which was smoothed over. set final=1)
data<-read.csv("~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/final_datasets_copied_from_phdfolder/cross_sectional_data_Feb2016.csv")
bruc<-read.csv("~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/final_datasets_copied_from_phdfolder/Brucellosis_ELISA_summary_Oct2012.csv")
data$SP<-bruc$S.P.[match(data$capid, bruc$Capt.ID)]
data$status<-bruc$status_sm[match(data$capid, bruc$Capt.ID)]
data$status_conservative<-NA
data$brucconvert<-NA
write.csv(data, "~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/final_datasets_copied_from_phdfolder/cross_sectional_data_Feb2016.csv")

####################
# In spreadsheet changes
####################
# removed redarts
# added values to columns for brucellosis= 
	#1) SP= raw test value (>159 is positive)
	#2) status_conservative: negative or positive. need to be >159 for two captures- based on diagnostic test paper.  Or >159 if only one capture before death/removal
	#3) convert= 0, 1 for whole animal ID.  Set value=1 if converted, 0 ow (postive from 1st captures and non-converters both 0, converters=1 for whole time period)
	# 4) check= needs follow up
	#5) Before vs. After: pfc, nc, 0, 1. To get incidence data! Subset by beforeafter = 0 or 1.  Before=0, after=1
	# 6) WHAT TO DO ABOUT THE CULL CAPTURE with missing data?  
	#Now both positives and negatives smoothed (5 Bs negative forward; 3 Os). Added indicatro variable for final capture. For example B41!!!!!!!! is really bad!

#sometimes had some positive, negative, removed? (O9)
# or postive then missing (O31, O10)
# missing start times: R20
# R43= positive, negative, restpositive
# what about the cull capture (Feb12 and July12)- 5 bs smoothed. 

# will have to exclude a few with NA values for brucellosis status 

########################################
# (1) add bTB data and immune data
setwd("~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/final_datasets_copied_from_phdfolder")
data<-read.csv("cross_sectional_data_Feb2016.csv")
bTB<-read.csv("bTB data for Erin.csv")
data$tborig<-bTB$TB_OPT.long[match(data$capid, bTB$capture.ID)]
data$tb<-bTB$TB.test..seropositivity..t.[match(data$capid, bTB$capture.ID)]
data$incid<-bTB$TB.Incidenc[match(data$capid, bTB$capture.ID)]
data$tb_beforeafter<-data$incid  # filled in as 0 before incidence, 1 after. pfc and nc
# assumed if positive, that final capture btb status same as previous.  
# Exclude these with by removing final capture. 
### testerror results were rechecked and added !!!  
###they largely match those in my cross-sectional patterns thesis chapter. 
### But are different than the final_data_with_dz_Feb2016 file. 

# immune data
data$ifngpokeconc<-bTB$IFNg.poke_ng.ml[match(data$capid, bTB$capture.ID)]
data$ifngplate<-bTB$IFNg.Plate.Number[match(data$capid, bTB$capture.ID)]
data$ifngdelay<-bTB$lag.time..days.[match(data$capid, bTB$capture.ID)]
data$bka_cont_ecoli<-bTB$BKA.Control.Ecoli[match(data$capid, bTB$capture.ID)]
data$bka_exp_ecoli<-bTB$BKA.Exp.Ecoli[match(data$capid, bTB$capture.ID)]
data$bka_killed_ecoli<-bTB$BKA.Ecoli.Control.Experimental[match(data$capid, bTB$capture.ID)]
data$hapto<-bTB$Haptoglobin.ng.ml[match(data$capid, bTB$capture.ID)]
data$hapto_plate<-bTB$Haptoglobin.Plate[match(data$capid, bTB$capture.ID)]

write.csv(data, "~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/final_datasets_copied_from_phdfolder/cross_sectional_data_withdz_Feb2016_mine.csv")


###################################################
# final groom for analysis
setwd("~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/final_datasets_copied_from_phdfolder")
data<-read.csv("cross_sectional_data_withdz_Feb2016_mine.csv")
cross <-data[data$missing_bruc_problem==0,]  # 
length(data[,1])-length(cross[,1]) # 11 excluded for missing brucellosis data at culls
length(unique(data$id))-length(unique(cross$id))  # but 0 animals
cross<-cross[cross$bruc!="",]  # removed 4 NA bruc results, as expected (92 total excluded)

# should be able to remove testerrors (for bTB)
table(cross$tb_beforeafter)   # 30 test errors to remove, all with NA bruc statuses
cross2<-cross[cross$tb_beforeafter!="testerror", ]    
length(cross[,1])-length(cross2[,1])     # 29 removed.
length(unique(cross$id))-length(unique(cross2$id))   # from 5 individuals with flip flop statuses. 

cross3<-cross2[!is.na(cross2$tb),]   # removes cull animals without data
length(cross2[,1])-length(cross3[,1])    # 7 removed- this is primarily removing the uncertain final captures. 
length(unique(cross2$id))-length(unique(cross3$id))  # from 0 animal
cross3$incid[is.na(cross3$incid)]<-0
cross3$convert[is.na(cross3$convert)]<-0


length(cross[,1])-length(cross3[,1]) # 36 (89 previously) sample times removed from bTB test errors
length(unique(cross$id))-length(unique(cross3$id))  #5 (12 previously) indiv lost from excluding tb errors
# missing conversions assumed to occur at the time of first positive test;
# if missing test result at begining before first positive capture, called pfc at first known test result. 

summary(as.data.frame(table(cross3$id)))
length(cross3[,1]); length(unique(cross3$id))   # 145 animals sampled a total of 834 times, with a median of 6 recaptures (range from 1-9 recaptures).

write.csv(cross3, "cross_sectional_data_withdz_cleandisease_withfinal_Feb2016.csv")

###################################################
# final groom for incidence analysis
data<-read.csv("cross_sectional_data_withdz_cleandisease_withfinal_Feb2016.csv")
brconverters<-data[!(data$bruc_beforeafter=="nc"),]
brconverters<-brconverters[!(brconverters$bruc_beforeafter=="pfc"),]
length(brconverters$id)  
length(unique(brconverters$id))  # 33 animals became seropositive for brucellosis; 226 time points 

tbconverters<-data[!(data$tb_beforeafter=="nc"),]
tbconverters<-tbconverters[!(tbconverters$tb_beforeafter=="pfc"),]
length(tbconverters$id)
length(unique(tbconverters$id))  # 44 animals became seropositive for bBT; 291 time points. 

incidtb<-tbconverters[tbconverters$incid=="1",]
summary(incidtb$age_sel)  # median=58.5= 4.87; mean=61.28= 5.1
quantile(incidtb$age_sel, c(0.05, 0.95)) 41.75 to 90.05
incidbr<-brconverters[brconverters$bruc_incid=="1",] # 3.4 to 7.5 yrs old by percentile
summary(incidbr$age_sel)  # med=58, mean=61.1
quantile(incidbr$age_sel, c(0.05, 0.95)) # 43.35= 3.6 to 86.55=7.2 yrs old by the 5th and 95th percentile. 

# season/month- when are they converting
incidbr$capid<-as.character(incidbr$capid)
for (i in 1:length(incidbr[,1])){
incidbr$month[i]<-substr(strsplit(incidbr$capid[i], "-")[[1]][2], 1, 2)
}
hist(as.numeric(incidbr$month), breaks=seq(1,12,1), xlab="Month", ylab="Number of brucellosis serovonversions", main="")

incidtb$capid<-as.character(incidtb$capid)
for (i in 1:length(incidtb[,1])){
  incidtb$month[i]<-substr(strsplit(incidtb$capid[i], "-")[[1]][2], 1, 2)
}
hist(as.numeric(incidtb$month), breaks=seq(1,12,1), xlab="Month", ylab="Number of tb serovonversions", main="")

# Make converters only dataset
bothconvert<-tbconverters[tbconverters$brucconvert=="1",]
length(bothconvert$id); length(unique(bothconvert$id))  # still 102, 13
write.csv(bothconvert, "convertersonly_Feb2016_mine.csv")

######################################
######################################
# Run GLMM
######################################
######################################
setwd("~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/final_datasets_copied_from_phdfolder/")
#data<-read.csv("cross_sectional_data_withdz_cleandisease_Feb2016.csv")
data<-read.csv("cross_sectional_data_withdz_cleandisease_withfinal_Feb2016.csv")
library(lme4)
data$age_sel<-as.numeric(data$age_sel)
data$age_sel_sq<-data$age_sel^2
data$age_sel_sd<-NA
data$age_yrsel<-data$age_sel/12; data_yrsel_sd<-NA
data$agecat<-NA
for (i in 1:length(data[,1])){
	data$age_sel_sd[i]<-(data$age_sel[i]-mean(data$age_sel))/sd(data$age_sel)
	data$age_sel_sq_sd[i]<-(data$age_sel_sq[i]-mean(data$age_sel_sq))/sd(data$age_sel_sq)
	data$age_yrsel_sd[i]<- (data$age_yrsel[i]- mean(data$age_yrsel))/sd(data$age_yrsel)
	if(data$age_yrsel[i]<4){
		data$agecat[i]<-"aayoungadult"}
		else if(data$age_yrsel[i]>=4 & data$age_yrsel[i]<7){
			data$agecat[i]<-"primeadult"}
		else if (data$age_yrsel[i]>=7 & data$age_yrsel[i]<10){
			data$agecat[i]<-"adult" }
			else if (data$age_yrsel[i]>=10) {data$agecat[i]<-"old"}
	}


full.mod<-glmer(bruc~age_sel_sd+ herdorig+ tb+ (1|id), family=binomial(link="logit"), data=data)
data$yr<-as.factor(data$yr)
full.mod<-glmer(bruc~age_sel_sd+ age_sel_sq_sd+ herdorig+ tb+ (1|yr/id), family=binomial(link="logit"), data=data)
red.mod<-glmer(bruc~age_sel_sd+ age_sel_sq_sd+ tb+ (1|yr/id), family=binomial(link="logit"), data=data)

table(data$bruc, data$tb)
# final capture included
188/(188+397)  # 32.1  # in bTB negatives
106/(106+143)  # 42.5  # in bTB positives


# with initial dataset
155/(351+155)  # prevalence in bTB negatives = 30%   # this is 50 less data points for bTB negatives. 
77/(77+121)  # prevalence in bTB positives= 39%
# final capture included
160/(358+160)  # 41.4
94/(133+94)  # 30.8


test<-read.csv("~/Documents/phd research/co-infection paper/cross-sectional patterns/Demogall_Jan2013dbupdate_annaages_TBaddedMay2013_noredartscull_suspos_control.csv")
red.mod<-glmer(Status_sm_conservative ~ TB+age_yrsel+ (1|Animal.ID), family=binomial(link="logit"), data=test)
# 33.7% vs. 46.8%
red.mod<-glmer(bruc~age_yrsel+ tb+ (1|id), family=binomial(link="logit"), data=data)
red.mod<-glmer(bruc~agecat+ tb+ (1|id), family=binomial(link="logit"), data=data)

data2<-data[data$final_capture =="0",]
red.mod<-glmer(bruc~agecat+ tb+ (1|id), family=binomial(link="logit"), data=data2)


table(data$bruc, data$tb, data$agecat)
# in <4: TB+= 6/19= 32%; TB-= 30/190= 15.8% 
# in 4-10: TB+= 95/(128+95)= 42.6; TB- = 135/(135+220)= 38%

table(data2$bruc, data2$tb, data2$agecat)


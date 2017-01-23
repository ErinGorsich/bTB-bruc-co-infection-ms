#############################################################
#############################################################
# Make dataset to look at number of new infections over time
#############################################################
#############################################################
data<-read.csv("~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/final_datasets_copied_from_phdfolder/cross_sectional_data_withdz_cleandisease_nofinal_Feb2016_capturetime_forsurv.csv")  # 5 animals added from last time

get_dz_status = function(data, capture, herd){
	# Input= data = data, capture= number specifiying capturetime, herd= "LS" or "CB"
	# for cLS$id, return a vector for each animal with... "NA", "S" "TB", "BR", "CO"
	
	if(herd == "LS"){
		id <- cLS$id
	}
	if(herd == "CB"){
		id <- cCB$id
	}
	dzdata <- data[data$capturetime %in% capture, c(2, 3, 8, 21, 30)]
	
	# fill status with disease status for each animal id
	status <- c(NA)
	for(i in 1:length(id)){
		tempdz <- dzdata[dzdata$id == id[i],]
		if(length(tempdz[,1]) > 1){
			print(paste("errror with ", id[i], "more than one animal with this id"))
			status[i] <- "error"
		} else if (length(tempdz[,1]) == 0){
			status[i] <- "NA"
		} else if(tempdz$bruc == "negative"& tempdz$tb == 0){
			status[i] <- "S"
		} else if(tempdz$bruc == "negative"& tempdz$tb == 1){
			status[i] <- "T"
		} else if(tempdz$bruc == "positive"& tempdz$tb == 0){
			status[i] <- "B"
		} else {
			status[i] <- "C"
		}
	}
	return(status)
}


################################################################
################################################################
# Lower Sabie dataset
################################################################
################################################################
cLS <- data.frame(id = unique(data$id[data$herdorig=="LS"]))
table(data$tb, data$bruc, data$capturetime, data$herdorig)

# fill in disease status from appropriate capture times	
################################################################
cLS$June08 <- get_dz_status(data = data, capture=c(0, 3), "LS")
cLS$Nov08 <- get_dz_status(data = data, capture=c(6), "LS")
cLS$June09 <- get_dz_status(data = data, capture=c(12), "LS")  # B16 = T on June09, 
cLS$Nov09 <- get_dz_status(data = data, capture=c(18), "LS")
cLS$June10 <- get_dz_status(data = data, capture=c(24), "LS")
cLS$Nov10 <- get_dz_status(data = data, capture=c(30), "LS")
cLS$June11 <- get_dz_status(data = data, capture=c(36), "LS")
cLS$Nov11 <- get_dz_status(data = data, capture=c(42), "LS")

# need to catch the few captured on off times and update
################################################################
data[data$herdorig=="LS" & data$capturetime== 9,]
# B16 always bruc negative, TB convert on June 2009 (error in capturetime, should be 12 not 9)
# O29b, started on April 29 (9); always bruc and bTB negative; 

data[data$herdorig=="LS" & data$capturetime== 15,]
# O32b, started Sept 09 (15); always bruc and bTB negative; 

data[data$herdorig=="LS" & data$capturetime== 33,]
#O1b- (capture time incorrect, should be 30, 36, 42 (T, C, C)
data[data$herdorig=="LS" & data$capturetime== 39,]
#B14 (capture time incorrect, 39 should be 42)
data[data$herdorig=="LS" & data$capturetime== 48,]
#B11- really captured august 12th, ok. 

cLS$June09[cLS$id =="B16"] <- "T"
cLS$Nov11[cLS$id =="O1b"] <- "C"
cLS$Nov11[cLS$id == "B14"] <- "B"

# check mortalities by double checking rows with missing values at the end 
################################################################
# ADD B DEATH EVENTS (ALL B DEATHS IN SURVIVAL SPREADSHEET MATCH MAIN DATA): 
cLS$Nov11[cLS$id =="B1"] <- "C"  # died on Feb 12 at end of this capture period
cLS$Nov11[cLS$id =="B10"] <- "C"  # died, retrieved on June 12
cLS[cLS$id =="B13", c(5:9)] <- "M"  # died, lions on July 09
#B14 #died, last capture Nov11 = end
cLS[cLS$id =="B2", c(5:9)] <- "M"  # died after 21 June 09
cLS[cLS$id =="B22", c(6:9)] <- "M"  # died, retreived, 21- June 2010; nov 09 last capture
cLS[cLS$id =="B26", c(8:9)] <- "M" # died, collar retrieved July-11B26
cLS[cLS$id =="B37", c(9)] <- "M"   # dead, Aug 11
cLS[cLS$id =="B42", c(3:9)] <- "M" # dead, unkn timeB42
# ADD O DEATH EVENTS 
#O1b  # died after Nov 2011
cLS[cLS$id =="O23", c(5:9)] <- "M"   # Died July2009
cLS[cLS$id =="O26", c(6:9)] <- "M"	# Died Last Nov 2009
cLS[cLS$id =="O29", c(3:9)] <- "M" 	# just after first cap (not clear, maybe capture induced)
cLS[cLS$id =="O32", c(3:9)] <- "M"	# died as above ??
cLS[cLS$id =="O32b", c(7:9)] <- "M"	# died, last July 2010
cLS[cLS$id =="O47", c(3:9)] <- "M"	# date not specified
cLS[cLS$id =="O51", c(5:9)] <- "M"	# Nov 09 is last
cLS[cLS$id =="O8", c(5:9)] <- "M"		# Dec 2009

# MISSING AND CAPTURE DEATHS (NAs)
# B21 went missing, Feb2011 last data point
# B22b missing, last data point
# B28- missing
# B37b went to Moz (Aug 2012)
# B42b (June 11 = last capture = last seen)
cLS[cLS$id =="O1", c(7:9)] <- "X"  # O1 capture death, ok not to exclude, Jan2011
# O1b died Feb2012
# 035 missing after Feb 2011
# O43 missing Nov 09 is last
#O10- went missing during conversion to bTB for 2009
write.csv(cLS, "~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/final_datasets_copied_from_phdfolder/LS_conversion_summary.csv")



################################################################
################################################################
# Crocodile Bridge dataset
################################################################
################################################################
cCB <- data.frame(id = unique(data$id[data$herdorig=="CB"]))
table(data$tb, data$bruc, data$capturetime, data$herdorig)

# fill in disease status from appropriate capture times	
################################################################
cCB$Oct08 <- get_dz_status(data = data, capture=c(0, 3), "CB")  # R46-1008 has error on capturetime so include t=0
cCB$April08 <- get_dz_status(data = data, capture=c(9), "CB")
cCB$Oct09 <- get_dz_status(data = data, capture=c(15), "CB")
cCB$April09 <- get_dz_status(data = data, capture=c(21), "CB")
cCB$Oct10 <- get_dz_status(data = data, capture=c(27), "CB")
cCB$April10 <- get_dz_status(data = data, capture=c(33), "CB")
cCB$Oct11 <- get_dz_status(data = data, capture=c(39), "CB")
cCB$April11 <- get_dz_status(data = data, capture=c(45), "CB")
# R20 captured at t=45 (0312) and t= 48 (0712), drop final capture

# only two buffalo specified above were captured on off times, no changes necessary 
################################################################
# strange errors. 
cCB$April11[cCB$id == "R11"] <- "B"
cCB$April11[cCB$id == "R13"] <- "B"
cCB$Oct10[cCB$id == "R15b"] <- "S"  # missing from cross_sectional datasheet, maybe because noted as redart in Nov
cCB$April10[cCB$id == "R15b"] <- "S" 
cCB$April11[cCB$id == "R43"] <- "B" # was missing last one
cCB$April11[cCB$id == "R7b"] <- "T" # was missing last one
cCB[cCB$id == "R9", c(5:9)] <- "M" # NOT SURE WHY NOT IN SURVIVAL DATASHEET, short time between capture and lions?
cCB[cCB$id == "Y17", c(9)] <- "S" 
cCB[cCB$id == "Y36", c(9)] <- "B" 
cCB[cCB$id == "Y40", c(9)] <- "C" 
cCB[cCB$id == "Y42", c(9)] <- "T"  # skipped this capture, but remained Br negative at culls, should have smoothed in datasheet 


# check mortalities by double checking rows with missing values at the end 
################################################################
#mortalities in survival datasheet (final_fixed)= 
Rmorts_survdatasheet <- data.frame(
	id = c("R15c", "R24", "R27", "R31", "R32", "R34b", "R34c", "R35", "R39", "R44", 
	"R45", "R45b", "R46", "R5", "R50", "R6", "R8", "W1", "Y1", "Y10", 
	"Y16", "Y20", "Y26", "Y30", "Y30b", "Y32", "Y34", "Y38", "Y39", "Y43", 
	"Y44"), 
	deathtime = c(27, 27, 15, 18, 15, 27, 45, 39, 30, 15, 
	9, 12, 0, 15, 15, 15, 3, 9, 15, 27,
	 27, 3, 15, 3, 33, 12, 21, 21, 21, 27, 
	 3))
	 
cCB[cCB$id == "R15c", c(7:9)] <- "M"	# darted on Oct10, found dead by pred Nov10
#cCB[cCB$id == "R24", c(7:9)] <- "M"	# bolus, not in datasheet, removed from survival later. 
#cCB[cCB$id == "R27", c(7:9)] <- "M"	# bolus, not in datasheet, removed from survival later- retrieved March2010, last picked up Nov09
cCB[cCB$id == "R31", c(5:9)] <- "M"		# died Nov09
cCB[cCB$id == "R32", c(4:9)] <- "M"		# R32- retrieved Oct2009 
cCB[cCB$id == "R34b", c(7:9)] <- "M"	# died March 2011, date set as Mar2010
# R34c- mortality after times here		# retrieved and last picked up on July 2012
#cCB[cCB$id == "R35", c(7:9)] <- "M"  	# BOLUS
cCB[cCB$id == "R39", c(7:9)] <- "M"	  	# last picked up Oct10, retrieved MArch 2011
cCB[cCB$id == "R44", c(5:9)] <- "M" 		# date set as Oct 09, picked up by Neels Dec 2009
#cCB[cCB$id == "R45", c(7:9)] <- "M"  	# BOLUS - no date
#cCB[cCB$id == "R45b", c(7:9)] <- "M"  	# BOLUS
cCB[cCB$id == "R46", c(3:9)] <- "M"		# picked up June2009, (darted Oct 08)
#cCB[cCB$id == "R5", c(7:9)] <- "M"		# BOLUS
cCB <- cCB[!cCB$id=="R5",] 				# shouldn't be included here, should remove from data. 
#cCB[cCB$id == "R50", c(7:9)] <- "M"	# BOLUS, not in datasheet, removed from survival later.
#cCB[cCB$id == "R6", c(7:9)] <- "M"		# BOLUS, not in datasheet, removed from survival later.
cCB[cCB$id == "R8", c(3:9)] <- "M"  	# ony one capture, unclear death date

# W & Ys
cCB[cCB$id == "W1", c(4:9)] <- "M"		# retreived 22Oct2009, last capture ?
cCB[cCB$id == "Y1", c(5:9)] <- "M"		# died 11Oct2009, lions (died close to 6)
cCB[cCB$id == "Y10", c(8:9)] <- "M"		# last observed 21-Oct-2010, picked up March 2011 (not sure about in between)
cCB[cCB$id == "Y16", c(9)] <- "M"		# last Oct2010, dead 11 March 2011 (so missing in between)
#cCB[cCB$id == "Y20", c(:9)] <- "M"		# BOLUS, not in datasheet, removed from survival later.
cCB[cCB$id == "Y26", c(5:9)] <- "M"		# died- Oct-Nov09
cCB[cCB$id == "Y30", c(3:9)] <- "M"		# unspecified death date, retrieved, May 2009
cCB[cCB$id == "Y30b", c(9)] <- "M"		# last picked up May 2011, retrieved June 2011 (not sure about in between)
cCB[cCB$id == "Y32", c(4:9)] <- "M"		# unspecified date, retrieved Sept 2009
cCB[cCB$id == "Y34", c(6:9)] <- "M"		# April 2010 (died close to 7 but set to 6)
cCB[cCB$id == "Y38", c(6:9)] <- "M"		# May 2010
#cCB[cCB$id == "Y39", c(:9)] <- "M"		# BOLUS, not in datasheet, removed from survival later.
#cCB[cCB$id == "Y43", c(:9)] <- "M"		# BOLUS
#cCB[cCB$id == "Y44", c(:9)] <- "M"		# BOLUS, not in datasheet, removed from survival later.


# MISSING AND CAPTURE DEATHS- KEEP NA
# R15- capture death- kept NA
# R21- capture death- kept NA
# R21b retrieved Oct 2009 - not in because unclear TB test result
# R40 missing- kept NA
# R7, unknown death time.  disease status for one capture (S) ... ?????? not sure if/when missing
cCB[cCB$id == "R7", c(3:9)] <- "M"
cCB[cCB$id == "R34", c(4:9)] <- "M"  # NOT SURE WHY NOT IN DEATH SPREADSHEET? unknown death time?, good bTB 
# R7b had its collar removed. 

#W3 went missing April 2012
# Y2 died 11-Oct-2009, captured on 24, Sept 2009, excluded b/c suspect test result (avian bTB and bovine bTB reactor...) at capture before death.
cCB[cCB$id == "Y2", c(5:9)] <- "M"	
# Y21 is missing	
# Y23b is missing
# Y23 died June 2011, good bTB status, not sure why not in capture datasheet? 
cCB[cCB$id == "Y23", c(7:9)] <- "M"
# Y23b = missing
# Y24 died, unknown date...bTB neg, NOT SURE WHY NOT IN SURVIVAL DATASET? maybe just unknown date? 
cCB[cCB$id == "Y24", c(3:9)] <- "M"

# Y26b = missing Nov 2011
# Y30c went missing March 2012
#Y33, death, unknown 
cCB[cCB$id == "Y33", c(3:9)] <- "M" #NOT SURE WHY NOT IN SURVIVAL DATASET? maybe just unknown date? 
# Y45 missing  March 2012
# Y46 collar removed April 2009
# Y7- unspecified death date- neg bTB, NOT SURE WHY NOT IN SURVIVAL DATASHEET. 
cCB[cCB$id == "Y7", c(3:9)] <- "M"
# Y7b missing Sept 2011
# Y8 missing

write.csv(cCB, "~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/final_datasets_copied_from_phdfolder/conversion_summary_CB.csv")

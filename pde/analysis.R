#############################################################
#############################################################
# Erin Gorsich
# This Code analyzes the pde co-infection model
# model is defined in rhs; parameters in fixed_parameters
#############################################################
#############################################################
#############################################################

#############################################################
#############################################################
#1) Load fixed parameters, model
#############################################################
#############################################################
rm(list = ls())
require("deSolve")
library("plyr")
library("ggplot2")
library("lattice") # for levelplots
library("gridExtra") # layout for lattice
library("RColorBrewer")
set.seed(5)

setwd("~/GitHub/bTB-bruc-co-infection-ms/pde")

# parameters
source('fixed_parameters.R', chdir = TRUE)

# rhs function - model structure
source('rhs.R', chdir = TRUE)

# age divisions in rhs function
agemax <- 20
agestep <- 0.1
N <- agemax / agestep
ages <- seq(1, agemax + 1, by = agestep)[-(N)]
N == length(ages)

# generate parameters with correct agebins
f.params <- gen_fixed_params(agemax, agestep, p = p, recovery = FALSE)
f.params.recov <- gen_fixed_params(agemax, agestep, p = p, recovery = TRUE)

# Functions for plotting (and define indecies based on ages, N): 
source('~/GitHub/bTB-bruc-co-infection-ms/pde/plotting_functions.R', chdir = TRUE)

# Starting agestructure (Jolles 2007; Caron et al. 2001)
juv <- rep(0.137 / length(ages[ages < 2]), length(ages[ages < 2]))
sa <- rep(0.368 / length(ages[ages >= 2 & ages < 6]), length(ages[ages >= 2 & ages < 6]))
a <- rep(0.185 / length(ages[ages >= 6 & ages < 9]), length(ages[ages >= 6 & ages < 9]))
ma <- rep(0.235 / length(ages[ages >= 9 & ages < 14]), length(ages[ages >= 9 & ages < 14]))
sen <- rep(0.075 / length(ages[ages >= 14 ]), length(ages[ages >= 14]))

relage <- c(juv, sa, a, ma, sen); length(relage) == N									
plot.relage <- c(0.137, rep(0.368/4, 4), rep(0.185/3, 3), rep(0.235/5, 5), rep(0.075/7, 7)) 	

#############################################################
#############################################################
#2) Analyses plotting Ro and EE using MC simulation from stats
#############################################################
#############################################################

# Ro calculations
#############################################################
# generate n random samples 
n = 1000
set.seed(1)
Ro_bTB_single_list2 <- NA; Ro_bTB_co_list2 <- NA
Ro_brucellosis_single_list2 <- NA; Ro_brucellosis_co2 <- NA

for(i in 1:n){
	params$rhoB<- exp(rnorm(n = 1, mean = 0.75579, sd = 0.40714))
	params$muB <- params$muS * exp(rnorm(n = 1, mean = 1.1060, sd = 0.3505))
	params$muT <- params$muS * exp(rnorm(n = 1, mean = 1.0370, sd = 0.3483 ))
	params$muC <- params$muS * exp(rnorm(n = 1, mean = 2.1430, sd = 0.5004329))
	params$muR <- params$muT
	params$muRC <- params$muC
	Ro_bTB_single_list2[i] <- Ro_bTB_single_age(params)
	Ro_bTB_co_list2[i] <- Ro_bTB_co_age(params)
	Ro_brucellosis_single_list2[i] <- Ro_brucellosis_single_age(params)
	Ro_brucellosis_co2[i] <- Ro_brucellosis_co_age(params)
}

quantile(Ro_bTB_single_list2, c(0.025, 0.05, 0.25, 0.5, 0.75, 0.90, 0.975))
quantile(Ro_bTB_co_list2, c(0.025, 0.05, 0.25, 0.5, 0.75, 0.90, 0.975))
quantile(Ro_brucellosis_single_list2, c(0.025, 0.05, 0.25, 0.5, 0.75, 0.90, 0.975))
quantile(Ro_brucellosis_co2, c(0.025, 0.05, 0.25, 0.5, 0.75, 0.90, 0.975))

data <- list(Ro_bTB_single_list2 = Ro_bTB_single_list2,  
	Ro_bTB_co_list2 = Ro_bTB_co_list2, 
	Ro_brucellosis_single_list2 = Ro_brucellosis_single_list2, 
	Ro_brucellosis_co2 = Ro_brucellosis_co2)

saveRDS(data, file = "Ro_confidence_interval_simulation_results.rds")


# For endemic equibrilium
#############################################################
# generate n random samples 
n = 1000
set.seed(1)
EE_bTB_single_list <- NA; EE_bTB_co_list <- NA
EE_brucellosis_single_list <- NA; EE_brucellosis_co <- NA
rhoB <- NA; mortB <- NA; mortT <- NA; mortC <- NA
pb <- txtProgressBar(min = 0, max = n, style = 3, char = "=")

for(i in 1:n){
	if(i %% 20 == 0){
		setTxtProgressBar(pb, i)
	}
	params$rhoB<- exp(rnorm(n = 1, mean = 0.75579, sd = 0.40714))
	rhoB[i] <- params$rhoB
	mortB[i] <- exp(rnorm(n = 1, mean = 1.1060, sd = 0.3505))
	mortT[i] <- exp(rnorm(n = 1, mean = 1.0370, sd = 0.3483 ))
	mortC[i] <- exp(rnorm(n = 1, mean = 2.1430, sd = 0.5004329))
	params$muB <- params$muS * mortB[i]
	params$muT <- params$muS * mortT[i]
	params$muC <- params$muS * mortC[i]
	params$muB[params$muB < 0] <- 0
	params$muT[params$muT < 0] <- 0
	params$muC[params$muC < 0] <- 0
	params$muR <- params$muT
	params$muRC <- params$muC
	val <- getEE(params)
	EE_bTB_single_list[i] <- val[1]
	EE_bTB_co_list[i] <- val[2]
	EE_brucellosis_single_list[i] <- val[3]
	EE_brucellosis_co[i] <- val[4]
}

quantile(EE_bTB_single_list, c(0.025, 0.05, 0.25, 0.5, 0.75, 0.90, 0.975))
quantile(EE_bTB_co_list, c(0.025, 0.05, 0.25, 0.5, 0.75, 0.90, 0.975))
quantile(EE_brucellosis_single_list, c(0.025, 0.05, 0.25, 0.5, 0.75, 0.90, 0.975))
quantile(EE_brucellosis_co, c(0.025, 0.05, 0.25, 0.5, 0.75, 0.90, 0.975))

data2 <- list(EE_bTB_single_list = EE_bTB_single_list,  
	EE_bTB_co_list = EE_bTB_co_list, 
	EE_brucellosis_single_list = EE_brucellosis_single_list, 
	EE_brucellosis_co = EE_brucellosis_co, 
	rhoB = rhoB, 
	mortB = mortB, 
	mortT = mortT,
	mortC = mortC)

saveRDS(data2, file = "EE_confidence_interval_simulation_results.rds")




#############################################################
#############################################################
#3) Analyses varying mortality and transmission rates
#############################################################
#############################################################

# Set up dataframe to hold the results
#############################################################
rhoB_test <- seq(0, 8, length.out = 101)
rhoT_test <- seq(0, 8, length.out = 101)
mort_test <- seq(0, 15, length.out = 101) 

# Data frame to hold results of changing bruc effects on bTB
epiTB <- data.frame(
	rhoT= rep(rhoB_test, length(rhoB_test)), 
	mort = rep(mort_test, each = length(rhoB_test)), 
	bTBprev = NA, brucprev = NA, 	finalN = NA, 
	bTB_inS = NA, bTB_inB = NA, bruc_inS = NA, bruc_inTB = NA)

# Data frame to hold results of changing bTB effects on bruc
epiB <- data.frame(
	rhoB= rep(rhoB_test, length(rhoB_test)), 
	mort = rep(mort_test, each = length(rhoB_test)),
	bTBprev = NA, brucprev = NA, finalN = NA, 
	bTB_inS = NA, bTB_inB = NA, bruc_inS = NA, bruc_inTB = NA)

# For BTB analyses set x0 as endemic prevalence of bruc in absence of BTB
#############################################################
S0 <- 400 * relage; It0 <- 0 * relage; Ib0 <- 10 * relage
Ic0 <- 0 * relage; R0 <- 10 * relage; Rc0 <- 0 * relage
times <- seq(1, 300, 1)
x0 <- c(S0, It0, Ib0, Ic0, R0, Rc0)

params <- c(f.params, list(gamma = 1/2, betaB = 0.6087396, 
		betaT = 0.0012974553, rhoT = 1, rhoB = 2.1))
sol <- as.data.frame(ode.1D(x0, times, rhs, params, 
		nspec = 6, dimens = N, method = "ode45"))
plot_raw_numbers(sol)
get_prevalence(sol)
endemic_agestructure <- unname(unlist( sol[500, c(2:length(colnames(sol)))] ))
xB <- endemic_agestructure
xB[(min(it.index)+50):(min(it.index)+53)] <- 1

imax <- c(length(epiTB[,1]))
pb <- txtProgressBar(min = 0, max = imax, style = 3)	
for (i in 1:length(epiTB[,1])){
	params.test <- c(f.params, list(gamma = 1/2, betaB = 0.6087396, 
		betaT = 0.0012974553, rhoT = epiTB$rhoT[i], rhoB = 2.1))
	params.test$muC <- epiTB$mort[i] * params.test$muS
	params.test$muC[params.test$muC > 1] <- 1
	params.test$muRC <- epiTB$mort[i] * params.test$muS
	params.test$muRC[params.test$muRC > 1] <- 1
		
	sol <- as.data.frame(ode(xB, times, rhs, params.test))
	temp <- get_prevalence(sol)
	
	epiTB$bTBprev[i] = temp$prevTB
	epiTB$brucprev[i] = temp$prevB 	
	epiTB$finalN[i] = sum(sol[length(sol), c(2:121)])
	epiTB$bTB_inS[i] = temp$prevTinS 
	epiTB$bTB_inB[i] = temp$prevTinB
	epiTB$bruc_inS[i] = temp$prevBinS 
	epiTB$bruc_inTB[i] = temp$prevBinT
	rm(params.test, sol, temp)
	setTxtProgressBar(pb, i)
}
cat("\n")
summary(epiTB)

write.csv(epiTB, "epiT.csv")


# For Brucellosis analyses set x0 as endemic prevalence BTB
# NOT UPDATED!
#############################################################

imax <- c(length(epiB[,1]))
pb <- txtProgressBar(min = 0, max = imax, style = 3)	
for (i in 1:length(epiB[,1])){
	params.test = c(fixed.params, list(gamma=1/2, theta = 4, K = 433, 
		betaB = 0.6087, betaT = 0.0012974553, rhoT = 1, rhoB = epiB$rhoB[i]))
	params.test$muC <- epiB$mort[i] * params.test$muS
	params.test$muC[params.test$muC > 1] <- 1
	params.test$muRC <- epiB$mort[i] * params.test$muS
	params.test$muRC[params.test$muRC > 1] <- 1
		
	sol <- as.data.frame(ode(x0, times, rhs_age_matrix, params.test))
	temp <- get_prevalence(sol)
	
	epiB$bTBprev[i] = temp$prevTB
	epiB$brucprev[i] = temp$prevB 	
	epiB$finalN[i] = sum(sol[length(sol[,1]), c(2:length(colnames(sol)))])
	epiB$bTB_inS[i] = temp$prevTinS 
	epiB$bTB_inB[i] = temp$prevTinB
	epiB$bruc_inS[i] = temp$prevBinS 
	epiB$bruc_inTB[i] = temp$prevBinT
	rm(params.test, sol, temp)
	setTxtProgressBar(pb, i)
}
cat("\n")

write.csv(epiB, "~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/draft2/post-labmeeting/epiB.csv")

################################################
# Again for epiB - but for when the rhoT = 1
rhoB_test <- seq(0, 8, length.out = 30)
rhoT_test <- seq(0, 8, length.out = 30)
mort_test <- seq(0, 7.2, length.out = 30) 

epib <- data.frame(
	rhoB= rep(rhoB_test, length(rhoB_test)), 
	mort = rep(mort_test, each = length(rhoB_test)),
	bTBprev = NA, brucprev = NA, 	finalN = NA, 
	bTB_inS = NA, bTB_inB = NA, bruc_inS = NA, bruc_inTB = NA)

imax <- c(length(epib[,1]))
pb <- txtProgressBar(min = 0, max = imax, style = 3)	
for (i in 1:length(epib[,1])){
	params.test_log = c(fixed.params.olddz, list(gamma=1/2, theta = 4, K = 433, 
		betaB = 1.025, betaT = 12.833531/10000, rhoT = 1, rhoB = epib$rhoB[i]))
	params.test_log$muC <- epib$mort[i] * params.test_log$muS
	#params.test_log$muC[params.test_log$muC > 1] <- 1
	params.test_log$muRC <- epib$mort[i] * params.test_log$muS
	#params.test_log$muRC[params.test_log$muRC > 1] <- 1
		
	sol <- as.data.frame(ode(x0, times, rhs_age_matrix, params.test_log))
	temp <- get_prevalence(sol)
	
	epib$bTBprev[i] = temp$prevTB
	epib$brucprev[i] = temp$prevB 	
	epib$finalN[i] = sum(sol[length(sol[,1]), c(2:length(colnames(sol)))])
	epib$bTB_inS[i] = temp$prevTinS 
	epib$bTB_inB[i] = temp$prevTinB
	epib$bruc_inS[i] = temp$prevBinS 
	epib$bruc_inTB[i] = temp$prevBinT
	rm(params.test, sol, temp)
	setTxtProgressBar(pb, i)
}
write.csv(epib, "~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/draft2/post-labmeeting/epib_rhoT1.csv")



# Brucellosis analyses
#############################################################


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
library("doParallel")
library("foreach")
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
source('plotting_functions.R', chdir = TRUE)

# Starting agestructure (Jolles 2007; Caron et al. 2001)
juv <- rep(0.137 / length(ages[ages < 2]), length(ages[ages < 2]))
sa <- rep(0.368 / length(ages[ages >= 2 & ages < 6]), length(ages[ages >= 2 & ages < 6]))
a <- rep(0.185 / length(ages[ages >= 6 & ages < 9]), length(ages[ages >= 6 & ages < 9]))
ma <- rep(0.235 / length(ages[ages >= 9 & ages < 14]), length(ages[ages >= 9 & ages < 14]))
sen <- rep(0.075 / length(ages[ages >= 14 ]), length(ages[ages >= 14]))

relage <- c(juv, sa, a, ma, sen); length(relage) == N									
plot.relage <- c(0.137, rep(0.368/4, 4), rep(0.185/3, 3), rep(0.235/5, 5), rep(0.075/7, 7)) 	
# Define x0, parameter for getEE function, disease free values
#############################################################
S0 <- 400 * relage; It0 <- 0 * relage; Ib0 <- 0 * relage
Ic0 <- 0 * relage; R0 <- 0 * relage; Rc0 <- 0 * relage
times <- seq(1, 300, 1)
x0 <- c(S0, It0, Ib0, Ic0, R0, Rc0)
params <- c(f.params, list(gamma = 1/2, betaB = 0, 
	betaT = 0, rhoT = 1, rhoB = 2.1))
test <- as.data.frame(ode.1D(x0, times, rhs, params, 
	nspec = 6, dimens = N, method = "ode45"))
x0 <- unname(unlist(test[length(test[,1]), c(2:(length(colnames(test))))])) 

#############################################################
#############################################################
#2) Analyses plotting Ro and EE using MC simulation from stats
#############################################################
#############################################################
# Functions to calculate EE and Ro
source('get_Ro.R', chdir = TRUE)
source('get_EE.R', chdir = TRUE)

# At ML parameters, find Ro and EE alone and with both pathogens
#############################################################
params <- c(f.params, list(gamma = 1/2, betaB = 0.5764065, 
		betaT = 1.3305462/1000, rhoT = 1, rhoB = 2.1))
#get_EE(params, x0, "beverton-holt")
# 65.855 TB alone
# 27.861 TB with co-infection
# 21.072 BRUC alone
# 31.493 BRUC with co-infection

# make sure get ok age-prev
xtest <- x0


# MCMC Ro calculations
# WARNING Ro reads in x0 frame global environment
#############################################################
n = 1000
set.seed(1)
binsize <- N/agemax

# set xB
xstart <- x0
xstart[min(ib.index)+1 + 3*binsize] <- 1
xstart[min(ib.index)+1 + 3*binsize] <- 1
params <- c(f.params, list(gamma = 1/2, betaB = 0.5764065, 
		betaT = 1.3305462/1000, rhoT = 1, rhoB = 2.1))
test <- as.data.frame(ode.1D(xstart, times, rhs, params, 
	nspec = 6, dimens = N, method = "ode45"))
xB <- unname(unlist(test[length(test[,1]), c(2:(length(colnames(test))))])) 
plot_raw_numbers(test)

# set xT
xstart <- x0
xstart[min(it.index)+1 + 3*binsize] <- 1
xstart[min(it.index)+1 + 4*binsize] <- 1
test <- as.data.frame(ode.1D(xstart, times, rhs, params, 
	nspec = 6, dimens = N, method = "ode45"))
xT <- unname(unlist(test[length(test[,1]), c(2:(length(colnames(test))))])) 
plot_raw_numbers(test)

Ro_bTB_single(params, x0)  # 4.025996
Ro_brucellosis_single(params, x0) #3.3862
Ro_bTB_co(params, xB)   # 1.596
Ro_brucellosis_co(params, xT) # 1.483981

# Generate 1000 samples
cl <- makeCluster(6)
registerDoParallel(cl)
d <- foreach(icount(n), .combine = rbind) %dopar% {
	params <- c(f.params, list(gamma = 1/2, betaB = 0.5764065, 
		betaT = 1.3305462/1000, rhoT = 1, rhoB = 2.1))
	params$rhoB<- exp(rnorm(n = 1, mean = 0.75579, sd = 0.40714))
	B <- rnorm(n = 1, mean = 1.1060, sd = 0.3505)
	TB <- rnorm(n = 1, mean = 1.0370, sd = 0.3483)
	params$muB <- params$muS * exp(B)
	params$muT <- params$muS * exp(TB)
	params$muC <- params$muS * exp(B + TB)
	params$muR <- params$muT
	params$muRC <- params$muC
	data <- data.frame(
		rhoB = params$rhoB, 
		dB = params$muB[1]/ params$muS[1],
		dT = params$muT[1]/ params$muS[1],
		Ro_bTB_single = Ro_bTB_single(params, x0),
		Ro_bTB_co = Ro_bTB_co(params, xB), 
		Ro_brucellosis_single = Ro_brucellosis_single(params, x0),
		Ro_brucellosis_co = Ro_brucellosis_co(params, xT))
}
stopCluster(cl)

quantile(d$Ro_bTB_single, c(0.025, 0.05, 0.25, 0.5, 0.75, 0.90, 0.975))
quantile(d$Ro_bTB_co, c(0.025, 0.05, 0.25, 0.5, 0.75, 0.90, 0.975))
quantile(d$Ro_brucellosis_single, c(0.025, 0.05, 0.25, 0.5, 0.75, 0.90, 0.975))
quantile(d$Ro_brucellosis_co, c(0.025, 0.05, 0.25, 0.5, 0.75, 0.90, 0.975))

saveRDS(d, file = "Ro_confidence_interval_simulation_results.rds")


# MCMC for endemic equibrilium
#############################################################
# generate n random samples 
n = 1000
set.seed(1)

# Generate 1000 samples
cl <- makeCluster(6)
registerDoParallel(cl)
d2 <- foreach(icount(n), .combine = rbind) %dopar% {
	params <- c(f.params, list(gamma = 1/2, betaB = 0.5764065, 
		betaT = 1.3305462/1000, rhoT = 1, rhoB = 2.1))
	params$rhoB<- exp(rnorm(n = 1, mean = 0.75579, sd = 0.40714))
	B <- rnorm(n = 1, mean = 1.1060, sd = 0.3505)
	TB <- rnorm(n = 1, mean = 1.0370, sd = 0.3483)
	params$muB <- params$muS * exp(B)
	params$muT <- params$muS * exp(TB)
	params$muC <- params$muS * exp(B + TB)
	params$muR <- params$muT
	params$muRC <- params$muC
	# think about runs with too high mortality!
	val <- get_EE(params, x0, "beverton-holt")
	data <- data.frame(
		rhoB = params$rhoB, 
		dB = params$muB[1]/ params$muS[1],
		dT = params$muT[1]/ params$muS[1],
		EE_bTB_single = val[1],  
		EE_bTB_co = val[2], 
		EE_brucellosis_single = val[3],
		EE_brucellosis_co = val[4])
}
stopCluster(cl)

saveRDS(d2, file = "EE_confidence_interval_simulation_results.rds")


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
# xB, endemic age structure for brucellosis only set above
xB.test <- xB[(min(it.index)+50):(min(it.index)+53)] <- 1

imax <- c(length(epiTB[,1]))
pb <- txtProgressBar(min = 0, max = imax, style = 3)	
for (i in 1:length(epiTB[,1])){
	params.test <- c(f.params, list(gamma = 1/2, betaB = 0.5764065, 
		betaT = 1.3305462/1000, rhoT = epiTB$rhoT[i], rhoB = 2.1))
	params.test$muC <- epiTB$mort[i] * params.test$muS
	params.test$muRC <- epiTB$mort[i] * params.test$muS
		
	sol <- as.data.frame(ode(xB.test, times, rhs, params.test))
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
xT.test <- xT[(min(ib.index)+50):(min(ib.index)+53)] <- 1

imax <- c(length(epiB[,1]))
pb <- txtProgressBar(min = 0, max = imax, style = 3)	
for (i in 1:length(epiB[,1])){
	params.test <- c(f.params, list(gamma = 1/2, betaB = 0.5764065, 
		betaT = 1.3305462/1000, rhoT = 1, rhoB = epiB$rhoB[i]))
	params.test$muC <- epiB$mort[i] * params.test$muS
	params.test$muRC <- epiB$mort[i] * params.test$muS
		
	sol <- as.data.frame(ode(xT.test, times, rhs, params.test))
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


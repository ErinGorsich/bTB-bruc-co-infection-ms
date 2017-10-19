#############################################################
#############################################################
# Erin Gorsich
# This Code runs and checks the age-structured co-infection model,
# it uses discrete aging, defined in rhs
#############################################################
#############################################################
#############################################################
# Outline:
# 1) Load fixed parameters, model
# 2) Set-up features of aging; Functions for plotting
# 3) Run model checks (no disease, single infections, together)
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


setwd("~/GitHub/bTB-bruc-co-infection-ms/discrete_aging_continuousdynamics")
# get fixed.params (assuming no recovery)
source('fixed_parameters.R', chdir = TRUE)
# rhs function, determinitic model
source('rhs_age.R', chdir = TRUE)
# function to run the model once, with aging
source('run_one.R', chdir = TRUE)

#############################################################
#############################################################
#2) Set-up features of aging, Functions for plotting
#############################################################
#############################################################
# age divisions in rhs function
s.index <- 1:20
it.index <- 21:40
ib.index <- 41:60
ic.index <- 61:80
r.index <- 81:100
rc.index <- 101:120

# initial age structure (Jolles 2007, Caron et al. 2001, max age at 20)
relage <- c(0.137, rep(0.368/4, 4), rep(0.185/4, 4),  	
	rep(0.235/6, 6), rep(0.075/5, 5))									

plot_agestructure = function(x){ # true based on Jolles et al. 2007
	# Input = row in the output of rhs, matrix form
	 xcounts <- NA
	 if (length(x) != 120){
	 	print("The age structure should include 20 ages,
	 		for 6 disease classes, giving 120 columns")
	 }
	 for(i in 1:20){
	 	xcounts[i] <- (x[i] + x[20+i] + x[40+i] + x[60+i] + x[80+i] + x[100+i])/sum(x)
	 }
	 d<-matrix(c(relage, xcounts), nrow=2, byrow=TRUE, 
	 	dimnames=list(c("Observed", "Predicted"), c(seq(1:20))) )
	print(d)
	barplot(d, beside = TRUE, col = c("light gray", "dark gray"),
		ylab = "Frequency (%)", ylim = c(0, max(d) + 0.01))
	box(	)
	legend("topright", legend = c("Observed", "Predicted"), fill = c("light gray", "dark gray"))
}

plot_raw_numbers = function(sol){
	plot(sol$time / 365, apply(sol[s.index+1], 1, sum), col= "black",
		type= 'l', ylim = c(0, 800), ylab = "Number of animals", 
		xlab = "Time (in years)")
	lines(sol$time / 365, apply(sol[it.index+1], 1, sum), col= "red")
	lines(sol$time / 365, apply(sol[ib.index+1], 1, sum), col= "blue")
	lines(sol$time / 365, apply(sol[ic.index+1], 1, sum), col= "green")
	lines(sol$time / 365, apply(sol[r.index+1], 1, sum), col = "orange")
	lines(sol$time / 365, apply(sol[rc.index+1], 1, sum), col = "pink")
	legend("topright", legend = c("S", "It", "Ib", "Ic", "R", "Rc"),
		col = c("black", "red", "blue", "green", "orange", "pink"), 
		bty="n", lty = 1)
}

get_prevalence = function(sol){
	S <- sum(sol[length(sol[,1]), s.index+1])
	It <- sum(sol[length(sol[,1]), it.index +1])
	Ib <- sum(sol[length(sol[,1]) , ib.index +1])
	Ic <- sum(sol[length(sol[,1]) , ic.index +1])
	R <- sum(sol[length(sol[,1]) , r.index +1])
	Rc <- sum(sol[length(sol[,1]) , rc.index +1])
	N <- sum(sol[length(sol[,1]), 2:121])
	prevTB <- (It + Ic + Rc) / N 
	prevB <- (Ib + Ic + R + Rc) / N
	prevBinS <- (Ib + R) / (S + Ib + R)
	prevBinT <- (Ic + Rc) / (It + Ic + Rc)
	prevTinS <- (It) / (S + It)
	prevTinB <- (Ic + Rc) / (Ib + Ic + R + Rc)
	return(list(prevTB = prevTB, prevB = prevB,
		prevBinS = prevBinS, prevBinT = prevBinT, 
		prevTinS = prevTinS, prevTinB = prevTinB))
}

#############################################################
#############################################################
# 3) Run model
#############################################################
#############################################################
# params and inits (no dz check, none takes off!)
S0 <- 400*relage; It0 <- 0*relage; Ib0 <- 0*relage; 
Ic0 <- 0*relage; R0 <- 0 * relage; Rc0 <- 0 * relage
x0 <- c(S0, It0, Ib0, Ic0, R0, Rc0)

params <- c(fixed.params, list(gamma=0.5/365.25, betaB = 0.6087396/365.25,
	betaT = 0.0012974553/365.25, rhoT = 1, rhoB = 2.1))

testout <- run_one(tmax = 500*365, x0 = x0, params = params)


par(mfrow = c(2,2))
#plot(x = testout$time, y = apply(testout[c(2:121)], 1, sum), 
#	pch = 19, main = "Annual aging")
plot_raw_numbers(testout)
plot_agestructure(as.matrix(testout[50,c(2:121)]))
plot_agestructure(as.matrix(testout[100,c(2:121)]))
plot_agestructure(as.matrix(testout[500,c(2:121)]))

# age structure stableizes by 50 (but not yet by 10) years
stable_age <- unname(unlist( testout[500, c(2:21)]/sum(testout[500, c(2:21)]) ))


# Test 2: Add brucellosis, only get brucellosis 
#############################################################
S0 <- 400* stable_age; It0 <- 0 * stable_age; Ib0 <- 20* stable_age; 
Ic0 <- 0* stable_age; R0 <- 30 * stable_age; Rc0 <- 0 * stable_age
x0 <- c(S0, It0, Ib0, Ic0, R0, Rc0)

betaB_temp <- seq(0.4, 1, 0.01) 
prevB <- NA; prevBrecov <- NA
for (i in 1:length(betaB_temp)){
	params.test = c(fixed.params, list(gamma=0.5/365.25, 
		betaB = betaB_temp[i]/365.25, betaT = 0.0012974553/365.25, rhoT = 1, rhoB = 2.1))
	sol <- run_one(tmax = 100*365, x0 = x0, params = params.test)
	prevB[i]<- get_prevalence(sol)$prevB
	get_prevalence(sol)$prevB
	rm(sol, params.test)
}

plot(x = betaB_temp, y = prevB, type = "l", xlab = expression(beta), 
	ylab = "Brucellosis prevalence", main = "Single disease")
abline(h = c(0.1, 0.2, 0.3, 0.4), col = "dark red")
# Beta value at:
betaB_temp[which(prevB < 0.06 & prevB > 0.04)]  		
betaB_temp[which(prevB < 0.11 & prevB > 0.09)]		
betaB_temp[which(prevB < 0.21 & prevB > 0.19)]		
betaB_temp[which(prevB < 0.31 & prevB > 0.29)]		
betaB_temp[which(prevB < 0.405 & prevB > 0.395)]	
# 30% prev alone at 0.58/365.25

# set endemic age-structure
params.test = c(fixed.params, list(gamma=0.5/365.25, 
	betaB = 0.6087396/365.25, betaT = 0.0012974553/365.25, rhoT = 1, rhoB = 2.1))
sol <- run_one(tmax = 500*365, x0 = x0, params = params.test)
endemic_agestructure <- unname(unlist( sol[500, c(2:121)] ))


# Test 3: Add bTB, only get bTB
#############################################################
S0 = 400*stable_age; It0 = 20 * stable_age; Ib0 = 0* stable_age; 
Ic0 = 0* stable_age; R0 = 0 * stable_age; Rc0 = 0 * stable_age
x0 = c(S0, It0, Ib0, Ic0, R0, Rc0)

betaT_temp <- seq(0.0001, 0.005, 0.0001) 
prevT <- NA; prevTrecov <- NA
for (i in 1:length(betaT_temp)){
	params.test = c(fixed.params, list(gamma=0.5/365.25, 
		betaB = 0.6087396/365.25, betaT = betaT_temp[i]/365.25, rhoT = 1, rhoB = 2.1))
	sol <- run_one(tmax = 100*365, x0 = x0, params = params.test)
	prevT[i]<- get_prevalence(sol)$prevTB
	rm(sol, params.test)
}

plot(x = betaT_temp, y = prevT, type = "l", xlab = expression(beta), 
	ylab = "BTB prevalence", main = "Single disease")
abline(h = c(0.1, 0.2, 0.3, 0.4), col = "dark red")
# Beta value at:
betaT_temp[which(prevT < 0.06 & prevT > 0.04)]  		
betaT_temp[which(prevT < 0.11 & prevT > 0.09)]		
betaT_temp[which(prevT < 0.21 & prevT > 0.19)]		
betaT_temp[which(prevT < 0.32 & prevT > 0.28)]		#0.0006
betaT_temp[which(prevT < 0.405 & prevT > 0.395)]


# Test 4: Add co-infection at set levels of brucellosis (vary change in susceptibility)
#############################################################
x0 = endemic_agestructure
x0[28] <- 5; x0[8] <- x0[8] - 5
params <- c(fixed.params, list(gamma=0.5/365.25, betaB = 0.6087396/365.25,
	betaT = 0.0012974553/365.25, rhoT = 1, rhoB = 2.1))
sol <- run_one(tmax = 200*365, x0 = x0, params = params.test)
plot_raw_numbers(sol)
get_prevalence(sol)


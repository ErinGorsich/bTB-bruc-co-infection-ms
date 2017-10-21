#############################################################
#############################################################
# Erin Gorsich
# This Code Reads in and runs the pde co-infection model
# model is defined in rhs; parameters in fixed_parameters
#############################################################
#############################################################
#############################################################
# Outline:
# 1) Load fixed parameters, model
# 2) Set-up features of aging; Functions for plotting
# 3) Test plots, with no Disease, none takes off (Generalized Beverton-Holt)
# 3) Test plots, BTB only 
# 4) Test plots Bruc only
# 5) Co-infection 
#############################################################
#############################################################
#############################################################
# This code owes much to:
# Age Structured Models by A. King & H. Wearing, &
# Structured models for host heterogenieties by J. Drake and P. Rohani
#  http://ms.mcmaster.ca/~bolker/eeid/2011_eco/waifw.pdf
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
library("pracma")
set.seed(5)
setwd("~/GitHub/bTB-bruc-co-infection-ms/pde")

# parameters
source('fixed_parameters.R', chdir = TRUE)

# rhs function - model structure
source('rhs.R', chdir = TRUE)

#############################################################
#############################################################
#2) Set-up features of aging, functions for plotting
#############################################################
#############################################################
agemax <- 20
agestep <- 0.1
N <- agemax / agestep
ages <- seq(1, agemax + 1, by = agestep)[-(N)]
N == length(ages)
s.index <- 1:N
it.index <- seq(N+1, 2*N)
ib.index <- seq((2*N+1), 3*N)
ic.index <- seq((3*N+1), 4*N)
r.index <- seq((4*N+1), 5*N)
rc.index <- seq((5*N+1), 6*N)

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
# 3) Test plots, with no Disease, none takes off
#############################################################
#############################################################
# Initial conditions
S0 <- 400 * relage; It0 <- 0 * relage; Ib0 <- 0 * relage
Ic0 <- 0 * relage; R0 <- 0 * relage; Rc0 <- 0 * relage
times <- seq(1, 500, 1)
x0 <- c(S0, It0, Ib0, Ic0, R0, Rc0)

params <- c(f.params, list(gamma = 1/2, betaB = 0.6087396, 
	betaT = 0.0012974553, rhoT = 1, rhoB = 2.1))
#params$muC[params$muC > 1] <- 0.9
#params$muRC[params$muRC > 1] <- 0.9

params.recov <- c(f.params.recov, list(gamma = 1/2, betaB = 0.6087396, 
	betaT = 0.0012974553, rhoT = 1, rhoB = 2.1))

system.time(
	test <- as.data.frame(ode(x0, times, rhs, params, method = "ode45"))  )
#test2 <- ode.1D(x0, times, rhs, params, nspec = 6, dimens = N, method = "ode45")
system.time(
	test3 <- as.data.frame(ode.1D(x0, times, rhs, params, 
		nspec = 6, dimens = N, method = "ode45")) )
system.time(
	test4 <- as.data.frame(ode.1D(x0, times, rhs, params,
		nspec = 6, dimens = N, method = "lsoda")) )

# plot age structure and disease-free dynamics
par(mfrow = c(2,2))
plot_agestructure(test, t = 100)
plot_agestructure(test, t = 500)
plot(x = test$time, y = apply(test[c(2:length(colnames(test)))], 1, sum),
	pch = 19, main = "sum")
plot_raw_numbers(test)

par(mfrow = c(2,2))
plot_agestructure(test3, t = 100)
plot_agestructure(test3, t = 500)
plot(x = test3$time, y = apply(test3[c(2:length(colnames(test3)))], 1, sum),
	pch = 19, main = "sum")
plot_raw_numbers(test3)

# set stable age structure in disease free context (1200 long N*6)
stable_age <- unname(unlist( test[length(test[,1]), c(2:(length(ages)+1))] / 
	sum(test[length(test[,1]), c(2:(length(ages)+1))]) ))

# choose which is faster test or test 3 and go...

#############################################################
#############################################################
# 4) Test plots, bTB only 
#############################################################
#############################################################
S0 <- 400 * stable_age; It0 <- 10 * stable_age; Ib0 <- 0 * stable_age
Ic0 <- 0 * stable_age; R0 <- 0 * stable_age; Rc0 <- 0 * stable_age
times <- seq(1, 500, 1)
x0 <- c(S0, It0, Ib0, Ic0, R0, Rc0)

betaT_temp <- seq(0.0001, 0.001, 0.00002) 
prevT <- NA; prevTrecov <- NA
prevB <- NA
for (i in 1:length(betaT_temp)){
	params.test = c(f.params, list(gamma=0.5, 
		betaB = 0.6087396, betaT = betaT_temp[i], rhoT = 1, rhoB = 2.1))
	sol <- as.data.frame(ode.1D(x0, times, rhs, params.test, 
		nspec = 6, dimens = N, method = "ode45"))
	prevT[i]<- get_prevalence(sol)$prevTB
	rm(sol, params.test)
}

plot(x = betaT_temp, y = prevT, type = "l", xlab = expression(beta), 
	ylab = "BTB prevalence", main = "Single disease")
abline(h = c(0.1, 0.2, 0.3, 0.4), col = "dark red")
# Beta value at:
betaT_temp[which(prevT < 0.12 & prevT > 0.08)]		
betaT_temp[which(prevT < 0.21 & prevT > 0.19)]		
betaT_temp[which(prevT < 0.32 & prevT > 0.28)]
betaT_temp[which(prevT < 0.405 & prevT > 0.395)]

# 30% TB around 0.0006

#############################################################
#############################################################
# 5) Test plots, Bruc only 
#############################################################
#############################################################
S0 <- 400* stable_age; It0 <- 0 * stable_age; Ib0 <- 20* stable_age; 
Ic0 <- 0* stable_age; R0 <- 30 * stable_age; Rc0 <- 0 * stable_age
x0 <- c(S0, It0, Ib0, Ic0, R0, Rc0)

betaB_temp <- seq(0.4, 1, 0.01) 
prevB <- NA; prevBrecov <- NA; prevT <- NA
for (i in 1:length(betaB_temp)){
	params.test = c(f.params, list(gamma=0.5, 
		betaB = betaB_temp[i], betaT = 0.00129, rhoT = 1, rhoB = 2.1))
	sol <- as.data.frame(ode.1D(x0, times, rhs, params.test, 
		nspec = 6, dimens = N, method = "ode45"))
	prevB[i]<- get_prevalence(sol)$prevB
	prevT[i]<- get_prevalence(sol)$prevTB
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


params.test = c(f.params, list(gamma=0.5, 
	betaB = 0.63, betaT = 0.00129, rhoT = 1, rhoB = 2.1))
sol <- as.data.frame(ode.1D(x0, times, rhs, params.test, 
		nspec = 6, dimens = N, method = "ode45"))
plot_raw_numbers(sol)

# Set endemic bruc age structure
endemic_agestructure <- unname(unlist( sol[500, c(2:length(colnames(sol)))] ))


#############################################################
#############################################################
# 6) Test plots, co-infection
#############################################################
#############################################################
x0 = endemic_agestructure
x0[min(it.index) + 100] <- 5
params <- c(fixed.params, list(gamma=0.5/365.25, betaB = 0.6087396/365.25,
	betaT = 0.0012974553/365.25, rhoT = 1, rhoB = 2.1))
sol <- as.data.frame(ode.1D(x0, times, rhs, params.test, 
		nspec = 6, dimens = N, method = "ode45"))
plot_raw_numbers(sol)
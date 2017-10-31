#############################################################
#############################################################
# Erin Gorsich
# This code fits the pde co-infection model to prevalence data
# model is defined in rhs; parameters in fixed_parameters
#############################################################
#############################################################
#############################################################
# Outline:
# 1) Load fixed parameters, model
# 2) Set-up features of aging
# 3) Define true prevalence and objective functions
# 4) Fit 
#############################################################
#############################################################
#############################################################


##################################################
##################################################
# 1) Load fixed parameters, model
###################################################
###################################################
# Accessory functions for fitting: 
rm(list = ls())
require("deSolve")
library("gridExtra")
library("ggplot2")
library("lattice") # for levelplots
setwd("~/GitHub/bTB-bruc-co-infection-ms/pde")

# parameters
source('fixed_parameters.R', chdir = TRUE)
# rhs function - model structure
source('rhs.R', chdir = TRUE)


#############################################################
#############################################################
#2) Set-up features of aging
#############################################################
#############################################################
agemax <- 20
agestep <- 0.1
N <- agemax / agestep
ages <- seq(1, agemax + 1, by = agestep)[-(N)]
binsize <- N / agemax
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

# Initial conditions
S0 <- 400 * relage; It0 <- 0 * relage; Ib0 <- 0 * relage
Ic0 <- 0 * relage; R0 <- 0 * relage; Rc0 <- 0 * relage
times <- seq(1, 500, 1)
x0 <- c(S0, It0, Ib0, Ic0, R0, Rc0)

# set stable age structure in disease free context (1200 long N*6)
params <- c(f.params, list(gamma = 1/2, betaB = 0.6087396, 
	betaT = 0.0012974553, rhoT = 1, rhoB = 2.1))
#params$muC > 1

params.recov <- c(f.params.recov, list(gamma = 1/2, betaB = 0.6087396, 
	betaT = 0.0012974553, rhoT = 1, rhoB = 2.1))
	
# use based on speed
test <- as.data.frame(ode.1D(x0, times, rhs, params, nspec = 6, dimens = N, method = "ode45"))
stable_age <- unname(unlist( test[length(test[,1]), c(2:(length(ages)+1))] / 
	sum(test[length(test[,1]), c(2:(length(ages)+1))]) ))

# no infections...
plot_raw_numbers(test)
get_prevalence(test)
	
##################################################
##################################################
# 3) Define true prevalence and objective functions
###################################################
###################################################

# Data: age structure and prevalence
#############################################################
prevTBobs <- 0.27  # for test- bootstrap estimate of overall prevalence
prevBobs <- 0.34
data <- read.csv("~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/final_datasets_copied_from_phdfolder/cross_sectional_data_withdz_cleandisease_nofinal_Feb2016_capturetime_forsurv.csv")
counts<- hist(data$age_sel/12, plot = FALSE)$counts  # youngest = 1.4 so aged 1-2
agestructure<- counts/sum(counts)
agestructure_yr = c(agestructure, 0, 0, 0, 0, 0)

# put in terms of aging features: 
data_agestructure <- agestructure_yr
#binsize <- 1 / agestep
#data_agestructure <- rep(agestructure_yr, each = binsize)

# Set starting prevalence for brucellosis in simulations
S0 <- 400 * stable_age; It0 <- 0 * stable_age; Ib0 <- 10 * stable_age
Ic0 <- 0 * stable_age; R0 <- 0 * stable_age; Rc0 <- 0 * stable_age
x0 <- c(S0, It0, Ib0, Ic0, R0, Rc0)
sol <- as.data.frame(ode.1D(x0, times, rhs, params, nspec = 6, dimens = N, method = "ode45"))

# Bruc prev to introduct BTB to in the optimizer
xB <- unname(unlist( sol[length(sol[,1]), c(2:length(colnames(sol)))] )) # 27% bruc prev

# functions used for fitting
#############################################################
objective = function(params.est){
	# params.est = 2 long = c(betaB, betaT)
	params <- c(f.params, list(gamma = 1/2, betaB = params.est[1], 
		betaT = params.est[2]/1000, rhoT = 1, rhoB = 2.1))
	
	# seed from endemic brucellosis conditions, 10 bTB positive buffalo
	x0 = xB
	x0[[min(it.index) + 5*binsize]] <- 5
	x0[[min(it.index) + 5*binsize + 1]] <- 5
	sol <- as.data.frame(ode.1D(x0, times, rhs, params, 
		nspec = 6, dimens = N, method = "ode45"))
	df <- get_structured_prevalence(sol)
	error <- sqrt(((prevTBobs - df$prevTB)^2 + (prevBobs - df$prevB)^2))
	return (error)
}

objective_recov = function(params.est){
	# params.est = 2 long = c(betaB, betaT)
	params <- c(f.params.recov, list(gamma = 1/2, betaB = params.est[1], 
		betaT = params.est[2]/1000, rhoT = 1, rhoB = 2.1))
	
	# seed from endemic brucellosis conditions, 10 bTB positive buffalo
	x0 = xB
	x0[[min(it.index) + 5*binsize]] <- 5
	x0[[min(it.index) + 5*binsize + 1]] <- 5

	sol <- as.data.frame(ode.1D(x0, times, rhs, params, 
		nspec = 6, dimens = N, method = "ode45"))
	df <- get_structured_prevalence(sol)
	error <- sqrt(((prevTBobs - df$prevTB)^2 + (prevBobs - df$prevB)^2))
	return (error)
}


#############################################################
#############################################################
#4) Fit
#############################################################
#############################################################
# No recovery (betaB = 0.5764065, betaT = 1.3305462)
par <- optim(c(1.025, 0.00094*1000), objective); par
par2 <- optim(c(1.5, 0.0007*1000), objective); par2
par3 <- optim(c(0.8, 0.001*1000), objective); par3

# Recovery (betaB = 0.6282482, betaT = 0.7848374)
par.r <- optim(c(1.025, 0.00094*1000), objective_recov); par.r


# Check optimization
#x <- seq(0.1, 1.5, 0.1)
#yval1 <- NA; yval2 <- NA
#for(i in 1:length(x)){
#	yval1[i] <- objective(c(x[i], 0.001382876))
#	yval2[i] <- objective(c(0.530219, x[i]))
#}

#par(mfrow = c(1, 2))
#plot(x = x, y = yval1, 
#	xlab = "betaB/365.25", ylab = "objective")
#plot(x = x, y = yval2, 
#	xlab = "betaT/(1000*365.25)", ylab = "objective")

# Plot best fit, no recovery
params <- c(f.params, list(gamma = 1/2, betaB = 0.5764065, 
		betaT = 1.3305462/1000, rhoT = 1, rhoB = 2.1))
x0 <- xB
x0[[min(it.index) + 5*binsize]] <- 5
x0[[min(it.index) + 5*binsize + 1]] <- 5
sol <- as.data.frame(ode.1D(x0, times, rhs, params, 
	nspec = 6, dimens = N, method = "ode45"))
get_structured_prevalence(sol)
plot_raw_numbers(sol)

# Plot best fit, recovery
params <- c(f.params.recov, list(gamma = 1/2, betaB = 0.6282482, 
		betaT = 0.7848374/1000, rhoT = 1, rhoB = 2.1))
x0 <- xB
x0[[min(it.index) + 5*binsize]] <- 1
x0[[min(it.index) + 5*binsize + 1]] <- 1
sol <- as.data.frame(ode.1D(x0, times, rhs, params, 
	nspec = 6, dimens = N, method = "ode45"))
get_structured_prevalence(sol)
plot_raw_numbers(sol)

#######################################
# TB Prev without brucellosis: 65.855
# TB prev in populations with co: 27.861
# BRUC prev without TB: 21.072
# BRUC prev in populations with co: 31.072
#######################################

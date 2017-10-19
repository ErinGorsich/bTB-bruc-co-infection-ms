#############################################################
#############################################################
# Erin Gorsich
# This Code estimates the transmission parameters for the age-structured
# co-infection model with discrete aging, defined in rhs
#############################################################
#############################################################
#############################################################
# Outline:
# 1) Load fixed parameters, model
# 2) Set-up features of aging; Functions for plotting
# 3) Define true prevalence and objective functions
# 4) Fit 
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
library("gridExtra")
library("ggplot2")

# get fixed.params (assuming no recovery)
source('~/GitHub/bTB-bruc-co-infection-ms/fixed_parameters.R', chdir = TRUE)
# rhs function, determinitic model
source('~/GitHub/bTB-bruc-co-infection-ms/rhs_age.R', chdir = TRUE)
# function to run the model once, with aging
source('~/GitHub/bTB-bruc-co-infection-ms/run_one.R', chdir = TRUE)

#############################################################
#############################################################
#2) Set-up features of aging and plotting functions
#############################################################
#############################################################
# set starting age structure as stable age distribution without disease
s.index <- 1:20
it.index <- 21:40
ib.index <- 41:60
ic.index <- 61:80
r.index <- 81:100
rc.index <- 101:120

relage <- c(0.137, rep(0.368/4, 4), rep(0.185/4, 4),  	
	rep(0.235/6, 6), rep(0.075/5, 5))			
S0 <- 400*relage; It0 <- 0*relage; Ib0 <- 0*relage; 
Ic0 <- 0*relage; R0 <- 0 * relage; Rc0 <- 0 * relage
x0 <- c(S0, It0, Ib0, Ic0, R0, Rc0)
params <- c(fixed.params, list(gamma=0.5/365.25, betaB = 0.6087396/365.25,
	betaT = 0.0012974553/365.25, rhoT = 1, rhoB = 2.1))	
sol <- run_one(tmax = 500*365, x0 = x0, params = params)
stable_age <- unname(unlist( sol[length(sol[,1]), c(2:21)]/sum(sol[length(sol[,1]), c(2:21)]) ))
	
# Functions for plotting
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
# 3) Define true prevalence and objective functions
#############################################################
#############################################################

# Data: age structure and prevalence
#############################################################
prevTBobs <- 0.27  # for test- bootstrap estimate of overall prevalence
prevBobs <- 0.34
data <- read.csv("~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/final_datasets_copied_from_phdfolder/cross_sectional_data_withdz_cleandisease_nofinal_Feb2016_capturetime_forsurv.csv")
counts<- hist(data$age_sel/12, plot = FALSE)$counts  # youngest = 1.4 so aged 1-2
agestructure<- counts/sum(counts)
data_agestructure = c(agestructure, 0, 0, 0, 0, 0)

# functions used for fitting
#############################################################
get_starting_eqbruc = function(params){
	S0 = 400* stable_age; It0 = 0 * stable_age; Ib0 = 20* stable_age; 
	Ic0 = 0* stable_age; R0 = 30 * stable_age; Rc0 = 0 * stable_age
	x0 = c(S0, It0, Ib0, Ic0, R0, Rc0)
	tmax <- 500*365
	sol <- run_one(tmax = tmax, x0 = x0, params = params)
	out <- unname(unlist( sol[length(sol[,1]), c(2:121)] ))
	return(out)
}

#make_summary_plots = function(sol){
#	df <- get_prevalence(sol)
#	df2 <- data.frame(Evaluation = c("Model", "Model", "Data", "Data"), 
#		Infection = c("Single", "Co", "Single", "Co"), 
#		BrucellosisPrevalence = c(df$prevBinS, df$prevBinT, 0.3035, 0.4524), 
#		TBPrevalence = c(df$prevTinS, df$prevTinB, 0.227, 0.3585))
#	p1 <- ggplot(df2, aes(x = Evaluation, y = BrucellosisPrevalence, fill = Infection)) + 
#		geom_bar(position = position_dodge(), stat = "identity") + ylim(0, 0.8)
#	p2 <- ggplot(df2, aes(x = Evaluation, y = TBPrevalence, fill = Infection)) + 
#		geom_bar(position = position_dodge(), stat = "identity") + ylim(0, 0.8)
#	grid.arrange(p1, p2, ncol = 2)
#}

make_structured_summary_plots = function(sol){
	df <- get_structured_prevalence(sol)
	df2 <- data.frame(Evaluation = c("Model", "Model", "Data", "Data"), 
		Infection = c("Single", "Co", "Single", "Co"), 
		BrucellosisPrevalence = c(df$prevBinS, df$prevBinT, 0.3035, 0.4524), 
		TBPrevalence = c(df$prevTinS, df$prevTinB, 0.227, 0.3585))
	p1 <- ggplot(df2, aes(x = Evaluation, y = BrucellosisPrevalence, fill = Infection)) + 
		geom_bar(position = position_dodge(), stat = "identity") + ylim(0, 0.8)
	p2 <- ggplot(df2, aes(x = Evaluation, y = TBPrevalence, fill = Infection)) + 
		geom_bar(position = position_dodge(), stat = "identity") + ylim(0, 0.8)
	grid.arrange(p1, p2, ncol = 2)
}

get_structured_prevalence = function(sol){
	S <-sum(sol[length(sol[,1]) , s.index+1] * data_agestructure)  # should give a scalar
	It <- sum(sol[length(sol[,1]) , it.index +1] * data_agestructure)
	Ib <- sum(sol[length(sol[,1]) , ib.index +1] * data_agestructure)
	Ic <- sum(sol[length(sol[,1]) , ic.index +1] * data_agestructure)
	R <- sum(sol[length(sol[,1]) , r.index +1] * data_agestructure)
	Rc <-sum(sol[length(sol[,1]) , rc.index +1] * data_agestructure)
	N <- sum(S + It + Ib + Ic + R + Rc)
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

# Objective function
#############################################################
objective = function(params.est){
	# params.est = 2 long = c(betaB, betaT)
	params <- c(fixed.params, list(gamma=0.5/365.25, betaB = params.est[1]/365.25,
		betaT = params.est[2]/(1000*365.25), rhoT = 1, rhoB = 2.1))
	
	# seed from endemic brucellosis conditions, 10 bTB positive buffalo
	x0 = get_starting_eqbruc(params = c(params))
	x0[[25]] <- x0[[25]] + 2
	if(x0[[5]] > 2){
		x0[[5]] <- x0[[5]] - 2} 
	tmax <- 1000*365
	sol <- run_one(tmax = tmax, x0 = x0, params = params)
 	df <- get_structured_prevalence(sol)
	error <- sqrt(((prevTBobs - df$prevTB)^2 + (prevBobs - df$prevB)^2))
	return (error)
}


#############################################################
#############################################################
#4) Fit
#############################################################
#############################################################
par <- optim(c(1.025, 0.00094*1000), objective); par 
par2 <- optim(c(1.5, 0.0007*1000), objective); par2
par3 <- optim(c(0.8, 0.001*1000), objective); par3
#0.530219, 1.382876

x <- seq(0.1, 1.5, 0.1)
yval1 <- NA; yval2 <- NA
for(i in 1:length(x)){
	yval1[i] <- objective(c(x[i], 0.001382876))
	yval2[i] <- objective(c(0.530219, x[i]))
}

par(mfrow = c(1, 2))
plot(x = x, y = yval1, 
	xlab = "betaB/365.25", ylab = "objective")
plot(x = x, y = yval2, 
	xlab = "betaT/(1000*365.25)", ylab = "objective")


# Add BTB to brucellosis system
#############################################################
params <- c(fixed.params, list(gamma=0.5/365.25, betaB = 0.530219/365.25,
	betaT = 0.001382876/365.25, rhoT = 1, rhoB = 2.1))	
x0 <- get_starting_eqbruc(params)  #19%prev alone
x0[[25]] <- x0[[25]] + 2
x0[[5]] <- x0[[5]] - 2

sol <- run_one(tmax = tmax, x0 = x0, params = params)
par(mfrow = c(1, 2))
plot_raw_numbers(sol)
plot_ageprevalence(sol)

# TB = 27.7 (23.7 in S vs. 36.7 in co)
# Bruc = 30.7 (26.8 in S vs. 40.7 in co)
get_prevalence(sol) 

# TB = 26.999 (23.69 vs 33.43)
# Bruc = 34 (31.006 vs. 42.0945)
get_structured_prevalence(sol)
make_structured_summary_plots(sol)
sum(sol[length(sol),]) # final population size

# Add brucellosis to BTB system
#############################################################
params <- c(fixed.params, list(gamma=0.5/365.25, betaB = 0.530219/365.25,
	betaT = 0.001382876/365.25, rhoT = 1, rhoB = 2.1))	
S0 = 400* stable_age; It0 = 20 * stable_age; Ib0 = 0* stable_age; 
Ic0 = 0* stable_age; R0 = 0 * stable_age; Rc0 = 0 * stable_age
x0 = c(S0, It0, Ib0, Ic0, R0, Rc0)
tmax <- 500*365
sol <- run_one(tmax = tmax, x0 = x0, params = params)  # 64.57% prevalence alone
x0 <- unname(unlist( sol[length(sol[,1]), c(2:121)] ))
get_prevalence(sol)

# add Bruc
x0[[45]] <- x0[[45]] + 2
x0[[5]] <- x0[[5]] - 2
sol <- run_one(tmax = tmax, x0 = x0, params = params)
plot_raw_numbers(sol)
get_prevalence(sol)  # 27.7 with bruc present


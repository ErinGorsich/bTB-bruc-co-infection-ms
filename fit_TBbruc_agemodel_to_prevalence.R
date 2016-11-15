#############################################################
#############################################################
#############################################################
#############################################################
# Code to fit age specific model to prevalence data 
# 1- Oct- 2016
#############################################################
#############################################################
#############################################################
#############################################################
# Outline
#############################################################
# 1) Load fixed parameters and model; diffeqs; stochastic version
# 2) Test Plot, grooming functions
# 3) Fit Deterministic Version to overall prevalence
#############################################################

#############################################################
#############################################################
#1) Load fixed parameters, model
#############################################################
#############################################################
require("deSolve")
library("plyr")
library("ggplot2")
set.seed(5)
# get fixed.params & fixed.params.recov
source('~/GitHub/bTB-bruc-co-infection-ms/fixed_parameters_norecovery_age.R', chdir = TRUE)
source('~/GitHub/bTB-bruc-co-infection-ms/fixed_parameters_recovery_age.R', chdir = TRUE)
# rhs function, determinitic model, age structure
source('~/GitHub/bTB-bruc-co-infection-ms/rhs_age.R', chdir = TRUE)

# Age structure information, used to calculate mortality rates in susceptibles. 
relageall = c(0.137, rep(0.368/4, 4), rep(0.185/4, 4),  # Jolles 2005, set max age at 18
	rep(0.235/6, 6), rep(0.075/3, 3))					# Also in Caron et al. from 2001 KNP
relage = c(sum(relageall[1:3]), relageall[4],  
	sum(relageall[5:14]), sum(relageall[15:length(relageall)]) ) 
	
#############################################################
#############################################################
#2) Test plots and grooming functions
#############################################################
#############################################################
# Test plots
#############################################################
S0 = 500*relage; It0 = 0*relage; Ib0 = 0*relage; 
Ic0 = 0*relage; R0 = 0 * relage; Rc0 = 0 * relage
x0 = c(S0, It0, Ib0, Ic0, R0, Rc0)
times <- seq(0, 100, 1)
params.test = c(fixed.params,list(gamma=1/2, betaB = c(0.001, 0.002, 0.001, 0.001), betaT = rep(0.002, 4), rhoT = 1.2, rhoB = 4))

sol <- as.data.frame(ode(x0, times, rhs_age, params.test))
soldd <- as.data.frame(ode(x0, times, rhs_age_logistic, params.test), method = ode45)

#sol.recov<- as.data.frame(ode(x0, times, rhs, params.test.recov))
par(mfrow(1,2))
plot(x = sol$time, y = (sol[,2] + sol[,3] + sol[,4] + sol[,5]), 
	pch = 19, main = "Density Independent")
plot(x = soldd$time, y = (soldd[,2] + soldd[,3] + soldd[,4] + soldd[,5]), 
	pch = 19, main = "Density Dependent")
(sol[99,2] + sol[99,3] + sol[99,4] + sol[99,5])/(sol[98,2] + sol[98,3] + sol[98,4] + sol[98,5])


# want 1.10- 1.2 growth rate with brucellosis
S0 = 500*relage; It0 = 0*relage; Ib0 = 100*relage; 
Ic0 = 0*relage; R0 = 0 * relage; Rc0 = 0 * relage
x0 = c(S0, It0, Ib0, Ic0, R0, Rc0)











# functions required to interpret results
#############################################################
plot_raw_numbers = function(sol){
	plot(sol$time, sol$S, col= "black", type= 'l', ylim = c(0, 1200), ylab = "Number of animals", xlab = "Time (in years)")
	lines(sol$time, sol$It, col= "red")
	lines(sol$time, sol$Ib, col= "blue")
	lines(sol$time, sol$Ic, col= "green")
	lines(sol$time, sol$R, col = "orange")
	lines(sol$time, sol$Rc, col = "pink")
	legend("topright", legend = c("S", "It", "Ib", "Ic", "R", "Rc"), col = c("black", "red", "blue", "green", "orange", "pink"), bty="n", lty = 1)
}

groom_sol = function(sol){
	colnames(sol) <- c("times", "S", "It", "Ib", "Ic", "R", "Rc") 
	sol$N <- sol$S + sol$It + sol$Ib + sol$Ic + sol$R + sol$Rc
	sol$TBprev <- (sol$It + sol$Ic + sol$Rc) / sol$N
	sol$Brucprev <- (sol$Ib + sol$R + sol$Ic + sol$Rc) / sol$N
	sol$propTB_co <- (sol$Ic + sol$Rc) / (sol$It + sol$Ic + sol$Rc)
	sol$propBruc_co <- (sol$Ic + sol$Rc) / (sol$Ib + sol$R + sol$Ic + sol$Rc)
	return(sol)
}









par(mfrow=c(1,2))
plot_raw_numbers(sol)
plot_raw_numbers(sol.recov)


get_starting_eqbruc = function(params){
	x0 <- c(S = params['K'][[1]]-100, It = 0, Ib = 50, Ic = 0, R = 50, Rc = 0)
	times <- seq(0, 1000, 1)
	sol <- as.data.frame(ode(x0, times, rhs, params))
	out <- c(S = sol$S[length(times)], 
		It = 0, Ib = sol$Ib[length(times)], Ic = 0,
		R = sol$R[length(times)], Rc = 0)
	return(out)
}

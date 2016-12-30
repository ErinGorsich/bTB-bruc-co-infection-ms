#############################################################
#############################################################
# Erin Gorsich
# This Code Reads in and runs an age-structured co-infection 
# model (defined in rhs_age) for a range of parameters 
#############################################################
#############################################################
#############################################################
# Outline: 
# 1) Load fixed parameters, model
# 2) Set-up features of aging; Functions for plotting
# 3) Getz / Generalized Beverton-Holt form of density dependence
# 4) Ricker form of density dependence
#############################################################
#############################################################
#############################################################
# This code owes much to:
#  http://ms.mcmaster.ca/~bolker/eeid/2011_eco/waifw.pdf
# King & Wearing, Age Structured Models
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
set.seed(5)
# get fixed.params & fixed.params.recov
source('~/GitHub/bTB-bruc-co-infection-ms/fixed_parameters_norecovery_agematrix.R', chdir = TRUE)
source('~/GitHub/bTB-bruc-co-infection-ms/fixed_parameters_recovery_agematrix.R', chdir = TRUE)
# rhs function, determinitic model, age structure
source('~/GitHub/bTB-bruc-co-infection-ms/rhs_age.R', chdir = TRUE)

#############################################################
#############################################################
#2) Set-up features of aging, Functions for plotting
#############################################################
#############################################################
# age divisions in rhs function
#(age= 1-3.9, 4-4.9, 5-14.9, 15+)..subsume calf mortality in births
s_index <- 1:20
it_index <- 21:40
ib_index <- 41:60
ic_index <- 61:80
r_index <- 81:100
rc_index <- 101:120
juveniles <- 1:3
subadult<- 4
adult <- 5:14
mature <- 15:20

# Age structure information, used to calculate mortality rates in susceptibles. 
relageall = c(0.137, rep(0.368/4, 4), rep(0.185/4, 4),  # Jolles 2007, set max age at 20
	rep(0.235/6, 6), rep(0.075/5, 5))					# Also in Caron et al. from 2001 KNP
	
relage = relageall

plot_agestructure = function(x){ # true based on Jolles et al. 2007
	 xcounts <- NA
	 if (length(x) != 120){
	 	print("The age structure should include 20 ages,
	 		for 6 disease classes, giving 120 columns")
	 }
	 for(i in 1:20){
	 	xcounts[i] <- (x[i] + x[20+i] + x[40+i] + x[60+i] + x[80+i] + x[100+i])/sum(x)
	 }
	 d<-matrix(c(relageall, 
	 	xcounts), nrow=2, byrow=TRUE, 
	 	dimnames=list(c("Observed", "Predicted"), 
	 	c(seq(1:20))))
	barplot(d, beside = TRUE, col = c("light gray", "dark gray"), ylab = "Frequency (%)")
	box(	)
	legend("topleft", legend = c("Observed", "Predicted"), fill = c("light gray", "dark gray"))
}
#plot_agestructure(x = seq(1, 20))

plot_dz_agestructure = function(x, dz){
	# Disease options =
	# "bruc" = brucellosis only
	# "tb" = bTB only
	# "co" = all subtypes
	# x = a row in sol
	if (length(x) != 120){
	 	print("The age structure should include 20 ages, 
	 	for 6 disease classes, giving 120 columns")
	}
	
	Sp <- x[1:20]/apply(x, 1, sum); colnames(Sp)<- seq(1:20)
	Itp <- x[21:40]/apply(x, 1, sum); colnames(Itp)<- seq(1:20)
	Ibp <- x[41:60]/apply(x, 1, sum); colnames(Ibp)<- seq(1:20)
	Icp <- x[61:80]/apply(x, 1, sum); colnames(Icp)<- seq(1:20)
	Rp <- x[81:100]/apply(x, 1, sum); colnames(Rp)<- seq(1:20)
	Rcp <- x[101:120]/apply(x, 1, sum); colnames(Rcp)<- seq(1:20)
	
	mat <- as.matrix(rbind(Sp, Itp, Ibp, Icp, Rp, Rcp))
	mat[is.na(mat)]<- 0
	barplot(mat, # columns = age, rows = proportions
	xlab = 'age', main = "Population structure")	

}

plot_raw_numbers = function(sol){
	plot(sol$time, apply(sol[s_index+1], 1, sum), col= "black",
		type= 'l', ylim = c(0, 800), ylab = "Number of animals", 
		xlab = "Time (in years)")
	lines(sol$time, apply(sol[it_index+1], 1, sum), col= "red")
	lines(sol$time, apply(sol[ib_index+1], 1, sum), col= "blue")
	lines(sol$time, apply(sol[ic_index+1], 1, sum), col= "green")
	lines(sol$time, apply(sol[r_index+1], 1, sum), col = "orange")
	lines(sol$time, apply(sol[rc_index+1], 1, sum), col = "pink")
	legend("topright", legend = c("S", "It", "Ib", "Ic", "R", "Rc"),
		col = c("black", "red", "blue", "green", "orange", "pink"), 
		bty="n", lty = 1)
}

get_prevalence = function(sol){
	S <- sum(sol[length(sol) , s_index+1])
	It <- sum(sol[length(sol) , it_index +1])
	Ib <- sum(sol[length(sol) , ib_index +1])
	Ic <- sum(sol[length(sol) , ic_index +1])
	R <- sum(sol[length(sol) , r_index +1])
	Rc <- sum(sol[length(sol) , rc_index +1])
	N <- sum(sol[length(sol), 2:121])
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

plot_ageprevalence = function(sol){
	S <- rep(0, 20); It <- rep(0, 20); Ib <- rep(0, 20);
	Ic <- rep(0, 20); R <- rep(0, 20); Rc <- rep(0, 20)
	N <- rep(0, 20)
	for (i in 1:20){
		S[i] <- sum(sol[length(sol) , i + 1])
		It[i] <- sum(sol[length(sol) , i + 21])
		Ib[i] <- sum(sol[length(sol) , 41 +i])
		Ic[i] <- sum(sol[length(sol) , 61 +i])
		R[i] <- sum(sol[length(sol) , 81 +i])
		Rc[i] <- sum(sol[length(sol) , 101 +i])
		N[i] <- sum(sol[length(sol), c(1+i, 21+i, 41+i, 61 + i, 81+i, 101+i)])
	}
	
	prevB <- (Ib + Ic + R + Rc) / N
	prevT <- (It + Ic + Rc)/ N
	prevBinS <- (Ib + R) / (S + Ib + R)
	prevBinT <- (Ic + Rc) / (It + Ic + Rc) 
	
	overall_prevB <- (sum(Ib) + sum(Ic) + sum(R) + sum(Rc)) / sum(N)
	overall_prevT <- (sum(It) + sum(Ic) + sum(Rc)) / sum(N)
	overall_prevBinS <- (sum(Ib) + sum(R)) / (sum(S) + sum(Ib) + sum(R))
	overall_prevBinT <- (sum(Ic) + sum(Rc)) / (sum(It) + sum(Ic) + sum(Rc))
	overallN = sum(N) 
	
	par(mfrow = c(1,2))
	plot(y = prevB, x= seq(1, 20, 1), type = "b", col = "dark blue", ylim = c(0, 0.8),
		ylab = "Prevalence", xlab = "Age", pch = 19, 
		main = paste("Overall prevalences, Br =", round(overall_prevB, 3), 
		" TB = ", round(overall_prevT, 3) ))
	points(y = prevT, x = seq(1,20,1), type = "b", col = "dark red", pch = 19)
	legend("bottomright", bty = "n", legend = c("Bruc", "TB"), 
		pch = c(19, 19), col = c("dark blue", "dark red"))
	plot(y = prevBinS, x = seq(1,20,1), type = "b", col = "dark blue", pch = 19, 
		ylab = "Brucellosis prevalence", xlab = "Age", ylim = c(0, 0.8), 
		main = paste("Br|S =", round(overall_prevBinS, 3), 
		" Br|Co = ", round(overall_prevBinT, 3) )
		)
	text(x = 10, y = 0.55, labels = paste("Final N = ", round(overallN, 2)))
	points(y = prevBinT, x = seq(1,20,1), type = "b", col = "dark red", pch = 19)
}

#############################################################



#############################################################
#############################################################
# 3) Getz / Generalized Beverton-Holt form of density dependence
#############################################################
#############################################################

# Figure out parameters that give reasonable age structure with Getz density dependence
# Test plots, with no Disease, none takes off
#############################################################
thetaL = seq(0.1, 0.9, by = 0.1)
thetaH = seq(1.1, 1.9, by = 0.1)
N = seq(1, 2000, 1)
f_N = function(N, theta){
	0.5 / (1 + ((N/1000)^theta))
}
plot(x = N, y = f_N(N, theta = 1), type = "l", ylab = "R(N)")
abline(v = 1000, col = "dark red")
for (i in 1:length(thetaL)){
	lines(x = N, y = f_N(N, theta = thetaL[i]), type = "l", col = "light gray", lty = 3)
	lines(x = N, y = f_N(N, theta = thetaH[i]), type = "l", col = "light gray", lty = 5)
}


# params and inits
S0 = 400*relage; It0 = 0*relage; Ib0 = 0*relage; 
Ic0 = 0*relage; R0 = 0 * relage; Rc0 = 0 * relage
x0 = c(S0, It0, Ib0, Ic0, R0, Rc0)
times <- seq(0, 500, 1)


params.test_log = c(fixed.params, list(gamma=1/2, betaB = 0.01,
	betaT = 0.0001, rhoT = 1.2, rhoB = 4, theta= 4, K = 433))
params.test.recov_log = c(fixed.params.recov, list(gamma=1/2, 
	betaB = 0.0001, betaT = 0.001, rhoT = 1.2, rhoB = 4, theta = 4, K = 433))

sol <- as.data.frame(ode(x0, times, rhs_age_matrix, params.test_log))
sol.recov<- as.data.frame(ode(x0, times, rhs_age_matrix, params.test.recov_log))

par(mfrow = c(2,2))
plot(x = sol$time, y = apply(sol[c(2:120)], 1, sum), 
	pch = 19, main = "Density Dependent, no Recovery")
plot(x = sol.recov$time, y = apply(sol.recov[c(2:120)], 1, sum), 
	pch = 19, main = "Density Dependent, Recovery")
plot_raw_numbers(sol)
plot_agestructure(as.matrix(sol[101,c(2:121)]))
stable_age <- unname(unlist( sol[500, c(2:21)]/sum(sol[500, c(2:21)]) ))


# get an idea of the abruptness parameter: 
theta_temp = seq(1, 10, 1)
plot(NA, ylim = c(600, 1500), xlim = c(0, 500), xlab = "Time", ylab = "N")
for (i in 1:length(theta_temp)){
	params.test_log = c(fixed.params, list(gamma=1/2, betaB = 0.01,
		betaT = 0.0001, rhoT = 1.2, rhoB = 4, theta= theta_temp[i], K = 1000))
	sol <- as.data.frame(ode(x0, times, rhs_age_matrix, params.test_log))
	lines(x = sol$time, y = apply(sol[c(2:120)], 1, sum), type = "l",
		 col = i, lty = 3)
}


# Test 2: Add brucellosis, only get brucellosis 
#############################################################
S0 = 400* stable_age; It0 = 0 * stable_age; Ib0 = 20* stable_age; 
Ic0 = 0* stable_age; R0 = 30 * stable_age; Rc0 = 0 * stable_age
x0 = c(S0, It0, Ib0, Ic0, R0, Rc0)
betaB_temp <- seq(0.5, 1.5, 0.01) # slow
betaB_temp <- betaB_temp[1:75]
prevB <- NA; prevBrecov <- NA
for (i in 1:length(betaB_temp)){
	params.test_log = c(fixed.params, list(gamma=1/2, theta = 4, K = 433,
		betaB = betaB_temp[i], betaT = 0.001, rhoT = 1, rhoB = 4))
	params.test.recov_log = c(fixed.params.recov, list(gamma=1/2, theta = 4, K = 433,
		betaB = betaB_temp[i], betaT = 0.001, rhoT = 1.2, rhoB = 4))
	sol <- as.data.frame(ode(x0, times, rhs_age_matrix, params.test_log))
	sol.recov<- as.data.frame(ode(x0, times, rhs_age_matrix, params.test.recov_log))
	prevB[i]<- get_prevalence(sol)$prevB
	prevBrecov[i] <- get_prevalence(sol.recov)$prevB
	rm(sol, sol.recov, params.test_log, params.test.recov_log)
}
plot(x = betaB_temp, y = prevB, type = "l", xlab = expression(beta), ylab = "Brucellosis prevalence", main = "Single disease")
#lines(x = betaB_temp, y = prevB, type = "l", col = "dark blue") # same pattern
abline(h = c(0.1, 0.2, 0.3, 0.4), col = "dark red")
# Beta value at:
betaB_temp[which(prevB < 0.06 & prevB > 0.04)]  		# 0.79
betaB_temp[which(prevB < 0.11 & prevB > 0.09)]		# 0.825
betaB_temp[which(prevB < 0.21 & prevB > 0.19)]		# 0.915
betaB_temp[which(prevB < 0.31 & prevB > 0.29)]		# 1.025
betaB_temp[which(prevB < 0.405 & prevB > 0.395)]	# 1.185

# set at 30% Brucellosis prevalence wihtout bTB, betaB = 1.025
params.test_log = c(fixed.params, list(gamma=1/2, theta = 4, K = 433,
	betaB = 1.025, betaT = 0.001, rhoT = 1, rhoB = 4))
params.test.recov_log = c(fixed.params.recov, list(gamma=1/2, theta = 4, K = 433,
	betaB = 1.025, betaT = 0.001, rhoT = 1, rhoB = 4))

times <- seq(0, 500, 1)
sol <- as.data.frame(ode(x0, times, rhs_age_matrix, params.test_log))
sol.recov<- as.data.frame(ode(x0, times, rhs_age_matrix, params.test.recov_log))

par(mfrow = c(2,2))
plot(x = sol$time, y = apply(sol[c(2:120)], 1, sum), 
	pch = 19, main = "Density Dependent, no Recovery")
plot(x = sol.recov$time, y = apply(sol.recov[c(2:120)], 1, sum), 
	pch = 19, main = "Density Dependent, Recovery")
#plot_age_prev_by_coinfection(sol)
plot_raw_numbers(sol)
plot_raw_numbers(sol.recov)
get_prevalence(sol); 
################
get_prevalence(sol.recov)
plot_agestructure(as.matrix(sol[101,c(2:121)])) 
plot_ageprevalence(sol)

# proportions
endemic_agestructure_p <- unname(unlist( sol[500, c(2:121)]/sum(sol[500, c(2:121)]) ))
endemic_agestructure_recov_p <- unname(unlist( sol.recov[500, c(2:121)]/sum(sol.recov[500, c(2:121)]) ))
#raw numbers
endemic_agestructure <- unname(unlist( sol[500, c(2:121)] ))
endemic_agestructure_recov <- unname(unlist( sol.recov[500, c(2:121)]))


# Test 3: Add bTB, only get bTB --> works! 
#############################################################
S0 = 400* stable_age; It0 = 20 * stable_age; Ib0 = 0* stable_age; 
Ic0 = 0* stable_age; R0 = 0 * stable_age; Rc0 = 0 * stable_age
x0 = c(S0, It0, Ib0, Ic0, R0, Rc0)

betaT_temp <- seq(0.0001, 0.0009, 0.00001) 
prevT <- NA; prevTrecov <- NA
for (i in 1:length(betaT_temp)){
	params.test_log = c(fixed.params, list(gamma=1/2, theta = 4, K = 433,
		betaB = 1.025, betaT = betaT_temp[i], rhoT = 1, rhoB = 4))
	sol <- as.data.frame(ode(x0, times, rhs_age_matrix, params.test_log))
	prevT[i]<- get_prevalence(sol)$prevTB
	rm(sol,  params.test_log)
}
plot(x = betaT_temp, y = prevT, type = "l", xlab = expression(beta), ylab = "bTB prevalence", main = "Single disease")
abline(h = c(0.1, 0.2, 0.3, 0.4), col = "dark red")
# BetaT value at:
betaT_temp[which(prevT < 0.12 & prevT > 0.08)]		# 0.000395
betaT_temp[which(prevT < 0.21 & prevT > 0.19)]	 	# 0.00044
betaT_temp[which(prevT < 0.31 & prevT > 0.29)] 	# 0.000505
betaT_temp[which(prevT < 0.405 & prevT > 0.395)]	# 0.00059

# set to give 30% ??????
params.test_log = c(fixed.params, list(gamma=1/2, theta = 4, K = 433,
	betaB = 1.025, betaT = 0.000223, rhoT = 1, rhoB = 4))
sol <- as.data.frame(ode(x0, times, rhs_age_matrix, params.test_log))
par(mfrow = c(1,3))
plot(x = sol$time, y = apply(sol[c(2:120)], 1, sum), 
	pch = 19, main = "Density Dependent, no Recovery")
#plot_age_prev_by_coinfection(sol)
plot_raw_numbers(sol)
get_prevalence(sol); 
################
plot_agestructure(as.matrix(sol[101,c(2:121)])) 



# Test 4: Add co-infection (Introduce bTB after set levels of brucellosis)
#############################################################
# endemic_agestructure is set to final prevalence/age structure in bruc only model
x0 = endemic_agestructure
x0recov = endemic_agestructure_recov
x0[28] <- 5; x0[8] <- x0[8] - 5
x0recov[28] <- 5; x0recov[8] <- x0recov[8] - 5
times <- seq(0, 500, 1)

# Gives 30% bruc prev without bTB; Gives 20% bTB prevalence without brucellosis
params.test_log = c(fixed.params, list(gamma=1/2, theta = 4, K = 433,
	betaB = 1.025, betaT = 0.00044, rhoT = 1, rhoB = 4))
params.test.recov_log = c(fixed.params.recov, list(gamma=1/2, theta = 4, K = 433,
	betaB = 1.025, betaT = 0.00044, rhoT = 1, rhoB = 4))

sol <- as.data.frame(ode(x0, times, rhs_age_matrix, params.test_log))
sol.recov<- as.data.frame(ode(x0recov, times, rhs_age_matrix, params.test.recov_log))

par(mfrow = c(1,2))
plot_raw_numbers(sol)
plot_raw_numbers(sol.recov)
get_prevalence(sol); get_prevalence(sol.recov);   # almost no change in Bruc prev, no bTB!

# Make Summary Plots
x0 = endemic_agestructure
x0recov = endemic_agestructure_recov
x0[28] <- 5; x0[8] <- x0[8] - 5
x0recov[28] <- 5; x0recov[8] <- x0recov[8] - 5

rhoB_test <- seq(1, 10, 0.5)
rhoT_test <- seq(1, 10, 0.5)
epi <- data.frame(
	rhoB= rep(rhoB_test, length(rhoB_test)), 
	rhoT = rep(rhoT_test, each = length(rhoB_test)), 
	bTBprev = NA, brucprev = NA, 	finalN = NA, 
	bTB_inS = NA, bTB_inB = NA, bruc_inS = NA, bruc_inTB = NA,
	rbTBprev = NA, rbrucprev = NA, rfinalN = NA, 
	rbTB_inS = NA, rbTB_inB = NA, rbruc_inS = NA, rbruc_inTB = NA)
	
for (i in 1:length(epi[,1])){
	params.test_log = c(fixed.params, list(gamma=1/2, theta = 4, K = 433,
		betaB = 1.025, betaT = 0.00044, rhoT = epi$rhoT[i], rhoB = epi$rhoB[i]))
	params.test.recov_log = c(fixed.params.recov, list(gamma=1/2, theta = 4, K = 433,
		betaB = 1.025, betaT = 0.00044, rhoT = epi$rhoT[i], rhoB = epi$rhoB[i]))
	sol <- as.data.frame(ode(x0, times, rhs_age_matrix, params.test_log))
	sol.recov<- as.data.frame(ode(x0recov, times, rhs_age_matrix, params.test.recov_log))
	temp <- get_prevalence(sol); rtemp <- get_prevalence(sol.recov)
	
	epi$bTBprev[i] = temp$prevTB
	epi$brucprev[i] = temp$prevB 	
	epi$finalN[i] = sum(sol[length(sol), c(2:121)])
	epi$bTB_inS[i] = temp$prevTinS 
	epi$bTB_inB[i] = temp$prevTinB
	epi$bruc_inS[i] = temp$prevBinS 
	epi$bruc_inTB[i] = temp$prevBinT

	epi$rbTBprev[i] = rtemp$prevTB 
	epi$rbrucprev[i] = rtemp$prevB 	
	epi$rfinalN[i] = sum(sol.recov[length(sol.recov), c(2:121)])
	epi$rbTB_inS[i] = rtemp$prevTinS  
	epi$rbTB_inB[i] = rtemp$prevTinB 
	epi$rbruc_inS[i] = rtemp$prevBinS
	epi$rbruc_inTB[i] = rtemp$prevBinT
	rm(params.test_log, params.test.recov_log, sol, sol.recov, temp, rtemp)
}

write.csv(epi, "~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/draft2/post-labmeeting/vary_rho.csv")
# Add lines for what it looks like without brucellosis around
# choose diverging palette (... but needs turned into range...col.region = brewer.pal(8, "RdYlBu")), centered on non-co-infected values
epi <- read.csv("~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/draft2/post-labmeeting/vary_rho/vary_rho.csv") 

# bTB prevalence (vs. 30% alone)
p1 <- levelplot(bTBprev~rhoB*rhoT, data = epi, main = "bTB prevalence, no recovery", 
	xlab = expression(paste(beta[B]^{"'"}, "/", beta[B]) ), 
	ylab = expression(paste(beta[TB]^{"'"}, "/", beta[TB])), 
	at = seq(0, max(epi$bTBprev, epi$rbTBprev), 0.02 ))
# better than: frac(beta[B], beta[b]) )) )

p2 <- levelplot(rbTBprev~rhoB*rhoT, data = epi,  main = "bTB prevalence, recovery", 
	xlab = expression(paste(beta[B]^{"'"}, "/", beta[B]) ), 
	ylab = expression(paste(beta[TB]^{"'"}, "/", beta[TB])), 
	at = seq(0, max(epi$bTBprev, epi$rbTBprev), 0.02 ))
grid.arrange(p1, p2, ncol = 2)

# brucellosis prevalence (vs 30% alone)
p3 <- levelplot(brucprev~rhoB*rhoT, data = epi, main = "bruc prevalence, no recovery", 
	xlab = expression(paste(beta[B]^{"'"}, "/", beta[B]) ), 
	ylab = expression(paste(beta[TB]^{"'"}, "/", beta[TB])), 
	at = seq(0, max(epi$brucprev, epi$rbrucprev), 0.02 ))
	
p4 <- levelplot(rbrucprev ~rhoB*rhoT, data = epi,  main = "bruc prevalence, recovery", 
	xlab = expression(paste(beta[B]^{"'"}, "/", beta[B]) ), 
	ylab = expression(paste(beta[TB]^{"'"}, "/", beta[TB])), 
	at = seq(0, max(epi$brucprev, epi$rbrucprev), 0.02))
grid.arrange(p3, p4, ncol = 2)


# final N  ( )
p5 <- levelplot(finalN~rhoB*rhoT, data = epi, main = "Final N, no recovery", 
	xlab = expression(paste(beta[B]^{"'"}, "/", beta[B]) ), 
	ylab = expression(paste(beta[TB]^{"'"}, "/", beta[TB])), 
	at = seq(350, max(epi$finalN, epi$rfinalN), 5) )
p6<- levelplot(rfinalN~rhoB*rhoT, data = epi,  main = "Final N, recovery", 
	xlab = expression(paste(beta[B]^{"'"}, "/", beta[B]) ), 
	ylab = expression(paste(beta[TB]^{"'"}, "/", beta[TB])), 
	at = seq(350, max(epi$finalN, epi$rfinalN), 5) )
grid.arrange(p5, p6, ncol = 2)




library(RColorBrewer)





#############################################################
#############################################################
# 4) Ricker Model  (DOUBLE CHECK FREQ/DENSITY DEPENDENCE ASSUMPTIONS!)
#############################################################
#############################################################


#############################################################
# Figure out parameters that give reasonable age structure with Ricker model
# Test plots, with no Disease, none takes off
# STILL NEED TO CLARIFY ASSUMPTIONS ON AGE STRUCTURE WITH AND WITHOUT DZ
#############################################################
S0 = 1000*relage; It0 = 0*relage; Ib0 = 0*relage; 
Ic0 = 0*relage; R0 = 0 * relage; Rc0 = 0 * relage
x0 = c(S0, It0, Ib0, Ic0, R0, Rc0)
times <- seq(0, 500, 1)
params.test = c(fixed.params, list(gamma=1/2, betaB = 0.01,
	betaT = 0.0001, rhoT = 1.2, rhoB = 4))
params.test.recov = c(fixed.params.recov, list(gamma=1/2, 
	betaB = 0.0001, betaT = 0.001, rhoT = 1.2, rhoB = 4))
sol <- as.data.frame(ode(x0, times, rhs_age_matrix_ricker, params.test))
sol.recov<- as.data.frame(ode(x0, times, rhs_age_matrix_ricker, params.test.recov))

par(mfrow = c(2,2))
plot(x = sol$time, y = apply(sol[c(2:120)], 1, sum), 
	pch = 19, main = "Density Dependent, no Recovery")
plot(x = sol.recov$time, y = apply(sol.recov[c(2:120)], 1, sum), 
	pch = 19, main = "Density Dependent, Recovery")

plot_raw_numbers(sol)
plot_agestructure(as.matrix(sol[101,c(2:121)]))  # better, no density dependence...

stable_age <- unname(unlist( sol[500, c(2:21)]/sum(sol[500, c(2:21)]) ))

# Figure out birth rates taht give age structure above... 
f = function(b, x){
	1 * exp(-b * x)  # x = Nall, b = scale constant
}
# Number of offspring produced (assuming stable equal births/age) at K
300* f(1/100, 300) # 14.93612
300* f(1/500, 300)  # 164.6435
300* f(1/1000, 300)  # 222.2455



# Test 2: Add brucellosis, only get brucellosis --> works! 
#############################################################
params.test = c(fixed.params, list(gamma=1/2, 
	betaB = 0.0015, betaT = 0.001, rhoT = 1.2, rhoB = 4))
params.test.recov = c(fixed.params.recov, list(gamma=1/2, 
	betaB = 0.0015, betaT = 0.001, rhoT = 1.2, rhoB = 4))
S0 = 980* stable_age; It0 = 0* stable_age; Ib0 = 20* stable_age; 
Ic0 = 0* stable_age; R0 = 30 * stable_age; Rc0 = 0 * stable_age
x0 = c(S0, It0, Ib0, Ic0, R0, Rc0)
times <- seq(0, 500, 1)
sol <- as.data.frame(ode(x0, times, rhs_age_matrix_ricker, params.test))
sol.recov<- as.data.frame(ode(x0, times, rhs_age_matrix_ricker, params.test.recov))

par(mfrow = c(2,2))
plot(x = sol$time, y = apply(sol[c(2:120)], 1, sum), 
	pch = 19, main = "Density Dependent, no Recovery")
plot(x = sol.recov$time, y = apply(sol.recov[c(2:120)], 1, sum), 
	pch = 19, main = "Density Dependent, Recovery")
#plot_age_prev_by_coinfection(sol)
plot_raw_numbers(sol)
plot_raw_numbers(sol.recov)
get_prevalence(sol); 
################
get_prevalence(sol.recov)
plot_agestructure(as.matrix(sol[101,c(2:121)])) 
plot_ageprevalence(sol)

endemic_agestructure <- unname(unlist( sol[500, c(2:121)]/sum(sol[500, c(2:121)]) ))
endemic_agestructure_recov <- unname(unlist( sol.recov[500, c(2:121)]/sum(sol.recov[500, c(2:121)]) ))


# Test 3: Add bTB, only get bTB --> Check! 
#############################################################	
params.test = c(fixed.params, list(gamma=1/2,
	betaB = 0.002, betaT = 0.0005, rhoT = 1.2, rhoB = 4))
params.test.recov = c(fixed.params.recov, list(gamma=1/2, 
	betaB = 0.002, betaT = 0.0005, rhoT = 1.2, rhoB = 4))
S0 = 980* stable_age; It0 = 2* stable_age; Ib0 = 0* stable_age; 
Ic0 = 0* stable_age; R0 = 0 * stable_age; Rc0 = 0 * stable_age
x0 = c(S0, It0, Ib0, Ic0, R0, Rc0)
times <- seq(0, 500, 1)
sol <- as.data.frame(ode(x0, times, rhs_age_matrix_ricker, params.test))
sol.recov<- as.data.frame(ode(x0, times, rhs_age_matrix_ricker, params.test.recov))

par(mfrow = c(2,2))
plot(x = sol$time, y = apply(sol[c(2:120)], 1, sum), 
	pch = 19, main = "Density Dependent, no Recovery")
plot(x = sol.recov$time, y = apply(sol.recov[c(2:120)], 1, sum), 
	pch = 19, main = "Density Dependent, Recovery")
#plot_age_prev_by_coinfection(sol)
plot_raw_numbers(sol)
plot_raw_numbers(sol.recov)
get_prevalence(sol); 
################
get_prevalence(sol.recov)

plot_agestructure(as.matrix(sol[101,c(2:121)])) 
plot_ageprevalence(sol)

	
	

#############################################################
#############################################################
#2) Add bTB to EE population... 
#############################################################
#############################################################
params.test = c(fixed.params, list(gamma=1/2, 
	betaB = 0.002, betaT = 0.0005, rhoT = 1.2, rhoB = 2))  
params.test.recov = c(fixed.params.recov, list(gamma=1/2, 
	betaB = 0.002, betaT = 0.0005, rhoT = 1.2, rhoB = 2))
# betaT = 0.005 gives bTB prev ~25% without bruc; essentially no prevalence with bruc.


x0 = 1000 * endemic_agestructure
x0recov = 1000 * endemic_agestructure_recov
x0[21:40] <- 5
x0recov[21:40] <- 5
times <- seq(0, 500, 1)
sol <- as.data.frame(ode(x0, times, rhs_age_matrix_ricker, params.test))
sol.recov<- as.data.frame(ode(x0recov, times, rhs_age_matrix_ricker, params.test.recov))

par(mfrow = c(2,2))
plot(x = sol$time, y = apply(sol[c(2:120)], 1, sum), 
	pch = 19, main = "Density Dependent, no Recovery")
plot(x = sol.recov$time, y = apply(sol.recov[c(2:120)], 1, sum), 
	pch = 19, main = "Density Dependent, Recovery")
#plot_age_prev_by_coinfection(sol)
plot_raw_numbers(sol)
plot_raw_numbers(sol.recov)
get_prevalence(sol); 
################
get_prevalence(sol.recov)

par(mfrow = c(1,2))
plot_agestructure(as.matrix(sol[101,c(2:121)])) 
plot_ageprevalence(sol)
plot_ageprevalence(sol.recov)

# only works without recovery: But without recovery, we get waay less bTB than without brucellosis
params.test = c(fixed.params, list(gamma=1/2, 
	betaB = 0.002, betaT = 0.0006, rhoT = 1, rhoB = 1))  
sol <- as.data.frame(ode(x0, times, rhs_age_matrix_ricker, params.test))
plot_ageprevalence(sol)
get_prevalence(sol)


#############################################################
#############################################################
#2) Looop through transmission parameters to see...
#############################################################
#############################################################
x0 = 1000 * endemic_agestructure
x0recov = 1000 * endemic_agestructure_recov
x0[21:40] <- 5
x0recov[21:40] <- 5
times <- seq(0, 500, 1)
modeloutput <- data.frame(
	betaTvec = rep(c(0.0005, 0.001, 0.0015, 0.002, 0.0025, 0.003, 0.0035), 8),
	rhoB = rep(c(1, 2, 3, 4), each = 14),
	rhoT = rep(c(rep(1, 7), rep(1.3, 7)), 4), 
	brucprev = NA, tbprev = NA, 
	brucprevins = NA, brucprevintb = NA, 
	brucprevR = NA, tbprevR = NA, 
	brucprevinsR = NA, brucprevintbR = NA)

setwd("~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/draft2/figures/model_output_agestructure")
for (i in 1:length(modeloutput[,1])){
	print (modeloutput$betaTvec[i])
	params.test = c(fixed.params, list(gamma=1/2, 
		betaB = 0.0025, betaT = modeloutput$betaTvec[i], 
		rhoT = modeloutput$rhoT[i], rhoB = modeloutput$rhoB[i]))  
	params.test.recov = c(fixed.params.recov, list(gamma=1/2, 
		betaB = 0.0025, betaT = modeloutput$betaTvec[i],
		 rhoT = modeloutput$rhoT[i], rhoB = modeloutput$rhoB[i]))

	sol <- as.data.frame(ode(x0, times, rhs_age_matrix_ricker, params.test))
	sol.recov<- as.data.frame(ode(x0recov, times, rhs_age_matrix_ricker, 
		params.test.recov))
	
	temp <- get_prevalence(sol)
	tempR <- get_prevalence(sol.recov)
	
	modeloutput$brucprev[i] <- temp$prevB
	modeloutput$brucprevR[i] <- tempR$prevB
	modeloutput$tbprev[i] <- temp$prevTB
	modeloutput$tbprevR[i] <- tempR$prevTB
	modeloutput$brucprevins[i] <- temp$prevBinS
	modeloutput$brucprevinsR[i] <- tempR$prevBinS
	modeloutput$brucprevintb[i] <- temp$prevBinT
	modeloutput$brucprevintbR[i] <- tempR$prevBinT

	jpeg(paste("Norecov_betaT=", modeloutput$betaTvec[i], 
		"_rho_", modeloutput$rhoB[i], "_", modeloutput$rhoT[i], 
		".jpg", sep = ""), width = 650, height = 480, units = "px")		
	plot_ageprevalence(sol)
	dev.off()
	
	jpeg(paste("Recov_betaT=", modeloutput$betaTvec[i], 
		"_rho_", modeloutput$rhoB[i], "_", modeloutput$rhoT[i],
		".jpg", sep = ""), , width = 650, height = 480, units = "px")
	plot_ageprevalence(sol.recov)
	dev.off()
}



# Age structured model owes much to:  http://ms.mcmaster.ca/~bolker/eeid/2011_eco/waifw.pdf

#############################################################
#############################################################
#1) Load fixed parameters, model
#############################################################
#############################################################
rm(list = ls())
require("deSolve")
library("plyr")
library("ggplot2")
set.seed(5)
# get fixed.params & fixed.params.recov
source('~/GitHub/bTB-bruc-co-infection-ms/fixed_parameters_norecovery_agematrix.R', chdir = TRUE)
source('~/GitHub/bTB-bruc-co-infection-ms/fixed_parameters_recovery_agematrix.R', chdir = TRUE)
# rhs function, determinitic model, age structure
source('~/GitHub/bTB-bruc-co-infection-ms/rhs_age.R', chdir = TRUE)


#############################################################
#############################################################
#2) Set-up features of aging
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
	 	print("The age structure should include 20 ages, for 6 disease classes, giving 120 columns")
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
plot_agestructure(x = seq(1, 20))

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
		type= 'l', ylim = c(0, 1100), ylab = "Number of animals", 
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
	return(list(prevTB = prevTB, prevB = prevB,
		prevBinS = prevBinS, prevBinT = prevBinT))
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
# Figure out parameters that give reasonable age structure 
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
# Need to make that = deaths at stable age distribution...


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



##########################################################
##########################################################
##########################################################
# Model sensitivity and evaluation 
##########################################################
##########################################################
##########################################################
# Outline
# 1) Run with other form of density dependence (logistic)
# 2) Do whole process with longer infection duraiton (1/gamma = 3 yr; 4 yr)
# 2) Latin hypercube PRCC
##########################################################

##########################################################
##########################################################
# 1) Run with other form of density dependence (logistic)
##########################################################
##########################################################
rm(list = ls())
require("deSolve")

source('~/GitHub/bTB-bruc-co-infection-ms/fixed_parameters_norecovery_agematrix.R', chdir = TRUE)
source('~/GitHub/bTB-bruc-co-infection-ms/rhs_age.R', chdir = TRUE)

# Run B-H model with no disease (endemic population size = 609)
s_index <- 1:20
it_index <- 21:40
ib_index <- 41:60
ic_index <- 61:80
r_index <- 81:100
rc_index <- 101:120
relageall = c(0.137, rep(0.368/4, 4), rep(0.185/4, 4),  
	rep(0.235/6, 6), rep(0.075/5, 5))					
relage = relageall
S0 = 400*relage; It0 = 0*relage; Ib0 = 0*relage; 
Ic0 = 0*relage; R0 = 0 * relage; Rc0 = 0 * relage
x0 = c(S0, It0, Ib0, Ic0, R0, Rc0)
times <- seq(0, 1000, 1)  
params.test <- c(fixed.params, list(gamma=1/2, theta = 4, K = 433,
	betaB = 0.6087396, betaT = 1.2974554/1000, rhoT = 1, rhoB = 2.1))
sol <- as.data.frame(ode(x0, times, rhs_age_matrix, params.test))
stable_age <- unname(unlist( sol[1000, c(2:21)]/sum(sol[500, c(2:21)]) ))

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

plot_raw_numbers(sol)
sum(sol[1001,c(2:length(sol))])

# Set K s.t. get same popuation size/age structure (endemic population size = 609.0097)
K <- seq(1030, 1038, 0.02)
S0 = 400*relage; It0 = 0*relage; Ib0 = 0*relage; 
Ic0 = 0*relage; R0 = 0 * relage; Rc0 = 0 * relage
x0 = c(S0, It0, Ib0, Ic0, R0, Rc0)
times <- seq(0, 1000, 1)  
val <- NA
for (i in 1:length(K)){
	params.test <- c(fixed.params, list(gamma=1/2, theta = 4, K = K[i],
		betaB = 0.6087396, betaT = 1.2974554/1000, rhoT = 1, rhoB = 2.1))
	sol <- as.data.frame(ode(x0, times,  rhs_age_logistic, params.test))
	val[i] <- sum(sol[1001,c(2:length(sol))])	
}
plot(x = K, y = val, type = "b", pch = 19)
abline(h = 609.0097) # K = 1034.99

val <- NA
K <- seq(382, 383, 0.05)
for (i in 1:length(K)){
	params.test <- c(fixed.params, list(gamma=1/2, theta = 4, K = K[i],
		betaB = 0.6087396, betaT = 1.2974554/1000, rhoT = 1, rhoB = 2.1))
	sol <- as.data.frame(ode(x0, times,  rhs_age_matrix_ricker, params.test))
	val[i] <- sum(sol[1001,c(2:length(sol))])	
}
plot(x = K, y = val[1:length(K)], type = "b", pch = 19)
abline(h = 609.0097) # K = 382.6 

# Make plots
###############################
n = 2
source('~/GitHub/bTB-bruc-co-infection-ms/get_EE.R', chdir = TRUE)
set.seed(1)

EE_bTB_single_list <- NA; EE_bTB_co_list <- NA
EE_brucellosis_single_list <- NA; EE_brucellosis_co <- NA
rEE_bTB_single_list <- NA; rEE_bTB_co_list <- NA
rEE_brucellosis_single_list <- NA; rEE_brucellosis_co <- NA
lEE_bTB_single_list <- NA; lEE_bTB_co_list <- NA
lEE_brucellosis_single_list <- NA; lEE_brucellosis_co <- NA

rhoB <- NA; mortB <- NA; mortT <- NA; mortC <- NA
pb <- txtProgressBar(min = 0, max = n, style = 3, char = "=")

for(i in 1:n){
	if(i %% 5 == 0){
		setTxtProgressBar(pb, i)
	}
	params <- c(fixed.params, list(gamma=1/2, theta = 4, K = 433,
		betaB = 0.6087396, betaT = 1.2974554/1000, rhoT = 1, rhoB = 2.1))
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

	paramslog <- params; paramslog$K <- 1034.99
	vallog <- getEE(paramslog)
	
	paramsricker <- params; paramsricker$K <- 382.6
	valricker <- getEE(paramsricker)

	EE_bTB_single_list[i] <- val[1]; EE_bTB_co_list[i] <- val[2]
	EE_brucellosis_single_list[i] <- val[3]; 	EE_brucellosis_co[i] <- val[4]
	rEE_bTB_single_list[i] <- val[1]; rEE_bTB_co_list[i] <- val[2]
	rEE_brucellosis_single_list[i] <- val[3]; rEE_brucellosis_co[i] <- val[4]
	lEE_bTB_single_list[i] <- val[1]; lEE_bTB_co_list[i] <- val[2]
	lEE_brucellosis_single_list[i] <- val[3]; lEE_brucellosis_co[i] <- val[4]
}
senslist <- list(
	EE_bTB_single_list = EE_bTB_single_list, 
	EE_bTB_co_list = EE_bTB_co_list,
	EE_brucellosis_single_list = EE_brucellosis_single_list, 
	EE_brucellosis_co = EE_brucellosis_co,
	rEE_bTB_single_list = rEE_bTB_single_list, 
	rEE_bTB_co_list = rEE_bTB_co_list,
	rEE_brucellosis_single_list = rEE_brucellosis_single_list, 
	rEE_brucellosis_co = rEE_brucellosis_co,	
	lEE_bTB_single_list = lEE_bTB_single_list, 
	lEE_bTB_co_list = lEE_bTB_co_list,
	lEE_brucellosis_single_list = lEE_brucellosis_single_list, 
	lEE_brucellosis_co = lEE_brucellosis_co)
str(senslist)

saveRDS(senslist, "~/GitHub/bTB-bruc-co-infection-ms/sensitivity_densitydependence_simulation_results.rds")
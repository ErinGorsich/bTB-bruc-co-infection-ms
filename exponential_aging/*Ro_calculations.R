rm(list = ls())
library(MASS)
library(deSolve)
#library(gmp)  # precision on matrix inverse calculation, not yet usnig
setwd("~/GitHub/bTB-bruc-co-infection-ms")
###############################################
# Ro of bTB in the absence of brucellosis, with age structure
###############################################
source('~/GitHub/bTB-bruc-co-infection-ms/fixed_parameters_norecovery_agematrix.R', chdir = TRUE)
source('~/GitHub/bTB-bruc-co-infection-ms/rhs_age.R', chdir = TRUE)
params <- c(fixed.params, list(gamma=1/2, betaB = 0.6087,
	betaT = 0.0012974, rhoT = 1, rhoB = 2.1, theta= 4, K = 433))

Ro_bTB_single_age = function(params){
	###################################
	# Input: x = c(S = final, 1*20 S vector at params)
	# Output: Ro
	###################################
	# Get stable age distribution in dz free conditions
	relage = c(0.137, rep(0.368/4, 4), rep(0.185/4, 4),
		rep(0.235/6, 6), rep(0.075/5, 5))
	S0 = 400*relage; It0 = 0*relage; Ib0 = 0*relage; 
	Ic0 = 0*relage; R0 = 0 * relage; Rc0 = 0 * relage
	x0 = c(S0, It0, Ib0, Ic0, R0, Rc0)
	times <- seq(0, 1000, 1)
	sol <- as.data.frame(ode(x0, times, rhs_age_matrix, params))
	S <- unname(unlist( sol[length(sol[,1]), c(2:21)]))

	# Calculate next generation matrix
	Fmat = params$betaT* matrix(rep(S, each = 20), nrow = 20, byrow = T)
	Vmat = diag(x = params$muT + 1)
	Vmat[row(Vmat) - col(Vmat) == 1] <- -1
	Vinv <- solve(Vmat)
	vals <- eigen(Fmat %*% Vinv)$values

	Vmat2 <- diag(x = params$muT)
	V2inv <- solve(Vmat2)
	vals2 <- eigen(Fmat %*% V2inv)$values
	return(max(Re(vals2)))
	#return(list(Ro =vals, Ro2 =vals2, Vmat = Vmat, Vmat2 = Vmat2, 
	#	check = Vmat%*%Vinv, check2 = Vmat2 %*% V2inv, N = sum(S)))
}


###############################################
# Ro of bTB in the absence of brucellosis, with age structure (and age sp. FOI, brucellosis)
###############################################

Ro_bTB_co_age = function(params){
	###################################
	# Input: x = c(S = final, 1*20 S vector at params)
	# Output: Ro
	###################################
	# Get stable age distribution in dz free conditions
	relage = c(0.137, rep(0.368/4, 4), rep(0.185/4, 4),
		rep(0.235/6, 6), rep(0.075/5, 5))
	S0 = 400*relage; It0 = 0*relage; Ib0 = 2*relage; 
	Ic0 = 0*relage; R0 = 2 * relage; Rc0 = 0 * relage
	x0 = c(S0, It0, Ib0, Ic0, R0, Rc0)
	times <- seq(0, 1000, 1)
	sol <- as.data.frame(ode(x0, times, rhs_age_matrix, params))
	S <- unname(unlist( sol[length(sol[,1]), c(2:21)]))
	Ib <- unname(unlist( sol[length(sol[,1]), c(42:61)]))
	R <- unname(unlist( sol[length(sol[,1]), c(82:101)]))
	N <- sum(S) + sum(Ib) + sum(R)
	Sall <- sum(S)
	Iball <- sum(Ib)
	Rall <- sum(R)
	
	#  age specific FOI
	params$betaBm <- rep(params$betaB, length.out = 20)
	params$betaBm[2:5] <- exp(0.885) * params$betaB 

	# Calculate next generation matrix, 1:60 columns, 1:60 rows
	Fmat = rbind(params$betaT * matrix(rep(S, each = 60), nrow = 20, byrow = T),
		params$rhoT * params$betaT * matrix(rep(Ib, each = 60), nrow = 20, byrow = T),
		params$rhoT * params$betaT * matrix(rep(R, each = 60), nrow = 20, byrow = T))
	
	# rows 1:20, cols 1:60
	diag <- (params$rhoB * params$betaBm * Iball / N) + params$muT + 1
	M1 <- diag(x = diag)
	M1[row(M1) - col(M1) == 1] <- -1
	Vmat1 = cbind(M1, 		
		matrix(0, nrow = 20, ncol = 40) )

	# rows 21:40, cols 1:60		 
	diag <- - params$rhoB * params$betaBm * Iball / N   
	M4 <- diag(x =diag)
	
	diag <- params$gamma + params$muC + 1
	M5 <- diag(x = diag)
	M5[row(M5) - col(M5) == 1] <- -1
	M6 <- diag(x = rep(- params$epsilon), 20)
	Vmat2 <- cbind(M4, M5, M6)
	
	# rows 41:60, cols 1:60		 
	M7 <- matrix(0, nrow = 20, ncol = 20)
	M8 <- diag(x = rep(-params$gamma, 20))
	M9 <- diag(x = params$epsilon + params$muRC + 1)
	M9[row(M9 - col(M9) == 1)] <- 1
	Vmat3 <- cbind(M7, M8, M9)
	
	Vmat <- rbind(Vmat1, Vmat2, Vmat3)
	Vinv <- solve(Vmat)
	vals <- eigen(Fmat %*% Vinv)$values
	return(max(Re(vals)))   # 0.9634113
}


###############################################
# Ro of brucellosis in the absence of bTB, age structure
###############################################

Ro_brucellosis_single_age = function(params){
	###################################
	# Input: paramter file
	###################################
	# Get stable age distribution in dz free conditions
	relage = c(0.137, rep(0.368/4, 4), rep(0.185/4, 4),
		rep(0.235/6, 6), rep(0.075/5, 5))
	S0 = 400*relage; It0 = 0*relage; Ib0 = 0*relage; 
	Ic0 = 0*relage; R0 = 0 * relage; Rc0 = 0 * relage
	x0 = c(S0, It0, Ib0, Ic0, R0, Rc0)
	times <- seq(0, 1000, 1)
	sol <- as.data.frame(ode(x0, times, rhs_age_matrix, params))
	S <- unname(unlist( sol[length(sol[,1]), c(2:21)]))

    betaB_age = c(params$betaB, rep(exp(0.885) * params$betaB, 4), rep(params$betaB, 15))
    
	# Calculate next generation matrix
	Fmat = params$betaB * matrix(rep(S/sum(S), each = 20), nrow = 20, byrow = T)
	Fmat[2:5,] <- exp(0.885) * Fmat[2:5,] # increased FOI in young
	
	Vmat = diag(x = params$muB + params$gamma + 1)
	Vmat[row(Vmat) - col(Vmat) == 1] <- -1
	Vinv <- solve(Vmat)
	vals <- eigen(Fmat %*% Vinv)$values
	return(max(Re(vals)))
	#return(list(Ro =vals, Ro2 =vals2, Vmat = Vmat, Vmat2 = Vmat2, 
	#	check = Vmat%*%Vinv, check2 = Vmat2 %*% V2inv, N = sum(S)))
}


###############################################
# Ro of brucellosis in the presence of bTB, age structure
###############################################

Ro_brucellosis_co_age = function(params){
	###################################
	# Input: paramter file
	###################################
	# Get stable age distribution in TB only
	relage = c(0.137, rep(0.368/4, 4), rep(0.185/4, 4),
		rep(0.235/6, 6), rep(0.075/5, 5))
	S0 = 400*relage; It0 = 2*relage; Ib0 = 0*relage; 
	Ic0 = 0*relage; R0 = 0 * relage; Rc0 = 0 * relage
	x0 = c(S0, It0, Ib0, Ic0, R0, Rc0)
	times <- seq(0, 1000, 1)
	sol <- as.data.frame(ode(x0, times, rhs_age_matrix, params))
	S <- unname(unlist( sol[length(sol[,1]), c(2:21)]))
	It <- unname(unlist( sol[length(sol[,1]), c(22:41)]))
	N <- sum(S) + sum(It)
	Sall <- sum(S)
	Itall <- sum(It)
	
	# age specific FOI (a vector)
	params$betaBm <- rep(params$betaB, length.out = 20)
	params$betaBm[2:5] <- exp(0.885) * params$betaB 
	params$betaT <- rep(params$betaT, length.out = 20)

	# Calculate next generation matrix, 1:40 columns, 1:40 rows
	Fmat = rbind(matrix(rep(params$betaBm * S / N, each = 40), nrow = 20, byrow = T),
		matrix(rep(params$rhoB * params$betaBm * It / N, each = 40), nrow = 20, byrow = T))

	# Rows 1:20, cols 1:40
	diag <- params$rhoT * params$betaT * Itall + params$muB + params$gamma + 1
	M1 <- diag(x = diag)
	M1[row(M1) - col(M1) == 1] <- -1
	Vmat1 = cbind(M1, 		
		matrix(0, nrow = 20, ncol = 20) )

	# Rows 21:40, cols 1:40		 
	diag <- - params$rhoT * params$betaT * Itall  
	M4 <- diag(x = diag)
	
	diag <- params$gamma + params$muC + 1
	M5 <- diag(x = diag)
	M5[row(M5) - col(M5) == 1] <- -1
	Vmat2 <- cbind(M4, M5)
	
	Vmat <- rbind(Vmat1, Vmat2)
	Vinv <- solve(Vmat)
	vals <- eigen(Fmat %*% Vinv)$values
	return(max(Re(vals))) 
}

Ro_bTB_single_age(params)  # 3.903531
Ro_bTB_co_age(params)  # 1.9156

Ro_brucellosis_single_age(params)  # 1.165354 
Ro_brucellosis_co_age(params)  #1.584019

###############################################
# Add error bars for Ro
###############################################
par(mfrow = c(1,2))
hist(exp(rnorm(n = 100, mean = 1.1060, sd = 0.3505)), main = mean( exp(rnorm(n = 100, mean = 1.1060, sd = 0.3505)) ) )
hist(rnorm(n = 100, mean = exp(1.1060), sd = exp(0.3505)), main = mean( rnorm(n = 100, mean = exp(1.1060), sd = exp(0.3505)) ))

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

saveRDS(data, file = "~/GitHub/bTB-bruc-co-infection-ms/Ro_confidence_interval_simulation_results.rds")




#######################################################################
#######################################################################
# Add errorbars for endemic equilibrium
#######################################################################
#######################################################################
params <- c(fixed.params, list(gamma=1/2, betaB = 0.6087,
	betaT = 0.0012974, rhoT = 1, rhoB = 2.1, theta= 4, K = 433))
s_index <- 1:20
it_index <- 21:40
ib_index <- 41:60
ic_index <- 61:80
r_index <- 81:100
rc_index <- 101:120
times <- seq(0, 1000, 1)  

source('~/GitHub/bTB-bruc-co-infection-ms/get_EE.R', chdir = TRUE)

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

saveRDS(data2, file = "~/GitHub/bTB-bruc-co-infection-ms/EE_confidence_interval_simulation_results2.rds")

# 7 Errors occur on the lower end of rhoB and upper end of mortB
# these have no brucellosis in single infection case





#######################################################################
#######################################################################
#######################################################################
# Checks for math with no age structure... 
#######################################################################
#######################################################################
#######################################################################

###############################################
# Ro of bTB in the absence of brucellosis, no age structure
###############################################
params = list(theta = 4, K = 400, betaT = 0.0012974, b = 0.5, bT = 0.86, muS = 0.06)
x = c(S = 541)   # S set in age age structured model

Ro_bTB_single = function(params, x){
	###################################
	# Input: 
	# x = c(S = final S value at params)
	# params = parameter values()
	# Output: Ro for bTB alone
	###################################
	with(as.list(c(x, params)), {
		v21 <- (-b*bT - b*bT*((S/K)^theta) + 
			b*S*theta*((S/K)^(theta-1))/K)/(1+((S/K)^theta)) + betaT*S
		v22 <- (-b - b*((S/K)^theta) + b*S*theta*((S/K)^(theta-1))/K)/(1+((S/K)^theta)) - muS
		Fmat <- matrix(c(betaT*S,0,0,0), nrow = 2, byrow = T)
		Vmat <- matrix(c(muT, 0, v21, v22), nrow = 2, byrow = T)
		Vinv <- solve(Vmat)
		vals1 <- eigen(Fmat%*%Vinv)$values #??? matrix algebra?
		vals2 <- eigen(betaT*S/muT)$values
		return(list(Ro_longway = max(vals1), Ro_shortway = vals2, 
			Vmatrix =Vmat, check1 = Vmat%*%Vinv))
		}
	)
}
Ro_bTB_single(params, x)  # 11.5

###############################################
# Ro of bTB in the presence of brucellosis, no age structure
###############################################
noage_params = list(theta = 4, K = 400, betaT = 0.0012, b = 0.5, bT = 0.86, muS = 0.06,
	muC = 0.348, rhoB = 2.1, betaB = 1.004592, epsilon = 0.01, gamma = 0.5, muT = 3*0.06, rhoT = 1, rhoB = 2.1)
x = c(S = 319.6701, Ib = 45.69181, R = 93.95181)

Ro_bTB_co = function(params, x){
	###################################
	# Input: 
	# x = c(S = final S value at params, Ib= final Ib value, R = final value)
	# params = parameter values()
	# Output: Ro for bTB alone
	###################################
	with(as.list(c(x, params)), {
		# Calculate next generation matrix
		#Fmat <- matrix(c(betaT*S, betaT*S, betaT*S, 0, 0, 0, 0, 0, 0), nrow = 3, byrow = TRUE)
		Fmat <- matrix(c(betaT*S, betaT*S, betaT*S, 
			rhoT*betaT*Ib, rhoT*betaT*Ib, rhoT*betaT*Ib, 
			rhoT*betaT*R, rhoT*betaT*R, rhoT*betaT*R), nrow = 3, byrow = TRUE)
		v11 <- ((S + Ib + R) * rhoB * betaB * Ib / ((S+ Ib + R)^2)) + muT
		v21 <- (-(S + Ib + R) * rhoB * betaB * Ib / ((S+ Ib + R)^2)) 
		
		Vmat = matrix(c(v11, 0, 0, v21, muC + gamma, -epsilon, 0, 0, epsilon + muC), nrow = 3, byrow = TRUE)
		Vinv <- solve(Vmat)
		vals <- eigen(Fmat %*% Vinv)$values
		return(max(Re(vals)))
		}
	)
}
Ro_bTB_co(noage_params, x)   # 1.6


###############################################
# Ro of brucellosis in the presence of bTB, no age structure
###############################################
source('~/GitHub/bTB-bruc-co-infection-ms/fixed_parameters_norecovery_agematrix.R', chdir = TRUE)
source('~/GitHub/bTB-bruc-co-infection-ms/rhs_age.R', chdir = TRUE)
params <- c(fixed.params.olddz, list(gamma=1/2, betaB = 1.004592,
	betaT = 12.833531/10000, rhoT = 1, rhoB = 2.1, theta= 4, K = 433))


Ro_brucellosis_co  = function(params) {
	# Generate guess at endemic conditions
	relage = c(0.137, rep(0.368/4, 4), rep(0.185/4, 4),
		rep(0.235/6, 6), rep(0.075/5, 5))
	S0 = 400*relage; It0 = 2*relage; Ib0 = 0*relage; 
	Ic0 = 0*relage; R0 = 0 * relage; Rc0 = 0 * relage
	x0 = c(S0, It0, Ib0, Ic0, R0, Rc0)
	times <- seq(0, 500, 1)
	sol <- as.data.frame(ode(x0, times, rhs_age_matrix, params))
	S <- sum(unname(unlist( sol[500, c(2:21)])))
	It <- sum(unname(unlist( sol[500, c(22:41)])))
	N <- sum(S) + sum(It)
	
	Fmat = matrix(c(params$betaB * S / N, params$betaB * S / N, 
		params$rhoB * params$betaB * It / N, params$rhoB * It / N), byrow = T, nrow = 2)	
	Vmat = matrix(c(params$rhoT * params$betaT * It + params$gamma + mean(params$muB), 0,
		-params$rhoT * params$betaT * It, params$gamma + mean(params$muB)	), byrow = T, nrow = 2)
	Vinv <- solve(Vmat)
	vals <- eigen(Fmat %*% Vinv)$values
	return(max(Re(vals)))
}
Ro_brucellosis_co(params)  # 2.15


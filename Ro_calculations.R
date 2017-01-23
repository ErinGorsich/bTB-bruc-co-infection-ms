library(MASS)
#library(gmp)  # precision on matrix inverse calculation, not yet usnig

###############################################
# Ro of bTB in the absence of brucellosis, no age structure
###############################################
params = list(theta = 4, K = 400, betaT = 0.001283353, b = 0.5, bT = 0.86, muS = 0.06)
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
# Ro of bTB in the absence of brucellosis, with age structure
###############################################
source('~/GitHub/bTB-bruc-co-infection-ms/fixed_parameters_norecovery_agematrix.R', chdir = TRUE)
source('~/GitHub/bTB-bruc-co-infection-ms/rhs_age.R', chdir = TRUE)
params <- c(fixed.params.olddz, list(gamma=1/2, betaB = 1.004592,
	betaT = 12.833531/10000, rhoT = 1, rhoB = 2.1, theta= 4, K = 433))

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
	times <- seq(0, 500, 1)
	sol <- as.data.frame(ode(x0, times, rhs_age_matrix, params))
	S <- unname(unlist( sol[500, c(2:21)]))

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

Ro_bTB_single_age(params)  # 3.431935


###############################################
# Ro of bTB in the presence of brucellosis, no age structure
###############################################

params = list(theta = 4, K = 400, betaT = 0.0012, b = 0.5, bT = 0.86, muS = 0.06,
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
Ro_bTB_co(params, x)   # 1.6



###############################################
# Ro of bTB in the absence of brucellosis, with age structure (and age sp. FOI, brucellosis)
###############################################
source('~/GitHub/bTB-bruc-co-infection-ms/fixed_parameters_norecovery_agematrix.R', chdir = TRUE)
source('~/GitHub/bTB-bruc-co-infection-ms/rhs_age.R', chdir = TRUE)
params <- c(fixed.params.olddz, list(gamma=1/2, betaB = 1.004592,
	betaT = 12.833531/10000, rhoT = 1, rhoB = 2.1, theta= 4, K = 433))

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
	times <- seq(0, 500, 1)
	sol <- as.data.frame(ode(x0, times, rhs_age_matrix, params))
	S <- unname(unlist( sol[500, c(2:21)]))
	Ib <- unname(unlist( sol[500, c(42:61)]))
	R <- unname(unlist( sol[500, c(82:101)]))
	N <- sum(S) + sum(Ib) + sum(R)
	Sall <- sum(S)
	Iball <- sum(Ib)
	Rall <- sum(R)
	
	#  age specific FOI
	params$betaBm <- rep(params$betaB, length.out = 20)
	params$betaBm[2:3] <- exp(0.0885) * params$betaB 

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

Ro_bTB_co_age(params)  # 0.96
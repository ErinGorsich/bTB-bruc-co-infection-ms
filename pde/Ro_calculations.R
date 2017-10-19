# NEEDS UPDATED

Ro_bTB_single = function(params, agemax, agestep){
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

Ro_bTB_co = function(params, agemax, agestep){
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

Ro_brucellosis_single = function(params, agemax, agestep){
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

Ro_brucellosis_co = function(params, agemax, agestep){
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
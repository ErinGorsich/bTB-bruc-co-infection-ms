# WARNING: reads in x0, xB, and xT from global environment
# agemax and N need to be global parameters
# indicies as well

Ro_bTB_single = function(params, x0){
	###################################
	# Input: x0 = c(S = final, 1*N S vector at params)
	# Output: Ro
	###################################
	S <- x0[s.index]
	ageflux <- params$aging[1,1]

	# Calculate next generation matrix
	Fmat = params$betaT* matrix(rep(S, each = N), nrow = N, byrow = T)
	Vmat = diag(x = params$muT - ageflux)
	Vmat[row(Vmat) - col(Vmat) == 1] <- ageflux
	Vinv <- solve(Vmat)
	vals <- eigen(Fmat %*% Vinv)$values
	
	#Vmat <- diag(x = params$muT)
	#Vinv <- solve(Vmat)
	#vals <- eigen(Fmat %*% Vinv)$values
	return(max(Re(vals)))
}

###############################################
# Ro of brucellosis in the absence of bTB, age structure
###############################################

Ro_brucellosis_single = function(params, x0){
	###################################
	# Input: paramter file
	###################################
	# Get stable age distribution in dz free conditions
	S <- x0[s.index]
	S.fd <- S/sum(S)
	ageflux <- params$aging[1,1]

	Fmat <- params$betaB* matrix(rep(S.fd, each = N), nrow = N, byrow = T)
	Fmat[which(ages >= 2 & ages <=5), ] <- exp(0.885) * Fmat[which(ages >= 2 & ages <=5), ]
	
	Vmat <- diag(x = params$muB + params$gamma - ageflux)
	Vmat[row(Vmat) - col(Vmat) == 1] <- ageflux
	Vinv <- solve(Vmat)
	vals <- eigen(Fmat %*% Vinv)$values
	return(max(Re(vals2)))
}


###############################################
# Ro of bTB in the absence of brucellosis, with age structure (and age sp. FOI, brucellosis)
###############################################
Ro_bTB_co = function(params, xB){
	###################################
	# Input: x = c(S = final, 1*20 S vector at params)
	# Output: Ro
	###################################
	# Get stable age distribution in dz free conditions
	S <- xB[s.index]
	Ib <- xB[ib.index]
	R <- xB[r.index]

	Tot <- sum(S) + sum(Ib) + sum(R)
	Sall <- sum(S)
	Iball <- sum(Ib)
	Rall <- sum(R)
	
	ageflux <- params$aging[1,1]
	
	#  age specific FOI (vector)
	params$betaB <- rep(params$betaB, length.out = N)
	params$betaB[which(ages >= 2 & ages <=5)] <- exp(0.885) * params$betaB[which(ages >= 2 & ages <=5)] 

	# Calculate next generation matrix, 1:60 columns, 1:60 rows
	Fmat = rbind(params$betaT * matrix(rep(S, each = 3*N), nrow = N, byrow = T),
		params$rhoT * params$betaT * matrix(rep(Ib, each = 3*N), nrow = N, byrow = T),
		params$rhoT * params$betaT * matrix(rep(R, each = 3*N), nrow = N, byrow = T))
	
	# rows 1:N, cols 1:3*N
	diag <- (params$rhoB * params$betaB * Iball / Tot) + params$muT - ageflux
	M1 <- diag(x = diag)
	M1[row(M1) - col(M1) == 1] <- ageflux
	Vmat1 = cbind(M1, 		
		matrix(0, nrow = N, ncol = 2*N) )

	# rows 21:40, cols 1:60		 
	diag <- - params$rhoB * params$betaB * Iball / Tot 
	M4 <- diag(x = diag)
	
	diag <- params$gamma + params$muC - ageflux
	M5 <- diag(x = diag)
	M5[row(M5) - col(M5) == 1] <- ageflux
	
	M6 <- diag(x = rep(- params$epsilon), N)
	Vmat2 <- cbind(M4, M5, M6)
	
	# rows 41:60, cols 1:60		 
	M7 <- matrix(0, nrow = N, ncol = N)
	M8 <- diag(x = rep(-params$gamma, N))
	M9 <- diag(x = params$epsilon + params$muRC - ageflux)
	M9[row(M9 - col(M9) == 1)] <- ageflux
	Vmat3 <- cbind(M7, M8, M9)
	
	Vmat <- rbind(Vmat1, Vmat2, Vmat3)
	Vinv <- solve(Vmat)
	vals <- eigen(Fmat %*% Vinv)$values
	return(max(Re(vals)))   # 0.9634113
}



###############################################
# Ro of brucellosis in the presence of bTB, age structure
###############################################
Ro_brucellosis_co = function(params, xT){
	###################################
	# Input: paramter file
	###################################
	# Get stable age distribution in TB only
	S <- xT[s.index]
	It <- xT[it.index]

	Tot <- sum(S) + sum(It) 
	Sall <- sum(S)
	Itall <- sum(It)

	ageflux <- params$aging[1,1]

	# age specific FOI (a vector)
	betaBm <- rep(params$betaB, length.out = N)
	betaBm[which(ages >= 2 & ages <=5)] <- exp(0.885) * params$betaB
	betaTm <- rep(params$betaT, length.out = N)

	# Calculate next generation matrix, 1:40 columns, 1:40 rows
	Fmat = rbind(
		matrix(rep(betaBm * S / Tot, each = 2* N), nrow = N, byrow = T),
		matrix(rep(params$rhoB * betaBm * It / Tot, each = 2*N), nrow = N, 
		byrow = T))

	# Rows 1:20, cols 1:40
	diag <- params$rhoT * params$betaT * Itall + params$muB + 
		params$gamma - ageflux
	M1 <- diag(x = diag)
	M1[row(M1) - col(M1) == 1] <- ageflux
	Vmat1 = cbind(M1, 		
		matrix(0, nrow = N, ncol = N) )

	# Rows 21:40, cols 1:40		 
	diag <- - params$rhoT * betaTm * Itall  
	M4 <- diag(x = diag)
	
	diag <- params$gamma + params$muC - ageflux
	M5 <- diag(x = diag)
	M5[row(M5) - col(M5) == 1] <- ageflux
	Vmat2 <- cbind(M4, M5)
	
	Vmat <- rbind(Vmat1, Vmat2)
	Vinv <- solve(Vmat)
	vals <- eigen(Fmat %*% Vinv)$values
	return(max(Re(vals))) 
}
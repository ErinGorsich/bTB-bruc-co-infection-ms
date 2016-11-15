#############################################################
#############################################################
# Write ODEs, density dependence... age structure, gamma ageing & wafwi
#############################################################
#############################################################

# logistic
rhs_age_matrix = function(times, x, params){
	##########################
	# Inputs: t = time sequence; 
	# x = initial conditions, vector(length=24), 4 age categories
	# params= list(
	# b1,b2,b3,b4,b5: age sp. prop reduction in fecundity (20 long)
	# with TB/bruc/chronic bruc/co/co-chronic 
	# betaT, betaB: transmission rates, age independent
	# rhoB & rhoT: prop increase in transmission 
	# with co-infection (scalor=age independent)
	# gamma, epsilon: scalors	 
	# b = a vector (0, small no., birth rate set to get growth = 1.2, same rate)
	# muS, muT, muB, muC, muR, muRC age specific mortality rates.  
	##########################
	# Output: differences, for 20 ages
	##########################
	with(as.list(c(x, params)), {

		# Assign state variables, each 20 long, 
		# 4 categories: 1-3.9, 4-4.9, 5-14.9, 15+)
		# from incidence anlaysis, use [0-3),[3 +) 
		# from mortality analysis, use [0-3),[3 +)
		# from birth analyis, use, [3-5),[5+)
		s_index <- 1:20
		it_index <- 21:40
		ib_index <- 41:60
		ic_index <- 61:80
		r_index <- 81:100
		rc_index <- 101:120

		S = x[s_index] 				
		It = x[it_index] 			
		Ib = x[ib_index]
		Ic = x[ic_index]
		R = x[r_index]
		Rc = x[rc_index]

		# Population size (N)
		Nall <- sum(S + It + Ib + Ic + R + Rc)  # overall
		N <- S + It + Ib +Ic + R + Rc 			# by age category

		# Frequency dependent force of infection is age dependent 
		# turn betaB into a matrix
		betaBm <- matrix(nrow = 20, ncol = 20)
		betaBm[1:20, 1:20] <- betaB
		betaBm[4:6,] <- exp(0.0885) * betaB  #4-6yr have higher suscept.
		betaTm <- matrix(c(rep(betaT, 400)), nrow = 20, ncol = 20)
	
		# lambdaB = age specific vectors!
		lambdaT <- betaTm %*% (It + Ic + Rc) / Nall 
		lambdaB <- betaBm %*% (Ib + Ic) / Nall
		lambdapT <- rhoT * betaTm %*% (It + Ic + Rc) / Nall
		lambdapB <- rhoB * betaBm %*% (Ib + Ic) / Nall

		# Age spcific pop contributing to births (Nb = vector); 
		# reduced births(b1, b2...=vectors; Nb = vector, birth = vector) 
		Nb <- S + b1 * It + b2 * Ib + 
			b3 * R + b4 * Ic + b5 * Rc
		birth <- c(b %*% Nb)
		r <- 5
		logistic_births <- c(birth * ((1- (r/birth)*(Nall/K))), rep(0, 19))
		dS <- logistic_births + aging %*% S
			- (lambdaT + lambdaB) * S - muS * S
		dIt <- lambdaT * S - (lambdapB + muT) * It + aging %*% It
		dIb <- lambdaB * S + aging %*% Ib +
			epsilon * R - (gamma + lambdapT + muB) * Ib
		dIc <- lambdapT * Ib + lambdapB * It + aging %*% Ic +
			epsilon * Rc - (gamma + muC) * Ic
		dR <- gamma * Ib - (epsilon + muR + lambdapT) * R + aging %*% R
		dRc <- lambdapT * R + gamma * Ic + aging %*% Rc -
			(epsilon + muRC) * Rc
		
		out = list(c(dS, dIt, dIb, dIc, dR, dRc))
		return(out)	
		}
	)	
}


# Ricker
rhs_age_matrix_ricker = function(times, x, params){
	##########################
	# Inputs: t = time sequence; 
	# x = initial conditions, vector(length=24), 4 age categories
	# params= list(
	# b1,b2,b3,b4,b5: age sp. prop reduction in fecundity (20 long)
	# with TB/bruc/chronic bruc/co/co-chronic 
	# betaT, betaB: transmission rates, age independent
	# rhoB & rhoT: prop increase in transmission 
	# with co-infection (scalor=age independent)
	# gamma, epsilon: scalors	 
	# b = a vector (0, small no., birth rate set to get growth = 1.2, same rate)
	# muS, muT, muB, muC, muR, muRC age specific mortality rates.  
	##########################
	# Output: differences, for 20 ages
	##########################
	with(as.list(c(x, params)), {

		# Assign state variables, each 20 long, 
		# 4 categories: 1-3.9, 4-4.9, 5-14.9, 15+)
		# from incidence anlaysis, use [0-3),[3 +) 
		# from mortality analysis, use [0-3),[3 +)
		# from birth analyis, use, [3-5),[5+)
		s_index <- 1:20
		it_index <- 21:40
		ib_index <- 41:60
		ic_index <- 61:80
		r_index <- 81:100
		rc_index <- 101:120

		S = x[s_index];  S[S<0] <- 0 				
		It = x[it_index]; It[It<0] <- 0 			
		Ib = x[ib_index]; Ib[Ib <0] <- 0
		Ic = x[ic_index]; Ic[Ic <0] <- 0
		R = x[r_index];   R[R<0] <- 0
		Rc = x[rc_index]; Rc[Rc < 0] <- 0
	
		# Population size (N)
		Nall <- sum(S + It + Ib + Ic + R + Rc)   # overall
		N <- S + It + Ib +Ic + R + Rc 			# by age category

		# Frequency dependent force of infection is age dependent 
		# turn betaB into a matrix
		betaBm <- matrix(nrow = 20, ncol = 20)
		betaBm[1:20, 1:20] <- betaB
		betaBm[4:6,] <- exp(0.0885) * betaB  #4-6yr have higher suscept.
		betaTm <- matrix(c(rep(betaT, 400)), nrow = 20, ncol = 20)
	
		# lambdaB = age specific vectors!
		lambdaT <- betaTm %*% (It + Ic + Rc) #/ Nall
		lambdaB <- betaBm %*% (Ib + Ic) #/ Nall  # extra 100 from freqdep
		lambdapT <- rhoT * betaTm %*% (It + Ic + Rc) #/ Nall
		lambdapB <- rhoB * betaBm %*% (Ib + Ic) #/ Nall

		# Age spcific pop contributing to births (Nb = vector); 
		# reduced births(b1, b2...=vectors; Nb = vector, birth = vector) 
		Nb <- S + b1 * It + b2 * Ib + 
			b3 * R + b4 * Ic + b5 * Rc
		birth <- c(b %*% Nb, rep(0, 19))
#		dS <- birth * (exp(- Nall/800)) + aging %*% S
#			- (lambdaT + lambdaB) * S - muS * S
		dS <- birth  + aging %*% S
			- (lambdaT + lambdaB) * S - muS * S		
		dIt <- lambdaT * S - (lambdapB + muT) * It + aging %*% It 
		dIb <- lambdaB * S + aging %*% Ib +
			epsilon * R - (gamma + lambdapT + muB) * Ib
		dIc <- lambdapT * Ib + lambdapB * It + aging %*% Ic +
			epsilon * Rc - (gamma + muC) * Ic
		dR <- gamma * Ib - (epsilon + muR + lambdapT) * R + aging %*% R
		dRc <- lambdapT * R + gamma * Ic + aging %*% Rc -
			(epsilon + muRC) * Rc
		
		out = list(c(dS, dIt, dIb, dIc, dR, dRc))
		return(out)	
		}
	)	
}


#############################################################
#############################################################
# Write ODEs, NOT density dependence... with age structure, exponential aging
#############################################################
#############################################################
rhs_age = function(times, x, params){
	##########################
	# Inputs: t = time sequence; 
	# x = initial conditions, vector(length=24), 4 age categories
	# params= list(
	# b1,b2,b3,b4,b5: age sp. prop reduction in fecundity
	# with TB/bruc/chronic bruc/co/co-chronic 
	# betaT, betaB: age sp, transmission rates
	# rhoB & rhoT: prop increase in transmission 
	# with co-infection (scalor=age independent)
	# gamma, epsilon: scalors	 
	# b = a vector (0, small no., birth rate set to get growth = 1.2, same rate)
	# muS, muT, muB, muC, muR, muRC age specific mortality rates.  
	# l = aging! = vector 
	##########################
		
		# Assign state variables, each 4 long, 
		# 4 categories 
		#(age= 1-3.9, 4-4.9, 5-14.9, 15+)..subsume calf mortality in births
		# from incidence anlaysis, use [0-3),[3 +) 
		# from mortality analysis, use [0-3),[3 +)
		# from birth analyis, use, [3-5),[5+)
		S = x[1:4] 			# S = x[1:nage]
		It = x[5:8] 			# It = x[(nage+1): (2*nage)]
		Ib = x[9:12]
		Ic = x[13:16]
		R = x[17:20]
		Rc = x[21:24]
							
		# Population size (N)
		Nall <- sum(S + It + Ib + Ic + R + Rc)  # overall
		N <- S + It + Ib +Ic + R + Rc 			# by age category
		
		# Age spcific population contributing to births (Nb = vector); 
		# reduced births with dz (b1, b2...=vectors; Nb = vector, birth = vector) 
		Nb <- S + params$b1 * It + params$b2 * Ib + 
			params$b3 * R + params$b4 * Ic + params$b5 * Rc
		birth <- c(params$b %*% Nb, 0, 0, 0) 
		# Maturation into or out of a stage: 1/stage duration
		# 4 categories (age= 1-3.9, 4-4.9, 5-14.9, 15+)... subsume calf mortality in births...
		lin = c(0, 1/3, 1, 1/10)  # birth in first age only!
		lout = c(1/3, 1, 1/10, 0)
					
		# Frequency dependent force of infection is age dependent 
		# All lambdas = scalors; all betas = age specific vectors!
		lambdaT <- params$betaT %*% (It + Ic + Rc) 
		lambdaB <- params$betaB %*% (Ib + Ic) 
		lambdapT <- params$rhoT * params$betaT %*% (It + Ic + Rc)
		lambdapB <- params$rhoB * params$betaB %*% (Ib + Ic) 
		
		# differential equations
		dx <- vector(length=24)
		for (i in 1:4){
			dx[i] <-  birth[i] + (lambdaT + lambdaB) * S[i] -				# dS
				 params$muS[i] * S[i] + (lin[i] - lout[i]) * S[i]
			dx[i+4] <- lambdaT * S[i] - (lambdapB + params$muT[i]) * It[i]+ # dIt
				(lin[i] - lout[i]) * It[i] 
			dx[i+8] <- lambdaB * S[i] + params$epsilon * R[i] - 				# dIb
				(params$gamma + lambdapT + params$muB[i]) * Ib[i] + 
				(lin[i] - lout[i]) * Ib[i]  
			dx[i+12] <- lambdapT * Ib[i] + lambdapB * It[i] + 				# dIc
				params$epsilon * Rc[i] - 
				(params$gamma + params$muC[i]) * Ic[i] + 
				(lin[i] - lout[i]) * Ic[i]  
			dx[i+16] <-  params$gamma * Ib[i] - 								# dR
				(params$epsilon + params$muR[i] + 			
				lambdapT) * R[i] +			
				(lin[i] - lout[i]) * R[i]
			dx[i+20] <- lambdapT * R[i] + params$gamma * Ic[i] - 			#dRc
				(params$epsilon + params$muRC[i]) * Rc[i] + 
				(lin[i] - lout[i]) * Rc[i]
		}
		return(list(dx))
}



#############################################################
#############################################################
# Write ODEs, density dependence... with age structure (not updated), exponential ageing
#############################################################
#############################################################
rhs_age_logistic_exp = function(times, x, params){
	##########################
	# Inputs: t = time sequence; 
	# x = initial conditions
	# params= list(
	# b1, b2, b3, b4, b5-prop reduction in fecundity with TB, bruc, chronic bruc, co, coinfect/chronic 
	# K 
	# betaT, betapT, betaB, betapB- transmission rates, p fo co-inf first
	# gamma, epsilon	 
	# b = max birth rate in susceptibles(atN=0) = logistic growth, see pg 376 b with beta = 1/K
	# muS, muT, muB, muC, muR, muRC mortality rates.  
	# l = aging! = vector 
	##########################
	# Assign state variables, each 4 long, 
		# 4 categories (age= 1-3.9, 4-4.9, 5-14.9, 15+)... subsume calf mortality in births...
		# from incidence anlaysis, use [0-3),[3 +) 
		# from mortality analysis, use [0-3),[3 +)
		# from birth analyis, use, [3-5),[5+)
		S = x[1:4] 			# S = x[1:nage]
		It = x[5:8] 			# It = x[(nage+1): (2*nage)]
		Ib = x[9:12]
		Ic = x[13:16]
		R = x[17:20]
		Rc = x[21:24]
							
		# Population size (N)
		Nall <- sum(S + It + Ib + Ic + R + Rc)  # overall
		N <- S + It + Ib +Ic + R + Rc 			# by age category
		
		# Stage durations
		lin = c(0, 1/3, 1, 1/10)  # birth in first age only!
		lout = c(1/3, 1, 1/10, 0)
		
		# Age spcific population contributing to births (Nb = vector); 
		# account for reduced births with infection (b2, b2... vectors; Nb scalor; 
		Nb <- S + params$b1 * It + params$b2 * Ib + 
			params$b3 * R + params$b4 * Ic + params$b5 * Rc
		birth_noDD <- c(params$b %*% Nb, 0, 0, 0)
		birth <- birth_noDD * (1 - (Nall/K))				  # 
					
		# Frequency dependent force of infection is age dependent 
		# All lambdas = scalors; all betas = age specific vectors!
		lambdaT <- params$betaT %*% (It + Ic + Rc) 
		lambdaB <- params$betaB %*% (Ib + Ic) 
		lambdapT <- params$rhoT * params$betaT %*% (It + Ic + Rc)
		lambdapB <- params$rhoB * params$betaB %*% (Ib + Ic) 
		
		# differential equations
		dx <- vector(length=24)
		for (i in 1:4){
			dx[i] <-  birth[i] + (lambdaT + lambdaB) * S[i] -					# dS
				 params$muS[i] * S[i] + (lin[i] - lout[i]) * S[i]
			dx[i+4] <- lambdaT * S[i] - (lambdapB + params$muT[i]) * It[i] +   	# dIt
				(lin[i] - lout[i]) * It[i] 
			dx[i+8] <- lambdaB * S[i] + params$epsilon * R[i] - 				# dIb
				(params$gamma + lambdapT + params$muB[i]) * Ib[i] + 
				(lin[i] - lout[i]) * Ib[i]  
			dx[i+12] <- lambdapT * Ib[i] + lambdapB * It[i] + 					# dIc
				params$epsilon * Rc[i] - 
				(params$gamma + params$muC[i]) * Ic[i] + 
				(lin[i] - lout[i]) * Ic[i]  
			dx[i+16] <-  params$gamma * Ib[i] - 								# dR
				(params$epsilon + params$muR[i] + 			
				lambdapT) * R[i] +			
				(lin[i] - lout[i]) * R[i]
			dx[i+20] <- lambdapT * R[i] + params$gamma * Ic[i] - 				# dRc
				(params$epsilon + params$muRC[i]) * Rc[i] + 
				(lin[i] - lout[i]) * Rc[i]
		}
	return(list(dx))
}

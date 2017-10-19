#############################################################
#############################################################
# Write ODEs, density dependence... age structure
#############################################################
#############################################################

#############################################################
# Beverton & Holt, in main text 
#############################################################
rhs = function(times, x, params){
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
		# 4 categories: 1-2.9, 3-4.9, 5-15.9, 16+)
		# from incidence anlaysis, use [0-4.9),[5 +) 
		# from mortality analysis, use [1-2.9),[3 +)
		# from birth analyis, use, 5+
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
		betaBv <- rep(betaB, 20)
		betaBv[2:5] <- exp(0.885) * betaBv[2:5]   #2-3yr have higher suscept.
		betaTv <- rep(betaT, 20)
	
		# lambdaB = age specific vectors!
		# bTB = density dependent, bruc = freq dependent
		lambdaT <- betaTv * sum((It + Ic + Rc))
		lambdaB <- betaBv * (sum(Ib) + sum(Ic)) / Nall
		lambdapT <- rhoT * betaTv * sum((It + Ic + Rc))
		lambdapB <- rhoB * betaBv * (sum(Ib) + sum(Ic)) / Nall

		# Age spcific pop contributing to births (Nb = vector); 
		# reduced births(b1, b2...=vectors; Nb = vector, birth = vector) 
		Nb <- S + b1 * It + b2 * Ib + 
			b3 * R + b4 * Ic + b5 * Rc
			
		birth <- c(b %*% Nb)
		recruitment <- c(birth / ( 1 + (Nall/K)^theta), rep(0, 19))
		dS <- recruitment - (lambdaT + lambdaB) * S - muS * S
		dIt <- lambdaT * S - (lambdapB + muT) * It
		dIb <- lambdaB * S +	epsilon * R - (gamma + lambdapT + muB) * Ib
		dIc <- lambdapT*Ib + lambdapB*It + epsilon * Rc - (gamma + muC)*Ic
		dR <- gamma * Ib - (epsilon + muR + lambdapT) * R
		dRc <- lambdapT * R + gamma * Ic - (epsilon + muRC) * Rc
		out = list(c(dS, dIt, dIb, dIc, dR, dRc))
		return(out)	
		}
	)
}

birthplay = function(r, Nall){
	0.5*800- r * 800 * (Nall/1000)
}


#############################################################
# Logistic growth - no aging
#############################################################
rhs_logistic = function(times, x, params){
	##########################
	# Inputs: t = time sequence; 
	# x = initial conditions, vector(length=24), 4 age categories
	# params= list(
	# b1,b2,b3,b4,b5: age sp. prop reduction in fecundity (20 long)
	# with TB/bruc/chronic bruc/co/co-chronic 
	# betaT, betaB: transmission rates, age independent
	# rhoB & rhoT: prop increase in transmission (age independent)
	# gamma, epsilon: scalors	 
	# b = a vector (0, small no., birth rate set to get growth = 1.2, same rate)
	# muS, muT, muB, muC, muR, muRC age specific mortality rates.  
	##########################
	# Output: differences, for 20 ages
	##########################
	with(as.list(c(x, params)), {
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
		betaBv <- rep(betaB, 20)
		betaBv[2:5] <- exp(0.885) * betaBv[2:5]   #2-3yr have higher suscept.
		betaTv <- rep(betaT, 20)
	
		# lambdaB = age specific vectors!
		# bTB = density dependent, bruc = freq dependent
		lambdaT <- betaTv * (It + Ic + Rc)
		lambdaB <- betaBv * (Ib + Ic) / Nall
		lambdapT <- rhoT * betaTv * (It + Ic + Rc)
		lambdapB <- rhoB * betaBv * (Ib + Ic) / Nall

		# Age spcific pop contributing to births (Nb = vector); 
		# reduced births(b1, b2...=vectors; Nb = vector, birth = vector) 
		Nb <- S + b1 * It + b2 * Ib + 
		b3 * R + b4 * Ic + b5 * Rc
			
		birth <- c(b %*% Nb) # scalar
		recruitment <- c(birth * ( 1 - (1/b[7]) * (Nall/K)), rep(0, 19))
		dS <- recruitment - (lambdaT + lambdaB) * S - muS * S
		dIt <- lambdaT * S - (lambdapB + muT) * It 
		dIb <- lambdaB * S +	epsilon * R - (gamma + lambdapT + muB) * Ib
		dIc <- lambdapT*Ib + lambdapB*It + epsilon * Rc - (gamma + muC)*Ic
		dR <- gamma * Ib - (epsilon + muR + lambdapT) * R
		dRc <- lambdapT * R + gamma * Ic - (epsilon + muRC) * Rc
		out = list(c(dS, dIt, dIb, dIc, dR, dRc))
		return(out)	
		}
	)
}



#############################################################
# Ricker
#############################################################
rhs_ricker = function(times, x, params){
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
		betaBv <- rep(betaB, 20)
		betaBv[2:5] <- exp(0.885) * betaBv[2:5]   #2-3yr have higher suscept.
		betaTv <- rep(betaT, 20)
	
		# lambdaB = age specific vectors!
		# bTB = density dependent, bruc = freq dependent
		lambdaT <- betaTv * (It + Ic + Rc)
		lambdaB <- betaBv * (Ib + Ic) / Nall
		lambdapT <- rhoT * betaTv * (It + Ic + Rc)
		lambdapB <- rhoB * betaBv * (Ib + Ic) / Nall

		# Age spcific pop contributing to births (Nb = vector); 
		# reduced births(b1, b2...=vectors; Nb = vector, birth = vector) 
		Nb <- S + b1 * It + b2 * Ib +  b3 * R + b4 * Ic + b5 * Rc
		birth <- c(b %*% Nb, rep(0, 19))
		dS <- birth * (exp(- Nall/K)) - (lambdaT + lambdaB) * S - muS * S
		dIt <- lambdaT * S - (lambdapB + muT) * It
		dIb <- lambdaB * S + 	epsilon * R - (gamma + lambdapT + muB) * Ib
		dIc <- lambdapT * Ib + lambdapB * It + epsilon * Rc - (gamma + muC) * Ic
		dR <- gamma * Ib - (epsilon + muR + lambdapT) * R
		dRc <- lambdapT * R + gamma * Ic + (epsilon + muRC) * Rc
		
		out = list(c(dS, dIt, dIb, dIc, dR, dRc))
		return(out)	
		}
	)
}
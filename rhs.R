#############################################################
#############################################################
# Write ODEs, density dependence... no age structure
#############################################################
#############################################################
rhs = function(times, x, params){
	##########################
	# Inputs: t = time sequence; 
	# x = initial conditions
	# params= list(
	# b1, b2, b3, b4, b5-prop reduction in fecundity with TB, bruc, chronic bruc, co, coinfect/chronic 
	# K 
	# betaT, betapT, betaB, betapB- transmission rates, p fo co-inf first
	# gamma, epsilon	 
	# b = max birth rate in susceptibles (when N = 0) 
	# muS, muT, muB, muC, muR, muRC mortality rates.  
	##########################
	with(as.list(c(x, params)), {
					
		# Overall population size (N)
		N <- S + It + Ib + Ic + R + Rc
		
		# Population contributing to births (Nb); 
		# account for reduced births with infection
		Nb <- S + b1 * It + b2 * Ib + b3 *R + b4 * Ic + b5 * Rc
					
		# Frequency dependent force of infection is independent of age
		lambdaT <- betaT * (It + Ic + Rc) 
		lambdaB <- betaB * (Ib + Ic) 
		lambdapT <- betapT * (It + Ic + Rc)
		lambdapB <- betapB * (Ib + Ic) 
		
		# differential equations
		# assume mortality in chronics is same as active infection
		dS <- b * Nb * (1 - (r/b) * (N/K) ) - (lambdaT + lambdaB) * S - muS * S 
		dIt <- lambdaT * S - (lambdapB + muT) * It 
		dIb <- lambdaB * S + epsilon * R - (gamma + lambdapT + muB) * Ib
		dIc <- lambdapT * Ib + lambdapB * It + epsilon * Rc - (gamma + muC) * Ic
		dR <-  gamma * Ib - (epsilon + muR + lambdapT) * R
		dRc <- lambdapT * R + gamma * Ic - (epsilon + muRC) * Rc
		
		out = list(c(dS, dIt, dIb, dIc, dR, dRc))
		return(out)
		}
	)
}
#############################################################
#############################################################
# Code defining the three functional forms of the model: 
# rhs = Beverton & Holt form of DD in births
# rhs_logistic = Logistic form of DD in births
# rhs_ricker = Ricker form of DD in births

#############################################################
# Beverton & Holt descretized PDE
#############################################################
rhs = function(times, x, params){
	##########################
	# co-infection model pde - 
	# x = 6*length(ages) long initial conditions
	# params list containing features of aging
	##########################
	with(as.list(c(x, params)), {
		n.ages <- length(ages) # total number of bins

		# Define states
		S = x[1:n.ages] 				
		It = x[seq(n.ages + 1, 2 * n.ages)] 			
		Ib = x[seq(2 * n.ages + 1, 3 * n.ages)]
		Ic = x[seq(3 * n.ages + 1, 4 * n.ages)]
		R = x[seq(4 * n.ages + 1, 5 * n.ages)]
		Rc = x[seq(5 * n.ages + 1, 6 * n.ages)]
		
		# Population size (N)
		Nall <- sum(S + It + Ib + Ic + R + Rc)  # overall
		N <- S + It + Ib +Ic + R + Rc 			# by age category

		# Frequency dependent force of infection is age dependent 
		betaBm <- matrix(nrow = n.ages, ncol = n.ages)
		betaBm[1:n.ages, 1:n.ages] <- betaB
		dims <- which(ages >= 2 & ages <=5)
		betaBm[dims,] <- exp(0.885) * betaB
		betaTm <- matrix(nrow = n.ages, ncol = n.ages)
		betaTm[1:n.ages, 1:n.ages] <- betaT
			
		lambdaT <- betaTm %*% (It + Ic + Rc)
		lambdaB <- betaBm %*% (Ib + Ic) / Nall
		lambdapT <- rhoT * betaTm %*% (It + Ic + Rc)
		lambdapB <- rhoB * betaBm %*% (Ib + Ic) / Nall

		Nb <- N
			
		birth <- c(b %*% Nb)
		recruitment <- c(birth / ( 1 + (Nall/K)^theta), rep(0, times = n.ages - 1))
		dS <- recruitment + aging %*% S - (lambdaT + lambdaB) * S - muS * S
		dIt <- lambdaT * S - (lambdapB + muT) * It + aging %*% It
		dIb <- lambdaB * S + aging %*% Ib +	epsilon * R - (gamma + lambdapT + muB) * Ib
		dIc <- lambdapT*Ib + lambdapB*It + aging %*% Ic + epsilon * Rc - (gamma + muC)*Ic
		dR <- gamma * Ib - (epsilon + muR + lambdapT) * R + aging %*% R
		dRc <- lambdapT * R + gamma * Ic + aging %*% Rc - (epsilon + muRC) * Rc
		out = list(c(dS, dIt, dIb, dIc, dR, dRc))
		return(out)	
		}
	)	
}

#############################################################
# Logistic growth (needs param K- set equilibrium population size with no dz)
#############################################################
rhs_logistic = function(times, x, params){
	##########################
	with(as.list(c(x, params)), {
		n.ages <- length(ages) # total number of bins

		# Define states
		S = x[1:n.ages] 				
		It = x[seq(n.ages + 1, 2 * n.ages)] 			
		Ib = x[seq(2 * n.ages + 1, 3 * n.ages)]
		Ic = x[seq(3 * n.ages + 1, 4 * n.ages)]
		R = x[seq(4 * n.ages + 1, 5 * n.ages)]
		Rc = x[seq(5 * n.ages + 1, 6 * n.ages)]
		
		# Population size (N)
		Nall <- sum(S + It + Ib + Ic + R + Rc)  # overall
		N <- S + It + Ib +Ic + R + Rc 			# by age category

		# Frequency dependent force of infection is age dependent 
		# turn betaB into a matrix
		betaBm <- matrix(nrow = n.ages, ncol = n.ages)
		betaBm[1:n.ages, 1:n.ages] <- betaB
		dims <- which(ages >= 2 & ages <=5)
		betaBm[dims,] <- exp(0.885) * betaB
		betaTm <- matrix(nrow = n.ages, ncol = n.ages)
		betaTm[1:n.ages, 1:n.ages] <- betaT
						
		lambdaT <- betaTm %*% (It + Ic + Rc)
		lambdaB <- betaBm %*% (Ib + Ic) / Nall
		lambdapT <- rhoT * betaTm %*% (It + Ic + Rc)
		lambdapB <- rhoB * betaBm %*% (Ib + Ic) / Nall

		Nb <- N
			
		birth <- c(b %*% Nb)
		recruitment <- c(birth * ( 1 - (1/max(b)) * (Nall/K)), rep(0, times = n.ages - 1))
		dS <- recruitment + aging %*% S - (lambdaT + lambdaB) * S - muS * S
		dIt <- lambdaT * S - (lambdapB + muT) * It + aging %*% It
		dIb <- lambdaB * S + aging %*% Ib +	epsilon * R - (gamma + lambdapT + muB) * Ib
		dIc <- lambdapT*Ib + lambdapB*It + aging %*% Ic + epsilon * Rc - (gamma + muC)*Ic
		dR <- gamma * Ib - (epsilon + muR + lambdapT) * R + aging %*% R
		dRc <- lambdapT * R + gamma * Ic + aging %*% Rc - (epsilon + muRC) * Rc
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
	# co-infection model pde - 
	# x = 6*length(ages) long initial conditions
	# params list containing features of aging
	##########################
	with(as.list(c(x, params)), {
		n.ages <- length(ages) # total number of bins

		# Define states
		S = x[1:n.ages] 				
		It = x[seq(n.ages + 1, 2 * n.ages)] 			
		Ib = x[seq(2 * n.ages + 1, 3 * n.ages)]
		Ic = x[seq(3 * n.ages + 1, 4 * n.ages)]
		R = x[seq(4 * n.ages + 1, 5 * n.ages)]
		Rc = x[seq(5 * n.ages + 1, 6 * n.ages)]
		
		# Population size (N)
		Nall <- sum(S + It + Ib + Ic + R + Rc)  # overall
		N <- S + It + Ib +Ic + R + Rc 			# by age category

		# Frequency dependent force of infection is age dependent 
		# turn betaB into a matrix
		betaBm <- matrix(nrow = n.ages, ncol = n.ages)
		betaBm[1:n.ages, 1:n.ages] <- betaB
		dims <- which(ages >= 2 & ages <=5)
		betaBm[dims,] <- exp(0.885) * betaB
		betaTm <- matrix(nrow = n.ages, ncol = n.ages)
		betaTm[1:n.ages, 1:n.ages] <- betaT
			
		lambdaT <- betaTm %*% (It + Ic + Rc)
		lambdaB <- betaBm %*% (Ib + Ic) / Nall
		lambdapT <- rhoT * betaTm %*% (It + Ic + Rc)
		lambdapB <- rhoB * betaBm %*% (Ib + Ic) / Nall

		Nb <- N
			
		birth <- c(b %*% Nb)
		recruitment <- c(birth * exp(- Nall / K) , rep(0, times = n.ages - 1))
		dS <- recruitment + aging %*% S - (lambdaT + lambdaB) * S - muS * S
		dIt <- lambdaT * S - (lambdapB + muT) * It + aging %*% It
		dIb <- lambdaB * S + aging %*% Ib +	epsilon * R - (gamma + lambdapT + muB) * Ib
		dIc <- lambdapT*Ib + lambdapB*It + aging %*% Ic + epsilon * Rc - (gamma + muC)*Ic
		dR <- gamma * Ib - (epsilon + muR + lambdapT) * R + aging %*% R
		dRc <- lambdapT * R + gamma * Ic + aging %*% Rc - (epsilon + muRC) * Rc
		out = list(c(dS, dIt, dIb, dIc, dR, dRc))
		return(out)	
		}
	)	
}







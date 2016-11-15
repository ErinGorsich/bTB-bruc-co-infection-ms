#############################################################
#############################################################
# Fixed parameters, Animals in recoverd class have similar birth/death rates to infecteds 
# No age structure model 
#############################################################
#############################################################
get_new_prop_birth = function(logOR, p1){
	OR <- exp(logOR)
	x <- OR/(p1/(1-p1))
	b <- x/(p1*(1-x))
	return(b)
}

K = 1000

#############################################################
# Mortality, susceptible females
#############################################################
muS <- 1 - 0.96
muB <- 1 - 0.89		# = 3.02 * muS
muT <- 1 - 0.899  	# = 2.82 * muS
muC <- 1 - 0.724		# = (3.02 + 2.82) * muS
muRC <- muC
muR <- muB


#############################################################
# births
#############################################################
b <- 0.41/2 # proportion in LS right before calving that were pregs (or had milk-> assume, can check this)- divide by two for females
# In interdrought periods, buffalo populations grew at rates from 5-15%. 
r <- b - muS

# these represent prob_birth_Infected/probability_birth_Susceptible
b1 <- 0.65 # proportional reduction with bTB
b2 <- 0.68   # proportional reduction with bruc
b3 <- b2 # ASSUME animals with chronic brucellosis have same fecundity as active brucellosis 
# FIX ME- SET TO NO CHANGE
b4 <- 0.8 #get_new_prop_birth(-0.839, b) # proportional reduction if co; -1.619 - 1.50 + 2.28
b5 <- b4 # 

#############################################################
# transmission and recovery rates to vary later. 
#############################################################
epsilon = 0.01


fixed.params = c(b1= b1, b2 = b2, b3 = b3, b4= b4, b5 = b5, 
	b = b, r = r, K = K,
	muS = muS, muB = muB, muT = muT, muC = muC, muR = muR, muRC = muRC, 
	epsilon = epsilon)
# missing: gamma, betaT = betaT, betaB = betaB, betapT = betapT, betapB = betapB

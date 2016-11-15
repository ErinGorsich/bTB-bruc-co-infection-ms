#############################################################
#############################################################
# Fixed parameters, Animals in recoverd class have similar birth/death rates to uninfecteds 
# Model does not include age structure
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
muB <- 1 - 0.89
muT <- 1 - 0.899
muC <- 1 - 0.724
muRC <- muT
muR <- muS
#############################################################
# births
#############################################################
b <- 0.41/2 # proportion in LS right before calving that were pregs (or had milk-> assume, can check this)
# In interdrought periods, buffalo populations grew at rates from 5-15%. 
r <- b - muS

# these represent prob_birth_Infected/probability_birth_Susceptible
b1 <- 0.65 #get_new_prop_birth(-1.619, b) # proportional reduction with bTB
b2 <- 0.68 #get_new_prop_birth(-1.5, b)   # proportional reduction with bruc
b3 <- 1 # ASSUME animals with chronic brucellosis have same fecundity as susceptibles
b4 <- 0.8 
# get_new_prop_birth((-1.619 - 1.50 + 2.28), b) # proportional reduction if co-infected
b5 <- b1 # 

#############################################################
# transmission and recovery rates to vary later. 
#############################################################
epsilon = 0.01


fixed.params.recov = c(b1= b1, b2 = b2, b3 = b3, b4= b4, b5 = b5, 
	b = b, r = r, K = K,
	muS = muS, muB = muB, muT = muT, muC = muC, muR = muR, muRC = muRC, 
	epsilon = epsilon)
# missing: gamma, betaT = betaT, betaB = betaB, betapT = betapT, betapB = betapB

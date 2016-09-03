#############################################################
#############################################################
# Fixed parameters, Animals in recoverd class have similar birth/death rates to infecteds  
#############################################################
#############################################################
get_new_prop_birth = function(logOR, p1){
	OR <- exp(logOR)
	x <- OR/(p1/(1-p1))
	b <- x/(p1*(1-x))
	return(b)
}

# age structure information, used to calculate mortality rates in susceptibles. 
relageall = c(0.137, rep(0.368/4, 4), rep(0.185/4, 4),  # Jolles 2005, set max age at 18
	rep(0.235/6, 6), rep(0.075/3, 3))					# Also in Caron et al. from 2001 KNP
relage = c(relageall[1], sum(relageall[2:3]), relageall[4],  
	relageall[5], sum(relageall[6:length(relageall)]) )  # sums to 1
K = 1000


#############################################################
# Mortality, susceptible females
#############################################################
muSa <- NA; muTa <- NA; muBa <- NA; muCa <- NA
muSa[1] <- 1- 0.7 # mortality rate in calves age [0-1) (1/yr)
muSa[2] <- 1- 0.884 # mortality rate in yearlings [1-3) (1/yr)
muSa[3] <- 1- 0.884 # mortality rate in juveniles [3-4) (1/yr)
muSa[4] <-  1- 0.963 # mortality rate in sub-adults [4-5)  (1/yr)
muSa[5] <-  1- 0.963 # mortality rate in adults 5+  (1/yr)
muS <- sum(muSa* relage) 

# mortality, TB, and Brucellosis positive animals
muT <- 2.82 * muS
muB <- 3.02 * muS
muC <- (2.82 + 3.02) * muS
muRC <- muC
muR <- muB

#############################################################
# births
#############################################################
b <- 0.41 # proportion in LS right before calving that were pregs (or had milk-> assume, can check this)
# In interdrought periods, buffalo populations grew at rates from 5-15%. 
r <- b - muS

# these represent prob_birth_Infected/probability_birth_Susceptible
b1 <- get_new_prop_birth(-1.619, b) # proportional reduction with bTB
b2 <- get_new_prop_birth(-1.5, b)   # proportional reduction with bruc
b3 <- b2 # ASSUME animals with chronic brucellosis have same fecundity as active brucellosis 
# FIX ME- SET TO NO CHANGE
b4 <- 1 #get_new_prop_birth(-0.839, b) # proportional reduction if co; -1.619 - 1.50 + 2.28
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

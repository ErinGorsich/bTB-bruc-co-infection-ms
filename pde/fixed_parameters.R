#############################################################
#############################################################
# Fixed parameters
#############################################################
#############################################################
# Continuous aging, 
# but functionally age categories but we represent as:
# 1-2; 3-4; 5-15; 16 +

# Mortality
#############################################################
muS <- NA
muS <- c(1- 0.9,  1- 0.94, 1- 0.94, 1- 0.9) 
muS <- muS
dmuT <- 2.8
dmuB <- 3.03
dmuC <- 8.56 

# Births
#############################################################
# now b is the maximum possible birth rate in S for age >=5
b <- 0.5
b1 <- 1; b2 <- 1; b3 <- 1; b4 <- 1; b5 <- 1

# Density dependence
#############################################################
theta= 4
K = 433

# Disease
#############################################################
epsilon <- 0.03

# fixed.params
#############################################################
p <- list(b = b, muS = muS, dmuB = dmuB, dmuT = dmuT, 
	dmuC =dmuC, epsilon = epsilon, K = K, theta = theta)


gen_fixed_params = function(agemax, agestep, p, recovery = FALSE) {
	N <- agemax / agestep
	ages <- seq(1, agemax + 1, by = agestep)[-(N+1)]
	juv.index <- which(ages <= 2)
	adult.index <- which(ages > 2 & ages <= 16)
	sens.index <- which(ages > 16)
	fecund.index <- which(ages >= 5)	
	
	# define age fluxes used in rhs
	da <- diff(c(min(ages) - agestep, ages))
	aging <- diag(-1/da)
	aging[row(aging) - col(aging) == 1] <- 1 / head(da, -1)
	
	#Mortality vector
	muS <- NA; muT <- NA; muB <- NA; muC <- NA
	muS[juv.index] <- p$muS[1]
	muS[adult.index] <- p$muS[2]
	muS[sens.index] <- p$muS[4]
	muT <- p$dmuT * muS
	muB <- p$dmuB * muS
	muC <- p$dmuC * muS
	if (recovery == FALSE) {
		muR <- muB
		muRC <- muC
	}
	if (recovery == TRUE) {
		muR <- muS
		muRC <- muT
	}
	
	# Fecundity
	b <- rep(0, length(ages))
	b[fecund.index] <- p$b
	
	fixed.params = list(b = b, aging = aging, ages = ages, 
		muS = muS, muB = muB, muT = muT, muC = muC, muR = muR, 
		muRC = muRC, epsilon = epsilon, K = K, theta = theta)
}

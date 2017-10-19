#############################################################
#############################################################
# Fixed parameters
# Animals in recoverd class have similar birth/death rates to infecteds 
# Units in 1/days
#############################################################
#############################################################
# Age categories but we represent with each age
# 1-2; 3-4; 5-15; 16 +

#############################################################
# Mortality, susceptible females
#############################################################
muS <- NA; muT <- NA; muB <- NA; muC <- NA
muS[1:2] <- 1- 0.86  # annual mortality rate in yearlings 
muS[3:16] <-  1- 0.94 # annual mortality rate in adults 
muS[17:20] <-  1- 0.86 # annual mortality rate in adults 15+  (1/yr)
muS <- muS/365.25
muT <- 2.8 * muS
muB <- 3.03 * muS
muC <- 8.56 * muS
muRC <- muC
muR <- muB

# if back to categories:
#muS <- c(1- 0.86,  1- 0.94, 1- 0.94, 1- 0.86)
#muC[muC > 0.002737] <- 0.002737


#############################################################
# births
#############################################################
# now b is the maximum possible birth rate in S
b <- c(0,0,0,0, rep(0.5/365.25, 16))
#b <- c(0, 0, 0.6/365.25, 0.6/365.25)
b1 <- 1; b2 <- 1; b3 <- 1; b4 <- 1; b5 <- 1

#############################################################
# Density dependence
#############################################################
theta= 4
K = 433

epsilon <- 0.03/ 365.25

#############################################################
fixed.params = list( 
	b1= b1, b2 = b2, b3 = b3, b4= b4, b5 = b5, b = b, 
	muS = muS, muB = muB, muT = muT, muC = muC, 
	muR = muR, muRC = muRC, epsilon = epsilon, 
	K = K, theta = theta)


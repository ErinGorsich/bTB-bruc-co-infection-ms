#############################################################
#############################################################
# Fixed parameters
# Animals in recoverd class have similar birth/death rates to infecteds 
# Units in years
#############################################################
#############################################################

#############################################################
# Aging
#############################################################
ages <- c(seq(1, 20, by = 1)) # upper end of age classes
da <- diff(c(0, ages))
aging <- diag(-1/da)
aging[row(aging) - col(aging) == 1] <- 1/head(da, -1)

#dim(betaB)
#filled.contour(betaBm,plot.title=title(main="WAIFW matrix"))
#############################################################
# Mortality, susceptible females
#############################################################
muS <- NA; muT <- NA; muB <- NA; muC <- NA
muS[1:2] <- 1- 0.86 # mortality rate in yearlings 
muS[3:16] <-  1- 0.94 # mortality rate in adults 
muS[17:20] <-  1- 0.86 # mortality rate in adults 15+  (1/yr)
muT <- 2.8 * muS
muB <- 3.03 * muS
muC <- 8.56 * muS
muC[muC > 1] <- 1
muRC <- muC
muR <- muB

#############################################################
# births
#############################################################
# NOTES: data informing births are from ages 5-10.  Age 4 is the youngest sucessful mom with calf.
# birth rate in uninfected buffalo of each age category
# now b is the maximum possible birth rate 

# Birth rate in uninfecteds
b <- c(0,0,0,0, rep(0.74, 16))
b1 <- 1; b2 <- 1; b3 <- 1; b4 <- 1; b5 <- 1
epsilon = 0.03

fixed.params = list(aging = aging, 
	b1= b1, b2 = b2, b3 = b3, b4= b4, b5 = b5, b = b, 
	muS = muS, muB = muB, muT = muT, muC = muC, 
	muR = muR, muRC = muRC, 	epsilon = epsilon)

# missing: gamma, betaT = betaT, betaB = betaB, betapT = betapT, betapB = betapB



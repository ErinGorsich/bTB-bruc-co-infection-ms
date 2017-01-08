#############################################################
#############################################################
# Fixed parameters
# Animals in recoverd class have similar birth/death rates to infecteds 
# Model contains 4 age classes, 20*6 compartments
#############################################################
#############################################################

#############################################################
# Aging!!!!
#############################################################
ages <- c(seq(1,20,by=1)) # upper end of age classes
da <- diff(c(0,ages))
aging <- diag(-1/da)
aging[row(aging)- col(aging)==1] <- 1/head(da, -1)

#dim(betaB)
#filled.contour(betaBm,plot.title=title(main="WAIFW matrix"))
#############################################################
# Mortality, susceptible females
#############################################################
muS <- NA; muT <- NA; muB <- NA; muC <- NA
muS[1:2]<- 1- 0.86 # mortality rate in yearlings 
muS[3:15]<-  1- 0.94 # mortality rate in adults 
muS[16:20]<-  1- 0.86 # mortality rate in adults 15+  (1/yr)

muT <- 3.12 * muS
muB <- 0
muB[1:2] <- 4.9 * muS[1:2]
muB[3:6] <- muS[3:6]
muB[7:20] <- 4.5 * muS[7:20]

muC <- 0
muC[1:2] <- (3.12 + 4.9) * muS[1:2]
muC[3:6] <- (3.12) * muS[3:6]
muC[7:20] <- (3.12 + 4.5) * muS[7:20]
muRC <- muC
muR <- muB

#############################################################
# births
#############################################################
# NOTES: data informing birts are from ages 4-10.  Age 4 is the youngest sucessful mom with calf.
# birth rate in uninfected buffalo of each age category
# now b is the maximum possible birth rate 
b <- NA; b1 <- NA; b2<- NA; b3 <- NA; b4 <- NA; b5 <- NA
# Birth rate in uninfecteds
b[1:3] <- 0						# [1-4)
b[5:15] <- 0.5 #0.56/2			# age [5-15); raw data (14/25)
b[4] <- 0.45						# age [4-5); raw data (2/26)
b[16:20] <- max(b[5] * 0.7, 0)	# age 15 +; raw data (14/25)
#b <- b/2						
# Proportional reductions
b1 <- c(rep(1, 4), rep(0.86, 20-4)) # c(rep(1, 4), rep(0.65, 20-4))		# with bTB (from raw data)
b2 <- c(1, 1, 1, 1, rep(0.8, 20-4))		# brucellosis
b3<- b2  								# chronic/recovered = active
b4<- rep(1, 20) 		# coinfected
b5 <- b4 								# chronic-coinf = coinfected

epsilon = 0.01

fixed.params = list(aging = aging, 
	b1= b1, b2 = b2, b3 = b3, b4= b4, b5 = b5, b = b, 
	muS = muS, muB = muB, muT = muT, muC = muC, 
	muR = muR, muRC = muRC, 	epsilon = epsilon)
# missing: gamma, betaT = betaT, betaB = betaB, betapT = betapT, betapB = betapB




################################################################
# TEST- DELETE LATER (a little lower)
################################################################
# same rates, higher baseline
muS <- NA; muT <- NA; muB <- NA; muC <- NA
muS[1:2]<- 1- 0.8 # mortality rate in yearlings 
muS[3:15]<-  1- 0.9 # mortality rate in adults 
muS[16:20]<-  1- 0.8 # mortality rate in adults 15+  (1/yr)

muT <- 3.12 * muS
muB <- 0
muB[1:2] <- 4.9 * muS[1:2]
muB[3:6] <- muS[3:6]
muB[7:20] <- 4.5 * muS[7:20]

muC <- 0
muC[1:2] <- (3.12 + 4.9) * muS[1:2]
muC[3:6] <- (3.12) * muS[3:6]
muC[7:20] <- (3.12 + 4.5) * muS[7:20]
muRC <- muC
muR <- muB

test.fixed.params = list(aging = aging, 
	b1= b1, b2 = b2, b3 = b3, b4= b4, b5 = b5, b = b, 
	muS = muS, muB = muB, muT = muT, muC = muC, 
	muR = muR, muRC = muRC, 	epsilon = epsilon)

################################

# higher baseline (as above), old disease rates
muS <- NA; muT <- NA; muB <- NA; muC <- NA
muS[1:2]<- 1- 0.8 # mortality rate in yearlings 
muS[3:15]<-  1- 0.9 # mortality rate in adults 
muS[16:20]<-  1- 0.8 # mortality rate in adults 15+  (1/yr)

muT <- 2.8 * muS
muB <- 3 * muS
muC <- (2.8 + 3) * muS
muRC <- muC
muR <- muB


test.fixed.params.olddz =  list(aging = aging, 
	b1= b1, b2 = b2, b3 = b3, b4= b4, b5 = b5, b = b, 
	muS = muS, muB = muB, muT = muT, muC = muC, 
	muR = muR, muRC = muRC, 	epsilon = epsilon)

################################



muS <- NA; muT <- NA; muB <- NA; muC <- NA
muS[1:2]<- 1- 0.86 # mortality rate in yearlings 
muS[3:15]<-  1- 0.94 # mortality rate in adults 
muS[16:20]<-  1- 0.86 # mortality rate in adults 15+  (1/yr)

muT <- 2.8 * muS
muB <- 3 * muS
muC <- (2.8 + 3) * muS
muRC <- muC
muR <- muB
fixed.params.olddz =  list(aging = aging, 
	b1= b1, b2 = b2, b3 = b3, b4= b4, b5 = b5, b = b, 
	muS = muS, muB = muB, muT = muT, muC = muC, 
	muR = muR, muRC = muRC, 	epsilon = epsilon)
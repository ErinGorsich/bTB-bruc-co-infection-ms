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
muS[1:3]<- 1- 0.884 # mortality rate in yearlings [1-4) (1/yr)
muS[4]<-  1- 0.963 # mortality rate in sub-adults [4-5)  (1/yr)
muS[4:14]<-  1- 0.963 # mortality rate in adults 5+  (1/yr)
muS[15:20]<-  1- 0.7 # mortality rate in adults 15+  (1/yr)

# mortality, TB, and Brucellosis positive animals
muT <- 2.82 * muS
muB <- 3.02 * muS
muC <- (2.82 + 3.02) * muS
muC[15:20] <- 1  # because otherwise goes over
muRC <- muC
muR <- muB

K = 1000
# only consider female mortality... 
# Mortality in males is 2-8% higher than mortality in females. (?)
#mum<- muS + 0.05

#############################################################
# births
#############################################################
# NOTES: data informing birts are from ages 4-10.  Age 4 is the youngest sucessful mom with calf.
# birth rate in uninfected buffalo of each age category
# now b is the maximum possible birth rate 
b <- NA; b1 <- NA; b2<- NA; b3 <- NA; b4 <- NA; b5 <- NA
# Birth rate in uninfecteds
b[1:3] <- 0						# [1-4)
b[5:14] <- 0.56/2				# age [5-15); raw data (14/25)
b[4] <- max(b[5] * 0.15, 0)		# age [4-5); raw data (2/26)
b[15:20] <- max(b[5] * 0.7/2, 0)	# age 15 +; raw data (14/25)
#b <- b/2						# now max births...
# Proportional reductions
b1 <- c(rep(1, 4), rep(0.65, 20-4))		# with bTB (from raw data)
b2 <- c(1, 1, 1, 1, rep(0.68, 20-4))		# brucellosis
b3<- b2  								# chronic/recovered = active
b4<- c(1, 1, 1, 1, rep(0.8, 20-4)) 		# coinfected
b5 <- b4 								# chronic-coinf = coinfected

epsilon = 0.01

fixed.params = list(aging = aging, 
	b1= b1, b2 = b2, b3 = b3, b4= b4, b5 = b5, b = b, 
	K = K, muS = muS, muB = muB, muT = muT, muC = muC, 
	muR = muR, muRC = muRC, 	epsilon = epsilon)
# missing: gamma, betaT = betaT, betaB = betaB, betapT = betapT, betapB = betapB





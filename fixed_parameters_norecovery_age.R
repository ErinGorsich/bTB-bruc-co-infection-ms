#############################################################
#############################################################
# Fixed parameters, Animals in recoverd class have similar birth/death rates to infecteds 
# Model contains 4 age classes
#############################################################
#############################################################

# translate odds ratios and inital birth probability to new probability
#get_new_prop_birth = function(logOR, p1){
#	OR <- exp(logOR)
#	x <- OR/(p1/(1-p1))
#	b <- x/(p1*(1-x))
#	return(b)
#}

#get_prop_reduction = function(logOR, p1){
#	# input = log odds ratio from model, birth rate in susceptibles
#	OR <- exp(logOR)
#	x <- OR/(p1/(1-p1))  # x = newp/(1-newp)
#	newp <- 1 / ((1/x) + 1)
#	prop <- newp/p1
#	return(prop)
#} # CHECK THIS returned value * p1 should give proportion giveing birth with dz

#get_new_probability = function(logOR, p1){
#	OR <- exp(logOR)
#	p2 <- p1 * OR / (1 - p1 + p1*OR)
#	return(p2)
#}

#test = function(logOR, p1){
#	oldodds <- p1/(1-p1)
#	OR <- exp(logOR)
#	inv <- (1/(OR*oldodds)) + 1
#	pn <- 1/inv
#	proportion <- pn/p1
#	return(proportion)
#}

# use raw data to get proportional reductions

#############################################################
# Mortality, susceptible females
#############################################################
muS <- NA; muT <- NA; muB <- NA; muC <- NA
muS[1]<- 1- 0.884 # mortality rate in yearlings [1-4) (1/yr)
muS[2]<-  1- 0.963 # mortality rate in sub-adults [4-5)  (1/yr)
muS[3]<-  1- 0.963 # mortality rate in adults 5+  (1/yr)
muS[4]<-  1- 0.7 # mortality rate in adults 15+  (1/yr)


# mortality, TB, and Brucellosis positive animals
muT <- 2.82 * muS
muB <- 3.02 * muS
muC <- (2.82 + 3.02) * muS
muC[4] <- 1  # because otherwise goes over
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
b[1] <- 0					# [1-4)
b[3] <- 0.56			    		# Birth rate in uninifected at age [5-15) ;  from raw data (14/25)
b[2] <- max(b[3] * 0.15, 0)	# Birth rate in uninifected at age [4-5); from raw data (2/26)
b[4] <- max(b[3] * 0.7, 0)			# Birth rate in uninifected at age 15 +;  from raw data (14/25)

b1 <- c(1, 1, 1, 0.65)			# prop reduction in fecundity bTB (from raw data)
b2 <- c(1, 1, 1, 0.68)			# prop reduction with brucellosis
b3<- b2  						# chronic/recovered = active
b4<- c(1, 1, 1, 0.8) 			# prop reduction in fecundity if coinfected
b5 <- b4 						# chronic-coinfected = coinfected

#############################################################
# transmission and recovery rates are not fixed #############################################################
# transmission parameters (constant with age)
# betaT =  0.1    		# transmission rate of TB
# rhoT
# betaB
# rhoB
# gamma = 1/2  # recovery rate = 1/2 years

epsilon = 0.01

#############################################################
# Age distribution information
#############################################################
# Age distribution information (? expected values)
#npop = 1000
#f = c(0.08, 0.14, 0.31, 0.29, 0.18) # Age structure based on Northern herds in Caron et al. 2003 (low bTB prevalence)
#N = npop * f
#nu = c(1, 1/2, 1, 1, 0) # 1 / duration of time spent in each age category
#############################################################


fixed.params = list(b1= b1, b2 = b2, b3 = b3, b4= b4, b5 = b5, 
	b = b, K = K, 
	muS = muS, muB = muB, muT = muT, muC = muC, muR = muR, muRC = muRC, 
	epsilon = epsilon)
# missing: gamma, betaT = betaT, betaB = betaB, betapT = betapT, betapB = betapB





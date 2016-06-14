#############################################################
#############################################################
#############################################################
# Age structured SIR model for TB and Brucellosis v.1
# Wednesday 8, June 2016
#############################################################
#############################################################
#############################################################
require("deSolve")

get_new_probability = function(logOR, p1){
	OR <- exp(logOR)
	p2 <- p1 * OR / (1 - p1 + p1*OR)
	return(p2)
}

#############################################################
#############################################################
# parameters
#############################################################
#############################################################
K = 1000 # carrying capacity
# Mortality, susceptible females
muS <- NA; muT <- NA; muB <- NA; muC <- NA
muS[1] <- 1- 0.7 # mortality rate in calves age 0-1 (1/yr)
muS[2]<- 1- 0.884 # mortality rate in juveniles 1-3 (1/yr)
muS[3]<-  1- 0.963 # mortality rate in adults 4-5  (1/yr)
muS[4]<-  1- 0.963 # mortality rate in adults 5+  (1/yr)

# mortality, TB mad Brucellosis positive animals
muT[1] <- 1 - 0.7; muT[2] <- 1- 0.689; muT[3] <- 0.892; muT[4] <- 0.892  
muB[1] <- 1 - 0.7; muB[2] <- 1- 0.706; muB[3] <- 0.899; muB[4] <- 0.899  
muC[1] <- 1 - 0.7; muC[2] <- 1- 0.35; muC[3] <- 0.724; muC[4] <- 0.724  

# Mortality in males is 2-8% higher than mortality in females. 
#mum<- muS + 0.05
#############################################################

# births
# NOTES: data informing birts are from ages 4-10.  Age 4 is the youngest sucessful mom with calf.
b0 <- NA
b0[1] <- 0
b0[2] <- 0
b0[3] <- 0.2 #get_new_probability(-2.684, 0.56)
b0[4] <- 0.56  # Birth rate in uninifected buffalo 
# without age structure r = b-mu; set b to get r=10% growth rate w/out density dependence, e.g. Cross et al. 2009.

#b1 <- get_new_probability(-1.619, b0) # proportional reduction in births with TB. Convert OR to proportions ratio: exp(-1.619)/(1+ exp(-1.619))
#b2 <- get_new_probability(-1.50, b0) # proportional reduction in births with brucellosis
#b3 <- get_new_probability((-1.619 - 1.50 + 2.28), b0) # proportional reduction in births if co-infected
# initial run
b1 <- c(0, 0, 1, 30/56)
b2 <- c(0, 0, 1, 30/56)
b3 <- c(0, 0, 1, 43/56)

# transmission parameters (constant with age)
betaT =  0.02    		# transmission rate of TB
betapT <- betaT  		# transnmission rate of TB in brucellosis infected buffalo
betaB = 0.02      		# transmission rate of brucellosis
betapB = 1.4 * betaB	# tranmsission rate of bruceillosis in TB + buffalo

gamma = 1/2  # recovery rate = 1/2 years
epsilon = 0.01
rho <- 0.05  # rate of vertical transmission

# Age distribution information 
npop = 1000
f = c(0.1, 0.25, 0.05, 0.7) # later make this the stable age distribution...
N = npop * f
nu = c(1, 1/2, 1, 0) # 1/duration of time spent in each age category
#############################################################


#############################################################
#############################################################
# Write ODEs
#############################################################
#############################################################
# 6 variables, each with four age classes

# NOW ALL FEMALES, WILL LATER MAYBE CONSIDER BOTH MALES AND FEMALES...
SIRfunc_age = function(t, x, params){
	# Inputs: t = time sequence; x = initial conditions in one vector in order of variable
	
	# Assign state variables, each 4 long, to capture the 4 categories (age= 0-1, 1-3, 4-5, 5+)
	nage = length(x)/4
	S = x[[1]]
	It = x[[2]]
	Ib = x[[3]]
	Ic = x[[4]]
	Rb = x[[5]]
	Rc = x[[6]]
	
	# make sure things make sense, (NEED TO ADD WARNING ERROR MESSAGES)
	Ib[Ib<0] = 0
	It[It<0] = 0
	Ic[Ic<0] = 0
	
	# Overall population size, and population contributing to births
	N <- S + It + Ib + Ic + Rb + Rc
	Nb <- S + b1 * (It+Rc) + b2 * Ib + b3 * Ic + Rb  # assuming recovereds don't suffer birth consequence of brucellosis. 
	
	# Frequency dependent force of infection is independent of age
	lambdaT <- betaT * (Ib + Ic + Rc) / N
	lambdaB <- betaB * (Ib + Ic) / N 
	lambdapT <- betapT * (Ib + Ic + Rc) / N
	lambdapB <- betapB * (Ib + Ic) / N 
	
	# differential equations, assume mortality in chronics is same as active infection
	for (i in 1:nage){
		dS[i] <- b0[i] * (1- rho) * Nb * (1 - (r/b) * (N/K)) - (lambdaT + lambdaB) * S[i] - muS[i] * S[i] - nu[i] * S[i]
		dIt[i] <- lambdaT * S[i] - (lambdapB + muT[i] + nu[i]) * It[i] 
		dIb[i] <- b0[i] * rho * Nb * (1 - (r/b) * (N/K)) + lambdaB * S[i] + epsilon * Rb[i] - gamma * Ib[i] - (lambdapT + muB[i] -nu[i]) * Ib[i]
		dIc[i] <- lambdapT * Ib[i] + lambdapB * It[i] - gamma * Ic[i] + epsilon * Rc[i] - (muC[i] + nu[i]) * Ic[i]
		dRb[i] <- gamma * Ib[i] - epsilon * Rb[i] - (muB[i] + nu[i]) * Rb[i] 
		dRc[i] <- lambdapT * Rb[i] + gamma * Ic[i] - epsilon * Rc[i] - (muC[i] + nu[i]) * Rc[i]
	}
}

#############################################################
#############################################################
# Test Run
#############################################################
#############################################################
#Initial Conditions
S0 = N-20
It0 = c(10, 10, 10, 10)
Ib0 = c(10, 10, 10, 10)
Ic0 = c(0, 0, 0, 0)
Rb0 = c(0, 0, 0, 0)
Rc0 = c(0, 0, 0, 0)
x = list(S0, It0, Ib0, Ic0, Rb0, Rc0)




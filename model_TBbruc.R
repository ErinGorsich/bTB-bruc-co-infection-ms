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

#############################################################
# Mortality, susceptible females
#############################################################
muS <- NA; muT <- NA; muB <- NA; muC <- NA
muS[1] <- 1- 0.7 # mortality rate in calves age [0-1) (1/yr)
muS[2]<- 1- 0.884 # mortality rate in yearlings [1-3) (1/yr)
muS[3] <- 1- 0.884 # mortality rate in juveniles [3-4) (1/yr)
muS[4]<-  1- 0.963 # mortality rate in sub-adults [4-5)  (1/yr)
muS[5]<-  1- 0.963 # mortality rate in adults 5+  (1/yr)

# mortality, TB and Brucellosis positive animals
muT[1] <- 1 - 0.7; muT[2] <- 1- 0.689; muT[3] <- 1-0.689; muT[4] <- 1- 0.892; muT[5] <- 1- 0.892; 
muB[1] <- 1 - 0.7; muB[2] <- 1- 0.706; muB[3] <- 1- 0.706; muB[4] <- 1- 0.899; muB[5] <- 1- 0.899;  
muC[1] <- 1 - 0.7; muC[2] <- 1- 0.35; muC[3] <- 1- 0.35; muC[4] <- 1- 0.724; muC[5] <- 1- 0.724
muRB <- muB
muRC <- muC

# Mortality in males is 2-8% higher than mortality in females. (?)
#mum<- muS + 0.05
#############################################################
# births
#############################################################
# NOTES: data informing birts are from ages 4-10.  Age 4 is the youngest sucessful mom with calf.

# birth rate in uninfected buffalo of each age category
b0 <- NA
b0[1] <- 0
b0[2] <- 0
b0[3] <- 0   
b0[4] <- 0.08   # Birth rate in uninifected buffalo at age [4-5); from raw data (2/26)
b0[5] <- 0.56   # Birth rate in uninifected buffalo at age 5 +;  from raw data (14/25)

b1 <- get_new_probability(-1.619, b0[5])   	# proportional reduction in fecundity with bTB, statistical model output
b2 <- get_new_probability(-1.5, b0[5])     	# proportional reduction in fecundity with brucellosis
b3<- b2  									# proportional reduction in fecundity with chronic brucellosis (Rb), assume ????????
b4<- get_new_probability((-1.619 - 1.50 + 2.28), b0[5]) # proportional reduction in fecundity if co-infected (Ic)
b5 <- b4 									# proportional reduction in fecundity if have chronic brucellosis and bTB (Rc)

#############################################################
# transmission and recovery rates to vary later. 
#############################################################
# transmission parameters (constant with age)
betaT =  0.1    		# transmission rate of TB
betapT <- betaT  		# transnmission rate of TB in brucellosis infected buffalo
betaB = 0.02      		# transmission rate of brucellosis
betapB = 3.92 * betaB	# tranmsission rate of bruceillosis in TB + buffalo

gamma = 1/2  # recovery rate = 1/2 years
epsilon = 0.01
rho <- 0.05  # rate of vertical transmission

#############################################################
# Age distribution information
#############################################################
# Age distribution information (? expected values)
npop = 1000
f = c(0.08, 0.14, 0.31, 0.29, 0.18) # Age structure based on Northern herds in Caron et al. 2003 (low bTB prevalence)
N = npop * f
nu = c(1, 1/2, 1, 1, 0) # 1 / duration of time spent in each age category
#############################################################


#############################################################
#############################################################
# Write ODEs, density dependence... 
#############################################################
#############################################################
rhs = function(times, x, params){
	##########################
	# Inputs: t = time sequence; x = initial conditions in one vector in order of variable orderd juv, adult for each category
	# params= list(
	# b1, b2, b3, b4, b5,  proportional reduction in fecundity with TB, bruc, chronic bruc, co-infect, coinfect/chronic (1 value each);
	# betaT, betapT, betaB, betapB, transmission rate of TB, without/ with co-infection, trans Bruc without/with co-infection (1 value each);
	# rho, proportion of buffalo born from infectious mothers that are infectious
	# gamma, epsilon
	# b = max birth rate (when N = 0) 
	# r = growth rate (set = b-d)???? with ????????
	# Age- specific rates (vectors of length 5)
	# muSa, muBa, muTa, muCa; muRBa; muRCa, muSj, muBj, muTj, muCj; muRBj; muRCj mortality rates.  Assume chronic stages have mortality rates equal to acute stage
	# nu, inverse of the duration of time spent in each age class
	##########################
	with(as.list(c(x, params)), {
		
		nage = 5
		# Assign state variables, each 4 long, to capture the 4 categories (age= 0-1, 1-3, 4-5, 5+)
		Sj = x[1] 
		Sa = x[2]
		Itj = x[3]
		Ita = x[4]
		Ibj = x[5]
		Iba = x[6]
		Icj = x[7]
		Ica = x[8]
		Rbj = x[9]
		Rba = x[10]
		Rcj = x[11]
		Rca = x[12]
			
		# Overall population size (N), and population contributing to births (Nb) summd over age classes
		N <- sum(Sj + Sa + Itj + Ita + Ibj + Iba + Icj + Ica + Rbj + Rba + Rcj+ Rca)
		Nbirthu <- Sa + b1 * Ita + b3 * Rba + b5 * Rca  # 5 long vector with age specific numbers contributing to births
		Nbirthb <- b2 * Iba + b4 * Ica
	
		# Frequency dependent force of infection is independent of age
		lambdaT <- betaT * sum(Ita + Ica + Rca + Itj + Icj + Rcj) / N
		lambdaB <- betaB * sum(Iba + Ibj + Ica + Icj) / N 
		lambdapT <- betapT * sum(Ita + Ica + Rca + Itj + Icj + Rcj) / N
		lambdapB <- betapB * sum(Iba + Ibj + Ica + Icj) / N 
	
		#### r = b - d
	
		# differential equations, assume mortality in chronics is same as active infection
		
		# juveniles
		dSj <- (b - r * (N/K)) * Nbirthu - (lambdaT + lambdaB) * Sj - muSj * Sj - nu * Sj
		dItj <- lambdaT * Sj - (lambdapB + muTj + nu) * Itj 
		dIbj <- (b - r * (N/K)) * Nbirthb	 + lambdaB * Sj + epsilon * Rbj - (gamma + lambdapT + muBj -nu) * Ibj
		dIcj <- lambdapT * Ibj + lambdapB * Itj + epsilon * Rcj - (gamma + muCj + nu) * Icj
		dRbj <- gamma * Ibj - (epsilon + muRBj + nu) * Rbj
		dRcj <- lambdapT * Rbj + gamma * Icj - (epsilon + muRCj + nu) * Rcj
		
		# adults
		dSa <-  nu * Sj - (lambdaT + lambdaB) * Sa - muSa * Sa 
		dIta <- nu * Itj + lambdaT * Sa - (lambdapB + muTa) * Ita 
		dIba <- nu * Ibj + lambdaB * Sa + epsilon * Rba - (gamma + lambdapT + muBa) * Iba
		dIca <- nu * Icj + lambdapT * Iba + lambdapB * Ita + epsilon * Rca - (gamma + muCa) * Ica
		dRba <- nu * Rbj + gamma * Iba - (epsilon + muRBa) * Rba
		dRca <- nu * Rcj + lambdapT * Rba + gamma * Ica - (epsilon + muRCa) * Rca

		
		out = list(c(dSj, dSa, dItj, dIta, dIbj, dIba, dIcj, dIca, dRbj, dRba, dRcj, dRca))
		return(out)
		}
	)
}
	

#############################################################



	
#############################################################
#############################################################
# Test Run
#############################################################
#############################################################
params = c(b1= b1, b2 = b2, b3 = b3, b4= b4, b5 = b5, 
	muS = muS, muB = muB, muT = muT, muRC = muRC, muRB = muRB, 
	nu = nu, gamma = gamma, epsilon = epsilon, rho = rho,
	betaT = betaT, betaB = betaB, betapT = betapT, betapB = betapB)

#Start with brucellosis only, make sure no bTB takes off
#############################################################
S0 = N-10
It0 = c(0, 0, 0, 0, 0)
Ib0 = c(10, 10, 10, 10, 10)
Ic0 = c(0, 0, 0, 0, 0)
Rb0 = c(0, 0, 0, 0, 0)
Rc0 = c(0, 0, 0, 0, 0)
x0 = c(S0, It0, Ib0, Ic0, Rb0, Rc0)
times <- seq(0, 100, 1)
sol <- as.data.frame(ode(x0, times, rhs, params))
groom_sol <- function(sol){
	colnames(sol) <- c("times", "S_calf", "S_yearling", "S_juv", "S_subadult", "S_adult",
		"It_calf", "It_yearling", "It_juv", "It_subadult", "It_adult",
		"Ib_calf", "Ib_yearling", "Ib_juv", "Ib_subadult", "Ib_adult",
		"Ic_calf", "Ic_yearling", "Ic_juv", "Ic_subadult", "Ic_adult",
		"Rb_calf", "Rb_yearling", "Rb_juv", "Rb_subadult", "Rb_adult",
		"Rc_calf", "Rc_yearling", "Rc_juv", "Rc_subadult", "Rc_adult") 
	sol$S <- sol$S_calf + sol$S_yearling + sol$S_juv + sol$S_subadult + sol$S_adult
	sol$It <- sol$It_calf + sol$It_yearling + sol$It_juv + sol$It_subadult + sol$It_adult
	sol$Ib <- sol$Ib_calf + sol$Ib_yearling + sol$Ib_juv + sol$Ib_subadult + sol$Ib_adult
	sol$Ic <- sol$Ic_calf + sol$Ic_yearling + sol$Ic_juv + sol$Ic_subadult + sol$Ic_adult
	sol$Rb <- sol$Rb_calf + sol$Rb_yearling + sol$Rb_juv + sol$Rb_subadult + sol$Rb_adult
	sol$Rc <- sol$Rc_calf + sol$Rc_yearling + sol$Rc_juv + sol$Rc_subadult + sol$Rc_adult
	return(sol)
}
sol<- groom_sol(sol)

plot_raw_numbers = function(sol){
	plot(sol$time, sol$S, col= "black", type= 'l', ylim = c(0, 5000))
	lines(sol$time, sol$It, col= "red")
	lines(sol$time, sol$Ib, col= "dark blue")
	lines(sol$time, sol$Ic, col= "dark green")
	lines(sol$time, sol$Rb, col = "blue")
	lines(sol$time, sol$Rc, col = "green")
	legend("topright", legend = c("S", "It", "Ib", "Ic", "Rb", "Rc"), col = c("black", "red", "dark blue", "dark green", "blue", "green"), bty="n", lty = 1)
}
plot_raw_numbers(sol)

#Start with bTB only, make sure no brucellosis takes off
#############################################################
Ib0 = c(0, 0, 0, 0, 0)
It0 = c(10, 10, 10, 10, 10)
x0 = c(S0, It0, Ib0, Ic0, Rb0, Rc0)
sol <- as.data.frame(ode(x0, times, rhs, params))
sol<- groom_sol(sol)
plot_raw_numbers(sol)


# Get realistic parameter values, feel for what is going on... 
#############################################################













#############################################################
#############################################################
# Write ODEs, non-density dependence... 
#############################################################
#############################################################
# 6 variables, each with five age classes

# NOW ALL FEMALES, WILL LATER MAYBE CONSIDER BOTH MALES AND FEMALES...
#rhs = function(times, x, params){
	##########################
	# Inputs: t = time sequence; x = initial conditions in one vector in order of variable
	# params= list(
	# b1, b2, b3, b4, b5,  proportional reduction in fecundity with TB, bruc, chronic bruc, co-infect, coinfect/chronic (1 value each);
	# betaT, betapT, betaB, betapB, transmission rate of TB, without/ with co-infection, trans Bruc without/with co-infection (1 value each);
	# rho, proportion of buffalo born from infectious mothers that are infectious
	# gamma, epsilon
	
	# Age- specific rates (vectors of length 5)
	# muS, muB, muT, muC; muRB; muRC mortality rates.  Assume chronic stages have mortality rates equal to acute stage
	# nu, inverse of the duration of time spent in each age class
	##########################
#	with(as.list(c(x, params)), {
#		
#		nage = 5
#		# Assign state variables, each 4 long, to capture the 4 categories (age= 0-1, 1-3, 4-5, 5+)
#		S = x[1:nage]
#		It = x[(nage+1): (2*nage)]
#		Ib = x[(2*nage+1): (3*nage)]
#		Ic = x[(3*nage+1): (4*nage)]
#		Rb = x[(4*nage+1): (5*nage)]
#		Rc = x[(5*nage+1): (6*nage)]
#	
#		# make sure things make sense (these should be unnecessary, because if you are at 0, nothing should be comming out... )
#		#if (Ib[Ib<0]){Ib[Ib<0]<- 0}
#		#if (It[It<0]){It[It<0]<- 0}
#		#if (Ic[Ic<0]){Ic[Ic<0]<- 0}
#		
#		# Overall population size (N), and population contributing to births (Nb) summd over age classes
#		N <- sum(S + It + Ib + Ic + Rb + Rc)
#		Nbirth <- S + b1 * It + b3 * Rb + b5 * Rc  # 5 long vector with age specific numbers contributing to births
#		Nbruc <- b2 * Ib + b4 * Ic
#		bu <- sum(b0 * Nbirth + b0 * Nbruc * (1-rho))
#		bb <- sum(b0 * rho * Nbruc)
#	
#		# Frequency dependent force of infection is independent of age
#		lambdaT <- betaT * sum(It + Ic + Rc) / N
#		lambdaB <- betaB * sum(Ib + Ic) / N 
#		lambdapT <- betapT * sum(It + Ic + Rc) / N
#		lambdapB <- betapB * sum(Ib + Ic) / N 
#		
#		dS <- vector(length=nage); dIt <- vector(length=nage); dIb<- vector(length=nage);
#		dIc <- vector(length=nage); dRb <- vector(length=nage); dRc <- vector(length=nage)
#	
#		# differential equations, assume mortality in chronics is same as active infection
#		for (i in 1:length(S0)){
#			dS[i] <- bu - (lambdaT + lambdaB) * S[i] - muS[i] * S[i] - nu[i] * S[i]
#			dIt[i] <- lambdaT * S[i] - (lambdapB + muT[i] + nu[i]) * It[i] 
#			dIb[i] <- bb + lambdaB * S[i] + epsilon * Rb[i] - (gamma + lambdapT + muB[i] -nu[i]) * Ib[i]
#			dIc[i] <- lambdapT * Ib[i] + lambdapB * It[i] + epsilon * Rc[i] - (gamma + muC[i] + nu[i]) * Ic[i]
#			dRb[i] <- gamma * Ib[i] - (epsilon + muRB[i] + nu[i]) * Rb[i] 
#			dRc[i] <- lambdapT * Rb[i] + gamma * Ic[i] - (epsilon + muRC[i] + nu[i]) * Rc[i]
#		}
#		out = list(c(	dS, dIt, dIb, dIc, dRb, dRc))
#		return(out)
#		}
#	)
#}
	



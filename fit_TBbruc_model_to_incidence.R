#############################################################
#############################################################
# Data, assume poisson distributed likelihood
#############################################################
#############################################################
new_TB_LS <- c(0, 1, 10, 3, 3, 2, 2, 1)
new_TB_CB <- c(0, 0, 4, 1, 5, 2, 2, 1)
new_B_LS <- c(0, 0, 0, 1, 4, 4, 3, 2)
new_B_CB <- c(0, 4, 0, 4, 0, 3, 1, 1)
#############################################################
#############################################################
# Write ODEs, density dependence... no age structure; \
# note betapT and betapB not here
# Keeping track of more than in model_TBbruc_derivs
#############################################################
#############################################################
require("deSolve")
source('~/GitHub/bTB-bruc-co-infection-ms/fixed_parameters.R', chdir = TRUE)
agedata <- read.csv("~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/final_datasets_copied_from_phdfolder/cross_sectional_data_withdz_cleandisease_nofinal_Feb2016_capturetime.csv")
LSage <- agedata$age_yr[agedata$herdorig == "LS" & agedata$capturetime == 0]
CBage <- agedata$age_yr[agedata$herdorig == "CB" & agedata$capturetime == 3]

LS<- agedata[agedata$herdorig=="LS",]
LSmass <- LS[LS$capturetime==0,]
LSmass[,c(3, 22, 30)]



#############################################################
rhs = function(times, x, params){
	##########################
	# Inputs: t = time sequence; x = initial conditions in one vector in order of variable orderd juv, adult for each category
	# params= list(
	# b1, b2, b3, b4, b5,  prop reduction in fecundity with TB, bruc, chronic bruc, co, coinfect/chronic 
	# K 
	# betaT, betapT, betaB, betapB, transmission rates, p indicates co-infected first
	# gamma, epsilon, ignore rho for now (proportion of buffalo born from infectious mothers that are infectious).
	# b = max birth rate (when N = 0) 
	# muS, muT, muB, muC, muR, muRC mortality rates.  
	# NOW Assume chronic stages have mortality rates equal to acute stage!!!
	##########################
	with(as.list(c(x, params)), {
					
		# Overall population size (N)
		N <- S + It + Ib + Ic + R + Rc
		
		# Population contributing to births (Nb); account for reduced births with infection
		Nb <- S + b1 * It + b2 * Ib + b3 * Ic + b4 * R + b5 * Rc
			
		# Frequency dependent force of infection is independent of age
		betapB <- 3.92 * betaB
		betapT <- betaT
		lambdaT <- betaT * (It + Ic + Rc) / N
		lambdaB <- betaB * (Ib + Ic) / N 
		lambdapT <- betapT * (It + Ic + Rc) / N
		lambdapB <- betapB * (Ib + Ic) / N
		
		# differential equations, assume mortality in chronics is same as active infection
		dS <- b * Nb * (1 - (r/b) * (N/K) ) - (lambdaT + lambdaB) * S - muS * S 
		dIt <- lambdaT * S - (lambdapB + muT) * It 
		dIb <- lambdaB * S + epsilon * R - (gamma + lambdapT + muB) * Ib
		dIc <- lambdapT * Ib + lambdapB * It + epsilon * Rc - (gamma + muC) * Ic
		dR <-  gamma * Ib - (epsilon + muR + lambdapT) * R
		dRc <- lambdapT * R + gamma * Ic - (epsilon + muRC) * Rc
		
		# derived quantities (number of co-infecteds with Bruc first/ TB first)
		ifelse (Ic1 >= 0, 
			dIc1 <- lambdapT * Ib + epsilon * Rc1 - (gamma + muC) * Ic1, 
			dIc1 <- 0) # ???????
		ifelse (Ic2 >= 0, 
			dIc2 <- lambdapB * It + epsilon * Rc2 - (gamma + muC) * Ic2, 
			dIc2 <- 0)
		ifelse (Rc1 >= 0,
			dRc1 <- lambdapT * R  + gamma * Ic1 - (epsilon+ muRC) * Rc1, 
			dRc1 <- 0)
		ifelse(Rc2 >= 0, 
			dRc2 <- gamma * Ic2 - (epsilon+ muRC) * Rc2, 
			dRc2 <- 0)

		# cumulative number of bTB cases from Susceptibles that did not die
		dTBfS <- lambdaT * S - muT * It	
		# cumulative number of brucellosis cases from Susceptibles that did not die
		dBfS <- lambdaB * S - muB * Ib
		# cumulative number of bTB cases from buffalo with brucellosis first
		dTBfB <- dIc1 + dRc1
		# cumulative number of brucellosis cases from buffalo with bTB first
		dBfTB <- dIc2 + dRc2
		
		
		out = list(c(dS, dIt, dIb, dIc, dR, dRc, dIc1, dIc2, dRc1, dRc2, dTBfS, dBfS, dTBfB, dBfTB))
		return(out)
		}
	)
}

plot_raw_numbers = function(sol){
	plot(sol$time, sol$S, col= "black", type= 'l', ylim = c(0, 1200), ylab = "Number of animals", xlab = "Time (in years)")
	lines(sol$time, sol$It, col= "red")
	lines(sol$time, sol$Ib, col= "blue")
	lines(sol$time, sol$Ic, col= "green")
	lines(sol$time, sol$R, col = "orange")
	lines(sol$time, sol$Rc, col = "pink")
	legend("topright", legend = c("S", "It", "Ib", "Ic", "R", "Rc"), col = c("black", "red", "blue", "green", "orange", "pink"), bty="n", lty = 1)
}


#############################################################
# Test Runs
par(mfrow = c(1, 2))
S0 = 500
It0 = 50
Ib0 = 50
Ic0 = 0
R0 = 0
Rc0 = 0
times <- seq(0, 100, 1)
betaT <- 0.4
betaB <- 0.6
params = c(b1= b1, b2 = b2, b3 = b3, b4= b4, b5 = b5, 
	b = b, r = r, K = K,
	muS = muS, muB = muB, muT = muT, muR = muR, muRC = muRC, 
	gamma = gamma, epsilon = epsilon,
	betaT = betaT, betaB = betaB)  # betapT & betapB specified in function
x0 = c(S = S0, It = It0, Ib = Ib0, Ic = Ic0, R = R0, Rc = Rc0,
	Ic1 = 0, Ic2 = 0, Rc1 = 0, Rc2 = 0, 
	TBfS = 0, BfS = 0, TBfB = 0, BfTB = 0)
sol <- as.data.frame(ode(x0, times, rhs, params))  # returns as many columns as elements in x0
plot_raw_numbers(sol)

#############################################################
# with longer lower mortality in buffalo recovered from brucellosis
params = c(b1= b1, b2 = b2, b3 = b3, b4= b4, b5 = b5, 
	b = b, r = r, K = K,
	muS = muS, muB = muB, muT = muT, muR = muS, muRC = muT, 
	gamma = gamma, epsilon = epsilon,
	betaT = betaT, betaB = betaB) 
sol <- as.data.frame(ode(x0, times, rhs, params))  # returns as many columns as elements in x0
plot_raw_numbers(sol)

#############################################################
#############################################################
# Add the simulation/sampling parts
#############################################################
#############################################################

get_num_new_infections = function(betavals, params.fixed, herd){
	# Run model at parameters
    params0 = c(betaT = betavals[1], betaB = betavals[2])
    Y0 <- c(S = 450, It = 20, Ib = 20, Ic = 0, R = 0, Rc = 10)
	sol <- as.data.frame(ode(Y0, times, rhs, params= c(params0, params.fixed)))

	# Randomly sample buffalo of age specified in first capture of LS/CB
	#############################################################
	# First get ages
	if (herd == "LS"){
		Page = LSage
	} else if (herd == "CB"){
		Page = CBage
	}
	# Define P(infection|age) based on model parameters... (now end state)
	
	
}


#############################################################
#############################################################
# Fit to Data by maximizing poisson likelihood, deterministic
#############################################################
#############################################################
# ASSUMTION 1: FIXED RECOVERY, REACTIVATION RATE; muR = muS; muRC = muT
#############################################################
ML.sir = function(betavals, params.fixed, data){
	##############################################
	# INPUT: 
	# betavals = c(betat, betab)
	# params.fixed = all other parameters
	# data = number of new cases per capture (0.5/yr)= a list of length ?
	# Returns: NLL value to be optimized. 
    	##############################################

	# Data
    obs_num_new_tb <- data$tb  # calculated as change in n frm t -1 to t
    obs_num_new_bruc <- data$br
    # obs_num_new_tb_co <- data$tbco
    # obs_num_new_bruc_co <- data$brucco

	# Run model at parameters
    params0 = c(betaT = betavals[1], betaB = betavals[2])
    Y0 <- c(S = 450, It = 20, Ib = 20, Ic = 0, R = 0, Rc = 10)
	sol <- as.data.frame(ode(Y0, times, rhs, params= c(params0, params.fixed)))

	# Randomly sample 200 buffalo, with age distribution in data...
	
	# Expected number of new cases
    pred_num_new_tb <- data$It
    pred_num_new_tb <- data$It + data$

    times <- data$Time
    
    pred_num_new_tb <- sol$
    pred_num_new_bruc <- sol$
    
    NLL_TB <- - sum(dpois(x=?? , lambda = ???, log=TRUE)) 
    NLL_BR <- - sum(dpois(x=?? , lambda = ???, log=TRUE)) 
	NLL <- NLL_TB + NLL_BR
	return(NLL)
}


params.fixed = 
times <- seq(0, 500, 0.5)  # now want data in half year time intervals.



# ASSUMTION 2: FIXED RECOVERY, REACTIVATION RATE; muR = muR; muRC = muRC
#############################################################



# ASSUMTION 3: ESTIMATE RECOVERY, REACTIVATION RATE; muR = muS; muRC = muT
#############################################################
#############################################################
#############################################################
#############################################################
# Age structured SIR model for TB and Brucellosis v.1
# Last update: 1 September 2016
#############################################################
#############################################################
#############################################################
require("deSolve")
require(plyr)

get_new_prop_birth = function(logOR, p1){
	OR <- exp(logOR)
	x <- OR/(p1/(1-p1))
	b <- x/(p1*(1-x))
	return(b)
}

#############################################################
#############################################################
# parameters- set values in analyses * age structure
#############################################################
#############################################################
relageall = c(0.137, rep(0.368/4, 4), rep(0.185/4, 4),  # Jolles 2005, set max age at 18
	rep(0.235/6, 6), rep(0.075/3, 3))					# Also in Caron et al. from 2001 KNP
relage = c(relageall[1], sum(relageall[2:3]), relageall[4],  
	relageall[5], sum(relageall[6:length(relageall)]) )  # sums to 1
N = 1000
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
# Births
#############################################################
b <- 0.41 # proportion in LS right before calving that were pregs (or had milk-> assume, can check this)
# In interdrought periods, buffalo populations grew at rates from 5-15%. 
# This is in populations with brucellosis but not bTB.
# set r = 1.1... in susceptible poulations.
b1 <- get_new_prop_birth(-1.619, b) # proportional reduction with bTB
b2 <- get_new_prop_birth(-1.5, b)   # proportional reduction with bruc
b3 <- b2 # ASSUME animals with chronic brucellosis have same fecundity as active brucellosis 
b4 <- get_new_prop_birth((-1.619 - 1.50 + 2.28), b) # proportional reduction if co-infected
b5 <- b4 # 
r <- b - muS   # Gao & Hethcote; r = birth - death

#############################################################
# transmission and recovery rates to vary later. 
#############################################################
# transmission parameters (constant with age)
betaT =  0.001           # 0.043 transmission rate of TB (with N = 250)
betapT <- betaT  		# transnmission rate of TB in brucellosis infected buffalo
betaB = 0.001      		# transmission rate of brucellosis
betapB = 3.92 * betaB	# tranmsission rate of bruceillosis in TB + buffalo

gamma = 1/3  # recovery rate = 1/2 years
epsilon = 0.01
#rho <- 0.05  # rate of vertical transmission


# Jolles et al. 2009 modeled bTB as density dependent
# Cross & Getz 2008 model bTB as frequency dependent (but with a rel constat population size)
#############################################################
#############################################################
# Write ODEs, density dependence... no age structure
#############################################################
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
		Nb <- S + b1 * It + b2 * Ib + b3 *R + b4 * Ic + b5 * Rc
					
		# Frequency dependent force of infection is independent of age
		lambdaT <- betaT * (It + Ic + Rc) 
		lambdaB <- betaB * (Ib + Ic) 
		lambdapT <- betapT * (It + Ic + Rc)
		lambdapB <- betapB * (Ib + Ic) 
		
		# differential equations, assume mortality in chronics is same as active infection
		dS <- b * Nb * (1 - (r/b) * (N/K) ) - (lambdaT + lambdaB) * S - muS * S 
		dIt <- lambdaT * S - (lambdapB + muT) * It 
		dIb <- lambdaB * S + epsilon * R - (gamma + lambdapT + muB) * Ib
		dIc <- lambdapT * Ib + lambdapB * It + epsilon * Rc - (gamma + muC) * Ic
		dR <-  gamma * Ib - (epsilon + muR + lambdapT) * R
		dRc <- lambdapT * R + gamma * Ic - (epsilon + muRC) * Rc
		
		out = list(c(dS, dIt, dIb, dIc, dR, dRc))
		return(out)
		}
	)
}


groom_sol = function(sol){
	colnames(sol) <- c("times", "S", "It", "Ib", "Ic", "R", "Rc") 
	sol$N <- sol$S + sol$It + sol$Ib + sol$Ic + sol$R + sol$Rc
	sol$TBprev <- (sol$It + sol$Ic + sol$Rc) / sol$N
	sol$Brucprev <- (sol$Ib + sol$R + sol$Ic + sol$Rc) / sol$N
	sol$propTB_co <- (sol$Ic + sol$Rc) / (sol$It + sol$Ic + sol$Rc)
	sol$propBruc_co <- (sol$Ic + sol$Rc) / (sol$Ib + sol$R + sol$Ic + sol$Rc)
	return(sol)
}

plot_raw_numbers = function(sol){
	plot(sol$time, sol$S, col= "black", type= 'l', ylim = c(0, 1200), ylab = "Number of animals", xlab = "Time (in years")
	lines(sol$time, sol$It, col= "red")
	lines(sol$time, sol$Ib, col= "blue")
	lines(sol$time, sol$Ic, col= "green")
	lines(sol$time, sol$R, col = "orange")
	lines(sol$time, sol$Rc, col = "pink")
	legend("topright", legend = c("S", "It", "Ib", "Ic", "R", "Rc"), col = c("black", "red", "blue", "green", "orange", "pink"), bty="n", lty = 1)
}
#############################################################

	
#############################################################
#############################################################
# Test Run- 
#############################################################
#############################################################
params = c(b1= b1, b2 = b2, b3 = b3, b4= b4, b5 = b5, 
	b = b, r = r, K = K,
	muS = muS, muB = muB, muT = muT, muC = muC, muR = muR, muRC = muRC, 
	gamma = gamma, epsilon = epsilon,
	betaT = betaT, betaB = betaB, betapT = betapT, betapB = betapB)
# params assuming fecundity/survival effects of brucellosis are limited to infected class
params_recov = c(b1= b1, b2 = b2, b3 = b3, b4= 1, b5 = b2, # assuming b4 =1, b5 = b2 
	b = b, r = r, K = K,
	muS = muS, muB = muB, muT = muT, muC = muC, muR = muS, muRC = muT, 
	gamma = gamma, epsilon = epsilon,
	betaT = betaT, betaB = betaB, betapT = betapT, betapB = betapB)

#Start no disease, make sure ends with none, pop goes to K
#############################################################
S0 = 500
It0 = 0
Ib0 = 0
Ic0 = 0
R0 = 0
Rc0 = 0
x0 = c(S = S0, It = It0, Ib = Ib0, Ic = Ic0, R = R0, Rc = Rc0)
times <- seq(0, 100, 1)
sol <- as.data.frame(ode(x0, times, rhs, params))
sol_recov <- as.data.frame(ode(x0, times, rhs, params_recov))
par(mfrow=c(1,2))
plot_raw_numbers(sol)
plot_raw_numbers(sol_recov)

#Start with brucellosis only, make sure no bTB takes off
#############################################################
S0 = N-100
It0 = 0
Ib0 = 100
Ic0 = 0
R0 = 0
Rc0 = 0
x0 = c(S = S0, It = It0, Ib = Ib0, Ic = Ic0, R = R0, Rc = Rc0)
times <- seq(0, 100, 1)
sol <- as.data.frame(ode(x0, times, rhs, params))
sol_recov <- as.data.frame(ode(x0, times, rhs, params_recov))
sol<- groom_sol(sol)
par(mfrow = c(1,2))
plot_raw_numbers(sol)
plot_raw_numbers(sol_recov) # with lax in mortality effects in recovered, much more recoverds... 


#Start with bTB only, make sure no brucellosis takes off
#############################################################
Ib0 = N - 100
It0 = 100
Ib0 = 0
Ic0 = 0
R0 = 0
Rc0 = 0
x0 = c(S = S0, It = It0, Ib = Ib0, Ic = Ic0, R = R0, Rc = Rc0)
sol <- as.data.frame(ode(x0, times, rhs, params))
sol<- groom_sol(sol)
plot_raw_numbers(sol)



#############################################################
#############################################################
# Explore realistic parameter values, feel for what is going on... 
#############################################################
#############################################################

# TB only 
#############################################################
K = 1000; N= 1000
beta_t <- seq(0.0001, 0.001, 0.0003)
Ib0 = N - 100
It0 = 100
Ib0 = 0
Ic0 = 0
R0 = 0
Rc0 = 0
x0 = c(S = S0, It = It0, Ib = Ib0, Ic = Ic0, R = R0, Rc = Rc0)
times <- seq(0, 1000, 1)
out <- list(NA); 

par(mfrow = c(1, 2))
plot(NULL, xlim = c(0, 500), xlab= "Time (in years)", ylim= c(0, 1300), 
	ylab = "Number of bTB positive animals")
legend("topright", 	title= "beta_t", legend=c("0.0001", "0.0004", "0.0007","0.001"),
	col = seq(1:length(beta_t)), lty = 1, bty = "n") 	
for(i in 1:length(beta_t)){
	params <- c(b1= b1, b2 = b2, b3 = b3, b4= b4, b5 = b5, 
		b = b, r = r, K = K,
		muS = muS, muB = muB, muT = muT, muC = muC, muR = muR, muRC = muRC, 
		gamma = gamma, epsilon = epsilon,
		betaT = beta_t[i], betaB = betaB, betapT = betapT, betapB = betapB)
	sol <- as.data.frame(ode(x0, times, rhs, params))
	sol<- groom_sol(sol)
	out[[i]] <- sol
	lines(y = out[[i]]$It, x = out[[i]]$times, col = i)
}


# For a given beta_t, find the endemic bTB prevalence
#############################################################
endprev<- c(NA)
beta_t <- seq(0.0001, 0.001, 0.00001)  # make more finely spaced
for (i in 1:length(beta_t)){
	params <- c(b1= b1, b2 = b2, b3 = b3, b4= b4, b5 = b5, 
		b = b, r = r, K = K,
		muS = muS, muB = muB, muT = muT, muC = muC, muR = muR, muRC = muRC, 
		gamma = gamma, epsilon = epsilon,
		betaT = beta_t[i], betaB = betaB, betapT = beta_t[i], betapB = 3.92 * betaB)	
	sol <- as.data.frame(ode(x0, times, rhs, params))
	sol<- groom_sol(sol)
	endprev[i] <- sol$TBprev[101]
}
plot(x=beta_t, y=endprev, xlab = "Density dependent transmission rate for bTB", ylab = "bTB prevalence", type="l")
abline(h = 0.05, col = "red"); abline(h = 0.1, col = "red")
abline(h = 0.2, col = "red"); abline(h = 0.3, col = "red")
want= c(0.05, 0.1, 0.2, 0.3)  # want beta_t values associated with each of these end prevalences
approx(x = endprev, y = beta_t, xout = want)
# 0.0002890386 0.0003149813 0.0003764191 0.0004592954
legend("topleft", legend= c("0.00029", "0.00031", "0.00038", "0.00046"), bty= "n")


# Brucellosis only
#############################################################
betaT <- 0.0002909  # re-set just in case
beta_b <- c(seq(0.0001, 0.001, 0.0003), seq(0.001, 0.003, 0.001))
Ib0 = N - 100
It0 = 0
Ib0 = 100
Ic0 = 0
R0 = 0
Rc0 = 0
x0 = c(S = S0, It = It0, Ib = Ib0, Ic = Ic0, R = R0, Rc = Rc0)
times <- seq(0, 1000, 1)
out <- list(NA); outr <- list(NA)

par(mfrow=c(1, 2))
plot(NULL, xlim = c(0, 100), xlab= "Time (in years)", ylim= c(0, 500), 
	ylab = "Number of Ib and Rb animals")
legend("topright", title= "beta_b", legend=c("0.0001", "0.0004", "0.0007", "0.001"), 
	col = seq(1:length(beta_b)-2), lty = 1, bty = "n")
title(main = list("Fitness effects in Brucellosis Recoverds", cex = 0.7)) 
for(i in 1:5){
	params <- c(b1= b1, b2 = b2, b3 = b3, b4= b4, b5 = b5, 
		b = b, r = r, K = K,
		muS = muS, muB = muB, muT = muT, muC = muC, muR = muR, muRC = muRC, 
		gamma = gamma, epsilon = epsilon,
		betaT = betaT, betaB = beta_b[i], betapT = betapT, betapB = 3.92 * beta_b[i])	
	params_recov <- c(b1= b1, b2 = b2, b3 = 1, b4= b4, b5 = b1, 
		b = b, r = r, K = K,
		muS = muS, muB = muB, muT = muT, muC = muC, muR = muS, muRC = muT, 
		gamma = gamma, epsilon = epsilon,
		betaT = betaT, betaB = beta_b[i], betapT = betapT, betapB = 3.92 * beta_b[i])	
	sol <- as.data.frame(ode(x0, times, rhs, params))
	sol <- groom_sol(sol)
	sol_recov <- as.data.frame(ode(x0, times, rhs, params_recov))
	#out[[i]] <- groom_sol(sol)
	outr[[i]] <- groom_sol(sol_recov)
	lines(y = sol$Ib, x = sol$times, col = i)
	lines(y = sol$R, x = sol$times, col = i, lty = 2)
	rm(sol, sol_recov, params, params_recov)
}

plot(NULL, xlim = c(0, 100), xlab= "Time (in years)", ylim= c(0, 500), 
	ylab = "Number of animals")
title(main = list("No fitness effects in Brucellosis Recoverds", cex = 0.7)) 
for(i in 1:5){
	lines(y = outr[[i]]$Ib, x = outr[[i]]$times, col = i)
	lines(y = outr[[i]]$R, x = outr[[i]]$times, col = i, lty = 2)
}

################################################
endprev<- c(NA); endprev_recov <- c(NA)
beta_b <- seq(0.0005, 0.002, 0.00005)
for (i in 1:length(beta_b)){
	params <- c(b1= b1, b2 = b2, b3 = b3, b4= b4, b5 = b5, 
		b = b, r = r, K = K,
		muS = muS, muB = muB, muT = muT, muC = muC, muR = muR, muRC = muRC, 
		gamma = gamma, epsilon = epsilon,
		betaT = betaT, betaB = beta_b[i], betapT = betapT, betapB = 3.92 * beta_b[i])	
	params_recov <- c(b1= b1, b2 = b2, b3 = 1, b4= b4, b5 = b1, 
		b = b, r = r, K = K,
		muS = muS, muB = muB, muT = muT, muC = muC, muR = muS, muRC = muT, 
		gamma = gamma, epsilon = epsilon,
		betaT = betaT, betaB = beta_b[i], betapT = betapT, betapB = 3.92 * beta_b[i])	
	sol <- as.data.frame(ode(x0, times, rhs, params))
	sol<- groom_sol(sol)
	sol_recov <- as.data.frame(ode(x0, times, rhs, params_recov))
	sol_recov <- groom_sol(sol_recov)
	endprev[i] <- sol$Brucprev[101]
	endprev_recov[i] <- sol_recov$Brucprev[101]
}
par(mfrow = c(1,2))
plot(x=beta_b, y=endprev, xlab = "Density dependent transmission rate for brucellosis", ylab = "Brucellosis prevalence", type="l")
abline(h = 0.05, col = "red"); abline(h = 0.1, col = "red")
abline(h = 0.2, col = "red"); abline(h = 0.3, col = "red")
want= c(0.05, 0.1, 0.2, 0.3)  # want beta_t values associated with each of these end prevalences
approx(x = endprev, y = beta_b, xout = want)
# 0.0006592174 0.0007162451 0.0008537897 0.0010365370

plot(x = beta_b, y = endprev_recov, type = "l")
approx(x = endprev_recov, y = beta_b, xout = want)
# 0.0006200171 0.0006627512 0.0007553400 0.0008754510


# TB and Brucellosis, loops through transmission rate values for each set of hypotheses!
#############################################################
#bTB# 0.0002875412 0.0003147644 0.0003764948 0.0004593076
# bruc 0.0006200171 0.0006627512 0.0007553400 0.0008754510  # recov
#0.0006414904 0.0006830101 0.0007736913 0.0008886672 # no recov
K = 1000; N= 1000
# beta_b values that cover range giving 5, 10, 20, 30% prev
beta_b <- c(0.0005, 0.0006200171, 0.0006627512, 0.0007553400, 0.0008754510,
	0.0006414904, 0.0006830101, 0.0007736913, 0.0008886672,
	0.0009, 0.001, 0.0015, 0.002, 0.0025, 0.003) 
beta_t <- c(0.0001, 0.0002, 0.0002875412, 0.0003147644, 0.0003764948, 0.0004593076,
	seq(0.0005, 0.0015,0.0001)) 
# beta_t values giving 5, 10, 20, 30% prev: 0.0002728621 0.0002998335 0.0003412031 0.0003923223
Ib0 = N - 100
It0 = 50
Ib0 = 50
Ic0 = 0
R0 = 0
Rc0 = 0
x0 = c(S = S0, It = It0, Ib = Ib0, Ic = Ic0, R = R0, Rc = Rc0)
times <- seq(0, 100, 1)

# parameters set as mortality in recoverd = mortality in infected; sigma, epsilon= as above
out_allmort <- data.frame(betaB = rep(beta_b, each= length(beta_t)), 
	betaT = rep(beta_t, length(beta_b)), 
	N = NA, TBprev = NA, Brucprev = NA, BrucprecentR = NA, 
	prevTBinBrneg = NA, prevTBinCo = NA, prevBinTBneg= NA, prevBinCo = NA)

# parameters set as mortality in recoverd = mortality in susceptibles; sigma, epsilon= as above
out_redmort <- out_allmort

# change sigma/ep also later
val <- 1
for(i in 1:length(out_allmort[,1])){
	# for out_allmort
	params <- c(b1= b1, b2 = b2, b3 = b3, b4= b4, b5 = b5, 
		b = b, r = r, K = K,
		muS = muS, muB = muB, muT = muT, muR = muR, muRC = muRC, 
		gamma = gamma, epsilon = epsilon,
		betaT = out_allmort$betaT[i], betaB = out_allmort$betaB[i], 
		betapT = out_allmort$betaT[i], betapB = 3.92 * out_allmort$betaB[i])
			
	params_recov = c(b1= b1, b2 = b2, b3 = b3, b4= 1, b5 = b2, # assuming b4 =1, b5 = b2 
		b = b, r = r, K = K,
		muS = muS, muB = muB, muT = muT, muC = muC, muR = muS, muRC = muT, 
		gamma = gamma, epsilon = epsilon,
		betaT = out_redmort$betaT[i], betaB = out_redmort$betaB[i], 
		betapT = out_redmort$betaT[i], betapB = 3.92 * out_redmort$betaB[i])		
						
	sol <- as.data.frame(ode(x0, times, rhs, params))
	sol<- groom_sol(sol)
	d <- sol[101,]	
	out_allmort$N[i] <- d$N
	out_allmort$TBprev[i] <- d$TBprev 
	out_allmort$Brucprev[i] <-  d$Brucprev
	out_allmort$BrucprecentR[i] = d$R / d$N
	out_allmort$prevTBinBrneg[i] = d$It / (d$S + d$It)
	out_allmort$prevTBinCo[i] = (d$Ic + d$Rc)/ (d$Ib + d$Ic + d$R + d$Rc)
	out_allmort$prevBinTBneg[i]= (d$Ib+ d$R) / (d$S + d$Ib + d$R)
	out_allmort$prevBinCo[i] = (d$Ic + d$Rc) / (d$It + d$Ic + d$Rc)
	
	solr <- as.data.frame(ode(x0, times, rhs, params_recov))
	solr<- groom_sol(sol)
	d <- solr[101,]	
	out_redmort$N[i] <- d$N
	out_redmort$TBprev[i] <- d$TBprev 
	out_redmort$Brucprev[i] <-  d$Brucprev
	out_redmort$BrucprecentR[i] = d$R / d$N
	out_redmort$prevTBinBrneg[i] = d$It / (d$S + d$It)
	out_redmort$prevTBinCo[i] = (d$Ic + d$Rc)/ (d$Ib + d$Ic + d$R + d$Rc)
	out_redmort$prevBinTBneg[i]= (d$Ib+ d$R) / (d$S + d$Ib + d$R)
	out_redmort$prevBinCo[i] = (d$Ic + d$Rc) / (d$It + d$Ic + d$Rc)
	val <- val + 1
	rm(d, sol, solr)
}

############################################################
############################################################
############################################################
#bTB# 0.0002875412 0.0003147644 0.0003764948 0.0004593076
# bruc 0.0006200171 0.0006627512 0.0007553400 0.0008754510  # recov
#0.0006414904 0.0006830101 0.0007736913 0.0008886672 # no recov


############################################################
############################################################
############################################################
# plot- consequences of bTB invaison at default brucellosis prevalence at 5%
par(mfrow = c(2,2), mar =c(4, 4, 2, 1))
plot(y= rep(0.05, 5),  
	x = c(0.0001, 0.0002728621, 0.0002998335, 0.0003412031, 0.001), 
	xlab = "Transmission rate for bTB", 
	ylab = "Brucellosis prevalence", type = "l", col = 'red', 
	ylim = c(0, 1), xlim = c(0.0001, 0.001), cex.lab = 0.8, cex.axis = 0.8)
#df <- out_allmort[out_allmort$betaB == 0.0006414904,]
df <- out_redmort[out_redmort$betaB == 0.0006200171,]
lines(y = df$Brucprev, x = df$betaT, col = "blue")
points(x = df$betaT, y = df$prevBinTBneg, col = "light blue", pch = 19)
points(x = df$betaT, y = df$prevBinCo, col = "dark blue", pch = 19)
legend("topright", bty="n", 
	legend = c("before bTB", "after bTB, overall", "after bTB, prev in bTB -", 
		"after bTB, prev in bTB +"), 
	cex = 0.8, col = c("red", "blue", "light blue", "dark blue"),
	lty = c(1, 1, 0, 0), pch = c(NA, NA, 19, 19))

# plot- consequences of bTB invaison at default brucellosis prevalence at 10%
plot(y= rep(0.10, 5),  # default brucellosis prevalence at 10%
	x = c(0.0001, 0.0002728621, 0.0002998335, 0.0003412031, 0.001), 
	xlab = "Transmission rate for bTB", 
	ylab = "Brucellosis prevalence", type = "l", col = 'red', 
	ylim = c(0, 1), xlim = c(0.0001, 0.001), cex.lab = 0.8, cex.axis = 0.8)
#df <- out_allmort[out_allmort$betaB ==0.0006830101,]
df <- out_redmort[out_redmort$betaB == 0.0006627512,]
lines(y = df$Brucprev, x = df$betaT, col = "blue")
points(x = df$betaT, y = df$prevBinTBneg, col = "light blue", pch = 19)
points(x = df$betaT, y = df$prevBinCo, col = "dark blue", pch = 19)
legend("topright", bty="n", 
	legend = c("before bTB", "after bTB, overall", "after bTB, prev in bTB -", 
		"after bTB, prev in bTB +"), 
	cex = 0.8, col = c("red", "blue", "light blue", "dark blue"),
	lty = c(1, 1, 0, 0), pch = c(NA, NA, 19, 19))

# plot- consequences of bTB invaison at default brucellosis prevalence at 20%
plot(y= rep(0.20, 5),  # default brucellosis prevalence at 10%
	x = c(0.0001, 0.0002728621, 0.0002998335, 0.0003412031, 0.001), 
	xlab = "Transmission rate for bTB", 
	ylab = "Brucellosis prevalence", type = "l", col = 'red', 
	ylim = c(0, 1), xlim = c(0.0001, 0.001), cex.lab = 0.8, cex.axis = 0.8)
#df <- out_allmort[out_allmort$betaB == 0.0007736913,]
df <- out_redmort[out_redmort$betaB == 0.0007553400,]

lines(y = df$Brucprev, x = df$betaT, col = "blue")
points(x = df$betaT, y = df$prevBinTBneg, col = "light blue", pch = 19)
points(x = df$betaT, y = df$prevBinCo, col = "dark blue", pch = 19)
legend("topright", bty="n", 
	legend = c("before bTB", "after bTB, overall", "after bTB, prev in bTB -", 
		"after bTB, prev in bTB +"), 
	cex = 0.8, col = c("red", "blue", "light blue", "dark blue"),
	lty = c(1, 1, 0, 0), pch = c(NA, NA, 19, 19))

# plot- consequences of bTB invaison at default brucellosis prevalence at 30%
plot(y= rep(0.30, 5),  # default brucellosis prevalence at 10%
	x = c(0.0001, 0.0002728621, 0.0002998335, 0.0003412031, 0.001), 
	xlab = "Transmission rate for bTB", 
	ylab = "Brucellosis prevalence", type = "l", col = 'red', 
	ylim = c(0, 1), xlim = c(0.0001, 0.001), cex.lab = 0.8, cex.axis = 0.8)
#df <- out_allmort[out_allmort$betaB == 0.0008886672,]
df <- out_redmort[out_redmort$betaB == 0.0008754510,]
lines(y = df$Brucprev, x = df$betaT, col = "blue")
points(x = df$betaT, y = df$prevBinTBneg, col = "light blue", pch = 19)
points(x = df$betaT, y = df$prevBinCo, col = "dark blue", pch = 19)
legend("topright", bty="n", 
	legend = c("before bTB", "after bTB, overall", "after bTB, prev in bTB -", 
		"after bTB, prev in bTB +"), 
	cex = 0.8, col = c("red", "blue", "light blue", "dark blue"),
	lty = c(1, 1, 0, 0), pch = c(NA, NA, 19, 19))

######################################################



######################################################
NOT USED BELOW HERE
######################################################




make_three_plots = function(outlist, beta_t_val){
	# outlist = list of length(beta_b) = 6 containing ode solver output
	par(mfrow= c(1, 3))
	title <- paste("beta_t = ", as.character(beta_t_val))
	beta_b <- seq(0.05, 0.55, 0.1)
	
	# plot Brucellosis prevalence
	plot(NULL, xlim = c(0, 100), xlab= "Time (in years)", ylim= c(0, 0.5), 
		ylab = "Brucellosis prevalence", main = title)
	for (i in 1:length(beta_b)){
		lines(x = times, y = outlist[[i]]$Brucprev, col = i)
	}
	
	# plot TB prevalence
	plot(NULL, xlim = c(0, 100), xlab= "Time (in years)", ylim= c(0, 0.5), 
		ylab = "bTB prevalence", main = title)
	for (i in 1:length(beta_b)){
		lines(x = times, y = outlist[[i]]$TBprev, col = i)
	}

	# plot proportion brucellosis infections that are co-infected with bTB
	plot(NULL, xlim = c(0, 100), xlab= "Time (in years)", ylim= c(0, 0.5), 
		ylab = "Proportion Brucellosis with co-infection", main = title)
	for (i in 1:length(beta_b)){
		lines(x = times, y = outlist[[i]]$propBruc_co, col = i)
	}
	legend("topright", title= "beta_b", 
		legend=c("0.05", "0.15", "0.25", "0.35", "0.45", "0.55"), 
		col = seq(1:length(beta_t)), lty = 1, bty = "n")
}

# generates < 5%, ~20% and >30% bTB prevalence in bTB only model
make_three_plots(outlist = out[13:18], beta_t_val = 0.25)  # all die out, regardless of brucellosis transmission rate
make_three_plots(outlist = out[19:24], beta_t_val = 0.35)
make_three_plots(outlist = out[25:30], beta_t_val = 0.45)





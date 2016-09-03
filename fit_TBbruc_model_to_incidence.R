#############################################################
#############################################################
#############################################################
#############################################################
# Code to fit model (defined in run_stochastic_coinfection_model), to data incidence rate data
# 1- Sept- 2016
#############################################################
#############################################################
#############################################################
#############################################################
# Outline
#############################################################
# 1) Load fixed parameters and model
# 2) Load Data to fit to and data for evaluation
# 3) Impliment sto
# 4) Fit to Data by maximizing binomial likelihood, deterministic



#############################################################
#############################################################


#############################################################
#############################################################
1) Load fixed parameters, model
#############################################################
#############################################################
require("deSolve")
# get fixed.params & fixed.params.recov
source('~/GitHub/bTB-bruc-co-infection-ms/fixed_parameters_recovery.R', chdir = TRUE)
source('~/GitHub/bTB-bruc-co-infection-ms/fixed_parameters_norecovery.R', chdir = TRUE)
# rhs function, determinitic model, no agestr
source('~/GitHub/bTB-bruc-co-infection-ms/rhs.R', chdir = TRUE) 
# functions to run the stochastic co-infection model
source('~/GitHub/bTB-bruc-co-infection-ms/run_stochastic_coinfection_model.R', chdir = TRUE)
#############################################################
#############################################################
2) Load Data to fit and data for evaluation
#############################################################
#############################################################
# incidence rates, Lower Sabie (do not use...)
pS_TB <- c(0.0303, 0.28125, 0.074, 0.0435, 0.0555, 0.0667)
pS_B <- c(0, 0, 0.0363, 0.0390, 0.0740, 0.0000)
pB_TB <- c(0, 0, 0, 0.134, 0.095, 0.08667, 0.0139)
pTB_B <- c(0, 0.083, 0.091, 0.091, 0.083, 0.071, 0)
hist(c(pS_TB, pS_B, pB_TB, pTB_B))
# Total number avaliable to become infected
NS <- c(33, 32, 27, 25, 23, 18, 15)
nB <- c(15, 12, 11, 11, 12, 14, 15)
nTB <- c(2, 3, 12, 11, 10, 11, 9, 8)
# Number of new infections
x_S_TB <- c(0, 1, 9, 2, 2, 1, 1, 1) # total: 0, 1, 10, 3, 3, 2, 2, 1
x_B_TB <- c(0, 0, 0, 0, 1, 3, 1, 2)
x_S_B <- c(0, 0, 0, 1, 3, 1, 2, 0) # total brucellosis: 0, 0, 0, 1, 4, 4, 3, 2
x_TB_B <- c(0, 0, 1, 1, 1, 1, 1, 0)


# totals and new infections for Crocodile Bridge
NS_CB <- c(35, 31, 28, 26, 26, 23, 20)
nB_CB <- c(8, 10, 8, 10, 6, 9, 10)
nTB_CB <- c(3, 3, 10, 9, 10, 10, 12, 12)

x_S_TB_CB <- c(0, 0, 4, 1, 2, 1, 2, 1) #total TB: c(0, 0, 4, 1, 5, 2, 2, 1)
x_B_TB_CB <- c(0, 0, 0, 0, 3, 1, 0, 0)
x_S_B_CB  <- c(0, 3, 0, 4, 0, 3, 1, 1) #total Bruc: c(0, 4, 0, 4, 0, 3, 1, 1)
x_TB_B_CB <- c(0, 1, 0, 0, 0, 0, 0, 0)

# Assume the probability that the number of new bTB infections = 2 on first chunk
# use data on the number sampled; model predicted incidence rate... 
# ... calculate probability of observing the number of new infections observed
dbinom(x = 2, size = NS[1], prob = pS_TB[1], log = FALSE) 

# Evaluation data- Bootstrap prevalence estimates!
#############################################################
data <- read.csv("~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/		final_datasets_copied_from_phdfolder/cross_sectional_data_withdz_cleandisease_nofinal_Feb2016_capturetime_forsurv.csv")

overallbruc <- NA; overallbtb <- NA 
brucintbneg <- NA; brucintbpos <- NA
tbinbrucneg <- NA; tbinbrucpos <- NA
for (i in 1:1000){
	# sample one time point for each individual id (should be 151)
	ssdata <- ddply(data, .(id), function(id) {id[sample(nrow(id), size = 1),]})
	# calculate prevalence overall and infection specific
	overallbruc[i] <- length(ssdata$bruc[ssdata$bruc=="positive"])/ length(ssdata$bruc)
	overallbtb[i] <- length(ssdata$tb[ssdata$tb==1])/ length(ssdata$tb)
	brucintbneg[i] <- length(ssdata$tb[ssdata$bruc=="positive" & ssdata$tb == 0]) / 
		length(ssdata$tb[ssdata$tb==0])
	brucintbpos[i] <- length(ssdata$tb[ssdata$bruc=="positive" & ssdata$tb == 1]) / 
		length(ssdata$tb[ssdata$tb==1])
	tbinbrucneg[i] <- length(ssdata$tb[ssdata$bruc=="negative" & ssdata$tb == 1]) / 
		length(ssdata$tb[ssdata$bruc=="negative"])
	tbinbrucpos[i] <- length(ssdata$tb[ssdata$bruc=="positive" & ssdata$tb == 1]) / 
		length(ssdata$tb[ssdata$bruc=="positive"])
}
quantile(overallbruc, c(0.025, 0.25, 0.5, 0.75, 0.975)) 
# 0.3112583 0.3311258 0.3443709 0.3576159 0.3774834 
quantile(overallbtb, c(0.025, 0.25, 0.5, 0.75, 0.975))
# 0.2317881 0.2582781 0.2715232 0.2847682 0.3112583
quantile(brucintbneg, c(0.025, 0.25, 0.5, 0.75, 0.975))
# 0.2654405 0.2892206 0.3035714 0.3181818 0.3451408
quantile(brucintbpos, c(0.025, 0.25, 0.5, 0.75, 0.975))
# 0.3657982 0.4222222 0.4523810 0.4864865 0.5405985
quantile(tbinbrucneg, c(0.025, 0.25, 0.5, 0.75, 0.975))
# 0.1781895 0.2103947 0.2268041 0.2448980 0.2783505
quantile(tbinbrucpos, c(0.025, 0.25, 0.5, 0.75, 0.975))
# 0.2857143 0.3333333 0.3584906 0.3818182 0.4313725



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
#############################################################
3) Impliment Stochastic Version
#############################################################
#############################################################
# no recovery
params.fixed = c(fixed.params, gamma=1/2)
params <- c(params.fixed, betaB = 0.003, betapB = 0.009, betaT = 0.0006, betapT = 0.0006000)

# recovery assumption
params.recov <- c(fixed.params.recov, gamma = 1/2, betaB = 0.003, betapB = 0.009, betaT = 0.0006, betapT = 0.0006000)

nsims <- 10	
nstep <- 10000
data<- run_stochastic_coinfection_model(params, nstep, nsims)
data.recov <- run_stochastic_coinfection_model(params.recov, nstep, nsims)
#betaB=0.001, betapB=0.003,betaT=0.0006,betapT=0.0006000)-> bruc dies out 
#betaB=0.003, betapB=0.009,betaT=0.0006,betapT=0.0006000)-> tb dies out
make_stochastic_plots = function(data){
	par(mfrow = c(2, 2))
	# It
	plot(It~ cumtime, data= data[[1]], 
		xlab = "Time", ylab = "It", type = 'o', cex = 0.3)
	for (k in 1:nsims){
		lines(It~cumtime, data = data[[k]], col = k, 
		type = 'o', cex = 0.3)
	}
	plot(Ib~ cumtime, data= data[[1]], 
		xlab = "Time", ylab = "Ib", type = 'o', cex = 0.3)
	for (k in 1:nsims){
		lines(Ib~cumtime, data = data[[k]], col = k,
		type = 'o', cex = 0.3)
	}
	plot(R~ cumtime, data= data[[1]], 
		xlab = "Time", ylab = "Rb", type = 'o', cex = 0.3)
	for (k in 1:nsims){
		lines(R~cumtime, data = data[[k]], col = k, 
		type = 'o', cex = 0.3)
	}
	plot(Ic~ cumtime, data= data[[1]],
		xlab = "Time", ylab = "Ic", type = 'o', cex = 0.3)
	for (k in 1:nsims){
		lines(Ic~cumtime, data = data[[k]], col = k, 
		type = 'o', cex = 0.3)
	}
}
make_stochastic_plots(data)
make_stochastic_plots(data.recov)  # both persist!
#############################################################
#############################################################
# 4) Fit to Data by maximizing binomial likelihood, stochastic
#############################################################
#############################################################
# ASSUMTION 1: FIXED RECOVERY (2yrs) & REACTIVATION; muR = muI; muRC = muC
#############################################################
ML.sir.norecovery = function(betavals, params.fixed, data){
	##############################################
	# INPUT: 
	# betavals = c(betat, betab)
	# params.fixed = all other parameters
	# data = number of new cases per capture (0.5/yr)= a list of length ?
	# Returns: NLL value to be optimized. 
    ##############################################

	# model, expected incidence rates (per 0.5 year)
	times <- seq(0, 500, 0.5)  

	Y0 <- c(S = 950, It = 20, Ib = 20, Ic = 0, R = 0, Rc = 10)
    params0 = c(betaT = betavals[1], betapT = betavals[1], 
    		betaB = betavals[2], betapB = 3.92 * betavals[2])
	sol <- as.data.frame(ode(Y0, times, rhs, parms= c(params0, params.fixed)))
	
	# Data
    obs_num_new_tb <- data$tb  # calculated as change in n frm t -1 to t
    obs_num_new_bruc <- data$br
    # obs_num_new_tb_co <- data$tbco
    # obs_num_new_bruc_co <- data$brucco



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


params.fixed = c(fixed.params, gamma=1/2)
times <- seq(0, 500, 0.5)  # now want data in half year time intervals.



# ASSUMTION 2: FIXED RECOVERY, REACTIVATION RATE; muR = muR; muRC = muRC
#############################################################



# ASSUMTION 3: ESTIMATE RECOVERY, REACTIVATION RATE; muR = muS; muRC = muT
#############################################################






# CUT STUFF: 
#############################################################
#############################################################
# Data- assume binomially distributed likelihood
# Data plots in prims
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
agedata <- read.csv("~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/final_datasets_copied_from_phdfolder/cross_sectional_data_withdz_cleandisease_nofinal_Feb2016_capturetime.csv")
LSage <- agedata$age_yr[agedata$herdorig == "LS" & agedata$capturetime == 0]
CBage <- agedata$age_yr[agedata$herdorig == "CB" & agedata$capturetime == 3]

LS<- agedata[agedata$herdorig=="LS",]
LSmass <- LS[LS$capturetime==0,]
LSmass[,c(3, 22, 30)]

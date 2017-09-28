##########################################################
##########################################################
##########################################################
# Model sensitivity and evaluation 
##########################################################
##########################################################
##########################################################
# Outline
# 1) Run with other form of density dependence (logistic)
# 2) Do whole process with longer infection duraiton (1/gamma = 3 yr; 4 yr)
# 2) Latin hypercube PRCC
##########################################################
library('sensitivity')
library('deSolve')
library('lhs')


##########################################################
##########################################################
# 1) Run with other form of density dependence (logistic)
##########################################################
##########################################################
rm(list = ls())
source('~/GitHub/bTB-bruc-co-infection-ms/fixed_parameters_norecovery_agematrix.R', chdir = TRUE)
source('~/GitHub/bTB-bruc-co-infection-ms/rhs_age.R', chdir = TRUE)

# Run B-H model with no disease (endemic population size = 609)
s_index <- 1:20
it_index <- 21:40
ib_index <- 41:60
ic_index <- 61:80
r_index <- 81:100
rc_index <- 101:120
relageall = c(0.137, rep(0.368/4, 4), rep(0.185/4, 4),  
	rep(0.235/6, 6), rep(0.075/5, 5))					
relage = relageall
S0 = 400*relage; It0 = 0*relage; Ib0 = 0*relage; 
Ic0 = 0*relage; R0 = 0 * relage; Rc0 = 0 * relage
x0 = c(S0, It0, Ib0, Ic0, R0, Rc0)
times <- seq(0, 1000, 1)  
params.test <- c(fixed.params, list(gamma=1/2, theta = 4, K = 433,
	betaB = 0.6087396, betaT = 1.2974554/1000, rhoT = 1, rhoB = 2.1))
sol <- as.data.frame(ode(x0, times, rhs_age_matrix, params.test))
stable_age <- unname(unlist( sol[1000, c(2:21)]/sum(sol[500, c(2:21)]) ))

plot_raw_numbers = function(sol){
	plot(sol$time, apply(sol[s_index+1], 1, sum), col= "black",
		type= 'l', ylim = c(0, 800), ylab = "Number of animals", 
		xlab = "Time (in years)")
	lines(sol$time, apply(sol[it_index+1], 1, sum), col= "red")
	lines(sol$time, apply(sol[ib_index+1], 1, sum), col= "blue")
	lines(sol$time, apply(sol[ic_index+1], 1, sum), col= "green")
	lines(sol$time, apply(sol[r_index+1], 1, sum), col = "orange")
	lines(sol$time, apply(sol[rc_index+1], 1, sum), col = "pink")
	legend("topright", legend = c("S", "It", "Ib", "Ic", "R", "Rc"),
		col = c("black", "red", "blue", "green", "orange", "pink"), 
		bty="n", lty = 1)
}
plot_raw_numbers(sol)
sum(sol[1001,c(2:length(sol))])

# Set K s.t. get same popuation size/age structure (endemic population size = 609.0097)
K <- seq(1033.2, 1033.35, 0.01)
S0 = 400*relage; It0 = 0*relage; Ib0 = 0*relage; 
Ic0 = 0*relage; R0 = 0 * relage; Rc0 = 0 * relage
x0 = c(S0, It0, Ib0, Ic0, R0, Rc0)
times <- seq(0, 1000, 1)  

# test (run with bruc only ,hten introduce bTB)
params.test <- c(fixed.params, list(gamma=1/2, theta = 4, K = 1033.29,
		betaB = 0.6087396, betaT = 1.2974554/1000, rhoT = 1, rhoB = 2.1))
S0 = 398*relage; It0 = 0*relage; Ib0 = 2*relage; 
Ic0 = 0*relage; R0 = 0 * relage; Rc0 = 0 * relage
x0 = c(S0, It0, Ib0, Ic0, R0, Rc0)
sol <- as.data.frame(ode(x0, times,  rhs_age_logistic, params.test))
x0new <- unlist(sol[length(sol[,1]), 2:length(sol)])
x0new[27] <- 2
sol <- as.data.frame(ode(x0new, times,  rhs_age_logistic, params.test))
plot_raw_numbers(sol)

val <- NA
for (i in 1:length(K)){
	params.test <- c(fixed.params, list(gamma=1/2, theta = 4, K = K[i],
		betaB = 0.6087396, betaT = 1.2974554/1000, rhoT = 1, rhoB = 2.1))
	sol <- as.data.frame(ode(x0, times,  rhs_age_logistic, params.test))
	val[i] <- sum(sol[1001,c(2:length(sol))])	
}
plot(x = K, y = val, type = "b", pch = 19)
abline(h = 609.0097) # K = 1033.29

val <- NA
K <- seq(382.5, 382.56, 0.002)
for (i in 1:length(K)){
	params.test <- c(fixed.params, list(gamma=1/2, theta = 4, K = K[i],
		betaB = 0.6087396, betaT = 1.2974554/1000, rhoT = 1, rhoB = 2.1))
	sol <- as.data.frame(ode(x0, times,  rhs_age_matrix_ricker, params.test))
	val[i] <- sum(sol[1001,c(2:length(sol))])	
}
plot(x = K, y = val[1:length(K)], type = "b", pch = 19)
abline(h = 609.0097) # K = 382.556 

# Make plots
###############################
n = 500
source('~/GitHub/bTB-bruc-co-infection-ms/get_EE.R', chdir = TRUE)
set.seed(1)

EE_bTB_single_list <- NA; EE_bTB_co_list <- NA
EE_brucellosis_single_list <- NA; EE_brucellosis_co <- NA
rEE_bTB_single_list <- NA; rEE_bTB_co_list <- NA
rEE_brucellosis_single_list <- NA; rEE_brucellosis_co <- NA
lEE_bTB_single_list <- NA; lEE_bTB_co_list <- NA
lEE_brucellosis_single_list <- NA; lEE_brucellosis_co <- NA

rhoB <- NA; mortB <- NA; mortT <- NA; mortC <- NA
pb <- txtProgressBar(min = 0, max = n, style = 3, char = "=")

for(i in 1:n){
	if(i %% 5 == 0){
		setTxtProgressBar(pb, i)
	}
	params <- c(fixed.params, list(gamma=1/2, theta = 4, K = 433,
		betaB = 0.6087396, betaT = 1.2974554/1000, rhoT = 1, rhoB = 2.1))
	params$rhoB<- exp(rnorm(n = 1, mean = 0.75579, sd = 0.40714))
	rhoB[i] <- params$rhoB
	mortB[i] <- exp(rnorm(n = 1, mean = 1.1060, sd = 0.3505))
	mortT[i] <- exp(rnorm(n = 1, mean = 1.0370, sd = 0.3483 ))
	mortC[i] <- exp(rnorm(n = 1, mean = 2.1430, sd = 0.5004329))
	params$muB <- params$muS * mortB[i]
	params$muT <- params$muS * mortT[i]
	params$muC <- params$muS * mortC[i]
	params$muB[params$muB > 1] <- 1
	params$muT[params$muT > 1] <- 1
	params$muC[params$muC > 1] <- 1
	params$muR <- params$muT
	params$muRC <- params$muC
	val <- getEE(params, method = "beverton-holt")

	paramslog <- params; paramslog$K <- 1033.29
	vallog <- getEE(paramslog, method = "logistic")
	
	paramsricker <- params; paramsricker$K <- 382.556
	valricker <- getEE(paramsricker, method = "ricker")

	EE_bTB_single_list[i] <- val[1]; EE_bTB_co_list[i] <- val[2]
	EE_brucellosis_single_list[i] <- val[3]; 	EE_brucellosis_co[i] <- val[4]
	
	rEE_bTB_single_list[i] <- valricker[1]; rEE_bTB_co_list[i] <- valricker[2]
	rEE_brucellosis_single_list[i] <- valricker[3]; rEE_brucellosis_co[i] <- valricker[4]
	lEE_bTB_single_list[i] <- vallog[1]; lEE_bTB_co_list[i] <- vallog[2]
	lEE_brucellosis_single_list[i] <- vallog[3]; lEE_brucellosis_co[i] <- vallog[4]
}

senslist <- list(
	EE_bTB_single_list = EE_bTB_single_list, 
	EE_bTB_co_list = EE_bTB_co_list,
	EE_brucellosis_single_list = EE_brucellosis_single_list, 
	EE_brucellosis_co = EE_brucellosis_co,
	rEE_bTB_single_list = rEE_bTB_single_list, 
	rEE_bTB_co_list = rEE_bTB_co_list,
	rEE_brucellosis_single_list = rEE_brucellosis_single_list, 
	rEE_brucellosis_co = rEE_brucellosis_co,	
	lEE_bTB_single_list = lEE_bTB_single_list, 
	lEE_bTB_co_list = lEE_bTB_co_list,
	lEE_brucellosis_single_list = lEE_brucellosis_single_list, 
	lEE_brucellosis_co = lEE_brucellosis_co)
str(senslist)

saveRDS(senslist, "~/GitHub/bTB-bruc-co-infection-ms/sensitivity_densitydependence_simulation_results.rds")

##########################################################
##########################################################
# 3) LHS Sensitivity
##########################################################
##########################################################
rm(list = ls())

source('~/GitHub/bTB-bruc-co-infection-ms/fixed_parameters_norecovery_agematrix.R', chdir = TRUE)
source('~/GitHub/bTB-bruc-co-infection-ms/rhs_age.R', chdir = TRUE)
source('~/GitHub/bTB-bruc-co-infection-ms/get_EE.R', chdir = TRUE)
source('~/GitHub/bTB-bruc-co-infection-ms/get_Ro.R', chdir = TRUE)

# draw 100 lhs samples from 13 parameters (13 rows, 100 columns)
##########################################################
set.seed(1)
X <- randomLHS(13, 100)  
params <- c(fixed.params, list(gamma=1/2, betaB = 0.6087396,
	betaT = 1.2974554/1000, rhoT = 1, rhoB = 2.1, theta= 4, K = 433))

df <- data.frame(
	gamma = qunif(X[1, ], min = params$gamma/2, max = params$gamma*2),
	betaB = qunif(X[2, ], min = params$betaB/2, params$betaB*2),
	betaT = qunif(X[3, ], min = params$betaT/2, params$betaT*2), 
	epsilon = qunif(X[4, ], min = params$epsilon/2, params$epsilon*2), 
	rhoT = qunif(X[5, ], min = params$rhoT/2, params$rhoT*2), 
	rhoB = qunif(X[6, ], min = params$rhoB/2, params$rhoB*2), 
	theta =  qunif(X[7, ], min = params$theta/2, params$theta*2), 
	K = qunif(X[8, ], min = params$K/2, params$K*2), 
	pchange_muB = qunif(X[9, ], min = 3.02/2, max = 3.02*2), 
	pchange_muT = qunif(X[10, ], min = 2.82/2, max = 2.82*2), 
	pchange_muC = qunif(X[11, ], min = 8.58/2, max = 8.58),
	pchange_muR = qunif(X[12, ], min = 3.02/2, max = 3.02*2), 
	pchange_muRC = qunif(X[13, ], min = 8.58/2, max = 8.58), 
	# outputs of interest
	RoTB = NA, RoTBco = NA, change_RoTB = NA, 
	RoB = NA, RoBco = NA, change_RoB = NA, 
	prevTB = NA, prevTBco = NA, 	change_prevTB = NA, 
	prevB = NA, prevBco = NA, 	change_prevB = NA)
		
# run base model with no disease
##########################################################
s_index <- 1:20
it_index <- 21:40
ib_index <- 41:60
ic_index <- 61:80
r_index <- 81:100
rc_index <- 101:120
relageall = c(0.137, rep(0.368/4, 4), rep(0.185/4, 4),  
	rep(0.235/6, 6), rep(0.075/5, 5))					
relage = relageall
S0 = 400*relage; It0 = 0*relage; Ib0 = 0*relage; 
Ic0 = 0*relage; R0 = 0 * relage; Rc0 = 0 * relage
x0 = c(S0, It0, Ib0, Ic0, R0, Rc0)
times <- seq(0, 1000, 1)  
sol <- as.data.frame(ode(x0, times, rhs_age_matrix, params))
stable_age <- unname(unlist( sol[1000, c(2:21)]/sum(sol[500, c(2:21)]) ))

##########################################################
pb <- txtProgressBar(min = 0, max = length(df[,1]), style = 3, char = "=")
for (i in 1:length(df[ ,1])){
	if(i %% 5 == 0){
		setTxtProgressBar(pb, i)
	}
	params <- c(fixed.params, list(gamma=1/2, betaB = 0.6087396,
		betaT = 1.2974554/1000, rhoT = 1, rhoB = 2.1, theta= 4, K = 433))
	params$muB <- params$muS * df$pchange_muB[i]
	params$muT <- params$muS * df$pchange_muT[i]
	params$muC <- params$muS * df$pchange_muC[i]
	params$muR <- params$muS * df$pchange_muR[i]
	params$muRC <- params$muS * df$pchange_muRC[i]
	params$K <- df$K[i]
	params$theta <- df$theta[i]
	params$rhoT <- df$rhoT[i]
	params$rhoB <- df$rhoB[i]
	params$epsilon <- df$epsilon[i]
	params$betaT <- df$betaT[i]
	params$betaB <- df$betaB[i]
	params$gamma <- df$gamma[i]

	# check no mortality rates > 1
	params$muB[params$muB > 1] <- 1
	params$muT[params$muT > 1] <- 1
	params$muC[params$muC > 1] <- 1
	params$muR[params$muR > 1] <- 1
	params$muRC[params$muRC > 1] <- 1

	# fill in solutions for Ro
	df$RoTB[i] <- Ro_bTB_single_age(params)
	df$RoTBco[i] <- Ro_bTB_co_age(params)
	df$RoB[i] <- Ro_brucellosis_single_age(params)
	df$RoBco[i] <- Ro_brucellosis_co_age(params)

	# fill in solutions for endemic prevalence
	val <- getEE(params, method = "beverton-holt")
	df$prevTB[i] <- val[1]
	df$prevTBco[i] <- val[2]
	df$change_prevTB[i] <- df$prevTBco[i] - df$prevTB[i]
	df$prevB[i] <- val[3]
	df$prevBco[i] <- val[4]
	df$change_prevB[i] <- df$prevBco[i] - df$prevB[i]
	rm(val); rm(params)
}

saveRDS(df, "~/GitHub/bTB-bruc-co-infection-ms/LHS_sensitivity_results.rds")

bonferroni.alpha <- 0.05/13

# Ro TB
prcc_RoTB <- pcc(df[ , 1:13], df[ ,15], nboot = 1000, rank = TRUE, conf = 1- bonferroni.alpha)
prcc_RoB <- pcc(df[ , 1:13], df[ ,18], nboot = 1000, rank = TRUE, conf = 1- bonferroni.alpha)
prcc_EETB <- pcc(df[ , 1:13], df[ ,21], nboot = 1000, rank = TRUE, conf = 1- bonferroni.alpha)
prcc_EEB <- pcc(df[ , 1:13], df[ ,24], nboot = 1000, rank = TRUE, conf = 1- bonferroni.alpha)

pt <- print(prcc_RoTB)
pb <- print(prcc_RoB)
dfRo <- data.frame(
	param = rep(colnames(df)[1:13], 2),
	infection = c(rep("TB", 13), rep("brucellosis", 13)),
	Ro = c(pt$original, pb$original),
	cilow = c(pt[ , 4], pb[ , 4]),
	cihigh = c(pt[ , 5], pb[ , 5])
)

pt <- print(prcc_EETB)
pb <- print(prcc_EEB)
dfEE <- data.frame(
	param = rep(colnames(df)[1:13], 2),
	infection = c(rep("TB", 13), rep("brucellosis", 13)),
	EE = c(pt$original, pb$original),
	cilow = c(pt[ , 4], pb[ , 4]),
	cihigh = c(pt[ , 5], pb[ , 5])
)

write.csv(dfRo, "~/GitHub/bTB-bruc-co-infection-ms/Ro_sensitivity.rds")
write.csv(dfEE, "~/GitHub/bTB-bruc-co-infection-ms/EE_sensitivity.rds")

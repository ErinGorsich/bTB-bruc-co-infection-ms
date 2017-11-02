##########################################################
##########################################################
##########################################################
# Erin Gorsich
# This code conducts sensitivity analsyes for the co-infection model
##########################################################
##########################################################
##########################################################
# Outline
# 1) Run with other form of density dependence (logistic, ricker)
# 2) Do whole process with longer infection duraiton (1/gamma = 3 yr; 4 yr)
# 2) Latin hypercube, PRCC sensitivity
##########################################################
##########################################################
library('sensitivity')
library('deSolve')
library('lhs')
library("doParallel")
library("foreach")

setwd("~/GitHub/bTB-bruc-co-infection-ms/pde")

##########################################################
##########################################################
# 1) Run with other form of density dependence (logistic)
##########################################################
##########################################################
rm(list = ls())
source('fixed_parameters.R', chdir = TRUE)
source('rhs.R', chdir = TRUE)

# age divisions in rhs function
agemax <- 20
agestep <- 0.1
N <- agemax / agestep
ages <- seq(1, agemax + 1, by = agestep)[-(N)]
N == length(ages)

# generate parameters with correct agebins
f.params <- gen_fixed_params(agemax, agestep,
    p = p, recovery = FALSE)
f.params.recov <- gen_fixed_params(agemax, agestep, 
    p = p, recovery = TRUE)

source('plotting_functions.R', chdir = TRUE)

# Starting agestructure (Jolles 2007; Caron et al. 2001)
juv <- rep(0.137 / length(ages[ages < 2]), length(ages[ages < 2]))
sa <- rep(0.368 / length(ages[ages >= 2 & ages < 6]), 
    length(ages[ages >= 2 & ages < 6]))
a <- rep(0.185 / length(ages[ages >= 6 & ages < 9]), 
    length(ages[ages >= 6 & ages < 9]))
ma <- rep(0.235 / length(ages[ages >= 9 & ages < 14]), 
    length(ages[ages >= 9 & ages < 14]))
sen <- rep(0.075 / length(ages[ages >= 14 ]), length(ages[ages >= 14]))

relage <- c(juv, sa, a, ma, sen); length(relage) == N									
plot.relage <- c(0.137, rep(0.368/4, 4), rep(0.185/3, 3), 
    rep(0.235/5, 5), rep(0.075/7, 7))


# Run B-H model with no disease (endemic population size = 559.8903)
#############################################################
S0 <- 400 * relage; It0 <- 0 * relage; Ib0 <- 0 * relage
Ic0 <- 0 * relage; R0 <- 0 * relage; Rc0 <- 0 * relage
times <- seq(1, 300, 1)
x0 <- c(S0, It0, Ib0, Ic0, R0, Rc0)
params <- c(f.params, list(gamma = 1/2, betaB = 0.5764065, 
	betaT = 1.3305462/1000, rhoT = 1, rhoB = 2.1))
test <- as.data.frame(ode.1D(x0, times, rhs, params, 
	nspec = 6, dimens = N, method = "ode45"))
x0 <- unname(unlist(test[length(test[,1]), c(2:(length(colnames(test))))])) 
sum(x0)
plot_agestructure(test, t = 300)

# Set K for logistic to get the same popuation size/age structure 
# (endemic population size = 559.8903)
#############################################################
#K <- seq(1033.2, 1033.35, 0.01)
times <- seq(0, 300, 1)  

# test run for logistic model 
# at K = 1033-> 380; K = 2000 -> 736; K = 1500 -> 552.4; 1500 -> 589
# at 1550 --> 570.81 # at 1525 -> 561
params.test <- params
params.test$K <- 1520 
sol <- as.data.frame(ode.1D(x0, times, rhs_logistic, params.test, 
	nspec = 6, dimens = N, method = "ode45"))
plot_raw_numbers(sol)
sum(sol[length(sol[,1]), c(2:(length(colnames(sol))))])
plot_agestructure(sol, t = 200)

K <- seq(1520.2, 1520.4, 0.05)
val <- NA
params <- c(f.params, list(gamma = 1/2, betaB = 0.5764065, 
	betaT = 1.3305462305462/1000, rhoT = 1, rhoB = 2.1))
for (i in 1:length(K)){
	params.test <- params
	params.test$K <- K[i]
	sol <-  as.data.frame(ode.1D(x0, times, rhs_logistic, params.test, 
		nspec = 6, dimens = N, method = "ode45"))
	val[i] <- sum(sol[length(sol[,1]),c(2:length(colnames(sol)))])	
	rm(params.test, sol)
}
plot(x = K, y = val, type = "b", pch = 19)
abline(h = 559.8903) 
# K = 1520.35


# Set K for ricker to get the same popuation size (size = 559.8903)
#############################################################
# old -> K <- seq(382.5, 382.56, 0.002)
# K = 382 -> 509 K = 400 -> 533.5; K = 450-> 500
# 410 -> 546; 415-> 553; 420 -> 560
params.test <- params
params.test$K <- 417 
sol <- as.data.frame(ode.1D(x0, times, rhs_ricker, params.test, 
	nspec = 6, dimens = N, method = "ode45"))
plot_raw_numbers(sol)
sum(sol[length(sol[,1]), c(2:(length(colnames(sol))))])
plot_agestructure(sol, t = 200)

#K2 <- seq(417, 420, 0.2)
K2 <- seq(419, 420, 0.05)
val2 <- NA
params <- c(f.params, list(gamma = 1/2, betaB = 0.5764065, 
	betaT = 1.3305462305462/1000, rhoT = 1, rhoB = 2.1))
for (i in 1:length(K2)){
	params.test <- params
	params.test$K <- K2[i]
	sol <-  as.data.frame(ode.1D(x0, times, rhs_ricker, params.test, 
		nspec = 6, dimens = N, method = "ode45"))
	val2[i] <- sum(sol[length(sol[,1]),c(2:length(colnames(sol)))])	
}
plot(x = K2, y = val2, type = "b", pch = 19)
abline(h = 559.8903) # K = 419.75


# Make plots
#############################################################
set.seed(1)
n = 500
source('get_EE.R', chdir = TRUE)

# Generate 500 samples
cl <- makeCluster(6)
# cl <- makeCluster(10) # work computer
registerDoParallel(cl)

densitydependence <- foreach(icount(n), .combine = rbind, .packages = "deSolve") %dopar% {
	params <- c(f.params, list(gamma = 1/2, betaB = 0.5764065, 
		betaT = 1.3305462/1000, rhoT = 1, rhoB = 2.1))
	params$rhoB<- exp(rnorm(n = 1, mean = 0.75579, sd = 0.40714))
	B <- rnorm(n = 1, mean = 1.1060, sd = 0.3505)
	TB <- rnorm(n = 1, mean = 1.0370, sd = 0.3483)
	params$muB <- params$muS * exp(B)
	params$muT <- params$muS * exp(TB)
	params$muC <- params$muS * exp(B + TB)
	params$muR <- params$muT
	params$muRC <- params$muC

	val <- get_EE(params, x0, method = "beverton-holt")

	paramslog <- params; paramslog$K <- 1520.35
	vallog <- get_EE(paramslog, x0, method = "logistic")
	
	paramsricker <- params; paramsricker$K <- 382.556
	valricker <- get_EE(paramsricker, x0, method = "ricker")

	data <- data.frame(
		rhoB = params$rhoB,
		dB = B,
		dT = TB,
		EE_bTB_single = val[1], EE_bTB_co = val[2], 
		EE_brucellosis_single = val[3], EE_brucellosis_co = val[4], 
		rEE_bTB_single = valricker[1], rEE_bTB_co = valricker[2], 
		rEE_brucellosis_single = valricker[3], rEE_brucellosis_co = valricker[4], 
		lEE_bTB_single = vallog[1], lEE_bTB_co = vallog[2], 
		lEE_brucellosis_single = vallog[3], lEE_brucellosis_co = vallog[4])
}
stopCluster(cl)

str(densitydependence)

saveRDS(densitydependence, "sensitivity_densitydependence_results.rds")

##########################################################
##########################################################
# 3) LHS Sensitivity
##########################################################
##########################################################
source('get_Ro.R', chdir = TRUE)
source('get_EE.R', chdir = TRUE)

# draw 100 lhs samples from 13 parameters (13 rows, 100 columns)
##########################################################
set.seed(1)
X <- randomLHS(13, 100)  
params <- c(f.params, list(gamma = 1/2, betaB = 0.5764065, 
	betaT = 1.3305462/1000, rhoT = 1, rhoB = 2.1))

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
	pchange_muRC = qunif(X[13, ], min = 8.58/2, max = 8.58))

cl <- makeCluster(6)
# cl <- makeCluster(10)
registerDoParallel(cl)
lhs <- foreach(d = iter(df, by = "row"), .combine = rbind, .packages = "deSolve") %dopar% {
	params <- c(f.params, list(gamma = 1/2, betaB = 0.5764065, 
		betaT = 1.3305462/1000, rhoT = d$rhoT, rhoB = 2.1))

	# update parameters to those specified in LHS sample
	params$muB <- params$muS * d$pchange_muB
	params$muT <- params$muS * d$pchange_muT
	params$muC <- params$muS * d$pchange_muC
	params$muR <- params$muS * d$pchange_muR
	params$muRC <- params$muS * d$pchange_muRC
	params$K <- d$K
	params$theta <- d$theta
	params$rhoT <- d$rhoT
	params$rhoB <- d$rhoB
	params$epsilon <- d$epsilon
	params$betaT <- d$betaT
	params$betaB <- d$betaB
	params$gamma <- d$gamma

	val <- get_EE_lhs(params, x0, method = "beverton-holt")

	# set xB and xT (endemic conditions when only )
	xB <- val$xB
	xT <- val$xT
	
	RoTB <- Ro_bTB_single(params, x0)
	RoTBco <- Ro_bTB_co(params, xB)
	RoB <- Ro_brucellosis_single(params, x0)
	RoBco <- Ro_brucellosis_co(params, xT)

	# fill in solutions
	data <- data.frame(
		RoTB = RoTB,
		RoTBco = RoTBco,
		change_RoTB = RoTBco - RoTB,
		RoB = RoB,
		RoBco = RoBco,
		change_RoB = RoBco - RoB,
		prevTB = val$EE[1],
		prevTBco = val$EE[2],
		change_prevTB = val$EE[2] - val$EE[1],
		prevB = val$EE[3],
		prevBco = val$EE[4],
		change_prevB = val$EE[4] - val$EE[3])
}
stopCluster(cl)

summary(lhs)
sensitivity <- as.data.frame(cbind(df, lhs))

saveRDS(sensitivity, "LHS_sensitivity_results.rds")

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

write.csv(dfRo, "Ro_sensitivity.rds")
write.csv(dfEE, "EE_sensitivity.rds")

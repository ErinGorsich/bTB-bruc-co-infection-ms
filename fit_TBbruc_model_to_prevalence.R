#############################################################
#############################################################
#############################################################
#############################################################
# Code to fit model to prevalence data (BETAB = 3.8* BETAB)
# 27- Sept- 2016
#############################################################
#############################################################
#############################################################
#############################################################
# Outline
#############################################################
# 1) Load fixed parameters and model; diffeqs; stochastic version
# 2) Load Data to fit to and data for evaluation
# 3) Fit Deterministic Version to overall prevalence
### Four models fit, one for each assumption about recovery, mortality
# 4) Fit Determinisit Version to overall prevlalence after systematically varying epsilon, gamma
# 5) Plot deterministic version of model predictions
# 6) Plot stochastic version of model predictions
#############################################################
#############################################################


#############################################################
#############################################################
1) Load fixed parameters, model
#############################################################
#############################################################
require("deSolve")
library("plyr")
library("ggplot2")
set.seed(5)
# get fixed.params & fixed.params.recov
source('~/GitHub/bTB-bruc-co-infection-ms/fixed_parameters_recovery.R', chdir = TRUE)
source('~/GitHub/bTB-bruc-co-infection-ms/fixed_parameters_norecovery.R', chdir = TRUE)

# rhs function, determinitic model, no agestr
source('~/GitHub/bTB-bruc-co-infection-ms/rhs.R', chdir = TRUE) 
# functions to run the stochastic co-infection model
source('~/GitHub/bTB-bruc-co-infection-ms/run_stochastic_coinfection_model.R', chdir = TRUE)
source('~/GitHub/bTB-bruc-co-infection-ms/rhs_optim.R', chdir = TRUE)   

#############################################################
#############################################################
2) Load Data to fit and data for evaluation
#############################################################
#############################################################

# Evaluation data= Bootstrap prevalence estimates!
#############################################################
data <- read.csv("~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/final_datasets_copied_from_phdfolder/cross_sectional_data_withdz_cleandisease_nofinal_Feb2016_capturetime_forsurv.csv")

overallbruc <- NA; overallbtb <- NA 
brucintbneg <- NA; brucintbpos <- NA
tbinbrucneg <- NA; tbinbrucpos <- NA
overallbrucLS <- NA; overallbtbLS <- NA 
brucintbnegLS <- NA; brucintbposLS <- NA
tbinbrucnegLS <- NA; tbinbrucposLS <- NA
overallbrucCB <- NA; overallbtbCB <- NA 
brucintbnegCB <- NA; brucintbposCB <- NA
tbinbrucnegCB <- NA; tbinbrucposCB <- NA
dataLS <- data[data$herdorig == "LS",]
dataCB <- data[data$herdorig == "CB",]
# REALLY SLOW!
for (i in 1:1000){
	# sample one time point for each individual id (should be 151)
	ssdata <- ddply(data, .(id), function(id) {id[sample(nrow(id), size = 1),]})
	ssdataLS <- ddply(dataLS, .(id), function(id) {id[sample(nrow(id), size = 1),]})
	ssdataCB <- ddply(dataCB, .(id), function(id) {id[sample(nrow(id), size = 1),]})
	
	# calculate prevalence overall and infection specific, all data
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

	# calculate prevalence overall and infection specific, Lower Sabie
	overallbrucLS[i] <- length(ssdataLS$bruc[ssdataLS$bruc=="positive"])/ length(ssdataLS$bruc)
	overallbtbLS[i] <- length(ssdataLS$tb[ssdataLS$tb==1])/ length(ssdataLS$tb)
	brucintbnegLS[i] <- length(ssdataLS$tb[ssdataLS$bruc=="positive" & ssdataLS$tb == 0]) / 
		length(ssdataLS$tb[ssdataLS$tb==0])
	brucintbposLS[i] <- length(ssdataLS$tb[ssdataLS$bruc=="positive" & ssdataLS$tb == 1]) / 
		length(ssdataLS$tb[ssdataLS$tb==1])
	tbinbrucnegLS[i] <- length(ssdataLS$tb[ssdataLS$bruc=="negative" & ssdataLS$tb == 1]) / 
		length(ssdataLS$tb[ssdataLS$bruc=="negative"])
	tbinbrucposLS[i] <- length(ssdataLS$tb[ssdataLS$bruc=="positive" & ssdataLS$tb == 1]) / 
		length(ssdataLS$tb[ssdataLS$bruc=="positive"])

	# calculate prevalence overall and infection specific, Crocodile Bridge
	overallbrucCB[i] <- length(ssdataCB$bruc[ssdataCB$bruc=="positive"])/ length(ssdataCB$bruc)
	overallbtbCB[i] <- length(ssdataCB$tb[ssdataCB$tb==1])/ length(ssdataCB$tb)
	brucintbnegCB[i] <- length(ssdataCB$tb[ssdataCB$bruc=="positive" & ssdataCB$tb == 0]) / 
		length(ssdataCB$tb[ssdataCB$tb==0])
	brucintbposCB[i] <- length(ssdataCB$tb[ssdataCB$bruc=="positive" & ssdataCB$tb == 1]) / 
		length(ssdataCB$tb[ssdataCB$tb==1])
	tbinbrucnegCB[i] <- length(ssdataCB$tb[ssdataCB$bruc=="negative" & ssdataCB$tb == 1]) / 
		length(ssdataCB$tb[ssdataCB$bruc=="negative"])
	tbinbrucposCB[i] <- length(ssdataCB$tb[ssdataCB$bruc=="positive" & ssdataCB$tb == 1]) / 
		length(ssdataCB$tb[ssdataCB$bruc=="positive"])
}

list <- list(overallbrucLS, overallbtbLS, brucintbnegLS, brucintbposLS,
	tbinbrucnegLS, tbinbrucposLS, overallbrucCB, overallbtbCB, 
	brucintbnegCB, brucintbposCB, tbinbrucnegCB, tbinbrucposCB)
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

 lapply(list, quantile, c(0.025, 0.25, 0.5, 0.75, 0.975))
# LOWER SABIE
#     2.5%       25%       50%       75%     97.5% 
#0.3000000 0.3428571 0.3571429 0.3714286 0.4142857  # bruc overall
#0.2000000 0.2428571 0.2714286 0.2857143 0.3285714  # tb overall
#0.2745098 0.3125000 0.3333333 0.3529412 0.3922488  # bruc in TB-
#0.2725758 0.3750000 0.4210526 0.4736842 0.5791796  # bruc in TB + 
#0.1590909 0.2119000 0.2380952 0.2666667 0.3182060 
#0.2000000 0.2758621 0.3200000 0.3571429 0.4347826  
# CROC BRIDGE
#0.2962963 0.3209877 0.3333333 0.3456790 0.3703704  # bruc overall
#0.2345679 0.2592593 0.2839506 0.2962963 0.3209877  # tb overall
#0.2280702 0.2631579 0.2807018 0.2982456 0.3333333 
#0.3913043 0.4400000 0.4761905 0.5000000 0.5652174 
#0.1666667 0.2000000 0.2222222 0.2407407 0.2777778 
#0.3214286 0.3666667 0.3928571 0.4230769 0.4814815 

#############################################################
#############################################################
3) Fit Deterministic Version to overall prevalence
#############################################################
#############################################################

# functions required to interpret results
plot_raw_numbers = function(sol){
	plot(sol$time, sol$S, col= "black", type= 'l', ylim = c(0, 1200), ylab = "Number of animals", xlab = "Time (in years)")
	lines(sol$time, sol$It, col= "red")
	lines(sol$time, sol$Ib, col= "blue")
	lines(sol$time, sol$Ic, col= "green")
	lines(sol$time, sol$R, col = "orange")
	lines(sol$time, sol$Rc, col = "pink")
	legend("topright", legend = c("S", "It", "Ib", "Ic", "R", "Rc"), col = c("black", "red", "blue", "green", "orange", "pink"), bty="n", lty = 1)
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

get_starting_eqbruc = function(params){
	x0 <- c(S = params['K'][[1]]-100, It = 0, Ib = 50, Ic = 0, R = 50, Rc = 0)
	times <- seq(0, 1000, 1)
	sol <- as.data.frame(ode(x0, times, rhs, params))
	out <- c(S = sol$S[length(times)], 
		It = 0, Ib = sol$Ib[length(times)], Ic = 0,
		R = sol$R[length(times)], Rc = 0)
	return(out)
}


# data for fitting, currently using overall, switch to herd specific later
prevTBobs <- 0.27; prevBobs <- 0.34
prevTBobsLS <- 0.27; prevBobsLS <- 0.35
prevTBobsCB <- 0.28; prevBobsCB <- 0.33

# H1:  fit just transmission rates
#############################################################
objective = function(params.est){
	# params.est = 2 long = c(betaB, betaT)
	params <- c(params.fixed, betaB = params.est[1], betaT = params.est[2], rhoB = 2.1, rhoT = 1.3)
	# seed from endemic brucellosis conditions, 10 bTB positive buffalo
	x0 = get_starting_eqbruc(params = c(params))
	x0[[2]] <- x0[[2]] + 10
	x0[[1]] <- x0[[1]] - 10  
	sol <- as.data.frame(ode(x0, times, rhs_optim, params)) # rhs_optim shouldn't get betap's
	df <- groom_sol(sol)
	prevTB <- df$TBprev[length(df[,1])]
	prevB <- df$Brucprev[length(df[,1])]
	error <- sqrt(((prevTBobs - prevTB)^2 + (prevBobs - prevB)^2))
	return (error)
}
betaB = 0.003; betaT = 0.0006
times <- seq(0, 100, 1)
params.fixed = c(fixed.params, gamma=1/2)
# make sure optimizer is smooth
betatest<- seq(0.00001, 0.001, 0.00002); 
df <- data.frame(betaB = rep(betatest, 3), 
	betaT = rep(c(0.0001, 0.005, 0.001), each = length(betatest)), 
	out = NA)
for(i in 1:length(df[,1])){
	df$out[i] <- objective(c(df$betaB[i], df$betaT[i]))
}
par(mfrow = c(1,3))
plot(df$out[df$betaT == 0.0001]~df$betaB[df$betaT == 0.0001])
plot(df$out[df$betaT == 0.005]~df$betaB[df$betaT == 0.005])
plot(df$out[df$betaT == 0.001]~df$betaB[df$betaT == 0.001])

# H1_1: no recovery from mortality rates!,
params.fixed = c(fixed.params, gamma=1/2)
parH1_1 <- optim(c(0.001, 0.001), objective)
paramsH1.1 <- c(params.fixed, betab = parH1_1$par[1], betat = parH1_1$par[2])  # run with rhs_optim

# H1_2: no recovery from mortality rates!, 
params.fixed = c(fixed.params.recov, gamma = 1/2)
parH1_2 <- optim(c(0.001, 0.001), objective) 
paramsH1.2 <- c(params.fixed, betab = parH1_2$par[1], betat = parH1_2$par[2]) 


# H2: fit beta and gamma (DOES NOT WORK WITHONLY PREVALENCE)
#############################################################
objective = function(params.est){
	# params.est = 2 long = c(betaB, betaT, gamma)
	params <- c(params.fixed, betaB = params.est[1], betaT = params.est[2], gamma = params.est[3], rhoB = 4.05, rhoT = 1.3)
	# seed from endemic brucellosis conditions, 10 bTB positive buffalo
	x0 = get_starting_eqbruc(params = c(params))  
	x0[[2]] <- x0[[2]] + 10
	x0[[1]] <- x0[[1]] - 10  
	sol <- as.data.frame(ode(x0, times, rhs_optim, params)) # rhs_optim shouldn't get betap's
	df <- groom_sol(sol)
	prevTB <- df$TBprev[length(df[,1])]
	prevB <- df$Brucprev[length(df[,1])]
	error <- sqrt(((prevTBobs - prevTB)^2 + (prevBobs - prevB)^2))
	return (error)
}
betaB = 0.0003; betaT = 0.0006; gamma = 1/2
# H1: no recovery from mortality rates!,
params.fixed = c(fixed.params)
parH2_1 <- optim(c(0.0003, 0.0001, 0.5), objective)  
parH2_1.1 <- optim(c(0.0001, 0.0001, 0.1), objective) # get different values... based on initials.  Don't estimate both!

#parH2_1LS <- optim(c(0.001, 0.001, 0.5), objective)  
#parH2_1CB <- optim(c(0.001, 0.001, 0.5), objective) 
paramsH2.1 <- c(params.fixed, gamma = parH2_1.1$par[3], betab = parH2_1.1$par[1], betat = parH2_1$par[2])  

# H2: recovery from mortality rates!, 
params.fixed = c(fixed.params.recov)
parH2_2 <- optim(c(0.0001, 0.001, 0.1), objective)
paramsH2.2 <- c(params.fixed, gamma = parH2_2$par[3], betab = parH2_2$par[1], betat = parH2_2$par[2])  


df <- data.frame(paramsestimated = c(rep("transmission", 2), rep("tranmssion&recovery", 2)), 
	recoveryassumption = c("none", "recovery", "none", "recovery"),
	betaB = c(paramsH1.1[[17]], paramsH1.2[[17]], paramsH2.1[[17]], paramsH2.2[[17]]),
	betaT = c(paramsH1.1[[18]], paramsH1.2[[18]], paramsH2.1[[18]], paramsH2.2[[18]]), 
	gamma = c(0.5, 0.5, paramsH2.1[[16]], paramsH2.2[[16]]),
	TBprev = NA, Brucprev = NA, TBprevinS = NA, TBprevinCo = NA, BrucprevinS = NA, BrucprevinCo = NA)
	
get_sum_stats = function(recoveryassumption, params.variable){
	if(recoveryassumption == "none"){
		params = c(fixed.params, params.variable)
	}
	if(recoveryassumption == "recovery"){
		params = c(fixed.params.recov, params.variable)
	}
	times <- seq(0, 1000, 1)
	x0 = get_starting_eqbruc(params = params)
	x0[[2]] <- x0[[2]] + 10
	x0[[1]] <- x0[[1]] - 10  
	sol <- as.data.frame(ode(x0, times, rhs, params))
	df <- groom_sol(sol)
	df$TBprevinS <- df$It / (df$S + df$It)
	df$TBprevinCo <- (df$Ic + df$Rc) / (df$Ib + df$R + df$Ic + df$Rc)
	df$BrucprevinS <- (df$Ib + df$R) / (df$S + df$Ib + df$R)
	df$BrucprevinCo <- (df$Ic + df$Rc) / (df$Ic + df$Rc + df$It)
	
	sumstats <- c(
	df$TBprev[length(df$TBprev)], #TB prev overall
	df$Brucprev[length(df$Brucprev)], #Bruc prev overall
	df$TBprevinS[length(df$TBprevinS)],
	df$TBprevinCo[length(df$TBprevinCo)],
	df$BrucprevinS[length(df$BrucprevinS)],
	df$BrucprevinCo[length(df$BrucprevinCo)]	)

	return(sumstats)
}
	
# fill with predictions!
df$recoveryassumption <- as.character(df$recoveryassumption)
#for (i in 1:length(df[,1])){
for (i in 1:2){
	ra <- df$recoveryassumption[i]
	params.variable <- c(gamma = df$gamma[i], betaB = df$betaB[i], betaT = df$betaT[i],
		rhoB = 3.8, rhoT = 1)
	df[i,c(6:11)] <- get_sum_stats(recoveryassumption = ra, 
		params.variable = params.variable)
	rm(ra, params.variable)
}

# WITH RHOB = 3.8, RHOT = 1
# paramsestimated recoveryassumption        betaB        betaT     gamma    		TBprev  Brucprev TBprevinS 	TBprevinCo BrucprevinS BrucprevinCo
#1        transmission               none 0.0007804677 0.0003206634 0.5000000 0.2699995 0.3399999 0.1777097  0.4491502   0.2565600    0.5655975
#2        transmission           recovery 0.0006768438 0.0001704345 0.5000000 0.2885072 0.3435391 0.1750562  0.5052979   0.2388633    0.6016821
#3 tranmssion&recovery               none 0.0007816310 0.0003206634 0.5010354 0.2699995 0.3399999 0.1777097  0.4491502   0.2565600    0.5655975
#4 tranmssion&recovery           recovery 0.0003010961 0.0002121311 0.1008721 0.2704056 0.3400351 0.1674170  0.4702934   0.2468753    0.5913941

# WITH RhoB = 2.1, RHOT = 1


# WITH RhoB = 4.05, RhoT = 1.3
#     paramsestimated recoveryassumption         betaB         betaT    gamma       TBprev     Brucprev    TBprevinS    TBprevinCo  BrucprevinS BrucprevinCo
#1        transmission               none  0.0007853113  0.0003092429 0.500000 2.497377e-01 3.377416e-01 1.639872e-01  4.178810e-01 2.620494e-01    0.5651361
#2        transmission           recovery  0.0006815121  0.0001591978 0.500000 2.482541e-01 3.347141e-01 1.503945e-01  4.427623e-01 2.481096e-01    0.5969641
#3 tranmssion&recovery               none  0.0003377280  0.0003092429 0.101361 2.413435e-01 3.464417e-01 1.555824e-01  4.031309e-01 2.725612e-01    0.5786829
#4 tranmssion&recovery           recovery -0.0032333333 -0.0056666667 0.110000 1.840844e-72 2.424844e-44 1.722342e-69 -7.095308e-26 2.424844e-44 -934.6264147


# WITH RhoB = 2.1, RhoT = 1.3
#   paramsestimated recoveryassumption         betaB         betaT    gamma    TBprev  Brucprev  TBprevinS TBprevinCo BrucprevinS BrucprevinCo
#1        transmission               none  0.0009316195  0.0002797068 0.500000 0.1381098 0.3895812 0.08004476  0.2290894   0.3484577    0.6462174
#2        transmission           recovery  0.0008050712  0.0001533393 0.500000 0.2175334 0.4089502 0.11139154  0.3709386   0.3287742    0.6973432





# Brucellosis prevalence before bTB (no recovery, RhoB=4, rhoT=1.3 looks best!)- same as after bTB
params = c(params.fixed, betaB= 0.0007853113, betaT = 0.0003092429, rhoT = 1.3, rhoB = 4.05)#check,gamma=0.5
get_starting_eqbruc(params = params)  
#       S        It        Ib        Ic         R        Rc 
#649.42399   0.00000  31.37793   0.00000 313.77927   0.00000 
(31.377+313.779) / (649.424 + 31.377 + 313.779) #0.3470369

# Effect of brucellosis on bTB invasion: 
RoTB = params['betaT']* params['K']/ params['muT']  # 3.061811
RoB = (params[['betaB']] * params[['K']]* (params[['epsilon']] + params[['muR']]) ) /
	( (params[['gamma']] + params[['muB']]) * (params[['epsilon']] + params[['muR']]) + params[['epsilon']] * params[['gamma']])  # 1.106

RoTBwithC = function(params, x){
	# Input = parameters and x = equlibrium brucellosis conditions
	with(as.list(c(x, params)), {
		betapT = rhoT * betaT
		betapB = rhoB * betaB
		Ro = (betapT * R) / (I * betapB + muT) + 
		(betapT * R * betapB * I * (muRC + epsilon)) / ((I * betapB + muT) * (muC * muRC + muC * epsilon + muRC * gamma)) +
		(betapT * R * betapB * I * gamma) / ( (I * betapB + muT) * (muC * muRC + muC * epsilon + muRC * gamma) ) + 
		((betapT * R * (muRC + epsilon)) / (muC * muRC + muC * epsilon + muRC * gamma)) *
		((betapT * R * gamma) / (muC * muRC + muC * epsilon + muRC * gamma)) + 
		((betapT * R * gamma) / (muC * muRC + muC * epsilon + muRC * gamma)) * 
		((betapT * R * epsilon) / (muC * muRC + muC * epsilon + muRC * gamma)) + 
		(betapT * R * (gamma + muC)) / (muC * muRC + muC * epsilon + muRC * gamma)
		return(Ro)
	}
	)
}

get_x_analytic_endemic_brucellosis = function(params, RoB= RoB){ 
	with(as.list(c(params)), {
		S = K / RoB
		
		# solve polynomial to get I 
		delta = gamma / (epsilon + muR)
		ap = (-b2 - b3 * delta - b3 * delta - b3 * (delta^2)) * r / K
		bp = b * b2 + b * b3 * delta - betaB * S - 
			(r / K) * (S + S * delta + b2 *S + b3 * delta * S)
		cp = b * S - muS * S - (r * (S^2) / K)
		I1 = 0.5 * (- bp + sqrt(bp^2 - 4 * ap * cp )) / ap
		I2 = 0.5 * (- bp - sqrt(bp^2 - 4 * ap * cp )) / ap 
		I = max(I1, I2)
		
		# R is a ratio of I
		R = I * delta

		x = c(S= S, I = I, R = R)
	return(x)
	}
	)
}

x = get_x_analytic_endemic_brucellosis(params)
RoTBwithC(params = params, x = x)  # 0.2684736

#############################################################
#############################################################
# 4) Repeat and systematically change gamma and epsilon
#############################################################
#############################################################
df_norecov <- data.frame( 
	recoveryassumption = rep("none", 36),
	betaB = NA,
	betaT = NA, 
	gamma = c(rep(seq(0.1, 0.9, by = 0.1), 4)),
	epsilon = rep(c(0.01, 0.03, 0.06, 0.09), each = 9),
	TBprev = NA, Brucprev = NA, TBprevinS = NA, TBprevinCo = NA, BrucprevinS = NA, BrucprevinCo = NA)
df_recov <- df_norecov; df_recov$recoveryassumption <- "recovery"
df_all_norecov <- df_norecov; df_all_norecov$recoveryassumption <- "none_all"
df_all_recov <- df_norecov; df_all_recov $recoveryassumption <- "recovery_all"

objective_recov = function(params.est, epsilon, gamma){
	params <- c(fixed.params.recov[-15], epsilon = epsilon, gamma = gamma, betaB = params.est[1], 
		betaT = params.est[2])
	# seed from endemic brucellosis conditions, 10 bTB positive buffalo
	x0 = get_starting_eqbruc(params = c(params, betapT = params.est[2], betapB = 2.1 * params.est[1]))  #3.92
	x0[[2]] <- x0[[2]] + 10
	x0[[1]] <- x0[[1]] - 10  
	sol <- as.data.frame(ode(x0, times, rhs_optim, params)) # rhs_optim shouldn't get betap's
	df <- groom_sol(sol)
	prevTB <- df$TBprev[length(df[,1])]
	prevB <- df$Brucprev[length(df[,1])]
	error <- sqrt(((prevTBobs - prevTB)^2 + (prevBobs - prevB)^2))
	return (error)
}

objective_norecov = function(params.est, epsilon, gamma){
	params <- c(fixed.params[-15], epsilon = epsilon, gamma = gamma, betaB = params.est[1], 
		betaT = params.est[2])
	# seed from endemic brucellosis conditions, 10 bTB positive buffalo
	x0 = get_starting_eqbruc(params = c(params, betapT = params.est[2], betapB = 2.1 * params.est[1]))   # 3.92
	x0[[2]] <- x0[[2]] + 10
	x0[[1]] <- x0[[1]] - 10  
	sol <- as.data.frame(ode(x0, times, rhs_optim, params)) # rhs_optim shouldn't get betap's
	df <- groom_sol(sol)
	prevTB <- df$TBprev[length(df[,1])]
	prevB <- df$Brucprev[length(df[,1])]
	error <- sqrt(((prevTBobs - prevTB)^2 + (prevBobs - prevB)^2))
	return (error)
}

prevTBinSobs <- 0.24  # these are Lower Sabie values
prevTBinCoobs <-0.32
prevBinSobs <- 0.33
prevBinCoobs <- 0.42
objective_all_norecov = function(params.est, epsilon, gamma){
	params <- c(fixed.params[-15], epsilon = epsilon, gamma = gamma, betaB = params.est[1], 
		betaT = params.est[2])
	# seed from endemic brucellosis conditions, 10 bTB positive buffalo
	x0 = get_starting_eqbruc(params = c(params, betapT = params.est[2], betapB = 2.1 * params.est[1]))   # 3.92
	x0[[2]] <- x0[[2]] + 10
	x0[[1]] <- x0[[1]] - 10  
	sol <- as.data.frame(ode(x0, times, rhs_optim, params)) # rhs_optim shouldn't get betap's
	df <- groom_sol(sol)
	df$TBprevinS <- df$It / (df$S + df$It)
	df$TBprevinCo <- (df$Ic + df$Rc) / (df$Ib + df$R + df$Ic + df$Rc)
	df$BrucprevinS <- (df$Ib + df$R) / (df$S + df$Ib + df$R)
	df$BrucprevinCo <- (df$Ic + df$Rc) / (df$Ic + df$Rc + df$It)

	prevTBinS <- df$TBprevinS[length(df[,1])]
	prevTBinCo <- df$TBprevinCo[length(df[,1])]
	prevBinS <- df$BrucprevinS[length(df[,1])]
	prevBinCo <- df$BrucprevinCo[length(df[,1])]
	
	error <- sqrt(((prevTBinSobs - prevTBinS)^2 + (prevTBinCoobs - prevTBinCo)^2 +
		(prevBinSobs - prevBinS)^2 + (prevBinCoobs - prevBinCo)^2))
	return (error)
}

objective_all_recov = function(params.est, epsilon, gamma){
	params <- c(fixed.params.recov[-15], epsilon = epsilon, gamma = gamma, betaB = params.est[1], 
		betaT = params.est[2])
	# seed from endemic brucellosis conditions, 10 bTB positive buffalo
	x0 = get_starting_eqbruc(params = c(params, betapT = params.est[2], betapB = 2.1 * params.est[1]))  # 3.92
	x0[[2]] <- x0[[2]] + 10
	x0[[1]] <- x0[[1]] - 10  
	sol <- as.data.frame(ode(x0, times, rhs_optim, params)) # rhs_optim shouldn't get betap's
	df <- groom_sol(sol)
	df$TBprevinS <- df$It / (df$S + df$It)
	df$TBprevinCo <- (df$Ic + df$Rc) / (df$Ib + df$R + df$Ic + df$Rc)
	df$BrucprevinS <- (df$Ib + df$R) / (df$S + df$Ib + df$R)
	df$BrucprevinCo <- (df$Ic + df$Rc) / (df$Ic + df$Rc + df$It)

	prevTBinS <- df$TBprevinS[length(df[,1])]
	prevTBinCo <- df$TBprevinCo[length(df[,1])]
	prevBinS <- df$BrucprevinS[length(df[,1])]
	prevBinCo <- df$BrucprevinCo[length(df[,1])]
	
	error <- sqrt(((prevTBinSobs - prevTBinS)^2 + (prevTBinCoobs - prevTBinCo)^2 +
		(prevBinSobs - prevBinS)^2 + (prevBinCoobs - prevBinCo)^2))
	return (error)
}

get_params = function(recoveryassumption, epsilon, gamma){
	# Input: recovery assumption (none or recovery)
	# epsilon and gamma values. 
	# Output: vector containing c(beta_b, beta_t) 
	if (recoveryassumption == "none"){
		par <- optim(c(0.003, 0.003), objective_norecov, epsilon = epsilon, gamma = gamma) 
	}
	if (recoveryassumption == "recovery"){
		par <- optim(c(0.003, 0.003), objective_recov, epsilon = epsilon, gamma = gamma)
	}
	if (recoveryassumption == "none_all"){
		par <- optim(c(0.003, 0.003), objective_all_norecov, epsilon = epsilon, gamma = gamma) 
	}
	if (recoveryassumption == "recovery_all"){
		par <- optim(c(0.003, 0.003), objective_all_recov, epsilon = epsilon, gamma = gamma) 
	}
	return(c(par$par[1], par$par[2]))  # betab, betat
}
	
get_sum_stats_ge = function(params){
	times <- seq(0, 1000, 1)
	x0 = get_starting_eqbruc(params = params)
	x0[[2]] <- x0[[2]] + 10
	x0[[1]] <- x0[[1]] - 10  
	sol <- as.data.frame(ode(x0, times, rhs, params))
	df <- groom_sol(sol)
	df$TBprevinS <- df$It / (df$S + df$It)
	df$TBprevinCo <- (df$Ic + df$Rc) / (df$Ib + df$R + df$Ic + df$Rc)
	df$BrucprevinS <- (df$Ib + df$R) / (df$S + df$Ib + df$R)
	df$BrucprevinCo <- (df$Ic + df$Rc) / (df$Ic + df$Rc + df$It)
	
	sumstats <- c(
	df$TBprev[length(df$TBprev)], #TB prev overall
	df$Brucprev[length(df$Brucprev)], #Bruc prev overall
	df$TBprevinS[length(df$TBprevinS)],
	df$TBprevinCo[length(df$TBprevinCo)],
	df$BrucprevinS[length(df$BrucprevinS)],
	df$BrucprevinCo[length(df$BrucprevinCo)]	)

	return(sumstats)
}	
	
for (i in 1:length(df_norecov[,1])){
	params <- get_params("none", epsilon = df_norecov$epsilon[i], gamma = df_norecov$gamma[i])
	df_norecov$betaB[i] <- params[1]
	df_norecov$betaT[i] <- params[2]
	# get summary statistics
	temp_params <- c(params.fixed[-15], epsilon = df_norecov$epsilon[i], gamma = df_norecov$gamma[i],
		 betaB = df_norecov$betaB[i], betaT = df_norecov$betaT[i],
		betapT = df_norecov$betaT[i], betapB = 3.92 * df_norecov$betaB[i])
	df_norecov[i,c(6:11)] <- get_sum_stats_ge(params = temp_params)
	rm(temp_params)
}
	
for (i in 1:length(df_recov[,1])){
	params <- get_params("recovery", epsilon = df_recov$epsilon[i], gamma = df_recov$gamma[i])
	df_recov$betaB[i] <- params[1]
	df_recov$betaT[i] <- params[2]
	# get summary statistics
	temp_params <- c(fixed.params.recov[-15], epsilon = df_recov$epsilon[i], gamma = df_recov$gamma[i],
		betaB = df_recov$betaB[i], betaT = df_recov$betaT[i],
		betapT = df_recov$betaT[i], betapB = 3.92 * df_recov$betaB[i])
	print(temp_params)
	df_recov[i,c(6:11)] <- get_sum_stats_ge(params = temp_params)
	rm(temp_params)
}
write.csv(df_recov, "~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/modelpredictions_varygammaepsilon_recovery.csv")
write.csv(df_norecov, "~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/modelpredictions_varygammaepsilon_norecovery.csv")
# LESSON: When only fitting to overall prevalence, we note that epsilon, gamma, and the betas just move to get same values. 

for (i in 1:length(df_all_norecov[,1])){
	# get parameters betaB, betaT
	params <- get_params("none_all", epsilon = df_all_norecov$epsilon[i], gamma = df_all_norecov$gamma[i])
	df_all_norecov$betaB[i] <- params[1]
	df_all_norecov$betaT[i] <- params[2]
	# get summary statistics
	temp_params <- c(params.fixed[-15], epsilon = df_all_norecov$epsilon[i], gamma = df_all_norecov$gamma[i],
		betaB = df_all_norecov$betaB[i], betaT = df_all_norecov$betaT[i],
		betapT = df_all_norecov$betaT[i], betapB = 3.92 * df_all_norecov$betaB[i])
	df_all_norecov[i,c(6:11)] <- get_sum_stats_ge(params = temp_params)
	rm(temp_params)
	print(i)
}
for (i in 1:length(df_all_recov[,1])){
	# get parameters betaB, betaT
	params <- get_params("recovery_all", epsilon = df_all_recov$epsilon[i], gamma = df_all_recov$gamma[i])
	df_all_recov$betaB[i] <- params[1]
	df_all_recov$betaT[i] <- params[2]
	# get summary statistics
	temp_params <- c(params.fixed[-15], epsilon = df_all_recov$epsilon[i], gamma = df_all_recov$gamma[i],
		betaB = df_all_recov$betaB[i], betaT = df_all_recov$betaT[i],
		betapT = df_all_recov$betaT[i], betapB = 3.92 * df_all_recov$betaB[i])
	df_all_recov[i,c(6:11)] <- get_sum_stats_ge(params = temp_params)
	rm(temp_params)
	print(i)
}
write.csv(df_all)norecov, "~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/modelpredictions_varygammaepsilon_all_norecovery.csv")

df_recov$epsilon = as.character(df_recov$epsilon)
source('~/GitHub/bTB-bruc-co-infection-ms/multiplot.R', chdir = TRUE)

p1 <- ggplot(data= df_all_recov, aes(x = gamma, y = TBprevinS, group = epsilon)) + 
	geom_line(aes(colour = epsilon)) + geom_point() + 
	geom_hline(yintercept = 0.23) + geom_hline(yintercept = 0.22)
p2 <- ggplot(data= df_all_recov, aes(x = gamma, y = TBprevinCo, color = epsilon)) + 
	geom_line(aes(colour = epsilon)) + geom_point() +
	geom_hline(yintercept = 0.32) + geom_hline(yintercept = 0.39)
p3 <- ggplot(data= df_all_recov, aes(x = gamma, y = BrucprevinS, group = epsilon)) + 
	geom_line(aes(colour=epsilon)) + geom_point() + 
	geom_hline(yintercept = 0.28) + geom_hline(yintercept = 0.33)
p4 <- ggplot(data= df_all_recov, aes(x = gamma, y = BrucprevinCo, group = epsilon)) + 
	geom_point() + geom_line(aes(colour=epsilon)) + geom_hline(yintercept = 0.42)+ geom_hline(yintercept = 0.47) 
multiplot(p1, p2, p3, p4, cols = 2)

p1 <- ggplot(data= df_all_norecov, aes(x = gamma, y = TBprevinS, group = epsilon)) + 
	geom_line(aes(colour = epsilon)) + geom_point() + 
	geom_hline(yintercept = 0.23) + geom_hline(yintercept = 0.22)
p2 <- ggplot(data= df_all_norecov, aes(x = gamma, y = TBprevinCo, color = epsilon)) + 
	geom_line(aes(colour = epsilon)) + geom_point() +
	geom_hline(yintercept = 0.32) + geom_hline(yintercept = 0.39)
p3 <- ggplot(data= df_all_norecov, aes(x = gamma, y = BrucprevinS, group = epsilon)) + 
	geom_line(aes(colour=epsilon)) + geom_point() + 
	geom_hline(yintercept = 0.28) + geom_hline(yintercept = 0.33)
p4 <- ggplot(data= df_all_norecov, aes(x = gamma, y = BrucprevinCo, group = epsilon)) + 
	geom_point() + geom_line(aes(colour=epsilon)) + geom_hline(yintercept = 0.42)+ geom_hline(yintercept = 0.47) 
multiplot(p1, p2, p3, p4, cols = 2)
#############################################################
#############################################################
5) Plot deterministic version
#############################################################
#############################################################

















#############################################################
#############################################################
6) Impliment Stochastic Version
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




# MORE CUT STUFF!
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

##################################################################################################
##################################################################################################
#1) Load fixed parameters, model
##################################################################################################
##################################################################################################
# Accessory functions for fitting: 
rm(list = ls())
require("deSolve")
library("gridExtra")
library("ggplot2")
# get fixed.params & fixed.params.recov
source('~/GitHub/bTB-bruc-co-infection-ms/fixed_parameters_norecovery_agematrix.R', chdir = TRUE)
source('~/GitHub/bTB-bruc-co-infection-ms/fixed_parameters_recovery_agematrix.R', chdir = TRUE)
# rhs function, determinitic model, age structure
source('~/GitHub/bTB-bruc-co-infection-ms/rhs_age.R', chdir = TRUE)

# age divisions in rhs function (age= 1-3.9, 4-4.9, 5-14.9, 15+)
s_index <- 1:20
it_index <- 21:40
ib_index <- 41:60
ic_index <- 61:80
r_index <- 81:100
rc_index <- 101:120
juveniles <- 1:3
subadult<- 4
adult <- 5:14
mature <- 15:20

# Stating age structure information 
relageall = c(0.137, rep(0.368/4, 4), rep(0.185/4, 4),  
	rep(0.235/6, 6), rep(0.075/5, 5))					
relage = relageall
S0 = 400*relage; It0 = 0*relage; Ib0 = 0*relage; 
Ic0 = 0*relage; R0 = 0 * relage; Rc0 = 0 * relage
x0 = c(S0, It0, Ib0, Ic0, R0, Rc0)
times <- seq(0, 1000, 1)
params.test_log <- c(fixed.params, list(gamma=1/2, theta = 4, K = 433,
	betaB = 1.025, betaT = 0.000223, rhoT = 1, rhoB = 2.1))
sol <- as.data.frame(ode(x0, times, rhs_age_matrix, params.test_log))
stable_age <- unname(unlist( sol[1000, c(2:21)]/sum(sol[500, c(2:21)]) ))


plot_agestructure = function(x){ # true based on Jolles et al. 2007
	 xcounts <- NA
	 if (length(x) != 120){
	 	print("The age structure should include 20 ages,
	 		for 6 disease classes, giving 120 columns")
	 }
	 for(i in 1:20){
	 	xcounts[i] <- (x[i] + x[20+i] + x[40+i] + x[60+i] + x[80+i] + x[100+i])/sum(x)
	 }
	 d<-matrix(c(relageall, 
	 	xcounts), nrow=2, byrow=TRUE, 
	 	dimnames=list(c("Observed", "Predicted"), 
	 	c(seq(1:20))))
	barplot(d, beside = TRUE, col = c("light gray", "dark gray"), ylab = "Frequency (%)")
	box(	)
	legend("topleft", legend = c("Observed", "Predicted"), fill = c("light gray", "dark gray"))
}

plot_ageprevalence = function(sol){
	S <- rep(0, 20); It <- rep(0, 20); Ib <- rep(0, 20);
	Ic <- rep(0, 20); R <- rep(0, 20); Rc <- rep(0, 20)
	N <- rep(0, 20)
	for (i in 1:20){
		S[i] <- sum(sol[length(sol) , i + 1])
		It[i] <- sum(sol[length(sol) , i + 21])
		Ib[i] <- sum(sol[length(sol) , 41 +i])
		Ic[i] <- sum(sol[length(sol) , 61 +i])
		R[i] <- sum(sol[length(sol) , 81 +i])
		Rc[i] <- sum(sol[length(sol) , 101 +i])
		N[i] <- sum(sol[length(sol), c(1+i, 21+i, 41+i, 61 + i, 81+i, 101+i)])
	}
	
	prevB <- (Ib + Ic + R + Rc) / N
	prevT <- (It + Ic + Rc)/ N
	prevBinS <- (Ib + R) / (S + Ib + R)
	prevBinT <- (Ic + Rc) / (It + Ic + Rc) 
	
	overall_prevB <- (sum(Ib) + sum(Ic) + sum(R) + sum(Rc)) / sum(N)
	overall_prevT <- (sum(It) + sum(Ic) + sum(Rc)) / sum(N)
	overall_prevBinS <- (sum(Ib) + sum(R)) / (sum(S) + sum(Ib) + sum(R))
	overall_prevBinT <- (sum(Ic) + sum(Rc)) / (sum(It) + sum(Ic) + sum(Rc))
	overallN = sum(N) 
	
	#par(mfrow = c(1,2))
	#plot(y = prevB, x= seq(1, 20, 1), type = "b", col = "dark blue", ylim = c(0, 0.8),
	#	ylab = "Prevalence", xlab = "Age", pch = 19, 
	#	main = paste("Overall prevalences, Br =", round(overall_prevB, 3), 
	#	" TB = ", round(overall_prevT, 3) ))
	#points(y = prevT, x = seq(1,20,1), type = "b", col = "dark red", pch = 19)
	#legend("bottomright", bty = "n", legend = c("Bruc", "TB"), 
	#	pch = c(19, 19), col = c("dark blue", "dark red"))
	plot(y = prevBinS, x = seq(1,20,1), type = "b", col = "dark blue", pch = 19, 
		ylab = "Brucellosis prevalence", xlab = "Age", ylim = c(0, 0.8), 
		main = paste("Br|S =", round(overall_prevBinS, 3), 
		" Br|Co = ", round(overall_prevBinT, 3) )
		)
	text(x = 10, y = 0.55, labels = paste("Final N = ", round(overallN, 2)))
	points(y = prevBinT, x = seq(1,20,1), type = "b", col = "dark red", pch = 19)
}

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

# use params or params_recov based on assumptions
get_starting_eqbruc = function(params){
	S0 = 400* stable_age; It0 = 0 * stable_age; Ib0 = 20* stable_age; 
	Ic0 = 0* stable_age; R0 = 30 * stable_age; Rc0 = 0 * stable_age
	x0 = c(S0, It0, Ib0, Ic0, R0, Rc0)
	times <- seq(0, 1000, 1)
	sol <- as.data.frame(ode(x0, times, rhs_age_matrix, params))
	out <- unname(unlist( sol[1000, c(2:121)] ))
	return(out)
}
test <- get_starting_eqbruc(params.test_log); sum(test[21:120])/sum(test)

get_prevalence = function(sol){
	S <- sum(sol[length(sol[,1]) , s_index+1])
	It <- sum(sol[length(sol[,1]) , it_index +1])
	Ib <- sum(sol[length(sol[,1]) , ib_index +1])
	Ic <- sum(sol[length(sol[,1]) , ic_index +1])
	R <- sum(sol[length(sol[,1]) , r_index +1])
	Rc <- sum(sol[length(sol[,1]) , rc_index +1])
	N <- sum(sol[length(sol[,1]), 2:121])
	prevTB <- (It + Ic + Rc) / N 
	prevB <- (Ib + Ic + R + Rc) / N
	prevBinS <- (Ib + R) / (S + Ib + R)
	prevBinT <- (Ic + Rc) / (It + Ic + Rc)
	prevTinS <- (It) / (S + It)
	prevTinB <- (Ic + Rc) / (Ib + Ic + R + Rc)
	return(list(prevTB = prevTB, prevB = prevB,
		prevBinS = prevBinS, prevBinT = prevBinT, 
		prevTinS = prevTinS, prevTinB = prevTinB))
}

make_summary_plots = function(sol){
	df <- get_prevalence(sol)
	df2 <- data.frame(Evaluation = c("Model", "Model", "Data", "Data"), 
		Infection = c("Single", "Co", "Single", "Co"), 
		BrucellosisPrevalence = c(df$prevBinS, df$prevBinT, 0.3035, 0.4524), 
		TBPrevalence = c(df$prevTinS, df$prevTinB, 0.227, 0.3585))
	p1 <- ggplot(df2, aes(x = Evaluation, y = BrucellosisPrevalence, fill = Infection)) + 
		geom_bar(position = position_dodge(), stat = "identity") + ylim(0, 0.8)
	p2 <- ggplot(df2, aes(x = Evaluation, y = TBPrevalence, fill = Infection)) + 
		geom_bar(position = position_dodge(), stat = "identity") + ylim(0, 0.8)
	grid.arrange(p1, p2, ncol = 2)
}

make_structured_summary_plots = function(sol){
	df <- get_structured_prevalence(sol)
	df2 <- data.frame(Evaluation = c("Model", "Model", "Data", "Data"), 
		Infection = c("Single", "Co", "Single", "Co"), 
		BrucellosisPrevalence = c(df$prevBinS, df$prevBinT, 0.3035, 0.4524), 
		TBPrevalence = c(df$prevTinS, df$prevTinB, 0.227, 0.3585))
	p1 <- ggplot(df2, aes(x = Evaluation, y = BrucellosisPrevalence, fill = Infection)) + 
		geom_bar(position = position_dodge(), stat = "identity") + ylim(0, 0.8)
	p2 <- ggplot(df2, aes(x = Evaluation, y = TBPrevalence, fill = Infection)) + 
		geom_bar(position = position_dodge(), stat = "identity") + ylim(0, 0.8)
	grid.arrange(p1, p2, ncol = 2)
}


	
##################################################################################################
##################################################################################################
# Test 1: Fit model to overall bTB and bruc prevalences (estimate betaB, betaT)
# do not subsample model output to match age structure
##################################################################################################
##################################################################################################	
prevTBobs <- 0.27  # for test- bootstrap estimate of overall prevalence
prevBobs <- 0.34
objective = function(params.est){
	# params.est = 2 long = c(betaB, betaT)
	params <- c(fixed.params, list(gamma=1/2, betaB = params.est[1],
	betaT = params.est[2]/10000, rhoT = 1, rhoB = 2.1, theta= 4, K = 433))
	
	# seed from endemic brucellosis conditions, 10 bTB positive buffalo
	x0 = get_starting_eqbruc(params = c(params))
	x0[[25]] <- x0[[25]] + 2
	if(x0[[5]] > 2){
		x0[[5]] <- x0[[5]] - 2}  
	sol <- as.data.frame(ode(x0, times, rhs_age_matrix, params))
	df <- get_prevalence(sol)
	error <- sqrt(((prevTBobs - df$prevTB)^2 + (prevBobs - df$prevB)^2))
	return (error)
}
par <- optim(c(1.025, 0.00054*10000), objective)  # starting vals=30% prev w/out co-infection
# 1.082371 4.546418 # Value: 0.27 at rho = 4
#params <- c(fixed.params, list(gamma=1/2, betaB = 1.082371,
#	betaT = 4.546418/10000, rhoT = 1, rhoB = 2.1, theta= 4, K = 433))

# 0.9353968, 8.6858712 # at rho = 2.1
params <- c(fixed.params, list(gamma=1/2, betaB = 0.9353968,
	betaT = 8.6858712/10000, rhoT = 1, rhoB = 2.1, theta= 4, K = 433))
x0 = get_starting_eqbruc(params = c(params))
x0[[25]] <- x0[[25]] + 2
x0[[5]] <- x0[[5]] - 2
sol <- as.data.frame(ode(x0, times, rhs_age_matrix, params))
get_prevalence(sol)  # 27, 34!!!  (26.6, 53.8, 28.8, 42.8)
make_summary_plots(sol)
par(mfrow = c(1, 2))
plot_raw_numbers(sol)
plot_ageprevalence(sol)

objective_recov = function(params.est){
	# params.est = 2 long = c(betaB, betaT)
	params <- c(fixed.params.recov, list(gamma=1/2, betaB = params.est[1],
	betaT = params.est[2]/10000, rhoT = 1, rhoB = 2.1, theta= 4, K = 433))
	
	# seed from endemic brucellosis conditions, 10 bTB positive buffalo
	x0 = get_starting_eqbruc(params = c(params))
	x0[[25]] <- x0[[25]] + 2
	if(x0[[5]] > 2){
		x0[[5]] <- x0[[5]] - 2}  
	sol <- as.data.frame(ode(x0, times, rhs_age_matrix, params))
	df <- get_prevalence(sol)
	error <- sqrt(((prevTBobs - df$prevTB)^2 + (prevBobs - df$prevB)^2))
	return (error)
}

par_recov <- optim(c(0.8, 0.00044*10000), objective_recov)  # starting vals=30% prev w/out co-infection
# with rho = 4
# 0.7812941 5.6520986# Value: 2.68*10^-9  
#params <- c(fixed.params.recov, list(gamma=1/2, betaB = 0.7812941,
#	betaT = 5.6520986/10000, rhoT = 1, rhoB = 4, theta= 4, K = 433))

#with rho = 2: 0.9327971 5.5166360
params <- c(fixed.params.recov, list(gamma=1/2, betaB = 0.9327971,
	betaT = 5.5166360/10000, rhoT = 1, rhoB = 4, theta= 4, K = 433))
x0 = get_starting_eqbruc(params = c(params))
x0[[25]] <- x0[[25]] + 2
x0[[5]] <- x0[[5]] - 2
sol <- as.data.frame(ode(x0, times, rhs_age_matrix, params))
get_prevalence(sol)
make_summary_plots(sol)
par(mfrow = c(1, 2))
plot_raw_numbers(sol)
plot_ageprevalence(sol)
#plot_agestructure(sol[length(sol[,1]), c(2:121)])


# ASIDE:
brucinS <- 0.3035 (0.265-0.345)
brucinT <- 0.4524 (0.365- 0.5406)
tbinS <- 0.227 (0.1782-0.2784)
tbinB <- 0.3585 (0.1857 - 0.4314)

#prevTBinSobs <- 0.24  # these are Lower Sabie values
#prevTBinCoobs <-0.32
#prevBinSobs <- 0.33
#prevBinCoobs <- 0.42

# to plot objective function over rage of betaB, betaT (for no recov, rho=4)
bB <- seq(0.05, 1.05, 0.01)
bT <- seq(1, 10, 1)

bB <- seq(0.7, 1.05, 0.01)
bT <- seq(2, 10, 1)

val <- NA
for(b in bB){
	for(t in bT){
		i <- 1
		val[i] <- objective(c(b, t))
		i <- i + 1
	}
}
val
df <- matrix(val, nrow = length(bT) byrow = TRUE)
##################################################################################################
##################################################################################################
# Test 2: Fit model to overall bTB and bruc prevalences (estimate betaB, betaT)
# subsample model output to match age structure
##################################################################################################
##################################################################################################
data <- read.csv("~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/final_datasets_copied_from_phdfolder/cross_sectional_data_withdz_cleandisease_nofinal_Feb2016_capturetime_forsurv.csv")
counts<- hist(data$age_sel/12, plot = FALSE)$counts  # youngest = 1.4 so aged 1-2
agestructure<- counts/sum(counts)
data_agestructure = c(agestructure, 0, 0, 0, 0, 0)

get_structured_prevalence = function(sol){
	S <-sum(sol[length(sol[,1]) , s_index+1] * data_agestructure)  # should give a scalar
	It <- sum(sol[length(sol[,1]) , it_index +1] * data_agestructure)
	Ib <- sum(sol[length(sol[,1]) , ib_index +1] * data_agestructure)
	Ic <- sum(sol[length(sol[,1]) , ic_index +1] * data_agestructure)
	R <- sum(sol[length(sol[,1]) , r_index +1] * data_agestructure)
	Rc <-sum(sol[length(sol[,1]) , rc_index +1] * data_agestructure)
	N <- sum(S + It + Ib + Ic + R + Rc)
	prevTB <- (It + Ic + Rc) / N 
	prevB <- (Ib + Ic + R + Rc) / N
	prevBinS <- (Ib + R) / (S + Ib + R)
	prevBinT <- (Ic + Rc) / (It + Ic + Rc)
	prevTinS <- (It) / (S + It)
	prevTinB <- (Ic + Rc) / (Ib + Ic + R + Rc)
	return(list(prevTB = prevTB, prevB = prevB,
		prevBinS = prevBinS, prevBinT = prevBinT, 
		prevTinS = prevTinS, prevTinB = prevTinB))
}


objective = function(params.est){
	# params.est = 2 long = c(betaB, betaT)
	#params <- c(fixed.params, list(gamma=1/2, betaB = params.est[1],
	#betaT = params.est[2]/10000, rhoT = 1, rhoB = 2.1, theta= 4, K = 433))
	params <- c(fixed.params, list(gamma=1/2, betaB = params.est[1],
	betaT = params.est[2]/10000, rhoT = 1, rhoB = 1.2, theta= 4, K = 433))
	
	# seed from endemic brucellosis conditions, 10 bTB positive buffalo
	x0 = get_starting_eqbruc(params = c(params))
	x0[[25]] <- x0[[25]] + 2
	if(x0[[5]] > 2){
		x0[[5]] <- x0[[5]] - 2}  
	sol <- as.data.frame(ode(x0, times, rhs_age_matrix, params))
	df <- get_structured_prevalence(sol)
	error <- sqrt(((prevTBobs - df$prevTB)^2 + (prevBobs - df$prevB)^2))
	return (error)
}

objective_recov = function(params.est){
	# params.est = 2 long = c(betaB, betaT)
	params <- c(fixed.params.recov, list(gamma=1/2, betaB = params.est[1],
	betaT = params.est[2]/10000, rhoT = 1, rhoB = 2.1, theta= 4, K = 433))
	
	# seed from endemic brucellosis conditions, 10 bTB positive buffalo
	x0 = get_starting_eqbruc(params = c(params))
	x0[[25]] <- x0[[25]] + 2
	if(x0[[5]] > 2){
		x0[[5]] <- x0[[5]] - 2}  
	sol <- as.data.frame(ode(x0, times, rhs_age_matrix, params))
	df <- get_structured_prevalence(sol)
	error <- sqrt(((prevTBobs - df$prevTB)^2 + (prevBobs - df$prevB)^2))
	return (error)
}

# Rho = 4, no recovery
par <- optim(c(1.025, 0.00094*10000), objective) 
#0.797551 9.640622
params <- c(fixed.params, list(gamma=1/2, betaB = 0.797551,
	betaT = 9.640622/10000, rhoT = 1, rhoB = 4, theta= 4, K = 433))
x0 = get_starting_eqbruc(params = c(params))
x0[[25]] <- x0[[25]] + 2
x0[[5]] <- x0[[5]] - 2
sol <- as.data.frame(ode(x0, times, rhs_age_matrix, params))
get_prevalence(sol)  # 27, 34!!!  (26.6, 53.8, 28.8, 42.8)
make_structured_summary_plots(sol)
par(mfrow = c(1, 2))
plot_raw_numbers(sol)
plot_ageprevalence(sol)

# Rho = 2.1, no recovery
par <- optim(c(1.025, 0.00094*10000), objective)  # starting vals=30% prev w/out co-infection
#0.9636965 9.1928187
params <- c(fixed.params, list(gamma=1/2, betaB = 0.9636965,
	betaT = 9.1928187/10000, rhoT = 1, rhoB = 2.1, theta= 4, K = 433))
x0 = get_starting_eqbruc(params = c(params))
x0[[25]] <- x0[[25]] + 2
x0[[5]] <- x0[[5]] - 2
sol <- as.data.frame(ode(x0, times, rhs_age_matrix, params))
get_prevalence(sol)  # 27, 34!!!  (26.6, 53.8, 28.8, 42.8)
get_structured_prevalence(sol)
make_structured_summary_plots(sol)
par(mfrow = c(1, 2))
plot_raw_numbers(sol)
plot_ageprevalence(sol)

# rho = 1.5, no recovery
par <- optim(c(1.025, 0.00094*10000), objective)  # starting vals=30% prev w/out co-infection
params <- c(fixed.params, list(gamma=1/2, betaB = 1.038374,
	betaT = 9.001541/10000, rhoT = 1, rhoB = 1.2, theta= 4, K = 433))

# Rho = 4, recovery

# Rho = 2.1, recovery





##################################################################################################
##################################################################################################
# Test 3: Fit model to minimize co-infection patterns at of the first capture (estimate betaB, betaT)
##################################################################################################
##################################################################################################
# in first capture
data <- read.csv("~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/final_datasets_copied_from_phdfolder/cross_sectional_data_withdz_cleandisease_nofinal_Feb2016_capturetime_forsurv.csv")
c1<- data[data$capturetime < 6,]
counts<- hist(c1$age_sel/12, plot = FALSE)$counts
agestructure<- counts/sum(counts)

prevBinSobs <- NA; prevBinTobs <- NA; 
prevTinSobs <- NA; prevTinBobs <- NA
objective = function(params.est){
	# params.est = 2 long = c(betaB, betaT)
	params <- c(fixed.params, list(gamma=1/2, betaB = params.est[1],
	betaT = params.est[2], rhoT = 1, rhoB = 4, theta= 4, K = 433))
	
	# seed from endemic brucellosis conditions, 10 bTB positive buffalo
	x0 = get_starting_eqbruc(params = c(params))
	x0[[25]] <- x0[[25]] + 2
	if(x0[[5]] > 2){
		x0[[5]] <- x0[[5]] - 2}  
	sol <- as.data.frame(ode(x0, times, rhs_age_matrix, params))
	df <- get_prevalence(sol)
	error <- sqrt(((prevBinSobs - df$prevBinS)^2 + (prevBinTobs - df$prevBinT)^2 +
		(prevTinSobs -df$prevTinS)^2 + (prevTinBobs- df$prevTinB)^2 ))
	return (error)
}

objective_recov = function(params.est){
	# params.est = 2 long = c(betaB, betaT)
	params <- c(fixed.params.recov, list(gamma=1/2, betaB = params.est[1],
	betaT = params.est[2], rhoT = 1, rhoB = 4, theta= 4, K = 433))
	
	# seed from endemic brucellosis conditions, 10 bTB positive buffalo
	x0 = get_starting_eqbruc(params = c(params))
	x0[[25]] <- x0[[25]] + 2
	if(x0[[5]] > 2){
		x0[[5]] <- x0[[5]] - 2}  
	sol <- as.data.frame(ode(x0, times, rhs_age_matrix, params))
	df <- get_prevalence(sol)
	error <- sqrt(((prevBinSobs - df$prevBinS)^2 + (prevBinTobs - df$prevBinT)^2 +
		(prevTinSobs -df$prevTinS)^2 + (prevTinBobs- df$prevTinB)^2 ))
	return (error)
}

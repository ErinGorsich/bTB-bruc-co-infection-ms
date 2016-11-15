# Age structured model owes much to:  http://ms.mcmaster.ca/~bolker/eeid/2011_eco/waifw.pdf

#############################################################
#############################################################
#1) Load fixed parameters, model
#############################################################
#############################################################
rm(list = ls())
require("deSolve")
library("plyr")
library("ggplot2")
set.seed(5)
# get fixed.params & fixed.params.recov
source('~/GitHub/bTB-bruc-co-infection-ms/fixed_parameters_norecovery_agematrix.R', chdir = TRUE)
source('~/GitHub/bTB-bruc-co-infection-ms/fixed_parameters_recovery_agematrix.R', chdir = TRUE)
# rhs function, determinitic model, age structure
source('~/GitHub/bTB-bruc-co-infection-ms/rhs_age.R', chdir = TRUE)


#############################################################
#############################################################
#2) Set-up features of aging
#############################################################
#############################################################
# age divisions in rhs function
#(age= 1-3.9, 4-4.9, 5-14.9, 15+)..subsume calf mortality in births
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

# Age structure information, used to calculate mortality rates in susceptibles. 
relageall = c(0.137, rep(0.368/4, 4), rep(0.185/4, 4),  # Jolles 2007, set max age at 20
	rep(0.235/6, 6), rep(0.075/5, 5))					# Also in Caron et al. from 2001 KNP
	
relage = relageall

plot_agestructure = function(x){ # true based on Jolles et al. 2007
	 xcounts <- NA
	 if (length(x) != 120){
	 	print("The age structure should include 20 ages, for 6 disease classes, giving 120 columns")
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
plot_agestructure(x = seq(1, 20))

plot_dz_agestructure = function(x, dz){
	# Disease options =
	# "bruc" = brucellosis only
	# "tb" = bTB only
	# "co" = all subtypes
	# x = a row in sol
	if (length(x) != 120){
	 	print("The age structure should include 20 ages,for 6 disease classes, giving 120 columns")
	}
	Sp <- x[1:20]/apply(x, 1, sum); colnames(Sp)<- seq(1:20)
	Itp <- x[21:40]/apply(x, 1, sum); colnames(Itp)<- seq(1:20)
	Ibp <- x[41:60]/apply(x, 1, sum); colnames(Ibp)<- seq(1:20)
	Icp <- x[61:80]/apply(x, 1, sum); colnames(Icp)<- seq(1:20)
	Rp <- x[81:100]/apply(x, 1, sum); colnames(Rp)<- seq(1:20)
	Rcp <- x[101:120]/apply(x, 1, sum); colnames(Rcp)<- seq(1:20)
	
	mat <- as.matrix(rbind(Sp, Itp, Ibp, Icp, Rp, Rcp))
	mat[is.na(mat)]<- 0
	barplot(mat, # columns = age, rows = proportions
	xlab = 'age', main = "Population structure")	
}


#############################################################
#############################################################

# Figure out parameters that give reasonable age structure 
# Test plots, No Disease
#############################################################
S0 = 500*relage; It0 = 0*relage; Ib0 = 0*relage; 
Ic0 = 0*relage; R0 = 0 * relage; Rc0 = 0 * relage
x0 = c(S0, It0, Ib0, Ic0, R0, Rc0)
times <- seq(0, 100, 1)
params.test = c(fixed.params, list(gamma=1/2, betaB = 0.001,
	betaT = 0.001, rhoT = 1.2, rhoB = 4))
params.test.recov = c(fixed.params.recov, list(gamma=1/2, 
	betaB = 0.001, betaT = 0.001, rhoT = 1.2, rhoB = 4))


sol <- as.data.frame(ode(x0, times, rhs_age_matrix_ricker, params.test))
sol.recov<- as.data.frame(ode(x0, times, rhs_age_matrix_ricker, params.test.recov))

par(mfrow = c(1,2))
plot(x = sol$time, y = apply(sol[c(2:120)], 1, sum), 
	pch = 19, main = "Density Dependent, no Recovery")
plot(x = sol.recov$time, y = apply(sol.recov[c(2:120)], 1, sum), 
	pch = 19, main = "Density Dependent, no Recovery")

# visualize age distribution: 
as.matrix(sol[101,c(2:121)], dimnames = NULL)

par(mfrow = c(2, 2))
plot_agestructure(as.matrix(sol[5,c(2:121)]))
plot_agestructure(as.matrix(sol[25,c(2:121)]))
plot_agestructure(as.matrix(sol[50,c(2:121)]))
plot_agestructure(as.matrix(sol[101,c(2:121)]))  # better, no density dependence...


# Test with disease... make freq dependent later...
sol <- as.data.frame(ode(x0dz, times, rhs_age_matrix_ricker, params.test))
plot_dz_agestructure(sol[101, c(2:121)])
x0dz <- x0; x0dz[ib_index] <- 10
soldz <- as.data.frame(ode(x0dz, times, rhs_age_matrix_ricker, params.test))
plot_dz_agestructure(soldz[101, c(2:121)])


# Figure ure out birth rates taht give age structure above... 
f = function(b, x){
	1 * exp(-b * x)  # x = Nall, b = scale constant
}
# Number of offspring produced (assuming stable equal births/age) at K
300* f(1/100, 300) # 14.93612
300* f(1/500, 300)  # 164.6435
300* f(1/1000, 300)  # 222.2455
# Need to make that = deaths at stable age distribution...

sol <- as.data.frame(ode(x0, times, rhs_age_matrix_ricker, params.test))
sol.recov<- as.data.frame(ode(x0, times, rhs_age_matrix_ricker, params.test.recov))

par(mfrow = c(1,2))
plot(x = sol$time, y = apply(sol[c(2:120)], 1, sum), 
	pch = 19, main = "Density Dependent, no Recovery")
plot(x = sol.recov$time, y = apply(sol.recov[c(2:120)], 1, sum), 
	pch = 19, main = "Density Dependent, no Recovery")

# visualize age distribution: 
plot_agestructure(sol[101,c(2:120)])
# get eignevectors for disease free system!





	
#############################################################
#############################################################
#2) Test plots and grooming functions
#############################################################
#############################################################




# functions required to interpret results
#############################################################
plot_raw_numbers_mat = function(sol){
	plot(sol$time, apply(sol$S[s_index], 1, sum), col= "black",
		type= 'l', ylim = c(0, 1200), ylab = "Number of animals", 
		xlab = "Time (in years)")
	lines(sol$time, apply(sol[it_index], 1, sum), col= "red")
	lines(sol$time, apply(sol[ib_index], 1, sum), col= "blue")
	lines(sol$time, apply(sol[ic_index], 1, sum), col= "green")
	lines(sol$time, apply(sol[r_index], 1, sum), col = "orange")
	lines(sol$time, apply(sol[rc_index], 1, sum), col = "pink")
	legend("topright", legend = c("S", "It", "Ib", "Ic", "R", "Rc"),
		col = c("black", "red", "blue", "green", "orange", "pink"), 
		bty="n", lty = 1)
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




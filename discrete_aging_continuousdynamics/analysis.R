#############################################################
#############################################################
# Erin Gorsich
# This Code analyzes the stage-structured co-infection model,
# it uses discrete aging, defined in rhs, parameters estimated in fit_TBbruc
#############################################################
#############################################################
#############################################################

#############################################################
#############################################################
#1) Load fixed parameters, model
#############################################################
#############################################################
rm(list = ls())
require("deSolve")
library("plyr")
library("ggplot2")
library("lattice") # for levelplots
library("gridExtra") # layout for lattice
library("RColorBrewer")
set.seed(5)

setwd("~/GitHub/bTB-bruc-co-infection-ms/discrete_aging_continuousdynamics")
# get fixed.params (assuming no recovery)
source('fixed_parameters.R', chdir = TRUE)
# rhs function, determinitic model, age structure
source('rhs_age.R', chdir = TRUE)

# age divisions in rhs function
s.index <- 1:20
it.index <- 21:40
ib.index <- 41:60
ic.index <- 61:80
r.index <- 81:100
rc.index <- 101:120

#############################################################
#############################################################
#2) Analyses plotting Ro and EE using MC simulation from stats
#############################################################
#############################################################
#(previously in Ro_calculations)



#############################################################
#############################################################
#3) Analyses varying mortality and transmission rates
#############################################################
#############################################################
#(previously in run_TBbruc_agemodel)
rhoB_test <- seq(0, 8, length.out = 101)
rhoT_test <- seq(0, 8, length.out = 101)
mort_test <- seq(0, 15, length.out = 101) 

# Data frame to hold results of changing bruc effects on bTB
epiTB <- data.frame(
	rhoT= rep(rhoB_test, length(rhoB_test)), 
	mort = rep(mort_test, each = length(rhoB_test)), 
	bTBprev = NA, brucprev = NA, 	finalN = NA, 
	bTB_inS = NA, bTB_inB = NA, bruc_inS = NA, bruc_inTB = NA)

# Data frame to hold results of changing bTB effects on bruc
epiB <- data.frame(
	rhoB= rep(rhoB_test, length(rhoB_test)), 
	mort = rep(mort_test, each = length(rhoB_test)),
	bTBprev = NA, brucprev = NA, finalN = NA, 
	bTB_inS = NA, bTB_inB = NA, bruc_inS = NA, bruc_inTB = NA)

imax <- c(length(epiTB[,1]))
pb <- txtProgressBar(min = 0, max = imax, style = 3)	
for (i in 1:length(epiTB[,1])){
	params.test = c(fixed.params, list(gamma=1/2, theta = 4, K = 433,
		betaB = 0.6087396, betaT = 0.0012974553, rhoT = epiTB$rhoT[i], rhoB = 2.1))
	params.test$muC <- epiTB$mort[i] * params.test$muS
	params.test$muC[params.test$muC > 1] <- 1
	params.test$muRC <- epiTB$mort[i] * params.test$muS
	params.test$muRC[params.test$muRC > 1] <- 1
		
	sol <- as.data.frame(ode(x0, times, rhs_age_matrix, params.test))
	temp <- get_prevalence(sol)
	
	epiTB$bTBprev[i] = temp$prevTB
	epiTB$brucprev[i] = temp$prevB 	
	epiTB$finalN[i] = sum(sol[length(sol), c(2:121)])
	epiTB$bTB_inS[i] = temp$prevTinS 
	epiTB$bTB_inB[i] = temp$prevTinB
	epiTB$bruc_inS[i] = temp$prevBinS 
	epiTB$bruc_inTB[i] = temp$prevBinT
	rm(params.test, sol, temp)
	setTxtProgressBar(pb, i)
}
cat("\n")
summary(epiTB)

write.csv(epiTB, "~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/draft2/post-labmeeting/epiT.csv")


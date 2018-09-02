#############################################################
#############################################################
# Erin Gorsich
# This code runs a pde co-infection model and fits it to prevalence data
# It assumes density dependent brucellosis transmission
#############################################################
#############################################################
#############################################################
# Outline:
# 1) Load fixed parameters, model, set-up features of aging
# 2) Test plots - no disease / one disease ect. 
# 3) Fit - Define true prevalence and objective functions
# 4) Recreate a quick version of figure S6
#############################################################
#############################################################
#############################################################

rm(list = ls())
require("deSolve")
library("plyr")
library("pracma")
library("doParallel")
library("foreach")

setwd("~/GitHub/bTB-bruc-co-infection-ms/pde")


#############################################################
#############################################################
# 1) Parameters, model definition
#############################################################
#############################################################

# Parameters
#############################################################
# mortality
muS <- NA
muS <- c(1- 0.9,  1- 0.94, 1- 0.94, 1- 0.9) 
muS <- muS
dmuT <- 2.8
dmuB <- 3.03
dmuC <- 8.56 

# births- b is the maximum possible birth rate in S for age >=5
b <- 0.5
b1 <- 1; b2 <- 1; b3 <- 1; b4 <- 1; b5 <- 1

# density dependence
theta= 4
K = 433

# disease
epsilon <- 0.03

# fixed.params
p <- list(b = b, muS = muS, dmuB = dmuB, dmuT = dmuT, 
          dmuC =dmuC, epsilon = epsilon, K = K, theta = theta)

gen_fixed_params = function(agemax, agestep, p, recovery = FALSE) {
    N <- agemax / agestep
    ages <- seq(1, agemax + 1, by = agestep)[-(N+1)]
    juv.index <- which(ages <= 2)
    adult.index <- which(ages > 2 & ages <= 16)
    sens.index <- which(ages > 16)
    fecund.index <- which(ages >= 5)	
    
    # define age fluxes used in rhs
    da <- diff(c(min(ages) - agestep, ages))
    aging <- diag(-1/da)
    aging[row(aging) - col(aging) == 1] <- 1 / head(da, -1)
    
    #Mortality vector
    muS <- NA; muT <- NA; muB <- NA; muC <- NA
    muS[juv.index] <- p$muS[1]
    muS[adult.index] <- p$muS[2]
    muS[sens.index] <- p$muS[4]
    muT <- p$dmuT * muS
    muB <- p$dmuB * muS
    muC <- p$dmuC * muS
    if (recovery == FALSE) {
        muR <- muB
        muRC <- muC
    }
    if (recovery == TRUE) {
        muR <- muS
        muRC <- muT
    }
    
    # Fecundity
    b <- rep(0, length(ages))
    b[fecund.index] <- p$b
    
    fixed.params = list(b = b, aging = aging, ages = ages, 
        muS = muS, muB = muB, muT = muT, muC = muC, muR = muR, 
        muRC = muRC, epsilon = epsilon, K = K, theta = theta)
}

# Model definition - beverton & holt descritized PDE 
# only difference from rhs in main file is in lambdaB; lambdapB
#############################################################
rhs = function(times, x, params){
    ##########################
    # co-infection model pde - 
    # x = 6*length(ages) long initial conditions
    # params list containing features of aging
    ##########################
    with(as.list(c(x, params)), {
        n.ages <- length(ages) # total number of bins
        
        # Define states
        S = x[1:n.ages] 				
        It = x[seq(n.ages + 1, 2 * n.ages)] 			
        Ib = x[seq(2 * n.ages + 1, 3 * n.ages)]
        Ic = x[seq(3 * n.ages + 1, 4 * n.ages)]
        R = x[seq(4 * n.ages + 1, 5 * n.ages)]
        Rc = x[seq(5 * n.ages + 1, 6 * n.ages)]
        
        # Population size (N)
        Nall <- sum(S + It + Ib + Ic + R + Rc)  # overall
        N <- S + It + Ib +Ic + R + Rc 			# by age category
        
        # Frequency dependent force of infection is age dependent 
        betaBm <- matrix(nrow = n.ages, ncol = n.ages)
        betaBm[1:n.ages, 1:n.ages] <- betaB
        dims <- which(ages >= 2 & ages <=5)
        betaBm[dims,] <- exp(0.885) * betaB
        betaTm <- matrix(nrow = n.ages, ncol = n.ages)
        betaTm[1:n.ages, 1:n.ages] <- betaT
        
        lambdaT <- betaTm %*% (It + Ic + Rc) 
        lambdaB <- betaBm %*% (Ib + Ic)
        lambdapT <- rhoT * betaTm %*% (It + Ic + Rc)
        lambdapB <- rhoB * betaBm %*% (Ib + Ic)
        
        Nb <- N
        
        birth <- c(b %*% Nb)
        recruitment <- c(birth / ( 1 + (Nall/K)^theta), rep(0, times = n.ages - 1))
        dS <- recruitment + aging %*% S - (lambdaT + lambdaB) * S - muS * S
        dIt <- lambdaT * S - (lambdapB + muT) * It + aging %*% It
        dIb <- lambdaB * S + aging %*% Ib +	epsilon * R - (gamma + lambdapT + muB) * Ib
        dIc <- lambdapT*Ib + lambdapB*It + aging %*% Ic + epsilon * Rc - (gamma + muC)*Ic
        dR <- gamma * Ib - (epsilon + muR + lambdapT) * R + aging %*% R
        dRc <- lambdapT * R + gamma * Ic + aging %*% Rc - (epsilon + muRC) * Rc
        out = list(c(dS, dIt, dIb, dIc, dR, dRc))
        return(out)	
    }
    )	
}

# Set up features of aging for plotting
#############################################################
agemax <- 20
agestep <- 0.1
N <- agemax / agestep
ages <- seq(1, agemax + 1, by = agestep)[-(N)]
binsize <- N / agemax
N == length(ages)
s.index <- 1:N
it.index <- seq(N+1, 2*N)
ib.index <- seq((2*N+1), 3*N)
ic.index <- seq((3*N+1), 4*N)
r.index <- seq((4*N+1), 5*N)
rc.index <- seq((5*N+1), 6*N)

# generate parameters with correct agebins
f.params <- gen_fixed_params(agemax, agestep, p = p, recovery = FALSE)

# Functions for plotting (and define indecies based on ages, N): 
source('~/GitHub/bTB-bruc-co-infection-ms/pde/plotting_functions.R', chdir = TRUE)

# Starting agestructure (Jolles 2007; Caron et al. 2001)
juv <- rep(0.137 / length(ages[ages < 2]), length(ages[ages < 2]))
sa <- rep(0.368 / length(ages[ages >= 2 & ages < 6]), length(ages[ages >= 2 & ages < 6]))
a <- rep(0.185 / length(ages[ages >= 6 & ages < 9]), length(ages[ages >= 6 & ages < 9]))
ma <- rep(0.235 / length(ages[ages >= 9 & ages < 14]), length(ages[ages >= 9 & ages < 14]))
sen <- rep(0.075 / length(ages[ages >= 14 ]), length(ages[ages >= 14]))

relage <- c(juv, sa, a, ma, sen); length(relage) == N									
plot.relage <- c(0.137, rep(0.368/4, 4), rep(0.185/3, 3), rep(0.235/5, 5), rep(0.075/7, 7)) 	

#############################################################
#############################################################
# 2) Test plots - no disease/single infection ect
#############################################################
#############################################################

# No disease (as a guess, start by divinding betaB by K)
#############################################################
S0 <- 400 * relage; It0 <- 0 * relage; Ib0 <- 0 * relage
Ic0 <- 0 * relage; R0 <- 0 * relage; Rc0 <- 0 * relage
times <- seq(1, 500, 1)
x0 <- c(S0, It0, Ib0, Ic0, R0, Rc0)
params <- c(f.params, list(gamma = 1/2, betaB = 0.6087396/433, 
    betaT = 0.0012974553, rhoT = 1, rhoB = 2.1))
test <- as.data.frame(ode.1D(x0, times, rhs, params, 
    nspec = 6, dimens = N, method = "ode45"))

par(mfrow = c(2,2))
plot_agestructure(test, t = 100)
plot_agestructure(test, t = 500)
plot(x = test$time, y = apply(test[c(2:length(colnames(test)))], 1, sum),
    pch = 19, main = "sum")
plot_raw_numbers(test)

stable_age <- unname(unlist( test[length(test[,1]), c(2:(length(ages)+1))] / 
    sum(test[length(test[,1]), c(2:(length(ages)+1))]) ))

# bTB only
#############################################################
S0 <- 400 * stable_age; It0 <- 10 * stable_age; Ib0 <- 0 * stable_age
Ic0 <- 0 * stable_age; R0 <- 0 * stable_age; Rc0 <- 0 * stable_age
times <- seq(1, 500, 1)
x0 <- c(S0, It0, Ib0, Ic0, R0, Rc0)

betaT_temp <- seq(0.0001, 0.001, 0.00005) 
prevT <- NA; prevTrecov <- NA
prevB <- NA
for (i in 1:length(betaT_temp)){
    params.test = c(f.params, list(gamma=0.5, 
        betaB = 0.6087396, betaT = betaT_temp[i], rhoT = 1, rhoB = 2.1))
    sol <- as.data.frame(ode.1D(x0, times, rhs, params.test, 
        nspec = 6, dimens = N, method = "ode45"))
    prevT[i]<- get_prevalence(sol)$prevTB
    rm(sol, params.test)
}

plot(x = betaT_temp, y = prevT, type = "l", xlab = expression(beta), 
     ylab = "BTB prevalence", main = "Single disease")
abline(h = c(0.1, 0.2, 0.3, 0.4), col = "dark red")
# Beta value at:
betaT_temp[which(prevT < 0.12 & prevT > 0.08)]		
betaT_temp[which(prevT < 0.21 & prevT > 0.19)]		
betaT_temp[which(prevT < 0.32 & prevT > 0.28)]
betaT_temp[which(prevT < 0.405 & prevT > 0.395)]

# 30% prev at ~0.0006

# brucellosis only
#############################################################
S0 <- 400* stable_age; It0 <- 0 * stable_age; Ib0 <- 20* stable_age; 
Ic0 <- 0* stable_age; R0 <- 30 * stable_age; Rc0 <- 0 * stable_age
x0 <- c(S0, It0, Ib0, Ic0, R0, Rc0)

betaB_temp <- seq(0.4, 1, 0.1) / 433 # changed for fd case
prevB <- NA; prevBrecov <- NA; prevT <- NA
for (i in 1:length(betaB_temp)){
    params.test = c(f.params, list(gamma=0.5, 
        betaB = betaB_temp[i], betaT = 0.00129, rhoT = 1, rhoB = 2.1))
    sol <- as.data.frame(ode.1D(x0, times, rhs, params.test, 
        nspec = 6, dimens = N, method = "ode45"))
    prevB[i]<- get_prevalence(sol)$prevB
    prevT[i]<- get_prevalence(sol)$prevTB
    get_prevalence(sol)$prevB
    rm(sol, params.test)
}

plot(x = betaB_temp, y = prevB, type = "l", xlab = expression(beta), 
     ylab = "Brucellosis prevalence", main = "Single disease")
abline(h = c(0.1, 0.2, 0.3, 0.4), col = "dark red")
# Beta value at:
# 30% around 0.0013

params.test = c(f.params, list(gamma=0.5, 
    betaB = 0.0013, betaT = 0.00129, rhoT = 1, rhoB = 2.1))
sol <- as.data.frame(ode.1D(x0, times, rhs, params.test, 
    nspec = 6, dimens = N, method = "ode45"))
xB <- unname(unlist( sol[length(sol[,1]), c(2:length(sol))] / 
    sum(test[length(sol[,1]), c(2:length(sol))]) ))

#############################################################
#############################################################
# 3) Fit - define true prevalence and objective function
#############################################################
#############################################################
prevTBobs <- 0.27  # for test- bootstrap estimate of overall prevalence
prevBobs <- 0.34
data <- read.csv("~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/final_datasets_copied_from_phdfolder/cross_sectional_data_withdz_cleandisease_nofinal_Feb2016_capturetime_forsurv.csv")
counts<- hist(data$age_sel/12, plot = FALSE)$counts  # youngest = 1.4 so aged 1-2
agestructure<- counts/sum(counts)
agestructure_yr = c(agestructure, 0, 0, 0, 0, 0)
data_agestructure <- agestructure_yr

objective = function(params.est){
    # params.est = 2 long = c(betaB, betaT)
    # now both betaT and betaB are rescaled
    params <- c(f.params, list(gamma = 1/2, betaB = params.est[1]/1000, 
        betaT = params.est[2]/1000, rhoT = 1, rhoB = 2.1))
    
    # seed from endemic brucellosis conditions, 10 bTB positive buffalo
    x0 = xB
    x0[[min(it.index) + 5*binsize]] <- 5
    x0[[min(it.index) + 5*binsize + 1]] <- 5
    sol <- as.data.frame(ode.1D(x0, times, rhs, params, 
                                nspec = 6, dimens = N, method = "ode45"))
    df <- get_structured_prevalence(sol)
    error <- sqrt(((prevTBobs - df$prevTB)^2 + (prevBobs - df$prevB)^2))
    return (error)
}

par <- optim(c(0.0013*1000, 0.00094*1000), objective); par
par2 <- optim(c(0.01*1000, 0.0007*1000), objective); par2
par3 <- optim(c(0.0001*1000, 0.001*1000), objective); par3
#$par
#[1] 1.678073 1.330546
# par2
#1.678073 1.330546
# par3
#0.2127101 0.6208829

# test
S0 <- 400* stable_age; It0 <- 10 * stable_age; Ib0 <- 10* stable_age; 
Ic0 <- 0* stable_age; R0 <- 0 * stable_age; Rc0 <- 0 * stable_age
x0 <- c(S0, It0, Ib0, Ic0, R0, Rc0)
params <- c(f.params, list(gamma = 1/2, betaB = 1.678073/1000, 
    betaT = 1.330546/1000, rhoT = 1, rhoB = 2.1))
sol <- as.data.frame(ode.1D(x0, times, rhs, params, 
    nspec = 6, dimens = N, method = "ode45"))

get_structured_prevalence(sol)
plot_raw_numbers(sol)


#############################################################
#############################################################
# 4) Ugly version of figure 3
#############################################################
#############################################################
get_EE = function(params, x0, method){
    ###################################
    # Input: x = c(S = final, 1*20 S vector at params)
    # method = "logistic", "ricker", "beverton-holt"
    # Output: Ro
    ###################################
    # stable age in disease free conditions now read in
    ###################################
    binsize <- N / agemax #number of delta age bins to get one age bin...
    
    # run single infection EE, bTB
    ###################################
    x_singleTB <- x0
    x_singleTB[(N+1)  +  3*binsize] <- 1
    x_singleTB[(N+1) +  4*binsize] <- 1
    if (method == "ricker"){
        sol <- as.data.frame(ode.1D(x_singleTB, times, 
                                    rhs_ricker, params, nspec = 6, dimens = N, 
                                    method = "ode45"))}
    if (method == "logistic"){
        sol <- as.data.frame(ode.1D(x_singleTB, times, 
                                    rhs_logistic, params, nspec = 6, dimens = N,
                                    method = "ode45"))}
    if (method == "beverton-holt"){
        sol <- as.data.frame(ode.1D(x_singleTB, times, 
                                    rhs, params, nspec = 6, dimens = N, 
                                    method = "ode45"))}
    S <-sum(sol[length(sol[,1]) , s.index+1])
    It <- sum(sol[length(sol[,1]) , it.index +1])
    Ic <- sum(sol[length(sol[,1]) , ic.index +1])
    Ib <- sum(sol[length(sol[,1]) , ib.index +1])
    Ic <- sum(sol[length(sol[,1]) , ic.index +1])
    R <- sum(sol[length(sol[,1]) , r.index +1])
    Rc <-sum(sol[length(sol[,1]) , rc.index +1])
    Tot <- sum(sol[length(sol[,1]), 2:length(colnames(sol))])	
    EE_bTB_alone <- (It + Ic + Rc) / Tot
    rm(sol)
    
    # run single infection EE, bruc
    ###################################
    x_singleBruc <- x0
    x_singleBruc[min(ib.index) + 1 + 3*binsize] <- 1
    x_singleBruc[min(ib.index) + 1 + 4*binsize] <- 1
    
    if (method == "ricker"){
        sol <- as.data.frame(ode.1D(x_singleBruc, times, rhs_ricker, 
                                    params, nspec = 6, dimens = N, method = "ode45"))}
    if (method == "logistic"){
        sol <- as.data.frame(ode.1D(x_singleBruc, times, rhs_logistic, 
                                    params, nspec = 6, dimens = N, method = "ode45"))}
    if (method == "beverton-holt"){
        sol <- as.data.frame(ode.1D(x_singleBruc, times, rhs,
                                    params, nspec = 6, dimens = N, method = "ode45"))}	
    S <-sum(sol[length(sol[,1]) , s.index+1])
    It <- sum(sol[length(sol[,1]) , it.index +1])
    Ic <- sum(sol[length(sol[,1]) , ic.index +1])
    Ib <- sum(sol[length(sol[,1]) , ib.index +1])
    Ic <- sum(sol[length(sol[,1]) , ic.index +1])
    R <- sum(sol[length(sol[,1]) , r.index +1])
    Rc <-sum(sol[length(sol[,1]) , rc.index +1])
    Tot <- sum(sol[length(sol[,1]), 2:length(colnames(sol))])	
    EE_bruc_alone <- (Ib + Ic + R + Rc) / Tot
    
    x_endB <- unname(unlist(sol[length(sol[,1]), 2:length(colnames(sol))] ))
    rm(sol)
    
    # introduce both...
    ###################################
    x_endB[min(it.index) + 1  +  3*binsize] <- 1
    x_endB[min(it.index) + 1  +  4*binsize] <- 1
    
    if (method == "ricker"){
        sol <- as.data.frame(ode.1D(x_endB, times, rhs_ricker, 
                                    params, nspec = 6, dimens = N, method = "ode45"))}
    if (method == "logistic"){
        sol <- as.data.frame(ode.1D(x_endB, times, rhs_logistic, 
                                    params, nspec = 6, dimens = N, method = "ode45"))}
    if (method == "beverton-holt"){
        sol <- as.data.frame(ode.1D(x_endB, times, rhs, 
                                    params, nspec = 6, dimens = N, method = "ode45"))}	
    
    S <-sum(sol[length(sol[,1]) , s.index+1])
    It <- sum(sol[length(sol[,1]) , it.index +1])
    Ic <- sum(sol[length(sol[,1]) , ic.index +1])
    Ib <- sum(sol[length(sol[,1]) , ib.index +1])
    Ic <- sum(sol[length(sol[,1]) , ic.index +1])
    R <- sum(sol[length(sol[,1]) , r.index +1])
    Rc <-sum(sol[length(sol[,1]) , rc.index +1])
    Tot <- sum(sol[length(sol[,1]), 2:length(colnames(sol))])	
    
    EE_bTB_co <- (It + Ic + Rc) / Tot
    EE_bruc_co <- (Ib + Ic + R + Rc) / Tot
    rm(sol)
    
    return(c(EE_bTB_alone, EE_bTB_co, EE_bruc_alone, EE_bruc_co))
}


# MCMC for endemic equibrilium
#############################################################
# generate n random samples 
n = 500
set.seed(1)

S0 <- 400* stable_age; It0 <- 0 * stable_age; Ib0 <- 0* stable_age; 
Ic0 <- 0* stable_age; R0 <- 0 * stable_age; Rc0 <- 0 * stable_age
x0 <- c(S0, It0, Ib0, Ic0, R0, Rc0)

# Generate 1000 samples
cl <- makeCluster(6)
# cl <- makeCluster(10) # work computer
registerDoParallel(cl)
d3 <- foreach(icount(n), .combine = rbind, .packages = "deSolve") %dopar% {
    params <- c(f.params, list(gamma = 1/2, betaB = 1.678073/1000,
        betaT = 1.330546/1000, rhoT = 1, rhoB = 2.1))
    params$rhoB<- exp(rnorm(n = 1, mean = 0.75579, sd = 0.40714))
    B <- rnorm(n = 1, mean = 1.1060, sd = 0.3505)
    TB <- rnorm(n = 1, mean = 1.0370, sd = 0.3483)
    params$muB <- params$muS * exp(B)
    params$muT <- params$muS * exp(TB)
    params$muC <- params$muS * exp(B + TB)
    params$muR <- params$muT
    params$muRC <- params$muC
    # think about runs with too high mortality!
    val <- get_EE(params, x0, "beverton-holt")
    data <- data.frame(
        rhoB = params$rhoB, 
        dB = params$muB[1]/ params$muS[1],
        dT = params$muT[1]/ params$muS[1],
        EE_bTB_single = val[1],  
        EE_bTB_co = val[2], 
        EE_brucellosis_single = val[3],
        EE_brucellosis_co = val[4])
}
stopCluster(cl)

saveRDS(d3, file = "EE_confidence_interval_simulation_results_revisiondensitydep.rds")


# Figure
#############################################################
d <- readRDS("~/GitHub/bTB-bruc-co-infection-ms/pde/sensitivity_densitydependence_results.rds")
d <- d[!is.na(d$rEE_bTB_co), ] # one NA

# test figure with violin plots
n <- length(d$EE_bTB_single)
n3 <- length(d3$EE_bTB_single)
df <- data.frame(E = c(d$EE_bTB_single, d$EE_bTB_co, 
    d$EE_brucellosis_single, d$EE_brucellosis_co, 
    d3$EE_bTB_single, d3$EE_bTB_co, 
    d3$EE_brucellosis_single, d3$EE_brucellosis_co),
    meanE = c(rep(median(d$EE_bTB_single),n), 
        rep(median(d$EE_bTB_co),n),
        rep(median(d$EE_brucellosis_single),n), 
        rep(median(d$EE_brucellosis_co),n),
        rep(median(d3$EE_bTB_single),n3), 
        rep(median(d3$EE_bTB_co),n3),
        rep(median(d3$EE_brucellosis_single),n3), 
        rep(median(d3$EE_brucellosis_co),n3)
        ),
    sdE = c(rep(sd(d$EE_bTB_single),n), 
        rep(sd(d$EE_bTB_co),n),
        rep(sd(d$EE_brucellosis_single),n), 
        rep(sd(d$EE_brucellosis_co),n),
        rep(sd(d3$EE_bTB_single),n3), 
        rep(sd(d3$EE_bTB_co),n3),
        rep(sd(d3$EE_brucellosis_single),n3), 
        rep(sd(d3$EE_brucellosis_co),n3)),
    infection = c(rep("bTB", 2*n), rep("brucellosis", 2*n),
        rep("bTB", 2*n3), rep("brucellosis", 2*n3)), 
    model = c(rep("frequency-dependent", 4*n), rep("density-dependent", 4*n3)),
    singleco = c(
        rep(c("single", "co-infection", "single","co-infection"), each = n), 
        rep(c("single", "co-infection", "single","co-infection"), each = n3)),
    X = c(
        rep(c("bTB-single-F", "bTB-co-infection-F", "bruc-single-F", 
            "bruc-co-infection-F"), each = n), 
        rep(c("bTB-single-D", "bTB-co-infection-D", "bruc-single-D", 
            "bruc-co-infection-D"), each = n3)) )
df$X <- factor(df$X, 
    levels = c("bTB-single-F", "bTB-co-infection-F", "bruc-single-F", 
    "bruc-co-infection-F", "bTB-single-D", "bTB-co-infection-D", 
    "bruc-single-D", "bruc-co-infection-D"))
df$infection <- relevel(as.factor(df$infection), "bTB")
df$singleco <- relevel(as.factor(df$singleco), "single")
df <- df[!is.na(df$E),]
df <- df[!is.na(df$sdE),] # none removed!
df$infectionmodel <- paste(df$infection, df$model, sep = "_")

df$singleco <- as.character(df$singleco)
df$singleco[df$singleco == "single"] <- "one infection"
df$singleco[df$singleco == "co-infection"] <- "both infections"
df$singleco <- as.factor(df$singleco)
df$singleco <- relevel(as.factor(df$singleco), "one infection")

brucellosis <- df[df$infection  == "brucellosis",]
tb <- df[df$infection == "bTB",]


# Dot and error plots for Endemic Prevalence
pB <- ggplot(brucellosis, aes(x = model, y = E, shape = singleco, colour = singleco)) +
    geom_point(aes(x = model, y = meanE), size = 3,
               position= position_dodge(width = 0.5)) +
    geom_errorbar(aes(x = model, ymin = meanE - sdE, ymax = meanE + sdE),
                  width = 0, position= position_dodge(width = 0.5)) +
    ylim(0, 0.83) +
    xlab("") +
    ylab("Brucellosis prevalence") +
    labs(shape = "Populations with") +
    theme_bw() +
    scale_colour_manual(values = c("chartreuse4","chartreuse4"), guide = F) +
    scale_shape_manual(values = c(19, 17)) +
    theme(axis.line.x = element_line(colour= "black"),
          axis.line.y = element_line(colour= "black"),
          axis.title.x = element_text(size=16, vjust=-0.15),
          axis.title.y = element_text(size=16, vjust= 0.8, margin = margin(r = 8)),
          axis.text.x = element_text(size=16, vjust=-0.05, margin = margin(t = 7)),
          axis.text.y = element_text(size=14, margin = margin(r = 4)),
          panel.border = element_blank(), 
          legend.position= c(0.8, 0.9),  
          legend.text = element_text(size = 14),
          legend.background= element_rect(fill="white", colour="white"),
          legend.key= element_blank(),
          legend.title= element_text(size = 14))

pT <- ggplot(tb, aes(x = model, y = E, shape = singleco, colour = singleco)) +
    geom_point(aes(x = model, y = meanE), size = 3,
               position= position_dodge(width = 0.5)) +
    geom_errorbar(aes(x = model, ymin = meanE - sdE, ymax = meanE + sdE),
                  width = 0, position= position_dodge(width = 0.5)) +
    ylim(0, 0.83) +
    xlab("") +
    ylab("bTB prevalence") +
    theme_bw() +
    scale_colour_manual(values = c("slateblue3","slateblue3"), guide = FALSE) +
    scale_shape_manual(values = c(19, 17),  guide = FALSE) +
    theme(axis.line.x = element_line(colour= "black"),
          axis.line.y = element_line(colour= "black"),
          axis.title.x = element_text(size=16, vjust=-0.15),
          axis.title.y = element_text(size=16, vjust= 0.8, margin = margin(r = 8)),
          axis.text.x = element_text(size=16, vjust=-0.05, margin = margin(t = 7)),
          axis.text.y = element_text(size=14, margin = margin(r = 4)),
          panel.border = element_blank())
#legend.position=c(0.85, 0.9),  
#legend.text = element_text(size = 10),
#legend.background= element_rect(fill="white", colour="white"),
#legend.key= element_blank(),
#legend.title= element_blank())

source('~/GitHub/bTB-bruc-co-infection-ms/pde/multiplot.R', chdir = TRUE)
multiplot(pT,  pB, cols = 2)
# 800*400


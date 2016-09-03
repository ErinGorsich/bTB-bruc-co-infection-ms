#############################################################
#############################################################
3) Impliment Stochastic, co-infection model
#############################################################
#############################################################
coinfect.onestep <- function(x, params){
	S <- x[2]
	It <- x[3]
	Ib <- x[4]
	R <- x[5]
	Ic <- x[6]
	Rc <- x[7]
	N <- S + It + Ib + R + Ic +Rc
	with(as.list(params),
	{
	total.rate <- epsilon*Rc + gamma*Ic + epsilon*R + gamma*Ib +  
		betapT*(It+Ic+Rc)*R + betapT*(It+Ic+Rc)*Ib + betapB*(Ib+Ic)*It + 
		betaB*(Ib+Ic)*S + betaT*(It+Ic+Rc)*S + 
		muRC*Rc + muC*Ic + muR*R + muB*Ib + muT*It + muS*S +
	 	b*(S+b1*It+b2*Ib+b3*R+b4*Ic+b5*Rc)*(1-(r/b)*(N/K))  
	tau <- rexp(n=1, rate = total.rate)	 	# calculate inter-event time
	new.statevar <- c(S, It, Ib, R, Ic, Rc)	# local state variables 
	U <- runif(1)
	new.statevar <- c(S+1, It, Ib, R, Ic, Rc) 	# default = births
		
	# Deaths
	if (U <= (epsilon*Rc + gamma*Ic + epsilon*R + gamma*Ib +  
		betapT*(It+Ic+Rc)*R + betapT*(It+Ic+Rc)*Ib + 
		betapB*(Ib+Ic)*It + betaB*(Ib+Ic)*S + betaT*(It+Ic+Rc)*S + 
		muRC*Rc + muC*Ic + muR*R + muB*Ib + muT*It + muS*S)/ total.rate){
		new.statevar <- c(S-1, It, Ib, R, Ic, Rc)}
	if (U <= (epsilon*Rc + gamma*Ic + epsilon*R + gamma*Ib +  
		betapT*(It+Ic+Rc)*R + betapT*(It+Ic+Rc)*Ib + 
		betapB*(Ib+Ic)*It + betaB*(Ib+Ic)*S + betaT*(It+Ic+Rc)*S + 
		muRC*Rc + muC*Ic + muR*R + muB*Ib + muT*It)/ total.rate){
		new.statevar <- c(S, It-1, Ib, R, Ic, Rc)}
	if (U <= (epsilon*Rc + gamma*Ic + epsilon*R + gamma*Ib +  
		betapT*(It+Ic+Rc)*R + betapT*(It+Ic+Rc)*Ib + 
		betapB*(Ib+Ic)*It + betaB*(Ib+Ic)*S + betaT*(It+Ic+Rc)*S + 
		muRC*Rc + muC*Ic + muR*R + muB*Ib)/ total.rate){
		new.statevar <- c(S, It, Ib-1, R, Ic, Rc)}
	if (U <= (epsilon*Rc + gamma*Ic + epsilon*R + gamma*Ib +  
		betapT*(It+Ic+Rc)*R + betapT*(It+Ic+Rc)*Ib + 
		betapB*(Ib+Ic)*It + betaB*(Ib+Ic)*S + betaT*(It+Ic+Rc)*S + 
		muRC*Rc + muC*Ic + muR*R) / total.rate){
		new.statevar <- c(S, It, Ib, R-1, Ic, Rc)}
	if (U <= (epsilon*Rc + gamma*Ic + epsilon*R + gamma*Ib +  
		betapT*(It+Ic+Rc)*R + betapT*(It+Ic+Rc)*Ib + 
		betapB*(Ib+Ic)*It + betaB*(Ib+Ic)*S + betaT*(It+Ic+Rc)*S + 
		muRC*Rc + muC*Ic)/ total.rate){
		new.statevar <- c(S, It, Ib, R, Ic-1, Rc)}
	if (U <= (epsilon*Rc + gamma*Ic + epsilon*R + gamma*Ib +  
		betapT*(It+Ic+Rc)*R + betapT*(It+Ic+Rc)*Ib + 
		betapB*(Ib+Ic)*It + betaB*(Ib+Ic)*S + betaT*(It+Ic+Rc)*S + 
		muRC*Rc) / total.rate){
		new.statevar <- c(S, It, Ib, R, Ic, Rc-1)}  
			
	# Transmission from Susceptibles
	if (U <= (epsilon*Rc + gamma*Ic + epsilon*R + gamma*Ib +  
		betapT*(It+Ic+Rc)*R + betapT*(It+Ic+Rc)*Ib + 
		betapB*(Ib+Ic)*It + betaB*(Ib+Ic)*S + betaT*(It+Ic+Rc)*S) / 
		total.rate){
		new.statevar <- c(S-1, It+1, Ib, R, Ic, Rc)}
	if (U <= (epsilon*Rc + gamma*Ic + epsilon*R + gamma*Ib +  
		betapT*(It+Ic+Rc)*R + betapT*(It+Ic+Rc)*Ib + 
		betapB*(Ib+Ic)*It + betaB*(Ib+Ic)*S) / total.rate){
		new.statevar <- c(S-1, It, Ib+1, R, Ic, Rc)}
	
	# transmission from first infection, recoverds
	if (U <= (epsilon*Rc + gamma*Ic + epsilon*R + gamma*Ib +  
		betapT*(It+Ic+Rc)*R + betapT*(It+Ic+Rc)*Ib + 
		betapB*(Ib+Ic)*It) / total.rate){
		new.statevar <- c(S, It-1, Ib, R, Ic+1, Rc)}
	if (U <= (epsilon*Rc + gamma*Ic + epsilon*R + gamma*Ib +  
		betapT*(It+Ic+Rc)*R + betapT*(It+Ic+Rc)*Ib) / total.rate){
		new.statevar <- c(S, It, Ib-1, R, Ic+1, Rc)}
	if (U <= (epsilon*Rc + gamma*Ic + epsilon*R + gamma*Ib +  
		betapT*(It+Ic+Rc)*R) / total.rate){
		new.statevar <- c(S, It, Ib, R-1, Ic, Rc+1)}
	
	# recovery and recrudescence
	if (U <= (epsilon*Rc + gamma*Ic + epsilon*R + gamma*Ib) /
		total.rate){
		new.statevar <- c(S, It, Ib-1, R+1, Ic, Rc)}	
	if (U <= (epsilon*Rc + gamma*Ic + epsilon*R) / total.rate){
		new.statevar <- c(S, It, Ib+1, R-1, Ic, Rc)}
	if (U <= (epsilon*Rc + gamma*Ic) / total.rate){
		new.statevar <- c(S, It, Ib, R, Ic-1, Rc+1)}	
	if (U <= (epsilon*Rc) / total.rate){
		new.statevar <- c(S, It, Ib, R, Ic+1, Rc-1)}

	# return, store result
	c(tau, new.statevar) 
	}
	)
}

coinfect.model <- function(x, params, nstep){
	output <- array(dim = c(nstep + 1, 7))		# array to store results
	colnames(output) <- c("time", "S", "It", "Ib", "R", "Ic", "Rc")
	output[1,] <- x 			# first output is initial conditions
	for(k in 1:nstep){
		output[k+1,] <- x <- coinfect.onestep(x, params) 
	}
	return(output)
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

run_stochastic_coinfection_model= function(params, nstep, nsims){
	set.seed(2)
	brucendemic <- get_starting_eqbruc(params)
	xstart <- c(time = 0, S = brucendemic[1][[1]]-10, 
		It = 10, Ib = brucendemic[3][[1]], Ic = 0, 
		R = brucendemic[5][[1]], Rc = 0) 
	data <- list(NA, length = nsims)
	for (i in 1:nsims){
		data[[i]] <- as.data.frame(coinfect.model(xstart, params, nstep))
		data[[i]]$cumtime <- cumsum(data[[i]]$time)
	}	
	return(data)
}

# agemax and N need to be global parameters
# indicies as well

get_EE = function(params, x0, method){
	###################################
	# Input: x = c(S = final, 1*20 S vector at params)
	# method = "logistic", "ricker", "beverton-holt"
	# Output: Ro
	###################################
	# stable age in disease free conditions now read in
	###################################
	binsize <- N/agemax #number of delta age bins to get one age bin...
	
	# run single infection EE, bTB
	###################################
	x_singleTB <- x0
	x_singleTB[(N+1)  +  3*binsize] <- 1
	x_singleTB[(N+1) +  4*binsize] <- 1
	if (method == "ricker"){
		sol <- as.data.frame(ode.1D(x_singleTB, times, 
			rhs_ricker, params, nspec = 6, dimens = N,method = "ode45"))}
	if (method == "logistic"){
		sol <- as.data.frame(ode.1D(x_singleTB, times, rhs_logistic, 
			params, nspec = 6, dimens = N, method = "ode45"))}
	if (method == "beverton-holt"){
		sol <- as.data.frame(ode.1D(x_singleTB, times, rhs, 
			nspec = 6, dimens = N, params, method = "ode45"))}
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
		sol <- as.data.frame(ode.1D(x_endB, times, rhs_ricker, params,  
			nspec = 6, dimens = N, method = "ode45"))}
	if (method == "logistic"){
		sol <- as.data.frame(ode.1D(x_endB, times, rhs_logistic, params,
			nspec = 6, dimens = N, method = "ode45"))}
	if (method == "beverton-holt"){
		sol <- as.data.frame(ode.1D(x_endB, times, rhs, params, 
			nspec = 6, dimens = N, method = "ode45"))}	
	
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



getEE = function(params){
	###################################
	# Input: x = c(S = final, 1*20 S vector at params)
	# Output: Ro
	###################################
	# Get stable age distribution in dz free conditions
	relage = c(0.137, rep(0.368/4, 4), rep(0.185/4, 4),
		rep(0.235/6, 6), rep(0.075/5, 5))
	S0 = 400*relage; It0 = 0*relage; Ib0 = 0*relage; 
	Ic0 = 0*relage; R0 = 0 * relage; Rc0 = 0 * relage
	x0 = c(S0, It0, Ib0, Ic0, R0, Rc0)
	times <- seq(0, 1000, 1)
	sol <- as.data.frame(ode(x0, times, rhs_age_matrix, params))
	stable_age <- unname(unlist( sol[length(sol[,1]), c(2:21)]/sum(sol[length(sol[,1]), c(2:21)]) ))
	S0 = 400 * stable_age; It0 = 0 * stable_age; Ib0 = 0* stable_age; 
	Ic0 = 0* stable_age; R0 = 0 * stable_age; Rc0 = 0 * stable_age
	x0 = c(S0, It0, Ib0, Ic0, R0, Rc0)
	
	# run single infection EE, bTB
	###################################
	x_singleTB = x0
	x_singleTB[22:24] <- x_singleTB[22:24] + 1
	x_singleTB[2:4] <- x_singleTB[2:4] - 1
	sol <- as.data.frame(ode(x_singleTB, times, rhs_age_matrix, params, method = "ode45"))
	
	S <-sum(sol[length(sol[,1]) , s_index+1])
	It <- sum(sol[length(sol[,1]) , it_index +1])
	Ic <- sum(sol[length(sol[,1]) , ic_index +1])
	Ib <- sum(sol[length(sol[,1]) , ib_index +1])
	Ic <- sum(sol[length(sol[,1]) , ic_index +1])
	R <- sum(sol[length(sol[,1]) , r_index +1])
	Rc <-sum(sol[length(sol[,1]) , rc_index +1])
	N <- sum(sol[length(sol[,1]), 2:121])	
	EE_bTB_alone <- (It + Ic + Rc) / N

	# run single infection EE, bruc
	###################################
	x_singleBruc = x0
	x_singleBruc[42:44] <- x_singleBruc[42:44] + 1
	x_singleBruc[2:4] <- x_singleBruc[2:4] - 1
	sol <- as.data.frame(ode(x_singleBruc, times, rhs_age_matrix, params, method = "ode45"))
	
	S <-sum(sol[length(sol[,1]) , s_index+1])
	It <- sum(sol[length(sol[,1]) , it_index +1])
	Ic <- sum(sol[length(sol[,1]) , ic_index +1])
	Ib <- sum(sol[length(sol[,1]) , ib_index +1])
	Ic <- sum(sol[length(sol[,1]) , ic_index +1])
	R <- sum(sol[length(sol[,1]) , r_index +1])
	Rc <-sum(sol[length(sol[,1]) , rc_index +1])
	N <- sum(S + It + Ib + Ic + R + Rc)
	EE_bruc_alone <- (Ib + Ic + R + Rc) / N

	x_endB <- unname(unlist(sol[length(sol[,1]), 2:121] ))

	# introduce both...
	###################################
	x_endB[22:24] <- x_endB[22:24] + 1
	if(x_endB[2] > 1){x_endB[2] <- x_endB[2] - 1 }
	if(x_endB[3] > 1){x_endB[3] <- x_endB[3] - 1 }
	if(x_endB[4] > 1){x_endB[4] <- x_endB[4] - 1 }

	sol <- as.data.frame(ode(x_endB, times, rhs_age_matrix, params, method = "ode45"))
	
	S <-sum(sol[length(sol[,1]) , s_index+1])
	It <- sum(sol[length(sol[,1]) , it_index +1])
	Ic <- sum(sol[length(sol[,1]) , ic_index +1])
	Ib <- sum(sol[length(sol[,1]) , ib_index +1])
	Ic <- sum(sol[length(sol[,1]) , ic_index +1])
	R <- sum(sol[length(sol[,1]) , r_index +1])
	Rc <-sum(sol[length(sol[,1]) , rc_index +1])
	N <- sum(S + It + Ib + Ic + R + Rc)

	EE_bTB_co <- (It + Ic + Rc) / N
	EE_bruc_co <- (Ib + Ic + R + Rc) / N
	
	return(c(EE_bTB_alone, EE_bTB_co, EE_bruc_alone, EE_bruc_co))
}

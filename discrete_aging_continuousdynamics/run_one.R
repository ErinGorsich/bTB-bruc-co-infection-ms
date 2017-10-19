run_one = function(tmax, x0, params) {
	####################################
	# Function to run model for until tmax,
	# update ages every 365 day
	# input = tmax (when sim stops) in days
	#x0 = initial vector (20*6 in length), 
	# list of parameters to run rhs_age
	####################################
		
	# initialize vectors to hold time and variables
	T <- c(); 
	S <- c(); It <- c(); Ib <- c()
	Ic <- c(); R <- c(); Rc <- c()
	S0 <- c(); It0 <- c(); Ib0 <- c()
	Ic0 <- c(); R0 <- c(); Rc0 <- c()
	
	# loop over until tmax
	T0 <- 0
	x.init <- x0
	while(T0 < tmax) {
		# add one year conditions one row at a time
		x <- ode(x.init, c(T0, T0 + 365), rhs, params)
		T <- rbind(T, x[2, 1])
		S <- rbind(S, x[2, s.index + 1])
		It <- rbind(It, x[2, it.index + 1])
		Ib <- rbind(Ib, x[2, ib.index + 1])
		Ic <-  rbind(Ic, x[2, ic.index + 1])
		R <-  rbind(R, x[2, r.index + 1])
		Rc <-  rbind(Rc, x[2, rc.index + 1])

		# age everyone in ages 1:19 up one
		S0[2:20] <- tail(x, 1)[2:20]
		It0[2:20] <- tail(x, 1)[22:40] 
		Ib0[2:20] <- tail(x, 1)[42:60] 
		Ic0[2:20] <- tail(x, 1)[62:80] 
		R0[2:20] <- tail(x, 1)[82:100] 
		Rc0[2:20] <- tail(x, 1)[102:120] 
		T0 <- tail(T, 1)
		S0[1] <- 0; It0[1] <- 0; Ib0[1] <- 0
		Ic0[1] <- 0; R0[1] <- 0; Rc0[1] <- 0
		x.init <- c(S0, It0, Ib0, Ic0, R0, Rc0)
	}
	out <- data.frame(cbind(T, S, It, Ib, Ic, R, Rc) )
	return(out)
}
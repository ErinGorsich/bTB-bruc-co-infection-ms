# age/disease divisions in rhs function
juv.index <- which(ages <= 2)
adult.index <- which(ages > 2 & ages <= 16)
sens.index <- which(ages > 16)
fecund.index <- which(ages >= 5)	
s.index <- 1:N
it.index <- seq(N+1, 2*N)
ib.index <- seq((2*N+1), 3*N)
ic.index <- seq((3*N+1), 4*N)
r.index <- seq((4*N+1), 5*N)
rc.index <- seq((5*N+1), 6*N)


plot_agestructure = function(sol, t){ # true based on Jolles et al. 2007
	 ############################
	 # x = raw solver output; t = time/row wanted
	 # returns figure of age structure in 1 year bins
	 # Global: ages, plot.relage, indicies defined above
	 ############################
	 x <- as.matrix(sol[t, c(2:length(colnames(sol)))])
	 xcounts <- NA; xbins <- NA
	 if (length(x) != length(ages)*6){
	 	print(
	 		paste("The age structure should include ",  N, 
	 			" ages, for 6 disease classes, giving ", 6*N,  
	 			" columns!", sep = ""))
	 }
	 # Sum values accross each disease state
	 for (i in 1:length(ages)){
	 	xcounts[i] <- (x[i] + x[min(it.index) - 1 + i] + x[min(ib.index) - 1 + i] + 
	 		x[min(ic.index) - 1 + i] + x[min(r.index) - 1 + i] + 
	 		x[min(rc.index) - 1 + i]) / sum(x)
	 }
	 binsize <- 1 / agestep # number of agebins per 1-year age group
	 xbins <- sum(xcounts[1:binsize])
	 for (i in 1:19) { 
	 	minbin <- i * binsize + 1 # 1, 11, 21, 31, ...
	 	maxbin <- (i+1)*binsize # 10, 20, 30...
	 	xbins <- c(xbins, sum(xcounts[minbin:maxbin]))
	 }
	 
	 # Get all values in each, 1-year age bin
	 d<-matrix(c(plot.relage, xbins), nrow = 2, byrow = TRUE, 
	 	dimnames = list(c("Observed", "Predicted"), 
	 	c(seq(1:20))))
	barplot(d, beside = TRUE, col = c("light gray", "dark gray"), 
		ylab = "Frequency (%)")
	box(	)
	legend("topleft", legend = c("Observed", "Predicted"), fill = c("light gray", "dark gray"))
}


plot_raw_numbers = function(sol){	
	plot(sol$time, apply(sol[s.index+1], 1, sum), col= "black",
		type= 'l', ylim = c(0, 800), ylab = "Number of animals", 
		xlab = "Time (in years)")
	lines(sol$time, apply(sol[it.index+1], 1, sum), col= "red")
	lines(sol$time, apply(sol[ib.index+1], 1, sum), col= "blue")
	lines(sol$time, apply(sol[ic.index+1], 1, sum), col= "green")
	lines(sol$time, apply(sol[r.index+1], 1, sum), col = "orange")
	lines(sol$time, apply(sol[rc.index+1], 1, sum), col = "pink")
	legend("topright", legend = c("S", "It", "Ib", "Ic", "R", "Rc"),
		col = c("black", "red", "blue", "green", "orange", "pink"), 
		bty="n", lty = 1)
}

get_prevalence = function(sol){
	S <- sum(sol[length(sol[,1]), s.index+1])
	It <- sum(sol[length(sol[,1]), it.index +1])
	Ib <- sum(sol[length(sol[,1]) , ib.index +1])
	Ic <- sum(sol[length(sol[,1]) , ic.index +1])
	R <- sum(sol[length(sol[,1]) , r.index +1])
	Rc <- sum(sol[length(sol[,1]) , rc.index +1])
	N <- sum(sol[length(sol[,1]), 2:length(colnames(sol))])
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


###########################################
# Figures specific for fitting
###########################################
get_structured_prevalence = function(sol){
	S <- sol[length(sol[,1]) , s.index+1]
	It <- sol[length(sol[,1]) , it.index +1]
	Ib <- sol[length(sol[,1]) , ib.index +1]
	Ic <- sol[length(sol[,1]) , ic.index +1]
	R <- sol[length(sol[,1]) , r.index +1]
	Rc <- sol[length(sol[,1]) , rc.index +1]
	
	# bin disease stats by age
	binsize <- 1/agestep
	Sb <- sum(S[1:binsize]); Itb <- sum(It[1:binsize])
	Ibb <- sum(Ib[1:binsize]); Icb <- sum(Ic[1:binsize])
	Rb <- sum(R[1:binsize]); Rcb <- sum(Rc[1:binsize])
	for (i in 1:19) {
		minbin <- i * binsize + 1 # 1, 11, 21, 31, ...
	 	maxbin <- (i+1) * binsize # 10, 20, 30...
	 	Sb <- c(Sb, sum(S[minbin:maxbin]))
	 	Itb <- c(Itb, sum(It[minbin:maxbin]))
	 	Ibb <- c(Ibb, sum(Ib[minbin:maxbin]))
	 	Icb <- c(Icb, sum(Ic[minbin:maxbin]))
	 	Rb <- c(Rb, sum(R[minbin:maxbin]))
	 	Rcb <- c(Rcb, sum(Rc[minbin:maxbin]))
	}
	
	S <- sum(Sb * data_agestructure)
	It <- sum(Itb * data_agestructure)
	Ib <- sum(Ibb * data_agestructure)
	Ic <- sum(Icb * data_agestructure)
	R <- sum(Rb * data_agestructure)
	Rc <- sum(Rcb * data_agestructure)
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







# use params or params_recov based on assumptions
#get_starting_eqbruc = function(params){
#	S0 = 400* stable_age; It0 = 0 * stable_age; Ib0 = 20* stable_age; 
#	Ic0 = 0* stable_age; R0 = 30 * stable_age; Rc0 = 0 * stable_age
#	x0 = c(S0, It0, Ib0, Ic0, R0, Rc0)
#	times <- seq(0, 1000, 1)
#	sol <- as.data.frame(ode(x0, times, rhs_age_matrix, params))
#	out <- unname(unlist( sol[1000, c(2:121)] ))
#	return(out)
#}
#test <- get_starting_eqbruc(params.test_log); sum(test[21:120])/sum(test)


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






###########################################
# Not Updated
###########################################

plot_dz_agestructure = function(x, dz){
	# Disease options =
	# "bruc" = brucellosis only
	# "tb" = bTB only
	# "co" = all subtypes
	# x = a row in sol
	if (length(x) != 120){
	 	print("The age structure should include 20 ages, 
	 	for 6 disease classes, giving 120 columns")
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
	
	par(mfrow = c(1,2))
	plot(y = prevB, x= seq(1, 20, 1), type = "b", col = "dark blue", ylim = c(0, 0.8),
		ylab = "Prevalence", xlab = "Age", pch = 19, 
		main = paste("Overall prevalences, Br =", round(overall_prevB, 3), 
		" TB = ", round(overall_prevT, 3) ))
	points(y = prevT, x = seq(1,20,1), type = "b", col = "dark red", pch = 19)
	legend("bottomright", bty = "n", legend = c("Bruc", "TB"), 
		pch = c(19, 19), col = c("dark blue", "dark red"))
	plot(y = prevBinS, x = seq(1,20,1), type = "b", col = "dark blue", pch = 19, 
		ylab = "Brucellosis prevalence", xlab = "Age", ylim = c(0, 0.8), 
		main = paste("Br|S =", round(overall_prevBinS, 3), 
		" Br|Co = ", round(overall_prevBinT, 3) )
		)
	text(x = 10, y = 0.55, labels = paste("Final N = ", round(overallN, 2)))
	points(y = prevBinT, x = seq(1,20,1), type = "b", col = "dark red", pch = 19)
}

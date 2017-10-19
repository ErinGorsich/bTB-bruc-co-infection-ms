# ODE 1D play

Aphid <- function (t, APHIDS, parameters) {
	deltax <- c(0.5, rep(1, numboxes - 1), 0.5)
	Flux <- - D * diff(c(0, APHIDS, 0)) / deltax
	dAPHIDS <- -diff(Flux) / delx + APHIDS * r
	list(dAPHIDS)
}

# ODE 1D play
D <- 0.3 # diffusion rate in m2/day
r <- 0.01 # /day
delx <- 1 # thickness of boxes
numboxes <- 60
Distance <- seq(0.5, by = delx,  length.out = numboxes)

# initial conditions
APHIDS <- rep(0, times = numboxes)
APHIDS[30:31] <- 1
state <- c(APHIDS = APHIDS)
out <- ode.1D(state, times, Aphid, parms = 0, nspec = 1, names = "Aphid")



# 
-diff( - D * diff(c(0, APHIDS, 0)) / deltax ) / delx

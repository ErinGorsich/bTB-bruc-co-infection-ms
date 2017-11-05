#######################################################
#######################################################
# Figures for Supplement
# bTB-brucellosis co-infection manuscript
# Erin Gorsich
#######################################################
#######################################################
setwd("~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/draft2/post-labmeeting/post-labmeeting2/Vanessa&Rampalcomments/draft_with_pde/figures")
source('~/GitHub/bTB-bruc-co-infection-ms/multiplot.R', chdir = TRUE)
library(ggplot2)
library(tidyr)
library('grid')
library('gridExtra') # specifies layout
library(survival)
library(lattice)
library(RColorBrewer)
library(maptools)   # for geospatial services; also loads foreign and sp
library(rgdal)      # for map projection work; also loads sp
library(RgoogleMaps)
library(dismo)		# basemaps
library(deSolve)


#######################################################
#######################################################
# Figure A1- GIS map of park and sample locations
#######################################################
#######################################################

makeTransparent<-function(someColor, alpha=100)
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
    blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

# Load Park boundary & rivers
###########################################
# To plot any shp file, you read them and then project them, this is the common projection
# plot then the name plots them, you can use any of the normal commands to change things (fill, col to specify border and inside colors, lwd to specify line type ect.)
boundary<-readShapePoly("~/Documents/GIS/KNP_data/kruger_boundary.shp")
proj4string(boundary) <- "+proj=longlat +datum=WGS84"
rivers<-readShapeLines("~/Documents/GIS/KNP_data/rivers_main.shp")
proj4string(rivers) <- "+proj=longlat +datum=WGS84"
sections <- readShapePoly("~/Documents/GIS/KNP_data/sections.shp")
proj4string(sections) <- "+proj=longlat +datum=WGS84"
waterways <- readShapeLines("~/Documents/GIS/KNP_data/s-africa-lesetho-waterways-shape/waterways.shp")
proj4string(sections) <- "+proj=longlat +datum=WGS84"
countries <- readShapePoly("~/Documents/GIS/KNP_data/AfricanCountries/AfricanCountires.shp")
proj4string(sections) <- "+proj=longlat +datum=WGS84"

# Load capture locations
d <- read.csv("~/Documents/phd research/big datasheets_updated to March2013/Demogall_Jan2013dbupdate_annaages2.csv")
df <- d[, colnames(d) %in% c('Animal.ID','capture.ID', 
	'capture.termination', 'x.cood.S.', 'y.cood.E.', 'Herd.captured.with')]
colnames(df) <- c("ID", "capID", "Jocapt", "x", "y", "herd")
df$x <- as.numeric(as.character(df$x))
df$y <- as.numeric(as.character(df$y))
df <- df[!is.na(df$x),]
df <- df[!is.na(df$y),]
#df <- df[df$Jocapt == 2,]
df <- df[df$Jocapt == 3,]
coordinates(df) <- c("y", "x")
proj4string(df) <- CRS("+proj=longlat +datum=WGS84")

# outside/overview figure
plot(countries[countries@data$COUNTRY %in% c("South Africa", 
	"Mozambique", "Botswana", "Namibia", "Zimbabwe"),])
plot(countries[countries@data$COUNTRY %in% c("South Africa"),],
	 col = "gray", add = TRUE)
plot(boundary, axes=TRUE, border="black", col = "olivedrab", add = TRUE, 
	lwd = 1.2, main = "") 
text(25,-28,"South Africa")

# Load base image; reproject 
# Use RgoogleMaps AND the dismo package's basemap; specify Lat/Long range
x = c(3400000, 3588000); y = c(-2549398, -2940755)
xy = cbind(x, y)
base.map <- gmap(xy, type = "satellite", scale = 2) # type=hybrid, terrain, satellite, roadmap
reprojected.boundary <- spTransform(boundary, base.map@crs)
reprojected.sections <- spTransform(sections, base.map@crs)
reprojected.df <- spTransform(df, base.map@crs)

# inset google maps figure
plot(base.map, axes = FALSE)
plot(reprojected.boundary, add = TRUE, border = "black", 
	col = makeTransparent("gray", 110), lwd = 1.2)
plot(reprojected.sections[reprojected.sections@data$SECTION  
	%in% c("Lower Sabie"),],  border = "yellow", add = TRUE, lwd = 1.2)
plot(reprojected.sections[reprojected.sections@data$SECTION  
	%in% c("Crocodile Bridge"),],  border = "gold", add = TRUE, lwd = 1.2)
points(reprojected.df, pch = 21, col = "red", add = TRUE, cex = 0.2)



########################################################################
########################################################################
# Figure B1: Fecundity Figure
########################################################################
########################################################################
newdf<-data.frame(
	Calf=c(11/16, 7/24, 6/16, 4/7, 5/35, 3/17, 3/14, 3/14), 
	Agecategory = c(rep("Adult (age > 4)", 4), rep("Sub-adult (age = 4)", 4)),
	Infection = c("uninfected", "brucellosis", "BTB", "co-infected", "uninfected", "brucellosis", "BTB", "co-infected"),
	N = c(25, 33, 23, 14, 26, 8, 7, 7))
newdf$color <- as.factor(seq(1,8))
newdf$se<- sqrt(newdf$Calf * (1 - newdf$Calf) / newdf$N)
newdf$order <- c(1,3,2,4,1,3,2,4)
newdf$Infection <- as.factor(newdf$Infection)
newdf$Infection <- factor(newdf$Infection, levels = newdf$Infection[order(unique(newdf$order))])
        
# With Bree's color scheme & juvenile ones on! for adults only: newdf2 <- newdf[c(1:4),]
p9 <- ggplot(newdf, aes(x=Infection, y=Calf, colour = as.character(color), shape = Agecategory, alpha = as.character(color))) + 
	geom_point(position = position_dodge(0.3), size=4) + 
	geom_errorbar(aes(ymin= newdf$Calf-newdf$se, ymax=newdf$Calf+newdf$se), 
  		position = position_dodge(0.3), width=0.2) + 
	scale_colour_manual(values=c("goldenrod1", "slateblue3", "chartreuse4",
		"tomato3", "goldenrod1", "slateblue3", "chartreuse4","tomato3"), 
		guide = FALSE) + 
	xlab("") +
	scale_alpha_manual(values = c(1,1,1,1,0.5, 0.5, 0.5, 0.5), guide = FALSE)       
p10 <- p9 +  theme_bw() + # removes ugly gray.
  ylab("Fecundity") +  # Proportion of buffalo observed with a calf
  scale_y_continuous(limits=c(0,0.8))  + 
  theme(axis.line.x = element_line(colour= "black"),
  		axis.line.y = element_line(colour= "black"),
  		axis.title.x = element_text(size=16, vjust=- 0.15),
        axis.title.y = element_text(size=16, vjust= 0.8, margin = margin(r = 8)),
        axis.text.x = element_text(size=14, vjust=- 0.05, margin = margin(t = 7)),
        axis.text.y = element_text(size=14, margin = margin(r = 4)),
        panel.border = element_blank(),
        # legend information
        legend.position=c(0.85, 0.95),  
        legend.background= element_rect(fill="white", colour="white"),
        legend.key= element_blank(),
        legend.title= element_blank(),
        legend.text = element_text(size=13)) 

png("Figure_B1_fecundity.png", width = 600, height = 450, units = "px")
p10
dev.off()


#######################################################
#######################################################
# Figure C1- Beverton-Holt and stable age distribution
#######################################################
#######################################################
#######################################################
# Left Pannel, Beverton-Holt
# Right Pannel, stable age distribution with no disease

# parameters to vary and recruitment function
########################################################
thetaL = seq(0.1, 0.9, by = 0.1)
thetaH = seq(1.1, 1.9, by = 0.1)
N = seq(1, 800, 1)
f_N = function(N, theta){
	0.5 / (1 + ((N/443)^theta))
}

# set up features of aging, functions, parameters for pde
########################################################
source('~/GitHub/bTB-bruc-co-infection-ms/pde/rhs.R', chdir = TRUE)
source('~/GitHub/bTB-bruc-co-infection-ms/pde/fixed_parameters.R', chdir = TRUE)

agemax <- 20
agestep <- 0.1
N <- agemax / agestep
ages <- seq(1, agemax + 1, by = agestep)[-(N)]
N == length(ages)
s.index <- 1:N
it.index <- seq(N+1, 2*N)
ib.index <- seq((2*N+1), 3*N)
ic.index <- seq((3*N+1), 4*N)
r.index <- seq((4*N+1), 5*N)
rc.index <- seq((5*N+1), 6*N)

f.params <- gen_fixed_params(agemax, agestep, p = p, recovery = FALSE)
f.params.recov <- gen_fixed_params(agemax, agestep, p = p, recovery = TRUE)

# Starting agestructure (Jolles 2007; Caron et al. 2001)
juv <- rep(0.137 / length(ages[ages < 2]), length(ages[ages < 2]))
sa <- rep(0.368 / length(ages[ages >= 2 & ages < 6]), length(ages[ages >= 2 & ages < 6]))
a <- rep(0.185 / length(ages[ages >= 6 & ages < 9]), length(ages[ages >= 6 & ages < 9]))
ma <- rep(0.235 / length(ages[ages >= 9 & ages < 14]), length(ages[ages >= 9 & ages < 14]))
sen <- rep(0.075 / length(ages[ages >= 14 ]), length(ages[ages >= 14]))

relage <- c(juv, sa, a, ma, sen); length(relage) == N									
plot.relage <- c(0.137, rep(0.368/4, 4), rep(0.185/3, 3), rep(0.235/5, 5), rep(0.075/7, 7)) 	


S0 <- 400 * relage; It0 <- 0 * relage; Ib0 <- 0 * relage
Ic0 <- 0 * relage; R0 <- 0 * relage; Rc0 <- 0 * relage
times <- seq(1, 300, 1)
x0 <- c(S0, It0, Ib0, Ic0, R0, Rc0)

params <- c(f.params, list(gamma = 1/2, betaB = 0.5764065, 
	betaT =1.3305462/1000, rhoT = 1, rhoB = 2.1))

sol <- as.data.frame(ode.1D(x0, times, rhs, params, 
		nspec = 6, dimens = N, method = "ode45")) 



# Make plot
########################################################
png("FigureC1_BH_and_agestructure.png", width = 800, height = 400, units = "px")
par(mfrow = c(1,2), mar = c(5,6,2,2))
#par(mar = c(5,6,4,2))
plot(x = N, y = f_N(N, theta = 4), type = "l", ylab = "Per captia birth rate, R(a, N)", 
	xlab = "Population size, N", bty = "n", las = 1, cex.axis = 1.2, cex.lab = 1.4, 
	ylim = c(0, 0.5))
abline(v = 443, col = "dark red")
lines(x = N, y = f_N(N, theta = 2), type = "l", col = "black", lty = 3)
lines(x = N, y = f_N(N, theta = 6), type = "l", col = "black", lty = 5)
legend("topright", bty = "n", 
	legend = c(expression(phi == 2), expression(phi == 4), expression(phi == 6)),
	lty = c(3, 1, 5) )

# Plot final age structure at disease free equilibrium
x = seq(1:200)
plot(y = sol[300, 2:201] / (sum(sol[300, 2:201])* agestep ), x = seq(1:200) / 10, 
	pch = 19, ylab = "Density", xlab = "Age (years)", bty = "n", las = 1, 
	cex.axis = 1.2, cex.lab = 1.4 )
dev.off()

#######################################################
#######################################################
# Figure S5- Endemic prevalence figures for each form of density dependence
#######################################################
#######################################################
d <- readRDS("~/GitHub/bTB-bruc-co-infection-ms/pde/sensitivity_densitydependence_results.rds")
d <- d[!is.na(d$rEE_bTB_co), ] # one NA

# test figure with violin plots
n = length(d$EE_bTB_single)
df <- data.frame(E = c(d$EE_bTB_single, d$EE_bTB_co, 
	d$EE_brucellosis_single, d$EE_brucellosis_co, 
	d$rEE_bTB_single, d$rEE_bTB_co, 
	d$rEE_brucellosis_single, d$rEE_brucellosis_co, 
	d$lEE_bTB_single, d$lEE_bTB_co, 
	d$lEE_brucellosis_single, d$lEE_brucellosis_co),
	meanE = c(rep(mean(d$EE_bTB_single),n), 
		rep(mean(d$EE_bTB_co),n),
		rep(mean(d$EE_brucellosis_single),n), 
		rep(mean(d$EE_brucellosis_co),n),
		rep(mean(d$rEE_bTB_single),n), 
		rep(mean(d$rEE_bTB_co),n),
		rep(mean(d$rEE_brucellosis_single),n), 
		rep(mean(d$rEE_brucellosis_co),n), 
		rep(mean(d$lEE_bTB_single),n), 
		rep(mean(d$lEE_bTB_co),n),
		rep(mean(d$lEE_brucellosis_single),n), 
		rep(mean(d$lEE_brucellosis_co),n)),
	sdE = c(rep(sd(d$EE_bTB_single),n), 
		rep(sd(d$EE_bTB_co),n),
		rep(sd(d$EE_brucellosis_single),n), 
		rep(sd(d$EE_brucellosis_co),n),
		rep(sd(d$rEE_bTB_single),n), 
		rep(sd(d$rEE_bTB_co),n),
		rep(sd(d$rEE_brucellosis_single),n), 
		rep(sd(d$rEE_brucellosis_co),n),
		rep(sd(d$lEE_bTB_single),n), 
		rep(sd(d$lEE_bTB_co),n),
		rep(sd(d$lEE_brucellosis_single),n), 
		rep(sd(d$lEE_brucellosis_co),n)),
	infection = c(rep("BTB", 2*n), rep("brucellosis", 2*n),
		rep("BTB", 2*n), rep("brucellosis", 2*n), 
		rep("BTB", 2*n), rep("brucellosis", 2*n)), 
	model = rep(c("Beverton & Holt", "Ricker", "Logistic" ), each = 4*n),
	singleco = rep(c("single", "co-infection", "single","co-infection",
		"single", "co-infection", "single","co-infection",
		"single", "co-infection", "single","co-infection"), each = n),
	X = rep(c("bTB-single-BH", "bTB-co-infection-BH", "bruc-single-BH", "bruc-co-infection-BH", 
		"bTB-single-R", "bTB-co-infection-R", "bruc-single-R", "bruc-co-infection-R", 
		"bTB-single-L", "bTB-co-infection-L", "bruc-single-L", "bruc-co-infection-L"), each = n))
#	Xindex = rep(c(1.1, 1.4, 3.1, 3.4), each = n) )
df$X <- factor(df$X, levels = c("bTB-single-BH", "bTB-co-infection-BH", "bruc-single-BH", "bruc-co-infection-BH", 
		"bTB-single-R", "bTB-co-infection-R", "bruc-single-R", "bruc-co-infection-R", 
		"bTB-single-L", "bTB-co-infection-L", "bruc-single-L", "bruc-co-infection-L"))
df$infection <- relevel(as.factor(df$infection), "BTB")
df$singleco <- relevel(as.factor(df$singleco), "single")
df <- df[!is.na(df$E),]
df <- df[df$model != "logistic",]
df <- df[!is.na(df$sdE),]
df$infectionmodel <- paste(df$infection, df$model, sep = "_")

brucellosis <- df[df$infection  == "brucellosis",]
tb <- df[df$infection == "BTB",]

# Dot and error plots for Endemic Prevalence
pB <- ggplot(brucellosis, aes(x = model, y = E, shape = singleco, colour = singleco)) +
	geom_point(aes(x = model, y = meanE), size = 3,
		position= position_dodge(width = 0.5)) +
	geom_errorbar(aes(x = model, ymin = meanE - sdE, ymax = meanE + sdE),
		width = 0, position= position_dodge(width = 0.5)) +
	ylim(0, 0.8) +
	xlab("") +
	ylab("Brucellosis prevalence") +
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
		legend.position= c(0.85, 0.9),  
		legend.text = element_text(size = 14),
        legend.background= element_rect(fill="white", colour="white"),
        legend.key= element_blank(),
        legend.title= element_blank())

pT <- ggplot(tb, aes(x = model, y = E, shape = singleco, colour = singleco)) +
	geom_point(aes(x = model, y = meanE), size = 3,
		position= position_dodge(width = 0.5)) +
	geom_errorbar(aes(x = model, ymin = meanE - sdE, ymax = meanE + sdE),
		width = 0, position= position_dodge(width = 0.5)) +
	ylim(0, 0.8) +
	xlab("") +
	ylab("BTB prevalence") +
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


#######################################################
#######################################################
# Figure C - 6 in supplement- sensitivity
#######################################################
#######################################################
#######################################################
df <- read.csv("~/GitHub/bTB-bruc-co-infection-ms/pde/Ro_sensitivity.csv")
df$infection <- as.character(df$infection)
df$infection[df$infection == "TB"] <- "BTB"
df$infection <- as.factor(df$infection)
df$infection <- relevel(df$infection, "BTB")
df$order <- rep(c(5, 4, 3, 6, 7, 8, 2, 1, 10, 9, 12, 11, 13), 2)  # order by table 1
df$order <- as.factor(df$order)

# Ro
pB <- ggplot(df, aes(x = order, y = Ro, shape = infection, colour = infection)) +
	geom_point(aes(x = order, y = Ro), size = 3,
		position= position_dodge(width = 0.5)) +
	geom_errorbar(aes(x = order, ymin = cilow, ymax = cihigh),
		width = 0, position= position_dodge(width = 0.5)) +
	geom_hline(yintercept = 0, color = 'red', size = 0.5) +
	ylim(-1, 1) +
	xlab("") +
	ylab("PRCC, Ro") +
	theme_bw() +
	scale_x_discrete("", labels = c("K",
		expression(theta), expression(beta[T]), expression(beta[B]),  
		expression(gamma), expression(epsilon), 
		expression(beta[B]/beta[B]), expression(beta[T]/beta[T]),
		expression(paste(mu[T], "/", mu[S], sep = "")),
		expression(paste(mu[B], "/", mu[S], sep = "")), 
		expression(paste(mu[R], "/", mu[S], sep = "")), 
		expression(paste(mu[C], "/", mu[S], sep = "")), 
		expression(paste(mu[RC], "/", mu[S], sep = "")))) +
	scale_colour_manual(values = c("black","dark gray")) +
	scale_shape_manual(values = c(19, 17)) +
	theme(axis.line.x = element_line(colour= "black"),
  		axis.line.y = element_line(colour= "black"),
  		axis.title.x = element_text(size=16, vjust=-0.15),
        axis.title.y = element_text(size=16, vjust= 0.8, margin = margin(r = 8)),
        axis.text.x = element_text(size=16, vjust=-0.05, margin = margin(t = 7)),
        axis.text.y = element_text(size=14, margin = margin(r = 4)),
        panel.border = element_blank(), 
		legend.position= c(0.85, 0.9),  
		legend.text = element_text(size = 14),
        legend.background= element_rect(fill="white", colour="white"),
        legend.key= element_blank(),
        legend.title= element_blank())


df <- read.csv("~/GitHub/bTB-bruc-co-infection-ms/pde/EE_sensitivity.csv")
df$infection <- as.character(df$infection)
df$infection[df$infection == "TB"] <- "BTB"
df$infection <- as.factor(df$infection)
df$infection <- relevel(df$infection, "BTB")
df$order <- rep(c(5, 4, 3, 6, 7, 8, 2, 1, 10, 9, 12, 11, 13), 2)  # order by table 1
df$order <- as.factor(df$order)

# EE
pB <- ggplot(df, aes(x = order, y = EE, shape = infection, colour = infection)) +
	geom_point(aes(x = order, y = EE), size = 3,
		position= position_dodge(width = 0.5)) +
	geom_errorbar(aes(x = order, ymin = cilow, ymax = cihigh),
		width = 0, position= position_dodge(width = 0.5)) +
	geom_hline(yintercept = 0, color = 'red', size = 0.5) +
	ylim(-1, 1) +
	xlab("") +
	ylab("PRCC, Endemic prevalence") +
	theme_bw() +
	scale_x_discrete("", labels = c("K",
		expression(theta), expression(beta[T]), expression(beta[B]),  
		expression(gamma), expression(epsilon), 
		expression(beta[B]/beta[B]), expression(beta[T]/beta[T]),
		expression(paste(mu[T], "/", mu[S], sep = "")),
		expression(paste(mu[B], "/", mu[S], sep = "")), 
		expression(paste(mu[R], "/", mu[S], sep = "")), 
		expression(paste(mu[C], "/", mu[S], sep = "")), 
		expression(paste(mu[RC], "/", mu[S], sep = "")))) +
	scale_colour_manual(values = c("black","dark gray")) +
	scale_shape_manual(values = c(19, 17)) +
	theme(axis.line.x = element_line(colour= "black"),
  		axis.line.y = element_line(colour= "black"),
  		axis.title.x = element_text(size=16, vjust=-0.15),
        axis.title.y = element_text(size=16, vjust= 0.8, margin = margin(r = 8)),
        axis.text.x = element_text(size=16, vjust=-0.05, margin = margin(t = 7)),
        axis.text.y = element_text(size=14, margin = margin(r = 4)),
        panel.border = element_blank(), 
		legend.position= c(0.85, 0.9),  
		legend.text = element_text(size = 14),
        legend.background= element_rect(fill="white", colour="white"),
        legend.key= element_blank(),
        legend.title= element_blank())


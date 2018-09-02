#######################################################
setwd("~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/final_datasets_copied_from_phdfolder")
source('~/GitHub/bTB-bruc-co-infection-ms/multiplot.R', chdir = TRUE)
library('ggplot2')
library('grid')
library('gridExtra') # specifies layout
library('survival')
library('lattice')
library('RColorBrewer')

#######################################################
#######################################################
# Herd specific parameters (from fit_TB_bruc)
#######################################################
#######################################################
# Plot best fit, no recovery
params <- c(f.params, list(gamma = 1/2, betaB = 0.5764065, 
    betaT = 1.3305462/1000, rhoT = 1, rhoB = 2.1))
paramsLS <- c(f.params, list(gamma = 1/2, betaB = params.est[1], 
    betaT = params.est[2]/1000, rhoT = 1, rhoB = 1))
paramsCB <- c(f.params, list(gamma = 1/2, betaB = params.est[1], 
    betaT = params.est[2]/1000, rhoT = 1, rhoB = 1))

paramlist <- list(params, paramsLS, paramsCB)



#x0 <- xB
#x0[[min(it.index) + 5*binsize]] <- 5
#x0[[min(it.index) + 5*binsize + 1]] <- 5
#sol <- as.data.frame(ode.1D(x0, times, rhs, params, 
#                            nspec = 6, dimens = N, method = "ode45"))
#xTsol <- sol
#get_structured_prevalence(sol)
#plot_raw_numbers(sol)

#######################################################
#######################################################
# Figure 4- Levelplots, varying mort and rho for both pathogens
#######################################################
#######################################################
# individual level prevalence
epiT <- read.csv("~/GitHub/bTB-bruc-co-infection-ms/epiT.csv")
epiB <- read.csv("~/GitHub/bTB-bruc-co-infection-ms/epiB.csv")

epiTp <- read.csv("~/GitHub/bTB-bruc-co-infection-ms/epiTp.csv")
epiBp <- read.csv("~/GitHub/bTB-bruc-co-infection-ms/epiBp.csv")

epiT <- cbind(epiT, epiTp[ ,4:7])
epiB <- cbind(epiB, epiTp[ ,4:7])

#epiB <- epiB[epiB$mort < 10.01,]
#epiT <- epiT[epiT$mort < 10.01,]
cols <- brewer.pal(11, "RdBu")
cols2 <- colorRampPalette(brewer.pal(11, "RdBu"))

epiT$bTBplot <- epiT$bTBprev - 0.65855
epiB$bTBplot <- epiB$bTBprev - 0.65855
epiT$brucplot <- epiT$brucprev - 0.21072
epiB$brucplot <- epiB$brucprev - 0.21072

df <- data.frame(Difference = c(epiT$bTBplot, epiB$brucplot), 
                 infection = c(rep("bTB", length(epiT[,1])), 
                               rep("brucellosis", length(epiB[,1]))), 
                 rho = c(epiT$rhoT, epiB$rhoB), mort = c(epiT$mort, epiB$mort))
df <- df[df$mort < 14.2, ]
df <- df[df$infection == "brucellosis", ]

# interaction term regression
# df2 <- data.frame(rho = c(4.3, 0.4, 1.2), mort = c(8.5, 8.5, 8.5), 
#     infection = c("brucellosis", "brucellosis", "bTB"), 
#     index = c(1, 2, 3))  # 5.8 in old
# df2$rho_seup <- c(7.39, 1.23, 1.85)		# bTB 1.29 (CI 0.645 to 2.606)
# df2$rho_sedown <- c(2.528, 0.124, 0.91)     # bruc trans (CI = 0.91 to 4.98)
# df2$mort_selow <- c(5.2, 5.2, 5.2) # was 5.8
# df2$mort_seup <- c(14.06, 14.06, 14.06)

# fit to each herd separately
df2 <- data.frame(rho = c(3.8, 1, 2.1), mort = c(8.5, 8.5, 8.5), 
                  infection = c("herd", "herd", "overall"), 
                  colour = c("gray", "gray", "black"))  # 5.8 in old
df2$rho_seup <- c(6.67, 1, 3.28)		# bTB 1.29 (CI 0.645 to 2.606)
df2$rho_sedown <- c(2.17, 1, 1.38)     # bruc trans (CI = 0.91 to 4.98)
df2$mort_selow <- c(5.2, 5.2, 5.2) # was 5.8
df2$mort_seup <- c(14.06, 14.06, 14.06)
df2.5 <- df2[df2$infection == "overall", ]
df2 <- df2[df2$infection == "herd",]

# For publication final
p <- ggplot(data = df, aes(x = mort, y = rho))
p2 <- p + theme_bw() + 
    xlab("Proportional increase in mortality with co-infection") + 
    ylab("Proportional increase in transmission") +
    theme(axis.ticks = element_line(size = 0.1),
        panel.grid.major = element_blank(), 
        strip.text.x = element_text(size = 5.7),
        plot.title = element_text(size = 5.7),  # 14 in pdf format
        axis.title.x = element_text(size = 5.7, margin = margin(t = 2)), 
        axis.title.y = element_text(size = 5.7, margin = margin(r = 2)),
        axis.text.x = element_text(size=5.6, margin = margin(t = 2)),
        axis.text.y = element_text(size = 5.6, margin = margin(r = 2)), # 12 
        legend.text = element_text(size = 5.6),  # 12 
        legend.title = element_text(size = 5.6),
        legend.text.align = 1, 
        legend.key.size = unit(0.2, "cm"), 
        panel.grid.minor = element_blank()) + 
    geom_raster(aes(fill = Difference), interpolate = TRUE, show.legend = TRUE) +  # was geom_tile
    scale_fill_distiller(palette = "RdYlBu", direction = -1, 
        limits = c(-0.7, 0.7)) + 
    geom_point(data = df2, size = 0.6, pch = 19, colour = "gray") + 
    geom_errorbar(data = df2, aes(ymin = rho_sedown, ymax = rho_seup), 
        colour = "dark gray", size = 0.3, width = 0) + 
    geom_errorbarh(data = df2, aes(xmin = mort_selow, xmax = mort_seup), 
        colour = "dark gray", size = 0.3, height = 0) +     
    geom_point(data = df2.5, size = 0.6, pch = 19, colour = "black") + 
    geom_errorbar(data = df2.5, aes(ymin = rho_sedown, ymax = rho_seup), colour = "black", 
                  size = 0.3, width = 0) + 
    geom_errorbarh(data = df2.5, aes(xmin = mort_selow, xmax = mort_seup), colour = "black", 
                   size = 0.3, height = 0)
  

df$Difference2 <- df$Difference - 0.1  # want contours to span -5 to +5 not 0 to 10
df3 <- data.frame(Difference = c(0.10, -0.10, -0.3, 0.50, 0.30, 0.10, -0.10), 
    rho = c(7.8, 7.8, 3.75, 7.8, 7.8, 7.8, 1.2), 
    mort = c(1.6, 3.8, 13.25, 1.7, 4.3, 12.15, 13.3), 
    infection = c("bTB", "bTB", "bTB", "brucellosis", 
        "brucellosis", "brucellosis", "brucellosis"))
df3 <- df3[df3$infection == "brucellosis", ]
p3 <- p2  + 
    geom_contour(data = df, 
        aes(x = mort, y = rho, z = Difference2, weight = ..level..),
        binwidth = 0.2, color = "black", linetype = 3, size = 0.2)
p4 <- p3 + geom_text(data = df3, size = 2, 
    aes(z = NULL, label = Difference))


tiff("Figure_S9_herdcontext.tiff", width = 7.8, height = 5, 
     units = "cm", res = 300)
p4
dev.off()

# in Lower Sabie, see co-infection increasing brucellosis prevlaence - 24%
df[df$rho %in% c(3.76 , 3.84) & df$mort %in% c(5.7, 5.85),]
# in Crocodile Bridge, see co-infection decreasing brucellosis prevalence - 21%
df[df$rho ==0.4 & df$mort %in% c(5.7, 5.85),] 
df[df$rho %in% c(0.94, 1.04) & df$mort %in% c(5.7, 5.85),] # or no change







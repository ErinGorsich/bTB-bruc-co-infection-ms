#######################################################
#######################################################
# Figures
# bTB-brucellosis co-infection manuscript
# Erin Gorsich
#######################################################
#######################################################
setwd("~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/final_datasets_copied_from_phdfolder")
source('~/GitHub/bTB-bruc-co-infection-ms/multiplot.R', chdir = TRUE)
library(ggplot2)
library(tidyr)
library('grid')
library('gridExtra') # specifies layout
library(survival)
library(lattice)
library(RColorBrewer)
library(deSolve)

# read in data prepared in cross_sectional_dataprep, groomed for my bTB statuses
#data<-read.csv("cross_sectional_data_withdz_cleandisease_nofinal_Feb2016_capturetime_forsurv.csv")
data<-read.csv("cross_sectional_data_withdz_cleandisease_nofinal_Feb2016_capturetime.csv")

data_nofinal<-data[data$final_capture=="0",] 
d<-data.frame(btb=data_nofinal$tb , bruc=as.character(data_nofinal$bruc), age=data_nofinal$age_sel/12, 
              id=data_nofinal$id)
d<-d[d$age<14,]


#######################################################
#######################################################
# Figure 1- AGE HISTOGRAM FOR THE CONCEPTUAL FIGURE
#######################################################
#######################################################
# only inset plot
png("Age_histogram.png", width = 800, height = 800, units = "px")
ggplot(d, aes(x=age)) + 
	geom_histogram(data=d, fill= "gray", alpha=0.3, #steelblue, steelblue4
		colour = "dimgray", stat="bin", binwidth= 1, linetype = 1, size = 2) + 
	theme_bw()+ 
	xlab("Age (years)")+ ylab("Number of buffalo")+
	theme(panel.border = element_blank(), 
        panel.margin = element_blank(),
        panel.grid.major= element_blank(),
        panel.grid.minor= element_blank()) + 
	theme(axis.line.x = element_line(colour= "black", size = 2), 
		axis.line.y = element_line(colour= "black", size = 2)) + 
	theme(axis.title = element_blank(), axis.text = element_blank())
dev.off()


d <- data.frame(Prevalence = c(0.27, 0.34, 0.27, 0.34), 
	Infection = c("bTB", "brucellosis", "bTB", "brucellosis"), 
	Color = c("a", "b", "c", "d"), 
	Type = c("Model", "Model", "Data", "Data"))
d$Infection <- relevel(d$Infection, "bTB")
makeTransparent<-function(someColor, alpha=100)
{
  newColor<-col2rgb(someColor)
  apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2],
    blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}
blue2 <- makeTransparent("slateblue4", 130)
green2 <- makeTransparent("chartreuse4", 130)

png("Evaluation.png", width = 800, height = 600, units = "px")
ggplot(d, aes(x = Type, y = Prevalence, fill = Color)) + guides(fill = FALSE) +
	geom_bar(position = position_dodge(), stat = "identity",) + ylim(0, 0.5) + 
	theme_bw() +  
	scale_fill_manual(values=c(blue2, green2,  "slateblue4", "chartreuse4"))  + 
	scale_colour_manual(values=c("black"))  + 
	theme(panel.border = element_blank(), 
        panel.margin = element_blank(),
        panel.grid.major= element_blank(),
        panel.grid.minor= element_blank()) + 
	theme(axis.line.x = element_line(colour= "black", size = 2), 
		axis.line.y = element_line(colour= "black", size = 2)) + 
	theme(axis.title = element_blank(), axis.text = element_blank())
dev.off()


# Evaluation take 2: 
setwd("~/GitHub/bTB-bruc-co-infection-ms")
df<- readRDS("coinfection_pattern_results.rds")
df2$color <- c(1,2,0, 0, 3, 4, 0, 0)
f2 <- data.frame(Evaluation = c("Model", "Model", "Data", "Data"), 
	Infection = c("Single infection", "Co-infection", "Single infection", "Co-infection"), 
	BrucellosisPrevalence = c(df$prevBinS, df$prevBinT, 0.3035, 0.4524), 
	TBPrevalence = c(df$prevTinS, df$prevTinB, 0.227, 0.3585), 
	CILowB = c(NA, NA, 0.2654405, 0.3657982), 
	CIHighB = c(NA, NA, 0.3451408, 0.5405985), 
	CILowTB = c(NA, NA, 0.1781895, 0.2857143), 
	CIHighTB = c(NA, NA, 0.2783505, 0.4313725))
df2 <- reshape(f2, varying = c("BrucellosisPrevalence", "TBPrevalence"), 
	v.names = "Prevalence", timevar = "Name", times = c("brucellosis", "bTB"), 
	direction = "long")
df2$Infection <- relevel(df2$Infection, "Single infection")
df2$CILow <- c(NA, NA, 0.26544, 0.365798, NA, NA, 0.1781895, 0.2857)
df2$CIHigh <- c(NA, NA, 0.3451408, 0.540598, NA, NA, 0.27835, 0.43137)
df3 <- df2[df2$Evaluation == "Model",]
df3$color <- c(rgb(97, 139, 14, maxColorValue= 255))
df3$Name <- as.factor(df3$Name)
df3$Name <- relevel(df3$Name, "bTB")
df3$color <- c("a", "b", "c", "d")

df3$data <- c(0.3035, 0.4524, 0.227, 0.3585)
df3$CILow = c(0.2654405, 0.3657982, 0.1781895, 0.2857143) 
df3$CIHigh = c(0.3451408, 0.5405985, 0.2783505, 0.4313725)
 
 
png("Evaluation.png", width = 800, height = 600, units = "px")
limits <- aes(ymax = df3$CIHigh, ymin= df3$CILow)
p1 <- ggplot(df3, aes(x = Name, y = Prevalence, fill = color)) + #, colour = Name
	geom_bar(position = position_dodge(), stat = "identity", colour = "black", cex = 2) +
	geom_point(data = df3, aes(y = data, x = Name), 
		position = position_dodge(width = 0.9), size = 6, colour = "black") +
	geom_errorbar(limits, position = position_dodge(width = 0.9), width= 0.2, cex = 2 ) +
	ylim(0, 0.8) + xlab("") + guides(fill = FALSE) +
	theme_bw() +  
	scale_fill_manual(values=c("#458b00", "#295300", "#6959cd", "#38306f")) +
	theme(panel.border = element_blank(), 
        panel.margin = element_blank(),
        panel.grid.major= element_blank(),
        panel.grid.minor= element_blank()) + 
	theme(axis.line.x = element_line(colour= "black", size = 2), 
		axis.line.y = element_line(colour= "black", size = 2)) + 
	theme(axis.title = element_blank(), axis.text = element_blank())
p1
dev.off()



#SlateBlue3, Slateblue4


#######################################################
#######################################################
# SUPPLEMENTAL FIGURE OR CUT- AGE PREVALENCE FIGURES (Old figure 1)
#######################################################
#######################################################
# make a dataframe with age prevalence in tbneg and tbpos, long format for ggplot
val<-NA
get_prev<-function(dat) {
  if(length(dat)==0){
    val<-0
  }
  if(length(dat>0)){
    val<-length(dat[dat=="1"])/length(dat)
  }
  return(val)
}

binsize=2
agebins=c(seq(0, max(d$age), binsize), max(d$age+1))
agebins<- c(agebins[c(1:6)], 16.5)
d$bruc2<-NA
newdf<-data.frame(agebin=c(agebins, agebins), Brucprev=NA, 
                  TB=c(rep("bTB-", length(agebins)),rep("bTB+", length(agebins)) ),
                  N=NA)  


for (i in 1:length(d$bruc)){
  ifelse(d$bruc[i]=="negative", d$bruc2[i]<-0, d$bruc2[i]<-1)
}
d$bruc<-as.numeric(d$bruc2)
for (i in 1:(length(agebins)-1)){
  neg<-d[d$btb==0,]
  pos<-d[d$btb==1,]
  d_neg<-d[d$age>=agebins[i] & d$age<agebins[i+1] & d$btb=="0",]
  d_pos<-d[d$age>=agebins[i] & d$age<agebins[i+1] & d$btb=="1",]
  newdf$Brucprev[i]<-get_prev(d_neg$bruc)
  newdf$Brucprev[i+length(agebins)]<-get_prev(d_pos$bruc)
  newdf$N[i]<- length(d_neg$bruc)
  newdf$N[i+length(agebins)]<- length(d_pos$bruc)
}
newdf<-newdf[newdf$agebin<16,]
newdf$se<-sqrt(newdf$Brucprev*(1-newdf$Brucprev)/newdf$N)
newdf$se[is.na(newdf$se)]<-0

# remove 0-2 bin because no TB+
newdf2<-newdf[newdf$N>1,]; newdf<-newdf2

# plot 1: 
p<- ggplot(newdf, aes(x=agebin, y=Brucprev, group=TB, colour=TB)) + 
  geom_line()+
  geom_point(size=3, shape=19) + # colour="darkred", fill="darkred" +
  scale_colour_manual(values=c("steelblue", "orangered4"))
p2<- p+ geom_errorbar(aes(ymin= newdf$Brucprev-newdf$se, ymax=newdf$Brucprev+newdf$se ), width=0.3) + 
  xlab("Age (years)") + ylab("Brucellosis prevalence") +
  theme_bw() + # removes ugly gray.
  scale_x_discrete(breaks=c("0", "2", "4", "6", "8", "10"), limits=seq(0, 10, 1), 
                   labels=c("0-2", "2-4", "4-6", "6-8", "8-10", "10+")) +
  scale_y_continuous(limits=c(0,1)) + 
  theme(axis.line.x = element_line(colour= "black"),
  		axis.line.y = element_line(colour= "black"),
  		axis.title.x = element_text(size=16, vjust=-0.15),
        axis.title.y = element_text(size=16, vjust= 0.8),
        axis.text.x = element_text(size=14, vjust=-0.05),
        axis.text.y = element_text(size=14),
        panel.border = element_blank(), 
        # legend information
        legend.position=c(0.82, 0.9),  
        legend.background= element_rect(fill="white", colour="white"),
        legend.key= element_blank(),
        legend.title= element_blank(),
        legend.text = element_text(size=15)   )

# inset plot
d$btb2<-as.factor(d$btb)
p3<-ggplot(d, aes(x=age)) + 
  geom_histogram(data=subset(d, btb2=="0"), fill= "steelblue", alpha=0.3, stat="bin", binwidth= 1) + 
  geom_histogram(data=subset(d, btb2=="1"), fill= "orangered", alpha=0.3, stat="bin", binwidth= 1) + 
  theme_bw()+ 
  xlab("Age (years)")+ ylab("Count")+
  theme(panel.border = element_blank(), 
        panel.margin = element_blank(),
        panel.grid.major= element_blank(),
        panel.grid.minor= element_blank()) + 
  theme(axis.line.x = element_line(colour= "black"), 
		axis.line.y = element_line(colour= "black")  )

vp<- viewport(width=0.4, height=0.4, x=0.34, y=0.75)
print(p2)
print(p3, vp=vp)


#############################################
# By heard: 
d<-data.frame(btb=data_nofinal$tb , bruc=as.character(data_nofinal$bruc), age=data_nofinal$age_sel/12, 
              id=data_nofinal$id, herd=data_nofinal$herdorig)
d<-d[d$age<14,]
binsize=2
agebins=c(seq(0, max(d$age), binsize), max(d$age+1))
agebins<- c(agebins[c(1:6)], 16.5)
d$bruc2<-NA

# LS
d<-d[d$herd=="LS",]
newdf<-data.frame(agebin=c(agebins, agebins), Brucprev=NA, 
                  TB=c(rep("bTB-", length(agebins)),rep("bTB+", length(agebins)) ),
                  N=NA)  
for (i in 1:length(d$bruc)){
  ifelse(d$bruc[i]=="negative", d$bruc2[i]<-0, d$bruc2[i]<-1)
}
d$bruc<-as.numeric(d$bruc2)
for (i in 1:(length(agebins)-1)){
  neg<-d[d$btb==0,]
  pos<-d[d$btb==1,]
  d_neg<-d[d$age>=agebins[i] & d$age<agebins[i+1] & d$btb=="0",]
  d_pos<-d[d$age>=agebins[i] & d$age<agebins[i+1] & d$btb=="1",]
  newdf$Brucprev[i]<-get_prev(d_neg$bruc)
  newdf$Brucprev[i+length(agebins)]<-get_prev(d_pos$bruc)
  newdf$N[i]<- length(d_neg$bruc)
  newdf$N[i+length(agebins)]<- length(d_pos$bruc)
}
newdf<-newdf[newdf$agebin<16,]
newdf$se<-sqrt(newdf$Brucprev*(1-newdf$Brucprev)/newdf$N)
newdf$se[is.na(newdf$se)]<-0

# remove 0-2 bin because no TB+
newdf2<-newdf[newdf$N>0,]; newdf<-newdf2

# plot 1: 
p4<- ggplot(newdf, aes(x=agebin, y=Brucprev, group=TB, colour=TB)) + 
  geom_line()+
  geom_point(size=3, shape=19) + # colour="darkred", fill="darkred" +
  scale_colour_manual(values=c("steelblue", "orangered4"))
p5<- p4+ 
  geom_errorbar(limits=aes(ymin= newdf$Brucprev-newdf$se, ymax=newdf$Brucprev+newdf$se ), width=0.3) + 
  xlab("Age (years)") + ylab("Brucellosis prevalence") +
  theme_bw() + # removes ugly gray.
  scale_x_discrete(breaks=c("0", "2", "4", "6", "8", "10"), limits=seq(0, 10, 1), 
                   labels=c("0-2", "2-4", "4-6", "6-8", "8-10", "10+")) +
  scale_y_continuous(limits=c(0,1)) + 
  theme(axis.line = element_line(colour= "black"),
  		axis.title.x = element_text(size=16, vjust=-0.15),
        axis.title.y = element_text(size=16, vjust= 0.8),
        axis.text.x = element_text(size=14, vjust=-0.05),
        axis.text.y = element_text(size=14),
        panel.border = element_blank(), 
        # legend information
        legend.position=c(0.82, 0.9),  
        legend.background= element_rect(fill="white", colour="white"),
        legend.key= element_blank(),
        legend.title= element_blank(),
        legend.text = element_text(size=15)   )

# inset plot
d$btb2<-as.factor(d$btb)
p6<-ggplot(d, aes(x=age)) + 
  geom_histogram(data=subset(d, btb2=="0"), fill= "steelblue", alpha=0.3, stat="bin", binwidth= 1) + 
  geom_histogram(data=subset(d, btb2=="1"), fill= "orangered", alpha=0.3, stat="bin", binwidth= 1) + 
  theme_bw()+ xlab("Age (years)")+ ylab("Count")+
  theme(panel.border = element_blank(), 
        panel.margin = element_blank(),
        panel.grid.major= element_blank(),
        panel.grid.minor= element_blank(),
        axis.line = element_line(colour= "black"))

vp<- viewport(width=0.4, height=0.4, x=0.34, y=0.75)
print(p5)
print(p6, vp=vp)

#######################################################
#######################################################
# Figure S1
#######################################################
#######################################################
###### SAME FIGURE BUT WITH Capture period effects
newdf<-data.frame(captbin=c(seq(0.9,7.9,1), seq(1.1,8.1,1)),
	prev=c(0.32, 0.284, 0.247, 0.2660, 0.34, 0.3696, 0.4300, 0.442, 
			0.118, 0.1300, 0.3, 0.2980, 0.300, 0.3695, 0.387, 0.4155), 
	dz = c(rep("Brucellosis", 8), rep("Tuberculosis", 8)),
	N = c(59, 80, 93, 94, 50, 92, 92, 77, 59, 80, 93, 94, 50, 92, 92, 77), 
	param = c(-0.7261, 0.10051 , 0.01135 , -0.231094, -1.9331 ,-0.87186 ,-1.1061 , -0.959), 
	sep = c(0.689, 0.484 , 0.326 , 0.354, 0.81 , 0.455 ,0.485 , 0.424))  
newdf$se<- sqrt(newdf$prev * (1 - newdf$prev) / newdf$N)

p3<- ggplot(newdf, aes(x=captbin, y=prev, group=dz, colour=dz)) + 
  geom_line()+
  geom_point(size=3, shape=19) + # colour="darkred", fill="darkred" +
  geom_errorbar(aes(ymin= newdf$prev-newdf$se, ymax=newdf$prev+newdf$se), width=0.2) + 
  scale_colour_manual(values=c("darkslategray", "darkseagreen3"))
p4<- p3 +
  theme_bw() + # removes ugly gray.
  xlab("Capture period (6 months)") + ylab("Sero-prevalence") +
  scale_x_discrete(breaks=c("1", "2", "3", "4", "5", "6", "7", "8"), limits=seq(1, 8, 1), 
                   labels=c("1", "2", "3", "4", "5", "6", "7", "8")) +
  scale_y_continuous(limits=c(0,0.6)) + 
  theme(axis.line.x = element_line(colour= "black"),
  		axis.line.y = element_line(colour= "black"),
  		axis.title.x = element_text(size=16, vjust=-0.15),
        axis.title.y = element_text(size=16, vjust= 0.8),
        axis.text.x = element_text(size=14, vjust=-0.05),
        axis.text.y = element_text(size=14),
        panel.border = element_blank(), 
        # legend information
        legend.position=c(0.82, 0.9),  
        legend.background= element_rect(fill="white", colour="white"),
        legend.key= element_blank(),
        legend.title= element_blank(),
        legend.text = element_text(size=15)   )


p4<- p3 +
  theme_bw() + # removes ugly gray.
  xlab("Capture period (6 months)") + ylab("Sero-prevalence") +
  scale_x_discrete(breaks=c("1", "2", "3", "4", "5", "6", "7", "8"), limits=seq(1, 8, 1), 
                   labels=c("1", "2", "3", "4", "5", "6", "7", "8")) +
  scale_y_continuous(limits=c(0,0.6)) + 
  theme(axis.line.x = element_line(colour= "black"),
  		axis.line.y = element_line(colour= "black"),
  		axis.title.x = element_text(size=16, vjust=-0.15),
        axis.title.y = element_text(size=16, vjust= 0.8),
        axis.text.x = element_text(size=14, vjust=-0.05),
        axis.text.y = element_text(size=14),
        panel.border = element_blank(), 
        # legend information
        legend.position=c(0.82, 0.9),  
        legend.background= element_rect(fill="white", colour="white"),
        legend.key= element_blank(),
        legend.title= element_blank(),
        legend.text = element_text(size=15)   )


newdf2<- newdf[newdf$dz=="Brucellosis",]
p5<- ggplot(newdf2, aes(x=captbin, y= param, colour= dz, fill=dz)) + 
  	xlab("Capture period (6 months)") + ylab("Effect size (TB by age interaction)") +
	geom_bar(stat="identity", fill="darkseagreen3") + 
	geom_errorbar(aes(ymin= newdf2$param - newdf2$sep, ymax = newdf2$param + newdf2$sep), 
	width=0.2) + 
	scale_colour_manual(values=c("darkslategray")) + 
	scale_x_discrete(breaks=c("1", "2", "3", "4", "5", "6", "7", "8"), limits=seq(1, 8, 1), 
                   labels=c("1", "2", "3", "4", "5", "6", "7", "8")) +
    theme_bw() +
    guides(fill=FALSE, colour= FALSE) +
  	theme(axis.line.x = element_line(colour= "black"),
  		axis.line.y = element_line(colour= "black"),
  		axis.title.x = element_text(size=16, vjust=-0.15),
        axis.title.y = element_text(size=16, vjust= 0.8),
        axis.text.x = element_text(size=14, vjust=-0.05),
        axis.text.y = element_text(size=14),
        panel.border = element_blank() )

#######################################################
#######################################################
# Figure 2- Survival plots
#######################################################
#######################################################
data2<-read.csv("~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/final_datasets_copied_from_phdfolder/survival/brucsurvival_TB3controls_longresidnomissing_noerrors_season2_final_fixed.csv")
bolus <- c("R22", "R24", "R27", "R35", "R45", "R45b", "R50", "R6", "Y20", "Y31c", "Y31d", "Y39", "Y43", "Y44")
#"W1" is ok.
data3<- data2[!(data2$animal %in% bolus),]
data3<- data3[data3$animal!= "O10",]

# V1
#final.mod<-coxph(Surv(start, stop, death.time)~brucella+TB_3+herd2+ age6+ cluster(animal), data=data3); summary(final.mod)
#full.mod<-test.mod

# V2
# 0-2, 3-5, 6+
#data3$age2.4 <- data3$age2.3
#data3$age2.4[data3$age_yr2 == 6] <- "matureadult"
#data3$age2.4[data3$age_yr2 == 3] <- "subadult"
#data3$age2.4 <- as.character(data3$age2.4)
#final.mod<-coxph(Surv(start, stop, death.time)~brucella*age2.4+TB_3+herd2+ cluster(animal), data=data3); summary(final.mod)


#######################################################
# Age and infection specific survival rates 
#######################################################
#  Age-Infection specific survival estimates for our population
df<-data.frame(start=data3$start2, stop=data3$stop2, death.time=data3$death.time, 
	age=data3$age6, herd=data3$herd, bruc= data3$brucella, tb = data3$TB_3, animal = data3$animal)
mort<-with(df, data.frame(
	age=c("adult", "adult", "adult", "adult", "adult", "adult", "adult", "adult", 
		"juvenile", "juvenile", "juvenile", "juvenile", "juvenile", "juvenile", "juvenile", "juvenile"),
	herd = c("LS", "CB", "LS", "CB", "LS", "CB", "LS", "CB"),
	tb = c(0, 0, 1, 1, 0, 0, 1, 1), 
	bruc = c(0, 0, 0, 0, 1, 1, 1, 1)))
plot_add.mod<-coxph(Surv(start, stop, death.time)~herd+age+ tb+ bruc + cluster(animal), data=df)
m<-survfit(plot_add.mod, newdata=mort); summary(m)
m<-survfit(plot_add.mod, newdata=mort[mort$herd=="LS" & mort$age=="adult",]); summary(m)
m<-survfit(plot_add.mod, newdata=mort[mort$herd=="CB" & mort$age=="adult",]); summary(m)


plot_add.mod<-coxph(Surv(start, stop, death.time)~age, data=df)
m<-survfit(plot_add.mod, newdata=mort)
summary(m)


# Compile estimates above into a df with other age specific estimate: 
survivaldf <- data.frame(age = seq(1, 15, 1), 
	SurvUn = c(NA, 0.86, rep(0.94, 10), NA, NA, NA),
	SurvTB = c(NA, 0.608, rep(0.832, 10), NA, NA, NA), 
	SurvBR =  c(NA, 0.58, rep(0.82, 10), NA, NA, NA),
	SurvCo = c(NA, 0.188, rep(0.652, 10), NA, NA, NA)
)

# age specific effects of brucellosis... not using
#survivaldf <- data.frame(age = seq(1, 15, 1), 
#	SurvUn = c(NA, 0.86, rep(0.94, 10), NA, NA, NA),
#	SurvTB = c(NA, 0.7186, rep(0.8794, 10), NA, NA, NA), 
#	SurvBR =  c(NA, 0.314, rep(0.940, 3), rep(0.730, 7), NA, NA, NA),
#	SurvCo = c(NA, 0.0326, rep(0.8794, 3), rep(0.6094, 7), NA, NA, NA)
#)
#original dataframe with outside sources compiled
#survivaldf <- data.frame(age = seq(1, 15, 1), 
	#Cross2009female = c(0.9, 0.9, rep(0.95, 5), rep(0.85, 5), rep(0.9, 3)),
	#Cross2009male = c(0.82, 0.82, rep(0.77, 2), rep(0.97, 3), rep(0.65, 5), rep(0.1, 3)), 
	#CrossGetz2006male = c(NA, rep(0.84, 7), rep(0.59, 7)),
	#CrossGetz2006female = c(NA, rep(0.95, 7), rep(0.86, 7)),
	#Jolles2005TBneg = c(rep(0.85, 4), rep(0.97, 11)),
	#Jolles2005TBpos = c(rep(0.85, 4), rep(0.86, 11)),
	#FunstonMills = c(0.87, rep(0.92, 14)), 
#	SurvUn = c(NA, 0.884, 0.884, rep(0.963, 9), NA, NA, NA),
#	SurvTB = c(NA, rep(0.689, 2), rep(0.892, 9), NA, NA, NA), 
#	SurvBR =  c(NA, rep(0.706, 2), rep(0.899, 9), NA, NA, NA),
#	SurvCo = c(NA, rep(0.349, 2), rep(0.724, 9), NA, NA, NA)
#)
#survlong <- gather(survivaldf, key = dataset, value = estimate, Cross2009female:SurvCo)
#survlong$colour <- "Estimates from African buffalo in southern Africa"
survlong <- gather(survivaldf, key = dataset, value = estimate, SurvUn:SurvCo)
survlong$colour <- NA
survlong$colour[survlong$dataset %in% c("SurvUn")] <- "uninfected"
survlong$colour[survlong$dataset %in% c("SurvTB")] <- "bTB"
survlong$colour[survlong$dataset %in% c("SurvBR")] <- "brucellosis"
survlong$colour[survlong$dataset %in% c("SurvCo")] <- "co-infected"
survlong$order <- seq(1, length(survlong[,1]))
survlong$colour <- as.factor(survlong$colour)
survlong$colour <- factor(survlong$colour, levels = survlong$colour[order(unique(survlong$order))])


p6 <- ggplot(survlong, aes(x = age, y = estimate, colour = colour, 
	size = colour, shape = colour, group = dataset)) + 
	xlab("Age (years)") + ylab("Annual survival") +
	geom_point(aes(size = colour)) +  #shape = sex
	geom_line(aes(linetype = colour)) + 
	scale_x_continuous(breaks=seq(1, 15, 2)) + 
	theme_bw() + 
	scale_colour_manual(values = c("goldenrod1", "slateblue3", "chartreuse4","tomato3")) +
	scale_size_manual(values = c(2, 2, 2, 2)) +
	scale_shape_manual(values = c(18, 17, 17, 19)) +
	scale_linetype_manual(values = c('longdash', 'dotdash', 'dashed', 'twodash')) +
	theme(panel.border = element_blank(), 
		axis.title.x = element_text(size=16, vjust=-0.15),
        axis.title.y = element_text(size=16, vjust= 0.8),
        axis.text.x = element_text(size=14, vjust=-0.05),
        axis.text.y = element_text(size=14)) + 
        theme(axis.line.x = element_line(colour= "black"),
  			axis.line.y = element_line(colour= "black"),   
			#legend information
        	legend.position=c(0.85, 0.2),  
        	legend.background= element_rect(fill="white", colour="white"),
        	legend.key= element_blank(),
        	legend.title= element_blank(),
        	legend.text = element_text(size=14))

p6red <- ggplot(survlong[survlong$age < 13.2,], aes(x = age, y = estimate, colour = colour, 
	size = colour, shape = colour, group = dataset)) + 
	xlab("Age (years)") + ylab("Annual survival") +
	geom_point(aes(size = colour)) +  #shape = sex
	geom_line(aes(linetype = colour)) + 
	scale_x_continuous(breaks=seq(1, 13, 2)) + 
	theme_bw() + 
	scale_colour_manual(values = c("goldenrod1", "slateblue3", "chartreuse4","tomato3")) +
	scale_size_manual(values = c(2, 2, 2, 2)) +
	scale_shape_manual(values = c(18, 17, 17, 19)) +
	scale_linetype_manual(values = c('longdash', 'dotdash', 'dashed', 'twodash')) +
	theme(panel.border = element_blank(), 
		axis.title.x = element_text(size=16, vjust=-0.15),
        axis.title.y = element_text(size=16, vjust= 0.8),
        axis.text.x = element_text(size=14, vjust=-0.05),
        axis.text.y = element_text(size=14)) + 
        theme(axis.line.x = element_line(colour= "black"),
  			axis.line.y = element_line(colour= "black"),   
			#legend information
        	legend.position=c(0.85, 0.2),  
        	legend.background= element_rect(fill="white", colour="white"),
        	legend.key= element_blank(),
        	legend.title= element_blank(),
        	legend.text = element_text(size=14))



# 2b- CUT - Effect size
#######################################################
df<- data.frame(name=c("Brucellosis", "Tuberculosis", "Site (LS)", "Age (< 3 yr)"), 
	est= c(3.0, 2.81, 2.08, 3.26), 
	lower= c(1.52, 1.43, 1.1, 1.70), 
	upper= c(6.0, 5.58, 3.93, 6.28), 
	order=c(1,2,3,4))
df$name <- factor(df$name, levels= df$name[order(df$order, decreasing = TRUE)])

p7 <- ggplot(df, aes(x=df$name, y=df$est)) + 
	geom_point(df$estimate, size = 3, shape = 19)+ 
	geom_errorbar(aes(ymin=df$lower, ymax=df$upper, width=0.1)) +  				
    theme_bw() +
    theme(panel.border= element_blank(), 
          axis.title.x=element_text(size=14), axis.title.y=element_blank() ) + 
    theme(axis.line.x = element_line(color="black"), 
          axis.line.y = element_line(color="black")) + 
	coord_flip() + 
	ylab("Relative risk of mortality") +
	geom_segment(aes(x=0, xend=4.3, y=1, yend=1), linetype=2, colour = "dark red")+ ylim(-0.01,6.5) + 
	theme(
		axis.title.x = element_text(size=16, vjust=-0.15),
        axis.title.y = element_blank(), #element_text(size=16, vjust= 0.8),
        axis.text.x = element_text(size=14, vjust=-0.05),
        axis.text.y = element_text(size=14))





#######################################################
# Survival curves- not used but informative
#######################################################
# plot predicted survival time (same as before but good): 
df<-data.frame(start=data3$start2, stop=data3$stop2, death.time=data3$death.time, TB=data3$TB_3, Br=data3$brucella, age=data3$age6, herd=data3$herd)
# set at average herd value (0.53) so not LS or CB. 
mort<-with(df, data.frame(tb=0, bruc=0, age=rep(mean(age=="adult"), 2), herd=rep(mean(herd=="LS"), 2)))
mortTB<-with(df, data.frame(TB=1, Br=0, age=rep(mean(age=="adult"),2), herd=rep(mean(herd=="LS"), 2)))
mortBr<-with(df, data.frame(Br=1, TB=0, age=rp(mean(age=="adult"),2), herd=rep(mean(herd=="LS"), 2)))
mortco<-with(df, data.frame(Br=1, TB=1, age=rep(mean(age=="adult"),2), herd=rep(mean(herd=="LS"), 2)))
mort<-with(df, data.frame(tb=0, bruc=0, age=0, herd=1))
mortTB<-with(df, data.frame(tb=1, bruc=0, age=0, herd=1))
mortBr<-with(df, data.frame(bruc=1, tb=0, age=0, herd=1))# herds = 0 or 1 (LS)
mortco<-with(df, data.frame(bruc=1, tb=1, age=0, herd=1))

plot_add.mod<-coxph(Surv(start, stop, death.time)~bruc+herd+tb+herd+age, data=df)
m<-survfit(plot_add.mod, newdata=mort)  # adults, herd - CB
mt<-survfit(plot_add.mod, newdata=mortTB)
mb<-survfit(plot_add.mod, newdata=mortBr)
mco<-survfit(plot_add.mod, newdata=mortco)


plot(m, conf.int=FALSE, ylab="Survival", xlab="Time (months)", lty=c(1, 2),
	ylim=c(0.1, 1), cex.lab=1.4, col="dark blue", bty="n")
lines(mt, lty=c(5, 5), conf.int= FALSE, col="dark green")
lines(mb, lty=c(4, 4), conf.int= FALSE, col="purple")
lines(mco, lty=c(3,3), conf.int= FALSE, col="dark red")
legend("bottomleft", legend=c("Uninfected", "Tuberculosis+", "Brucellosis +", "Co-infected"),
	lty=c(1 ,5, 4, 3), inset=0.02, bty="n", col=c("dark blue", "dark green", "purple", "dark red"))


df<- data.frame(name=c("bruc", "bTB", "site", "age (< 3 yr)"), 
	est= c(3.0, 2.81, 2.08, 3.26), 
	lower= c(1.52, 1.43, 1.1, 1.70), 
	upper= c(6.0, 5.58, 3.93, 6.28), 
	order=c(1,2,3,4))
df$name <- factor(df$name, levels= df$name[order(df$order, decreasing = TRUE)])

p2 <-ggplot(df, aes(x=df$name, y=df$est)) + 
	geom_point(df$estimate)+ 
	geom_errorbar(aes(ymin=df$lower, ymax=df$upper, width=0.1)) +  				
	theme_bw() +
	theme(panel.border= element_blank(), 
	axis.title.x=element_text(size=14), axis.title.y=element_blank() ) + 
	theme(axis.line.x = element_line(color="black"), 
	axis.line.y = element_line(color="black")) + 
	coord_flip() + 
	ylab("Relative risk of mortality") +
	geom_segment(aes(x=0, xend=4.3, y=1, yend=1), linetype=2, colour = "dark red")+ ylim(-0.01,6.5) + 
	coord_flip() 
p2


########  OVERALL SURVIVAL ESTIMATES FOR OUR POPULATION
df<-data.frame(start=data3$start2, stop=data3$stop2, death.time=data3$death.time, age=data3$age6, herd=data3$herd)
mort<-with(df, data.frame(age=c("adult", "adult", "juvenile", "juvenile"), herd=c("LS", "CB", "LS", "CB")))
plot_add.mod<-coxph(Surv(start, stop, death.time)~herd+age, data=df)
m<-survfit(plot_add.mod, newdata=mort)
summary(m)

plot_add.mod<-coxph(Surv(start, stop, death.time)~age, data=df)
m<-survfit(plot_add.mod, newdata=mort)
summary(m)




########################################################################
########################################################################
# Figure 2: Fecundity Figure
########################################################################
########################################################################
newdf<-data.frame(
	Calf=c(14/25, 10/33, 7/23, 6/14, 2/26, 0/8, 2/7, 1/7), 
	Agecategory = c(rep("Adult (age > 4)", 4), rep("Juvenile (age = 4)", 4)),
	Infection = c("uninfected", "brucellosis", "bTB", "co-infected", "uninfected", "brucellosis", "bTB", "co-infected"),
	N = c(25, 33, 23, 14, 26, 8, 7, 7))
newdf$color <- as.factor(seq(1,8))
newdf$se<- sqrt(newdf$Calf * (1 - newdf$Calf) / newdf$N)
newdf$order <- c(1,3,2,4,1,3,2,4)
newdf$Infection <- as.factor(newdf$Infection)
newdf$Infection <- factor(newdf$Infection, levels = newdf$Infection[order(unique(newdf$order))])

p9<- ggplot(newdf, aes(x=Infection, y=Calf, group=Agecategory, colour=Agecategory)) + 
  #geom_line()+
  geom_point(size=3, shape=19) + # colour="darkred", fill="darkred" +
  geom_errorbar(aes(ymin= newdf$Calf-newdf$se, ymax=newdf$Calf+newdf$se), width=0.2) + 
  scale_colour_manual(values=c("darkslategray", "darkseagreen3"))
p10<- p9 +
  theme_bw() + # removes ugly gray.
  ylab("Fecundity") +  # Proportion of buffalo observed with a calf
  xlab("")+
  scale_y_continuous(limits=c(0,0.8)) + 
  theme(axis.line.x = element_line(colour= "black"),
  		axis.line.y = element_line(colour= "black"),
  		axis.title.x = element_text(size=16, vjust=-0.15),
        axis.title.y = element_text(size=16, vjust= 0.8),
        axis.text.x = element_text(size=14, vjust=-0.05),
        axis.text.y = element_text(size=14),
        panel.border = element_blank(), 
        # legend information
        legend.position=c(0.8, 0.9),  
        legend.background= element_rect(fill="white", colour="white"),
        legend.key= element_blank(),
        legend.title= element_blank(),
        legend.text = element_text(size=14)   )
        
# With Bree's color scheme:  
newdf2 <- newdf[c(1:4),]
p9<- ggplot(newdf2, aes(x=Infection, y=Calf, colour=Infection)) + 
	geom_point(position = position_dodge(0.2), size=4, shape=19) + 
	geom_errorbar(aes(ymin= newdf2$Calf-newdf2$se, ymax=newdf2$Calf+newdf2$se), 
  		position = position_dodge(0.2), width=0.2) + 
	scale_colour_manual(values=c("goldenrod1", "slateblue3", "chartreuse4","tomato3"), guide = FALSE)       
p10<- p9 +  theme_bw() + # removes ugly gray.
  ylab("Fecundity") +  # Proportion of buffalo observed with a calf
  scale_y_continuous(limits=c(0,0.8))  + 
  theme(axis.line.x = element_line(colour= "black"),
  		axis.line.y = element_line(colour= "black"),
  		axis.title.x = element_blank(),#(size=16, vjust=-0.15),
        axis.title.y = element_text(size=16, vjust= 0.8),
        axis.text.x = element_text(size=14, vjust=-0.05),
        axis.text.y = element_text(size=14),
        panel.border = element_blank(),
        # legend information
        legend.position=c(0.8, 0.9),  
        legend.background= element_rect(fill="white", colour="white"),
        legend.key= element_blank(),
        legend.title= element_blank(),
        legend.text = element_text(size=14)) 


# With juvenile ones on!
p9 <- ggplot(newdf, aes(x=Infection, y=Calf, colour = as.character(color), shape = Agecategory, alpha = as.character(color))) + 
	geom_point(position = position_dodge(0.3), size=4) + 
	geom_errorbar(aes(ymin= newdf$Calf-newdf$se, ymax=newdf$Calf+newdf$se), 
  		position = position_dodge(0.3), width=0.2) + 
	scale_colour_manual(values=c("goldenrod1", "slateblue3", "chartreuse4",
		"tomato3", "goldenrod1", "slateblue3", "chartreuse4","tomato3"), 
		guide = FALSE) + 
	xlab("") +
	scale_alpha_manual(values = c(1,1,1,1,0.5, 0.5, 0.5, 0.5), guide = FALSE)       
p10.2<- p9 +  theme_bw() + # removes ugly gray.
  ylab("Fecundity") +  # Proportion of buffalo observed with a calf
  scale_y_continuous(limits=c(0,0.8))  + 
  theme(axis.line.x = element_line(colour= "black"),
  		axis.line.y = element_line(colour= "black"),
  		axis.title.x = element_text(size=16, vjust=-0.15),
        axis.title.y = element_text(size=16, vjust= 0.8),
        axis.text.x = element_text(size=14, vjust=-0.05),
        axis.text.y = element_text(size=14),
        panel.border = element_blank(),
        # legend information
        legend.position=c(0.8, 0.9),  
        legend.background= element_rect(fill="white", colour="white"),
        legend.key= element_blank(),
        legend.title= element_blank(),
        legend.text = element_text(size=12)) 



# CUT - Fecundity Effect size
#######################################################
df2<- data.frame(name=c("Tuberculosis\n (LS)", "Tuberculosis\n (CB)", 
	"Site (LS)", "Age (< 3 yr)"), 
	est= c(4.32, 0.39, 0.49, 2.42), 
	lower= c(1.51, 0.04062, 0.18, 1.04), 
	upper= c(12.3, 3.69, 1.32, 5.62), 
	order=c(1,2,3,4))
	df2$name <- factor(df2$name, levels= df2$name[order(df2$order, decreasing = TRUE)])

p11 <- ggplot(df2, aes(x=df2$name, y=df2$est)) + 
	geom_point(df2$estimate, size = 3, shape = 19)+ 
	geom_errorbar(aes(ymin=df2$lower, ymax=df2$upper, width=0.1)) +  				
    theme_bw() +
    theme(panel.border= element_blank(), 
          axis.title.x=element_text(size=14), axis.title.y=element_blank() ) + 
    theme(axis.line.x = element_line(color="black"), 
          axis.line.y = element_line(color="black")) + 
	coord_flip() + 
	ylab("Relative risk of brucellosis infection") +
	geom_segment(aes(x=0, xend=4.3, y=1, yend=1), linetype=2, colour = "dark red") +
	ylim(-0.01,12.31) + 
	theme(
		axis.title.x = element_text(size=16, vjust=-0.15),
        axis.title.y = element_blank(), #element_text(size=16, vjust= 0.8),
        axis.text.x = element_text(size=14, vjust=-0.05),
        axis.text.y = element_text(size=14))
source('~/GitHub/bTB-bruc-co-infection-ms/multiplot.R', chdir = TRUE)



########################################################################
########################################################################
# Figure 2: Incidence Figure
########################################################################
########################################################################

# Number of new infections
########################################################################
df <- data.frame(incid = c(0,0, 1,4,14,0,4,5,8,4,4,7,4,4,2,3), captureperiod = rep(c("1", "2", "3", "4", "5", "6", "7", "8"), each=2), infection = c(rep(c("bTB", "brucellosis"), 8)))
df$infection <- relevel(df$infection, "bTB")

pincid <- ggplot(df, aes(x = captureperiod, y = incid, fill = infection, colour = infection)) + 
	geom_bar(position = position_dodge(), stat = "identity") + ylim(0, 14) + 
	theme_bw() +  scale_fill_manual(values=c("slateblue3", "chartreuse4")) +
	scale_colour_manual(values=c("slateblue4", "darkgreen")) +
	ylab("Number of new infections") + xlab("Capture period") + 
	theme(axis.line.x = element_line(colour= "black"),
  		axis.line.y = element_line(colour= "black"),
  		axis.title.x = element_text(size=16, vjust=-0.15),
        axis.title.y = element_text(size=16, vjust= 0.8),
        axis.text.x = element_text(size=14, vjust=-0.05),
        axis.text.y = element_text(size=14),
        panel.border = element_blank(), 
		legend.position=c(0.85, 0.9),  
        legend.background= element_rect(fill="white", colour="white"),
        legend.key= element_blank(),
        legend.title= element_blank(),
        legend.text = element_text(size=12))

#df <- data.frame(incid = c(1,2,3,4), infection = c("bTB", "bTB", "brucellosis", "brucellosis"),
#	status = c("Single", "Co-infected", "Single", "Co-infected"))
#df$infection <- relevel(df$infection, "bTB")
#df$status <- relevel(df$status, "Single")

# Incidence rate
########################################################################
df2 <- data.frame(Infection = c("brucellosis", "brucellosis", "bTB", "bTB"), 
	Coinfection = c("uninfected", "infected", "uninfected", "infected"),
	Incidence = c(0.052, 0.083, 0.081, 0.07),
	N = c(173+189-14-8-6, 58+57-6-5-7, 173+189-14-8-6, 90+61-10-1-4)) 
df2$SE <- sqrt(df2$Incidence*(1-df2$Incidence)/df2$N)
df2$color <- as.character(seq(1,4,1))
df2$Infection <- relevel(df2$Infection, "bTB")
df2$Coinfection <- relevel(df2$Coinfection, "uninfected")

p12 <- 
ggplot(df2, aes(x = Infection, y = Incidence, colour = color, shape = Coinfection)) + #colour = color,
	geom_point(position = position_dodge(0.3), stat = "identity", size =4) +
	geom_errorbar(aes(ymin=df2$Incidence - df2$SE, ymax= df2$Incidence + df2$SE, width=0.1), 
		position = position_dodge(0.3)) +
	#geom_line(aes(group =Infection), position = position_dodge(0.3), colour = "black") + 
	ylim(0, 0.13) + 
	xlab("") +
	theme_bw() +
	scale_colour_manual(values = c("goldenrod1","slateblue3","goldenrod1","chartreuse4"), guide = FALSE) +
	#scale_fill_manual(values=c("goldenrod1","slateblue3","goldenrod1","chartreuse4"), guide = FALSE) +  			#"slateblue4","darkgreen"
	#scale_colour_manual(values = c("goldenrod3","slateblue4","goldenrod3","darkgreen"), guide = FALSE) +
	theme(axis.line.x = element_line(colour= "black"),
  		axis.line.y = element_line(colour= "black"),
  		axis.title.x = element_text(size=16, vjust=-0.15),
         axis.title.y = element_text(size=16, vjust= 0.8),
         axis.text.x = element_text(size=14, vjust=-0.05),
         axis.text.y = element_text(size=14),
         panel.border = element_blank(), 
	legend.position=c(0.15, 0.9),  
        legend.background= element_rect(fill="white", colour="white"),
        legend.key= element_blank(),
        legend.title= element_blank(),
        legend.text = element_text(size=12)) 
        
# With effect size (95% CI)
df = data.frame(infection = c("bTB", "brucellosis"), y = c(1.3, 2.13), ymin = c(0.65, 0.91), ymax = c(2.6, 4.9))
# SE
df = data.frame(infection = c("bTB", "brucellosis"), y = c(1.3, 2.13), ymin = c(0.91, 1.38), ymax = c(1.85, 3.29))
df$infection<- relevel(df$infection, "bTB")
p13 <- ggplot(df, aes(x = infection, y = y, colour = infection)) +
	geom_point(stat = "identity", size = 4, shape = 15) + 
	geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.1)+
	scale_colour_manual(values = c("slateblue3","chartreuse4"), guide = FALSE) + 
	theme_bw() +
	xlab("") +
	ylim(0, 6)+
	geom_hline(yintercept = 0.98, linetype = 2) +
	#ylab("Proportional change in infection risk") +
	ylab("Change in infection risk") +
	theme(axis.line.x = element_line(colour= "black"),
  		axis.line.y = element_line(colour= "black"),
  		axis.title.x = element_text(size=16, vjust=-0.15),
         axis.title.y = element_text(size=16, vjust= 0.8),
         axis.text.x = element_text(size=14, vjust= -0.05),
         axis.text.y = element_text(size=14, vjust = 0.05),
         panel.border = element_blank())


# inset plot
inset <- ggplot(df2, aes(x = Infection, y = Incidence, colour = Infection, shape = Coinfection)) + #colour = color,
	geom_point(position = position_dodge(0.4), stat = "identity", size = 3) +
	ylim(0, 0.11) + 
	xlab("") +
	theme_bw() +
	scale_colour_manual(values = c("slateblue3","chartreuse4"), guide = FALSE) +
	theme(axis.line.x = element_line(colour= "black"),
  		axis.line.y = element_line(colour= "black"),
  		axis.title.x = element_blank(),
         	axis.title.y = element_text(size=16, vjust= 0.8),
         	axis.text.x = element_blank(),
         	axis.text.y = element_text(size=14),
         	panel.border = element_blank(), 
		panel.margin = element_blank(),
        	panel.grid.major= element_blank(),
        	panel.grid.minor= element_blank()) + 
        	guides(colour = FALSE, shape = FALSE)
		#legend.position=c(0.15, 0.9),  
        	#legend.background= element_rect(fill="white", colour="white"),
        	#legend.key= element_blank(),
        	#legend.title= element_blank(),
        	#legend.text = element_text(size=12)) 


	
2.13 (95% CI 0.91 - 4.98)  # bTB on brucellosis 
1.30  (95% CI 0.65 - 2.61)
 
 # Figure 2, take 2
 multiplot(p6red, pincid, p13, p10.2, cols =2)
 multiplot(p6red, p13, cols =2)

#vp<- viewport(width=0.21, height=0.17, x=0.66, y=0.89)
#print(inset, vp=vp)


        
# FIGURE 2 (Full)
#multiplot(p6,  pincid, p10.2, p12 , cols = 2)



#######################################################
#######################################################
# Figure 3b- co-infection patterns
#######################################################
#######################################################
setwd("~/GitHub/bTB-bruc-co-infection-ms")
df<- readRDS("coinfection_pattern_results.rds")
f2 <- data.frame(Evaluation = c("Model", "Model", "Data", "Data"), 
	Infection = c("Single infection", "Co-infection", "Single infection", "Co-infection"), 
	BrucellosisPrevalence = c(df$prevBinS, df$prevBinT, 0.3035, 0.4524), 
	TBPrevalence = c(df$prevTinS, df$prevTinB, 0.227, 0.3585), 
	CILowB = c(NA, NA, 0.2654405, 0.3657982), 
	CIHighB = c(NA, NA, 0.3451408, 0.5405985), 
	CILowTB = c(NA, NA, 0.1781895, 0.2857143), 
	CIHighTB = c(NA, NA, 0.2783505, 0.4313725))
df2 <- reshape(f2, varying = c("BrucellosisPrevalence", "TBPrevalence"), 
	v.names = "Prevalence", timevar = "Name", times = c("brucellosis", "bTB"), 
	direction = "long")
df2$color <- c(1,2,0, 0, 3, 4, 0, 0)
df2$Infection <- relevel(df2$Infection, "Single infection")
df2$CILow <- c(NA, NA, 0.26544, 0.365798, NA, NA, 0.1781895, 0.2857)
df2$CIHigh <- c(NA, NA, 0.3451408, 0.540598, NA, NA, 0.27835, 0.43137)
df3 <- df2[df2$Evaluation == "Model",]
df3$color <- c(rgb(97, 139, 14, maxColorValue= 255))
df3$Name <- as.factor(df3$Name)
df3$Name <- relevel(df3$Name, "bTB")

limits <- aes(ymax = df2$CIHigh[df2$Evaluation == "Data"], ymin= df2$CILow[df2$Evaluation == "Data"])
p1 <- 
ggplot(df3, aes(x = Name, y = Prevalence, fill = Infection)) + #, colour = Name
	geom_bar(position = position_dodge(), stat = "identity", colour = "black") +
	ylim(0, 0.8) + xlab("") +
	theme_bw() +  scale_fill_manual(values=c("gray83", "gray47")) + 
	theme(axis.line.x = element_line(colour= "black"),
  		axis.line.y = element_line(colour= "black"),
  		axis.title.x = element_text(size=16, vjust=-0.15),
        axis.title.y = element_text(size=16, vjust= 0.8),
        axis.text.x = element_text(size=14, vjust=-0.05),
        axis.text.y = element_text(size=14),
        panel.border = element_blank(), 
		legend.position=c(0.8, 0.9),  
        legend.background= element_rect(fill="white", colour="white"),
        legend.key= element_blank(),
        legend.title= element_blank(),
        legend.text = element_text(size=12))  + 
	geom_point(data = df2[df2$Evaluation == "Data",], colour = "black", 
		position = position_dodge(width = 0.9), size = 3)  +
	geom_errorbar(limits, position = position_dodge(width = 0.9), width= 0.2)      
# turn off dots in legend.  Add model/data to legend



p1


#######################################################
#######################################################
# Figure 3a- time sequence
#######################################################
#######################################################
source('~/GitHub/bTB-bruc-co-infection-ms/fixed_parameters_norecovery_agematrix.R', chdir = TRUE)
source('~/GitHub/bTB-bruc-co-infection-ms/rhs_age.R', chdir = TRUE)
s_index <- 1:20
it_index <- 21:40
ib_index <- 41:60
ic_index <- 61:80
r_index <- 81:100
rc_index <- 101:120
juveniles <- 1:3
subadult<- 4
adult <- 5:14
mature <- 15:20

# Age structure information, used to calculate mortality rates in susceptibles. 
relageall = c(0.137, rep(0.368/4, 4), rep(0.185/4, 4),  # Jolles 2007, set max age at 20
	rep(0.235/6, 6), rep(0.075/5, 5))					# Also in Caron et al. from 2001 KNP	
relage = relageall
# Get stable age structure without dz
S0 = 400*relage; It0 = 0*relage; Ib0 = 0*relage; 
Ic0 = 0*relage; R0 = 0 * relage; Rc0 = 0 * relage
x0 = c(S0, It0, Ib0, Ic0, R0, Rc0)
times <- seq(0, 500, 1)
params <- c(fixed.params.olddz, list(gamma=1/2, betaB = 1.004592,   
	betaT = 12.833531/10000, rhoT = 1, rhoB = 2.1, theta= 4, K = 433))  
sol <- as.data.frame(ode(x0, times, rhs_age_matrix, params))
stable_age <- unname(unlist( sol[500, c(2:21)]/sum(sol[500, c(2:21)]) ))
# Add Brucellosis, get endemic age structure
S0 = 400* stable_age; It0 = 0 * stable_age; Ib0 = 20* stable_age; 
Ic0 = 0* stable_age; R0 = 30 * stable_age; Rc0 = 0 * stable_age
x0 = c(S0, It0, Ib0, Ic0, R0, Rc0)
times <- seq(0, 500, 1)
sol <- as.data.frame(ode(x0, times, rhs_age_matrix, params))
endemic_agestructure <- unname(unlist( sol[500, c(2:121)] ))
# Introduce bTB
x0 = endemic_agestructure
x0[28] <- 5; x0[8] <- x0[8] - 5
sol <- as.data.frame(ode(x0, times, rhs_age_matrix, params))
sol$bTBoverall <- apply(sol[it_index+1], 1, sum)
sol$brucoverall <- apply(sol[ib_index+1], 1, sum) + apply(sol[r_index+1], 1, sum)
sol$co <- apply(sol[ic_index+1], 1, sum) + apply(sol[rc_index+1], 1, sum)
sol$Ib <- apply(sol[ib_index+1], 1, sum)
sol$It <- apply(sol[it_index+1], 1, sum)
sol$Ic <- apply(sol[ic_index+1], 1, sum)
sol$Rb <- apply(sol[r_index+1], 1, sum)
sol$Rc <- apply(sol[rc_index+1], 1, sum)

bTB <- data.frame(N = sol$bTBoverall[sol$time < 100], time = sol$time[sol$time < 100])
bruc <- data.frame(N = sol$bruc[sol$time < 100], time= sol$time[sol$time < 100])
co <- data.frame(N = sol$co[sol$time < 100], time = sol$time[sol$time < 100])
df <- data.frame(N = c(bTB$N, bruc$N, co$N), time = rep(bTB$time, 3), label = rep(c("bTB only", "brucellosis only", "co-infected"), each = length(bTB$N)) )
df$label <- relevel(df$label, "bTB")
p <- ggplot(df, aes(x = time, y = N, colour = label)) + 
	geom_line(size = 1.3) + ylim(0, 150) + #position = position_dodge(), stat = "identity"
	theme_bw() + xlab("Time (years)") + ylab("Number of buffalo") + 
	scale_colour_manual(values = c("slateblue3", "chartreuse4", "tomato3")) +
	theme(axis.line.x = element_line(colour= "black"),
		axis.line.y = element_line(colour= "black"),
		axis.title.x = element_text(size=16, vjust=-0.15),
		axis.title.y = element_text(size=16, vjust= 0.8),
		axis.text.x = element_text(size=14, vjust=-0.05),
       axis.text.y = element_text(size=14),
       panel.border = element_blank(), 
		legend.position=c(0.8, 0.9),  
       legend.background= element_rect(fill="white", colour="white"),
       legend.key= element_blank(),
       legend.title= element_blank(),
       legend.text = element_text(size=12)) 
p

multiplot(p, p1, cols = 2)

#######################################################
#######################################################
# Figure 4- Levelplots, varying mort and rho for both pathogens
#######################################################
#######################################################
epiT <- read.csv("~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/draft2/post-labmeeting/epiT.csv")
epiB <- read.csv("~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/draft2/post-labmeeting/epiB.csv")
epiB <- epiB[epiB$mort < 7,]
epiT <- epiT[epiT$mort < 7,]
cols <- brewer.pal(11, "RdBu")
cols2 <- colorRampPalette(brewer.pal(11, "RdBu"))

epiT$bTBplot <- epiT$bTBprev - 0.58
epiB$bTBplot <- epiB$bTBprev - 0.58
epiT$brucplot <- epiT$brucprev - 0.318
epiB$brucplot <- epiB$brucprev - 0.318

# see range of data by self overall
p1 <- ggplot(data = epiT, aes(x = mort, y = rhoT)) + geom_tile(aes(fill = bTBprev)) +
	scale_fill_distiller(palette = "RdYlBu", direction = -1, limits = c(0, 0.8))
#p2<- ggplot(data = epiT, aes(x = mort, y = rhoB)) + geom_tile(aes(fill = bTBplot)) +
#	scale_fill_distiller(palette = "RdYlBu", direction = -1, limits = c(-0.51, 0.51))	
p3 <- ggplot(data = epiB, aes(x =mort, y = rhoB)) + geom_tile(aes(fill = brucprev)) +
	scale_fill_distiller(palette = "RdYlBu", direction = -1, limits = c(0, 0.8))
#p4 <- ggplot(data = epiB, aes(x =mort, y = rhoT)) + geom_tile(aes(fill = bTBplot)) +
#	scale_fill_distiller(palette = "RdYlBu", direction = -1, limits = c(-0.51, 0.51))
multiplot(p1, p3, cols= 2)


# For publication
df <- data.frame(Difference = c(epiT$bTBplot, epiB$brucplot), 
	infection = c(rep("bTB", length(epiT[,1])), rep("brucellosis", length(epiB[,1]))), 
	rho = c(epiT$rhoT, epiB$rhoB), mort = c(epiT$mort, epiB$mort))
df2 <- rbind(df[df$infection == "brucellosis" & df$rho == 2.1 & df$mort == 5.8 , c("rho", "mort", "infection")], 
	df[df$infection == "bTB" & df$rho == 1 & df$mort == 5.8 , c("rho", "mort", "infection"),])
df$infection <- relevel(df$infection, "bTB")
p <- ggplot(data = df, aes(x = mort, y = rho))
p2 <- p + theme_bw() + facet_wrap(~ infection) + 
	xlab(expression(paste(Proportional~increase~"in"~mortality~with~"co-infection,"~ mu[C], "/", mu[S]))) + 
	ylab("Proportional increase in transmission with co-infection") +
	theme(panel.grid.major = element_blank(), 
	plot.title = element_text(size = 14),
	axis.title = element_text(size = 12), axis.text = element_text(size = 12), 
	legend.text = element_text(size=12), 
	legend.text.align = 1, 
	legend.key.size = unit(1, "cm") ) + 
	geom_tile(aes(fill = Difference)) +
	scale_fill_distiller(palette = "RdYlBu", direction = -1, limits = c(-0.51, 0.51)) +
	geom_point(data = df2, size = 2, pch = 8)
	#geom_rect(data = df2, aes(x = NULL, y = NULL, xmin = mort - 0.5, ymin = rho - 0.5, 
	#	xmax = mort + 0.5, ymax = rho + 0.5), size = 1, fill = NA, 
	#	colour = "black")
p2


# Just realistic mortality (or mort with same as single to greater)
df <- df[df$mort > 2.8,]
df$infection<- relevel(df$infection, "bTB")
p <- ggplot(data = df, aes(x = mort, y = rho))
p2 <- p + theme_bw() + facet_wrap(~ infection) + 
	xlab(expression(paste(Proportional~increase~"in"~mortality~with~"co-infection,"~ mu[C], "/", mu[S]))) + 
	ylab("Proportional increase in transmission with co-infection") +
	theme(panel.grid.major = element_blank(), 
	plot.title = element_text(size = 14),
	axis.title = element_text(size = 12), axis.text = element_text(size = 12), 
	legend.text = element_text(size=12), 
	legend.text.align = 1, 
	legend.key.size = unit(1, "cm") ) + 
	geom_tile(aes(fill = Difference)) +
	scale_fill_distiller(palette = "RdYlBu", direction = -1, limits = c(-0.35, 0.35)) +
	geom_point(data = df2, size = 1, pch = 8)
	#geom_rect(data = df2, aes(x = NULL, y = NULL, xmin = mort - 0.5, ymin = rho - 0.5, 
	#	xmax = mort + 0.5, ymax = rho + 0.5), size = 1, fill = NA, 
	#	colour = "black")
p2



# Just realistic mortality and prevalence
df <- data.frame(Prevalence = c(epiT$bTBprev, epiB$brucprev), 
	infection = c(rep("bTB", length(epiT[,1])), rep("brucellosis", length(epiB[,1]))), 
	rho = c(epiT$rhoB, epiB$rhoT), mort = c(epiT$mort, epiB$mort))
df2 <- rbind(df[df$infection == "brucellosis" & df$rho == 2.1 & df$mort == 5.8 , c("rho", "mort", "infection")], 
	df[df$infection == "bTB" & df$rho == 1 & df$mort == 5.8 , c("rho", "mort", "infection"),])
df <- df[df$mort > 2.8,]
p <- ggplot(data = df, aes(x = mort, y = rho))
p2 <- p + theme_bw() + facet_wrap(~ infection) + 
	xlab(expression(paste(Proportional~increase~"in"~mortality~with~"co-infection,"~ mu[C], "/", mu[S]))) + 
	ylab("Proportional increase in transmission with co-infection") +
	theme(panel.grid.major = element_blank(), 
	plot.title = element_text(size = 14),
	axis.title = element_text(size = 12), axis.text = element_text(size = 12), 
	legend.text = element_text(size=12), 
	legend.text.align = 1, 
	legend.key.size = unit(1, "cm") ) + 
	geom_tile(aes(fill = Prevalence)) +
	scale_fill_distiller(palette = "RdYlBu", direction = -1)   #, limits = c(-0.4, 0.06)
	#geom_point(data = df2, size = 1, pch = 8)
	#geom_rect(data = df2, aes(x = NULL, y = NULL, xmin = mort - 0.5, ymin = rho - 0.5, 
	#	xmax = mort + 0.5, ymax = rho + 0.5), size = 1, fill = NA, 
	#	colour = "black")
p2


p2<- ggplot(data = epiT, aes(x = mort, y = rhoB)) + geom_tile(aes(fill = bTBplot)) +
	scale_fill_distiller(palette = "RdYlBu", direction = -1, limits = c(-0.41, 0.41))	
p3 <- ggplot(data = epiB, aes(x =mort, y = rhoT)) + geom_tile(aes(fill = brucplot)) +
	scale_fill_distiller(palette = "RdYlBu", direction = -1, limits = c(-0.41, 0.41))
multiplot(p2, p3, cols= 2)


#######################################################
#######################################################
# Figure 4 old- Levelplots, varying rho for both pathogens
#######################################################
#######################################################
epi <- read.csv("~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/draft2/post-labmeeting/vary_rho/vary_rho.csv")
cols <- brewer.pal(11, "RdBu")
cols2 <- colorRampPalette(brewer.pal(11, "RdBu"))

#epi$bTBplot <- epi$bTBprev - 0.30
#epi$brucplot <- epi$brucprev - 0.30
epi$bTBplot <- epi$bTBprev - 0.58
epi$brucplot <- epi$brucprev - 0.318

df <- reshape(epi, direction = "long", varying = c("bTBplot", "brucplot"),
	v.names = "Difference", timevar = "infection", times = c("bTB", "brucellosis"))
p <- ggplot(data = df, aes(x = mort, y = rho))
p + theme_bw() + facet_wrap(~ infection) + 
	geom_tile(aes(fill = Difference)) +
	scale_fill_distiller(palette = "RdYlBu", direction = -1, limits = c(-0.40, 0.40)) + 
	xlab(expression(paste(Proportional~increase~"in"~brucellosis~
		transmission~with~"co-infection,"~ beta[B]^{"'"}, "/", beta[B]))) + 
	ylab(expression(paste(Proportional~increase~"in"~bTB~
		transmission~with~"co-infection,"~ beta[T]^{"'"}, "/", beta[T]) )) +
	theme(panel.grid.major = element_blank(), 
		plot.title = element_text(size = 14),
		axis.title = element_text(size = 12), axis.text = element_text(size = 12), 
		legend.text = element_text(size=12), 
		legend.text.align = 1, 
		legend.key.size = unit(1, "cm") )


	






########################################################################
########################################################################
# Fecundity Figure for supplememt with pfc in it. 
########################################################################
########################################################################
data<- read.csv("~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/final_datasets_copied_from_phdfolder/cross_sectional_data_withdz_cleandisease_nofinal_Feb2016_capturetime_forsurv.csv")
data2<- data[data$herdorig=="LS",]
calftime <- c(0, 12, 24, 36)
data4<- data2[data2$capturetime %in% calftime,]
data4$fec<- NA
table(data4$milk , data4$calf)
# With milk one 10 year old, two 12 year olds, no 9, 11, 0, or 1 yr olds (good)
# one 11 year old without milk
# called age <=10 yes
# called 10 unknowns with milk yes; 3 unknowns with no/unknown milks no
include = c( 10, 2, 3,4, 5, 6, 7, 8, 9, "yes")
data4$fec[data4$calf %in% include] <- 1 
data4$fec[data4$calf %in% c(11, 12, "no")] <- 0
data4$fec[data4$calf=="unknown" & data4$milk=="yes"] <- 1
data4$fec[data4$calf=="unknown" & data4$milk=="unknown"] <- 0
data4$fec[data4$calf=="unknown" & data4$milk=="no"] <- 0
d<- data4
d<- data4[data4$age_yr > 3 & data4$age_yr <= 9,]

# totals, denominator: 
table(d$tb[d$age1 == "adult"], d$bruc_beforeafter[d$age1 == "adult"])
table(d$tb[d$age1 == "juvenile"], d$bruc_beforeafter[d$age1 == "juvenile"])
# numerator
table(d$tb[d$age1 == "adult" & d$fec=="1"], d$bruc_beforeafter[d$age1 == "adult" & d$fec=="1"])
table(d$tb[d$age1 == "juvenile" & d$fec=="1"], d$bruc_beforeafter[d$age1 == "juvenile" & d$fec=="1"])


newdf<-data.frame(
	Calf=c(14/25, 1/5, 9/28, 7/23, 1/8, 5/6,  # adults
		 2/26, 0/3, 0/5, 2/7, 0/4, 1/3),  # juveniles
	Agecategory = c(rep("Adult (age > 4)", 6), rep("Juvenile (age = 4)", 6)),
	Infection = c("Uninfected", "BR Converter", "BR PFC", "bTB", "bTB & BR Converter", "bTB & BR PFC", "Uninfected", "BR Converter", "BR PFC", "bTB", "bTB & BR Converter", "bTB & BR PFC"),
	N = c(25, 5, 28, 23, 8, 6, 26, 3, 5, 7, 4, 3))
newdf$se<- sqrt(newdf$Calf * (1 - newdf$Calf) / newdf$N)
newdf$order <- c(seq(1, 6), seq(1,6))
newdf$Infection <- as.factor(newdf$Infection)
newdf$Infection <- factor(newdf$Infection, levels = newdf$Infection[order(unique(newdf$order))])


p11<- ggplot(newdf, aes(x=Infection, y=Calf, group=Agecategory, colour=Agecategory)) + 
  #geom_line()+
  geom_point(size=3, shape=19) + # colour="darkred", fill="darkred" +
  geom_errorbar(aes(ymin= newdf$Calf-newdf$se, ymax=newdf$Calf+newdf$se), width=0.2) + 
  scale_colour_manual(values=c("darkslategray", "darkseagreen3"))
p12<- p11 +
  theme_bw() + # removes ugly gray.
  ylab("Proportion of buffalo observed with a calf") +
  scale_y_continuous(limits=c(0,1.07)) + 
  theme(axis.line.x = element_line(colour= "black"),
  		axis.line.y = element_line(colour= "black"),
  		axis.title.x = element_text(size=16, vjust=-0.15),
        axis.title.y = element_text(size=16, vjust= 0.8),
        axis.text.x = element_text(size=14, vjust=-0.05),
        axis.text.y = element_text(size=14),
        panel.border = element_blank(), 
        # legend information
        legend.position=c(0.7, 0.9),  
        legend.background= element_rect(fill="white", colour="white"),
        legend.key= element_blank(),
        legend.title= element_blank(),
        legend.text = element_text(size=15)   ) + 
        annotate("text", x= "Uninfected", y = 0.77, label = "N = 25, 26") +
        annotate("text", x= "BR Converter", y = 0.55, label = "N = 5, 3") + 
        annotate("text", x= "BR PFC", y = 0.55, label = "N = 28, 5") + 
        annotate("text", x= "bTB", y = 0.55, label = "N = 23, 7") + 
        annotate("text", x= "bTB & BR Converter", y = 0.55, label = "N = 8, 4") + 
        annotate("text", x= "bTB & BR PFC", y = 1.05, label = "N = 6, 3")
p12




# CUT FIGURES & INFORMATION
#############################################
#############################################
Summary information, from cross sectional data
#############################################
#############################################
###############################################
data2<-read.csv("cross_sectional_data_withdz_cleandisease_nofinal_Feb2016_capturetime_forsurv.csv")
data_nofinal<-data2[data2$final_capture=="0",] 

length(data_nofinal[,1])
length(unique(data_nofinal$id))
length(data_nofinal$id[data_nofinal$tb=="1"])
length(data_nofinal$id[data_nofinal$bruc=="positive"])

# brucellosis prevalence in TB+, TB-
TBneg<- data_nofinal[data_nofinal$tb=="0",]
TBpos<- data_nofinal[data_nofinal$tb=="1",]
length(TBneg$bruc[TBneg$bruc=="positive"])/length(TBneg$bruc)
length(TBpos$bruc[TBpos$bruc=="positive"])/length(TBpos$bruc)


# number concverting
brconverters<-data_nofinal[!(data_nofinal$bruc_beforeafter=="nc"),]
brconverters<-brconverters[!(brconverters$bruc_beforeafter=="pfc"),]
length(brconverters$id); length(unique(brconverters$id))  # 29 buffalo became infected with brucellosis. 
    
tbconverters<-data_nofinal[!(data_nofinal$tb_beforeafter=="nc"),]
tbconverters<-tbconverters[!(tbconverters$tb_beforeafter=="pfc"),]
length(tbconverters$id); length(unique(tbconverters$id)) # 44 with bTB

bothconvert<-tbconverters[tbconverters$brucconvert=="1",]
length(bothconvert$id); length(unique(bothconvert$id))

# average age and month of first infection
incidtb<-tbconverters[tbconverters$incid=="1",]
summary(incidtb$age_sel)  # 31.00   50.00   56.00   59.62   71.00   91.00 
quantile(incidtb$age_sel, c(0.05, 0.95))  # 41  86 
incidbr<-brconverters[brconverters$incidbr =="1",] # 3.4 to 7.5 yrs old by percentile
summary(incidbr$age_sel)  # 60.56
quantile(incidbr$age_sel, c(0.05, 0.95)) #42.9 85.8 

# season/month- when are they converting (note incidbr added by hand.)
incidbr$capid<-as.character(incidbr$capid)
for (i in 1:length(incidbr[,1])){
	incidbr$month[i]<-substr(strsplit(incidbr$capid[i], "-")[[1]][2], 1, 2)
}

tiff(filename="Bruc_conversion_month.tiff", width=480, height=480, units="px")
hist(as.numeric(incidbr$month), breaks=seq(1,12,1), xlab="Month", ylab="Number of brucellosis seroconversions", main="", col="light gray")
dev.off()

incidtb$capid<-as.character(incidtb$capid)
for (i in 1:length(incidtb[,1])){
	incidtb$month[i]<-substr(strsplit(as.character(incidtb$capid)[i], "-")[[1]][2], 1, 2)
}
tiff(filename="bTB_conversion_month.tiff", width=480, height=480, units="px")
hist(as.numeric(incidtb$month), breaks=seq(1,12,1), xlab="Month", ylab="Number of bTB seroconversions", main="", col="light gray")
dev.off()




###############################################################
###############################################################
# FIGURE 1 UPDATES-  SAME FIGURE BUT WITHOUT PFC's for Brucellosis 
# (to avoid differential antibody loss)
datanpfc <- data[data$bruc_beforeafter != "pfc",]
data_nofinal<-datanpfc[datanpfc$final_capture=="0",] 
d2<-data.frame(btb=data_nofinal$tb , bruc=as.character(data_nofinal$bruc), age=data_nofinal$age_sel/12, 
              id=data_nofinal$id)
d2<-d2[d2$age<14,]
binsize=2
agebins=c(seq(0, max(d2$age), binsize), max(d2$age+1))
agebins<- c(agebins[c(1:6)], 16.5)
d2$bruc2<-NA
newdf<-data.frame(agebin=c(agebins, agebins), Brucprev=NA, 
                  TB=c(rep("bTB-", length(agebins)),rep("bTB+", length(agebins)) ),
                  N=NA)  


for (i in 1:length(d2$bruc)){
  ifelse(d2$bruc[i]=="negative", d2$bruc2[i]<-0, d2$bruc2[i]<-1)
}
d2$bruc<-as.numeric(d2$bruc2)
for (i in 1:(length(agebins)-1)){
  neg<-d2[d2$btb==0,]
  pos<-d2[d2$btb==1,]
  d_neg<-d2[d2$age>=agebins[i] & d2$age<agebins[i+1] & d2$btb=="0",]
  d_pos<-d2[d2$age>=agebins[i] & d2$age<agebins[i+1] & d2$btb=="1",]
  newdf$Brucprev[i]<-get_prev(d_neg$bruc)
  newdf$Brucprev[i+length(agebins)]<-get_prev(d_pos$bruc)
  newdf$N[i]<- length(d_neg$bruc)
  newdf$N[i+length(agebins)]<- length(d_pos$bruc)
}
newdf<-newdf[newdf$agebin<16,]
newdf$se<-sqrt(newdf$Brucprev*(1-newdf$Brucprev)/newdf$N)
newdf$se[is.na(newdf$se)]<-0

# remove 0-2 bin because no TB+
newdf2<-newdf[newdf$N>1,]; newdf<-newdf2














#############################################
#############################################
# Notes of GLMM
#############################################
#############################################
# AICs hanshed out ran on previous dataset.  IN repot results are presented for the dataset with capture time that matches the survival times.
data <- read.csv("~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/final_datasets_copied_from_phdfolder/cross_sectional_data_withdz_cleandisease_nofinal_Feb2016_capturetime.csv")
data_nofinal<-data[data$final_capture=="0",] 
temp<-data_nofinal[data_nofinal$age_yr<14,]
rescale= function(col){
  new=NA
  for (i in 1:length(col)){
    new[i]<-(col[i]-mean(col))/ 2*sd(col)  
  }
  return(new)
}
# try to rescale to get model lwm2.2 to converge
temp$herd<-NA
temp$tb2<-rescale(temp$tb)
temp$age_yr2<-rescale(temp$age_yr)
temp$age_yr_2sq<-rescale(temp$age_yr*temp$age_yr)
temp$herd[temp$herdorig=="LS"] <-1
temp$herd[temp$herdorig=="CB"]<-0
temp$herd2<-rescale(temp$herd)

# NOTE: NOTE AGE AND TB ARE COLINEAR- NO BTB POSTIIVE BUFFALO LESS THAN 2.5 yrs... so remove three data points
# AND OLDEST bTB positive buffalo is 10.5; so 
temp2<- temp[temp$age_yr>2.5 & temp$age_yr<10.5,]
temp3<- temp[temp$age_yr>2 & temp$age_yr<10,]
#temp3<- temp2[!(temp2$id %in% c("B14", "B32", "O33")),]
temp4<- temp2[!(temp2$capturetime %in% c(0,3,6,9,12,15, 18, 21)),]
temp5<- temp3[!(temp3$capturetime %in% c(0,3,6,9,12,15, 18, 21)),]
temp6<- temp[!(temp$capturetime %in% c(0,3,6,9,12,15, 18, 21)),]
temp6<- temp3[!(temp3$capturetime %in% c(0,3,6,9,12,15)),]

# Final in paper now:
 t1<-glmmPQL(bruc~ age_yr2*tb + I(age_yr2^2), correlation= corAR1(form=~capturetime2|id), random= ~ 1|id, data=temp2, family=binomial); summary(t1)

 t1<-glmmPQL(bruc~ age_yr2*tb + I(age_yr2^2), correlation= corAR1(form=~capturetime2|id), random= ~ 1|id, data=temp4, family=binomial); summary(t1)
 
t1<-glmmPQL(bruc~ age_yr2*tb, correlation= corAR1(form=~capturetime2|id), random= ~ 1|id, data=temp5, family=binomial); summary(t1)

#****** In legend t1<-glmmPQL(bruc~ age_yr2*tb+I(age_yr2^2) , correlation= corAR1(form=~capturetime2|id), random= ~ 1|id, data=temp6, family=binomial); summary(t1)



# For each capture separately. 
head(temp2)
t<- temp3[temp3$capturetime %in% c(0, 3),]
t1<- temp3[temp3$capturetime %in% c(6, 9),]
t2<- temp3[temp3$capturetime %in% c(12, 15),]
t3<- temp3[temp3$capturetime %in% c(18, 21),]
t4<- temp3[temp3$capturetime %in% c(24, 28),]
t5<- temp3[temp3$capturetime %in% c(30, 33),]
t6<- temp3[temp3$capturetime %in% c(36, 39),]
t7<- temp3[temp3$capturetime %in% c(42, 45),]

t<- temp[temp $capturetime %in% c(0, 3),]
t1<- temp[temp $capturetime %in% c(6, 9),]
t2<- temp[temp $capturetime %in% c(12, 15),]
t3<- temp[temp $capturetime %in% c(18, 21),]
t4<- temp[temp $capturetime %in% c(24, 28),]
t5<- temp[temp $capturetime %in% c(30, 33),]
t6<- temp[temp $capturetime %in% c(36, 39),]
t7<- temp[temp $capturetime %in% c(42, 45),]


test<-glmer(bruc~ age_yr*tb + (1|id), data=t, family=binomial); summary(test) 
test<-glmer(bruc~ age_yr2*tb+ I(age_yr2^2)+(1|id), data=t, family=binomial); summary(test) 
test<-glmer(bruc~ age_yr2*tb+ I(age_yr2^2)+(1|id), data=t1, family=binomial); summary(test) 
test<-glmer(bruc~ age_yr2*tb+ I(age_yr2^2)+(1|id), data=t2, family=binomial); summary(test) 
test<-glmer(bruc~ age_yr2*tb+ I(age_yr2^2)+(1|id), data=t3, family=binomial); summary(test) 
# age interaction starts to be significant
test<-glmer(bruc~ age_yr2*tb+ I(age_yr2^2)+(1|id), data=t4, family=binomial); summary(test) 
test<-glmer(bruc~ age_yr2*tb+ I(age_yr2^2)+(1|id), data=t5, family=binomial); summary(test) 
test<-glmer(bruc~ age_yr2*tb+ I(age_yr2^2)+(1|id), data=t6, family=binomial); summary(test) 
test<-glmer(bruc~ age_yr2*tb+ I(age_yr2^2)+(1|id), data=t7, family=binomial); summary(test) 


test<-glmer(bruc~ age_yr2+tb+(1|id), data=t3, family=binomial); summary(test) 

get_prev = function(data){
	temp <- as.data.frame(table(data$tb))
	prev1 <- temp$Freq[temp$Var1 == "1"]/sum(temp$Freq)
	temp <- as.data.frame(table(data$bruc))
	prev2 <- temp$Freq[temp$Var1 == "positive"]/sum(temp$Freq)
	prev<- c(prev1, prev2)
	return(prev)
	}


# AICs need changed
t<-glmer(bruc~ floor(age_yr)+ I(floor(age_yr^2))+(1|id), data=temp2, family=binomial(link="logit")); summary(t) #319.2
t<-glmer(bruc~ floor(age_yr)+(1|id), data=temp2, family=binomial); summary(t)
t<-glmer(bruc~ floor(age_yr)+ I(floor(age_yr)^2)+ tb+(1|id), data=temp2, family=binomial); summary(t) #349.4 # small error
t<-glmer(bruc~ floor(age_yr)+ I(floor(age_yr)^2)*tb+(1|id), data=temp2, family=binomial); summary(t) #345.2
t<-glmer(bruc~ floor(age_yr)*tb+ I(floor(age_yr)^2)+ herdorig+(1|id), data=temp2, family=binomial); summary(t) #351.2, small error
t<-glmer(bruc~ floor(age_yr)+ I(floor(age_yr)^2)*tb+ herdorig+(1|id), data=temp2, family=binomial); summary(t) #346.5
t<-glmer(bruc~ floor(age_yr)+ I(floor(age_yr)^2)*tb+tb*herdorig+(1|id), data=temp2, family=binomial); summary(t) # 354.9
t<-glmer(bruc~ floor(age_yr)+ I(floor(age_yr)^2)*tb+floor(age_yr)*herdorig+(1|id), data=temp2, family=binomial); summary(t) # 354.9
t<-glmer(bruc~age_yr2+ I(age_yr2^2)*tb2+(1|id), data=temp2, family=binomial); summary(t) # 354.9

# USE MODEL WITH STANDARDIZED AGE & TB*age+ TB*age^2 interaction!!!!!!
t<-glmer(bruc~ age_yr2+ I(age_yr2^2)+(1|id), data=temp2, family=binomial); summary(t) # 351.1
t<-glmer(bruc~ age_yr2+ I(age_yr2^2)*tb+(1|id), data=temp2, family=binomial); summary(t) # 355.1
t<-glmer(bruc~ age_yr2*tb+ I(age_yr2^2)+(1|id), data=temp2, family=binomial); summary(t) # AIC = 351.7
t<-glmer(bruc~ age_yr2*tb+ I(age_yr2^2)*tb+(1|id), data=temp2, family=binomial); summary(t) # AIC = 351.6

t2<-glmer(bruc~ yr+(1|yr/id), data=temp4, family=binomial); summary(t2) # AIC = 351.7
t2<-glmer(bruc~ yr+(capturetime|id), data=temp4, family=binomial); summary(t2) # AIC = 351.7
t2<-glmer(bruc~ yr+(yr|herd2/id), data=temp2, family=binomial); summary(t2) # AIC = 351.7

#http://rpsychologist.com/r-guide-longitudinal-lme-lmer
# need glmmPQL for binomial data
library(MASS)
t<- glmmPQL(bruc~ yr, random = ~1|id, family=binomial, data= temp4); summary(t)
t<- glmmPQL(bruc~ capturetime, random = ~capturetime|id, family=binomial, data= temp4); summary(t)
t<- glmmPQL(bruc~ yr, random = ~1+yr|id, family=binomial, data= temp4); summary(t)
t<- glmmPQL(bruc~ yr, random = ~0+yr|id, family=binomial, data= temp4); summary(t) # works!


t<- glmmPQL(bruc~ tb+ yr, random = ~0+yr|id, family=binomial, data= temp2); summary(t)
t<- glmmPQL(bruc~ tb*age_yr2+ I(age_yr2^2)+ yr, random = ~0+yr|id, family=binomial, data= temp2); summary(t)
#ns = t<- glmmPQL(bruc~ tb*age_yr2+ tb*I(age_yr2^2)+ yr, random = ~0+yr|id, family=binomial, data= temp2); summary(t)

# Can nest herd even 
t<- glmmPQL(bruc~ tb+ yr, random = ~0+yr|herd2/id, family=binomial, data= temp2); summary(t)
t<- glmmPQL(bruc~ tb*age_yr2+ I(age_yr2^2)+ yr, random = ~0+yr|herd2/id, family=binomial, data= temp2); summary(t)

# AR1 -> NO- autocorrelaiton does not decrease to 0 as time lags. 
t<- glmmPQL(bruc~ yr, random = ~1|id, family=binomial, data= temp2); summary(t)
E<- residuals(t, type="normalized")
acf(E)
t<- glmmPQL(bruc~ capturetime, random = ~1|id, family=binomial, data= temp2); summary(t)
E<- residuals(t, type="normalized")
acf(E)

t<- glmmPQL(bruc~ capturetime, random = ~1|id, correlation= corAR1(), family=binomial, data= temp4); summary(t)
t<- glmmPQL(bruc~ capturetime, random = ~0+capturetime|id, correlation= corAR1(), family=binomial, data= temp2); summary(t)
t<- glmmPQL(bruc~ 1, random = ~capturetime|id, correlation= corAR1(), family=binomial, data= temp4); summary(t)

t<- glmmPQL(bruc~ capturetime + tb, random = ~capturetime, family=binomial, data= temp2); summary(t)
t<- glmmPQL(bruc~age_yr2 * tb, random = ~ capturetime |id, family=binomial, data= temp4); summary(t)


t<-glmer(bruc~ age_yr2*tb+ I(age_yr2^2)*tb+(capturetime|herd2), data=temp2, family=binomial); summary(t) # AIC = 351.6


t<-glmer(bruc~ age_yr2*tb + I(age_yr2^2)+ (capturetime|id), data=temp4, family=binomial); summary(t) 
t<-glmer(bruc~ age_yr2*tb + I(age_yr2^2)+ (1|capturetime/id), data=temp4, family=binomial); summary(t) 

######################################################################
######################################################################
######################################################################
######################################################################
t0<-glmmPQL(bruc~ age_yr2*tb + I(age_yr2^2), random= ~ 1|id, data=temp4, family=binomial); summary(t) 
library(car)
t1<-glmmPQL(bruc~ age_yr2*tb + I(age_yr2^2), correlation= corAR1(), random= ~ 1|id, data=temp4, family=binomial); summary(t) 

t2<-glmmPQL(bruc~ age_yr2*tb + I(age_yr2^2), correlation= corARMA(c(0.2, 0.2, 0.2), form=~capturetime|id, p=1, q=2), random= ~ 1|id, data=temp4, family=binomial); summary(t) 

t1<-glmmPQL(bruc~ age_yr2*tb + I(age_yr2^2), correlation= corAR1(form=~capturetime2|id), random= ~ 1|id, data=temp4, family=binomial); summary(t) 

t<-glmmPQL(bruc~ age_yr2*tb + I(age_yr2^2), random= ~ 1|id, data=temp4, family=binomial); summary(t) 
E<- residuals(t, type="normalized")
par(mfrow=c(1,2))
acf(E)
pacf(E)

durbinWatsonTest(glm(bruc~ age_yr2*tb + I(age_yr2^2), data=temp4, family=binomial ))


###############
#However, with 4 time points you probably won't be able to fit the autocorrelation. So I would first fit a mixed effects model without autocorrelation structure (probably using package lme4, but you can also use lme) and test autocorrelation of the residuals using the Durbin-Watson test,
#http://stats.stackexchange.com/questions/71087/analysis-of-a-time-series-with-a-fixed-and-random-factor-in-r
######################################################################
######################################################################
######################################################################
######################################################################

# a 3 year old has age_yr2 of -1.07
# a 9 year old has age_yr2 of 1.68
# a 10 year old has age_yr2 of 2.144
temp2$age_yr3<- temp2$age_yr - 3  # 3 year olds
temp2$age_yr2 <- rescale(temp2$age_yr3)


temp2$age_yr4<- temp2$age_yr2 + 1.07
t<-glmer(bruc~ age_yr3*tb+ I(age_yr2^2)+(1|id), data=temp2, family=binomial); summary(t)
t<-glmer(bruc~ age_yr4*tb+ I(age_yr2^2)+(1|id), data=temp2, family=binomial); summary(t)



get_tb_increase = function(ageyrval){
	logodds = 15.562 - ageyrval * 0.4134
	odds = exp(logodds)
	prop = odds / (1 + odds)
	return(list(odds, prop))	
}

tage<- seq(-2, 4, 0.1)
get_odds = function(age, tb){
	logodds = -25.19 + age * 9.54 + tb * 3.51 + -4.255 * tb * age - 0.538 * age * age 
	odds = exp(logodds)
	return(logodds)
}

age2<- tage* 2* sd(temp2$age_yr) + mean(temp2$age_yr)

plot(x= age2, y= get_odds(tage, 0))
points(x= age2, y= get_odds(tage, 1), pch=19, ylab= "log odds Br+")

tage<- seq(-1, 2, 0.1)
age2<- tage* 2* sd(temp2$age_yr) + mean(temp2$age_yr)

plot(x= age2, y= exp(get_odds(tage, 0)), ylab = "odds Br+")
points(x= age2, y= exp(get_odds(tage, 1)), pch=19)





##############################################################################

# All play with random effects... 
t <- glmer(bruc~ floor(age_yr)+ I(floor(age_yr^2))+(1|id), data=temp3, family=binomial(link="logit")); summary(t) #319.2
t<-glmer(bruc~ floor(age_yr)+(1|id), data=temp3, family=binomial); summary(t)
t<-glmer(bruc~ floor(age_yr)+ I(floor(age_yr)^2)+ tb+(1|id), data=temp3, family=binomial); summary(t) #349.4 # small error
t<-glmer(bruc~ floor(age_yr)+ I(floor(age_yr)^2)*tb+(1|id), data=temp3, family=binomial); summary(t) #
t<-glmer(bruc~age_yr2+ I(age_yr2^2)*tb+(1|id), data=temp3, family=binomial); summary(t) #

t<-glmer(bruc~ age_yr2*tb+ I(age_yr2^2)*tb+(1|id), data=temp3, family=binomial); summary(t) #

t<-glmer(bruc~ floor(age_yr)*tb+ I(floor(age_yr)^2)*tb+(1|id), data=temp3, family=binomial); summary(t) #
t<-glmer(bruc~ floor(age_yr)*tb+ I(floor(age_yr)^2)+(1|id), data=temp3, family=binomial); summary(t) #

t<-glmer(bruc~ floor(age_yr)*tb+ I(floor(age_yr)^2)+ herdorig+(1|id), data=temp3, family=binomial); summary(t) #351.2, small error
t<-glmer(bruc~ floor(age_yr)+ I(floor(age_yr)^2)*tb+ herdorig+(1|id), data=temp3, family=binomial); summary(t) #346.5
t<-glmer(bruc~ floor(age_yr)+ I(floor(age_yr)^2)*tb+ herd+(1|id), data=temp3, family=binomial); summary(t) #346.5

t<-glmer(bruc~ floor(age_yr)+ I(floor(age_yr)^2)*tb+tb*herd+(1|id), data=temp3, family=binomial); summary(t) # 354.9
t<-glmer(bruc~ floor(age_yr)+ I(floor(age_yr)^2)*tb+floor(age_yr)*herdorig+(1|id), data=temp3, family=binomial); summary(t) # 354.9

t<-glmer(bruc~ floor(age_yr2)+ I(floor(age_yr2)^2)*tb+(1|id), data=temp3, family=binomial); summary(t) 
t<-glmer(bruc~ floor(age_yr2)+ I(floor(age_yr2)^2)+tb+(1|id), data=temp3, family=binomial); summary(t) 
t<-glmer(bruc~ floor(age_yr2)*tb+ I(floor(age_yr2)^2)*tb+(1|id), data=temp3, family=binomial); summary(t)
t<-glmer(bruc~ floor(age_yr2)*tb+ herd+ I(floor(age_yr2)^2)*tb+(1|id), data=temp3, family=binomial); summary(t) 
t<-glmer(bruc~ floor(age_yr2)*tb+ herd+ I(floor(age_yr2)^2)+(1|id), data=temp3, family=binomial); summary(t) 

temp3$time<- NA
temp3$time[temp3$capturetime %in% c(0, 3, 6, 9)]<- 0
temp3$time[temp3$capturetime %in% c(12, 15, 18, 21)]<- 1
temp3$time[temp3$capturetime %in% c(24, 27, 30, 33)]<- 2
temp3$time[temp3$capturetime %in% c(36, 39, 42, 45, 48)]<- 3

# convergence issues
t<-glmer(bruc~ floor(age_yr2)+ I(floor(age_yr2)^2)+(time|id), data=temp3, family=binomial); summary(t) 
t<-glmer(bruc~ floor(age_yr2)+ I(floor(age_yr2)^2)*tb+(time|id), data=temp3, family=binomial); summary(t) 
t<-glmer(bruc~ floor(age_yr2)+ I(floor(age_yr2)^2)+tb+(time|id), data=temp3, family=binomial); summary(t) 
t<-glmer(bruc~ floor(age_yr2)*tb+ I(floor(age_yr2)^2)*tb+(time|id), data=temp3, family=binomial); summary(t) # conv error
t<-glmer(bruc~ floor(age_yr2)*tb+ herd+ I(floor(age_yr2)^2)*tb+(time|id), data=temp3, family=binomial); summary(t) 
t<-glmer(bruc~ floor(age_yr2)*tb+ herd+ I(floor(age_yr2)^2)+(time|id), data=temp3, family=binomial); summary(t) 

t<-glmer(bruc~ floor(age_yr2)+time*tb+ I(floor(age_yr2)^2)+(time|id), data=temp3, family=binomial); summary(t) 
t<-glmer(bruc~ floor(age_yr2)*time+ tb+ I(floor(age_yr2)^2)*time+(time|id), data=temp3, family=binomial); summary(t) 

t<-glmer(bruc~ floor(age_yr2)*time+ tb*time + floor(age_yr2)*tb + (time|id), data=temp3, family=binomial); summary(t) 
#t<-glmer(bruc~ floor(age_yr2)*capturetime+ tb*capturetime + floor(age_yr2)*tb + (capturetime |id), data=temp3, family=binomial); summary(t) 
t<-glmer(bruc~ floor(age_yr2)+ tb*time + floor(age_yr2)*tb + (time|id), data=temp3, family=binomial); summary(t) 
t<-glmer(bruc~ floor(age_yr2)+ tb*time + floor(age_yr2)*tb + (time|id), data=temp3, family=binomial); summary(t) 


t<-glmer(bruc~ floor(age_yr2)*time+ tb+ I(floor(age_yr2)^2)*time+(time|id), data=temp3, family=binomial); summary(t) 


# same here
t<-glmer(bruc~ floor(age_yr2)+ I(floor(age_yr2)^2)*tb+(capturetime|herd/id), data=temp3, family=binomial); summary(t) 
t<-glmer(bruc~ floor(age_yr2)+ I(floor(age_yr2)^2)+tb+(capturetime |herd/id), data=temp3, family=binomial); summary(t) 
t<-glmer(bruc~ floor(age_yr2)*tb+ I(floor(age_yr2)^2)*tb+(capturetime |herd/id), data=temp3, family=binomial); summary(t)

t<-glmer(bruc~ floor(age_yr2)+ I(floor(age_yr2)^2)*tb+(1|herd/id), data=temp3, family=binomial); summary(t) 
t<-glmer(bruc~ floor(age_yr2)*tb+ I(floor(age_yr2)^2)+(1|herd/id), data=temp3, family=binomial); summary(t) 
t<-glmer(bruc~ floor(age_yr2)+ I(floor(age_yr2)^2)+tb+(1 |herd/id), data=temp3, family=binomial); summary(t) 
t<-glmer(bruc~ floor(age_yr2)*tb+ I(floor(age_yr2)^2)*tb+(1 |herd/id), data=temp3, family=binomial); summary(t)


t<-glmer(bruc~ 1+(capturetime|herd/id), data=temp3, family=binomial); summary(t) 
t<-glmer(bruc~ 1+(capturetime|capturetime*herd/id), data=temp3, family=binomial); summary(t) 






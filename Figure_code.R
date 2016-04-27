#######################################################
#######################################################
# Figure 1
#######################################################
#######################################################
setwd("~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/final_datasets_copied_from_phdfolder")
library(ggplot2)
library('grid')
library('gridExtra') # specifies layout
# read in data prepared in cross_sectional_dataprep, groomed for my bTB statuses
data<-read.csv("cross_sectional_data_withdz_cleandisease_nofinal_Feb2016_capturetime_forsurv.csv")
data_nofinal<-data[data$final_capture=="0",] 
d<-data.frame(btb=data_nofinal$tb , bruc=as.character(data_nofinal$bruc), age=data_nofinal$age_sel/12, 
              id=data_nofinal$id)
d<-d[d$age<14,]
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
newdf2<-newdf[newdf$N>0,]; newdf<-newdf2

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
  theme(axis.title.x = element_text(size=16, vjust=-0.15),
        axis.title.y = element_text(size=16, vjust= 0.8),
        axis.text.x = element_text(size=14, vjust=-0.05),
        axis.text.y = element_text(size=14),
        panel.border = element_blank(), 
        axis.line = element_line(colour= "black"),
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
  theme_bw()+ xlab("Age (years)")+ ylab("Count")+
  theme(panel.border = element_blank(), 
        panel.margin = element_blank(),
        panel.grid.major= element_blank(),
        panel.grid.minor= element_blank(),
        axis.line = element_line(colour= "black"))

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
  theme(axis.title.x = element_text(size=16, vjust=-0.15),
        axis.title.y = element_text(size=16, vjust= 0.8),
        axis.text.x = element_text(size=14, vjust=-0.05),
        axis.text.y = element_text(size=14),
        panel.border = element_blank(), 
        axis.line = element_line(colour= "black"),
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


#######################################################
#######################################################
# Figure 2
#######################################################
#######################################################
immunecross<- read.csv("~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/final_datasets_copied_from_phdfolder/cross_sectional_data_withdz_cleandisease_nofinal_Feb2016_capturetime_forsurv.csv")
gcross<-immunecross[!is.na(immunecross$ifng),]
gcross<-gcross[gcross$ifng_CV<14,]
immuneba<- read.csv("~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/final_datasets_copied_from_phdfolder/convertersonly_Feb2016.csv")
g1<-immuneba[!is.na(immuneba$ifng),]
g1<-g1[g1$ifng_CV<14,]

g<-data.frame(ifng=g1$ifng, infection=paste(g1$bruc, g1$tb, sep="_"), herd=g1$herdorig)
gc<-data.frame(ifng=gcross$ifng, infection=paste(gcross$bruc, gcross$tb, sep="_"), herd=gcross$herdorig)

p<- ggplot(g, aes(x=infection, y=log(ifng), fill=infection))+ geom_boxplot()+ 
	scale_fill_discrete(labels=c("Neg", "bTB+", "Br+", "Co"))+
	ylab("L(IFNg) concentration (converters only)")+ theme(legend.position="top", 
	legend.title=element_blank())
p1<- ggplot(g[g$herd=="LS",], aes(x=infection, y=log(ifng), fill=infection))+ geom_boxplot()+
	scale_fill_discrete(labels=c("Neg", "bTB+", "Br+", "Co"))+
	ylab("L(IFNg) concentration (LS converters only)")+ theme(legend.position="top", 
	legend.title=element_blank())
p2<- ggplot(g[g$herd=="CB",], aes(x=infection, y=log(ifng), fill=infection))+ geom_boxplot()+
	scale_fill_discrete(labels=c("Neg", "bTB+", "Br+", "Co"))+
	ylab("L(IFNg) concentration (CB converters only)")+ theme(legend.position="top", 
	legend.title=element_blank())+scale_y_continuous(limits=c(0,2))
grid.arrange(p, p1, p2, ncol=3)

p3<- ggplot(gc, aes(x=infection, y=log(ifng), fill=infection))+ geom_boxplot()+ 
	scale_fill_discrete(labels=c("Neg", "bTB+", "Br+", "Co"))+
	ylab("L(IFN) gamma concentration (all)")+ theme(legend.position="top", 
	legend.title=element_blank())
p4<- ggplot(gc[gc$herd=="LS",], aes(x=infection, y=log(ifng), fill=infection))+ geom_boxplot()+ 
	scale_fill_discrete(labels=c("Neg", "bTB+", "Br+", "Co"))+
	ylab("L(IFN) gamma concentration (all LS)")+ theme(legend.position="top", 
	legend.title=element_blank())
p5<- ggplot(gc[gc$herd=="CB",], aes(x=infection, y=log(ifng), fill=infection))+ geom_boxplot()+ 
	scale_fill_discrete(labels=c("Neg", "bTB+", "Br+", "Co"))+
	ylab("L(IFN) gamma concentration (all CB)")+ theme(legend.position="top", 
	legend.title=element_blank())
grid.arrange(p3, p4, p5, ncol=3)

# with age
g<-data.frame(ifng=g1$ifng, infection=paste(g1$bruc, g1$tb, sep="_"), herd=g1$herdorig, age= g1$age_yr)
gc<-data.frame(ifng=gcross$ifng, infection=paste(gcross$bruc, gcross$tb, sep="_"), 
	herd=gcross$herdorig, age=gcross$age_yr)
	
a<-ggplot(g, aes(x=age, y=log(ifng), colour=infection))+ geom_point() +guides(fill=FALSE, colour=FALSE)+ ylab("IFN gamma (converters only)") + scale_colour_discrete(labels=c("Neg", "TB", "Bruc", "Coinfected")) 
b<-ggplot(gc, aes(x=age, y=log(ifng), colour=infection))+ geom_point()+ 	ylab("IFN gamma (all data)")+ scale_colour_discrete(labels=c("Neg", "TB", "Bruc", "Co"))+ theme(legend.position=c(0.8,0.8))
grid.arrange(a, b, ncol=2)


# Plasma BKA: 
gcross<-immunecross[!is.na(immunecross$bka_killed),]
g1<-immuneba[!is.na(immuneba$bka_killed),]
g<-data.frame(bka=g1$bka_killed/g1$bka_control, infection=paste(g1$bruc, g1$tb, sep="_"), herd=g1$herdorig)
gc<-data.frame(bka=gcross$bka_killed/gcross$bka_control, infection=paste(gcross$bruc, gcross$tb, sep="_"), herd=gcross$herdorig)
p<- ggplot(g, aes(x=infection, y= bka, fill=infection))+ geom_boxplot()+ 
	scale_fill_discrete(labels=c("bTB-, Br-", "bTB+, Br-", "bTB-, Br+", "bTB+, Br-"))+
	ylab("Proportion killed (converters only)")+ theme(legend.position="top", 
	legend.title=element_blank())
p1<- ggplot(g[g$herd=="LS",], aes(x=infection, y= bka, fill=infection))+ geom_boxplot()+
	scale_fill_discrete(labels=c("bTB-, Br-", "bTB+, Br-", "bTB-, Br+", "bTB+, Br-"))+
	ylab("Proportion killed (LS converters only)")+ theme(legend.position="top", 
	legend.title=element_blank())
p2<- ggplot(g[g$herd=="CB",], aes(x=infection, y= bka, fill=infection))+ geom_boxplot()+
	scale_fill_discrete(labels=c("bTB-, Br-", "bTB+, Br-", "bTB-, Br+", "bTB+, Br-"))+
	ylab("Proportion killed (CB converters only)")+ theme(legend.position="top", 
	legend.title=element_blank())+scale_y_continuous(limits=c(0,2))
grid.arrange(p, p1, p2, ncol=3)

p3<- ggplot(gc, aes(x=infection, y=bka, fill=infection))+ geom_boxplot()+ 
	scale_fill_discrete(labels=c("bTB-, Br-", "bTB+, Br-", "bTB-, Br+", "bTB+, Br-"))+
	ylab("Proportion killed (all)")+ theme(legend.position="top", 
	legend.title=element_blank())
p4<- ggplot(gc[gc$herd=="LS",], aes(x=infection, y=bka, fill=infection))+ geom_boxplot()+ 
	scale_fill_discrete(labels=c("bTB-, Br-", "bTB+, Br-", "bTB-, Br+", "bTB+, Br-"))+
	ylab("Proportion killed (all LS)")+ theme(legend.position="top", 
	legend.title=element_blank())
p5<- ggplot(gc[gc$herd=="CB",], aes(x=infection, y=bka, fill=infection))+ geom_boxplot()+ 
	scale_fill_discrete(labels=c("bTB-, Br-", "bTB+, Br-", "bTB-, Br+", "bTB+, Br-"))+
	ylab("Proportion killed (all CB)")+ theme(legend.position="top", 
	legend.title=element_blank())
grid.arrange(p3, p4, p5, ncol=3)

# Whole Blood BKA: 
gcross<-immunecross[!is.na(immunecross$bka_wb_killed),]
gc<-data.frame(bka=gcross$bka_wb_killed/gcross$bka_wb_control, infection=paste(gcross$bruc, gcross$tb, sep="_"), herd=gcross$herdorig)
p3<- ggplot(gc, aes(x=infection, y=bka, fill=infection))+ geom_boxplot()+ 
	scale_fill_discrete(labels=c("bTB-, Br-", "bTB+, Br-", "bTB-, Br+", "bTB+, Br-"))+
	ylab("Proportion killed, whole blood (all)")+ theme(legend.position="top", 
	legend.title=element_blank())
p4<- ggplot(gc[gc$herd=="LS",], aes(x=infection, y=bka, fill=infection))+ geom_boxplot()+ 
	scale_fill_discrete(labels=c("bTB-, Br-", "bTB+, Br-", "bTB-, Br+", "bTB+, Br-"))+
	ylab("Proportion killed, whole blood (all LS)")+ theme(legend.position="top", 
	legend.title=element_blank())
p5<- ggplot(gc[gc$herd=="CB",], aes(x=infection, y=bka, fill=infection))+ geom_boxplot()+ 
	scale_fill_discrete(labels=c("bTB-, Br-", "bTB+, Br-", "bTB-, Br+", "bTB+, Br-"))+
	ylab("Proportion killed, whole blood (all CB)")+ theme(legend.position="top", 
	legend.title=element_blank())
grid.arrange(p3, p4, p5, ncol=3)

# xyplots
gcross<-immunecross[!is.na(immunecross$ifng),]
PBC.samp <- subset(gcross, id %in% c("B1",   "B10",  "B11",  "B13",  "B13b", "B14",
	"B14b", "B16",  "B19", "B2",   "B20",  "B22",  "B22b", "B25",  "B26",  "B26b") )
	#"B28",  "B28b", "B29",  "B2b",  "B30",  "B31",  "B32",  "B33"))

xyplot(ifng ~ yr | id, data = PBC.samp,
    type = c("p", "smooth"), lwd = 2, layout = c(4, 4),
    as.table = TRUE, ylab = "log IFNgamma",
    xlab = "Time (years)")
xyplot(log(ifng) ~ yr | id, data = PBC.samp,
    type = c("p", "smooth"), lwd = 2, layout = c(4, 4),
    as.table = TRUE, ylab = "log IFNgamma",
    xlab = "Time (years)")


# Haptoglobin: 
gcross<-immunecross[!is.na(immunecross$hapto),]
gc<-data.frame(hapto=log(gcross$hapto), infection=paste(gcross$bruc, gcross$tb, sep="_"), herd=gcross$herdorig)
p3<- ggplot(gc, aes(x=infection, y=hapto, fill=infection))+ geom_boxplot()+ 
	scale_fill_discrete(labels=c("Neg", "bTB+", "Br+", "Co"))+
	ylab("Log(Haptoglobin) (all)")+ theme(legend.position="top", 
	legend.title=element_blank())
p4<- ggplot(gc[gc$herd=="LS",], aes(x=infection, y=hapto, fill=infection))+ geom_boxplot()+ 
	scale_fill_discrete(labels=c("Neg", "bTB+", "Br+", "Co"))+
	ylab("Log(Haptoglobin) (all LS)")+ theme(legend.position="top", 
	legend.title=element_blank())
p5<- ggplot(gc[gc$herd=="CB",], aes(x=infection, y=hapto, fill=infection))+ geom_boxplot()+ 
	scale_fill_discrete(labels=c("Neg", "bTB+", "Br+", "Co"))+
	ylab("Log(Haptoglobin) (all CB)")+ theme(legend.position="top", 
	legend.title=element_blank())
grid.arrange(p3, p4, p5, ncol=3)




#######################################################
#######################################################
# Figure 3- Survival and incidence plots
#######################################################
#######################################################
# From book code: 
#AIDS.samp <- subset(aids, patient %in% c(82,152,213,236,332,
     335,353,407,410,452))
#KM <- survfit(Surv(Time, death) ~ 1, data = aids.id)

#par(mfrow = c(1, 2))
#plot(KM, mark.time = FALSE, ylab = "Survival Probability",
    xlab = "Time (months)")
 # USE STANDARDIZED ESTIMATES!!!
df<- data.frame(name=c("bruc", "bTB", "site", "age^2", "age"), est= c(3.0, 2.81, 1.99, 1.05,  0.51), 
	lower= c(1.50, 1.43, 1.05, 1.01 , 0.34), upper= c(6.01, 5.51, 3.78, 1.07, 0.76),  #78.83 
	order=c(1,2,3,4,5))
df$name <- factor(df$name, levels= df$name[order(df$order, decreasing = TRUE)])
	
p2 <-ggplot(df, aes(x=df$name, y=df$est)) + 
	geom_errorbar(aes(ymin=df$lower, ymax=df$upper, width=0)) +  				# for black error bars, add colour="black"
	geom_point(df$estimate)+
	theme_bw()+
	#theme(panel.grid.major.y=element_blank(), 
	#panel.border = element_blank(), 
 	#axis.line = element_line(colour= "black"),
 	#panel.grid.major.x=element_blank()) +  
	coord_flip()
	#ylim(-0.5, 0.5)+
	#ylim(min(out$sdbeta-out$sdse), max(out$sdbeta[out$sdbeta<100]+out$sdse[out$sdbeta<100]))   # change to this
	#theme(axis.text.x=element_text(size=14), axis.text.y=element_text(size=14))

p3<- p2 + 
	theme(axis.title.x=element_text(size=14), axis.title.y=element_blank()) +
	ylab("Relative risk compared to uninfected animals")+
	scale_colour_grey(start=0.6, end=0.7)+
	geom_segment(aes(x=0, xend=5, y=1, yend=1), linetype=2)+ ylim(-0.01,6.2) 

p3

p3 + theme(axis.line = element_line(colour= "black"))






#############################################
#############################################
# Notes of GLMM
#############################################
#############################################
# AICs hanshed out ran on previous dataset.  IN repot results are presented for the dataset with capture time that matches the survival times.
temp<-data_nofinal[data_nofinal$age_yr<14,]
rescale= function(col){
  new=NA
  for (i in 1:length(col)){
    new[i]<-(col[i]-mean(col))/sd(col)
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



t<-glmer(bruc~ floor(age_yr)+ I(floor(age_yr^2))+(1|id), data=temp2, family=binomial(link="logit")); summary(t) #319.2
t<-glmer(bruc~ floor(age_yr)+(1|id), data=temp2, family=binomial); summary(t)
t<-glmer(bruc~ floor(age_yr)+ I(floor(age_yr)^2)+ tb+(1|id), data=temp2, family=binomial); summary(t) #349.4 # small error
t<-glmer(bruc~ floor(age_yr)+ I(floor(age_yr)^2)*tb+(1|id), data=temp2, family=binomial); summary(t) #345.2
t<-glmer(bruc~ floor(age_yr)*tb+ I(floor(age_yr)^2)+ herdorig+(1|id), data=temp2, family=binomial); summary(t) #351.2, small error
t<-glmer(bruc~ floor(age_yr)+ I(floor(age_yr)^2)*tb+ herdorig+(1|id), data=temp2, family=binomial); summary(t) #346.5
t<-glmer(bruc~ floor(age_yr)+ I(floor(age_yr)^2)*tb+tb*herdorig+(1|id), data=temp2, family=binomial); summary(t) # 354.9
t<-glmer(bruc~ floor(age_yr)+ I(floor(age_yr)^2)*tb+floor(age_yr)*herdorig+(1|id), data=temp2, family=binomial); summary(t) # 354.9


t<-glmer(bruc~age_yr2+ I(age_yr2^2)*tb2+(1|id), data=temp2, family=binomial); summary(t) # 354.9




t<-glmer(bruc~ floor(age_yr)+ I(floor(age_yr^2))+(1|id), data=temp3, family=binomial(link="logit")); summary(t) #319.2
t<-glmer(bruc~ floor(age_yr)+(1|id), data=temp3, family=binomial); summary(t)
t<-glmer(bruc~ floor(age_yr)+ I(floor(age_yr)^2)+ tb+(1|id), data=temp3, family=binomial); summary(t) #349.4 # small error
t<-glmer(bruc~ floor(age_yr)+ I(floor(age_yr)^2)*tb+(1|id), data=temp3, family=binomial); summary(t) #345.2
t<-glmer(bruc~ floor(age_yr)*tb+ I(floor(age_yr)^2)+ herdorig+(1|id), data=temp3, family=binomial); summary(t) #351.2, small error
t<-glmer(bruc~ floor(age_yr)+ I(floor(age_yr)^2)*tb+ herdorig+(1|id), data=temp3, family=binomial); summary(t) #346.5
t<-glmer(bruc~ floor(age_yr)+ I(floor(age_yr)^2)*tb+ herd+(1|id), data=temp3, family=binomial); summary(t) #346.5

t<-glmer(bruc~ floor(age_yr)+ I(floor(age_yr)^2)*tb+tb*herd+(1|id), data=temp3, family=binomial); summary(t) # 354.9
t<-glmer(bruc~ floor(age_yr)+ I(floor(age_yr)^2)*tb+floor(age_yr)*herdorig+(1|id), data=temp3, family=binomial); summary(t) # 354.9















# Extra play when choosing random effects

# Use herd as a random effect
# not fitting
t0<-glmer(bruc~ tb+ (1|herdorig/id), data=temp2, family=binomial); summary(t0) # 481.0
t1<-glmer(bruc~ age_yr+ (1|herdorig/id), data=temp2, family=binomial); summary(t1) # 317.3 error?
t2<-glmer(bruc~ log(age_yr)+ (1|herdorig/id), data=temp2, family=binomial); summary(t2) # 320.2
t3<-glmer(bruc~ log(age_yr)*tb+ (1|herdorig/id), data=temp2, family=binomial); summary(t3) # 319.8 fit error!
t4<-glmer(bruc~ age_yr*tb+ (1|herdorig/id), data=temp2, family=binomial); summary(t4) # 322.9 # fit error!
t5<-glmer(bruc~ log(age_yr)+tb+ (1|herdorig/id), data=temp2, family=binomial); summary(t5) # 321.7 n.s.
t6<-glmer(bruc~ age_yr+tb+ (1|herdorig/id), data=temp2, family=binomial); summary(t6) # 328.4 n.s.
a<-sd(temp$age_yr)
t3<-glmer(bruc~ log(age_yr)*tb+ herdorig+ (1|id), data=temp2, family=binomial); summary(t3) 

# with herd as a main effect
t0<-glmer(bruc~ tb+ (1|id), data=temp2, family=binomial); summary(t0) # 479.0
t0.5<-glmer(bruc~ tb+ herdorig+ (1|id), data=temp2, family=binomial); summary(t0.5) # 480.8
t1<-glmer(bruc~ age_yr+ (1|id), data=temp2, family=binomial); summary(t1) # 317.6
t1.5<-glmer(bruc~ age_yr+ herdorig+ (1|id), data=temp2, family=binomial); summary(t1.5) #314.8
t2<-glmer(bruc~ log(age_yr)+ (1|id), data=temp2, family=binomial); summary(t2) # 310.0
t2.5<-glmer(bruc~ log(age_yr)+ herdorig+(1|id), data=temp2, family=binomial); summary(t2.5) # 311.7
t3<-glmer(bruc~ log(age_yr)*tb+ (1|id), data=temp2, family=binomial); summary(t3) # 309.7 !!!!!!!
t3.5<-glmer(bruc~ log(age_yr)*tb+herdorig+ (1|id), data=temp2, family=binomial); summary(t3.5) # fit error!
t4<-glmer(bruc~ age_yr*tb+ (1|id), data=temp2, family=binomial); summary(t4) # 314.3
t4.5<-glmer(bruc~ age_yr*tb+herdorig+(1|id), data=temp2, family=binomial); summary(t4.5) # fit error!
t5<-glmer(bruc~ log(age_yr)+tb+ (1|id), data=temp2, family=binomial); summary(t5) # 311.3
t5.5<-glmer(bruc~ log(age_yr)+tb+herdorig+ (1|id), data=temp2, family=binomial); summary(t5.5) # 312.6
t6<-glmer(bruc~ age_yr+tb+ (1|id), data=temp2, family=binomial); summary(t6) #319.0
t6.5<-glmer(bruc~ age_yr+tb+ herdorig+(1|id), data=temp2, family=binomial); summary(t6.5) #316.8

t7<-glmer(bruc~ tb+(1|id), data=temp2, family=binomial); summary(t7) #480.8
t7<-glmer(bruc~ tb+ herdorig+(1|id), data=temp2, family=binomial); summary(t7) #480.8
t8<-glmer(bruc~ tb*herdorig+(1|id), data=temp2, family=binomial); summary(t8) # 481.3
t9<-glmer(bruc~ age_yr+ herdorig+(1|id), data=temp2, family=binomial); summary(t9) #314.8
t10<-glmer(bruc~ age_yr*herdorig+(1|id), data=temp2, family=binomial); summary(t10) #313.2
t11<-glmer(bruc~ log(age_yr)+ herdorig+(1|id), data=temp2, family=binomial); summary(t11) #311.7
t12<-glmer(bruc~ log(age_yr)*herdorig+(1|id), data=temp2, family=binomial); summary(t12) # 313.2convergence error


t13<-glmer(bruc~ log(age_yr)*tb+herdorig*log(age_yr)+ 
             (1|id), data=temp2, family=binomial, REML=TRUE); summary(t13) # 309.5, fit error!
t14<-glmer(bruc~ age_yr*tb+herdorig*age_yr+ (1|id), data=temp2, family=binomial); summary(t14) # fit error!


t<-glmer(bruc~ age_yr*tb+tb*I(age_yr^2)+herdorig+(1|id), data=temp2, family=binomial); summary(t) 
t<-glmer(bruc~ age_yr*tb+I(age_yr^2)+herdorig+(1|id), data=temp2, family=binomial); summary(t) 
t<-glmer(bruc~ age_yr+tb+I(age_yr^2)+herdorig+(1|id), data=temp2, family=binomial); summary(t) 

# compare age curves
t<-glmer(bruc~ age_yr+ I(age_yr^2)+(1|id), data=temp2, family=binomial(link="logit")); summary(t) #298.5, error
t<-glmer(bruc~ age_yr+(1|id), data=temp2, family=binomial); summary(t) #305.8
t<-glmer(bruc~ log(age_yr)+(1|id), data=temp2, family=binomial); summary(t) #300.4

t<-glmer(bruc~ floor(age_yr)+ I(floor(age_yr^2))+(1|id), data=temp2, family=binomial(link="logit")); summary(t) #319.2
t<-glmer(bruc~ floor(age_yr)+(1|id), data=temp2, family=binomial); summary(t) # 357.7
t<-glmer(bruc~ log(floor(age_yr))+(1|id), data=temp2, family=binomial); summary(t) #350.6

t<-glmer(bruc~ floor(age_yr)+ I(floor(age_yr)^2)+ tb+(1|id), data=temp2, family=binomial); summary(t) #349.4 # small error
t<-glmer(bruc~ floor(age_yr)+ I(floor(age_yr)^2)*tb+(1|id), data=temp2, family=binomial); summary(t) #345.2

t<-glmer(bruc~ floor(age_yr)*tb+ I(floor(age_yr)^2)+ herdorig+(1|id), data=temp2, family=binomial); summary(t) #351.2, small error
t<-glmer(bruc~ floor(age_yr)+ I(floor(age_yr)^2)*tb+ herdorig+(1|id), data=temp2, family=binomial); summary(t) #346.5
t<-glmer(bruc~ floor(age_yr)+ I(floor(age_yr)^2)*tb+tb*herdorig+(1|id), data=temp2, family=binomial); summary(t) # 354.9
t<-glmer(bruc~ floor(age_yr)+ I(floor(age_yr)^2)*tb+tb*herdorig+(1|id), data=temp2, family=binomial); summary(t) #354.9
t<-glmer(bruc~ floor(age_yr)*tb+ I(floor(age_yr)^2)*tb+tb*herdorig+(1|id), data=temp2, family=binomial); summary(t) #377.4








# NOPE: 
t<-glmer(bruc~ floor(age_yr)+ I(floor(age_yr)^2)+ tb+(1|herdorig/id), data=temp2, family=binomial); summary(t) #351.4 
t<-glmer(bruc~ floor(age_yr)+ I(floor(age_yr)^2)*tb+(1|herdorig/id), data=temp2, family=binomial); summary(t) #345.2
t<-glmer(bruc~ floor(age_yr)*tb+ I(floor(age_yr)^2)*tb+(1|herdorig/id), data=temp2, family=binomial); summary(t) #345.2

a<-seq(1, 20, 0.1)
val=function(a, btb){14.954*a-0.5144*(a^2)-0.2772*(a^2)*btb + 9.592*btb }

par(mfrow= c(1,2))
plot(x=a, y=val(a, 1))
plot(x=a, y=val(a, 0))


# max at y[78]
a[78] = 8.7 in TB+




d<-data.frame(btb=temp2$tb , bruc=as.character(temp2$bruc), age=temp2$age_yr, 
              herd=temp2$herdorig, id=temp2$id)
make_age_odds_plot(d$btb, d$bruc, d$age, binsize=2)
make_age_odds_plot(d$btb[d$herd=="LS"], d$bruc[d$herd=="LS"], d$age[d$herd=="LS"], binsize=2)
make_age_odds_plot(d$btb[d$herd=="CB"], d$bruc[d$herd=="CB"], d$age[d$herd=="CB"], binsize=2)





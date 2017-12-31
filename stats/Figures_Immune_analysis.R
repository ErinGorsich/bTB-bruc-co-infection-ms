##################################################
##################################################
# Immune data analysis, cut
##################################################
##################################################
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




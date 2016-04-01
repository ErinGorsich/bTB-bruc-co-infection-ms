#######################################################
#######################################################
# Figure 1
#######################################################
#######################################################
setwd("~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/final_datasets_copied_from_phdfolder")
library(ggplot2)
library('grid')
# read in data prepared in cross_sectional_dataprep, groomed for my bTB statuses
data<-read.csv("cross_sectional_data_withdz_cleandisease_withfinal_Feb2016.csv")
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
agebins=c(seq(0, max(age), binsize), max(age+1))
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
#############################################
# Notes of GLMM
#############################################
#############################################
data<-read.csv("cross_sectional_data_withdz_cleandisease_withfinal_Feb2016.csv")
data_nofinal<-data[data$final_capture=="0",] 
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
val=function(a, btb){16.6487*a-0.48*(a^2)-0.48*(a^2)*btb}

plot(x=a, y=val(a, 1))
plot(x=a, y=val(a, 0))


# max at y[78]
a[78] = 8.7 in TB+




d<-data.frame(btb=temp2$tb , bruc=as.character(temp2$bruc), age=temp2$age_yr, 
              herd=temp2$herdorig, id=temp2$id)
make_age_odds_plot(d$btb, d$bruc, d$age, binsize=2)
make_age_odds_plot(d$btb[d$herd=="LS"], d$bruc[d$herd=="LS"], d$age[d$herd=="LS"], binsize=2)
make_age_odds_plot(d$btb[d$herd=="CB"], d$bruc[d$herd=="CB"], d$age[d$herd=="CB"], binsize=2)



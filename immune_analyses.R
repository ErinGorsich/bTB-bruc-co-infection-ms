#########################################
#########################################
# Immunological analyses
#########################################
#########################################
library(ggplots2)
library(lme4)
library(pbkrtest)
library("JM")
library("lattice")

############################################################
############################################################
############################################################
# IFNgamma- Cross sectional full data
############################################################
############################################################
############################################################
immune<-read.csv("~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/final_datasets_copied_from_phdfolder/cross_sectional_data_withdz_cleandisease_nofinal_Feb2016_capturetime_forsurv.csv")  # 5 animals added from last time
# added these five to cross_sectional_data_withdz_cleandisease_nofinal_Feb2016_capturetime.csv because in survival analyses.
immune2<-immune[!is.na(immune$ifng),]
immune3<-immune2[!is.na(immune2$ifng_plate),]  # 1 NA plate value!
immune3<-immune3[!is.na(immune3$ifng_delay),]
immune3<-immune3[immune3$ifng_CV<15,]
immune3$ifng_plate<-as.factor(immune3$ifng_plate)

# who is missing in the gamma data and why?
surv<-read.csv("~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/final_datasets_copied_from_phdfolder/survival/brucsurvival_TB3controls_longresidnomissing_noerrors_season2.csv")
unique(surv$animal[!(surv$animal %in% immune2$id)]) # 14 bolused animals in survival anlayses that need sorted. 


# Step 1- decide on random effects
#########################################
# note to self, order of listing in nlme implies nesting, wich we do not want
###http://stats.stackexchange.com/questions/58669/specifying-multiple-separate-random-effects-in-lme
###https://biostatmatt.com/archives/2718
# LME specification required for JM in joint_play code

test<-lmer(log(ifng)~1+ (1|ifng_plate) + (capturetime|id), data=immune3)
test<-lmer(log(ifng)~capturetime + (1|ifng_plate) + (capturetime|id), data=immune3)

red1<-test<-lmer(log(ifng)~1 + (capturetime|id), data=immune3)
anova(test, red1)
red2<-lmer(log(ifng)~1+ (1|ifng_plate) + (1|id), data=immune3) # crossed design is better
red3<-lmer(log(ifng)~1+ (1|ifng_plate/id), data=immune3) # nested

plot(resid(red2)~immune3$capturetime)


# Step 2- model selection
#########################################
full.mod<-lmer(log(ifng)~bruc+herdorig+tb+age_yr+I(age_yr^2)+ 
 	tb*age_yr+ tb*I(age_yr^2)+ tb*herdorig+ bruc*age_yr+ bruc*I(age_yr^2)+ 
 	bruc*herdorig+ tb*bruc+ifng_delay + (1|ifng_plate) + (capturetime|id), data=immune3)
 	
immune3$plate<-as.factor(immune3$ifng_plate)
full.mod<-lmer(log(ifng)~bruc+herdorig+tb+age_yr+I(age_yr^2)+ 
  	tb*age_yr+ tb*I(age_yr^2)+ tb*herdorig+ bruc*age_yr+ bruc*I(age_yr^2)+ 
  	bruc*herdorig+ tb*bruc+ifng_delay + (1|plate) + (capturetime|id), data=immune3)
lwm<-update(full.mod, .~. -tb:I(age_yr^2)); anova(full.mod, lwm) # close second
lwm<-update(full.mod, .~. -bruc:I(age_yr^2)); anova(full.mod, lwm)
lwm<-update(full.mod, .~. -tb:herdorig); anova(full.mod, lwm)  # highest p
lwm<-update(full.mod, .~. -bruc:herdorig); anova(full.mod, lwm)
lwm<-update(full.mod, .~. -bruc:age_yr); anova(full.mod, lwm)
lwm<-update(full.mod, .~. -bruc:age_yr); anova(full.mod, lwm)
lwm<-update(full.mod, .~. -tb:bruc); anova(full.mod, lwm)

full.mod<-update(full.mod, .~. -tb:herdorig)
lwm<-update(full.mod, .~. -tb:I(age_yr^2)); anova(full.mod, lwm) # highest p
lwm<-update(full.mod, .~. -bruc:I(age_yr^2)); anova(full.mod, lwm)
lwm<-update(full.mod, .~. -bruc:herdorig); anova(full.mod, lwm)
lwm<-update(full.mod, .~. -bruc:age_yr); anova(full.mod, lwm)
lwm<-update(full.mod, .~. -bruc:age_yr); anova(full.mod, lwm)
lwm<-update(full.mod, .~. -tb:bruc); anova(full.mod, lwm)

full.mod<-update(full.mod, .~. -tb:I(age_yr^2))
lwm<-update(full.mod, .~. -bruc:I(age_yr^2)); anova(full.mod, lwm)
lwm<-update(full.mod, .~. -bruc:herdorig); anova(full.mod, lwm)
lwm<-update(full.mod, .~. -tb:age_yr); anova(full.mod, lwm) # highest p.
lwm<-update(full.mod, .~. -bruc:age_yr); anova(full.mod, lwm)
lwm<-update(full.mod, .~. -tb:bruc); anova(full.mod, lwm)

full.mod<-update(full.mod, .~. -tb:age_yr)
lwm<-update(full.mod, .~. -bruc:I(age_yr^2)); anova(full.mod, lwm) # drop
lwm<-update(full.mod, .~. -bruc:herdorig); anova(full.mod, lwm)
lwm<-update(full.mod, .~. -bruc:age_yr); anova(full.mod, lwm)
lwm<-update(full.mod, .~. -tb:bruc); anova(full.mod, lwm)

full.mod<-update(full.mod, .~. -bruc:I(age_yr^2))
lwm<-update(full.mod, .~. -bruc:herdorig); anova(full.mod, lwm)
lwm<-update(full.mod, .~. -bruc:age_yr); anova(full.mod, lwm)
lwm<-update(full.mod, .~. -tb:bruc); anova(full.mod, lwm)

full.mod<-update(full.mod, .~. -bruc:age_yr)
lwm<-update(full.mod, .~. -bruc:herdorig); anova(full.mod, lwm)
lwm<-update(full.mod, .~. -tb:bruc); anova(full.mod, lwm)


full.mod<-update(full.mod, .~. -tb:bruc)
lwm<-update(full.mod, .~. -age_yr); anova(full.mod, lwm)
lwm<-update(full.mod, .~. -I(age_yr^2)); anova(full.mod, lwm)
lwm<-update(full.mod, .~. -ifng_delay); anova(full.mod, lwm)
lwm<-update(full.mod, .~. -tb); anova(full.mod, lwm)

full.mod<-update(full.mod, .~. -ifng_delay)
lwm<-update(full.mod, .~. -age_yr); anova(full.mod, lwm)
lwm<-update(full.mod, .~. -I(age_yr^2)); anova(full.mod, lwm)
lwm<-update(full.mod, .~. -bruc:herdorig); anova(full.mod, lwm)
lwm<-update(full.mod, .~. -tb); anova(full.mod, lwm)

full.mod<-update(full.mod, .~. -I(age_yr^2))
lwm<-update(full.mod, .~. -age_yr); anova(full.mod, lwm)
full.mod<-update(full.mod, .~. -age_yr)

immune3$herd<-relevel(immune3$herdorig, "LS")
final.mod<-lmer(log(ifng) ~ bruc + herd + tb + bruc:herd + (1 | ifng_plate) + (capturetime | id), data=immune3 )

# 3. Use the Kenward-Roger approximation to get approximate degrees of freedom and the t-distribution to get p-values, which is implemented in the pbkrtest package.
#http://mindingthebrain.blogspot.com/2014/02/three-ways-to-get-parameter-specific-p.html
# get p-values from the t-distribution using the t-values and approximated degrees of freedom
coefs <- data.frame(coef(summary(full.mod)))
df.KR <- get_ddf_Lb(full.mod, fixef(full.mod))
coefs$p.KR <- 2 * (1 - pt(abs(coefs$t.value), df.KR))
coefs


#Scaled residuals: 
    Min      1Q  Median      3Q     Max 
-5.5140 -0.5496  0.0918  0.5930  2.5391 

#Random effects:
# Groups   Name        Variance  Std.Dev. Corr 
# id       (Intercept) 2.038e-01 0.45145       
#          capturetime 5.084e-05 0.00713  -0.54
# plate    (Intercept) 5.134e-02 0.22659       
# Residual             3.129e-01 0.55938       
#Number of obs: 744, groups:  id, 146; plate, 42

#Fixed effects:
#                    Estimate Std. Error t value
#(Intercept)         -0.67489    0.08304  -8.128
#brucpositive         0.14937    0.09470   1.577
#herdCB               0.18303    0.11004   1.663
#tb                   0.19961    0.06995   2.853
#brucpositive:herdCB -0.39426    0.14011  -2.814

#Correlation of Fixed Effects:
#            (Intr) brcpst herdCB tb    
#brucpositiv -0.406                     
#herdCB      -0.682  0.327              
#tb          -0.208 -0.079  0.001       
#brcpstv:hCB  0.283 -0.670 -0.438 -0.001

# p-values
coefs <- data.frame(coef(summary(final.mod)))
df.KR <- get_ddf_Lb(lwm, fixef(final.mod))
coefs$p.KR <- 2 * (1 - pt(abs(coefs$t.value), df.KR))
coefs
#                      Estimate Std..Error   t.value         p.KR
#(Intercept)         -0.6667788 0.08281638 -8.051291 2.771117e-13
#brucpositive         0.1640323 0.09337660  1.756675 8.109329e-02
#herdCB               0.1679881 0.10868023  1.545710 1.243637e-01
#tb                   0.1795050 0.06709383  2.675432 8.325425e-03
#brucpositive:herdCB -0.4081314 0.13545922 -3.012947 3.056206e-03


#coefs
#                          Estimate Std..Error   t.value         p.KR
#(Intercept)             -0.4987900 0.07922996 -6.295472 3.467603e-09
#brucpositive            -0.2441006 0.09909454 -2.463310 1.494133e-02
#herdorigLS              -0.1679882 0.10868074 -1.545704 1.243653e-01
#tb                       0.1795036 0.06709389  2.675408 8.325984e-03
#brucpositive:herdorigLS  0.4081334 0.13545962  3.012952 3.056152e-03


immune3$TBcol[immune3$convert==1]<- "red"
immune3$TBcol[immune3$convert==0]<- "blue"
immune3$TBcol[immune3$tb_beforeafter=="pfc"]<- "dark green"

immune3$Brlty[immune3$brucconvert ==1]<- 1
immune3$Brlty[immune3$brucconvert ==0]<- 2
immune3$Brlty[immune3$bruc_beforeafter =="pfc"]<- 3

immune3$Brpch[immune3$brucconvert ==1]<- 15
immune3$Brpch[immune3$brucconvert ==0]<- 17
immune3$Brpch[immune3$bruc_beforeafter =="pfc"]<- 19


plot(log(ifng) ~ capturetime, data = immune3, type = "n",
    xlab = "Time (months)", ylab = expression(log("ifng concentration")))
for (i in unique(immune3$id)[1:50])
    points(log(ifng) ~ jitter(capturetime), data = immune3[immune3$id == i, ],
        pch = immune3$Brpch[immune3$id==i], col= immune3$TBcol[immune3$id==i])
for (i in unique(immune3$id)[1:50])
    lines(log(ifng) ~ jitter(capturetime), data = immune3[immune3$id == i, ],
        lty = immune3$Brlty[immune3$id==i], col= immune3$TBcol[immune3$id==i])
legend("bottomright", legend=c("TB non-converter", "TB, converter", "TB pfc"), lty= c(1,1,1), col=c("blue", "red","dark green"), bty="n")
legend("bottomleft", legend=c("Bruc, converter", "Bruc non-converter", "Bruc pfc"), pch= c(17,15,19), bty="n")


plot(log(ifng) ~ capturetime, data = immune3, type = "n",
    xlab = "Time (months)", ylab = expression(log("ifng concentration")))
for (i in unique(immune3$id))
    points(log(ifng) ~ jitter(capturetime), data = immune3[immune3$id == i, ],
        pch = immune3$Brpch[immune3$id==i], col= immune3$TBcol[immune3$id==i], cex=0.7)
for (i in unique(immune3$id))
    lines(log(ifng) ~ jitter(capturetime), data = immune3[immune3$id == i, ],
        lty = immune3$Brlty[immune3$id==i], col= immune3$TBcol[immune3$id==i])
legend("bottomright", legend=c("TB non-converter", "TB, converter", "TB pfc"), lty= c(1,1,1), col=c("red", "blue", "dark green"), bty="n")
legend("bottomleft", legend=c("Bruc non-converter", "Bruc, converter", "Bruc pfc"), pch= c(17,15,19), bty="n")





############################################################
# IFNgamma- Before vs. After infection
############################################################
# Time since bTB in only those that became infected with both:
temp<- read.csv("~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/final_datasets_copied_from_phdfolder/time_since_tb.csv")
ggplot(data=temp, aes(x=floor(tb_time), y=ifng, group=id, colour=bruc))+ geom_line()+ xlab("Time since TB")

# Time since bTB in all buffalo
all<- read.csv("~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/final_datasets_copied_from_phdfolder/time_since_tb_alltbconverters.csv")  # added time since tb by hand. 

ggplot(data=all, aes(x=floor(tsi_tb), y=ifng, group=id, colour=bruc))+ geom_line()+ 
	xlab("Time since bTB seroconversion") + ylab("IFNg concentration")+ theme_bw() +
	theme(axis.title.x= element_text(size=16), axis.title.y = element_text(size=16), 
	axis.text.x = element_text(size=14, vjust=-0.05), axis.text.y = element_text(size=14),
	axis.line= element_line(colour="black"), 
	legend.position=c(0.82, 0.9), 
	legend.background= element_rect(fill="white", colour="white"), 
	legend.key= element_blank(),
    legend.title= element_blank(),
    legend.text = element_text(size=15) ) + guides(fill=guide_legend(title="Brucellosis"))


# Model selection
data<-read.csv("~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/final_datasets_copied_from_phdfolder/convertersonly_Feb2016.csv")






############################################################
############################################################
############################################################
# Haptoglobin- Cross sectional full data
############################################################
############################################################
############################################################
immune<-read.csv("~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/final_datasets_copied_from_phdfolder/cross_sectional_data_withdz_cleandisease_nofinal_Feb2016_capturetime_forsurv.csv")  # 5 animals added from last time
# added these five to cross_sectional_data_withdz_cleandisease_nofinal_Feb2016_capturetime.csv because in survival analyses.
immune2<-immune[!is.na(immune$hapto),]
immune3<-immune2[!is.na(immune2$hapto_plate),]  # 1 NA plate value!
immune3$hapto_plate<-as.factor(immune3$hapto_plate)

table(immune3$id, immune3$hapto_plate)

# Step 1- decide on random effects
#########################################
# note to self, order of listing in nlme implies nesting, wich we do not want
###http://stats.stackexchange.com/questions/58669/specifying-multiple-separate-random-effects-in-lme
###https://biostatmatt.com/archives/2718
# LME specification required for JM in joint_play code

test<-lmer(log(hapto)~1+ (1|hapto_plate) + (capturetime|id), data=immune3)
test<-lmer(log(hapto)~capturetime + (1|ifng_plate) + (capturetime|id), data=immune3)

red1<-lmer(log(hapto)~1 + (capturetime|id), data=immune3)
anova(test, red1)
red2<-lmer(log(hapto)~1+ (1| hapto_plate) + (1|id), data=immune3) # crossed design is better
anova(test, red2) # time not necessary
red3<-lmer(log(hapto)~1+ (1| hapto_plate/id), data=immune3) # nested

plot(resid(red2)~immune3$capturetime)



# Step 2- model selection
#########################################
full.mod<-lmer(log(hapto)~bruc+herdorig+tb+age_yr+I(age_yr^2)+ 
 	tb*age_yr+ tb*I(age_yr^2)+ tb*herdorig+ bruc*age_yr+ bruc*I(age_yr^2)+ 
 	bruc*herdorig+ tb*bruc+bka_delay + (1|hapto_plate) + (capturetime|id), data=immune3)
 	
lwm<-update(full.mod, .~. -tb:I(age_yr^2)); anova(full.mod, lwm) 
lwm<-update(full.mod, .~. -bruc:I(age_yr^2)); anova(full.mod, lwm)
lwm<-update(full.mod, .~. -tb:herdorig); anova(full.mod, lwm)  
lwm<-update(full.mod, .~. -bruc:herdorig); anova(full.mod, lwm)
lwm<-update(full.mod, .~. -bruc:age_yr); anova(full.mod, lwm)
lwm<-update(full.mod, .~. -tb:age_yr); anova(full.mod, lwm)
lwm<-update(full.mod, .~. -tb:bruc); anova(full.mod, lwm)

full.mod<-update(full.mod, .~. -tb:I(age_yr^2))
lwm<-update(full.mod, .~. -tb:herdorig); anova(full.mod, lwm) 
lwm<-update(full.mod, .~. -bruc:I(age_yr^2)); anova(full.mod, lwm)
lwm<-update(full.mod, .~. -bruc:herdorig); anova(full.mod, lwm)
lwm<-update(full.mod, .~. -bruc:age_yr); anova(full.mod, lwm)
lwm<-update(full.mod, .~. -tb:age_yr); anova(full.mod, lwm)
lwm<-update(full.mod, .~. -tb:bruc); anova(full.mod, lwm)

full.mod<-update(full.mod, .~. -tb:bruc)
lwm<-update(full.mod, .~. -tb:herdorig); anova(full.mod, lwm) 
lwm<-update(full.mod, .~. -bruc:I(age_yr^2)); anova(full.mod, lwm)
lwm<-update(full.mod, .~. -bruc:herdorig); anova(full.mod, lwm)
lwm<-update(full.mod, .~. -bruc:age_yr); anova(full.mod, lwm)
lwm<-update(full.mod, .~. -tb:age_yr); anova(full.mod, lwm)

full.mod <-update(full.mod, .~. -bruc:herdorig)
lwm<-update(full.mod, .~. -tb:herdorig); anova(full.mod, lwm) 
lwm<-update(full.mod, .~. -bruc:I(age_yr^2)); anova(full.mod, lwm)
lwm<-update(full.mod, .~. -tb:age_yr); anova(full.mod, lwm)
lwm<-update(full.mod, .~. -bruc:age_yr); anova(full.mod, lwm)

full.mod <-update(full.mod, .~. -tb:age_yr)
lwm<-update(full.mod, .~. -tb:herdorig); anova(full.mod, lwm) 
lwm<-update(full.mod, .~. -bruc:I(age_yr^2)); anova(full.mod, lwm)
lwm<-update(full.mod, .~. -bruc:age_yr); anova(full.mod, lwm)

full.mod <-update(full.mod, .~. -bruc:I(age_yr^2))
lwm<-update(full.mod, .~. -tb:herdorig); anova(full.mod, lwm) 
lwm<-update(full.mod, .~. -bruc:age_yr); anova(full.mod, lwm)
lwm<-update(full.mod, .~. -bka_delay); anova(full.mod, lwm)


full.mod <-update(full.mod, .~. -bruc:age_yr)
lwm<-update(full.mod, .~. -tb:herdorig); anova(full.mod, lwm) 
lwm<-update(full.mod, .~. -bka_delay); anova(full.mod, lwm)
lwm<-update(full.mod, .~. -bruc); anova(full.mod, lwm)
lwm<-update(full.mod, .~. -age_yr); anova(full.mod, lwm)
lwm<-update(full.mod, .~. -I(age_yr^2)); anova(full.mod, lwm)

full.mod <-update(full.mod, .~. -tb:herdorig)
lwm<-update(full.mod, .~. -bka_delay); anova(full.mod, lwm)
lwm<-update(full.mod, .~. -bruc); anova(full.mod, lwm)
lwm<-update(full.mod, .~. -age_yr); anova(full.mod, lwm)
lwm<-update(full.mod, .~. -I(age_yr^2)); anova(full.mod, lwm)
lwm<-update(full.mod, .~. -tb); anova(full.mod, lwm)
lwm<-update(full.mod, .~. -herdorig); anova(full.mod, lwm)

full.mod <-update(full.mod, .~. -I(age_yr^2))
lwm<-update(full.mod, .~. -bka_delay); anova(full.mod, lwm)
lwm<-update(full.mod, .~. -bruc); anova(full.mod, lwm)
lwm<-update(full.mod, .~. -age_yr); anova(full.mod, lwm)
lwm<-update(full.mod, .~. -tb); anova(full.mod, lwm)
lwm<-update(full.mod, .~. -herdorig); anova(full.mod, lwm)


full.mod <-update(full.mod, .~. -bruc)
lwm<-update(full.mod, .~. -bka_delay); anova(full.mod, lwm)
lwm<-update(full.mod, .~. -age_yr); anova(full.mod, lwm)
lwm<-update(full.mod, .~. -tb); anova(full.mod, lwm)
lwm<-update(full.mod, .~. -herdorig); anova(full.mod, lwm)

full.mod <-update(full.mod, .~. -bka_delay)
lwm<-update(full.mod, .~. -age_yr); anova(full.mod, lwm)
lwm<-update(full.mod, .~. -tb); anova(full.mod, lwm)
lwm<-update(full.mod, .~. -herdorig); anova(full.mod, lwm)

full.mod <-update(full.mod, .~. -herdorig)
lwm<-update(full.mod, .~. -age_yr); anova(full.mod, lwm)
lwm<-update(full.mod, .~. -tb); anova(full.mod, lwm)
# keep full.mod


# summary 
#summary(full.mod)
#Linear mixed model fit by REML ['lmerMod']
#Formula: log(hapto) ~ tb + age_yr + (1 | hapto_plate) + (capturetime |      id)
#   Data: immune3

#REML criterion at convergence: 1909.2

#Scaled residuals: 
#    Min      1Q  Median      3Q     Max 
#-4.6200 -0.5121  0.0160  0.4731  3.2776 
#
#Random effects:
# Groups      Name        Variance  Std.Dev. Corr 
# id          (Intercept) 1.031e+00 1.015348      
#             capturetime 9.131e-05 0.009556 -0.96
# hapto_plate (Intercept) 5.187e-01 0.720193      
# Residual                1.145e+00 1.069815      
#Number of obs: 570, groups:  id, 131; hapto_plate, 33

#Fixed effects:
#            Estimate Std. Error t value
#(Intercept)  5.76870    0.27591  20.908
#tb           0.28581    0.14516   1.969
#age_yr       0.07113    0.03804   1.870

#Correlation of Fixed Effects:
#       (Intr) tb    
#tb     -0.108       
#age_yr -0.801 -0.101


coefs <- data.frame(coef(summary(full.mod)))
df.KR <- get_ddf_Lb(full.mod, fixef(full.mod))
coefs$p.KR <- 2 * (1 - pt(abs(coefs$t.value), df.KR))
coefs
#              Estimate Std..Error   t.value       p.KR
#(Intercept) 5.76870131 0.27590520 20.908273 0.00000000
#tb          0.28581390 0.14516258  1.968923 0.05087667
#age_yr      0.07113059 0.03804066  1.869857 0.06352866


test<-lmer(log(hapto)~tb*age_yr+ (1|hapto_plate) + (capturetime|id), data=immune3)
 	
 	
# Step 3- Repeat model selection without plate... because Bree says so. 
######################################### 	
full.mod<-lmer(log(hapto)~bruc+herdorig+tb+age_yr+I(age_yr^2)+ 
 	tb*age_yr+ tb*I(age_yr^2)+ tb*herdorig+ bruc*age_yr+ bruc*I(age_yr^2)+ 
 	bruc*herdorig+ tb*bruc+bka_delay + (capturetime|id), data=immune3)
plot(resid(full.mod), immune3$capturetime)

lwm<-update(full.mod, .~. -tb:I(age_yr^2)); anova(full.mod, lwm) 
lwm<-update(full.mod, .~. -bruc:I(age_yr^2)); anova(full.mod, lwm)
lwm<-update(full.mod, .~. -tb:herdorig); anova(full.mod, lwm)  
lwm<-update(full.mod, .~. -bruc:herdorig); anova(full.mod, lwm)
lwm<-update(full.mod, .~. -bruc:age_yr); anova(full.mod, lwm)
lwm<-update(full.mod, .~. -tb:age_yr); anova(full.mod, lwm)
lwm<-update(full.mod, .~. -tb:bruc); anova(full.mod, lwm)

full.mod<-update(full.mod, .~. -bruc:I(age_yr^2)) 	#  0.2589      1     0.6109
############################################################
lwm<-update(full.mod, .~. -tb:I(age_yr^2)); anova(full.mod, lwm)
lwm<-update(full.mod, .~. -tb:herdorig); anova(full.mod, lwm)  
lwm<-update(full.mod, .~. -bruc:herdorig); anova(full.mod, lwm)
lwm<-update(full.mod, .~. -bruc:age_yr); anova(full.mod, lwm)
lwm<-update(full.mod, .~. -tb:age_yr); anova(full.mod, lwm)
lwm<-update(full.mod, .~. -tb:bruc); anova(full.mod, lwm)

full.mod<-update(full.mod, .~. -tb:bruc) # 0.9442    1     0.3312
############################################################
lwm<-update(full.mod, .~. -tb:I(age_yr^2)); anova(full.mod, lwm)
lwm<-update(full.mod, .~. -tb:herdorig); anova(full.mod, lwm)  
lwm<-update(full.mod, .~. -bruc:herdorig); anova(full.mod, lwm)
lwm<-update(full.mod, .~. -bruc:age_yr); anova(full.mod, lwm)
lwm<-update(full.mod, .~. -tb:age_yr); anova(full.mod, lwm)

full.mod<-update(full.mod, .~. -tb:I(age_yr^2))# 1.6375   1     0.2007
############################################################
lwm<-update(full.mod, .~. -tb:herdorig); anova(full.mod, lwm)  
lwm<-update(full.mod, .~. -bruc:herdorig); anova(full.mod, lwm)
lwm<-update(full.mod, .~. -bruc:age_yr); anova(full.mod, lwm)
lwm<-update(full.mod, .~. -tb:age_yr); anova(full.mod, lwm)

full.mod <-update(full.mod, .~. -tb:age_yr) # 0.2365      1     0.6267
############################################################
lwm<-update(full.mod, .~. -tb:herdorig); anova(full.mod, lwm)  
lwm<-update(full.mod, .~. -bruc:herdorig); anova(full.mod, lwm)
lwm<-update(full.mod, .~. -bruc:age_yr); anova(full.mod, lwm)



full.mod<-update(full.mod, .~. -bruc:age_yr) # 1.4041      1      0.236
############################################################
lwm<-update(full.mod, .~. -tb:herdorig); anova(full.mod, lwm)  
lwm<-update(full.mod, .~. -bruc:herdorig); anova(full.mod, lwm)

full.mod<-update(full.mod, .~. -tb:herdorig) # 1.3375      1     0.2475
############################################################
lwm<-update(full.mod, .~. -bruc:herdorig); anova(full.mod, lwm) 
lwm<-update(full.mod, .~. -bka_delay); anova(full.mod, lwm) 
lwm<-update(full.mod, .~. -age_yr); anova(full.mod, lwm) 
lwm<-update(full.mod, .~. -I(age_yr^2)); anova(full.mod, lwm) 
lwm<-update(full.mod, .~. -tb); anova(full.mod, lwm) 

full.mod<-update(full.mod, .~. -bruc:herdorig) # 2.3986      1     0.1214
############################################################
lwm<-update(full.mod, .~. -bka_delay); anova(full.mod, lwm) 

full.mod<-update(full.mod, .~. -bka_delay) # 0.0447      1     0.8326
############################################################
lwm<-update(full.mod, .~. -age_yr); anova(full.mod, lwm) 
lwm<-update(full.mod, .~. -I(age_yr^2)); anova(full.mod, lwm) 
lwm<-update(full.mod, .~. -tb); anova(full.mod, lwm) 
lwm<-update(full.mod, .~. -bruc); anova(full.mod, lwm) 
lwm<-update(full.mod, .~. -herdorig); anova(full.mod, lwm) 

full.mod <-update(full.mod, .~. -herdorig) #0.4776      1     0.4895
############################################################
lwm<-update(full.mod, .~. -age_yr); anova(full.mod, lwm) 
lwm<-update(full.mod, .~. -I(age_yr^2)); anova(full.mod, lwm) # 1.7662      1     0.1839
lwm<-update(full.mod, .~. -tb); anova(full.mod, lwm) 
lwm<-update(full.mod, .~. -bruc); anova(full.mod, lwm)  # 1.3338      1     0.2481

t1<-update(full.mod, .~. -bruc)
t2<-update(full.mod, .~. -I(age_yr^2))  # doesn't matter which is dropped... 


full.mod <-update(full.mod, .~. -I(age_yr^2)) # 1.7662      1     0.1839
############################################################
lwm<-update(full.mod, .~. -age_yr); anova(full.mod, lwm) 
lwm<-update(full.mod, .~. -tb); anova(full.mod, lwm) 
lwm<-update(full.mod, .~. -bruc); anova(full.mod, lwm) 

full.mod <-update(full.mod, .~. -bruc)# 1.0339      1     0.3092
############################################################
lwm<-update(full.mod, .~. -age_yr); anova(full.mod, lwm) #3.3756      1    0.06617 
lwm<-update(full.mod, .~. -tb); anova(full.mod, lwm) 


coefs <- data.frame(coef(summary(full.mod)))
df.KR <- get_ddf_Lb(full.mod, fixef(full.mod))
coefs$p.KR <- 2 * (1 - pt(abs(coefs$t.value), df.KR))
coefs
#               Estimate Std..Error   t.value       p.KR
# (Intercept)   5.64262409 0.23323891 24.192464 0.00000000
# tb            0.30140289 0.16108493  1.871081 0.06286005
# age_yr        0.08749377 0.03663146  2.388487 0.01789021

plot(resid(full.mod)~immune3$capturetime)

#Random effects:
# Groups   Name        Variance  Std.Dev. Corr 
# id       (Intercept) 1.2911428 1.1363        
#          capturetime 0.0005904 0.0243   -0.74
# Residual             1.5078177 1.2279        
#Number of obs: 570, groups:  id, 131

#Fixed effects:
#            Estimate Std. Error t value
#(Intercept)  5.64262    0.23324  24.192
#tb           0.30140    0.16108   1.871
#age_yr       0.08749    0.03663   2.388
#
#Correlation of Fixed Effects:
#       (Intr) tb    
#tb     -0.108       
#age_yr -0.894 -0.147

immune.sample <-subset(immune3, id %in% c("B19", "Y13", "B5", "B41", "B31", "O47b", "Y38", "B1", "O12"))
xyplot(log(hapto)~ capturetime|id, data=immune.sample, type=c("p", "smooth"), lwd=2)

immune.sample <-subset(immune3, id %in% c("B19", "Y13", "B5", "B41", "O47b", "O12"))
xyplot(log(hapto)~capturetime|id, data=immune.sample, type=c("p", "smooth"), lwd=2)

plot(log(hapto) ~ capturetime, data = immune3, type = "n",
    xlab = "Time (months)", ylab = expression(log("haptoglobin concentration")))
for (i in unique(immune3$id))
    lines(log(hapto) ~ capturetime, data = immune3[immune3$id == i, ],
        lty = match(i, unique(immune3$id)))
        
plot(log(hapto) ~ capturetime, data = immune3, type = "n",
    xlab = "Time (months)", ylab = expression(log("haptoglobin concentration")))
for (i in unique(immune3$id)[1:30])
    lines(log(hapto) ~ capturetime, data = immune3[immune3$id == i, ],
        lty = match(i, unique(immune3$id)), )

immune3$TBcol[immune3$convert==1]<- "red"
immune3$TBcol[immune3$convert==0]<- "blue"
immune3$TBcol[immune3$tb_beforeafter=="pfc"]<- "dark green"

immune3$Brlty[immune3$brucconvert ==1]<- 1
immune3$Brlty[immune3$brucconvert ==0]<- 2
immune3$Brlty[immune3$bruc_beforeafter =="pfc"]<- 3

immune3$Brpch[immune3$brucconvert ==1]<- 15
immune3$Brpch[immune3$brucconvert ==0]<- 17
immune3$Brpch[immune3$bruc_beforeafter =="pfc"]<- 19

plot(log(hapto) ~ capturetime, data = immune3, type = "n",
    xlab = "Time (months)", ylab = expression(log("haptoglobin concentration")))
for (i in unique(immune3$id)[1:60])
    points(log(hapto) ~ jitter(capturetime), data = immune3[immune3$id == i, ],
        pch = immune3$Brpch[immune3$id==i], col= immune3$TBcol[immune3$id==i], cex=0.8)

for (i in unique(immune3$id)[1:60])
    lines(log(hapto) ~ jitter(capturetime), data = immune3[immune3$id == i, ],
        lty = immune3$Brlty[immune3$id==i], col= immune3$TBcol[immune3$id==i])
legend("bottomright", legend=c("TB non-converter", "TB, converter", "TB pfc"), lty= c(1,1,1), col=c("blue", "red", "dark green"), bty="n")
legend("bottomleft", legend=c("Bruc non-converter", "Bruc, converter", "Bruc pfc"), pch= c(17,15,19), bty="n")

# Smooth longitudinal profiles of 16 subjects from our dataset.  The solid lines represents the fit of the loess smoother
############################################################
# Other immunological measures- Before vs. After infection
############################################################
all$hapto<-as.numeric(as.character(all$hapto))
all2<-all[!is.na(all$hapto),]
ggplot(data=all2, aes(x=floor(tsi_tb), y=log(hapto), group=id, colour=bruc))+ geom_line()+ 
	xlab("Time since bTB seroconversion") + ylab("Log haptoglobin concentration")+ theme_bw() +
	theme(axis.title.x= element_text(size=16), axis.title.y = element_text(size=16), 
	axis.text.x = element_text(size=14, vjust=-0.05), axis.text.y = element_text(size=14),
	axis.line= element_line(colour="black"), 
	legend.position=c(0.82, 0.9), 
	legend.background= element_rect(fill="white", colour="white"), 
	legend.key= element_blank(),
    legend.title= element_blank(),
    legend.text = element_text(size=15) ) + guides(fill=guide_legend(title="Brucellosis"))



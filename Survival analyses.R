############################################
############################################
######       Survival Analyses        ######
############################################
############################################

# Started: 1 April 2016
library('JMbayes')
#http://www.r-bloggers.com/joint-models-for-longitudinal-and-survival-data/
#http://www.r-bloggers.com/dynamic-predictions-using-joint-models/
setwd("~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/final_datasets_copied_from_phdfolder/survival")
data<-read.csv("~/Documents/postdoc_buffology/Last-Thesis-Chapter!!!!!!/final_datasets_copied_from_phdfolder/survival/brucsurvival_TB3controls_longresidnomissing_noerrors_season2.csv")
data$age_yr2<-floor(data$age_yr)


# run cph model


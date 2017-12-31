#statistical tests
summary_stats= function(btb, bruc, age, id, ok){
  # overall chi-square test: 
  test<-chisq.test(x=btb, y=bruc)
  table(btb, bruc)

  # GLMER
  test1<-glmer(bruc~btb+(1|id), family=binomial(link="logit"))
  test2<-glmer(bruc~btb+age+(1|id), family=binomial(link="logit"))
  if (ok==TRUE){
  test3<-glmer(bruc~btb+age+I(age^2)+(1|id), family=binomial(link="logit"))
  test4<-glmer(bruc~btb*age+I(age^2)+(1|id), family=binomial(link="logit"))
  test5<-glmer(bruc~age+btb*I(age^2)+(1|id), family=binomial(link="logit"))
  test6<-glmer(bruc~btb*age+btb*I(age^2)+(1|id), family=binomial(link="logit"))
  
  summary<-data.frame(test=c("Chisquare", "GLMER tb", "GLMER tb+age", "GLMER tb+age+age^2", 
                             "GLMER+ tb*age+age^2", "GLMER age+tb*age^2", "GLMER tb*age+tb*age^2"), 
                      estimate=c(0, test1@beta[2], test2@beta[2], test3@beta[2], test4@beta[2], 
                                 test5@beta[2], test6@beta[2]), 
                      pval=c(test$p.value, summary(test1)$coefficients[[8]], summary(test2)$coefficients[[11]],
                             summary(test3)$coefficients[[14]], summary(test4)$coefficients[[17]], 
                             summary(test5)$coefficients[[17]], summary(test6)$coefficients[[20]]), 
                      AIC= c(0, summary(test1)$AICtab[[1]], summary(test2)$AICtab[[1]], 
                             summary(test3)$AICtab[[1]], summary(test4)$AICtab[[1]], 
                             summary(test5)$AICtab[[1]], summary(test6)$AICtab[[1]]))
  }
  if (ok==FALSE){
    summary<-data.frame(test=c("Chisquare", "GLMER tb", "GLMER tb+age"), 
                        estimate=c(0, test1@beta[2], test2@beta[2]), 
                        pval=c(test$p.value, summary(test1)$coefficients[[8]], summary(test2)$coefficients[[11]]), 
                        AIC= c(0, summary(test1)$AICtab[[1]], summary(test2)$AICtab[[1]]) )  
  }
  
  return(summary)
}
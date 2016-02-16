make_age_prev_plot= function(btb, bruc, age, binsize){
  require(gplots)
  ######################################################
  # Input: 
  ##btb= columns for bTB presence absence (with 0, 1 values)
  ##bruc= columns for brucellosis presence absence (positive, negative row values)
  ##binsize= how many years do you want to lump together
  # Plot and summary of prevalence 
  ######################################################
  agebins=c(seq(0, max(age), binsize), max(age+1))
  bruc2<-NA
  for (i in 1:length(bruc)){
    ifelse(bruc[i]=="negative", bruc2[i]<-0, bruc2[i]<-1)
  }
  bruc<-bruc2
  olddf<-data.frame(cbind(btb, bruc, age))
  newdf<-data.frame(agebins, prevneg=NA, prevpos=NA)   # one row longer than new df. 
  for (i in 1:(length(newdf[,1])-1)){
    neg<-olddf[olddf$btb==0,]
    pos<-olddf[olddf$btb==1,]
    olddf_neg<-olddf[olddf$age>=agebins[i] & olddf$age<agebins[i+1] & olddf$btb=="0",]
    olddf_pos<-olddf[olddf$age>=agebins[i] & olddf$age<agebins[i+1] & olddf$btb=="1",]
    newdf$prevneg[i]<-get_prev(olddf_neg$bruc)
    newdf$prevpos[i]<-get_prev(olddf_pos$bruc)
  }
  newdf<-newdf[!is.na(newdf$prevneg),]
  
  plotmat<-as.matrix(rbind(newdf$prevneg, newdf$prevpos), 
                     dimnames=list(c(agebins[-length(agebins)] ), c("btbneg", "btbpos") ) )
  rownames(plotmat)<-c("btbneg", "btbpos")
  colnames(plotmat)<-agebins[-length(agebins)]
  
  # make plot
  barplot(plotmat, beside=TRUE, col=c("light grey", "gray45"), ylim=c(0, max(plotmat)+0.1), 
          ylab="Brucellosis Prevalence")  
}


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
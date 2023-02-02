fitY$q025




#write.csv(fitY,"fitY.csv")

#  write.csv(x,"x.csv")
write.csv(bi_lst,"5weeksresults.csv")
# 提取 出  1 60 

cul <- function(x) {
   res = list()
   # n = length(x) / 60
   for (i in  1:366) {
     res[[i]] = as.numeric(x[seq(i,73200,366)])
   } 
   return(res)
  
}

cul <- function(x) {
   res = list()
   # n = length(x) / 60
   for (i in  1:420) {
      res[[i]] = as.numeric(x[seq(i,48000,420)])
   } 
   return(res)
   
}


mymean <- function(alist) {
   n = length(alist)
   res = c()
   for (i in 1:n) {
      res = c(res,mean(alist[[i]],na.rm = T))
   }
   return(res)
}

myq025 <- function(alist) {
   n = length(alist)
   res = c()
   for (i in 1:n) {
      res = c(res,quantile(alist[[i]],0.025))
   }
   return(res)
}


myq25 <- function(alist) {
   n = length(alist)
   res = c()
   for (i in 1:n) {
      res = c(res,quantile(alist[[i]],0.25))
   }
   return(res)
}




myq05 <- function(alist) {
   n = length(alist)
   res = c()
   for (i in 1:n) {
      res = c(res,quantile(alist[[i]],0.5))
   }
   return(res)
}



myq75 <- function(alist) {
   n = length(alist)
   res = c()
   for (i in 1:n) {
      res = c(res,quantile(alist[[i]],0.75))
   }
   return(res)
}



myq975 <- function(alist) {
   n = length(alist)
   res = c()
   for (i in 1:n) {
      res = c(res,quantile(alist[[i]],0.975))
   }
   return(res)
}



mymedian <- function(alist) {
   n = length(alist)
   res = c()
   for (i in 1:n) {
      res = c(res,quantile(alist[[i]],0.5))
   }
   return(res)
}



MSEIR <- read.csv("60wbeta1.csv", header=TRUE, stringsAsFactors=FALSE)
Evalue =cul(MSEIR$q50)
Eavgdata = mymean(Evalue)
Eavgdata = mymedian(Evalue)
Zmodel<-Eavgdata*0.2/5

Svalue=cul(SEIR$S.value)
Savgdata=mymean(Svalue)
Ivalue=cul(SEIR$I.value)
Iavgdata=mymean(Ivalue)
Rvalue=cul(SEIR$R.value)
Ravgdata=mymean(Rvalue)*2.2
Mvalue=cul(MSEIR$M.value)
Mavgdata=mymean(Mvalue)
avgdata025 = myq025(avgvalue)
avgdata25 = myq25(avgvalue)
avgdata975 = myq975(avgvalue)
avgdata05= myq05(avgvalue)
avgdata75= myq75(avgvalue)

copy_y = y
copy_y$mean = avgdata
copy_y$q025 = avgdata025
copy_y$q25=avgdata25
copy_y$q75=avgdata75
copy_y$q975=avgdata975


ggplot(data = copy_y) +
   geom_ribbon(aes(x = time, ymin = q25, ymax = q75), alpha = 0.3) +
   geom_ribbon(aes(x = time, ymin = q025, ymax = q975), alpha = 0.3) +
   geom_line(aes(x = time, y = mean)) +
   geom_point(aes(x = time, y = value), colour = "Red") +
   ylab("Incidence") +
   xlab("Time")  
#+  ylim(c(16000,30000))



#avgvalue025 =cul(fitY$q025)
#avgvalue25 =cul(fitY$q25)
#avgvalue75 =cul(fitY$q75)
#avgvalue975 =cul(fitY$q975)


library(tidyverse)
copy_y = y
copy_y$postmean = avgdata
#  
y_copy = copy_y %>%  pivot_longer(cols = -time)

ggplot(y_copy,aes(x = time, y = value,color = name,group = name)) + geom_point() + geom_line()+
   ylab("Daily new H1N1 clinical cases")+
   xlab("Time-Day")

#Reproduction number plot

fitE <- bi_lst$E %>% 
   group_by(time) %>%
   mutate(
      q025 = quantile(value, 0.025),
      q25 = quantile(value, 0.25),
      q50 = quantile(value, 0.5),
      q75 = quantile(value, 0.75),
      q975 = quantile(value, 0.975)
   ) %>% ungroup() 
extractE<-cul(fitE$value)
posteriormeanE<-mymean(extractE)

fitI <- bi_lst$I %>% 
   group_by(time) %>%
   mutate(
      q025 = quantile(value, 0.025),
      q25 = quantile(value, 0.25),
      q50 = quantile(value, 0.5),
      q75 = quantile(value, 0.75),
      q975 = quantile(value, 0.975)
   ) %>% ungroup() 
extractI<-cul(fitI$value)
posteriormeanI<-mymean(extractI)
ratioEI <-posteriormeanE/posteriormeanI

fitbeta <- bi_lst$x %>% mutate(value = exp(value)) %>%
   group_by(time) %>%
   mutate(
      q025 = quantile(value, 0.025),
      q25 = quantile(value, 0.25),
      q50 = quantile(value, 0.5),
      q75 = quantile(value, 0.75),
      q975 = quantile(value, 0.975)
   ) %>% ungroup()
extractbeta<-cul(fitbeta$value)
posteriormeanbeta<-mymean(extractbeta)

posteriormeank<-mymean(1/bi_lst$k)
posteriormeangamma<-mymean(1/bi_lst$gamma)
reprodnumI<-(posteriormeank)*ratioEI/(posteriormeangamma+0.0087)
#reprodnumI<-ts(reprodnumI)
reprodnumE<-posteriormeanbeta*(1+0.1/ratioEI)/(posteriormeangamma+posteriormeank)
reprodnumE<-ts(reprodnumE)
plot(reprodnumI,type = "b",col="Orange",ylab = TeX("Reproduction Num ($\\R_{t}^{I}$)"),xlab="Time-weekly")
abline(h=1)
#plot(reprodnumE,type = "l",col="Blue",ylab = TeX("Reproduction Num ($\\R_{t}^{E}$)"),xlab="Time-weekly")
#abline(h=1)

#Trace Plot
par(mfrow=c(2,1))
plot(1/bi_lst$k$value,type='l',main=TeX("Sampled onset London COVID-19 symptoms rate ($\\alpha$)"),xlab="PMCMC iterations after first 5000 burn-in and thinning by 5",ylab="alpha")
plot(1/bi_lst$gamma$value,type='l',main=TeX("Sampled London COVID-19 Recovery rate ($\\gamma$)"),xlab="PMCMC iterations after first 5000 burn-in and thinning by 5",ylab="gamma")


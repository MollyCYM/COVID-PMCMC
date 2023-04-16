rm(list=ls())
require(deSolve)
library(latex2exp)
set.seed(012)   #random seed 01 03 012
times <- 1:365
N=52196381
#Simulate a Brownian Motion Path
e <- rnorm(n = length(times) - 1, sd = sqrt(0.01))
e <- c(0, cumsum(e))

mu<-vector(length=365) 
for (t in 1:365) {
  if (t <= 120){
    mu[t]=-0.02      }                    #No lock-down policy
  else {mu[t]=-0.2-0.02 }                 #Lock-down policy
}
mu<-ts(mu)
# write.csv(mu,"poster2_mu1.csv")

#plot(mu,type='l')
#Simulate an O-U Process
ornstein_uhlenbeck <- function(n,theta,sigma,x0){
  dt  <- 1
  x<-vector(length=366)
  for (i in 1:(n+1)) {
    if (i==1){x[i]=x0}
    else{
      x[i]  <-  x[i-1] + theta*(mu[i-1]-x[i-1])*dt + sigma*e[i-1]}
  }
  return(x);
}

x<- ornstein_uhlenbeck(365,0.05,sqrt(0.004),log(0.8))
plot(x,type='l',xlab="time")
abline(v=121, col="red")
lines(mu, col="blue")
tt1 <-expression(mu==-0.02)
text(30,-0.1,tt1,col="green")
tt2 <-expression(mu==-0.02-0.2)
text(210,-0.3,tt2,col="green")
plot(exp(x),type='l',xlab="time",ylab=TeX("$\\e^{x}=beta$"),ylim=c(0.2,1.2))
abline(v=121, col="red")
lines(exp(mu+0.02),col="blue")
axis(1, at=121,labels=121, col.axis="red", las=2)
tt1 <-expression(mu[beta]==1)
text(50,1.1,tt1,col="green")
tt2 <-expression(mu[beta]==0.81)
text(220,0.9,tt2,col="green")
tt3 <-expression('lockdown policy start')
text(120,0.4,tt3,col="red")
# write.csv(x,"poster2_x1.csv")
#Main ODE Model
Covid_OU <- function(time, current_state, params){
  
  with(as.list(c(current_state, params)),{
    
    N <- S+E+I+R
    dt <- 1
    dS <- -exp(x[t+1])*S*(E+0.1*I)/N
    dE <- exp(x[t+1])*S*(E+0.1*I)/N - E*(1/k+1/gamma)
    dI <- E/k-I*(1/gamma+0.0087)
    dR <- (I+E)/gamma+0.0087*I
    
    return(list(c( dt, dS, dE, dI, dR)))
  })
}
params <- c(k=5, gamma=5)   #Trial: k=5 & gamma=5; k=7 & gamma=5

S <- N-1
E <- 1
R <- 0
I <- 0

initial_state<- c(t=1, S=S, E=E, I=I, R=R)
#initial_state <- c(t=1,S=52196380, E=1, I=0, R=0)
model5 <- ode(initial_state, times, Covid_OU, params)

summary(model5)

matplot(model5, type="l", lty=1, main="COVID—OU Model",ylab="counts", xlab="Time")
legend <- colnames(model5)[3:6]
legend("right", legend=legend, col=3:6, lty = 1)
Z2 <-model5[,4]/5

tau2 <- 0.8#runif(1,0,1)
Y2 <-vector(length = 365)
for (i in 1:365){
  #Y1[i]<- rlnorm(1,log(Z1[i,]/5),tau1)
  Y2[i]<- rlnorm(1,log(Z2[i]/5),tau2)
  
}
plot(Y2,type='l',ylab="Obs Y",xlab="time",col="blue")
# write.csv(Y2,"Covidou_Y5.csv")
# write.csv(model5,"Covidou5.csv")
# tryexp <- c(.3, 1.03, 2.67, 5, 8.91)
# round(tryexp, digits = 0)
###################################################################################
###################################################################################
#ODE solver version 1-put x inside ODE system
rm(list=ls())
require(deSolve)
library(latex2exp)
set.seed(012)   #random seed 01 03 012
times <- 1:365
N=52196381
#Simulate a Brownian Motion Path
e <- rnorm(n = length(times) - 1, sd = sqrt(0.01))
e <- c(0, cumsum(e))

mu<-vector(length=365) 
for (t in 1:365) {
  if (t <= 120){
    mu[t]=-0.02      }                    #No lock-down policy
  else {mu[t]=-0.2-0.02 }                 #Lock-down policy
}
mu<-ts(mu)
# write.csv(mu,"poster2_mu1.csv")
#plot(mu,type='l')

#Main ODE Model
Covid_OU <- function(time, current_state, params){
  
  with(as.list(c(current_state, params)),{
    
    N <- S+E+I+R
    dt <- 1
    dx <- 0.05*(mu[t]-x)*dt + sqrt(0.004)*e[t]
    dS <- -exp(x)*S*(E+0.1*I)/N
    dE <- exp(x)*S*(E+0.1*I)/N - E*(1/k+1/gamma)
    dI <- E/k-I*(1/gamma+0.0087)
    dR <- (I+E)/gamma+0.0087*I
    
    return(list(c( dt, dS, dE, dI, dR, dx)))
  })
}
params <- c(k=5, gamma=5)   #Trial: k=5 & gamma=5; k=7 & gamma=5

S <- N-1
E <- 1
R <- 0
I <- 0
x <- log(0.8)
initial_state<- c(t=1, S=S, E=E, I=I, R=R, x=x)
#initial_state <- c(t=1,S=52196380, E=1, I=0, R=0)
model6 <- ode(initial_state, times, Covid_OU, params)

summary(model6)
x<- model6[,7]
plot(x,type='l',xlab="time")
abline(v=121, col="red")
lines(mu, col="blue")
tt1 <-expression(mu==-0.02)
text(30,-0.1,tt1,col="green")
tt2 <-expression(mu==-0.02-0.2)
text(210,-0.3,tt2,col="green")
plot(exp(x),type='l',xlab="time",ylab=TeX("$\\e^{x}=beta$"),ylim=c(0.2,1.2))
abline(v=121, col="red")
lines(exp(mu+0.02),col="blue")
axis(1, at=121,labels=121, col.axis="red", las=2)
tt1 <-expression(mu[beta]==1)
text(50,1.1,tt1,col="green")
tt2 <-expression(mu[beta]==0.81)
text(220,0.9,tt2,col="green")
tt3 <-expression('lockdown policy start')
text(120,0.4,tt3,col="red")
write.csv(x,"poster3_x1.csv")
matplot(model6, type="l", lty=1, main="COVID—OU Model",ylab="counts", xlab="Time")
legend <- colnames(model6)[3:6]
legend("right", legend=legend, col=3:6, lty = 1)
Z3 <-model6[,4]/5

tau3 <- 0.8#runif(1,0,1)
Y3 <-vector(length = 365)
for (i in 1:365){
  #Y1[i]<- rlnorm(1,log(Z1[i,]/5),tau1)
  Y3[i]<- rlnorm(1,log(Z3[i]/5),tau3)
  
}
plot(Y3,type='l',ylab="Obs Y",xlab="time",col="blue")
write.csv(Y3,"Covidou_Y6.csv")
write.csv(model6,"Covidou6.csv")
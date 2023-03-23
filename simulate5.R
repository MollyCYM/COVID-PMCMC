rm(list=ls())
require(deSolve)
set.seed(34)   #random seed 32 34 36 37 38 39
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
text(30,0.1,tt1,col="green")
tt2 <-expression(mu==-0.02-0.2)
text(200,-0.1,tt2,col="green")
plot(exp(x),type='l',xlab="time",ylab=TeX("$e^{x}=beta$"))
abline(v=121, col="red")
lines(exp(mu+0.02),col="blue")
tt1 <-expression(mu==1)
text(30,1.1,tt1,col="green")
tt2 <-expression(mu==0.81)
text(200,0.9,tt2,col="green")
write.csv(exp(x),"Exp2_beta.csv")
#Main ODE Model
COVID_OU <- function(time, current_state, params){
  
  with(as.list(c(current_state, params)),{
    
    N <- S+E+I+R
    dt <- 1
    dS <- -exp(x[t+1])*S*I/N
    dE <- exp(x[t+1])*S*I/N - E/k
    dI <- E/k-I/gamma
    dR <- I/gamma
    
    return(list(c( dt, dS, dE, dI, dR)))
  })
}
params <- c(k=7, gamma=5)   #Trial: k=5 & gamma=5; k=7 & gamma=5

#library('truncnorm')
# R0 <- rtruncnorm(1, a=0, b=1, mean = 0.15, sd = 0.15)
# E0 <-runif(1,-16, -9)
# I0 <-runif(1,-16, -9)
# #x <-runif(1, -5,2)
# S <- N
# R <- R0*S
# S <- S - R
# 
# E <- exp(E0 + log(S))
# S <- S - E
# I <- exp(I0 + log(S))
# S <- S - I

R0 <-0.03
S <- N
R <- S*R0
S <- S - R

E0 <- -15
E <- exp(E0 + log(S))
S <- S - E

I0 <- -10
I <- exp(I0 + log(S))
S <- S - I
initial_state<- c(t=1, S=S, E=E, I=I, R=R)
#initial_state <- c(t=1,S=52196380, E=1, I=0, R=0)
model3 <- ode(initial_state, times, COVID_OU, params)

summary(model3)

matplot(model3, type="l", lty=1, main="COVID—OU Model",ylab="counts", xlab="Time")
legend <- colnames(model3)[3:6]
legend("right", legend=legend, col=3:6, lty = 1)
Z1 <-model3[,4]/7

tau1 <- 0.8#runif(1,0,1)
Y1 <-vector(length = 365)
for (i in 1:365){
  #Y1[i]<- rlnorm(1,log(Z1[i,]/5),tau1)
  Y1[i]<- rlnorm(1,log(Z1[i]/5),tau1)
}
plot(Y1,type='l',ylab="Obs Y",xlab="time",col="blue")
write.csv(Y1,"Exp2_Y.csv")
write.csv(model3,"Exp2_model.csv")


######################################################################

######################################################################
rm(list=ls())
require(deSolve)
times <- 1:365 #Simulate 1 year daily data
N=56000000     #England population size in 2023

#Simulate a Brownian Motion Path
e <- rnorm(n = length(times) - 1, sd = sqrt(0.01))
e <- c(0, cumsum(e))
sigma <- runif(1, min = 0, max = 1)
beta<- exp(sigma*e)

mu<-vector(length=365) 
for (t in 1:365) {
  if (t <= 150){
    mu[t]=exp(-0.92)      }                    #No lock-down policy
  else if (t <=280){                         
    mu[t]=exp(-0.92-1.38)     }                #Strict national lock-down 
  else {mu[t]=exp(-0.92-0.69) }                #Mitigate Tier System
}
mu<-ts(mu)


#Simulate an O-U Process
ornstein_uhlenbeck <- function(n,theta,sigma,beta0){
  #dw  <- rnorm(n, 0, sqrt(T/n))
  dt  <- 1
  beta<-vector(length=366)
  for (i in 1:(n+1)) {
    if (i==1){beta[i]=beta0}
    else{
    beta[i]  <-  beta[i-1] + theta*(mu[i-1]-beta[i-1])*dt + sigma*e[i-1]}
  }
  return(beta);
}

beta<- ornstein_uhlenbeck(365,0,0.02,0)

#Main ODE Model
COVID_OU <- function(time, current_state, params){
  
  with(as.list(c(current_state, params)),{
    
    N <- S+E+I+R
    for (t in 1:times){
    dS <- -beta[t+1]*S*(E+0.1*I)/N
    dE <- beta[t+1]*S*(E+0.1*I)/N - E*(1/alpha+1/gamma)
    dI <- E/alpha-I*(1/gamma+0.0087)
    dR <- (I+E)/gamma+0.0087*I}
    
    return(list(c( dS, dE, dI, dR)))
  })
}

###################################
H1N1 <- function(time, current_state, params){
  
  with(as.list(c(current_state, params)),{
    
    N <- S+E+I+R
    dt <- 1
    dS <- -beta[t]*S*I/N
    dE <- beta[t]*S*I/N - E/k
    dI <- E/k - I/gamma
    dR <- I/gamma
    
    return(list(c(dt, dS, dE, dI, dR)))
  })
}


params <- c(k=1.59, gamma=1.08)

COVID_BM <- function(time, current_state, params){
  
  with(as.list(c(current_state, params)),{
    
    N <- S+E+I+R
      dt <- 1
      dS <- -beta[t]*S*(E+0.1*I)/N
      dE <- beta[t]*S*(E+0.1*I)/N - E*(1/alpha+1/gamma)
      dI <- E/alpha-I*(1/gamma+0.0087)
      dR <- (I+E)/gamma+0.0087*I
    
    return(list(c(dt, dS, dE, dI, dR)))
  })
}
params <- c(alpha=1.59, gamma=1.08)
###################################

params <- c(alpha=0.125, gamma=0.07)

library('truncnorm')
R0 <- rtruncnorm(1, a=0, b=1, mean = 0.15, sd = 0.15)
E0 <-runif(1,-16, -9)
I0 <-runif(1,-16, -9)
#x <-runif(1, -5,2)
S <- N
R <- R0*S
S <- S - R

E <- exp(E0 + log(S))
S <- S - E
I <- exp(I0 + log(S))
S <- S - I

initial_state<- c(S=S, E=E, I=I, R=R)
#initial_state <- c(S=N-1, E=1, I=0, R=0)
model <- ode(initial_state, times, COVID_BM, params)
summary(model)

matplot(model, type="l", lty=1, main="SEIR model", xlab="Time")
legend <- colnames(model)[3:6]
legend("right", legend=legend, col=3:6, lty = 1)

Z <-model[,3]*0.125

tau=0.8
Y <-vector(length = 365)
for (i in 1:365){
  Y[i]<- rlnorm(1,log(Z[i]/5),tau)
}
plot(Y,type='l')
write.csv(Y,"simcovidY.csv")
write.csv(model,"simulatecovid.csv")

plot(beta,type='l')
abline(v=61, col="blue")
abline(v=281,col="blue")

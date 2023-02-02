rm(list=ls())
require(deSolve)
times <- 1:259
N=52196381
#Simulate a Brownian Motion Path
e <- rnorm(n = length(times) - 1, sd = sqrt(0.01))
e <- c(0, cumsum(e))

mu<-vector(length=259) 
for (t in 1:259) {
  if (t <= 60){
    mu[t]=exp(-0.92)      }                    #No lock-down policy
  else if (t <= 210){                         
    mu[t]=exp(-0.92-1.38)     }                #Strict national lock-down 
  else {mu[t]=exp(-0.92-0.69) }                #Mitigate Tier System
}
mu<-ts(mu)


#Simulate an O-U Process
ornstein_uhlenbeck <- function(n,theta,sigma,beta0){
  dt  <- 1
  beta<-vector(length=260)
  for (i in 1:(n+1)) {
    if (i==1){beta[i]=beta0}
    else{
      beta[i]  <-  beta[i-1] + theta*(mu[i-1]-beta[i-1])*dt + sigma*e[i-1]}
  }
  return(beta);
}

beta<- ornstein_uhlenbeck(259,0.3,0.02,0)
plot(beta,type='l')
abline(v=61, col="red")
abline(v=211,col="red")
lines(mu,col="blue")
#Main ODE Model
COVID_OU <- function(time, current_state, params){
  
  with(as.list(c(current_state, params)),{
    
    N <- S+E+I+R
    dt <- 1
    dS <- -beta[t]*S*(E+0.1*I)/N
    dE <- beta[t]*S*(E+0.1*I)/N - E*(1/k+1/gamma)
    dI <- E/k-I*(1/gamma+0.0087)
    dR <- (I+E)/gamma+0.0087*I
    
    return(list(c( dt, dS, dE, dI, dR)))
  })
}
params <- c(k=1.59, gamma=1.08)

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


initial_state<- c(t=1, S=S, E=E, I=I, R=R)
#initial_state <- c(t=1,S=52196380, E=1, I=0, R=0)
model1 <- ode(initial_state, times, COVID_OU, params)

summary(model1)

matplot(model1, type="l", lty=1, main="SEIR model", xlab="Time")
legend <- colnames(model1)[3:6]
legend("right", legend=legend, col=3:6, lty = 1)
Z1 <-model1[,4]/1.59

tau1 <- 0.8#runif(1,0,1)
Y1 <-vector(length = 259)
for (i in 1:259){
  #Y1[i]<- rlnorm(1,log(Z1[i,]/5),tau1)
  Y1[i]<- rlnorm(1,log(Z1[i]/5),tau1)
}
plot(Y1,type='l')
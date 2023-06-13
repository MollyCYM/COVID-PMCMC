rm(list=ls())
require(deSolve)

times <- 1:365
N=52196381

e <- rnorm(n = length(times) , sd = 1)

mu<-vector(length=365) 
for (t in 1:365) {
  if (t <= 120){
    mu[t]=-0.02      }                    #No lock-down policy
  else {mu[t]=-0.2-0.02 }                 #Lock-down policy
}
mu<-ts(mu)


#Simulate an O-U Process
ornstein_uhlenbeck <- function(n,theta,sigma,x0){
  dt  <- 1
  x<-vector(length=365)
  for (i in 1:n) {
    if (i==1){x[i]=x0}
    else{
      x[i]  <-  x[i-1] + theta*(mu[i]-x[i])*dt + sigma*e[i]}
  }
  return(x);
}

x<- ornstein_uhlenbeck(365,0.05,sqrt(0.004),log(0.8))

#Main ODE Model
Covid_OU <- function(time, current_state, params){
  
  with(as.list(c(current_state, params)),{
    
    # N <- S+E+I+R
    dt <- 1
    dS <- -exp(x[t])*S*(E+0.1*I)/N
    dE <- exp(x[t])*S*(E+0.1*I)/N - E*(1/k+1/gamma)
    dI <- E/k-I*(1/gamma+0.0087)
    dR <- (I+E)/gamma+0.0087*I
    
    return(list(c( dt, dS, dE, dI, dR)))
  })
}
params <- c(k=7, gamma=9)   
initial_state <- c(t=1,S=N-1, E=1, I=0, R=0)
model4 <- ode(initial_state, times, Covid_OU, params)

summary(model4)

Z1 <-model4[,4]/7

tau1 <- 0.1
Y1 <-vector(length = 365)
for (i in 1:365){
  Y1[i]<- rlnorm(1,log(Z1[i]/5),tau1)
}

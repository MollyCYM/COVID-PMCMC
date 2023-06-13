rm(list=ls())
require(deSolve)
times <- 1:365
N=52196381
sigma <- 0.07

e <- rnorm(n = length(times)-1 , sd = 1)
## now compute their cumulative sum
e <- c(0, cumsum(e))

beta<- exp(sigma*e)

H1N1 <- function(time, current_state, params){
  
  with(as.list(c(current_state, params)),{
    
    dt <- 1
    dS <- -beta[t]*S*I/N
    dE <- beta[t]*S*I/N - E/k
    dI <- E/k - I/gamma
    dR <- I/gamma
    
    return(list(c(dt, dS, dE, dI, dR)))
  })
}


params <- c(k=1.59, gamma=1.08)

library('truncnorm')
R0 <- rtruncnorm(1, a=0, b=1, mean = 0.15, sd = 0.15)
E0 <-runif(1,-16, -9)
I0 <-runif(1,-16, -9)
S <- N
R <- R0*S
S <- S - R

E <- exp(E0 + log(S))
S <- S - E
I <- exp(I0 + log(S))
S <- S - I


initial_state<- c(t=1, S=S, E=E, I=I, R=R)
model1 <- ode(initial_state, times, H1N1, params)

summary(model1)


Z1 <-model1[,4]/1.59

tau1 <- 0.1
Y1 <-vector(length = 365)
for (i in 1:365){
 
  Y1[i]<- rlnorm(1,log(Z1[i]/5),tau1)
}

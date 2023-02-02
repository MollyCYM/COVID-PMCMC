rm(list=ls())
require(deSolve)
times <- 1:259
N=52196381
sigma <- runif(1, min = 0, max = 1)
## first, simulate a set of random deviates
e <- rnorm(n = length(times) - 1, sd = sqrt(0.01))
## now compute their cumulative sum
e <- c(0, cumsum(e))

beta<- exp(sigma*e)
write.csv(beta,"simulatebeta.csv")


# H1N1 <- function(time, current_state, params){
# 
#   with(as.list(c(current_state, params)),{
#     for(t in 1:times){
#     N <- S+E+I+R
#     dS <- -(beta[t]*S*I)/N
#     dE <- (beta[t]*S*I)/N - E/k
#     dI <- E/k - I/gamma
#     dR <- I/gamma
#    }
#     return(list(c(dS, dE, dI, dR)))
#   })
# }
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
model1 <- ode(initial_state, times, H1N1, params)

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
write.csv(Y1,"simY11.csv")
write.csv(model1,"simulatestates1.csv")



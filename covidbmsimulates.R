rm(list=ls())
require(deSolve)
set.seed(001213)   #random seed 0012 00122 00124 00129 001211 001213 1215 1216 1219 1222
times <- 1:365     #1224
N=52196381
sigma <- 0.4
## first, simulate a set of random deviates
e <- rnorm(n = length(times) - 1, sd = sqrt(0.01))
## now compute their cumulative sum
e <- c(0, cumsum(e))


beta<- exp(sigma*e)
plot(beta,type='l')
write.csv(beta,"simulatecovidbeta.csv")


#Main ODE Model
Covid_OU <- function(time, current_state, params){
  
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
params <- c(k=7, gamma=5)   #Trial: k=5 & gamma=5; k=7 & gamma=5

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
model4 <- ode(initial_state, times, Covid_OU, params)

summary(model4)

matplot(model4, type="l", lty=1, main="COVID—BM Model",ylab="counts", xlab="Time")
legend <- colnames(model4)[3:6]
legend("right", legend=legend, col=3:6, lty = 1)
Z1 <-model4[,4]/7

tau1 <- 0.8#runif(1,0,1)
Y1 <-vector(length = 365)
for (i in 1:365){
  #Y1[i]<- rlnorm(1,log(Z1[i,]/5),tau1)
  Y1[i]<- rlnorm(1,log(Z1[i]/5),tau1)
}
plot(Y1,type='l',ylab="Obs Y",xlab="time",col="blue")
write.csv(Y1,"Covidbm_Y1.csv")
write.csv(model4,"Covidbm1.csv")
######################################################################################
#Resimulate Covid-BM data with fixed initial state values
rm(list=ls())
require(deSolve)
set.seed(001213)   #random seed 0012 00122 00124 00129 001211 001213 1215 1216 1219 1222
times <- 1:365     #1224
N=52196381
sigma <- sqrt(0.004)
## first, simulate a set of random deviates
e <- rnorm(n = length(times) - 1, sd = sqrt(0.01))
## now compute their cumulative sum
e <- c(0, cumsum(e))


beta<- exp(sigma*e)
plot(beta,type='l')
# write.csv(beta,"simulatecovidbeta1.csv")


#Main ODE Model
Covid_OU <- function(time, current_state, params){
  
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
params <- c(k=5, gamma=5)   #Trial: k=5 & gamma=5; k=7 & gamma=5

S <- N-1
E <- 1
R <- 0
I <- 0

initial_state<- c(t=1, S=S, E=E, I=I, R=R)
#initial_state <- c(t=1,S=52196380, E=1, I=0, R=0)
model11 <- ode(initial_state, times, Covid_OU, params)

summary(model11)

matplot(model11, type="l", lty=1, main="COVID—BM Model",ylab="counts", xlab="Time")
legend <- colnames(model11)[3:6]
legend("right", legend=legend, col=3:6, lty = 1)
Z11 <-model11[,4]/5

tau11 <- 0.8#runif(1,0,1)
Y11 <-vector(length = 365)
for (i in 1:365){
  #Y1[i]<- rlnorm(1,log(Z1[i,]/5),tau1)
  Y11[i]<- rlnorm(1,log(Z11[i]/5),tau11)
}
plot(Y11,type='l',ylab="Obs Y",xlab="time",col="blue")
# write.csv(Y11,"Covidbm_Y11.csv")
# write.csv(model11,"Covidbm11.csv")
######################################################################################
#Resimulate Covid-BM data with replaced sigma value
rm(list=ls())
require(deSolve)
set.seed(001213)   #random seed 0012 00122 00124 00129 001211 001213 1215 1216 1219 1222
times <- 1:365     #1224
N=52196381
sigma <- 0.4
## first, simulate a set of random deviates
e <- rnorm(n = length(times) - 1, sd = sqrt(0.01))
## now compute their cumulative sum
e <- c(0, cumsum(e))


beta<- exp(sigma*e)
plot(beta,type='l')
#write.csv(beta,"simulatecovidbeta2.csv")


#Main ODE Model
Covid_OU <- function(time, current_state, params){
  
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
params <- c(k=5, gamma=5)   #Trial: k=5 & gamma=5; k=7 & gamma=5

S <- N-1
E <- 1
R <- 0
I <- 0

initial_state<- c(t=1, S=S, E=E, I=I, R=R)
#initial_state <- c(t=1,S=52196380, E=1, I=0, R=0)
model12 <- ode(initial_state, times, Covid_OU, params)

summary(model12)

matplot(model12, type="l", lty=1, main="COVID—BM Model",ylab="counts", xlab="Time")
legend <- colnames(model12)[3:6]
legend("right", legend=legend, col=3:6, lty = 1)
Z12 <-model12[,4]/5

tau12 <- 0.8#runif(1,0,1)
Y12 <-vector(length = 365)
for (i in 1:365){
  #Y1[i]<- rlnorm(1,log(Z1[i,]/5),tau1)
  Y12[i]<- rlnorm(1,log(Z12[i]/5),tau12)
}
plot(Y12,type='l',ylab="Obs Y",xlab="time",col="blue")
# write.csv(Y12,"Covidbm_Y12.csv")
# write.csv(model12,"Covidbm12.csv")

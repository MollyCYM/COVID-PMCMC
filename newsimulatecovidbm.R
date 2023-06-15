rm(list=ls())
require(deSolve)
set.seed(1000020)   
#good seed 1000012 18 1000020
#general seed: 1000002 4 5 7 8 1000014 17
#question seed: beta is peaked on 200th day but pandemic die down already around 60th day
times <- 1:365     
N=52196381
sigma <- 0.07
## first, simulate a set of random deviates
e <- rnorm(n = length(times)-1, sd = 1)
## now compute their cumulative sum
e <- c(0, cumsum(e))
# write.csv(e,"rsimcovidbm_e1.csv")

beta<- exp(sigma*e)
plot(beta,type='l',)
plot(beta,type='p',ylab="Transmission rate beta",xlab="time",col="darkred",main="Covid_BM model generated transmission rate")
lines(beta,col="darkblue")
# write.csv(beta,"rsimcovidbm_beta1.csv")


#Main ODE Model
Covid_OU <- function(time, current_state, params){
  
  with(as.list(c(current_state, params)),{
    
    # N <- S+E+I+R
    dt <- 1
    dS <- -beta[t]*S*(E+0.1*I)/N
    dE <- beta[t]*S*(E+0.1*I)/N - E*(1/k+1/gamma)
    dI <- E/k-I*(1/gamma+0.0087)
    dR <- (I+E)/gamma+0.0087*I
    
    return(list(c( dt, dS, dE, dI, dR)))
  })
}
params <- c(k=2, gamma=9)   #Trial: k=5 & gamma=5; k=7 & gamma=5

S <- N-1
E <- 1
R <- 0
I <- 0

initial_state<- c(t=1, S=S, E=E, I=I, R=R)
model11 <- ode(initial_state, times, Covid_OU, params)

summary(model11)

matplot(model11, type="l", lty=1, main="Generated COVID_BM Model Trajectories",ylab="Counts", xlab="Time")
legend <- colnames(model11)[3:6]
legend("right", legend=legend, col=3:6, lty = 1)
Z11 <-model11[,4]/2

tau11 <- 0.1
Y11 <-vector(length = 365)
for (i in 1:365){
  Y11[i]<- rlnorm(1,log(Z11[i]/5),tau11)
}
plot(Y11,type='p',ylab="Daily incidence Y",xlab="time",col="darkred",main="Covid_BM model generated observation")
lines(Y11,col="darkblue")

# write.csv(Y11,"rsimcovidbm_Y1.csv")
# write.csv(model11,"rsimcovidbm_model1.csv")
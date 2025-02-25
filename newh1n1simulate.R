rm(list=ls())
require(deSolve)
set.seed(00001) #00033 #00036 #00043 #00048 previous:00033 #remove N= 00093

# rn<-rnorm(n=1000,m=0,sd=1)
# plot(rn, type=”l”)
# hist(rn,breaks=25)
# dp<-cumsum(rn)
# plot(dp,type=”l”)

times <- 1:365
N=52196381
sigma <- 0.07
## first, simulate a set of random deviates
e <- rnorm(n = length(times)-1 , sd = 1)
## now compute their cumulative sum
e <- c(0, cumsum(e))
plot(e,type='l')

# write.csv(e,"simulateh1n1x2.csv")
beta<- exp(sigma*e)
plot(beta,type='l')
plot(beta,type='p',ylab="Transmission rate beta",xlab="time",col="darkred",main="H1N1_BM model generated transmission rate")
lines(beta,col="darkblue")
# write.csv(beta,"simulateh1n1beta2.csv")

H1N1 <- function(time, current_state, params){
  
  with(as.list(c(current_state, params)),{
    
    # N <- S+E+I+R
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

matplot(model1, type="l", lty=1, main="Generated H1N1_BM Model Trajectories", xlab="Time",ylab = "Counts")
legend <- colnames(model1)[3:6]
legend("right", legend=legend, col=3:6, lty = 1)

Z1 <-model1[,4]/1.59

tau1 <- 0.1#runif(1,0,1)
Y1 <-vector(length = 365)
for (i in 1:365){
  #Y1[i]<- rlnorm(1,log(Z1[i,]/5),tau1)
  Y1[i]<- rlnorm(1,log(Z1[i]/5),tau1)
}
plot(Y1,type='l',ylab="Y obs", xlab = "Time")
plot(Y1,type='p',ylab="Daily incidence Y",xlab="time",col="darkred",main="H1N1_BM model generated observation")
lines(Y1,col="darkblue")
# write.csv(Y1,"simh1n1Y2.csv")
# write.csv(model1,"simulateh1n1states2.csv")
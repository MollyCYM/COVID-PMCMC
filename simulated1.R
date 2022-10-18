install.packages('deSolve')
library('deSolve')
rm(list=ls())
N = 55977178
times <- seq(1,420,by=1)
sig2 <- 0.01
sigma <- runif(1, min = 0, max = 1)
## first, simulate a set of random deviates
x <- rnorm(n = length(times) - 1, sd = sqrt(sig2))
## now compute their cumulative sum
x <- c(0, cumsum(x))
beta<- exp(sigma*x)

COVIDmod <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    for(t in 1:times){
    dS = -beta[t]*S*(0.1*I+E)/N
    dE = beta[t]*S*(0.1*I+E)/N - E*(1/k+1/gamma)
    dI = E/k-I*(1/gamma+0.0087)
    dR = (I+E)/gamma
    dZ = E/k}
    res <- c(dS, dE, dI, dR, dZ)
    list(res)
   # y <- rpois(1, Z) #????? 是1吗 还是total time steps
    
    #return(y)
  })
}

pars <- c(k = 6,   # /day: onset symptoms development period
          gamma = 14 #/day: recovery period
)

#x0<- runif(1,-5,2)
#install.packages('truncnorm')
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
Z <- 0

xstart<- c(S=S, E=E, I=I, R=R, Z=Z) 
#????都是0可以吗 尤其是I

out <-ode(xstart, times, COVIDmod,pars)
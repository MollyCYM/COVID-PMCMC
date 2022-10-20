


N = 52196381
times <- seq(1,259,by=1)
sig2 <- 0.01
sigma <- runif(1, min = 0, max = 1)
## first, simulate a set of random deviates
e <- rnorm(n = length(times) - 1, sd = sqrt(sig2))
## now compute their cumulative sum
e <- c(0, cumsum(e))
#beta<- exp(sigma*x)

H1N1mod <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    for(t in 1:times){
      dx = sigma*e[t]
      dS = -exp(x)*S*I/N
      dE = exp(x)*S*I/N - E/k
      dI = E/k-I/gamma
      dR = I/gamma
      dZ = E/k}
    res <- c(dx, dS, dE, dI, dR, dZ)
    list(res)
    # y <- rpois(1, Z) #????? 是1吗 还是total time steps
    
    #return(y)
  })
}
pars <- c(k = 1.6,   # /day: onset symptoms development period
          gamma = 1.1 #/day: recovery period
)
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
x0 <-runif(1,-5,2)
x=x0
xstart<- c(x=x, S=S, E=E, I=I, R=R, Z=Z) 
#????都是0可以吗 尤其是I

out <-ode(xstart, times, H1N1mod,pars)
write.csv(out,"simulatedH1N1.csv")
ZtransferH1N1 <- read.csv("simulatedH1N1Z.csv", header=FALSE, stringsAsFactors=FALSE) 
Ztransfer <- data.frame(ZtransferH1N1) 

tau <- runif(1,0,1)
y <-vector(length = 259)
for (i in 1:259){
  y[i]<- rlnorm(1,log(Ztransfer[i,]/10),tau)
}
write.csv(y,"simulatedyH1N1.csv")
v <- read.csv("simulatedyH1N1.csv", header=FALSE, stringsAsFactors=FALSE) 
v<-v[2:258,2] 
y <- data.frame(value = v) %>%
  mutate(time = seq(1,259, by = 1)) %>%
  dplyr::select(time, value)


plot(y,type='l')
###########################################

F<- rnorm(100000, -3,tau)
F=exp(F)
G <-rlnorm(100000,-3,tau)
hist(F,col='blue')
lines(G,col='red')
b <- min(c(G,F)) - 0.001 # Set the minimum for the breakpoints
e <- max(c(G,F)) # Set the maximum for the breakpoints
ax <- pretty(b:e, n = 12) # Make a neat vector for the breakpoints
hgF<- hist(F, plot = FALSE) # Save first histogram data
hgG <- hist(G, plot = FALSE) # Save 2nd histogram data

plot(hgF, col = 'red') # Plot 1st histogram using a transparent color
plot(hgG, col = 'blue', add = TRUE) # Add 2nd histogram using different color


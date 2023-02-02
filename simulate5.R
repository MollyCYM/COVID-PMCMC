rm(list=ls())
N=1000000
require(deSolve)
SEIR <- function(time, current_state, params){
  
  with(as.list(c(current_state, params)),{
    N <- S+E+I+R
    dS <- -(beta*S*I)/N
    dE <- (beta*S*I)/N - sigma*E
    dI <- sigma*E - gamma*I - mu*I
    dR <- gamma*I
    dM <- mu*I
    
    return(list(c(dS, dE, dI, dR, dM)))
  })
}
params <- c(beta=0.5, sigma=0.25, gamma=0.2, mu=0.001)
initial_state <- c(S=999999, E=1, I=0, R=0, M=0)
times <- 0:1000
model2 <- ode(initial_state, times, SEIR, params)

summary(model2)

matplot(model2, type="l", lty=1, main="SEIR model", xlab="Time")
legend <- colnames(model2)[2:6]
legend("right", legend=legend, col=2:6, lty = 1)


Z <-0.2*model2[,3]
S<-model2[,2]
I<-model2[,4]
R<-model2[,5]
M<-model2[,6]


tau <- runif(1,0,1)
Y <-vector(length = 1001)
for (i in 1:1001){
  Y[i]<- rlnorm(1,log(Z[i]/5),tau)
}
plot(Y,type='l')
write.csv(Y,"simYl.csv")
write.csv(model2,"simulatestatesl.csv")


require(deSolve)

SEIR <- function(time, current_state, params){
  
  with(as.list(c(current_state, params)),{
    N <- S+E+I+R
    dS <- -(beta*S*I)/N
    dE <- (beta*S*I)/N - sigma*E
    dI <- sigma*E - gamma*I - mu*I
    dR <- gamma*I+mu*I
    
    return(list(c(dS, dE, dI, dR)))
  })
}
params <- c(beta=0.5, sigma=0.25, gamma=0.2, mu=0.001)
initial_state <- c(S=999999, E=1, I=0, R=0)
times <- 0:365
model <- ode(initial_state, times, SEIR, params)

summary(model)

matplot(model, type="l", lty=1, main="SEIR model", xlab="Time")
legend <- colnames(model)[2:5]
legend("right", legend=legend, col=2:5, lty = 1)
write.csv(model,"simulatestates4.csv")


Z <-0.2*model[,3]




tau <- 0.8
Y <-vector(length = 366)
for (i in 1:366){
  Y[i]<- rlnorm(1,log(Z[i]/5),tau)
}
plot(Y,type='l')
write.csv(Y,"simY4.csv")

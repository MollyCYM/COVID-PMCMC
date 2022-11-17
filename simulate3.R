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
times <- 0:365
model <- ode(initial_state, times, SEIR, params)

summary(model)

matplot(model, type="l", lty=1, main="SEIR model", xlab="Time")
legend <- colnames(model)[2:6]
legend("right", legend=legend, col=2:6, lty = 1)

Z <-as.vector(model[,3])
Z <-0.2*Z
Z <-0.2*model[,3]
S<-model[,2]
I<-model[,4]
R<-model[,5]
M<-model[,6]

Mmodel <- read.csv("simulatestates.csv", header=TRUE, stringsAsFactors=FALSE)
S<-Mmodel[,3]
E<-Mmodel[,4]
I<-Mmodel[,5]
R<-Mmodel[,6]
M<-Mmodel[,7]
plot(Z/5,type='l',col='Red',main="Simulation states vs Model latent state estimations",ylab=TeX("State Z=($\\sigma E/5$)"))
lines(Zmodel,type='l',col='Blue')
plot(S,type='l',col='Red',main="Simulation states vs Model latent state estimations",ylab=TeX("State S"))
lines(Savgdata,type='l',col='Blue')
plot(I,type='l',col='Red',main="Simulation states vs Model latent state estimations",ylab=TeX("State I"))
lines(Iavgdata,type='l',col='Blue')
plot(R,type='l',col='Red',main="Simulation states vs Model latent state estimations",ylab=TeX("State R"))
lines(Ravgdata,type='l',col='Blue')
plot(M,type='l',col='Red',main="Simulation states vs Model latent state estimations",ylab=TeX("State M"))
lines(Mavgdata,type='l',col='Blue')
plot()
write.csv(Z,"simZ1.csv")
Z <- read.csv("simz_1.csv", header=FALSE, stringsAsFactors=FALSE) 
Z <- data.frame(Z) 

tau <- runif(1,0,1)
Y <-vector(length = 366)
for (i in 1:366){
  Y[i]<- rlnorm(1,log(Z[i,]/5),tau)
}
plot(Y,type='l')
write.csv(Y,"simY1.csv")
write.csv(model,"simulatestates.csv")


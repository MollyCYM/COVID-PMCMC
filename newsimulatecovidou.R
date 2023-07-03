rm(list=ls())
require(deSolve)
set.seed(10023) #00033 #00036 #00043 #00048 previous:00033 #remove N= 00093
#Good: 00001 00033 00048 10000 10007 10012-large outbreak 10014 10016 10018 10019 10020
#Pick 10023
#Gneneral: 00036 00043 10001 10003 10006 10009 10011 10013 10017
#Small outbreak: 10002 10010 
times <- 1:365
N=52196381
## first, simulate a set of random deviates
e <- rnorm(n = length(times) , sd = 1)
# write.csv(e,"rsimcovidou_e1.csv")
plot(e,type='l')
mu<-vector(length=365) 
for (t in 1:365) {
  if (t <= 120){
    mu[t]=-0.02      }                    #No lock-down policy
  else {mu[t]=-0.2-0.02 }                 #Lock-down policy
}
mu<-ts(mu)
# write.csv(mu,"rsimcovidou_mu1.csv")
plot(mu,type='l')

#Simulate an O-U Process
ornstein_uhlenbeck <- function(n,theta,sigma,x0){
  dt  <- 1
  x<-vector(length=365)
  for (i in 1:n) {
    if (i==1){x[i]=x0}
    else{
      x[i]  <-  x[i-1] + theta*(mu[i]-x[i])*dt + sigma*e[i]}
      # x[i]  <-  x[i-1] + sigma*e[i]}
  }
  return(x);
}

x<- ornstein_uhlenbeck(365,0.05,sqrt(0.004),0)
plot(x,type='l',xlab="time")
# abline(v=121, col="red")
# lines(mu, col="blue")
# tt1 <-expression(mu==-0.02)
# text(30,-0.1,tt1,col="green")
# tt2 <-expression(mu==-0.02-0.2)
# text(210,-0.3,tt2,col="green")
plot(exp(x),type='l',xlab="time",ylab=TeX("$\\e^{x}=beta$"))
abline(v=121, col="red")
lines(exp(mu+0.02),col="blue")
axis(1, at=121,labels=121, col.axis="red", las=2)
tt1 <-expression(mu[beta]==1)
text(50,1.1,tt1,col="green")
tt2 <-expression(mu[beta]==0.81)
text(220,0.9,tt2,col="green")
tt3 <-expression('lockdown policy start')
text(120,0.4,tt3,col="red")
# write.csv(x,"rsimcovidou_x1.csv")
#Main ODE Model
Covid_OU <- function(time, current_state, params){
  
  with(as.list(c(current_state, params)),{
    
    # N <- S+E+I+R
    dt <- 1
    dS <- -exp(x[t])*S*(E+0.1*I)/N
    dE <- exp(x[t])*S*(E+0.1*I)/N - E*(1/k+1/gamma)
    dI <- E/k-I*(1/gamma+0.0087)
    dR <- (I+E)/gamma+0.0087*I
    
    return(list(c( dt, dS, dE, dI, dR)))
  })
}
params <- c(k=2, gamma=9)   #Trial: k=5 & gamma=5; k=7 & gamma=5

# library('truncnorm')
# R0 <- rtruncnorm(1, a=0, b=1, mean = 0.15, sd = 0.15)
# E0 <-runif(1,-16, -9)
# I0 <-runif(1,-16, -9)
# 
# S <- N
# R <- R0*S
# S <- S - R
# 
# E <- exp(E0 + log(S))
# S <- S - E
# I <- exp(I0 + log(S))
# S <- S - I

# initial_state<- c(t=1, S=S, E=E, I=I, R=R)
initial_state <- c(t=1,S=N-1, E=1, I=0, R=0)
model4 <- ode(initial_state, times, Covid_OU, params)

summary(model4)


matplot(model4, type="l", lty=1, main="Generated COVID_OU Model Trajectories",ylab="Counts", xlab="Time")
legend <- colnames(model4)[3:6]
legend("right", legend=legend, col=3:6, lty = 1)
Z1 <-model4[,4]/2

tau1 <- 0.1
Y1 <-vector(length = 365)
for (i in 1:365){
  #Y1[i]<- rlnorm(1,log(Z1[i,]/5),tau1)
  Y1[i]<- rlnorm(1,log(Z1[i]/5),tau1)
}
plot(Y1,type='l',ylab="Obs Y",xlab="time",col="blue")
plot(Y1,type='p',ylab="Daily incidence Y",xlab="time",col="darkred",main="Covid_OU model generated observation")
lines(Y1,col="darkblue")
# write.csv(Y1,"rsimcovidou_Y1.csv")
# write.csv(model4,"rsimcovidou_model1.csv")
